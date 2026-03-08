# ==============================================================================
# Clean SSR Health Files with GTN ratios and make annual estimates
# (no name matching yet)
# msahu
# ==============================================================================

source('init.R')

# ==============================================================================
# Read and combine all files
# ==============================================================================

ssr_dir <- paste0(data_dir, 'raw/ssr_health/')
ssr_files <- list.files(ssr_dir, pattern = "\\.csv$", full.names = TRUE)

ssr_dt <- rbindlist(
  lapply(ssr_files, function(f) {
    dt <- as.data.table(read.delim(f, fileEncoding = "UTF-16LE"))
    setnames(
      dt,
      old = c("X", "X.1", "X.2"),
      new = c("product", "disease_area", "original_smoothing") # Rename to avoid confusion
    )
    dt[, source_file := basename(f)]
    # Extract payer and smoothing type from filename
    dt[, `:=`(
      payer = sub(
        "GTN_(4qMA_|Unsmoothed_)",
        "",
        sub("\\.csv$", "", source_file)
      ),
      smoothing_type = ifelse(
        grepl("4qMA", source_file),
        "4q moving average",
        "Unsmoothed"
      )
    )]
    dt
  }),
  use.names = T
)

# ==============================================================================
# Clean and reshape
# ==============================================================================

# Cleaning
ssr_dt[ssr_dt == ""] <- NA
ssr_dt[, product_clean := toupper(trimws(product))]
ssr_dt[, original_smoothing := NULL]

# Reshape to long
ssr_long <- melt(
  ssr_dt,
  id.vars = c(
    "product",
    "product_clean",
    "disease_area",
    "smoothing_type",
    "payer",
    "source_file"
  ),
  variable.name = "quarter",
  value.name = "gtn_scalar"
)

# Clean up the quarter column (remove "X" prefix)
ssr_long[, quarter := as.character(gsub("^X", "", quarter))]

# Extract year and quarter
ssr_long[, `:=`(
  year_id = as.integer(substr(quarter, 1, 4)),
  qtr = as.integer(substr(quarter, 7, 7))
)]

# Calculate mean of available quarters
ssr_annual <- ssr_long[,
  .(
    gtn_scalar_annual = mean(
      as.numeric(gsub("%", "", gtn_scalar)) / 100,
      na.rm = TRUE
    ),
    n_quarters = sum(!is.na(gtn_scalar))
  ),
  by = .(product, product_clean, disease_area, year_id, payer, smoothing_type)
]

# --- DIAGNOSTICS ---
cat("Total Rows:", nrow(ssr_annual), "\n")
cat(
  "Count of NAs (Missing years):",
  sum(is.na(ssr_annual$gtn_scalar_annual)),
  "\n"
)
cat(
  "Count of Negative Rebates (Price > WAC):",
  sum(ssr_annual$gtn_scalar_annual < 0, na.rm = T),
  "\n"
)
cat(
  "Count of Zero Rebates (Legacy Brands):",
  sum(ssr_annual$gtn_scalar_annual == 0, na.rm = T),
  "\n"
)
cat(
  "Count of Tiny Rebates (0 < x < 1%):",
  sum(
    ssr_annual$gtn_scalar_annual > 0 & ssr_annual$gtn_scalar_annual < 0.01,
    na.rm = T
  ),
  "\n"
)
cat(
  "Count of Valid Rebates (>= 1%):",
  sum(ssr_annual$gtn_scalar_annual >= 0.01, na.rm = T),
  "\n"
)

# ==============================================================================
# Save
# ==============================================================================

# Annual estimates
fwrite(
  ssr_annual,
  paste0(data_dir, "processed/ssr_cleaned/ssr_gtn_annual.csv")
)

# ==============================================================================
# Outlier removal + linear imputation
# ==============================================================================

# 1. OUTLIER REMOVAL -----------------------------------------------------------
# Logic based on diagnostics:
# - Drop NA: The 102k empty rows (drug didn't exist).
# - Drop < 0: Impossible prices.
# - Drop > 1: Accounting errors (>100% rebate).
# - KEEP == 0: Legacy brands (The 32 rows).
# - DROP (0, 0.01): The 186 rows of "noise" (likely inventory errors).
# - KEEP >= 0.01: Valid rebates.

ssr_annual_clean <- ssr_annual[
  !is.na(gtn_scalar_annual) &
    !is.nan(gtn_scalar_annual) &
    gtn_scalar_annual <= 1.00 &
    (gtn_scalar_annual == 0 | gtn_scalar_annual >= 0.01) # Keep 0 OR >1%
]

# 2. CREATE TIME GRID ----------------------------------------------------------
# We must create a full grid (Product x Year) so 'zoo' knows where the gaps are.

min_yr <- min(ssr_annual_clean$year_id)
max_yr <- max(ssr_annual_clean$year_id)
all_years <- data.table(year_id = min_yr:max_yr)

# Get unique keys from the CLEANED data
unique_prods <- unique(ssr_annual_clean[, .(
  product,
  product_clean,
  disease_area,
  payer,
  smoothing_type
)])

# Cross-join: Repeat every product for every year
template <- unique_prods[rep(1:.N, each = nrow(all_years))]
template[, year_id := rep(all_years$year_id, times = nrow(unique_prods))]

# Merge actual data onto the template
ssr_grid <- merge(
  template,
  ssr_annual_clean,
  by = c(
    "product",
    "product_clean",
    "disease_area",
    "payer",
    "smoothing_type",
    "year_id"
  ),
  all.x = TRUE
)

# 3. LINEAR INTERPOLATION ------------------------------------------------------
# maxgap = 1: Only interpolate if the gap is strictly 1 year (e.g., 2008, NA, 2010).
# rule = 1: Do NOT extrapolate. If 2008 is the first data point, 2007 remains NA.

ssr_grid[,
  gtn_imputed := na.approx(
    gtn_scalar_annual,
    maxgap = 1,
    na.rm = FALSE,
    rule = 1
  ),
  by = .(product, payer, smoothing_type)
]

# Flag imputed values for your records
ssr_grid[, is_imputed := is.na(gtn_scalar_annual) & !is.na(gtn_imputed)]

# 4. FINAL CLEANUP -------------------------------------------------------------
# Remove the empty leading/trailing years (NAs that weren't interpolated)
ssr_annual_imputed <- ssr_grid[!is.na(gtn_imputed)]

# Select final columns
ssr_annual_imputed <- ssr_annual_imputed[, .(
  product,
  product_clean,
  disease_area,
  payer,
  smoothing_type,
  year_id,
  gtn_final = gtn_imputed, # This contains original values + interpolations
  is_imputed
)]

# Final Validation
cat("Final Dataset Size:", nrow(ssr_annual_imputed), "\n")
cat("Total Imputed Years:", sum(ssr_annual_imputed$is_imputed), "\n")


# ==============================================================================
# Save
# ==============================================================================

# Annual estimates
fwrite(
  ssr_annual_imputed,
  paste0(data_dir, "processed/ssr_cleaned/ssr_gtn_annual_imputed.csv")
)

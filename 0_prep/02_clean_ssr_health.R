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

# ==============================================================================
# Save
# ==============================================================================

# Annual estimates
fwrite(
  ssr_annual,
  paste0(data_dir, "processed/ssr_cleaned/ssr_gtn_annual.csv")
)

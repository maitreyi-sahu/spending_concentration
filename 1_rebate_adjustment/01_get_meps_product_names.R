source('init.R')

# ==============================================================================
# Load MEPS data
# ==============================================================================

# Load MEPS NDCs
meps_ndcs <- fread(
  paste0(data_dir, "processed/meps_cleaned/ndc_lookup.csv"),
  colClasses = list(character = c("ndc_original", "ndc_clean"))
)

# Filter valid only (optional, but recommended for speed)
# We keep year_id to allow matching against FDB history
meps_ndcs <- meps_ndcs[ndc_quality_flag != 'missing_code'][order(-total_spend)]

# Create a unique list of NDC + Year to match against
# We create a "match_date" (July 1st) to check if the drug was active that year
meps_keys <- unique(meps_ndcs[, .(ndc_clean, year_id)])
meps_keys[, match_date := as.Date(paste0(year_id, "-07-01"))]

# ==============================================================================
# Process FDB
# ==============================================================================

# Load FDB file
fdb_path <- paste0(data_dir, "raw/fdb/master_2025_2.sas7bdat")
fdb_dt <- as.data.table(read_sas(fdb_path))
names(fdb_dt) <- toupper(names(fdb_dt))

# Filter immediately
fdb_dt[, ndc_clean := str_pad(NDC, 11, pad = "0")]
fdb_subset <- fdb_dt[ndc_clean %in% meps_ndcs$ndc_clean]

# --- HANDLE DATES ---
# SAS numeric dates to R Dates
if (is.numeric(fdb_subset$STARTDT)) {
  fdb_subset[, start_date := as.Date(STARTDT, origin = "1960-01-01")]
} else {
  fdb_subset[, start_date := as.Date(STARTDT)]
}

if (is.numeric(fdb_subset$ENDDT)) {
  fdb_subset[, end_date := as.Date(ENDDT, origin = "1960-01-01")]
} else {
  fdb_subset[, end_date := as.Date(ENDDT)]
}

# Handle active drugs
fdb_subset[is.na(end_date), end_date := as.Date("2099-12-31")]

# --- THE INTERVAL MATCH ---
# Join MEPS (Point in time) to FDB (Start/End Range)
fdb_matches <- fdb_subset[
  meps_keys,
  on = .(ndc_clean, start_date <= match_date, end_date >= match_date),
  # SELECT YOUR SPECIFIC COLUMNS HERE
  .(
    ndc_clean,
    year_id,
    BRAND,
    GENERIC,
    GI,
    CLASS,
    ATC_2025_2,
    start_date,
    end_date
  )
]

# Deduplicate (Take most recent start date if overlaps exist)
fdb_matches <- unique(
  fdb_matches,
  by = c("ndc_clean", "year_id"),
  fromLast = TRUE
)

# --- DEFINE STATUS & NAME ---
# 1. Product Name: Prefer BRAND, fall back to GENERIC if Brand is missing
fdb_matches[,
  fdb_product_name := ifelse(is.na(BRAND) | BRAND == "", GENERIC, BRAND)
]

# 2. Status: Check FDB Dictionary, but usually 1=Brand, 0=Generic
fdb_matches[,
  fdb_status := fcase(
    GI == "1" , "Brand"   ,
    GI == "0" , "Generic" ,
    default = "Unknown"
  )
]

# ==============================================================================
# Run checks
# ==============================================================================

# Check how many found a Product Name
match_count <- sum(!is.na(fdb_matches$BRAND))
total_count <- nrow(fdb_matches)

cat("Actual Match Rate:", round(match_count / total_count * 100, 2), "%\n")

# Overwrite meps_ndcs (or create a new object) to include the FDB data
meps_linked <- merge(
  meps_ndcs,
  fdb_matches,
  by = c("ndc_clean", "year_id"),
  all.x = TRUE
)

# 2. CREATE THE NAME COLUMN
# ==============================================================================
# We need to consolidate Brand and Generic names into one column to check coverage
meps_linked[,
  product_name_final := fcase(
    !is.na(BRAND) & BRAND != ""     , BRAND   , # Prefer Brand Name
    !is.na(GENERIC) & GENERIC != "" , GENERIC , # Fallback to Generic Name
    default = NA_character_
  )
]

# 3. CALCULATE SPENDING
# ==============================================================================
# Denominator: Total Spend in raw MEPS
total_spend_all <- sum(meps_linked$total_spend, na.rm = TRUE)

# Numerator: Spend associated with NDCs that successfully matched a Name
total_spend_matched <- sum(
  meps_linked[!is.na(product_name_final)]$total_spend,
  na.rm = TRUE
)

# 4. RESULTS
# ==============================================================================
cat("\n=== SPENDING COVERAGE CHECK ===\n")
cat(
  "Total Spend:    $",
  formatC(total_spend_all, format = "f", big.mark = ",", digits = 0),
  "\n"
)
cat(
  "Matched Spend:  $",
  formatC(total_spend_matched, format = "f", big.mark = ",", digits = 0),
  "\n"
)
cat(
  "Coverage Rate:  ",
  round((total_spend_matched / total_spend_all) * 100, 2),
  "%\n"
)

# 5. DIAGNOSTICS
# Check the top 10 most expensive NDCs that are MISSING a name
cat("\nTop 10 Unmatched by Spend (Likely Insulin supplies or bad codes):\n")
print(meps_linked[is.na(product_name_final)][order(-total_spend)][
  1:10,
  .(ndc_clean, year_id, total_spend, n_claims)
])

# ==============================================================================
# Process FDA (The Backup Source)
# ==============================================================================

# 3. Clean and Filter
# ------------------------------------------------------------------------------
# Standardize NDC
fdb_dt[, ndc_clean := str_pad(NDC, 11, pad = "0")]

# Filter: Only keep NDCs that exist in your MEPS dataset
# This shrinks the data from millions of rows to just your 77k NDCs
fdb_subset <- fdb_dt[ndc_clean %in% meps_ndcs$ndc_clean]

# 4. Handle Dates (The Jun Liu Method)
# ------------------------------------------------------------------------------
# Convert SAS dates (often just numbers) to R Dates if needed
# If they look like '15000', they are days since 1960.
# Check class first. If numeric:
# fdb_subset[, start_date := as.Date(STARTDT, origin = "1960-01-01")]
# fdb_subset[, end_date   := as.Date(ENDDT,   origin = "1960-01-01")]

# If they are already Date format, just rename:
setnames(fdb_subset, "STARTDT", "start_date", skip_absent = TRUE)
setnames(fdb_subset, "ENDDT", "end_date", skip_absent = TRUE)

# Handle Active Drugs (End Date often 2030 or NA)
fdb_subset[is.na(end_date), end_date := as.Date("2099-12-31")]

# 5. The Join (Interval Match)
# ------------------------------------------------------------------------------
# We need to link MEPS (NDC + Year) to FDB (NDC + Start/End Window)
# We use foverlaps or a simple non-equi join in data.table

# Create a "MEPS Date" (e.g., July 1st of the study year) for matching
meps_keys[, match_date := as.Date(paste0(year, "-07-01"))]

# Perform Non-Equi Join
matched <- fdb_subset[
  meps_keys,
  on = .(ndc_clean, start_date <= match_date, end_date >= match_date),
  nomatch = NULL, # Drops non-matches
  .(ndc_clean, year, LN, GI, GCN_SEQNO, GTC, start_date, end_date)
]

# 6. Deduplicate (If multiple rows still match)
# Sometimes ranges overlap slightly. Take the most recent start date.
final_matches <- unique(matched, by = c("ndc_clean", "year"), fromLast = TRUE)


# ==============================================================================
# Clean FDA data
# ==============================================================================

# Load FDA files
fda_product <- fread(paste0(data_dir, "raw/fda/product.txt"), sep = "\t")
fda_package <- fread(paste0(data_dir, "raw/fda/package.txt"), sep = "\t")

# Clean column names
setnames(fda_product, make.names(tolower(names(fda_product))))
setnames(fda_package, make.names(tolower(names(fda_package))))

# Standardize package NDCs to 11 digits
fda_package[, ndc_clean := gsub("-", "", ndcpackagecode)]
fda_package[, ndc_clean := sprintf("%011s", ndc_clean)]
fda_package[, ndc_clean := gsub(" ", "0", ndc_clean)]

# Merge package + product info
fda_full <- merge(
  fda_package[, .(ndc_clean, productndc, packagedescription)],
  fda_product[, .(
    productndc,
    proprietaryname,
    nonproprietaryname,
    active_numerator_strength,
    dosageformname,
    routename,
    labelername
  )],
  by = "productndc",
  all.x = TRUE
)

# ==============================================================================
# Match MEPS to FDA
# ==============================================================================

meps_fda <- merge(
  meps_ndcs,
  fda_full,
  by = "ndc_clean",
  all.x = TRUE
)

# Clean brand names
meps_fda[, brand_name := toupper(trimws(proprietaryname))]
meps_fda[, generic_name := toupper(trimws(nonproprietaryname))]
meps_fda[, ingredient := toupper(trimws(active_numerator_strength))]

# Brand vs Generic flag
meps_fda[,
  brand_generic := fcase(
    !is.na(proprietaryname) & proprietaryname != ""       , "Brand"   ,
    !is.na(nonproprietaryname) & nonproprietaryname != "" , "Generic" ,
    default = NA_character_
  )
]

# Extract product name for matching to SSR (use brand if available, else generic)
meps_fda[,
  product_name := fcase(
    !is.na(brand_name) & brand_name != ""     , brand_name   ,
    !is.na(generic_name) & generic_name != "" , generic_name ,
    default = NA_character_
  )
]

# ==============================================================================
# Save
# ==============================================================================

saveRDS(
  meps_fda,
  paste0(data_dir, "processed/meps_cleaned/meps_with_fda_names.rds")
)

# ==============================================================================
# Summary
# ==============================================================================

cat("\n=== FDA MATCHING SUMMARY ===\n")
cat("Total MEPS NDCs:", nrow(meps_fda), "\n")
cat(
  "Matched to FDA:",
  sum(!is.na(meps_fda$product_name)),
  "(",
  round(100 * sum(!is.na(meps_fda$product_name)) / nrow(meps_fda), 1),
  "%)\n"
)
cat(
  "  Brand drugs:",
  sum(meps_fda$brand_generic == "Brand", na.rm = TRUE),
  "\n"
)
cat(
  "  Generic drugs:",
  sum(meps_fda$brand_generic == "Generic", na.rm = TRUE),
  "\n"
)

# Spending coverage
matched_spend <- sum(meps_fda[!is.na(product_name)]$total_spend, na.rm = TRUE)
total_spend <- sum(meps_fda$total_spend, na.rm = TRUE)
cat("\nSpending coverage:", round(100 * matched_spend / total_spend, 1), "%\n")

# Top unmatched by spending
cat("\nTop 10 unmatched NDCs by spending:\n")
print(meps_fda[is.na(product_name)][order(-total_spend)][
  1:10,
  .(ndc_clean, total_spend, n_claims)
])

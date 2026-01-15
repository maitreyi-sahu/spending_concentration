# ==============================================================================
# MEPS NDC CLEANING
#
# Problem -- ndcs missing for 21.66% of spending -- Haley to try to keep RXNAME
# TODO create cleaned NDC file
#
# ==============================================================================

source('init.R')

# Load in MEPS data, which has NDC codes only
dt <- open_dataset(paste0(data_dir, "raw/meps/USA_MEPS_RX.parquet")) |>
  filter(toc == 'RX' & year_id >= START_YEAR) |>
  select(year_id, ndc, tot_pay_amt) |>
  as.data.table() |>
  collect()

# ==============================================================================
# Clean NDCs
# ==============================================================================

# Identify MEPS missing data codes (negative values)
dt[,
  is_missing := ifelse(is.na(ndc) | grepl("^-[0-9]+$", ndc) | ndc == "", T, F)
]

# Extract only numeric characters
dt[, ndc_numeric := gsub("[^0-9]", "", ndc)]
dt[, numeric_length := nchar(ndc_numeric)]

# Initialize cleaned NDC column
dt[, ndc_clean := NA_character_]

# For valid numeric NDCs: pad to 11 digits with leading zeros
# Use sprintf with %011s to pad with spaces, then replace spaces with zeros
dt[!is_missing & numeric_length > 0, ndc_clean := sprintf("%011s", ndc_numeric)]
dt[!is.na(ndc_clean), ndc_clean := gsub(" ", "0", ndc_clean)]

# Mark invalid NDCs (too long or unreasonably short)
dt[numeric_length > 11 & !is_missing, ndc_clean := NA_character_]
dt[numeric_length < 6 & !is_missing, ndc_clean := NA_character_]

# Create quality flags
dt[,
  suspicious_leading_zeros := grepl("^0{5,}", ndc_clean) &
    !is.na(ndc_clean) &
    !is_missing
]
dt[, was_modified := (ndc != ndc_clean) & !is.na(ndc_clean)]

dt[,
  ndc_quality_flag := fcase(
    is_missing                       , "missing_code"             ,
    numeric_length > 11              , "too_long_invalid"         ,
    numeric_length < 6 & !is_missing , "too_short_invalid"        ,
    grepl("^0+$", ndc_clean)         , "all_zeros"                ,
    suspicious_leading_zeros         , "suspicious_leading_zeros" ,
    !is.na(ndc_clean)                , "valid"                    ,
    default = "other_invalid"
  )
]
# ==============================================================================
# SUMMARIES
# ==============================================================================

cat("\n=== NDC CLEANING SUMMARY ===\n\n")
cat("Total NDCs processed:", nrow(dt), "\n\n")

cat("Quality Flag Distribution:\n")
print(dt[, .N, by = ndc_quality_flag][order(-N)])

# Overall quality distribution
cat("Quality Flag Distribution:\n")
quality_summary <- dt[,
  .(
    n_claims = .N,
    pct_claims = round(.N / nrow(dt) * 100, 2),
    total_spend = sum(tot_pay_amt, na.rm = TRUE),
    pct_spend = round(
      sum(tot_pay_amt, na.rm = TRUE) / sum(dt$tot_pay_amt, na.rm = TRUE) * 100,
      2
    )
  ),
  by = ndc_quality_flag
][order(-total_spend)]
print(quality_summary)

# ==============================================================================
# CREATE MEPS NDC map -- save in the processed dir
# ==============================================================================

# Create NDC lookup with unique NDCs and their quality flags
ndc_lookup <- dt[,
  .(
    n_claims = .N,
    total_spend = sum(tot_pay_amt, na.rm = TRUE),
    ndc_quality_flag = first(ndc_quality_flag),
    was_modified = first(was_modified),
    suspicious_leading_zeros = first(suspicious_leading_zeros)
  ),
  by = .(year_id, ndc_original = ndc, ndc_clean)
]

# Before saving ndc_lookup, ensure character type:
ndc_lookup[, ndc_original := as.character(ndc_original)]
ndc_lookup[, ndc_clean := as.character(ndc_clean)]

# Save the lookup
fwrite(ndc_lookup, paste0(data_dir, "processed/meps_cleaned/ndc_lookup.csv"))

# Save cleaned MEPS prescription data
# Keep only essential columns for matching
meps_prescriptions_clean <- dt[, .(
  year_id,
  ndc_original = ndc,
  ndc_clean,
  ndc_quality_flag,
  tot_pay_amt
)]

saveRDS(
  meps_prescriptions_clean,
  paste0(data_dir, "processed/meps_cleaned/meps_prescriptions_clean.rds")
)

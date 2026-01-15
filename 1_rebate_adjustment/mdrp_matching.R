# ==============================================================================
# MERGE STRATEGY: Local CSV (High Speed)
# ==============================================================================

library(data.table)
library(stringr)

# 1. Load your MEPS NDCs
# ------------------------------------------------------------------------------
# Assuming 'meps_ndcs' is already loaded from your previous code
# Ensure you have a clean 11-digit column
meps_ndcs[, ndc_clean := gsub("-", "", as.character(ndc_original))]
meps_ndcs[, ndc_clean := str_pad(ndc_clean, 11, pad = "0")]

# 2. Load the Medicaid CSV
# ------------------------------------------------------------------------------
# Medicaid Drug Rebate Program Data
medicaid_path <- paste0(
  data_dir,
  "raw/mdrp/drugproducts3q_2025Updated11132025.csv"
)

# Read all columns as character to preserve leading zeros
cms_data <- fread(medicaid_path, colClasses = "character")

# Clean column names (make lower case for easier reference)
names(cms_data) <- tolower(gsub(" ", "_", names(cms_data)))

# 3. Construct the 11-digit NDC from the Medicaid file
# ------------------------------------------------------------------------------
# The file has: labeler_code (5), product_code (4), package_size (2)
# We must pad them to ensure 00123-0004-01 format, then join.

cms_data[,
  ndc11 := paste0(
    str_pad(labeler_code, 5, pad = "0"),
    str_pad(product_code, 4, pad = "0"),
    str_pad(package_size_code, 2, pad = "0")
  )
]

# 4. Create the Lookup Table
# ------------------------------------------------------------------------------
# Keep the most relevant columns for your research
cms_lookup <- cms_data[, .(
  ndc11,
  fda_product_name, # Marketed Name
  drug_category, # The gold standard for Brand/Generic
  market_date # Helps verify if it was active 2000-2022
)]

# Deduplicate if necessary (rare, but good practice)
cms_lookup <- unique(cms_lookup, by = "ndc11")

# 5. Join
# ------------------------------------------------------------------------------
final_data <- merge(
  meps_ndcs,
  cms_lookup,
  by.x = "ndc_clean",
  by.y = "ndc11",
  all.x = TRUE
)

# 6. Interpret the Results
# ------------------------------------------------------------------------------
# N = Non-Innovator (Generic)
# I = Innovator (Brand)
# S = Single Source (Brand)

final_data[,
  brand_generic_status := fcase(
    drug_category == "N"           , "Generic" ,
    drug_category %in% c("I", "S") , "Brand"   ,
    default = "Unknown/Unmatched"
  )
]

# Check Coverage
print(final_data[, .N, by = brand_generic_status])

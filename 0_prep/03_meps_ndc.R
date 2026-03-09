source('init.R')

# ==============================================================================
# STEP 1: PREPARE FDB DICTIONARY (ndc to drug name / is_generic)
# ==============================================================================
# Load FDB Product/Name files
# You need columns: NDC, LABEL_NAME (Brand), GENERIC_NAME (Ingredient), GENERIC_INDICATOR

# Load data
fdb <- read_sas(paste0(data_dir, "raw/firstdatabank/master_2025_2.sas7bdat")) |>
  setDT()

# Clean column names
fdb[, `:=`(
  fdb_brand = toupper(trimws(Brand)),
  fdb_generic = toupper(trimws(Generic)),
  ndc = NDC
)]

# Create GENERIC INDICATOR using the GI column
fdb[,
  fdb_genind := fcase(
    fdb_brand == fdb_generic , "Generic" , # Explicit string match
    as.character(GI) == "1"  , "Generic" , # 1 = Multi-source (Generic)
    as.character(GI) == "2"  , "Brand"   , # 2 = Single-source (Brand)
    default = NA_character_ # Catch-all for any missing or unexpected values
  )
]

# DEDUPLICATE: Prioritize Brand if conflict exists
setorder(fdb, ndc, fdb_genind, na.last = TRUE)
fdb <- unique(fdb, by = "ndc")
fdb <- fdb[, .(ndc, fdb_brand, fdb_generic, fdb_genind)]

# BUILD NAME DICTIONARY (Crucial for the text matching step)
# This maps "LIPITOR" -> "Brand" and "ATORVASTATIN" -> "Generic"
dict_brands <- data.table(
  name = unique(fdb[fdb_genind == "Brand", fdb_brand]),
  type = "Brand"
)
dict_generics <- data.table(
  name = unique(fdb[fdb_genind == "Generic", fdb_generic]),
  type = "Generic"
)
name_dict <- rbind(dict_brands, dict_generics)
name_dict <- unique(name_dict, by = "name") # remove dups

# ==============================================================================
# STEP 2: PREPARE RXNORM (Rescue Lookup) & CLEANING FUNCTION
# ==============================================================================

# ENHANCED cleaning function for airtight text matching
clean_drug_text <- function(x) {
  x <- toupper(as.character(x))

  # 1. Remove anything in parentheses (often packaging/manufacturer noise)
  x <- str_remove_all(x, "\\(.*?\\)")

  # 2. Remove common MEPS survey temporal/dosing noise
  x <- str_remove_all(x, "\\b\\d+\\s*DAY\\b")
  x <- str_remove_all(x, "\\b(REFORMULATED|FIRST MONTH|STARTER|TITRATION).*")

  # 3. Remove numerical dosages, ratios, and percentages
  x <- str_remove_all(x, "\\b\\d+\\.?\\d*\\s*\\/\\s*\\d+\\.?\\d*\\b") # e.g. 50/50
  x <- str_remove_all(x, "\\b\\d+\\.?\\d*\\s*(MG|MCG|ML|G|L|MEQ|UNIT[S]?|%)\\b")

  # 4. Remove salts, forms, and delivery device suffixes
  # (Crucial for linking MEPS text to financial filings)
  forms_salts_devices <- "\\b(TAB|CAP|INJ|SOL|HCL|SODIUM|EQ|ER|XR|DR|SYR|SUSP|PEN|FLEXPEN|KWIKPEN|RESPIMAT|HANDIHALER|DISKUS|HFA|AUTOINJECTOR|VIAL|SOLOSTAR)\\b"
  x <- str_remove_all(x, forms_salts_devices)

  # 5. Final cleanup: remove punctuation and extra whitespace
  x <- str_remove_all(x, "[[:punct:]]")
  x <- str_squish(x)

  return(x)
}

# load file generated from rxnorm API
rxnorm_map <- readRDS(paste0(
  data_dir,
  "processed/rxnorm_mapped/meps_ndc_to_brand_full.rds"
)) |>
  setDT()

# Clean NDCs & Rename Columns
rxnorm_map[, ndc := stringr::str_pad(gsub("-", "", ndc), 11, pad = "0")]
rxnorm_map[, `:=`(
  rxnorm_rxname = toupper(trimws(product_name)),
  rxnorm_genind = brand_generic
)]

rxnorm_map <- unique(
  rxnorm_map[, .(ndc, rxnorm_rxname, rxnorm_genind)],
  by = "ndc"
)

# UPDATE DICTIONARY
rx_brands_raw <- unique(rxnorm_map[rxnorm_genind == "Brand", rxnorm_rxname])
rx_dict_brands <- data.table(
  name = unique(clean_drug_text(rx_brands_raw)),
  type = "Brand"
)

rx_generics_raw <- unique(rxnorm_map[rxnorm_genind == "Generic", rxnorm_rxname])
rx_dict_generics <- data.table(
  name = unique(clean_drug_text(rx_generics_raw)),
  type = "Generic"
)

name_dict <- rbind(name_dict, rx_dict_brands, rx_dict_generics)
name_dict <- name_dict[name != ""]
name_dict <- unique(name_dict, by = "name")

# ==============================================================================
# STEP 3: PREPARE SSR (imputed version)
# ==============================================================================

ssr <- fread(paste0(
  data_dir,
  "processed/ssr_cleaned/ssr_gtn_annual_imputed.csv"
))[
  payer == 'Total' & smoothing_type == "4q moving average"
]

ssr[, product_clean := toupper(trimws(product_clean))]

# Keep product_clean for the NDC merge, but create ssr_rxname for the text match
ssr <- ssr[, .(
  product_clean,
  ssr_rxname = product_clean, # Duplicate column for later text matching
  disease_area,
  year_id,
  gtn_final
)]

# ==============================================================================
# STEP 4: PREPARE SSR-NDC crosswalk
# ==============================================================================

ndc_xw_path <- paste0(data_dir, "raw/ssr_health_ndcs/NDCUOM.PrStr.DATA.csv")
ndc_xw <- as.data.table(read.delim(
  ndc_xw_path,
  fileEncoding = "UTF-16LE",
  colClasses = "character"
))
setnames(ndc_xw, new = c("product", "strength", "ndc", "unit"))

# STANDARDIZE THE NDCs
ndc_xw[, ndc := gsub("[^0-9]", "", ndc)]
ndc_xw[, ndc := str_pad(ndc, width = 11, side = "left", pad = "0")]
ndc_xw <- ndc_xw[!is.na(ndc) & nchar(ndc) == 11]

# STANDARDIZE THE PRODUCT NAMES
ndc_xw[, product_clean := toupper(trimws(product))]
ndc_xw <- unique(ndc_xw[, .(product_clean, ndc)])

# MERGE WITH SSR DATA
ssr_ndc_bridge <- merge(
  ssr,
  ndc_xw,
  by = "product_clean",
  all.x = TRUE,
  allow.cartesian = TRUE
)

# ==============================================================================
# STEP 5: LOAD & CLEAN MEPS
# ==============================================================================

# Load in MEPS data
meps <- open_dataset(paste0(data_dir, "raw/meps/USA_MEPS_RX.parquet")) |>
  filter(toc == 'RX' & year_id >= START_YEAR) |>
  select(year_id, ndc, rxname, tot_pay_amt) |>
  as.data.table() |>
  collect() |>
  mutate(year_id = as.integer(year_id))

# Clean missing codes
meps[, ndc := fifelse(nchar(trimws(ndc)) == 11, trimws(ndc), NA_character_)]
meps[, rxname := fifelse(trimws(rxname) == "-15", NA_character_, rxname)]

# Check valid NDC's
meps[, ndc_present := !is.na(ndc) & ndc != ""]
meps[, in_fdb := ndc_present & (ndc %in% fdb$ndc)]

# CHECK 1: What percent of spending has a valid NDC?
spend_valid_ndc <- meps[ndc_present == TRUE, sum(tot_pay_amt, na.rm = TRUE)] /
  sum(meps$tot_pay_amt)

print(paste0(
  "Spend with NDC: ",
  round(spend_valid_ndc * 100, 1),
  "%"
))

# ==============================================================================
# STEP 6: TIER 1 MATCH (DIRECT MEPS NDC TO SSR NDC)
# ==============================================================================
# This is Phase 2 & 3: We do this BEFORE text matching because NDCs are exact.

# 1. PREPARE THE BRIDGE FOR JOINING
# If an NDC accidentally maps to two SSR products in the same year, take the mean
# rebate to prevent blowing up the MEPS dataset with duplicate claims.
bridge_tier1 <- ssr_ndc_bridge[
  !is.na(ndc),
  .(
    gtn_tier1 = mean(gtn_final, na.rm = TRUE),
    ssr_product_tier1 = first(product_clean),
    ssr_disease_area_tier1 = first(disease_area)
  ),
  by = .(ndc, year_id)
]

# 2. PERFORM THE TIER 1 MERGE
# Join MEPS directly to the SSR Bridge using NDC and Year
meps <- merge(
  meps,
  bridge_tier1,
  by = c("ndc", "year_id"),
  all.x = TRUE
)

# Flag successful Tier 1 matches
meps[, matched_ssr_ndc := !is.na(gtn_tier1)]

# 3. TIER 1 VALIDATION
# How much MEPS spending mapped to SSR immediately via NDC?
meps_ndc_spend <- meps[matched_ssr_ndc == TRUE, sum(tot_pay_amt, na.rm = TRUE)]
meps_total_spend <- sum(meps$tot_pay_amt, na.rm = TRUE)

cat("\n=== TIER 1: SSR DIRECT NDC MATCH ===\n")
cat(
  "Direct NDC Spend Matched: $",
  formatC(meps_ndc_spend, format = "f", big.mark = ",", digits = 0),
  "\n"
)
cat(
  "Spend Match Rate (of Total MEPS Spend):",
  round((meps_ndc_spend / meps_total_spend) * 100, 1),
  "%\n\n"
)

# NOTE: This match rate is out of TOTAL spend (including generics).
# ~45.8%

# ==============================================================================
# STEP 7: MERGE MEPS WITH DICTIONARIES
# ==============================================================================

# 1. Merge FDB (Tier 1 - The Gold Standard)
meps_merged <- merge(meps, fdb, by = "ndc", all.x = TRUE)

# 2. Merge RxNorm (Tier 2 - The Rescue for Missing NDCs)
meps_merged <- merge(meps_merged, rxnorm_map, by = "ndc", all.x = TRUE)

# 3. Text Matching Prep (Tier 3 - The Rescue for Bad/Missing NDCs)
#    Apply the cleaning function to the raw MEPS name
meps_merged[, text_match_name := clean_drug_text(rxname)]

# 4. Merge Name Dictionary
#    This adds the 'type' column (Brand/Generic) based on the text name
meps_merged <- merge(
  meps_merged,
  name_dict,
  by.x = "text_match_name",
  by.y = "name",
  all.x = TRUE
)

# CHECK: How much spending is covered by each tier?
spend_total <- meps_merged[, sum(tot_pay_amt, na.rm = T)]
spend_fdb <- meps_merged[!is.na(fdb_genind), sum(tot_pay_amt, na.rm = T)]
spend_rxnorm <- meps_merged[
  is.na(fdb_genind) & !is.na(rxnorm_genind),
  sum(tot_pay_amt, na.rm = T)
]
spend_text <- meps_merged[
  is.na(fdb_genind) & is.na(rxnorm_genind) & !is.na(type),
  sum(tot_pay_amt, na.rm = T)
]

cat("\n=== IDENTIFICATION WATERFALL ===\n")
cat(
  "1. Identified by FDB NDC:   ",
  round(spend_fdb / spend_total * 100, 1),
  "%\n"
)
cat(
  "2. Rescued by RxNorm NDC:   ",
  round(spend_rxnorm / spend_total * 100, 1),
  "%\n"
)
cat(
  "3. Rescued by Text Match:   ",
  round(spend_text / spend_total * 100, 1),
  "%\n"
)


# ==============================================================================
# STEP 8: THE FUZZY RESCUE (Fixing the 19% Gap)
# ==============================================================================

# 1. Isolate the "Unidentified" pile (No NDC match, No Text match)
unidentified <- meps_merged[
  is.na(fdb_genind) & is.na(rxnorm_genind) & is.na(type)
]

# 2. Group by text name to save processing time (Don't match 1 million rows, just the unique names)
#    Focus on high-spend items first
unidentified_grouped <- unidentified[,
  .(spend = sum(tot_pay_amt, na.rm = T)),
  by = text_match_name
][order(-spend)]

# 3. Create a Lookup Vector from your Master Dictionary
#    (We want to match these messy names to the clean FDB/RxNorm names)
dict_names <- name_dict$name

# 4. RUN FUZZY MATCHING
#    Warning: This can be slow if 'unidentified_grouped' is huge.
#    We limit to items with > $0 spend to be efficient.
candidates <- unidentified_grouped[spend > 0]

#    'amatch' finds the index of the closest match in the dictionary
#    maxDist = 2 allows for 1-2 typos (e.g. "LIPTOR" -> "LIPITOR")
cat("Running fuzzy match on", nrow(candidates), "unique text strings...\n")
matches_index <- amatch(
  candidates$text_match_name,
  dict_names,
  method = "jw",
  maxDist = 0.2
)

# 5. Assign the Matched Name
candidates$fuzzy_name <- dict_names[matches_index]

# 6. Merge Fuzzy Results back to Main Data
#    We only keep matches where we found something
valid_fuzzy <- candidates[!is.na(fuzzy_name), .(text_match_name, fuzzy_name)]

#    Get the 'type' (Brand/Generic) from the dictionary for these new fuzzy matches
valid_fuzzy <- merge(
  valid_fuzzy,
  name_dict,
  by.x = "fuzzy_name",
  by.y = "name",
  all.x = TRUE
)

#    Update the main MEPS dataset
meps_merged <- merge(
  meps_merged,
  valid_fuzzy,
  by = "text_match_name",
  all.x = TRUE,
  suffixes = c("", "_fuzzy")
)

#    If we found a fuzzy match, update the 'type' column
meps_merged[!is.na(type_fuzzy), type := type_fuzzy]
meps_merged[!is.na(fuzzy_name), text_match_name := fuzzy_name] # Update the name to the clean version

meps_merged[,
  final_status := fcase(
    # 1. Tier 1: Proven Brand (FDB/RxNorm/Text said "Brand")
    (!is.na(fdb_genind) & fdb_genind == "Brand") | (!is.na(rxnorm_genind) & rxnorm_genind == "Brand") | (!is.na(type) & type == "Brand") | (!is.na(type_fuzzy) & type_fuzzy == "Brand")         , "Brand"   ,

    # 2. Tier 2: Proven Generic (FDB/RxNorm/Text said "Generic")
    (!is.na(fdb_genind) & fdb_genind == "Generic") | (!is.na(rxnorm_genind) & rxnorm_genind == "Generic") | (!is.na(type) & type == "Generic") | (!is.na(type_fuzzy) & type_fuzzy == "Generic") , "Generic" ,

    # 4. Tier 4: The Leftovers
    default = "Unidentified"
  )
]

# Before we assign 0% rebates, we MUST verify this pile is just devices/vague stuff.
unidentified_pile <- meps_merged[
  final_status == "Unidentified",
  .(spend = sum(tot_pay_amt, na.rm = T)),
  by = rxname
][order(-spend)]

print("Top 20 Unidentified Items (Check for Missed Brands):")
print(head(unidentified_pile, 20))

# ==============================================================================
# STEP 9: CREATE FINAL JOIN NAME
# ==============================================================================

# 1. Select the Best Available Name (The Hierarchy)
meps_merged[,
  raw_join_name := fcase(
    !is.na(fdb_brand)              , fdb_brand       ,
    !is.na(rxnorm_rxname)          , rxnorm_rxname   ,
    !is.na(type) & type == "Brand" , text_match_name ,
    default = text_match_name
  )
]

# 2. Clean it one last time using our new rigorous function
meps_merged[, final_join_name := clean_drug_text(raw_join_name)]

# 3. Fill NAs (If cleaning resulted in empty string, use original text)
meps_merged[
  final_join_name == "" | is.na(final_join_name),
  final_join_name := clean_drug_text(rxname)
]
# ==============================================================================
# STEP 12: THE FAKE BRANDS (Demote to Generic -> 0% Rebate)
# ==============================================================================

true_generics <- c(
  "ANTISEPTIC",
  "CARDIOVASCULAR SUPPORT",
  "PRAVASTATIN",
  "BUPROPION XL",
  "FENOFIBRATE",
  "AZITHROMYCIN",
  "METOPROLOL SUCCINATE",
  "DIVALPROEX",
  "TROPICAL STYLE",
  "CITALOPRAM HBR",
  "LATANOPROST",
  "FLOMAX",
  "CLOBETASOL PROPIONATE",
  "CYCLOBENZAPRINE",
  "MORPHINE SULFATE",
  "PROPRANOLOL",
  "TRIHIBIT",
  "FAMOTIDINE",
  "WIXELA INHUB",
  "THERA ANTIFUNGAL",
  "AMOXICILLINCLAVULANATE POTASS",
  "ENALAPRIL MALEATE",
  "CLARAVIS",
  "METHOTREXATE",
  "DIURETIC SOFTGEL",
  "HYDROXYZINE",
  "VITAMIN D2",
  "ANTICOAGULANT CITRATE",
  "VERAPAMIL",
  "TRIAMTERENEHYDROCHLOROTHIAZID",
  "ERYTHROMYCIN",
  "ARIMIDEX",
  "NEURONTIN",
  "CARTIA XT",
  "TRETINOIN",
  "NIACIN",
  "CIPROFLOXACIN",
  "ENDOCET",
  "CARBIDOPALEVODOPA",
  "VERAPAMIL SR",
  "PROPRANOLOLHYDROCHLOROTHIAZID",
  "VALIUM",
  "DORZOLAMIDE",
  "METHYLPREDNISOLONE",
  "SEROQUEL",
  "LEVOXYL",
  "VICODIN",
  "CLARINEX",
  "PROMETHAZINE",
  "DILANTIN",
  "DEPAKOTE",
  "PROTONIX",
  "NORVIR",
  "GLUCOPHAGE",
  "CIMETIDINE",
  "NASALCARE FOR KIDS",
  "HYDROCORTISONE",
  "ULTRAM",
  "TICALAST",
  "ALLEGRA",
  "NP THYROID",
  "LORTAB",
  "TRAVOPROST",
  "AMNESTEEM",
  "ANTIFUNGAL",
  "INDOMETHACIN",
  "ALTACE",
  "PROPECIA",
  "MILLIPRED",
  "VANCOMYCIN",
  "OCELLA",
  "BENZACLIN",
  "INSULIN ASPART",
  "CICLOPIROX",
  "LOSARTAN",
  "ZYRTEC",
  "ZANTAC",
  "SOMA",
  "DIGOX",
  "DILT",
  "ESTROGENS",
  "ARMOUR THYROID",
  "MYORISAN",
  "ZENATANE",
  "ANALGESIC"
)

meps_merged[final_join_name %in% true_generics, final_status := "Generic"]
meps_merged[final_join_name == "NOVOLIN N", final_join_name := "NOVOLIN"]


# ==============================================================================
# STEP 13: THE REBATE WATERFALL (Tiers 1, 2, 3)
# ==============================================================================

# Initialize everyone cleanly
meps_merged[, final_rebate_pct := 0.0]
meps_merged[, match_tier := "Unmatched"]
meps_merged[is.na(matched_ssr_ndc), matched_ssr_ndc := FALSE]

# ------------------------------------------------------------------------------
# TIER 1: EXACT NDC MATCH
# ------------------------------------------------------------------------------
meps_merged[
  matched_ssr_ndc == TRUE,
  `:=`(
    final_rebate_pct = gtn_tier1,
    match_tier = "Tier 1: NDC"
  )
]

# ------------------------------------------------------------------------------
# TIER 2: TEXT MATCH
# ------------------------------------------------------------------------------
meps_merged[, temp_row_id := .I]

brands_missing_ndc <- meps_merged[
  final_status == "Brand" & match_tier == "Unmatched"
]
ssr_lookup <- unique(ssr$ssr_rxname)

brands_missing_ndc[,
  ssr_match_key := fcase(
    final_join_name %in% ssr_lookup          , final_join_name          ,
    word(final_join_name, 1) %in% ssr_lookup , word(final_join_name, 1) ,
    default = NA_character_
  )
]

brands_text_matched <- merge(
  brands_missing_ndc[!is.na(ssr_match_key)],
  ssr[, .(ssr_rxname, year_id, gtn_final)],
  by.x = c("ssr_match_key", "year_id"),
  by.y = c("ssr_rxname", "year_id"),
  all.x = TRUE
)

meps_merged[
  brands_text_matched[!is.na(gtn_final)],
  on = "temp_row_id",
  `:=`(
    final_rebate_pct = i.gtn_final,
    match_tier = "Tier 2: Text"
  )
]
meps_merged[, temp_row_id := NULL]

# ------------------------------------------------------------------------------
# TIER 3: CLASS AVERAGE IMPUTATION (Private Brands & Vague Classes)
# ------------------------------------------------------------------------------

# 1. THE GHOST-BUSTER: Nuke lingering columns from previous runs
cols_to_remove <- c(
  "ssr_area",
  "ssr_area.x",
  "ssr_area.y",
  "avg_class_rebate",
  "avg_class_rebate.x",
  "avg_class_rebate.y"
)
for (col in cols_to_remove) {
  if (col %in% names(meps_merged)) meps_merged[, (col) := NULL]
}

# 2. Build the mappings
class_text_map <- rbind(
  data.table(name = "IMMUNOLOGIC AGENTS", ssr_area = "DMARDs (Anti-TNF)"),
  data.table(
    name = "SELECTIVE IMMUNOSUPPRESSANTS",
    ssr_area = "Immunosuppressants"
  ),
  data.table(name = "ANTINEOPLASTICS", ssr_area = "Oncology (MTI)"),
  data.table(name = "20 ANTINEOPLASTICS", ssr_area = "Oncology (MTI)"),
  data.table(
    name = "MISCELLANEOUS ANTINEOPLASTICS",
    ssr_area = "Oncology (MTI)"
  ),
  data.table(name = "ANTIVIRAL AGENTS", ssr_area = "Antiviral (Other)"),
  data.table(name = "PSYCHOTHERAPEUTIC AGENTS", ssr_area = "Atyp. antipsych"),
  data.table(
    name = "CENTRAL NERVOUS SYSTEM AGENTS",
    ssr_area = "Anticonvulsants"
  )
)

private_brand_map <- rbind(
  # Respiratory
  data.table(name = "SPIRIVA", ssr_area = "COPD (Bronchodilators)"),
  data.table(name = "ATROVENT", ssr_area = "COPD (Bronchodilators)"),
  data.table(name = "PERFOROMIST", ssr_area = "COPD (Bronchodilators)"),
  data.table(name = "BROVANA", ssr_area = "COPD (Bronchodilators)"),
  data.table(name = "COMBIVENT", ssr_area = "COPD (Combos)"),
  data.table(name = "STIOLTO", ssr_area = "COPD (Combos)"),
  data.table(name = "QVAR", ssr_area = "COPD (Glucocorticoids)"),
  data.table(name = "ASMANEX", ssr_area = "COPD (Glucocorticoids)"),
  data.table(name = "PROVENTIL", ssr_area = "COPD/Asthma (Other)"),
  data.table(name = "PROAIR", ssr_area = "COPD/Asthma (Other)"),
  data.table(name = "XOPENEX", ssr_area = "COPD/Asthma (Other)"),
  data.table(name = "AUVIQ", ssr_area = "COPD/Asthma (Other)"),

  # Pain / Opioids
  data.table(name = "OXYCONTIN", ssr_area = "Pain (Opioids)"),
  data.table(name = "HYSINGLA", ssr_area = "Pain (Opioids)"),
  data.table(name = "OPANA", ssr_area = "Pain (Opioids)"),
  data.table(name = "KADIAN", ssr_area = "Pain (Opioids)"),
  data.table(name = "AVINZA", ssr_area = "Pain (Opioids)"),
  data.table(name = "BUTRANS", ssr_area = "Pain (Opioids)"),
  data.table(name = "NORCO", ssr_area = "Pain (Opioids)"),

  # Pain / Non-Opioids
  data.table(name = "LIDODERM", ssr_area = "Pain, Inflammation (Non-Opioids)"),
  data.table(name = "FLECTOR", ssr_area = "Pain, Inflammation (Non-Opioids)"),
  data.table(name = "SKELAXIN", ssr_area = "Pain, Inflammation (Non-Opioids)"),
  data.table(name = "AMRIX", ssr_area = "Pain, Inflammation (Non-Opioids)"),
  data.table(name = "SOLARAZE", ssr_area = "Pain, Inflammation (Non-Opioids)"),

  # CNS - Anticonvulsants
  data.table(name = "GRALISE", ssr_area = "Anticonvulsants"),
  data.table(name = "HORIZANT", ssr_area = "Anticonvulsants"),
  data.table(name = "TRILEPTAL", ssr_area = "Anticonvulsants"),
  data.table(name = "CARBATROL", ssr_area = "Anticonvulsants"),

  # CNS - ADHD
  data.table(name = "VYVANSE", ssr_area = "ADHD"),
  data.table(name = "FOCALIN", ssr_area = "ADHD"),
  data.table(name = "ADDERALL", ssr_area = "ADHD"),
  data.table(name = "CONCERTA", ssr_area = "ADHD"),
  data.table(name = "DAYTRANA", ssr_area = "ADHD"),
  data.table(name = "QUILLIVANT", ssr_area = "ADHD"),
  data.table(name = "QUILLICHEW", ssr_area = "ADHD"),
  data.table(name = "INTUNIV", ssr_area = "ADHD"),
  data.table(name = "JORNAY PM", ssr_area = "ADHD"),

  # CNS - Depression
  data.table(name = "PROZAC", ssr_area = "Depression"),
  data.table(name = "LEXAPRO", ssr_area = "Depression"),

  # CNS - Other
  data.table(name = "PROVIGIL", ssr_area = "Narcolepsy"),
  data.table(name = "NUVIGIL", ssr_area = "Narcolepsy"),
  data.table(name = "MIRAPEX", ssr_area = "Parkinson's"),
  data.table(name = "AZILECT", ssr_area = "Parkinson's"),
  data.table(name = "IMITREX", ssr_area = "Migraine (Triptans)"),
  data.table(name = "VALIUM", ssr_area = "Benzodiazepines"),
  data.table(name = "KLONOPIN", ssr_area = "Benzodiazepines"),
  data.table(name = "SILENOR", ssr_area = "Insomnia"),

  # Cardiovascular
  data.table(name = "PRADAXA", ssr_area = "Oral Anticoagulants"),
  data.table(name = "COUMADIN", ssr_area = "Oral Anticoagulants"),
  data.table(name = "BENICAR", ssr_area = "Cardiovascular (ARBs, ARB combos)"),
  data.table(name = "MICARDIS", ssr_area = "Cardiovascular (ARBs, ARB combos)"),
  data.table(
    name = "MICARDIS HCT",
    ssr_area = "Cardiovascular (ARBs, ARB combos)"
  ),
  data.table(name = "BYSTOLIC", ssr_area = "Cardiovascular (Beta Blockers)"),
  data.table(name = "INDERAL LA", ssr_area = "Cardiovascular (Beta Blockers)"),
  data.table(
    name = "LOTREL",
    ssr_area = "Cardiovascular (Calcium Channel Blockers)"
  ),
  data.table(name = "LIVALO", ssr_area = "Cholesterol (Non-PCSK9s)"),
  data.table(name = "SIMCOR", ssr_area = "Cholesterol (Non-PCSK9s)"),
  data.table(name = "ADVICOR", ssr_area = "Cholesterol (Non-PCSK9s)"),
  data.table(name = "ANTARA", ssr_area = "Cholesterol (Non-PCSK9s)"),
  data.table(name = "NITROSTAT", ssr_area = "Cardiovascular (Other)"),
  data.table(name = "VASCULERA", ssr_area = "Cardiovascular (Other)"),

  # Diabetes / Metabolic
  data.table(
    name = "NOVOLOG",
    ssr_area = "Diabetes (Rapid-acting/mix insulins)"
  ),
  data.table(
    name = "NOVOLOG MIX 7030",
    ssr_area = "Diabetes (Rapid-acting/mix insulins)"
  ),
  data.table(
    name = "NOVOLIN 7030",
    ssr_area = "Diabetes (Rapid-acting/mix insulins)"
  ),
  data.table(name = "FIASP", ssr_area = "Diabetes (Rapid-acting/mix insulins)"),
  data.table(name = "NOVOLIN", ssr_area = "Diabetes (Short-acting insulins)"),
  data.table(name = "LEVEMIR", ssr_area = "Diabetes (Long-acting insulins)"),
  data.table(name = "SEMGLEE", ssr_area = "Diabetes (Long-acting insulins)"),
  data.table(name = "VICTOZA 3PAK", ssr_area = "Diabetes (GLP-1 Agonists)"),
  data.table(name = "ACTOS", ssr_area = "Diabetes (Other)"),
  data.table(name = "AVANDIA", ssr_area = "Diabetes (Other)"),
  data.table(name = "ACTOPLUS MET", ssr_area = "Diabetes (Other)"),
  data.table(name = "AVANDARYL", ssr_area = "Diabetes (Other)"),

  # Renal / Bone Metabolism
  data.table(name = "RENAGEL", ssr_area = "Hyperphosphatemia"),
  data.table(name = "RENVELA", ssr_area = "Hyperphosphatemia"),
  data.table(name = "FOSRENOL", ssr_area = "Hyperphosphatemia"),
  data.table(name = "ZEMPLAR", ssr_area = "Hyperparathyroidism"),
  data.table(name = "RAYALDEE", ssr_area = "Hyperparathyroidism"),
  data.table(name = "HECTOROL", ssr_area = "Hyperparathyroidism"),

  # Female Hormone Therapy / Contraceptives
  data.table(name = "VAGIFEM", ssr_area = "Female Hormone Therapy"),
  data.table(name = "ESTRING", ssr_area = "Female Hormone Therapy"),
  data.table(name = "VIVELLEDOT", ssr_area = "Female Hormone Therapy"),
  data.table(name = "CLIMARA", ssr_area = "Female Hormone Therapy"),
  data.table(name = "ESTROGEL", ssr_area = "Female Hormone Therapy"),
  data.table(name = "YUVAFEM", ssr_area = "Female Hormone Therapy"),
  data.table(name = "XULANE", ssr_area = "Contraceptives"),
  data.table(name = "TAYTULLA", ssr_area = "Contraceptives"),
  data.table(name = "JUNEL FE", ssr_area = "Contraceptives"),
  data.table(name = "ELURYNG", ssr_area = "Contraceptives"),

  # Thyroid
  data.table(name = "TIROSINT", ssr_area = "Hypothyroidism"),

  # GI
  data.table(name = "PENTASA", ssr_area = "IBD"),
  data.table(name = "LIALDA", ssr_area = "IBD"),
  data.table(name = "PREVACID", ssr_area = "Anti-Ulcer (PPIs)"),
  data.table(name = "DEXILANT", ssr_area = "Anti-Ulcer (PPIs)"),
  data.table(name = "AMITIZA", ssr_area = "Constipation"),
  data.table(name = "ZENPEP", ssr_area = "Exocrine Pancreatic Insufficiency"),
  data.table(name = "DICLEGIS", ssr_area = "Antiemetics"),

  # Ophthalmology
  data.table(name = "RESTASIS", ssr_area = "Ophthalmics (Other)"),
  data.table(name = "VIGAMOX", ssr_area = "Ophthalmics (Other)"),
  data.table(name = "DUREZOL", ssr_area = "Ophthalmics (Other)"),
  data.table(name = "CEQUA", ssr_area = "Ophthalmics (Other)"),
  data.table(name = "TOBRADEX", ssr_area = "Ophthalmics (Other)"),
  data.table(name = "PROLENSA", ssr_area = "Ophthalmics (Other)"),
  data.table(name = "LUMIGAN", ssr_area = "Glaucoma"),
  data.table(name = "ALPHAGAN P", ssr_area = "Glaucoma"),
  data.table(name = "COMBIGAN", ssr_area = "Glaucoma"),
  data.table(name = "AZOPT", ssr_area = "Glaucoma"),
  data.table(name = "TRAVATAN", ssr_area = "Glaucoma"),
  data.table(name = "BETIMOL", ssr_area = "Glaucoma"),

  # Dermatology
  data.table(name = "ORACEA", ssr_area = "Rosacea"),
  data.table(name = "ABSORICA", ssr_area = "Acne (Orals)"),
  data.table(name = "AMNESTEEM", ssr_area = "Acne (Orals)"),
  data.table(name = "EPIDUO", ssr_area = "Acne (Topicals)"),
  data.table(name = "DIFFERIN", ssr_area = "Acne (Topicals)"),
  data.table(name = "TAZORAC", ssr_area = "Acne (Topicals)"),
  data.table(name = "TACLONEX", ssr_area = "Psoriasis"),
  data.table(name = "CLOBEX", ssr_area = "Atopic dermatitis"),

  # Nasal
  data.table(name = "DYMISTA", ssr_area = "Nasal glucocorticoids"),
  data.table(name = "ASTEPRO", ssr_area = "Nasal glucocorticoids"),

  # Oncology
  data.table(name = "FEMARA", ssr_area = "Oncology (Hormonal / GnRH)"),

  # Infectious Disease
  data.table(name = "LEVAQUIN", ssr_area = "Antibacterials"),
  data.table(name = "CIPRODEX", ssr_area = "Antibacterials"),
  data.table(name = "CIPRO", ssr_area = "Antibacterials"),
  data.table(name = "DORYX", ssr_area = "Antibacterials"),
  data.table(name = "AUGMENTIN", ssr_area = "Antibacterials"),
  data.table(name = "AVELOX", ssr_area = "Antibacterials"),
  data.table(name = "VALCYTE", ssr_area = "Antiviral (Other)"),

  # Autoimmune
  data.table(name = "OTREXUP", ssr_area = "DMARDs (Other)"),
  data.table(name = "RASUVO", ssr_area = "DMARDs (Other)")
)

full_tier3_map <- rbind(private_brand_map, class_text_map)
full_tier3_map <- unique(full_tier3_map, by = "name")

# Name fixes (cleaning mangled these)
meps_merged[final_join_name == "NOVOLIN N", final_join_name := "NOVOLIN"]
meps_merged[
  final_join_name == "LOSARTAN POTASSIUM",
  final_join_name := "LOSARTAN"
]
meps_merged[final_join_name == "ESTROVEN", final_join_name := "ESTROGENS"] # fix fuzzy mismatch before demoting


# 3. FAST DATA.TABLE UPDATE JOIN: Attach the SSR Disease Area
meps_merged[
  full_tier3_map,
  on = c("final_join_name" = "name"),
  ssr_area := i.ssr_area
]


# 4. Calculate SSR Class Averages
ssr_class_avgs <- ssr[,
  .(avg_class_rebate = mean(gtn_final, na.rm = TRUE)),
  by = .(disease_area, year_id)
]

# 5. FAST DATA.TABLE UPDATE JOIN: Attach the Rebate
meps_merged[
  ssr_class_avgs,
  on = c("ssr_area" = "disease_area", "year_id"),
  avg_class_rebate := i.avg_class_rebate
]

# 6. Apply Tier 3 ONLY to rows that are still Unmatched
meps_merged[
  match_tier == "Unmatched" & !is.na(avg_class_rebate),
  `:=`(
    final_rebate_pct = avg_class_rebate,
    match_tier = "Tier 3: Class Imputed",
    final_status = "Class_Imputed" # We change status so the denominator captures these!
  )
]

# Clean up temporary columns
meps_merged[, c("ssr_area", "avg_class_rebate") := NULL]


# ==============================================================================
# STEP 14: CALCULATE FINAL NET SPENDING
# ==============================================================================

meps_merged[final_rebate_pct < 0, final_rebate_pct := 0]
meps_merged[final_rebate_pct > 0.95, final_rebate_pct := 0.95]
meps_merged[is.na(final_rebate_pct), final_rebate_pct := 0]

meps_merged[, net_pay_amt := tot_pay_amt * (1 - final_rebate_pct)]


# ==============================================================================
# STEP 15: FINAL WATERFALL DIAGNOSTICS
# ==============================================================================

cat("\n========================================\n")
cat(" FINAL REBATE WATERFALL MATCHING \n")
cat("========================================\n")

total_brand_universe <- meps_merged[
  final_status %in% c("Brand", "Class_Imputed"),
  sum(tot_pay_amt, na.rm = T)
]

ndc_spend <- meps_merged[
  match_tier == "Tier 1: NDC",
  sum(tot_pay_amt, na.rm = T)
]
text_spend <- meps_merged[
  match_tier == "Tier 2: Text",
  sum(tot_pay_amt, na.rm = T)
]
class_spend <- meps_merged[
  match_tier == "Tier 3: Class Imputed",
  sum(tot_pay_amt, na.rm = T)
]

cat(sprintf(
  "1. Exact NDC Match:     $%13s (%4.1f%% of Brand Spend)\n",
  formatC(ndc_spend, format = "d", big.mark = ","),
  ndc_spend / total_brand_universe * 100
))
cat(sprintf(
  "2. Exact Text Match:    $%13s (%4.1f%% of Brand Spend)\n",
  formatC(text_spend, format = "d", big.mark = ","),
  text_spend / total_brand_universe * 100
))
cat(sprintf(
  "3. Class Average Match: $%13s (%4.1f%% of Brand Spend)\n",
  formatC(class_spend, format = "d", big.mark = ","),
  class_spend / total_brand_universe * 100
))
cat("----------------------------------------\n")
cat(sprintf(
  "TOTAL BRAND MATCH RATE:               %4.1f%%\n",
  (ndc_spend + text_spend + class_spend) / total_brand_universe * 100
))
cat("========================================\n")


# ==============================================================================
# STEP 16: DIAGNOSTIC (THE NEW "MISSING")
# ==============================================================================

unmatched_brands <- meps_merged[
  final_status == "Brand" & match_tier == "Unmatched",
  .(spend = sum(tot_pay_amt, na.rm = TRUE)),
  by = .(rxname, final_join_name)
][order(-spend)]

cat("\n========================================\n")
cat(" TOP 15 UNMATCHED BRANDS (The Remaining Gap) \n")
cat("========================================\n")
print(head(unmatched_brands, 50))

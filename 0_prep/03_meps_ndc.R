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
# STEP 2: PREPARE RXNORM (Rescue Lookup)
# ==============================================================================

# cleaning function (to apply to both brand & generics)
clean_drug_text <- function(x) {
  x <- toupper(x)
  # Remove packaging/noise
  x <- str_remove_all(x, "\\(.*?\\)") # anything in parentheses
  x <- str_remove_all(x, "\\b\\d+\\s*DAY\\b")
  x <- str_remove_all(x, "\\bREFORMULATED.*")
  x <- str_remove_all(x, "\\b(FIRST MONTH|STARTER|TITRATION).*")
  # Remove dosages
  x <- str_remove_all(x, "\\b\\d+\\.?\\d*\\/\\d+\\.?\\d*") # dosage ratios
  x <- str_remove_all(x, "\\b\\d+\\.?\\d*\\s*(MG|MCG|ML|G|L|MEQ|UNIT[S]?)\\b")
  # Remove forms/salts
  x <- str_remove_all(
    x,
    "\\b(TAB|CAP|INJ|SOL|HCL|SODIUM|EQ|ER|XR|DR|SYR|SUSP)\\b"
  )
  # Cleanup
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

#  Clean NDCs & Rename Columns
rxnorm_map[, ndc := stringr::str_pad(gsub("-", "", ndc), 11, pad = "0")]
rxnorm_map[, `:=`(
  rxnorm_rxname = toupper(trimws(product_name)),
  rxnorm_genind = brand_generic
)]

# Deduplicate for the NDC Merge (One row per NDC)
rxnorm_map <- unique(
  rxnorm_map[, .(ndc, rxnorm_rxname, rxnorm_genind)],
  by = "ndc"
)

# UPDATE DICTIONARY

# brands
rx_brands_raw <- unique(rxnorm_map[rxnorm_genind == "Brand", rxnorm_rxname])
rx_dict_brands <- data.table(
  name = unique(clean_drug_text(rx_brands_raw)),
  type = "Brand"
)

# generics
rx_generics_raw <- unique(rxnorm_map[rxnorm_genind == "Generic", rxnorm_rxname])
rx_dict_generics <- data.table(
  name = unique(clean_drug_text(rx_generics_raw)),
  type = "Generic"
)

# Add to Master Dictionary
name_dict <- rbind(name_dict, rx_dict_brands, rx_dict_generics)
name_dict <- name_dict[name != ""]
name_dict <- unique(name_dict, by = "name")

# ==============================================================================
# STEP 3: PREPARE SSR (imputed version)
# ==============================================================================

# Update your Step 3 to keep the rebate column:
ssr <- fread(paste0(
  data_dir,
  "processed/ssr_cleaned/ssr_gtn_annual_imputed.csv"
))[
  payer == 'Total' & smoothing_type == "4q moving average",
][, ssr_rxname := toupper(trimws(product_clean))][, .(
  ssr_rxname,
  disease_area,
  year_id,
  gtn_final
)]

# ==============================================================================
# STEP 4: PREPARE SSR-NDC crosswalk
# ==============================================================================

# 1. READ THE NDC CROSSWALK ----------------------------------------------------

ndc_xw_path <- paste0(data_dir, "raw/ssr_health_ndcs/NDCUOM.PrStr.DATA.csv")
ndc_xw <- as.data.table(
  read.delim(ndc_xw_path, fileEncoding = "UTF-16LE", colClasses = "character")
)
setnames(ndc_xw, new = c("product", "strength", "ndc", "unit"))


# 2. STANDARDIZE THE NDCs ------------------------------------------------------
# A standard NDC must be exactly 11 numeric digits with no hyphens.
# - gsub("[^0-9]", "", ndc) drops hyphens, spaces, and letters.
# - str_pad left-pads with "0" to ensure it reaches 11 characters.

ndc_xw[, ndc := gsub("[^0-9]", "", ndc)]
ndc_xw[, ndc := str_pad(ndc, width = 11, side = "left", pad = "0")]

# Drop empty or blatantly invalid NDCs (must be exactly 11 characters now)
ndc_xw <- ndc_xw[!is.na(ndc) & nchar(ndc) == 11]


# 3. STANDARDIZE THE PRODUCT NAMES ---------------------------------------------
# Apply the exact same cleaning logic used on the SSR file so the strings match perfectly.
ndc_xw[, product_clean := toupper(trimws(product))]
ndc_xw <- unique(ndc_xw[, .(product_clean, ndc)])


# 4. MERGE WITH SSR DATA -------------------------------------------------------
# Join the crosswalk to the imputed SSR annual data.
# Note: allow.cartesian = TRUE is required because 1 Product = Many Years (in SSR)
# AND 1 Product = Many NDCs (in crosswalk). It will multiply the rows correctly.

ssr_ndc_bridge <- merge(
  ssr_annual_imputed,
  ndc_xw,
  by = "product_clean",
  all.x = TRUE, # Keep SSR products even if they don't have an NDC map
  allow.cartesian = TRUE
)


# 5. VALIDATION ----------------------------------------------------------------
cat("Rows in SSR Annual Data:", nrow(ssr_annual_imputed), "\n")
cat("Rows in the new NDC Bridge:", nrow(ssr_ndc_bridge), "\n\n")

# How successful was the crosswalk?
mapped_prods <- uniqueN(ssr_ndc_bridge[!is.na(ndc)]$product_clean)
total_prods <- uniqueN(ssr_ndc_bridge$product_clean)

cat("SSR Products successfully mapped to >= 1 NDC:", mapped_prods, "\n")
cat("SSR Products missing NDC mapping:", total_prods - mapped_prods, "\n")
cat(
  "Mapping Success Rate (by Product):",
  round((mapped_prods / total_prods) * 100, 1),
  "%\n"
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
# ~45.9%

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
# STEP 8.5: THE MISSING 20% RESCUE (Fixing FDB Errors & Private Brands)
# ==============================================================================

# 1. THE FAKE BRANDS (Demote to Generic -> 0% Rebate)
# FDB called these "Brands", but they are just generic chemicals or vague classes.
true_generics <- c(
  "ANTISEPTIC",
  "CARDIOVASCULAR SUPPORT",
  "PRAVASTATIN",
  "BUPROPION XL",
  "FENOFIBRATE",
  "AZITHROMYCIN"
)
meps_merged[final_join_name %in% true_generics, final_status := "Generic"]

# 2. THE PRIVATE BRANDS (Promote to Tier 3: Class Imputed)
# We map these to SSR Disease Areas to get the class average rebate.
private_brand_map <- rbind(
  data.table(name = "SPIRIVA HANDIHALER", ssr_area = "Respiratory"),
  data.table(name = "SPIRIVA RESPIMAT", ssr_area = "Respiratory"),
  data.table(name = "COMBIVENT RESPIMAT", ssr_area = "Respiratory"),
  data.table(name = "OXYCONTIN", ssr_area = "Central Nervous System"),
  data.table(name = "PRADAXA", ssr_area = "Cardiovascular"),
  data.table(name = "ACTOS", ssr_area = "Metabolic"),
  data.table(name = "BENICAR", ssr_area = "Cardiovascular"),
  data.table(name = "CONCERTA", ssr_area = "Central Nervous System"),
  data.table(name = "GRALISE", ssr_area = "Central Nervous System")
)

# Apply the private mapping
meps_merged <- merge(
  meps_merged,
  private_brand_map,
  by.x = "final_join_name",
  by.y = "name",
  all.x = TRUE
)

# Calculate the Class Average Rebate directly from SSR
ssr_class_avgs <- ssr[,
  .(
    avg_class_rebate = mean(gtn_final, na.rm = TRUE)
  ),
  by = .(disease_area, year_id)
]

# Merge the class averages onto our private brands
meps_merged <- merge(
  meps_merged,
  ssr_class_avgs,
  by.x = c("ssr_area", "year_id"),
  by.y = c("disease_area", "year_id"),
  all.x = TRUE
)

# If a private brand matched, officially change its status to Class_Imputed
meps_merged[
  !is.na(ssr_area),
  `:=`(
    final_status = "Class_Imputed",
    imputed_rebate = avg_class_rebate
  )
]

# Clean up temp columns
meps_merged[, c("ssr_area", "avg_class_rebate") := NULL]

# ==============================================================================
# STEP 9: FLAG THERAPEUTIC CLASSES (NON-MATCHABLE)
# ==============================================================================

# Define pattern of words that indicate a Class, not a Drug
# These are exactly what is showing in your screenshot
class_string <- "AGENTS|ANTINEOPLASTIC|ANTI-INFECTIVE|ANTIVIRAL|INSULIN|ATYPICAL ANTIPSYCHOTICS|
ANTIPSYCHOTIC|ANTICONVULSANT|MODIFIERS|COMBINATIONS|HORMONES|IMMUNOSUPPRESSANT|
ANALGESIC|INHIBITORS|ANTAGONIST|BLOCKER|DIURETIC|RELAXANT|STIMULANT|BENZODIAZEPINE|
ANTIBIOTIC|ANTICOAGUL|ANTIHYPERT|ANTIDIABET|ANTIHYPERLIPID|ANTIARRHYTHM|
ANTIPARKINSON|ANTIMIGRAINE|BRONCHODILATOR|DERMATOLOG|OPHTHALM|METABOLIC|
RESPIRATORY|URINARY|BOWEL|DOPAMINERGIC|NARCOTIC|CARDIOVASCULAR|PSYCHOTHERAP|
COAGULATION|MULTIKINASE|SKELETAL|SEX HORMONE|ANTIFUNGAL|ANTIDEPRESSANT|ANXIOLYTIC|SEDATIVE|HYPNOT|BIOLOGICAL|
BILE ACID|CORTICOSTEROID|ANTIEMETIC|MONOCLONAL|PHENOTHIAZINE|TRICYCLIC|
SSNRI|GLP-1|AZOLE|COAGULATION MEDICINE|5-AMINOSALICYLATES|BILE ACID SEQUESTRANTS|
TOPICAL STEROIDS|AMINOGLYCOSIDES|PURINE NUCLEOSIDES|GLUCOCORTICOIDS|NASAL STEROIDS|
TOPICAL ANESTHETICS|ANTICHOLINERGICS|ANTISPASMODICS|BISPHOSPHONATES|ESTROGENS|ANTIRHEUMATICS|
VASOPRESSORS|MISCELLANEOUS ANTIBIOTICS|THIRD GENERATION CEPHALOSPORINS|CHOLINERGIC AGONISTS|
TOPICAL ANTIBIOTICS|NUTRITIONAL PRODUCTS|NASAL ANTIHISTAMINES AND DECONGESTANTS|
RESPIRATORY INHALANT PRODUCTS|NASAL PREPARATIONS|TOPICAL EMOLLIENTS|TOPICAL ANTIPSORIATICS|
INCRETIN MIMETICS|VASODILATORS|TETRACYCLINES|VIRAL VACCINES|GLUCOCORTICOIDS|NASAL STEROIDS|
TOPICAL ANESTHETICS|ANTICHOLINERGICS|ANTISPASMODICS|BISPHOSPHONATES|ANTIRHEUMATICS|
ANTIHISTAMINES|ANTIMALARIAL|ANTIMETABOLITE|ANTISPASMODIC|ANTIDOTES|CEPHALOSPORINS|
FIBRIC ACID DERIVATIVES|MACROLIDES|MOUTH AND THROAT PRODUCTS|NUTRACEUTICAL PRODUCTS|
NUTRITIONAL SUPPLEMENT|OTIC PREPARATIONS|POLYENES|SULFONYLUR|VAGINAL PREPARATIONS|
ALTERNATIVE MEDICINES|DRUGS USED IN ALCOHOL DEPENDENCE|TOPICAL RUBEFACIENT"

# Remove newlines and extra spaces
class_pattern <- gsub("\n", "", class_string)
class_pattern <- gsub("\\s+", " ", class_pattern) # Squash spaces

device_string <- "MINIMED|ONETOUCH|SOFTCLIX|OT ULTRA|SET MMT|SYRINGE|NEEDLE|PRECISION Q-I-D|GLUCOMETER|
   LANCET|BLOOD GLUC|MONITOR|ACCU-CH|ASCENSIA|FREESTYLE|LIFESCAN|MEDISENSE|SURE STEP|SURESTEP|
   TRUTRACK|TRUETRACK|ACCU-CHECK|GLUCOSE STRIP|TEST STRIP|GLUCOSE KIT|GLUCOSE MONITOR|
   AEROCHAMBER|NEBULIZ|COMPRESSOR|SPACER|INHALER DEVICE|OPTICHAMBER|PEAK FLOW|
   INSULIN SYR|INS SYR|INSULIN PEN|PEN NEEDLE|BD |B-D |MONOJECT|NOVOFI|FLEXPEN|
   WHEELCHAIR|WALKER|CRUTCH|BRACE|SUPPORT|INSERT|INSOLE|SHOES|STOCKING|TED HOSE|
   BANDAGE|GAUZE|DRESSING|TAPE|PAD |WOUND|CATHETER|OSTOMY|COLOSTOMY|ILEOSTOMY|
   BREAST PUMP|BLOOD PRESSURE|SPHYGMO|THERMOMETER|SHARPS|DISPOSAL|CONTAINER|
   PILL BOX|PILL CUTTER|PILL SPLITTER|TABLET CUTTER|TABLET SPLITTER|
   BATTERY|TUBING|MASK|PUMP|RESERVOIR|INFUSION SET|ADMIN SET|
   DIABETIC SHOE|DIABETIC INSERT|ORTHO INSERT|CUSTOM INSERT|DIABETIC FOOT"

device_pattern <- gsub("\n", "", device_string)
device_pattern <- gsub("\\s+", " ", device_pattern)

# 2. Apply Flags
meps_merged[,
  is_drug_class := grepl(class_pattern, text_match_name, ignore.case = TRUE)
]
meps_merged[
  is.na(text_match_name) | text_match_name == "",
  is_drug_class := TRUE
] # Treat missing text as Class/Vague

meps_merged[,
  is_device := grepl(device_pattern, text_match_name, ignore.case = TRUE)
]


# ==============================================================================
# STEP 10: UPDATE FINAL STATUS (PRIORITIZE FLAGS)
# ==============================================================================
# We modify 'final_status' so these categories appear in the final breakdown.
# CRITICAL PRIORITY ORDER:
# 1. Class_Imputed (Already set by high-cost logic - KEEP THIS)
# 2. Devices (If Unidentified/Generic -> Device)
# 3. Vague Classes (If Unidentified/Generic -> Therapeutic_Class)

# Update "Unidentified" or "Generic" items to "Medical_Device" if flagged
meps_merged[
  final_status %in% c("Unidentified", "Generic") & is_device == TRUE,
  final_status := "Medical_Device"
]

# Update "Unidentified" or "Generic" items to "Therapeutic_Class" if flagged
# (Note: This only runs if it wasn't already caught by High-Cost Class Imputation)
meps_merged[
  final_status %in% c("Unidentified", "Generic") & is_drug_class == TRUE,
  final_status := "Therapeutic_Class"
]

# ==============================================================================
# STEP 11: Join name
# ==============================================================================

# Select the Best Available Name (The Hierarchy)
#    We prefer FDB names because they are usually cleaner than MEPS text.
meps_merged[,
  raw_join_name := fcase(
    !is.na(fdb_brand)              , fdb_brand       , # Priority 1: FDB Brand Name
    !is.na(rxnorm_rxname)          , rxnorm_rxname   , # Priority 2: RxNorm Name
    !is.na(type) & type == "Brand" , text_match_name , # Priority 3: Dictionary Match
    default = text_match_name # Fallback: Cleaned MEPS text
  )
]

# 3. Create the Final Join Name (Clean it one last time)
meps_merged[, final_join_name := clean_drug_text(raw_join_name)]

# 4. Fill NAs (If cleaning resulted in empty string, use original text)
meps_merged[
  final_join_name == "" | is.na(final_join_name),
  final_join_name := clean_drug_text(rxname)
]

# ==============================================================================
# STEP 12 & 13: TEXT MATCHING & FINAL REBATE ASSIGNMENT
# ==============================================================================

# We now execute the Waterfall to assign final_rebate_pct.

# 1. Initialize everyone at 0%
meps_merged[, final_rebate_pct := 0.0]

# ------------------------------------------------------------------------------
# TIER 1: EXACT NDC MATCH (From Step 4.5)
# ------------------------------------------------------------------------------
meps_merged[matched_ssr_ndc == TRUE, final_rebate_pct := gtn_tier1]

# ------------------------------------------------------------------------------
# TIER 2: TEXT MATCH (For Brands that failed the NDC match)
# ------------------------------------------------------------------------------
# Identify Brands that missed the NDC match
brands_missing_ndc <- meps_merged[
  final_status == "Brand" & matched_ssr_ndc == FALSE
]

# Create the match key (Exact name, or first word fallback)
ssr_lookup <- unique(ssr$ssr_rxname)
brands_missing_ndc[,
  ssr_match_key := fcase(
    final_join_name %in% ssr_lookup          , final_join_name          ,
    word(final_join_name, 1) %in% ssr_lookup , word(final_join_name, 1) ,
    default = NA_character_
  )
]

# Merge the text-matched Brands with SSR to get the rebate (gtn_final)
brands_text_matched <- merge(
  brands_missing_ndc[!is.na(ssr_match_key)],
  ssr[, .(ssr_rxname, year_id, gtn_final)], # Bring in the rebate!
  by.x = c("ssr_match_key", "year_id"),
  by.y = c("ssr_rxname", "year_id"),
  all.x = TRUE
)

# Bring the Text Match rebates back into the main MEPS dataset
# (We do a fast update join in data.table)
meps_merged[
  brands_text_matched[!is.na(gtn_final)],
  on = .(year_id, ndc, text_match_name), # Merge on unique row identifiers
  `:=`(
    final_rebate_pct = i.gtn_final,
    matched_ssr_text = TRUE
  )
]

# ------------------------------------------------------------------------------
# TIER 3: CLASS AVERAGE IMPUTATION (For Private/Missing Brands)
# ------------------------------------------------------------------------------
# If it is a Class_Imputed drug (from Step 6), apply the class average.
meps_merged[
  final_status == "Class_Imputed" & final_rebate_pct == 0.0,
  final_rebate_pct := imputed_rebate
]

# ------------------------------------------------------------------------------
# TIER 4: GENERICS / DEVICES / UNIDENTIFIED
# ------------------------------------------------------------------------------
# Already defaulted to 0.0!

# ==============================================================================
# CALCULATE FINAL NET SPENDING
# ==============================================================================

# Protect against weird data (cap rebates between 0 and 95%)
meps_merged[final_rebate_pct < 0, final_rebate_pct := 0]
meps_merged[final_rebate_pct > 0.95, final_rebate_pct := 0.95]
meps_merged[is.na(final_rebate_pct), final_rebate_pct := 0]

# The Final Calculation
meps_merged[, net_pay_amt := tot_pay_amt * (1 - final_rebate_pct)]

# ==============================================================================
# FINAL WATERFALL DIAGNOSTICS
# ==============================================================================

# ==============================================================================
# FINAL WATERFALL DIAGNOSTICS
# ==============================================================================

cat("\n========================================\n")
cat(" FINAL REBATE WATERFALL MATCHING \n")
cat("========================================\n")

# The true branded universe is anything that is a Brand OR was forced to a Class Imputed Brand
total_brand_universe <- meps_merged[
  final_status %in% c("Brand", "Class_Imputed"),
  sum(tot_pay_amt, na.rm = T)
]

ndc_spend <- meps_merged[matched_ssr_ndc == TRUE, sum(tot_pay_amt, na.rm = T)]
text_spend <- meps_merged[matched_ssr_text == TRUE, sum(tot_pay_amt, na.rm = T)]
class_spend <- meps_merged[
  final_status == "Class_Imputed",
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
cat("========================================\n")

# ==============================================================================
# DIAGNOSTIC: THE "MISSING 20%"
# ==============================================================================

# 1. Safely handle NAs in our match flags
meps_merged[is.na(matched_ssr_ndc), matched_ssr_ndc := FALSE]
meps_merged[is.na(matched_ssr_text), matched_ssr_text := FALSE]

# 2. Filter for drugs FDB confirmed are Brands, but failed BOTH SSR matches
unmatched_brands <- meps_merged[
  final_status == "Brand" &
    matched_ssr_ndc == FALSE &
    matched_ssr_text == FALSE,
  .(spend = sum(tot_pay_amt, na.rm = TRUE)),
  by = .(rxname, final_join_name)
][order(-spend)]

cat("\n========================================\n")
cat(" TOP 15 UNMATCHED BRANDS (The Missing 20%) \n")
cat("========================================\n")
print(head(unmatched_brands, 15))

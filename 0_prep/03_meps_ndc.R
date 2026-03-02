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
# STEP 3: PREPARE SSR (Target)
# ==============================================================================

# Load and filter to TOTAL payer (upon Will's advising to not split Medicaid)
ssr <- fread(paste0(data_dir, "processed/ssr_cleaned/ssr_gtn_annual.csv"))[
  payer == 'Total' & smoothing_type == "4q moving average",
][, ssr_rxname := toupper(trimws(product_clean))][, .(
  ssr_rxname,
  disease_area,
  year_id
)]

# ==============================================================================
# STEP 4: LOAD & CLEAN MEPS
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
# STEP 5: MERGE MEPS WITH DICTIONARIES
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
# STEP 5.5: THE FUZZY RESCUE (Fixing the 19% Gap)
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
# STEP 5.8: FLAG THERAPEUTIC CLASSES (NON-MATCHABLE)
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
# STEP 5.9: UPDATE FINAL STATUS (PRIORITIZE FLAGS)
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
# STEP 9: CALCULATE FINAL NET SPENDING
# ==============================================================================

# 1. Initialize Rebate Column
meps_merged[, final_rebate_pct := 0.0]

# 2. Apply Brand Rebates (Placeholder for your specific SSR merge)
#    meps_merged[final_status == "Brand" & !is.na(ssr_specific_rebate),
#                final_rebate_pct := ssr_specific_rebate]

# 3. Apply Class Average Rebates (Imputed High-Cost Classes)
meps_merged[final_status == "Class_Imputed", final_rebate_pct := avg_rebate]

# 4. Explicitly set 0.0 for the non-rebated categories
no_rebate_cats <- c(
  "Generic",
  "Unidentified",
  "Medical_Device",
  "Therapeutic_Class"
)
meps_merged[final_status %in% no_rebate_cats, final_rebate_pct := 0.0]

# 5. Calculate Net Pay
meps_merged[, net_pay_amt := tot_pay_amt * (1 - final_rebate_pct)]

# ==============================================================================
# FINAL SPENDING BREAKDOWN
# ==============================================================================

# Calculate Sums
total_gross <- meps_merged[, sum(tot_pay_amt, na.rm = T)]
total_net <- meps_merged[, sum(net_pay_amt, na.rm = T)]

# Group by Status (Now includes Devices and Classes!)
breakdown <- meps_merged[,
  .(
    gross_spend = sum(tot_pay_amt, na.rm = T),
    net_spend = sum(net_pay_amt, na.rm = T),
    n_records = .N
  ),
  by = final_status
]

# Add Percentages
breakdown[, pct_gross := round(gross_spend / total_gross * 100, 1)]
breakdown[, pct_net := round(net_spend / total_net * 100, 1)]

# Order for readability
setorder(breakdown, -gross_spend)

cat("\n=== FINAL DATASET COMPOSITION ===\n")
print(breakdown)

cat("\n=== TOTALS ===\n")
cat(
  "Total Gross Spend: $",
  formatC(total_gross, format = "d", big.mark = ","),
  "\n"
)
cat(
  "Total Net Spend:   $",
  formatC(total_net, format = "d", big.mark = ","),
  "\n"
)

# Verify Math
if (abs(total_gross - sum(breakdown$gross_spend)) < 100) {
  cat("SUCCESS: 100% of Gross Spending is accounted for.\n")
} else {
  cat("WARNING: Spending dropped during categorization.\n")
}


# ==============================================================================
# STEP 6: THERAPEUTIC CLASS AVERAGES FOR SELECT DRUGS
# ==============================================================================

# Load full SSR data (ensure you have the rebate/discount columns)
# Assuming you have columns: year_id, disease_area, gross_spend, net_spend
# OR a 'discount_pct' column.
# If you don't have gross/net columns, calculate from the discount %

ssr_full <- fread(paste0(data_dir, "processed/ssr_cleaned/ssr_gtn_annual.csv"))[
  payer == 'Total'
]

# Calculate Weighted Average Discount per Disease Area per Year
ssr_class_avgs <- ssr_full[,
  .(
    total_gross = sum(gross_sales, na.rm = T), # Adjust column name as needed
    total_net = sum(net_sales, na.rm = T) # Adjust column name as needed
  ),
  by = .(year_id, disease_area)
]

ssr_class_avgs[, avg_class_discount := 1 - (total_net / total_gross)]

# Handle weird data (caps at 0 and 1)
ssr_class_avgs[avg_class_discount < 0, avg_class_discount := 0]
ssr_class_avgs[avg_class_discount > 1, avg_class_discount := 1]
ssr_class_avgs[is.nan(avg_class_discount), avg_class_discount := 0]

print(head(ssr_class_avgs))

# ---------------------------

# Create a manual mapping table based on your screenshot
# MEPS Text -> SSR Disease Area
class_map <- rbind(
  data.table(meps_text = "IMMUNOLOGIC AGENTS", ssr_area = "Autoimmune"), # High Conf
  data.table(
    meps_text = "SELECTIVE IMMUNOSUPPRESSANTS",
    ssr_area = "Autoimmune"
  ), # High Conf
  data.table(meps_text = "ANTINEOPLASTICS", ssr_area = "Oncology"), # High Conf
  data.table(meps_text = "20 ANTINEOPLASTICS", ssr_area = "Oncology"),
  data.table(
    meps_text = "MISCELLANEOUS ANTINEOPLASTICS",
    ssr_area = "Oncology"
  ),
  data.table(meps_text = "ANTIVIRAL AGENTS", ssr_area = "Infectious Disease"), # Medium Conf
  data.table(
    meps_text = "PSYCHOTHERAPEUTIC AGENTS",
    ssr_area = "Central Nervous System"
  ),
  data.table(
    meps_text = "CENTRAL NERVOUS SYSTEM AGENTS",
    ssr_area = "Central Nervous System"
  )
)

# NOTE: We intentionally Exclude "MISCELLANEOUS AGENTS", "ANTI-INFECTIVES", etc.
# These will default to 0% rebate (Generic).

# ----------------------------

# 1. Identify rows that are Unmatched AND look like a Drug Class
# (You likely ran the 'is_drug_class' flag in previous steps)
unmatched_classes <- meps_merged[
  final_status == "Unidentified" & is_drug_class == TRUE
]

# 2. Match to our "Safe Map"
unmatched_classes <- merge(
  unmatched_classes,
  class_map,
  by.x = "text_match_name",
  by.y = "meps_text",
  all.x = TRUE
)

# 3. Bring in the Rebate % from SSR for the mapped classes
unmatched_classes <- merge(
  unmatched_classes,
  ssr_class_avgs,
  by.x = c("ssr_area", "year_id"),
  by.y = c("disease_area", "year_id"),
  all.x = TRUE
)

# 4. Update the Main Dataset
#    If we found a Class Match, we change status to "Class_Imputed"
#    and assign the discount factor.
meps_merged <- merge(
  meps_merged,
  unmatched_classes[, .(year_id, ndc, text_match_name, avg_class_discount)],
  by = c("year_id", "ndc", "text_match_name"),
  all.x = TRUE,
  suffixes = c("", "_class")
)

# Apply the logic
meps_merged[
  !is.na(avg_class_discount),
  `:=`(
    final_status = "Class_Imputed",
    imputed_rebate = avg_class_discount
  )
]

# 5. DIAGNOSTICS
recovered <- meps_merged[
  final_status == "Class_Imputed",
  sum(tot_pay_amt, na.rm = T)
]
cat(
  "Recovered via Class Averaging: $",
  formatC(recovered, big.mark = ","),
  "\n"
)

# ----------------------------------

# Define Rebate %
# 1. Brand Exact Match: Use SSR specific rebate (Need to merge that in)
# 2. Class Imputed: Use Class Average
# 3. Everything else (Generic/Unidentified): 0% Rebate

meps_merged[, final_rebate_pct := 0.0] # Default to 0

# (Assuming you merged specific SSR rebates earlier into 'ssr_specific_rebate')
# meps_merged[final_status == "Brand", final_rebate_pct := ssr_specific_rebate]

meps_merged[final_status == "Class_Imputed", final_rebate_pct := imputed_rebate]

# Calculate Final Net Price
meps_merged[, net_pay_amt := tot_pay_amt * (1 - final_rebate_pct)]

# ==============================================================================
# STEP 6: DETERMINE FINAL STATUS & JOIN NAME
# ==============================================================================

# 1. Determine Status (Brand vs Generic)
meps_merged[,
  final_status := fcase(
    # Priority 1: FDB NDC
    !is.na(fdb_genind)    , fdb_genind    ,

    # Priority 2: RxNorm NDC
    !is.na(rxnorm_genind) , rxnorm_genind ,

    # Priority 3: Text Dictionary Match
    !is.na(type)          , type          ,

    # Fallback: Assume Generic (Conservative for rebate analysis)
    default = "Generic"
  )
]

# 2. Select the Best Available Name (Raw)
meps_merged[,
  raw_join_name := fcase(
    !is.na(fdb_brand)              , fdb_brand       , # FDB Brand Name is usually best
    !is.na(rxnorm_rxname)          , rxnorm_rxname   , # RxNorm Name is second best
    !is.na(type) & type == "Brand" , text_match_name , # Use the matched text
    default = text_match_name # Fallback
  )
]

# 3. Final Cleaning for SSR Join
#    Apply the aggressive cleaner again to ensure the 'raw_join_name'
#    (which might be "TRI-SPRINTEC 28 DAY") matches SSR ("TRI-SPRINTEC")
meps_merged[, final_join_name := clean_drug_text(raw_join_name)]

# ==============================================================================
# STEP 7: MATCH TO SSR & REPORTING
# ==============================================================================

# 1. Filter: Define the Rebatable Universe (Brands Only)
target_brands <- meps_merged[final_status == "Brand"]

# 2. Prepare for Join
#    Ensure Year is integer in both
target_brands[, year_id := as.integer(year_id)]
ssr[, year_id := as.integer(year_id)]

# 3. Create Fallback Match Logic (The "First Word" Safety Net)
#    Get list of valid SSR names
ssr_lookup <- unique(ssr$ssr_rxname)

target_brands[,
  ssr_match_key := fcase(
    # A. Exact Match
    final_join_name %in% ssr_lookup          , final_join_name          ,

    # B. First Word Match (e.g. "TRI-SPRINTEC" from "TRI-SPRINTEC 28 DAY")
    #    Only use if the First Word is actually a valid drug in SSR
    word(final_join_name, 1) %in% ssr_lookup , word(final_join_name, 1) ,

    # C. No Match found
    default = NA_character_
  )
]

# 4. Perform the Merge
final_df <- merge(
  target_brands,
  ssr,
  by.x = c("ssr_match_key", "year_id"),
  by.y = c("ssr_rxname", "year_id"),
  all.x = TRUE
)

# Flag successful matches
final_df[, is_matched := !is.na(disease_area)]

# ==============================================================================
# FINAL CHECKS
# ==============================================================================

# Metric 1: Total Spending in MEPS
total_meps_spend <- meps_merged[, sum(tot_pay_amt, na.rm = T)]

# Metric 2: Spending identified as BRAND (The Denominator)
brand_spend <- target_brands[, sum(tot_pay_amt, na.rm = T)]

# Metric 3: Spending successfully matched to SSR (The Numerator)
matched_spend <- final_df[is_matched == TRUE, sum(tot_pay_amt, na.rm = T)]

cat("\n========================================\n")
cat(" FINAL MATCH RESULTS \n")
cat("========================================\n")
cat(paste0(
  "Total MEPS Spend:            $",
  formatC(total_meps_spend, format = "d", big.mark = ","),
  "\n"
))
cat(paste0(
  "Identified as Brand:         $",
  formatC(brand_spend, format = "d", big.mark = ","),
  " (",
  round(brand_spend / total_meps_spend * 100, 1),
  "% of Total)\n"
))
cat(paste0(
  "Matched to SSR Health:       $",
  formatC(matched_spend, format = "d", big.mark = ","),
  "\n"
))
cat("----------------------------------------\n")
cat(paste0(
  "MATCH RATE (Of Brand Spend): ",
  round(matched_spend / brand_spend * 100, 2),
  "%\n"
))
cat("========================================\n")

# Top Unmatched Brands (Manual Review List)
# This helps you spot if you need to add anything to your dictionary
top_misses <- final_df[
  is_matched == FALSE,
  .(spend = sum(tot_pay_amt, na.rm = T)),
  by = final_join_name
][order(-spend)]

print("Top 15 Unmatched Brands (Review for Naming Mismatches):")
print(head(top_misses, 15))

# Top "Unidentified" Items (Assumed Generic)
# Check this list to ensure no major brands are slipping through as Generics
top_unidentified <- meps_merged[
  final_status == "Generic" & is.na(fdb_genind) & is.na(rxnorm_genind),
  .(spend = sum(tot_pay_amt, na.rm = T)),
  by = text_match_name
][order(-spend)]

print("Top Unidentified Items (Classified as Generic):")
print(head(top_unidentified, 10))


# Define Rebate %
# 1. Brand Exact Match: Use SSR specific rebate (Need to merge that in)
# 2. Class Imputed: Use Class Average
# 3. Everything else (Generic/Unidentified): 0% Rebate

meps_merged[, final_rebate_pct := 0.0] # Default to 0

# (Assuming you merged specific SSR rebates earlier into 'ssr_specific_rebate')
# meps_merged[final_status == "Brand", final_rebate_pct := ssr_specific_rebate]

meps_merged[final_status == "Class_Imputed", final_rebate_pct := imputed_rebate]

# Calculate Final Net Price
meps_merged[, net_pay_amt := tot_pay_amt * (1 - final_rebate_pct)]


# JUNK FROM HERE DOWN

# ==============================================================================
# 5. MERGE MEPS & RUN CHECKS
# ==============================================================================

# Merge FDB
meps_merged <- merge(meps, fdb, by = 'ndc', all.x = TRUE)
# Merge RxNorm
meps_merged <- merge(meps_merged, rxnorm_map, by = "ndc", all.x = TRUE)

# CHECK 2: How much spending has a valid NDC match in FDB?
spend_total <- meps_merged[, sum(tot_pay_amt, na.rm = TRUE)]
spend_fdb <- meps_merged[!is.na(fdb_brand), sum(tot_pay_amt, na.rm = TRUE)]

print(paste0(
  "Total MEPS Spend: $",
  formatC(spend_total, format = "d", big.mark = ",")
))
print(paste0(
  "Spend successfully linked to FDB: ",
  round(spend_fdb / spend_total * 100, 1),
  "%"
))

# Strip dosages from MEPS names to improve matching (e.g. "LIPITOR 20MG" -> "LIPITOR")
meps_merged[,
  text_match_name := str_remove_all(
    rxname,
    "\\d+\\.?\\d*\\s*(MG|MCG|ML|G|L|TABS|CAPS|UNIT).*"
  )
]
meps_merged[,
  text_match_name := trimws(str_remove_all(text_match_name, "[[:punct:]]"))
]

#Apply the Dictionary to missing rows
meps_merged[
  is.na(fdb_genind) & text_match_name %in% known_generics,
  fdb_genind := "Generic"
]
meps_merged[
  is.na(fdb_genind) & text_match_name %in% known_brands,
  fdb_genind := "Brand"
]


# --- PILE 1: Missing NDC (Use Text Matching) ---

# what percent has valid rxname?
meps_missing_ndc[,
  .(
    n = .N,
    spend = sum(tot_pay_amt, na.rm = TRUE),
    pct_spend = sum(tot_pay_amt, na.rm = TRUE) /
      meps_missing_ndc[, sum(tot_pay_amt, na.rm = TRUE)] *
      100
  ),
  by = .(has_rxname = !is.na(rxname) & rxname != "" & rxname != "-15")
][order(-spend)]

# some are therapeutic classes!
# flag therapeutic category labels
drug_class_pattern <-
  meps_missing_ndc[,
    is_drug_class := grepl(
      drug_class_pattern,
      rxname,
      ignore.case = TRUE
    )
  ]


# check actual drug names still unmatched
x <- meps_missing_ndc[
  is_drug_class == FALSE & is_device == FALSE & is.na(type),
  .(
    spend = sum(tot_pay_amt, na.rm = TRUE)
  ),
  by = rxname
][order(-spend)]
fwrite(x, paste0(data_dir, 'name_matching.csv'))

# A. Clean the MEPS rxname
clean_drug_text <- function(x) {
  x <- toupper(x)
  x <- str_remove_all(x, "\\b\\d+\\.?\\d*\\s*(MG|MCG|ML|G|L|MEQ|UNIT[S]?)\\b")
  x <- str_remove_all(
    x,
    "\\b(TAB|CAP|INJ|SOL|HCL|SODIUM|EQ|ER|XR|DR|SYR|SUSP)\\b"
  )
  x <- str_remove_all(x, "[[:punct:]]")
  x <- str_squish(x)
  return(x)
}

meps_missing_ndc[, clean_rxname := clean_drug_text(rxname)]
meps_missing_ndc[, first_word := toupper(trimws(word(rxname, 1)))]

# B. Map against FDB Dictionary - try clean_rxname first, fall back to first_word
meps_missing_ndc <- merge(
  meps_missing_ndc,
  name_dict,
  by.x = "clean_rxname",
  by.y = "clean_name",
  all.x = TRUE
)

# For still-unmatched rows, try first_word
meps_missing_ndc <- merge(
  meps_missing_ndc,
  name_dict,
  by.x = "first_word",
  by.y = "clean_name",
  all.x = TRUE,
  suffixes = c("", "_fw")
)

# Combine: use clean_rxname match if available, else first_word match
meps_missing_ndc[,
  type := fcase(
    !is.na(type)    , type    ,
    !is.na(type_fw) , type_fw ,
    default = NA_character_
  )
]
meps_missing_ndc[, type_fw := NULL]

# C. Calculate Spending Breakdown
spend_generic <- meps_missing_ndc[
  type == "GENERIC",
  sum(tot_pay_amt, na.rm = TRUE)
]

# Check improvement
meps_missing_ndc[,
  .(
    n = .N,
    spend = sum(tot_pay_amt, na.rm = TRUE),
    pct_of_total = sum(tot_pay_amt, na.rm = TRUE) / total_spend * 100
  ),
  by = type
][order(-spend)]


#-----------------------------------------------------------------------------------

# 2. Identified as Brand (Match these to SSR!)
target_brands <- meps_missing_ndc[type == "BRAND"]

# Match to SSR
target_brands[, in_ssr := clean_rxname %in% ssr$ssr_join_name]

# ==============================================================================
# STEP 4: DIAGNOSE THE REMAINING UNMATCHED
# ==============================================================================

# Look at the top unmatched strings that FDB thinks are Brands
unmatched_brands <- target_brands[
  in_ssr == FALSE,
  .(spend = sum(tot_pay_amt, na.rm = T)),
  by = clean_rxname
][order(-spend)]

head(unmatched_brands, 20)

# If you still have high spending here, run fuzzy matching ONLY on this subset
ssr_names <- unique(ssr$ssr_join_name)
unmatched_brands[,
  best_match := ssr_names[amatch(clean_rxname, ssr_names, maxDist = 2)]
]

print("Top fuzzy candidates:")
print(unmatched_brands[!is.na(best_match), .(clean_rxname, best_match, spend)])

# ==============================================================================
# STEP 5: MATCH DIAGNOSTICS
# ==============================================================================

total_spend <- dt[, sum(tot_pay_amt, na.rm = TRUE)]

# --- PILE 2 overall breakdown ---
pile2_total <- meps_missing_ndc[, sum(tot_pay_amt, na.rm = TRUE)]

# How many rows/spend identified as Generic vs Brand vs Unidentified
meps_missing_ndc[,
  .(
    n = .N,
    spend = sum(tot_pay_amt, na.rm = TRUE),
    pct_of_pile2 = sum(tot_pay_amt, na.rm = TRUE) / pile2_total * 100,
    pct_of_total = sum(tot_pay_amt, na.rm = TRUE) / total_spend * 100
  ),
  by = type
][order(-spend)]

# --- Brand matching to SSR ---
target_brands[,
  .(
    n = .N,
    spend = sum(tot_pay_amt, na.rm = TRUE),
    pct_of_total = sum(tot_pay_amt, na.rm = TRUE) / total_spend * 100
  ),
  by = in_ssr
]

# --- Fuzzy match recovery ---
fuzzy_matches <- unmatched_brands[!is.na(best_match)]
cat("Fuzzy matches found:", nrow(fuzzy_matches), "\n")
cat(
  "Spend recovered via fuzzy:",
  sum(fuzzy_matches$spend, na.rm = TRUE) / total_spend * 100,
  "% of total spend\n"
)

# --- Summary waterfall ---
spend_generic_pct <- spend_generic / total_spend * 100
spend_brand_matched <- target_brands[
  in_ssr == TRUE,
  sum(tot_pay_amt, na.rm = TRUE)
]
spend_brand_matched_pct <- spend_brand_matched / total_spend * 100
spend_unmatched <- target_brands[
  in_ssr == FALSE,
  sum(tot_pay_amt, na.rm = TRUE)
]
spend_unmatched_pct <- spend_unmatched / total_spend * 100

cat("\n=== PILE 2 (Missing NDC) Spend Waterfall ===\n")
cat("Identified as Generic:       ", round(spend_generic_pct, 1), "%\n")
cat("Identified as Brand, in SSR: ", round(spend_brand_matched_pct, 1), "%\n")
cat("Brand, unmatched to SSR:     ", round(spend_unmatched_pct, 1), "%\n")
cat(
  "Unidentified (not in FDB):   ",
  round(
    (pile2_total - spend_generic - spend_brand_matched - spend_unmatched) /
      total_spend *
      100,
    1
  ),
  "%\n"
)


# ---------------------------------------

# Total spend unmatched to SSR across both piles

# Pile 1: has NDC but no SSR match
pile1_unmatched <- meps_valid_ndc[
  !(ndc %in% fdb$NDC),
  sum(tot_pay_amt, na.rm = TRUE)
]

# Pile 2: missing NDC, unidentified in FDB
pile2_unidentified <- meps_missing_ndc[
  is.na(type),
  sum(tot_pay_amt, na.rm = TRUE)
]

# Pile 2: identified as brand but not in SSR
pile2_brand_unmatched <- target_brands[
  in_ssr == FALSE,
  sum(tot_pay_amt, na.rm = TRUE)
]

cat("\n=== OVERALL MEPS SPEND UNMATCHED TO SSR ===\n")
cat(
  "Pile 1 (has NDC, no FDB/SSR match): ",
  round(pile1_unmatched / total_spend * 100, 1),
  "%\n"
)
cat(
  "Pile 2 (no NDC, unidentified):      ",
  round(pile2_unidentified / total_spend * 100, 1),
  "%\n"
)
cat(
  "Pile 2 (brand, not in SSR):         ",
  round(pile2_brand_unmatched / total_spend * 100, 1),
  "%\n"
)
cat(
  "Total unmatched:                    ",
  round(
    (pile1_unmatched + pile2_unidentified + pile2_brand_unmatched) /
      total_spend *
      100,
    1
  ),
  "%\n"
)

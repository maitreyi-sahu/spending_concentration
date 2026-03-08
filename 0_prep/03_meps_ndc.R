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
# STEP 6: Join name
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
# STEP 9.5: CLEANUP & BRIDGING (Fixing the Top 20 Mismatches)
# ==============================================================================

# 1. DEMOTE FAKE BRANDS (The "Garbage" Filter)
#    These are ingredients or vague terms that FDB called "Brand" but are actually Generics.
#    Removing these from 'final_status == Brand' shrinks your denominator and fixes the match rate.
fake_brands <- c(
  "ANTISEPTIC",
  "CARDIOVASCULAR SUPPORT",
  "METHYLPHENIDATE",
  "PRAVASTATIN",
  "AMOXICILLIN",
  "IBUPROFEN",
  "LUBRICANT EYE DROPS",
  "SODIUM CHLORIDE",
  "POTASSIUM CHLORIDE",
  "GABAPENTIN",
  "OMEPRAZOLE"
)

meps_merged[final_join_name %in% fake_brands, final_status := "Generic"]

# 2. THE SUFFIX STRIPPER (Programmatic Fix)
#    Removes "FLEXPEN", "KWIKPEN", "HANDIHALER", "DISCUS", "HFA"
#    This fixes NOVOLOG, HUMALOG, SPIRIVA, ADVAIR automatically.
suffix_pattern <- "\\b(FLEXPEN|KWIKPEN|HANDIHALER|DISKUS|HFA|SOLOSTAR|AUTOINJECTOR|PEN|VIAL|U-100|U-200|U-500)\\b"
meps_merged[
  final_status == "Brand",
  final_join_name := str_remove_all(final_join_name, suffix_pattern)
]
meps_merged[
  final_status == "Brand",
  final_join_name := str_squish(final_join_name)
] # Clean up spaces

# 3. MANUAL BRIDGE (The Specific Fixes)
#    Map the stubborn leftovers to their SSR names.
manual_bridge <- rbind(
  data.table(meps = "HUMULIN N", ssr = "HUMULIN"), # Check SSR for specific formulation
  data.table(meps = "NOVOLIN N", ssr = "NOVOLIN"),
  data.table(meps = "JANUMET", ssr = "JANUMET"), # Often mismatch on XR
  data.table(meps = "JANUMET XR", ssr = "JANUMET"),
  data.table(meps = "METOPROLOL SUCCINATE", ssr = "TOPROL XL"), # Common branded generic issue
  data.table(meps = "ADDERALL XR", ssr = "ADDERALL"),
  data.table(meps = "FARXIGA", ssr = "FARXIGA"), # Ensure spelling matches SSR exactly
  data.table(meps = "SEROQUEL XR", ssr = "SEROQUEL")
)

# Apply Manual Bridge
meps_merged <- merge(
  meps_merged,
  manual_bridge,
  by.x = "final_join_name",
  by.y = "meps",
  all.x = TRUE
)

# Overwrite name if bridged
meps_merged[!is.na(ssr), final_join_name := ssr]
meps_merged[, ssr := NULL]

# ==============================================================================
# PROCEED TO STEP 10 (MATCHING)
# ==============================================================================

# ==============================================================================
# STEP 10: MATCH "BRANDS" TO SSR HEALTH
# ==============================================================================

# 1. Filter for the candidates (The 62.3% of spending)
brand_candidates <- meps_merged[final_status == "Brand"]

# 2. Prepare SSR Data (if not already loaded/prepped)
#    Ensure years are integers and names are upper/trimmed
ssr[, year_id := as.integer(year_id)]
ssr_lookup <- unique(ssr$ssr_rxname) # Vector for fast checking

# 3. Create the Match Key (Exact + First Word Fallback)
brand_candidates[,
  ssr_match_key := fcase(
    # A. Exact Match (Best)
    final_join_name %in% ssr_lookup          , final_join_name          ,

    # B. First Word Match (Safety Net for "BRAND 28 DAY" vs "BRAND")
    word(final_join_name, 1) %in% ssr_lookup , word(final_join_name, 1) ,

    # C. No Match
    default = NA_character_
  )
]

# 4. Merge SSR Data
#    (We use the match key we just created)
brands_matched <- merge(
  brand_candidates,
  ssr,
  by.x = c("ssr_match_key", "year_id"),
  by.y = c("ssr_rxname", "year_id"),
  all.x = TRUE
)

# 5. Flag Success
brands_matched[, in_ssr := !is.na(disease_area)]

# ==============================================================================
# DIAGNOSTICS: HOW GOOD IS THE MATCH?
# ==============================================================================

brand_total_spend <- brands_matched[, sum(tot_pay_amt, na.rm = T)]
brand_match_spend <- brands_matched[in_ssr == TRUE, sum(tot_pay_amt, na.rm = T)]

cat("\n=== BRAND MATCHING PERFORMANCE ===\n")
cat(
  "Total Spending identified as 'Brand':      $",
  formatC(brand_total_spend, format = "d", big.mark = ","),
  "\n"
)
cat(
  "Spending successfully matched to SSR:      $",
  formatC(brand_match_spend, format = "d", big.mark = ","),
  "\n"
)
cat(
  "SSR MATCH RATE (Dollars):                  ",
  round(brand_match_spend / brand_total_spend * 100, 1),
  "%\n"
)

# ==============================================================================
# INSPECT THE LOSERS (Brands we couldn't find in SSR)
# ==============================================================================

# Look at the specific names that FDB called "Brand" but SSR didn't find
unmatched_brands <- brand_candidates[
  is.na(ssr_match_key),
  .(spend = sum(tot_pay_amt, na.rm = T)),
  by = final_join_name
][order(-spend)]

print("Top 20 Unmatched Brands:")
print(head(unmatched_brands, 20))


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

#  Goal: Match MEPS data on RX with Rebates data, using redbook
#
#



# Load in MEPS data, which has NDC codes only
MEPS <-  open_dataset("/mnt/share/dex/us_county/00_data_prep/data/MEPS/stage_3/USA_MEPS_RX.parquet") %>% filter() %>% as.data.table() %>% collect()


# Read in redbook, which has NDC code and drug names
redbook <- "/ihme/limited_use/IDENT/PROJECT_FOLDERS/USA/MARKETSCAN/DATA/USA_MARKETSCAN_2022_REDBOOK_Y2024M02D07.SAS7BDAT"  

marketscan_dir <- "/mnt/share/limited_use/IDENT/PROJECT_FOLDERS/USA/MARKETSCAN/DATA/PARQUET/"
redbook <- open_dataset(paste0(marketscan_dir, "redbook.parquet")) %>%  
  select(NDCNUM, GENNME, THRCLDS, GNINDDS, PRODNME) %>% 
  collect() %>% 
  select(NDCNUM, GENNME,PRODNME) %>%  arrange(PRODNME,GENNME, NDCNUM) %>% 
  distinct(NDCNUM, .keep_all = T) %>% setDT()

# Merge MEPS and redbook by NDC
MEPS <- merge(MEPS, redbook, by.x = 'ndc', by.y = 'NDCNUM', all.x = T)

# 91.8% of MEPS encounters get matched
nrow(MEPS[!is.na(GENNME)])/nrow(MEPS)
nrow(MEPS[is.na(GENNME)])/nrow(MEPS)

# Filter to just those matched
MEPS_matched <- MEPS[!is.na(GENNME) | !is.na(PRODNME)]

# Collapse by spend, n_claims, and drug name
frequency <- MEPS_matched[,.(n_claims = length(unique(claim_id)),
                     spend = sum(tot_pay_amt)), by = c('GENNME', 'PRODNME')]

frequency[rev(order(spend))]


# How many drugs are we working with now? 9194, but just 2321 generic
drugs_in_MEPS <- unique(MEPS_matched[,.(full_drug = PRODNME, generic = GENNME)])
drugs_in_MEPS[, drug := tstrsplit(full_drug, " ", fixed = TRUE)[[1]]]
length(unique(drugs_in_MEPS$drug))
length(unique(drugs_in_MEPS$generic))
#drugs_in_MEPS[drug %like% "-", drug := tstrsplit(drug, "-", fixed = TRUE)[[1]]] # this makes us lose some matches


# Rebates info
map1 <- fread('/mnt/share/dex/us_county/projects/meps_concentration_rebate/inputs/ssr_health_dummy/GrossToNet_TotalLessMedicaid_Smoothed_DUMMY.csv')
drugs_in_rebates <- map1[,.(drug = toupper(V1), have_rebates = 1)]

# 692 drug perfect match! 
test <- merge(drugs_in_MEPS, drugs_in_rebates, by = 'drug')

# Look at complete merge
all_matching <- merge(drugs_in_MEPS, drugs_in_rebates, by = 'drug', all.x = T)
all_matching[is.na(have_rebates)]


# need to match drug to drug name..........
drugs_in_rebates[drug %in% toupper(c('Adalimumab', 'Humira', 'Amgevita', 'Hyrimoz'))]


check_top_15_drugs <- frequency[rev(order(spend))]$PRODNME[1:50]
all_matching[full_drug %in% check_top_15_drugs]



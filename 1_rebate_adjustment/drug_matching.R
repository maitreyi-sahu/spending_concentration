#  Goal: Match MEPS data on RX with Rebates data, using redbook

# For MEPS NDCs, Use RxNorm to get brand names, etc
# Create brand flag in MEPS --> APPLY ONLY TO BRAND DRUGS
# Clean SSR Health
# Will only be able to apply the discounts to some drugs and not all

# ------------------------------------------------------------------

source('init.R')
library(httr)
library(jsonlite)

# -------------------------------------------------------------------

# Load in MEPS data, which has NDC codes only
meps <- open_dataset(paste0(data_dir, "meps/USA_MEPS_RX.parquet")) %>%
  filter(toc == 'RX' & year_id >= 2000) |>
  as.data.table() %>%
  collect()

# Clean NDCs
print(unique(frequency$ndc[!grepl("^\\d{11}$", frequency$ndc)]))
missing_codes <- -c(
  '-15' ~ `Cannot be computed`,
  0022.95
)


# Collapse by spend, n_claims, and drug name
frequency <- meps[,
  .(n_claims = length(unique(claim_id)), spend = sum(tot_pay_amt)),
  by = 'ndc'
]
frequency[rev(order(spend))]


# How many drugs are we working with now? 9194, but just 2321 generic
drugs_in_MEPS <- unique(MEPS_matched[, .(
  full_drug = PRODNME,
  generic = GENNME
)])
drugs_in_MEPS[, drug := tstrsplit(full_drug, " ", fixed = TRUE)[[1]]]
length(unique(drugs_in_MEPS$drug))
length(unique(drugs_in_MEPS$generic))
#drugs_in_MEPS[drug %like% "-", drug := tstrsplit(drug, "-", fixed = TRUE)[[1]]] # this makes us lose some matches

# Rebates info
map1 <- fread(
  '/mnt/share/dex/us_county/projects/meps_concentration_rebate/inputs/ssr_health_dummy/GrossToNet_TotalLessMedicaid_Smoothed_DUMMY.csv'
)
drugs_in_rebates <- map1[, .(drug = toupper(V1), have_rebates = 1)]

# 692 drug perfect match!
test <- merge(drugs_in_MEPS, drugs_in_rebates, by = 'drug')

# Look at complete merge
all_matching <- merge(drugs_in_MEPS, drugs_in_rebates, by = 'drug', all.x = T)
all_matching[is.na(have_rebates)]


# need to match drug to drug name..........
drugs_in_rebates[
  drug %in% toupper(c('Adalimumab', 'Humira', 'Amgevita', 'Hyrimoz'))
]


check_top_15_drugs <- frequency[rev(order(spend))]$PRODNME[1:50]
all_matching[full_drug %in% check_top_15_drugs]

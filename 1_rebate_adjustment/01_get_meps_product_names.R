# ==============================================================================
# RxNorm API Matching: Get proudct names for MEPS NDCs
# [TAKES ~XX MIN TO RUN]
#
# msahu
# ==============================================================================

source('init.R')

library(httr)
library(jsonlite)
library(fuzzyjoin)

# ==============================================================================
# Load cleaned data
# ==============================================================================

# Read MEPS NDCs and order by spend [starting list]
meps_ndcs <- fread(
  paste0(data_dir, "processed/meps_cleaned/ndc_lookup.csv"),
  colClasses = list(character = c("ndc_original", "ndc_clean"))
)
meps_ndcs <- meps_ndcs[ndc_quality_flag != 'missing_code'][order(-total_spend)]


# Read SSR Health data with Product Names
ssr_names <- data.table(unique(fread(paste0(
  data_dir,
  "processed/ssr_cleaned/ssr_gtn_annual.csv"
))[, .(product, product_clean)]))

# ==============================================================================
# RxNorm API Helper Function
# ==============================================================================

library(httr)
library(jsonlite)

get_rxnorm_from_ndc <- function(ndc) {
  Sys.sleep(0.05)

  base_url <- "https://rxnav.nlm.nih.gov/REST"
  ndc_clean <- gsub("-", "", as.character(ndc))

  tryCatch(
    {
      # Get RxCUI
      response <- GET(paste0(base_url, "/rxcui.json?idtype=NDC&id=", ndc_clean))
      content <- fromJSON(content(response, "text"))

      if (is.null(content$idGroup$rxnormId)) {
        return(data.table(
          ndc = ndc,
          rxcui = NA_character_,
          drug_name = NA_character_,
          tty = NA_character_,
          brand_generic = NA_character_,
          ingredient = NA_character_
        ))
      }

      rxcui <- content$idGroup$rxnormId[[1]]

      # Get properties
      props <- GET(paste0(
        base_url,
        "/rxcui/",
        rxcui,
        "/allProperties.json?prop=all"
      ))
      props_content <- fromJSON(content(props, "text"))
      props_df <- props_content$propConceptGroup$propConcept

      drug_name <- props_df[props_df$propName == "RxNorm Name", "propValue"][1]
      tty <- props_df[props_df$propName == "TTY", "propValue"][1]

      # Get ingredient
      ing <- GET(paste0(base_url, "/rxcui/", rxcui, "/related.json?tty=IN"))
      ing_content <- fromJSON(content(ing, "text"))
      ingredient <- NA_character_

      if (
        !is.null(ing_content$relatedGroup$conceptGroup$conceptProperties[[1]])
      ) {
        ingredient <- paste(
          ing_content$relatedGroup$conceptGroup$conceptProperties[[1]]$name,
          collapse = " / "
        )
      }

      brand_generic <- fcase(
        tty %in% c("SBD", "BPCK", "GPCK") , "Brand"   ,
        tty %in% c("SCD", "SCDC", "SCDF") , "Generic" ,
        default = NA_character_
      )

      return(data.table(
        ndc = ndc,
        rxcui = rxcui,
        drug_name = drug_name,
        tty = tty,
        brand_generic = brand_generic,
        ingredient = ingredient
      ))
    },
    error = function(e) {
      return(data.table(
        ndc = ndc,
        rxcui = NA_character_,
        drug_name = NA_character_,
        tty = NA_character_,
        brand_generic = NA_character_,
        ingredient = NA_character_
      ))
    }
  )
}

# ==============================================================================
# Process ALL NDCs with parallel
# ==============================================================================

library(parallel)
start_time = Sys.time()

n_cores <- detectCores() - 1
batch_size <- 5000
ndc_list <- meps_ndcs$ndc_clean
n_batches <- ceiling(length(ndc_list) / batch_size)

all_results <- list()

for (b in 1:n_batches) {
  batch_file <- paste0(data_dir, "processed/rxnorm_mapped/batch_", b, ".rds")

  if (file.exists(batch_file)) {
    all_results[[b]] <- readRDS(batch_file)
    next
  }

  start_idx <- (b - 1) * batch_size + 1
  end_idx <- min(b * batch_size, length(ndc_list))
  batch_ndcs <- ndc_list[start_idx:end_idx]

  cl <- makeCluster(n_cores)
  clusterEvalQ(cl, {
    library(httr)
    library(jsonlite)
    library(data.table)
  })
  clusterExport(cl, "get_rxnorm_from_ndc")

  batch_results <- parLapply(cl, batch_ndcs, get_rxnorm_from_ndc)
  stopCluster(cl)

  batch_results <- rbindlist(batch_results, fill = T)
  saveRDS(batch_results, batch_file)
  all_results[[b]] <- batch_results
}

# Bind and add column for product name
rxnorm_mapping <- rbindlist(all_results, fill = T)
rxnorm_mapping[, product_name := sub(".*\\[(.*)\\].*", "\\1", drug_name)]
rxnorm_mapping[!grepl("\\[", drug_name), product_name := NA_character_]

# Save
saveRDS(
  rxnorm_mapping,
  paste0(data_dir, "processed/rxnorm_mapped/meps_ndc_to_brand_full.rds")
)

cat("\n=== COMPLETE ===\n")
cat(
  "Total mapped:",
  sum(!is.na(rxnorm_mapping$rxcui)),
  "of",
  nrow(rxnorm_mapping),
  "(",
  round(100 * sum(!is.na(rxnorm_mapping$rxcui)) / nrow(rxnorm_mapping), 1),
  "%)\n"
)

time = Sys.time() - start_time
print(time)

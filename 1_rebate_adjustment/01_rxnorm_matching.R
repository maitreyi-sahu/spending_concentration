# ==============================================================================
# RxNorm API Matching: Link MEPS NDCs to SSR Health Brand Names
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
  base_url <- "https://rxnav.nlm.nih.gov/REST"

  Sys.sleep(0.05)

  ndc_clean <- gsub("-", "", as.character(ndc))

  tryCatch(
    {
      # Get RxCUI
      response <- GET(paste0(base_url, "/rxcui.json?idtype=NDC&id=", ndc_clean))

      if (status_code(response) != 200) {
        return(data.table(
          ndc = ndc,
          rxcui = NA_character_,
          drug_name = NA_character_,
          tty = NA_character_,
          brand_generic = NA_character_,
          ingredient = NA_character_
        ))
      }

      content <- fromJSON(content(response, "text", encoding = "UTF-8"))

      if (
        is.null(content$idGroup$rxnormId) ||
          length(content$idGroup$rxnormId) == 0
      ) {
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
      props_response <- GET(paste0(
        base_url,
        "/rxcui/",
        rxcui,
        "/allProperties.json?prop=all"
      ))

      drug_name <- NA_character_
      tty <- NA_character_

      if (status_code(props_response) == 200) {
        props_content <- fromJSON(content(
          props_response,
          "text",
          encoding = "UTF-8"
        ))

        if (!is.null(props_content$propConceptGroup$propConcept)) {
          props_df <- props_content$propConceptGroup$propConcept

          # Get drug name
          name_row <- props_df[props_df$propName == "RxNorm Name", ]
          if (nrow(name_row) > 0) {
            drug_name <- name_row$propValue[1]
          }

          # Get TTY
          tty_row <- props_df[props_df$propName == "TTY", ]
          if (nrow(tty_row) > 0) {
            tty <- tty_row$propValue[1]
          }
        }
      }

      # Get ingredients - FIX THIS PART
      ingredient <- NA_character_
      ing_response <- GET(paste0(
        base_url,
        "/rxcui/",
        rxcui,
        "/related.json?tty=IN"
      ))

      if (status_code(ing_response) == 200) {
        ing_content <- fromJSON(content(
          ing_response,
          "text",
          encoding = "UTF-8"
        ))

        # Access the conceptGroup data.frame
        if (!is.null(ing_content$relatedGroup$conceptGroup)) {
          concept_df <- ing_content$relatedGroup$conceptGroup

          # Check if conceptProperties exists and has data
          if (
            "conceptProperties" %in%
              names(concept_df) &&
              !is.null(concept_df$conceptProperties[[1]]) &&
              nrow(concept_df$conceptProperties[[1]]) > 0
          ) {
            ing_df <- concept_df$conceptProperties[[1]]
            ingredient <- paste(ing_df$name, collapse = " / ")
          }
        }
      }

      # Determine brand/generic
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
      message("Error processing NDC ", ndc, ": ", e$message)
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

n_cores <- detectCores() - 1
batch_size <- 5000
ndc_list <- meps_ndcs$ndc_clean[1:10] # ALL 77,680 NDCs
n_batches <- ceiling(length(ndc_list) / batch_size)

cat(
  "Processing",
  length(ndc_list),
  "NDCs in",
  n_batches,
  "batches using",
  n_cores,
  "cores\n\n"
)

all_results <- list()

for (b in 1:n_batches) {
  batch_file <- paste0(data_dir, "processed/rxnorm_mapped/batch_", b, ".rds")

  if (file.exists(batch_file)) {
    cat("Batch", b, "- loading existing\n")
    all_results[[b]] <- readRDS(batch_file)
    next
  }

  start_idx <- (b - 1) * batch_size + 1
  end_idx <- min(b * batch_size, length(ndc_list))
  batch_ndcs <- ndc_list[start_idx:end_idx]

  cat("Batch", b, "- processing", length(batch_ndcs), "NDCs...")
  batch_start <- Sys.time()

  cl <- makeCluster(n_cores)
  clusterEvalQ(cl, {
    library(httr)
    library(jsonlite)
    library(data.table)
  })
  clusterExport(cl, "get_rxnorm_from_ndc")

  batch_results <- parLapply(cl, batch_ndcs, get_rxnorm_from_ndc)
  stopCluster(cl)

  batch_results <- rbindlist(batch_results)
  saveRDS(batch_results, batch_file)
  all_results[[b]] <- batch_results

  cat(
    " done in",
    round(difftime(Sys.time(), batch_start, units = "mins"), 1),
    "min",
    "- mapped:",
    sum(!is.na(batch_results$rxcui)),
    "\n"
  )
}

rxnorm_mapping <- rbindlist(all_results)
saveRDS(
  rxnorm_mapping,
  paste0(data_dir, "processed/rxnorm_mapped/meps_ndc_to_brand_full.rds")
)

cat(
  "\nDONE! Mapped",
  sum(!is.na(rxnorm_mapping$rxcui)),
  "of",
  nrow(rxnorm_mapping),
  "NDCs\n"
)

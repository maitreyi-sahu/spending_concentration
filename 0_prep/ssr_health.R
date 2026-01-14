source('init.R')

ssr_dir <- paste0(data_dir, 'raw/ssr_health/')
ssr_files <- list.files(ssr_dir, pattern = "\\.csv$", full.names = TRUE)

# Read and combine all files
ssr_dt <- rbindlist(
  lapply(ssr_files, function(f) {
    dt <- as.data.table(read.delim(f, fileEncoding = "UTF-16LE"))
    setnames(
      dt,
      old = c("X", "X.1", "X.2"),
      new = c("product", "disease_area", "smoothing_type")
    )
    dt[, source_file := basename(f)]
    dt
  }),
  use.names = T
)
ssr_dt[ssr_dt == ""] <- NA

ssr_dt[, lapply(.SD, function(x) sum(is.na(x)))] / nrow(ssr_dt)

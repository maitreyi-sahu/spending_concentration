user = 'MS3390'
local = T

# Avoid scientific notation in all printed/saved outputs
options(scipen = 999)

# ---------------------- ARGUMENTS -------------------------------

START_YEAR = 2000
END_YEAR = 2023

# -------------------- DIRECTORIES -------------------------------

dir_head <- if (Sys.info()[["sysname"]] == "Darwin") "/Users/" else "C:/Users/"

if (local) {
  dir <- paste0(dir_head, user, '/projects_local/concentration_local/')
} else {
  dir <- paste0(
    dir_head,
    user,
    '/Partners HealthCare Dropbox/Maitreyi Sahu/projects/concentration_local/'
  )
}

data_dir <- paste0(dir, 'data/')
plot_dir <- paste0(dir, 'plots/')

# ---------------------- PACKAGES -------------------------------

pacman::p_load(
  data.table,
  readxl,
  openxlsx,
  arrow,
  dplyr,
  ggplot2,
  RColorBrewer,
  gridExtra
)

# ----------------- HELPER FUNCTIONS ----------------------------

ensure_dir <- function(path) {
  target <- if (grepl("\\.[^/]+$", path)) dirname(path) else path
  dir.create(target, recursive = TRUE, showWarnings = FALSE)
}

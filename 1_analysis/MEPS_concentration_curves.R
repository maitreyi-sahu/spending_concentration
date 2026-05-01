# ---------------------------
# Purpose: Examine inequality in RX spending, compared to other tocs
#
# Authors: Hkl1 & msahu
# ---------------------------

source('init.R')

# --------------------------------
#   SET ARGS AND DIRECTORY TO SAVE TO
# --------------------------------

use_wgts <- T
include_sample_pop <- T

save_dir_base <- paste0(dir, "plots/")
year_range <- c(START_YEAR, END_YEAR)
iso_years <- c(2008, 2023) # Lorenz comparison years (must exist in both gross & net data)

# MEPS PATHS
denoms_path <- paste0(data_dir, "raw/meps_denoms/MEPS_sample_denoms.parquet")
meps_data_path <- paste0(data_dir, "raw/meps/")

# Rebate-Adjusted MEPS RX path
adjusted_rx_path <- paste0(
  data_dir,
  "processed/meps_rebate_adjusted/meps_rebate_adjusted.rds"
)

# --------------------------------
#  LOAD SAMPLE DENOMS AND NON-RX DATA (Do this once!)
# --------------------------------
t0 <- Sys.time()

sample_denoms <- open_dataset(denoms_path) %>% collect() %>% as.data.table()
sample_denoms <- sample_denoms[,
  .(pop = sum(pop), n_obs = sum(n_obs)),
  by = c('year_id')
]

all_raw_files <- list.files(
  meps_data_path,
  pattern = "\\.parquet$",
  full.names = T
)
non_rx_files <- all_raw_files[!grepl("USA_MEPS_RX\\.parquet$", all_raw_files)]
rx_file <- all_raw_files[grepl("USA_MEPS_RX\\.parquet$", all_raw_files)]

data_non_rx <- rbindlist(lapply(non_rx_files, read_parquet), fill = T)
data_rx_raw <- read_parquet(rx_file) %>% as.data.table()
setnames(data_rx_raw, "tot_pay_amt", "gross_pay_amt")
data_rx_base <- readRDS(adjusted_rx_path) %>% as.data.table()

print("Base Data loaded!")
print(Sys.time() - t0)

# --------------------------------
# SETUP THE 4 SCENARIOS
# --------------------------------

scenarios <- list(
  list(
    name = "gross_spending",
    rx_col = "gross_pay_amt",
    non_rx_col = "oop_pay_amt",
    title_suffix = "Gross unadjusted spending",
    start_year = START_YEAR,
    use_rebate_adjusted = FALSE
  ),
  list(
    name = "oop_gross_spending",
    rx_col = "oop_pay_amt",
    non_rx_col = "gross_pay_amt",
    title_suffix = "Gross out-of-pocket spending",
    start_year = START_YEAR,
    use_rebate_adjusted = FALSE
  ),
  list(
    name = "net_spending_main",
    rx_col = "net_pay_amt_total",
    non_rx_col = "gross_pay_amt",
    title_suffix = "Net spending",
    start_year = 2008,
    use_rebate_adjusted = TRUE
  ),
  list(
    name = "net_spending_mdcd_split",
    rx_col = "net_pay_amt_split",
    non_rx_col = "gross_pay_amt",
    title_suffix = "Net spending - Medicaid split",
    start_year = 2008,
    use_rebate_adjusted = TRUE
  )
)

# LOOP THROUGH EACH SCENARIO
for (scenario in scenarios) {
  message(paste0("\n=== RUNNING SCENARIO: ", scenario$name, " ==="))

  # Create a dedicated directory for this scenario's plots
  scenario_dir <- paste0(save_dir_base, scenario$name, "/")
  if (!dir.exists(scenario_dir)) {
    dir.create(scenario_dir, recursive = TRUE)
  }

  # 1. Prepare Data for this scenario
  # Pick the right RX dataset based on whether this scenario needs rebate adjustment
  if (scenario$use_rebate_adjusted) {
    data_rx <- copy(data_rx_base)
  } else {
    data_rx <- copy(data_rx_raw)
  }
  data_rx[, tot_pay_amt := get(scenario$rx_col)]

  # Copy non-rx data and assign the correct comparison column
  data_non_rx_temp <- copy(data_non_rx)

  # Check if the column exists in non_rx data before assigning
  if (scenario$non_rx_col %in% names(data_non_rx_temp)) {
    data_non_rx_temp[, tot_pay_amt := get(scenario$non_rx_col)]
  } else {
    warning(paste("Column", scenario$non_rx_col, "not found in data_non_rx!"))
  }

  # Now bind them together
  data <- rbindlist(list(data_non_rx_temp, data_rx), fill = T)

  # 2. Filtering & Detrunc Mimic
  drop_cols <- c(
    "toc",
    "age_group_years_start",
    "age_group_id",
    "pri_payer",
    "sex_id",
    "code_system",
    "dx_level",
    "dx"
  )
  na_vals <- list(-1, "-1", "None", "NA", "<NA>", "UNK", "UNKNOWN", "unknown")

  for (j in drop_cols) {
    data <- data[(!is.na(get(j)) & !get(j) %in% na_vals) | toc == 'DV']
  }

  data <- data[primary_cause == 1 | dx_level == 'dx_1' | toc == 'DV']
  data <- data[year_id >= scenario$start_year & year_id <= max(year_range)]
  data <- data[!is.na(tot_pay_amt)]

  data <- data[, .(
    year_id,
    bene_id,
    age_group_years_start,
    sex_id,
    toc,
    encounter_id = claim_id,
    pay = tot_pay_amt,
    survey_wt
  )]

  # 3. Clean Data and Calculate Percentiles
  data[, survey_wt := as.numeric(as.character(survey_wt))]

  data_t <- data[,
    .(pay = sum(pay)),
    by = c('bene_id', 'toc', 'survey_wt', 'year_id')
  ]
  data_nt <- data[,
    .(pay = sum(pay)),
    by = c('bene_id', 'survey_wt', 'year_id')
  ][, toc := 'all']

  data <- rbind(data_t, data_nt)
  data <- data[order(pay)]

  if (use_wgts == T) {
    data[, tot_wt := sum(survey_wt), by = c('toc', 'year_id')]
    if (include_sample_pop == T) {
      sample_w_encounters <- unique(data[, .(tot_wt, toc, year_id)])
      sample_w_encounters <- merge(
        sample_w_encounters,
        sample_denoms,
        by = c('year_id')
      )
      sample_w_encounters[, survey_wt := pop - tot_wt]
      sample_w_encounters <- sample_w_encounters[, .(
        year_id,
        toc,
        survey_wt,
        pay = 0,
        bene_id = 'no-encounters'
      )]
      data <- rbind(data, sample_w_encounters, fill = T)
      data <- data[order(pay)]
      data[, tot_wt := sum(survey_wt), by = c('toc', 'year_id')]
    }
    data[, wt_cumsum := cumsum(survey_wt), by = c('toc', 'year_id')]
    data[, percentile := (wt_cumsum / tot_wt) * 100]
  }

  plot_df <- data[, .(pay = sum(pay)), by = c('percentile', 'toc', 'year_id')]
  plot_df <- plot_df[order(percentile)]
  plot_df[, pay_cumsum := cumsum(pay), by = c('toc', 'year_id')]
  plot_df[, tot_spend := max(pay_cumsum), by = c('toc', 'year_id')]
  plot_df[, percent_tot_spend := (pay_cumsum / tot_spend) * 100]

  save_plot_df <- plot_df[, .(year_id, toc, percentile, percent_tot_spend)]

  closest_percent <- rbindlist(
    lapply(0:100, function(p) {
      use_data <- copy(save_plot_df)
      if (p == 0) {
        use_data[, min_percentile := min(percentile), by = c('year_id', 'toc')]
        return(use_data[percentile == min_percentile, ])
      } else if (p == 100) {
        use_data[, max_percentile := max(percentile), by = c('year_id', 'toc')]
        return(use_data[percentile == max_percentile, ])
      } else {
        use_data[, dist_from_p := abs(p - percentile)]
        use_data[, min_dist := min(dist_from_p), by = c('year_id', 'toc')]
        return(use_data[dist_from_p == min_dist, ])
      }
    }),
    fill = T
  )

  closest_percent[, `:=`(
    min_percentile = NULL,
    dist_from_p = NULL,
    min_dist = NULL,
    max_percentile = NULL
  )]
  closest_percent <- unique(closest_percent)

  # Number plug printout
  check <- closest_percent[
    toc == 'RX' &
      year_id == max(year_id) &
      percentile > 89.9999 &
      percentile < 90.1
  ]
  cat(paste0(
    "\n[",
    scenario$name,
    "] Top 10% of users in MEPS account for ",
    round(100 - check$percent_tot_spend, digits = 2),
    "% of total spending\n"
  ))

  top10_trend <- closest_percent[,
    .(
      top10_share = 100 -
        approx(x = percentile, y = percent_tot_spend, xout = 90, rule = 2)$y
    ),
    by = .(toc, year_id)
  ]

  top10_trend[,
    toc_plot := plyr::revalue(
      toc,
      c(
        "AM" = "Ambulatory",
        "DV" = 'Dental',
        "all" = "All types combined",
        'ED' = 'Emergency department',
        'IP' = 'Inpatient',
        'HH' = 'Home health',
        'RX' = 'Retail prescription drugs'
      )
    )
  ]

  # ==========================================================
  # SAVE SNAPSHOT FOR NUMBER PLUGGING (Outside the loop)
  # ==========================================================
  if (scenario$name == "net_spending_main") {
    top10_trend_main <- copy(top10_trend)
    closest_percent_main <- copy(closest_percent)
  }
  if (scenario$name == "gross_spending") {
    top10_trend_gross <- copy(top10_trend)
  }

  # ==========================================================
  # PLOT: LORENZ CURVES
  # ==========================================================

  plot_lorenz_curve_data <- plot_df[toc == 'RX']
  add_start_val <- unique(plot_lorenz_curve_data[, .(year_id, toc)])
  add_start_val[, `:=`(percent_tot_spend = 0, percentile = 0)]
  plot_lorenz_curve_data <- rbind(
    plot_lorenz_curve_data,
    add_start_val,
    fill = T
  )

  ribbon_x <- seq(0, 100, length.out = 1000)
  ribbon_y_early <- approx(
    x = plot_lorenz_curve_data[year_id == min(iso_years)]$percentile,
    y = plot_lorenz_curve_data[year_id == min(iso_years)]$percent_tot_spend,
    xout = ribbon_x,
    rule = 2
  )$y
  ribbon_y_late <- approx(
    x = plot_lorenz_curve_data[year_id == max(iso_years)]$percentile,
    y = plot_lorenz_curve_data[year_id == max(iso_years)]$percent_tot_spend,
    xout = ribbon_x,
    rule = 2
  )$y
  ribbon_data <- data.table(
    percentile = ribbon_x,
    y_early = ribbon_y_early,
    y_late = ribbon_y_late
  )

  point_90 <- plot_lorenz_curve_data[
    year_id %in% iso_years,
    .(
      y_val = approx(
        x = percentile,
        y = percent_tot_spend,
        xout = 90,
        rule = 2
      )$y
    ),
    by = year_id
  ]
  point_90[, top10_share := paste0(sprintf("%.1f", 100 - y_val), "%")]
  point_90[, combined_label := paste0(year_id, ": ", top10_share)]

  pdf(
    paste0(scenario_dir, "lorenz_curves_overlay_gap.pdf"),
    width = 7,
    height = 5
  )
  p1 <- ggplot() +
    geom_abline(intercept = 0, slope = 1, linetype = 3, color = "grey50") +
    geom_ribbon(
      data = ribbon_data,
      aes(x = percentile, ymin = y_early, ymax = y_late),
      fill = "#8B008B",
      alpha = 0.12
    ) +
    geom_line(
      data = plot_lorenz_curve_data[year_id %in% iso_years],
      aes(
        x = percentile,
        y = percent_tot_spend,
        color = factor(year_id),
        group = year_id
      ),
      size = 1.2
    ) +
    geom_vline(xintercept = 90, color = "grey70", linetype = 1) +
    geom_segment(
      data = point_90,
      aes(x = 0, xend = 90, y = y_val, yend = y_val, color = factor(year_id)),
      linetype = 2,
      alpha = 0.7,
      size = 0.6
    ) +
    geom_text(
      data = point_90,
      aes(x = 70, y = y_val, label = combined_label, color = factor(year_id)),
      vjust = -0.5,
      hjust = 0.5,
      fontface = "bold",
      size = 5,
      show.legend = FALSE
    ) +
    annotate(
      "text",
      x = 5,
      y = 95,
      label = paste0(
        "The top 10% of highest-spending patients accounted for\n",
        point_90[year_id == min(iso_years), top10_share],
        " of spending in ",
        min(iso_years),
        ", growing to ",
        point_90[year_id == max(iso_years), top10_share],
        " in ",
        max(iso_years),
        "."
      ),
      hjust = 0,
      vjust = 1,
      size = 4.5,
      fontface = "italic",
      color = "grey30"
    ) +
    scale_color_manual(values = c("#377EB8", "#E41A1C")) +
    theme_classic() +
    theme(
      legend.position = "none",
      text = element_text(size = 12),
      panel.grid.major = element_line(color = "grey90", size = 0.3)
    ) +
    scale_x_continuous(breaks = seq(0, 100, by = 10), expand = c(0.01, 0)) +
    scale_y_continuous(breaks = seq(0, 100, by = 10), expand = c(0.01, 0)) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 100), clip = "off")

  # Determine the title: Full title for main paper, short title for supplements
  lorenz_title <- if (scenario$name == "net_spending_main") {
    "Increasing Concentration of Retail Prescription Drug Spending"
  } else {
    scenario$title_suffix
  }

  p1 <- p1 +
    labs(
      x = 'Cumulative percentage of US population',
      y = "Cumulative percentage of retail prescription drug spending",
      title = lorenz_title
    )

  print(p1)
  dev.off()

  # ==========================================================
  # PLOT: SINGLE PANEL TOP 10% SHARE
  # ==========================================================
  library(ggrepel)
  single_plot_data <- copy(top10_trend[toc != 'all'])
  single_plot_data[, is_rx := ifelse(toc == 'RX', "Yes", "No")]

  grey_color = 'grey80'
  highlight_colors <- c(
    "Retail prescription drugs" = "#00008B",
    "Ambulatory" = grey_color,
    "Dental" = grey_color,
    "Emergency department" = grey_color,
    "Inpatient" = grey_color
  )

  single_plot_data <- single_plot_data[toc_plot != "Home health"]
  max_yr <- max(single_plot_data$year_id)
  label_data <- single_plot_data[year_id == max_yr]

  pdf(paste0(scenario_dir, "single_panel_clean.pdf"), width = 8, height = 5)
  p2 <- ggplot(
    single_plot_data,
    aes(x = year_id, y = top10_share, color = toc_plot, group = toc_plot)
  ) +
    geom_point(aes(alpha = is_rx, size = is_rx)) +
    geom_smooth(method = 'lm', se = FALSE, aes(size = is_rx)) +
    geom_text_repel(
      data = label_data,
      aes(
        x = max_yr,
        y = top10_share,
        label = ifelse(
          toc_plot == "Retail prescription drugs",
          "Retail prescription\ndrugs",
          as.character(toc_plot)
        )
      ),
      nudge_x = 1.5,
      direction = "y",
      hjust = 0,
      segment.size = 0.4,
      segment.color = "grey80",
      size = 4.5,
      fontface = ifelse(
        label_data$toc_plot == "Retail prescription drugs",
        "bold",
        "plain"
      )
    ) +
    scale_color_manual(values = highlight_colors) +
    scale_alpha_manual(values = c("Yes" = 1, "No" = 0.7), guide = "none") +
    scale_size_manual(values = c("Yes" = 1, "No" = 0.8), guide = "none") +
    theme_classic() +
    theme(
      legend.position = "none",
      text = element_text(size = 12),
      plot.margin = margin(10, 20, 10, 10)
    ) +
    scale_x_continuous(
      breaks = seq(scenario$start_year, max(year_range), by = 5),
      limits = c(scenario$start_year, max(year_range) + 5)
    ) +
    coord_cartesian(ylim = c(60, 100), clip = "off") +
    labs(
      x = 'Year',
      y = 'Share of spending by top 10% of patients (%)',
      title = paste0(
        "Share of spending accounted for by the top 10% of patients\nby type of care, ",
        scenario$start_year,
        " to ",
        max(year_range),
        " — ",
        scenario$title_suffix
      )
    )
  print(p2)
  dev.off()
}
message("\n=== ALL SCENARIOS COMPLETE ===")

# ==========================================================
# OVERLAY PLOT: GROSS VS NET TOP-10% SHARE BY TOC
# ==========================================================
library(ggrepel)

top10_trend_gross[, spending_type := "Gross"]
top10_trend_main[, spending_type := "Net"]
overlay_data <- rbind(top10_trend_gross, top10_trend_main)
overlay_data <- overlay_data[toc != "all" & toc_plot != "Home health"]

# For non-RX tocs, gross and net are identical (same non_rx_col across scenarios),
# so drop the duplicates to avoid plotting them twice
overlay_data <- overlay_data[
  toc == "RX" | spending_type == "Gross"
]
# Relabel non-RX rows so they don't carry a spending_type distinction
overlay_data[toc != "RX", spending_type := "Other"]

# Build a single grouping variable for color/linetype mapping
overlay_data[,
  series := fcase(
    toc == "RX" & spending_type == "Gross" , "RX (Gross)" ,
    toc == "RX" & spending_type == "Net"   , "RX (Net)"   ,
    default = as.character(toc_plot)
  )
]

grey_color <- "grey80"
series_colors <- c(
  "RX (Gross)" = "#00008B",
  "RX (Net)" = "#2CA02C",
  "Ambulatory" = grey_color,
  "Dental" = grey_color,
  "Emergency department" = grey_color,
  "Inpatient" = grey_color
)
series_linetypes <- c(
  "RX (Gross)" = "solid",
  "RX (Net)" = "dashed",
  "Ambulatory" = "solid",
  "Dental" = "solid",
  "Emergency department" = "solid",
  "Inpatient" = "solid"
)
series_sizes <- c(
  "RX (Gross)" = 1,
  "RX (Net)" = 1,
  "Ambulatory" = 0.8,
  "Dental" = 0.8,
  "Emergency department" = 0.8,
  "Inpatient" = 0.8
)

max_yr <- max(overlay_data$year_id)
label_data <- overlay_data[year_id == max_yr]
label_data[,
  label := fcase(
    series == "RX (Gross)" , "Retail prescription\ndrugs (gross)" ,
    series == "RX (Net)"   , "Retail prescription\ndrugs (net)"   ,
    default = as.character(series)
  )
]

pdf(paste0(plot_dir, "single_panel_gross_vs_net.pdf"), width = 8, height = 5)
p_overlay <- ggplot(
  overlay_data,
  aes(
    x = year_id,
    y = top10_share,
    color = series,
    linetype = series,
    group = series
  )
) +
  geom_point(aes(size = series), alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, aes(size = series)) +
  geom_text_repel(
    data = label_data,
    aes(x = max_yr, y = top10_share, label = label, color = series),
    nudge_x = 1.5,
    direction = "y",
    hjust = 0,
    segment.size = 0.4,
    segment.color = "grey80",
    size = 4.5,
    fontface = ifelse(grepl("^RX", label_data$series), "bold", "plain"),
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = series_colors) +
  scale_linetype_manual(values = series_linetypes) +
  scale_size_manual(values = series_sizes, guide = "none") +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 12),
    plot.margin = margin(10, 20, 10, 10)
  ) +
  scale_x_continuous(
    breaks = seq(START_YEAR, END_YEAR, by = 5),
    limits = c(START_YEAR, END_YEAR + 5)
  ) +
  coord_cartesian(ylim = c(60, 100), clip = "off") +
  labs(
    x = "Year",
    y = "Share of spending by top 10% of patients (%)",
    title = paste0(
      "Share of spending accounted for by the top 10% of patients\n",
      "by type of care, ",
      START_YEAR,
      " to ",
      END_YEAR,
      " — gross vs net of rebates"
    )
  )
print(p_overlay)
dev.off()


# ==========================================================
# NUMBER PLUGGING FOR MANUSCRIPT
# ==========================================================

# 1. Total RX observations (Rows loaded from MEPS, full gross analysis window)
gross_scenario <- scenarios[[which(sapply(scenarios, function(s) {
  s$name == "gross_spending"
}))]]
total_rx_obs <- nrow(data_rx_raw[
  year_id >= gross_scenario$start_year & year_id <= max(year_range)
])
total_respondents <- sum(sample_denoms[
  year_id >= min(iso_years) & year_id <= max(iso_years),
  n_obs
])

# 2. Top 10% RX start/end
rx_top10_start <- top10_trend_main[
  toc == "RX" & year_id == min(iso_years),
  top10_share
]
rx_top10_end <- top10_trend_main[
  toc == "RX" & year_id == max(iso_years),
  top10_share
]
rx_top10_diff <- rx_top10_end - rx_top10_start

# 3. Bottom 50% RX
bottom50_trend <- closest_percent_main[,
  .(
    bottom50_share = approx(
      x = percentile,
      y = percent_tot_spend,
      xout = 50,
      rule = 2
    )$y
  ),
  by = .(toc, year_id)
]

rx_bottom50_max <- max(bottom50_trend[toc == "RX", bottom50_share])
rx_bottom50_ceiling <- ceiling(rx_bottom50_max * 10) / 10

# 4. Catastrophic Care ranges (IP and ED)
ip_ed_min <- min(top10_trend_main[toc %in% c("IP", "ED"), top10_share])
ip_ed_max <- max(top10_trend_main[toc %in% c("IP", "ED"), top10_share])

# 5. Dental Care ranges (Stable routine care)
dv_min <- min(top10_trend_main[toc == "DV", top10_share])
dv_max <- max(top10_trend_main[toc == "DV", top10_share])

# 6. Ambulatory Care start/end (To show the modest increase)
am_top10_start <- top10_trend_main[
  toc == "AM" & year_id == min(iso_years),
  top10_share
]
am_top10_end <- top10_trend_main[
  toc == "AM" & year_id == max(iso_years),
  top10_share
]
am_top10_diff <- am_top10_end - am_top10_start

cat(
  "\n========================================================================\n"
)
cat("   MANUSCRIPT NUMBER PLUG PARAGRAPH (RESULTS)\n")
cat(
  "========================================================================\n\n"
)
cat(sprintf(
  "Among %s prescription drug observations in MEPS (%s–%s), we assigned gross-to-net statuses to 97.9%% of identified branded drug expenditures [eTable 1]. Retail pharmaceutical spending became increasingly concentrated over the %d-year period; this pattern persisted after adjusting for rebates and other discounts[Figure 1; eFigure 2]. In %s, the top 10%% of patients accounted for %.1f%% of net retail prescription drug spending; by %s, this decile accounted for %.1f%%, a %.1f-percentage-point increase [Figure 2]. Over the same period, the bottom 50%% of patients accounted for less than %.1f%% of spending. By contrast, concentration remained highly stable for acute inpatient and emergency care (consistently between %.0f-%.0f%%) and dental care (%.0f-%.0f%%). While ambulatory care concentration increased modestly from %.1f%% to %.1f%% (a %.1f-percentage-point increase), prescription drugs were the only sector to exhibit a compounding shift toward a catastrophic financial distribution [Figure 1].\n\n",
  formatC(total_rx_obs, format = "d", big.mark = ","),
  min(iso_years),
  max(iso_years),
  min(iso_years),
  max(iso_years) - min(iso_years),
  rx_top10_start,
  max(iso_years),
  rx_top10_end,
  rx_top10_diff,
  rx_bottom50_ceiling,
  floor(ip_ed_min),
  ceiling(ip_ed_max),
  floor(dv_min),
  ceiling(dv_max),
  am_top10_start,
  am_top10_end,
  am_top10_diff
))
cat(
  "========================================================================\n\n"
)

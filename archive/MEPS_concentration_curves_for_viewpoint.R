# ---------------------------
# Purpose: Examine inequality in RX spending, compared to other tocs
#
# Authors: Hkl1 & msahu
# ---------------------------

source('init.R')

# --------------------------------
#   SET ARGS AND DIRECTORY TO SAVE TO
# --------------------------------

use_wgts <- T # use MEPS weights (instead of counts) to make nationally representative
include_sample_pop <- T # have denominator be anyone in MEPS not just people with at least one encounter (difference between # benes in spend data vs # benes in MEPS are all added on with 0 spend)

save_dir <- paste0(dir, "plots/")
iso_years <- c(2008, 2023) # update. max value dictates max of the timetrend range plotted

# MEPS PATHS
denoms_path <- paste0(data_dir, "raw/meps_denoms/MEPS_sample_denoms.parquet")
meps_data_path <- paste0(data_dir, "raw/meps/")

# Rebate-Adjusted MEPS RX path
adjusted_rx_path <- paste0(
  data_dir,
  "processed/meps_rebate_adjusted/meps_rebate_adjusted.rds"
)

# --------------------------------
#  LOAD SAMPLE DENOMS AND PROCESSED DATA
# --------------------------------
t0 <- Sys.time()

# 1. Sample denoms
sample_denoms <- open_dataset(denoms_path) %>% collect() %>% as.data.table()
sample_denoms <- sample_denoms[,
  .(pop = sum(pop), n_obs = sum(n_obs)),
  by = c('year_id')
]

# 2. Load NON-RX data
all_raw_files <- list.files(
  meps_data_path,
  pattern = "\\.parquet$",
  full.names = T
)
non_rx_files <- all_raw_files[!grepl("USA_MEPS_RX\\.parquet$", all_raw_files)]

data_non_rx <- rbindlist(lapply(non_rx_files, read_parquet), fill = T)

# 3. Load REBATE-ADJUSTED RX data (reading the .rds file)
data_rx <- readRDS(adjusted_rx_path) %>% as.data.table()

# CRITICAL: Overwrite the gross pay with the net pay for the RX data
data_rx[, tot_pay_amt := net_pay_amt]

# 4. Combine all types of care
data <- rbindlist(list(data_non_rx, data_rx), fill = T)

# --------------------------------
# FILTERING & DETRUNC MIMIC
# --------------------------------

# Mimic DETRUNC processing and drop rows with missing info
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

## drop NAs
for (j in drop_cols) {
  data <- data[(!is.na(get(j)) & !get(j) %in% na_vals) | toc == 'DV']
}

# select just one spend amount per claim
data <- data[primary_cause == 1 | dx_level == 'dx_1'] # note, dental doesn't include dx info
data <- data[year_id >= min(iso_years) & year_id <= max(iso_years)]

# Drop encounters with unknown payment; keep encounters with known payment that is 0
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


data[toc %in% c('OBGYN', 'OB', 'OS', 'PAD', 'PC', 'UC'), toc := 'AM']


# --------------------------------
#  CLEAN DATA AND CALCULATE PERCENTILES AND CUMSUM FOR EVERY ROW
# --------------------------------

# Make key in case we want to look within an age etc
data[, survey_wt := as.character(survey_wt)] # just to make sure trailing decimals don't effect this
data[, survey_wt := as.numeric(survey_wt)]

bene_key <- unique(data[, .(
  bene_id,
  age_group_years_start,
  sex_id,
  survey_wt,
  year_id
)])

# Aggregate spend by bene (so across encounters) by toc and all-toc
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

  # add benes not in sample
  if (include_sample_pop == T) {
    sample_w_encounters <- unique(data[, .(tot_wt, toc, year_id)])
    sample_w_encounters <- merge(
      sample_w_encounters,
      sample_denoms,
      by = c('year_id')
    )
    sample_w_encounters[, survey_wt := pop - tot_wt]
    sample_w_encounters[, prop := survey_wt / pop]

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

  # Calculate percentiles for every row
  data[, wt_cumsum := cumsum(survey_wt), by = c('toc', 'year_id')]
  data[, percentile := (wt_cumsum / tot_wt) * 100]
} else {
  data[, tot_n := .N, by = c('toc', 'year_id')]

  if (include_sample_pop == T) {
    sample_w_encounters <- unique(data[, .(tot_n, toc, year_id)])
    sample_w_encounters <- merge(
      sample_w_encounters,
      sample_denoms,
      by = c('year_id')
    )
    sample_w_encounters[, add_obs := n_obs - tot_n]
    sample_w_encounters[add_obs < 1, add_obs := 0]
    sample_w_encounters[, prop := add_obs / n_obs]

    sample_w_encounters <- sample_w_encounters[, .(
      year_id,
      toc,
      N = add_obs,
      pay = 0,
      bene_id = 'no-encounters'
    )]

    data[, N := 1]
    data <- rbind(data, sample_w_encounters, fill = T)
    data <- data[order(pay)]
    data[, tot_n := sum(N), by = c('toc', 'year_id')]
  } else {
    data[, N := 1]
  }

  # Calculate percentiles for every row
  data[, n_cumsum := cumsum(N), by = c('toc', 'year_id')]
  data[, percentile := (n_cumsum / tot_n) * 100]
}

# Sum across each bucket and toc
plot_df <- data[, .(pay = sum(pay)), by = c('percentile', 'toc', 'year_id')]

# add cumulative sum
plot_df <- plot_df[order(percentile)]
plot_df[, pay_cumsum := cumsum(pay), by = c('toc', 'year_id')]

# make cumsum actually % of total spend
plot_df[, tot_spend := max(pay_cumsum), by = c('toc', 'year_id')]
plot_df[, percent_tot_spend := (pay_cumsum / tot_spend) * 100]
print("Data loaded!")
print(Sys.time() - t0)

#-----------
# SIMPLIFY DATA TO JUST GET % spend for each percentile (substantially reduce data size and helps ensure integration is successful)
#-----------
save_plot_df <- plot_df[, .(year_id, toc, percentile, percent_tot_spend)]

# Goal = 100 rows, one for each percent. Find row that is the closest to each percentile, then take unique
closest_percent <- rbindlist(
  lapply(0:100, function(p) {
    #message(p)

    use_data <- copy(save_plot_df)

    if (p == 0) {
      # for p = 0 just find lowest
      use_data[, min_percentile := min(percentile), by = c('year_id', 'toc')]
      return(use_data[percentile == min_percentile, ])
    } else if (p == 100) {
      # for p = 100 just find highest
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

#-----------
# NUMBER PLUG - top 10% of users account for >X% of spending
#-----------
check <- closest_percent[
  toc == 'RX' &
    year_id == max(year_id) &
    percentile > 89.9999 &
    percentile < 90.1
]
paste0(
  "top 10% of users in MEPS account for ",
  round(100 - check$percent_tot_spend, digits = 2),
  "% of total spending"
)

# Extract the Top 10% share robustly using interpolation
# This safely handles highly concentrated care where >90% of people have $0 spend
top10_trend <- closest_percent[,
  .(
    top10_share = 100 -
      approx(x = percentile, y = percent_tot_spend, xout = 90, rule = 2)$y
  ),
  by = .(toc, year_id)
]

# Map the nice names for the plot facets
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

#-----------
# CALCULATE GINI COEFFICIENT
#-----------

gini_plots <- F # Option to make plot, keep F
plot_years <- iso_years

if (gini_plots == T) {
  pdf(
    paste0(
      save_dir,
      "/meps_concentration_curves_",
      ifelse(use_wgts, "", "notweighted_"),
      ifelse(include_sample_pop, "", "peoplewithenc"),
      paste0(plot_years, collapse = "_"),
      ".pdf"
    ),
    width = 11,
    height = 8
  )
}

# use simplified data
plot_df <- closest_percent
toc_gini_timetrend <- data.table()

for (t in sort(unique(plot_df$toc))) {
  message(t)

  years <- unique(plot_df[toc == t]$year_id)
  gini_by_year <- rbindlist(lapply(years, function(y) {
    #message(y)
    tmp_plot_data <- plot_df[toc == t & year_id == y]

    #tmp_plot_data <- rbind(tmp_plot_data, data.table(percentile = c(0,1), toc = t, year_id = y, pay = c(0,1), pay_cumsum = c(0,1), tot_spend = NA, percent_tot_spend = c(0,1)))
    tmp_plot_data <- rbind(
      tmp_plot_data,
      data.table(
        percentile = c(0, 100),
        toc = t,
        year_id = y,
        percent_tot_spend = c(0, 100)
      ),
      fill = T
    )

    # calculate gini:      A / (A+B)
    # C = area beneath line of equality
    # B = area beneath plotted line (percent total spend)
    # A = C-B (area between line of equality and plotted line)

    # use approx fun and integrate to get area under new curve
    line_fun <- approxfun(
      tmp_plot_data$percentile,
      tmp_plot_data$percent_tot_spend,
      yleft = 0,
      yright = 0
    )
    # tmp_plot_data$pred_percent <- line_fun(tmp_plot_data$percentile)

    get_gini <- tryCatch(
      expr = {
        integrate(
          line_fun,
          min(tmp_plot_data$percentile),
          max(tmp_plot_data$percentile)
        )$value
      },
      error = function(e) {
        return(NA)
      }
    )

    if (!is.na(get_gini)) {
      b <- integrate(
        line_fun,
        min(tmp_plot_data$percentile),
        max(tmp_plot_data$percentile)
      )$value
      c <- (max(tmp_plot_data$percentile) *
        max(tmp_plot_data$percent_tot_spend)) /
        2
      a <- c - b
      gini_coef_approx <- a / (a + b)
    } else {
      gini_coef_approx <- NA
    }

    return(data.table(
      'year_id' = y,
      'gini' = round(gini_coef_approx, digits = 2)
    ))
  }))

  toc_gini_timetrend <- rbind(toc_gini_timetrend, gini_by_year[, toc := t])
  #
  tmp_plot_data <- plot_df[toc == t]

  if (include_sample_pop == T) {
    x_axis_lab <- 'Percent of beneficaries in MEPS'
  } else {
    if (t == 'all') {
      x_axis_lab <- "Percent of people with a health system encounter"
    } else {
      x_axis_lab <- paste0(
        'Percent of people with at least one ',
        t,
        " encounter"
      )
    }
  }

  plot <- ggplot(
    tmp_plot_data[year_id %in% plot_years],
    aes(x = percentile, y = percent_tot_spend)
  ) +
    geom_line() +
    facet_wrap(~year_id) +
    geom_ribbon(
      aes(ymin = percent_tot_spend, ymax = percentile),
      fill = "grey",
      alpha = 0.5
    ) +
    geom_text(
      data = gini_by_year[year_id %in% plot_years],
      aes(x = 3, y = 15, label = gini),
      size = 3.5
    ) +
    labs(
      x = x_axis_lab,
      y = 'Percent of spending',
      title = paste0(
        'Spending equality ',
        ifelse(t == 'all', "across all TOCs", paste0("in ", t))
      )
    ) +
    #subtitle = paste0("Gini ", round(gini_coef_approx, digits = 4)))+
    geom_abline(aes(intercept = 0, slope = 1), linetype = 3) +
    #scale_y_continuous(sec.axis = sec_axis(~ . * tot_val/100, name = "Cumulative spending ($)", labels = scales::comma))+
    theme_bw()

  if (gini_plots == T) {
    print(plot)
  }
}

if (gini_plots == T) {
  dev.off()
}


#-----------
# FIGURE - overlap different years of RX gini time trends
#-----------

# Make sure there is a starting value everywhere
plot_lorenz_curve_data <- plot_df[toc == 'RX']
add_start_val <- unique(plot_lorenz_curve_data[, .(year_id, toc)])
add_start_val[, `:=`(percent_tot_spend = 0, percentile = 0)]
plot_lorenz_curve_data <- rbind(plot_lorenz_curve_data, add_start_val)

# 1. PREPARE THE "SHADED GAP" DATA (Interpolated for perfect ribbon matching)
# Create a sequence of 1000 points from 0 to 100
ribbon_x <- seq(0, 100, length.out = 1000)

# Interpolate the Y values for the early year (e.g. 2008)
ribbon_y_early <- approx(
  x = plot_lorenz_curve_data[year_id == min(iso_years)]$percentile,
  y = plot_lorenz_curve_data[year_id == min(iso_years)]$percent_tot_spend,
  xout = ribbon_x,
  rule = 2 # rule=2 ensures no NAs at the extreme edges
)$y

# Interpolate the Y values for the late year (e.g. 2023)
ribbon_y_late <- approx(
  x = plot_lorenz_curve_data[year_id == max(iso_years)]$percentile,
  y = plot_lorenz_curve_data[year_id == max(iso_years)]$percent_tot_spend,
  xout = ribbon_x,
  rule = 2
)$y

# Combine into a clean, NA-free dataframe for the ribbon
ribbon_data <- data.table(
  percentile = ribbon_x,
  y_early = ribbon_y_early,
  y_late = ribbon_y_late
)

# 2. CALCULATE EXACT 90TH PERCENTILE INTERSECTIONS
point_90 <- plot_lorenz_curve_data[
  year_id %in% iso_years,
  .(
    y_val = approx(x = percentile, y = percent_tot_spend, xout = 90, rule = 2)$y
  ),
  by = year_id
]

point_90[, top10_share := paste0(sprintf("%.1f", 100 - y_val), "%")]
point_90[, combined_label := paste0(year_id, ": ", top10_share)]


# 3. GENERATE THE PLOT
pdf(paste0(save_dir, "/lorenz_curves_overlay_gap.pdf"), width = 7, height = 5)

plot <- ggplot() +
  geom_abline(intercept = 0, slope = 1, linetype = 3, color = "grey50") +

  # THE SHADED GAP (Now using our perfectly interpolated NA-free dataset)
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
      "The top 10% of highest-cost patients accounted for\n",
      point_90[year_id == min(iso_years), top10_share],
      " of net spending in ",
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

  coord_cartesian(xlim = c(0, 100), ylim = c(0, 100), clip = "off") +

  labs(
    x = 'Cumulative percentage of US population',
    y = "Cumulative percentage of net retail prescription drug spending" #,
    #    title = "Figure 1: Increasing Concentration of Retail Prescription Drug Spending"
  )

print(plot)
dev.off()

#-----------
# FIGURE - time trends of gini coefficients for different types of care
#-----------

pdf_name <- paste0(
  save_dir,
  "/gini_coefficient_timetrend_",
  ifelse(use_wgts, "", "notweighted_"),
  ifelse(include_sample_pop, "", "peoplewithenc"),
  max(iso_years)
)


toc_gini_timetrend[,
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


base_colors2 <- brewer.pal(6, "Pastel2") # Or "Dark2", "Paired", etc.
base_colors2[5] <- "#D1E5F4" # Replace one with a pop color
base_colors2[6] <- "#00008C" # Replace one with a pop color


base_colors2 <- c(
  "#B3E2CD",
  "#FDCDAC",
  "#ACACE6",
  "#F4CAE4",
  "#AECCE4",
  "#00008C"
)

pdf(paste0(pdf_name, ".pdf"), width = 6, height = 4.5)
plot <- ggplot(
  toc_gini_timetrend[toc != 'all'],
  aes(x = year_id, y = gini, color = toc_plot, fill = toc_plot)
) +
  facet_wrap(~toc_plot, nrow = 2) +
  ylim(0.7, 1) +
  geom_point() +
  geom_smooth(method = 'lm', se = F) +
  theme_bw() +
  #theme(legend.position = 'bottom')+
  theme(legend.position = 'none') +
  scale_color_manual(values = base_colors2) +
  scale_fill_manual(values = base_colors2) +
  labs(
    x = 'Year',
    y = 'Gini coefficient',
    title = paste0(
      'Inequality index for health care spending over time \n by type of care, 2000 to 2022'
    ),
    color = '',
    fill = ''
  )
print(plot)
dev.off()

write.csv(
  toc_gini_timetrend,
  paste0(save_dir, "/gini_coef_timetrend_data.csv"),
  row.names = F
)

#-----------
# FIGURE - presented as top 10% share
#-----------

pdf(paste0(pdf_name, "_top10.pdf"), width = 6, height = 4.5)
plot <- ggplot(
  top10_trend[toc != 'all'],
  aes(x = year_id, y = top10_share, color = toc_plot, fill = toc_plot)
) +
  facet_wrap(~toc_plot, nrow = 2) +
  ylim(50, 100) + # Sets Y-axis from 50% to 100%
  geom_point() +
  geom_smooth(method = 'lm', se = F) +
  theme_bw() +
  theme(legend.position = 'none') +
  scale_color_manual(values = base_colors2) +
  scale_fill_manual(values = base_colors2) +
  labs(
    x = 'Year',
    y = 'Share of spending by top 10% of patients (%)',
    title = 'Panel B: Spending by the top 10% of patients over time \n by type of care, 2008 to 2022',
    color = '',
    fill = ''
  )
print(plot)
dev.off()

#-----------
# SINGLE PANEL FIGURE - presented as top 10% share
#-----------

library(ggrepel)

# 1. Prep data (Do not change the strings in the data table!)
single_plot_data <- single_plot_data[toc_plot != "Home health"]
max_yr <- max(single_plot_data$year_id)
label_data <- single_plot_data[year_id == max_yr]

# 2. Plot
pdf(paste0(pdf_name, "_single_panel_clean.pdf"), width = 8, height = 5)

plot <- ggplot(
  single_plot_data,
  aes(x = year_id, y = top10_share, color = toc_plot, group = toc_plot)
) +
  geom_point(aes(alpha = is_rx, size = is_rx)) +
  geom_smooth(method = 'lm', se = FALSE, aes(size = is_rx)) +

  # Add the labels
  geom_text_repel(
    data = label_data,
    # FIX: Change the label text directly in the aes() mapping!
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
    # Keep this looking for the original string name
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
    breaks = c(2008, 2013, 2018, 2023),
    limits = c(2008, 2028)
  ) +

  # FIX: Put clip = "off" back in so nothing ever gets chopped!
  coord_cartesian(ylim = c(60, 100), clip = "off") +

  labs(
    x = 'Year',
    y = 'Share of spending by top 10% of patients (%)',
    title = 'Figure 1: Share of spending accounted for by the top 10% of patients\nby type of care, 2008 to 2023'
  )

print(plot)
dev.off()

# ---------------------------
# Purpose: Examine inequality in RX spending, compared to other tocs
#
# Authors: Hkl1/ msahu
# ---------------------------

source('init.R')

# --------------------------------
#   SET ARGS AND DIRECTORY TO SAVE TO
# --------------------------------

use_wgts <- T # use MEPS weights (instead of counts) to make nationally representative
include_sample_pop <- T # have denominator be anyone in MEPS not just people with at least one encounter (difference between # benes in spend data vs # benes in MEPS are all added on with 0 spend)

save_dir <- paste0(plot_dir, 'concentration_curves')
ensure_dir(save_dir)
iso_years <- c(2000, 2022) # update. max value dictates max of the timetrend range plotted

# MEPS PATHS
denoms_path <- paste0(data_dir, 'meps_denoms')
meps_data_path <- paste0(data_dir, 'meps')

# --------------------------------
#  LOAD SAMPLE DENOMS AND PROCESSED DATA - takes 2.5 minutes to load and process all data
# --------------------------------
t0 <- Sys.time()

# Sample denoms
sample_denoms <- open_dataset(denoms_path) %>% collect() %>% as.data.table()
sample_denoms <- sample_denoms[,
  .(pop = sum(pop), n_obs = sum(n_obs)),
  by = c('year_id')
]

# Load data from data processing
data <- rbindlist(
  lapply(list.files(meps_data_path, full.names = T), read_parquet),
  fill = T
)

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
data <- data[primary_cause == 1 | dx_level == 'dx_1']
data <- data[year_id >= 2000 & year_id <= max(iso_years)]

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


data[toc %in% c('OBGYN', 'OS', 'PAD', 'PC', 'UC'), toc := 'AM']


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
#bene_key[,.N, by = c('bene_id', 'year_id')][N > 1]

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

# Make sure there is a starting value everywhere for the right shading
plot_lorenz_curve_data <- plot_df[toc == 'RX']
add_start_val <- unique(plot_lorenz_curve_data[, .(year_id, toc)])
add_start_val[, `:=`(percent_tot_spend = 0, percentile = 0)]
plot_lorenz_curve_data <- rbind(plot_lorenz_curve_data, add_start_val)

pdf(
  paste0(
    save_dir,
    "/lorenz_curves_overlay",
    ifelse(use_wgts, "", "_notweighted_"),
    ifelse(include_sample_pop, "", "_peoplewithenc"),
    max(iso_years),
    ".pdf"
  ),
  width = 6,
  height = 4.5
) # 8,6


# calculate the 90%
check <- plot_lorenz_curve_data[
  toc == 'RX' &
    percentile > 89.9999 &
    percentile < 90.1 &
    year_id %in% iso_years
]
percent_key <- check[,
  ten_percent_user_spend := paste0(
    round(100 - percent_tot_spend, digits = 0),
    "%"
  )
]
percent_key[, percentile := percentile + 3.5]
percent_key[, percent_tot_spend_plot := percent_tot_spend - 2]

# Make plot
plot <- ggplot(
  plot_lorenz_curve_data[year_id %in% iso_years],
  aes(x = percentile, y = percent_tot_spend, group = year_id)
) +
  #geom_hline(yintercept = 90, color = 'darkgrey')+
  geom_vline(xintercept = 90, color = 'darkgrey') +
  geom_line(aes(color = factor(year_id)), size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = 3) +
  #geom_hline(data = percent_key, aes(yintercept = percent_tot_spend, color = factor(year_id)), linetype = 2)+
  geom_segment(
    data = percent_key,
    aes(
      y = percent_tot_spend,
      x = 0,
      xend = 90,
      yend = percent_tot_spend,
      color = factor(year_id)
    ),
    linetype = 2
  ) +
  geom_text(
    data = percent_key,
    aes(
      label = ten_percent_user_spend,
      color = factor(year_id),
      y = percent_tot_spend_plot
    ),
    show.legend = FALSE
  ) +
  scale_color_brewer(palette = 'Set1') +
  scale_fill_brewer(palette = 'Set1') +
  geom_ribbon(
    aes(ymin = percent_tot_spend, ymax = percentile, fill = factor(year_id)),
    alpha = 0.3
  ) +
  labs(
    x = 'Percent of respondents in the Medical Expenditure Panel Survey [MEPS]',
    y = "Percent of retail prescription drug spending",
    title = paste0(
      "Panel A: Concentration of retail prescription drug \n spending, 2000 versus 2022"
    ),
    color = "Year",
    fill = "Year"
  ) +
  theme_bw()
print(plot)

dev.off()

write.csv(
  plot_lorenz_curve_data,
  paste0(save_dir, "/lorenz_curves_overlay_data.csv"),
  row.names = F
)


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
  ylim(0.5, 1) +
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
      'Panel B: Inequality index for health care spending over time \n by type of care, 2000 to 2022'
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

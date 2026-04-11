
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(readr)
  library(stringr)
  library(scales)
  library(ggplot2)
  library(mc2d) # for rpert()
})

options(scipen = 999)

# ============================================================
# WGS cost-effectiveness model for Hospital A / TARSS
# Price year: 2022 USD
# Includes:
#   - Base-case analysis
#   - One-way sensitivity analysis (OWSA)
#   - Probabilistic sensitivity analysis (PSA)
#   - Cost-effectiveness acceptability curve (CEAC)
#   - Budget impact analysis (BIA)
# ============================================================

# ------------------------------------------------------------
# 1) Input data: antibiotic cost context
# ------------------------------------------------------------
antibiotic_context <- tribble(
  ~year, ~mdr_pct, ~kpn_infected_uti_bsi, ~kpn_infected_all_spec, ~abx_spend_usd,
  2018,   0.65,     154,                    288,                    39469.51,
  2019,   0.58,     120,                    219,                    41951.74,
  2020,   0.62,      73,                    139,                    38119.77,
  2021,   0.53,     119,                    180,                    32433.68,
  2022,   0.52,      78,                    126,                    24547.08
)

get_antibiotic_bounds <- function(data, denominator = c("UTI_BSI", "ALL_SPECIMENS")) {
  denominator <- match.arg(denominator)
  
  data %>%
    mutate(
      denom_infections = if_else(
        denominator == "UTI_BSI",
        kpn_infected_uti_bsi,
        kpn_infected_all_spec
      ),
      est_3gcr_infections = denom_infections * mdr_pct,
      abx_cost_per_inf = abx_spend_usd / pmax(est_3gcr_infections, 1e-9)
    ) %>%
    summarise(
      abx_low = min(abx_cost_per_inf, na.rm = TRUE),
      abx_high = max(abx_cost_per_inf, na.rm = TRUE)
    )
}

abx_bounds <- get_antibiotic_bounds(antibiotic_context, denominator = "UTI_BSI")
abx_base <- 605.2

# ------------------------------------------------------------
# 2) Core parameters
# ------------------------------------------------------------
TAT_MIN <- 4
TAT_BASE <- 7
TAT_MAX <- 10

base_params <- list(
  # WGS micro-costing (pre-overhead per test)
  wgs_equip_preoh = 88.46,
  wgs_cons_preoh = 131.96,
  wgs_staff_preoh = 12.00,
  overhead_rate = 0.10,
  batch_eff = 1.00,
  
  # Epidemiology and operations
  clusters_4yrs = 11,
  cases_per_cluster = 4,
  tat_days = TAT_BASE,
  ipc_effect = 0.50,
  eligible_share = 0.60,
  pi_c = 1.00,
  
  # Savings per infection
  abx_cost_usd = abx_base,
  closed_bed_day_cost_usd = 246,
  excess_los_days = 5,
  cleaning_nurse_per_det = 122,
  ppe_cost_per_iso_day = 50,
  iso_days_per_prev_inf = 5,
  
  # Rule-out savings
  sro_override = 0,
  n_contacts = 8,
  d_isolation_avoided = 1,
  screening_rounds = 1,
  cscreen_unit_cost = 28.6
)

# ------------------------------------------------------------
# 3) Helper functions
# ------------------------------------------------------------
format_money <- function(x) {
  ifelse(is.finite(x), dollar(round(x, 0)), "∞")
}

format_money_cost_saving <- function(x) {
  out <- format_money(x)
  ifelse(is.finite(x) & x < 0, paste0(out, " (CS)"), out)
}

format_pct_change <- function(delta, base_value) {
  pct <- ifelse(
    is.finite(delta) & is.finite(base_value) & base_value != 0,
    100 * delta / base_value,
    NA_real_
  )
  ifelse(is.na(pct), "—", paste0(ifelse(pct < 0, "-", ""), round(abs(pct), 0), "%"))
}

tat_modifier <- function(tat_days, min_day = TAT_MIN, max_day = TAT_MAX, floor = 0.1) {
  x <- (max_day - tat_days) / (max_day - min_day)
  pmin(pmax(x, floor), 1)
}

wgs_cost_per_test <- function(params) {
  pre_oh <- params$wgs_equip_preoh +
    params$wgs_staff_preoh +
    (params$wgs_cons_preoh / params$batch_eff)
  
  pre_oh * (1 + params$overhead_rate)
}

ruleout_saving_per_investigation <- function(params) {
  if (!is.null(params$sro_override) && is.finite(params$sro_override)) {
    return(params$sro_override)
  }
  
  isolation_component <- params$n_contacts *
    params$d_isolation_avoided *
    params$ppe_cost_per_iso_day
  
  cleaning_component <- params$cleaning_nurse_per_det
  
  screening_component <- params$n_contacts *
    params$screening_rounds *
    params$cscreen_unit_cost
  
  isolation_component + cleaning_component + screening_component
}

savings_per_infection <- function(params) {
  params$abx_cost_usd +
    (params$closed_bed_day_cost_usd * params$excess_los_days) +
    params$cleaning_nurse_per_det +
    (params$ppe_cost_per_iso_day * params$iso_days_per_prev_inf)
}

sro_upper <- ruleout_saving_per_investigation(modifyList(base_params, list(sro_override = NA_real_)))

# ------------------------------------------------------------
# 4) Deterministic model
# ------------------------------------------------------------
run_model_full <- function(params) {
  per_test <- wgs_cost_per_test(params)
  tests <- params$clusters_4yrs * params$cases_per_cluster
  program_cost <- per_test * tests
  
  downstream_cases <- max(params$cases_per_cluster - 2, 0)
  infections_averted <- params$pi_c *
    params$clusters_4yrs *
    downstream_cases *
    params$ipc_effect *
    params$eligible_share *
    tat_modifier(params$tat_days)
  
  avoided_costs <- infections_averted * savings_per_infection(params)
  ruleout_savings <- (1 - params$pi_c) *
    params$clusters_4yrs *
    ruleout_saving_per_investigation(params)
  
  total_savings <- avoided_costs + ruleout_savings
  net_cost <- program_cost - total_savings
  beddays_avoided <- infections_averted * params$excess_los_days
  
  tibble(
    per_test = per_test,
    tests = tests,
    program_cost = program_cost,
    infections_averted = infections_averted,
    avoided_costs = avoided_costs,
    ruleout_savings = ruleout_savings,
    total_savings = total_savings,
    net_cost = net_cost,
    cpia = net_cost / pmax(infections_averted, 1e-9),
    beddays_avoided = beddays_avoided,
    cost_per_bedday = net_cost / pmax(beddays_avoided, 1e-9)
  )
}

run_model_cpia <- function(params) {
  run_model_full(params) %>% pull(cpia)
}

make_base_summary_row <- function(params, use_total_savings = FALSE) {
  out <- run_model_full(params)
  clusters <- params$clusters_4yrs
  tests <- clusters * params$cases_per_cluster
  avoided <- if (isTRUE(use_total_savings)) out$total_savings else out$avoided_costs
  
  tibble(
    `Number of clusters [number of tests]` = sprintf("%d [%d]", clusters, tests),
    `WGS program overall cost (USD)` = round(out$program_cost, 0),
    `Avoided costs` = round(avoided, 0),
    `Net cost [WGS – avoided]` = round(out$program_cost - avoided, 0),
    `Bed-days avoided (total)` = round(out$beddays_avoided, 1),
    `Cost per bed-day avoided` = round(out$cost_per_bedday, 0),
    `Number of infections averted (total)` = round(out$infections_averted, 1),
    `Cost per infection averted` = round(out$cpia, 0)
  )
}

# ------------------------------------------------------------
# 5) Base case
# ------------------------------------------------------------
base_output <- run_model_full(base_params)
base_summary <- make_base_summary_row(base_params)
write_csv(base_summary, "outputs/base_case_summary.csv")

# ------------------------------------------------------------
# 6) One-way sensitivity analysis (OWSA)
# ------------------------------------------------------------
param_specs <- tribble(
  ~param,                    ~low,                 ~base_val,                        ~high,              ~units,
  "overhead_rate",            0.05,                 base_params$overhead_rate,        0.30,              "proportion",
  "ipc_effect",               0.30,                 base_params$ipc_effect,           0.70,              "proportion",
  "clusters_4yrs",            8,                    base_params$clusters_4yrs,        20,                "count",
  "cases_per_cluster",        3,                    base_params$cases_per_cluster,    6,                 "count",
  "tat_days",                 TAT_MIN,              base_params$tat_days,             TAT_MAX,           "days",
  "eligible_share",           0.40,                 base_params$eligible_share,       0.80,              "proportion",
  "closed_bed_day_cost_usd",  151,                  base_params$closed_bed_day_cost_usd, 342,           "USD/day",
  "excess_los_days",          3,                    base_params$excess_los_days,      7,                 "days",
  "abx_cost_usd",             abx_bounds$abx_low,   base_params$abx_cost_usd,         abx_bounds$abx_high, "USD/infection",
  "cleaning_nurse_per_det",   90,                   base_params$cleaning_nurse_per_det, 155,            "USD/detection",
  "ppe_cost_per_iso_day",     35,                   base_params$ppe_cost_per_iso_day, 65,               "USD/day",
  "iso_days_per_prev_inf",    3,                    base_params$iso_days_per_prev_inf, 7,               "days",
  "pi_c",                     0.30,                 base_params$pi_c,                 0.70,              "proportion",
  "sro_override",             0.0,                  base_params$sro_override,         sro_upper,         "USD/investigation",
  "batch_eff",                0.50,                 base_params$batch_eff,            1.00,              "proportion"
)

base_cpia <- run_model_cpia(base_params)

owsa_raw <- param_specs %>%
  rowwise() %>%
  mutate(
    cost_low = {
      params <- base_params
      params[[param]] <- low
      run_model_cpia(params)
    },
    cost_high = {
      params <- base_params
      params[[param]] <- high
      run_model_cpia(params)
    }
  ) %>%
  ungroup() %>%
  mutate(
    cost_base = base_cpia,
    delta_low = cost_low - cost_base,
    delta_high = cost_high - cost_base
  )

owsa_table <- owsa_raw %>%
  transmute(
    Parameter = str_replace_all(param, "_", " "),
    `Base (unit)` = paste0(base_val, " ", units),
    `Lower bound` = low,
    `Higher bound` = high,
    `Cost at lower` = format_money_cost_saving(cost_low),
    `Cost at upper` = format_money_cost_saving(cost_high),
    `Delta vs base (lower, upper)` = paste0(
      format_pct_change(delta_low, cost_base), ", ",
      format_pct_change(delta_high, cost_base)
    )
  )

write_csv(owsa_raw, "outputs/owsa_raw.csv")
write_csv(owsa_table, "outputs/owsa_formatted.csv")

# ------------------------------------------------------------
# 7) Probabilistic sensitivity analysis (PSA)
# ------------------------------------------------------------
wilson_ci <- function(k, n, z = qnorm(0.975)) {
  p <- k / n
  denom <- 1 + (z^2 / n)
  center <- (p + (z^2 / (2 * n))) / denom
  half <- (z * sqrt(p * (1 - p) / n + (z^2 / (4 * n^2)))) / denom
  c(lower = max(0, center - half), upper = min(1, center + half))
}

sample_pert <- function(low, mode, high, units) {
  vals <- sort(c(low, mode, high))
  draw <- mc2d::rpert(1, min = vals[1], mode = vals[2], max = vals[3])
  
  units <- tolower(trimws(units))
  if (grepl("count", units)) draw <- round(pmax(draw, 0))
  if (grepl("proportion", units)) draw <- min(max(draw, 0), 1)
  
  draw
}

draw_psa_once <- function() {
  params <- base_params
  
  for (i in seq_len(nrow(param_specs))) {
    params[[param_specs$param[i]]] <- sample_pert(
      low = param_specs$low[i],
      mode = param_specs$base_val[i],
      high = param_specs$high[i],
      units = param_specs$units[i]
    )
  }
  
  out <- run_model_full(params)
  num_cols <- vapply(out, is.numeric, logical(1))
  out[num_cols] <- lapply(out[num_cols], function(x) {
    x[!is.finite(x)] <- NA_real_
    x
  })
  
  out
}

set.seed(123)
n_sim <- 1000
psa_draws <- map_dfr(seq_len(n_sim), ~ draw_psa_once())

psa_summary <- psa_draws %>%
  summarise(
    cpia_mean = mean(cpia, na.rm = TRUE),
    cpia_median = median(cpia, na.rm = TRUE),
    cpia_q025 = quantile(cpia, 0.025, na.rm = TRUE),
    cpia_q975 = quantile(cpia, 0.975, na.rm = TRUE),
    net_mean = mean(net_cost, na.rm = TRUE),
    net_median = median(net_cost, na.rm = TRUE),
    net_q025 = quantile(net_cost, 0.025, na.rm = TRUE),
    net_q975 = quantile(net_cost, 0.975, na.rm = TRUE),
    inf_mean = mean(infections_averted, na.rm = TRUE),
    inf_median = median(infections_averted, na.rm = TRUE),
    p_cost_saving = mean(net_cost < 0, na.rm = TRUE),
    p_cpia_negative = mean(cpia < 0, na.rm = TRUE)
  )

write_csv(psa_draws, "outputs/psa_draws_1000.csv")
write_csv(psa_summary, "outputs/psa_summary_1000.csv")

# ------------------------------------------------------------
# 8) CEAC
# ------------------------------------------------------------
wtp_grid <- seq(0, 10000, by = 100)

calc_ce_probability <- function(wtp, draws) {
  ce <- (wtp * draws$infections_averted - draws$net_cost) >= 0
  k <- sum(ce, na.rm = TRUE)
  p <- k / nrow(draws)
  ci <- wilson_ci(k, nrow(draws))
  
  tibble(
    WTP = wtp,
    p = p,
    p_low = ci["lower"],
    p_high = ci["upper"]
  )
}

ceac <- map_dfr(wtp_grid, calc_ce_probability, draws = psa_draws)
write_csv(ceac, "outputs/ceac_psa_1000.csv")

wtp_lines <- tribble(
  ~name,                            ~WTP,
  "0.5x GDP per capita (US$2,175)", 2175,
  "1x GDP per capita (US$4,350)",   4350,
  "2x GDP per capita (US$8,700)",   8700
)

p_ceac <- ggplot(ceac, aes(x = WTP, y = p)) +
  geom_ribbon(aes(ymin = p_low, ymax = p_high), fill = "grey80", alpha = 0.4) +
  geom_line(linewidth = 1) +
  geom_vline(data = wtp_lines, aes(xintercept = WTP, linetype = name), alpha = 0.9) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_x_continuous(labels = dollar_format(prefix = "$")) +
  labs(
    x = "Willingness-to-pay per infection averted (USD)",
    y = "Probability programme is cost-effective",
    linetype = "WTP thresholds"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

ggsave("outputs/ceac_psa_1000.png", p_ceac, width = 7.5, height = 4.5, dpi = 300)

# ------------------------------------------------------------
# 9) Budget impact analysis (BIA)
# ------------------------------------------------------------
bia_params <- list(
  wgs_cost_pre_oh = 232.42,
  overhead_rate = 0.10,
  cases_per_cluster = 4,
  ipc_effect = 0.50,
  eligible_share = 0.60,
  tat_base_days = 7,
  tat_min = 4,
  tat_max = 10,
  tat_floor = 0.10,
  abx_cost_usd = 605.2,
  closed_bed_day_cost_usd = 246,
  excess_los_days = 5,
  cleaning_nurse_per_det = 122,
  ppe_cost_per_iso_day = 50,
  iso_days_per_prev_inf = 5
)

bia_wgs_cost_per_test <- function(params) {
  params$wgs_cost_pre_oh * (1 + params$overhead_rate)
}

bia_tat_modifier <- function(tat_days, params) {
  x <- (params$tat_max - tat_days) / (params$tat_max - params$tat_min)
  pmin(pmax(x, params$tat_floor), 1)
}

bia_savings_per_infection <- function(params) {
  params$abx_cost_usd +
    (params$closed_bed_day_cost_usd * params$excess_los_days) +
    params$cleaning_nurse_per_det +
    (params$ppe_cost_per_iso_day * params$iso_days_per_prev_inf)
}

site_from_clusters <- function(clusters, tat_days, params) {
  tests <- clusters * params$cases_per_cluster
  program_cost <- tests * bia_wgs_cost_per_test(params)
  downstream_cases <- max(params$cases_per_cluster - 2, 0)
  infections_averted <- clusters * downstream_cases * params$ipc_effect * params$eligible_share * bia_tat_modifier(tat_days, params)
  avoided <- infections_averted * bia_savings_per_infection(params)
  net <- program_cost - avoided
  beddays <- infections_averted * params$excess_los_days
  
  tibble(
    clusters = clusters,
    tests = tests,
    program_cost = program_cost,
    avoided = avoided,
    net = net,
    infections_averted = infections_averted,
    beddays = beddays
  )
}

allocate_integer_total <- function(total, weights) {
  raw <- total * weights
  base <- floor(raw)
  remainder <- total - sum(base)
  
  if (remainder > 0) {
    idx <- order(raw - base, decreasing = TRUE)[seq_len(remainder)]
    base[idx] <- base[idx] + 1
  }
  
  as.integer(base)
}

clusters_vec <- c(4, 8, 12, 16, 20, 24)

hospital_table <- map_dfr(clusters_vec, function(clusters) {
  out <- site_from_clusters(clusters, tat_days = bia_params$tat_base_days, params = bia_params)
  
  tibble(
    Perspective = "I. Hospital perspective",
    clusters = clusters,
    tests = out$tests,
    `Number of clusters [number of tests]` = paste0(clusters, " [", out$tests, "]"),
    `Total program costs` = round(out$program_cost, 0),
    `Cost savings` = round(out$avoided, 0),
    `Net costs` = round(out$net, 0),
    `Number of infections averted` = round(out$infections_averted, 2),
    `Bed-days avoided (total)` = round(out$beddays, 1),
    `Cost per infection averted` = round(out$net / out$infections_averted, 0),
    `Cost per bed-day avoided` = round(out$net / out$beddays, 0)
  )
})

nrl_table <- map_dfr(clusters_vec, function(clusters) {
  tests <- clusters * bia_params$cases_per_cluster
  program_cost <- tests * bia_wgs_cost_per_test(bia_params)
  
  tibble(
    Perspective = "II. NRL perspective",
    clusters = clusters,
    tests = tests,
    `Number of clusters [number of tests]` = paste0(clusters, " [", tests, "]"),
    `Total program costs` = round(program_cost, 0),
    `Cost savings` = NA_real_,
    `Net costs` = round(program_cost, 0),
    `Number of infections averted` = NA_real_,
    `Bed-days avoided (total)` = NA_real_,
    `Cost per infection averted` = NA_real_,
    `Cost per bed-day avoided` = NA_real_
  )
})

sites_11 <- tribble(
  ~site,                                ~tat_add, ~weight,
  " (Ariana)",          1,     1/11,
  "(Ben Arous)",      1,     1/11,
  "(Manouba)",      1,     1/11,
  "(Monastir)",       2,     1/11,
  "(Sfax)",              2,     1/11,
  "(Sousse)",              2,     1/11,
  "(Tunis)",                0,     1/11,
  "(Tunis)",                       0,     1/11,
  "(Tunis)",                          0,     1/11,
  "(Tunis)",                           0,     1/11,
  "(Tunis)",               0,     1/11
)

system_table <- map_dfr(clusters_vec, function(clusters_per_hospital) {
  total_clusters <- clusters_per_hospital * 11
  clusters_by_site <- allocate_integer_total(total_clusters, sites_11$weight)
  total_tests <- total_clusters * bia_params$cases_per_cluster
  total_program_cost <- total_tests * bia_wgs_cost_per_test(bia_params)
  
  infections_averted_total <- map2_dbl(clusters_by_site, sites_11$tat_add, function(site_clusters, tat_add) {
    site_from_clusters(
      clusters = site_clusters,
      tat_days = bia_params$tat_base_days + tat_add,
      params = bia_params
    ) %>%
      pull(infections_averted)
  }) %>%
    sum()
  
  avoided_total <- infections_averted_total * bia_savings_per_infection(bia_params)
  net_total <- total_program_cost - avoided_total
  beddays_total <- infections_averted_total * bia_params$excess_los_days
  
  tibble(
    Perspective = "III. NRL + 11-hospital network (system)",
    clusters = total_clusters,
    tests = total_tests,
    `Number of clusters [number of tests]` = paste0(total_clusters, " [", total_tests, "]"),
    `Total program costs` = round(total_program_cost, 0),
    `Cost savings` = round(avoided_total, 0),
    `Net costs` = round(net_total, 0),
    `Number of infections averted` = round(infections_averted_total, 2),
    `Bed-days avoided (total)` = round(beddays_total, 1),
    `Cost per infection averted` = round(net_total / infections_averted_total, 0),
    `Cost per bed-day avoided` = round(net_total / beddays_total, 0)
  )
})

bia_raw <- bind_rows(hospital_table, nrl_table, system_table)
write_csv(bia_raw, "outputs/bia_all_perspectives_raw.csv")

# ------------------------------------------------------------
# 10) Plots
# ------------------------------------------------------------
ggsave(
  filename = "outputs/psa_cpia_hist.png",
  plot = ggplot(psa_draws, aes(cpia)) +
    geom_histogram(bins = 40) +
    labs(x = "USD per infection averted", y = "Count") +
    theme_minimal(),
  width = 7,
  height = 4,
  dpi = 300
)

ggsave(
  filename = "outputs/psa_net_cost_hist.png",
  plot = ggplot(psa_draws, aes(net_cost)) +
    geom_histogram(bins = 40) +
    labs(x = "USD (4-year)", y = "Count") +
    theme_minimal(),
  width = 7,
  height = 4,
  dpi = 300
)

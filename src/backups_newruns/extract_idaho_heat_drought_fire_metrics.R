# ============================================================
# Idaho heat + drought + fire metrics workflow
# Built to mirror the original fire_metrics2 structure as closely
# as possible, while extending it to the climate-comparable heat
# and drought questions for the five Idaho towns only.
# ============================================================

library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(lubridate)
library(zoo)
library(ggplot2)
library(trend)
library(forecast)

verbose <- FALSE

# -----------------------------
# Settings
# -----------------------------
project_root <- tryCatch(normalizePath(getwd(), mustWork = FALSE), error = function(e) getwd())

towns_keep <- c("Deary", "Bovill", "Kendrick", "Juliaetta", "ElkRiver")
towns_pretty <- c(
  Deary = "Deary",
  Bovill = "Bovill",
  Kendrick = "Kendrick",
  Juliaetta = "Juliaetta",
  ElkRiver = "Elk River"
)

baseline_years <- 1996:2025
early_years    <- 1996:2005
late_years     <- 2016:2025
wateryearselect <- 1980:2025
use_anomalies <- TRUE

# -----------------------------
# Helpers
# -----------------------------
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

safe_dir_first <- function(paths) {
  hit <- paths[file.exists(paths) & dir.exists(paths)]
  if (length(hit) == 0) stop("Could not find any candidate data directory.")
  normalizePath(hit[[1]], mustWork = TRUE)
}

safe_file_first <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) stop("Could not find any candidate file path.")
  normalizePath(hit[[1]], mustWork = TRUE)
}

std_town <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_replace_all("[[:space:]]+", "") %>%
    str_replace_all("[^A-Za-z0-9]", "")
}

to_num <- function(x) suppressWarnings(as.numeric(x))

max_run <- function(x) {
  x <- as.integer(x)
  if (length(x) == 0 || all(is.na(x))) return(NA_integer_)
  x[is.na(x)] <- 0L
  r <- rle(x)
  if (!any(r$values == 1L)) return(0L)
  max(r$lengths[r$values == 1L])
}

first_sustained_day <- function(high_vec, day_vec, k = 3) {
  hv <- as.integer(high_vec)
  hv[is.na(hv)] <- 0L
  if (length(hv) < k) return(NA_integer_)
  roll <- zoo::rollapply(hv, width = k, FUN = sum, align = "left", fill = NA)
  idx <- which(roll == k)
  if (length(idx) == 0) return(NA_integer_)
  as.integer(day_vec[min(idx)])
}

last_sustained_day <- function(high_vec, day_vec, k = 3) {
  hv <- as.integer(high_vec)
  hv[is.na(hv)] <- 0L
  if (length(hv) < k) return(NA_integer_)
  roll <- zoo::rollapply(hv, width = k, FUN = sum, align = "left", fill = NA)
  idx <- which(roll == k)
  if (length(idx) == 0) return(NA_integer_)
  as.integer(day_vec[max(idx)])
}

total_days_in_long_runs <- function(x, min_run = 6) {
  xv <- as.integer(x)
  xv[is.na(xv)] <- 0L
  r <- rle(xv)
  if (!any(r$values == 1L & r$lengths >= min_run)) return(0L)
  sum(r$lengths[r$values == 1L & r$lengths >= min_run])
}

safe_sen_slope_per_decade <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  as.numeric(tryCatch(trend::sens.slope(x)$estimates * 10, error = function(e) NA_real_))
}

complete_years_vec <- function(years, vals) {
  years <- as.integer(years)
  vals <- as.numeric(vals)
  ok <- is.finite(years)
  years <- years[ok]
  vals <- vals[ok]
  if (length(years) == 0) return(numeric(0))
  yrs <- seq(min(years), max(years))
  tibble(wateryear = years, v = vals) %>%
    distinct(wateryear, .keep_all = TRUE) %>%
    right_join(tibble(wateryear = yrs), by = "wateryear") %>%
    arrange(wateryear) %>%
    pull(v)
}

arima_no_drift_qcat_diag <- function(x, n_sim = 2000, seed = 2026) {
  set.seed(seed)
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 10) {
    return(list(qcat = NA_character_, qbreaks = rep(NA_real_, 8), obs_slope = NA_real_, slopes_sim = numeric(0)))
  }

  fit <- tryCatch(
    forecast::auto.arima(x, allowdrift = FALSE, allowmean = TRUE),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(list(qcat = NA_character_, qbreaks = rep(NA_real_, 8), obs_slope = NA_real_, slopes_sim = numeric(0)))
  }

  obs_slope <- safe_sen_slope_per_decade(x)
  slopes_sim <- numeric(n_sim)

  for (i in seq_len(n_sim)) {
    simx <- tryCatch(
      as.numeric(simulate(fit, nsim = length(x))),
      error = function(e) {
        r <- tryCatch(na.omit(residuals(fit)), error = function(e2) numeric(0))
        mu <- mean(x, na.rm = TRUE)
        if (length(r) > 5) mu + sample(r, length(x), replace = TRUE) else rnorm(length(x), mean = mu, sd = sd(x, na.rm = TRUE))
      }
    )
    slopes_sim[i] <- safe_sen_slope_per_decade(simx)
  }

  slopes_sim <- slopes_sim[is.finite(slopes_sim)]
  if (!is.finite(obs_slope) || length(slopes_sim) < 50) {
    return(list(qcat = NA_character_, qbreaks = rep(NA_real_, 8), obs_slope = obs_slope, slopes_sim = slopes_sim))
  }

  qbreaks <- as.numeric(quantile(slopes_sim, probs = seq(0, 1, length.out = 8), na.rm = TRUE, type = 7))

  # Guard against non-unique breaks for flat or near-flat simulated series
  if (length(unique(qbreaks)) < length(qbreaks)) {
    eps <- max(sd(slopes_sim, na.rm = TRUE), 1e-8)
    center <- median(slopes_sim, na.rm = TRUE)
    qbreaks <- seq(center - 3 * eps, center + 3 * eps, length.out = 8)
  }

  qcat <- cut(obs_slope,
              breaks = qbreaks,
              include.lowest = TRUE,
              labels = paste0("Q", 1:7)) %>%
    as.character()

  list(qcat = qcat, qbreaks = qbreaks, obs_slope = obs_slope, slopes_sim = slopes_sim)
}

compute_delta_one_series <- function(df_one, end_year, tenure_years, min_years = 8,
                                     delta_method = c("slope", "diff"), k_years = 5) {
  delta_method <- match.arg(delta_method)
  if (!is.finite(tenure_years) || tenure_years <= 0) return(NA_real_)

  df_one <- df_one %>% filter(!is.na(value), !is.na(year)) %>% arrange(year)
  if (nrow(df_one) == 0) return(NA_real_)

  max_year <- min(end_year, max(df_one$year, na.rm = TRUE))
  start_year <- max_year - floor(tenure_years) + 1

  df_w <- df_one %>% filter(year >= start_year, year <= max_year) %>% arrange(year)
  if (nrow(df_w) < min_years) return(NA_real_)

  if (delta_method == "slope") {
    fit <- try(stats::lm(value ~ year, data = df_w), silent = TRUE)
    if (inherits(fit, "try-error")) return(NA_real_)
    slope_per_year <- unname(stats::coef(fit)[["year"]])
    return(10 * slope_per_year)
  }

  k <- min(k_years, floor(nrow(df_w) / 2))
  if (k < 2) return(NA_real_)
  mean(df_w$value[(nrow(df_w) - k + 1):nrow(df_w)], na.rm = TRUE) - mean(df_w$value[1:k], na.rm = TRUE)
}

compute_pdi_tenure <- function(survey, annual_long, survey_q_map, q_metric_map,
                               town_col = "town", id_col = "record_id", tenure_col = "tenure_years",
                               end_year = 2025, delta_method = c("slope", "diff"), k_years = 5, min_years = 8) {
  delta_method <- match.arg(delta_method)

  perception_long <- survey %>%
    mutate(
      .town = .data[[town_col]],
      .id = .data[[id_col]],
      tenure_years = to_num(.data[[tenure_col]])
    ) %>%
    select(.id, .town, tenure_years, all_of(unique(survey_q_map$survey_col))) %>%
    pivot_longer(cols = all_of(unique(survey_q_map$survey_col)), names_to = "survey_item", values_to = "resp") %>%
    left_join(survey_q_map, by = c("survey_item" = "survey_col")) %>%
    mutate(
      resp = to_num(resp),
      resp = if_else(direction == "reverse" & !is.na(resp), 8 - resp, resp)
    ) %>%
    filter(!is.na(question)) %>%
    select(.id, .town, tenure_years, question, resp)

  person_q_metric <- perception_long %>%
    left_join(q_metric_map, by = "question") %>%
    filter(!is.na(metric))

  stopifnot(all(c("town", "year", "metric", "value") %in% names(annual_long)))

  deltas <- person_q_metric %>%
    left_join(annual_long, by = c(".town" = "town", "metric" = "metric")) %>%
    group_by(.id, .town, tenure_years, question, resp, metric) %>%
    summarise(
      delta = compute_delta_one_series(
        df_one = tibble(year = year, value = value),
        end_year = end_year,
        tenure_years = unique(tenure_years)[1],
        min_years = min_years,
        delta_method = delta_method,
        k_years = k_years
      ),
      .groups = "drop"
    )

  metric_z <- deltas %>%
    group_by(question, metric) %>%
    mutate(delta_z = ifelse(sd(delta, na.rm = TRUE) > 0, as.numeric(scale(delta)), NA_real_)) %>%
    ungroup()

  pdi_person <- metric_z %>%
    group_by(.id, .town, tenure_years, question, resp) %>%
    summarise(
      n_metrics = sum(!is.na(delta_z)),
      climate_delta_mean_z = mean(delta_z, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      pdi = ifelse(is.na(resp) | is.na(climate_delta_mean_z), NA_real_, resp - climate_delta_mean_z)
    )

  pdi_summary <- pdi_person %>%
    group_by(question, .town) %>%
    summarise(
      n = sum(!is.na(pdi)),
      mean_resp = mean(resp, na.rm = TRUE),
      mean_climate_delta_z = mean(climate_delta_mean_z, na.rm = TRUE),
      mean_pdi = mean(pdi, na.rm = TRUE),
      sd_pdi = sd(pdi, na.rm = TRUE),
      .groups = "drop"
    )

  list(person = pdi_person, summary = pdi_summary, deltas = deltas)
}

# -----------------------------
# Paths
# -----------------------------
town_dir <- safe_dir_first(c(
  file.path(project_root, "data", "ID_data", "ID"),
  "/mnt/ceph/erichs/git/WWL-climate-extraction/data/ID_data/ID",
  "/mnt/data"
))

survey_path <- safe_file_first(c(
  file.path(project_root, "data", "FW2L_DC.csv"),
  file.path(project_root, "FW2L_DC.csv"),
  "/mnt/data/FW2L_DC.csv"
))

plot_dir <- file.path(project_root, "plots", "climate_questions")
out_dir  <- file.path(project_root, "outputs", "climate_questions")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Climate file inventory + combined table
# -----------------------------
files <- list.files(town_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("No CSV files found in: ", town_dir)

file_index <- tibble(
  file = files,
  filename = basename(files),
  stem = str_remove(filename, "\\.csv$")
) %>%
  mutate(
    buffer = if_else(str_detect(stem, "10mibuffer"), "10mibuffer", "none"),
    stem_nobuf = str_remove(stem, "10mibuffer"),
    source = str_extract(stem_nobuf, "(?<=_)(era5land|gridmet|nclim)$"),
    town = str_remove(stem_nobuf, "_(era5land|nclim|gridmet)$"),
    town = std_town(town)
  ) %>%
  filter(!is.na(source), !is.na(town)) %>%
  filter(town %in% towns_keep)

if (nrow(file_index) == 0) stop("No matching Idaho town climate files were found for the five-town scope.")

if (verbose) print(file_index)

town_tbl <- file_index %>%
  filter(source == "gridmet") %>%
  mutate(data = map(file, ~ read_csv(.x, show_col_types = FALSE))) %>%
  select(town, buffer, source, data) %>%
  unnest(data)

if (nrow(town_tbl) == 0) stop("No GridMET rows were available after reading the Idaho town files.")

# -----------------------------
# Time fields
# -----------------------------
town_tbl <- town_tbl %>%
  mutate(
    date_str = str_extract(`system.index`, "^\\d{8}"),
    date = as.Date(date_str, format = "%Y%m%d"),
    year = year(date),
    month = month(date),
    wateryear = if_else(month >= 10, year + 1L, year),
    yday = yday(date),
    mday = mday(date),
    season2 = if_else(month %in% 4:9, "growing", "cool"),
    season4 = case_when(
      month %in% 3:5  ~ "MAM",
      month %in% 6:8  ~ "JJA",
      month %in% 9:11 ~ "SON",
      TRUE ~ "DJF"
    ),
    wyday = ifelse(
      year(date) == wateryear & leap_year(year(date)), yday + 92,
      ifelse(year(date) == wateryear, yday + 91,
             ifelse(leap_year(wateryear - 1), yday - 274, yday - 273))
    )
  )

stopifnot(!all(is.na(town_tbl$date)))

# Convert temperatures if they appear to be Kelvin
if (median(town_tbl$tmmx, na.rm = TRUE) > 100) town_tbl$tmmx <- town_tbl$tmmx - 273.15
if (median(town_tbl$tmmn, na.rm = TRUE) > 100) town_tbl$tmmn <- town_tbl$tmmn - 273.15

town_tbl <- town_tbl %>% mutate(drydays = if_else(pr <= 1, 1L, 0L))

# -----------------------------
# Daily climatology, anomalies, and smoothed series
# -----------------------------
climo_vars <- intersect(c("bi", "erc", "eto", "pr", "tmmn", "tmmx"), names(town_tbl))

# Mean daily climatology
climo_mean_long <- town_tbl %>%
  filter(wateryear %in% baseline_years) %>%
  group_by(town, buffer, month, mday) %>%
  summarise(across(all_of(climo_vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(cols = all_of(climo_vars), names_to = "name", values_to = "dailyclimo")

dfanoms <- town_tbl %>%
  pivot_longer(cols = all_of(climo_vars), names_to = "name", values_to = "daily") %>%
  left_join(climo_mean_long, by = c("town", "buffer", "month", "mday", "name")) %>%
  mutate(dailyanom = daily - dailyclimo) %>%
  select(town, buffer, date, year, month, wateryear, yday, wyday, mday, season2, season4, name, dailyanom) %>%
  pivot_wider(names_from = name, values_from = dailyanom)

quibble <- function(x, q = c(0.9, 0.95), dropNA = TRUE) tibble(x = quantile(x, q, na.rm = dropNA), q = q)

climo_quant_long <- town_tbl %>%
  filter(wateryear %in% baseline_years) %>%
  group_by(town, buffer, month, mday) %>%
  reframe(across(all_of(intersect(c("bi", "erc", "pr", "tmmn", "tmmx"), names(town_tbl))), ~ quibble(.x, q = c(0.9, 0.95), dropNA = TRUE))) %>%
  pivot_longer(cols = -c(town, buffer, month, mday), names_to = "var", values_to = "qb") %>%
  mutate(
    quantile = map_dbl(qb, "q"),
    value = map_dbl(qb, "x"),
    name = paste0("dailyclimo9625_", var, "_q", quantile)
  ) %>%
  select(town, buffer, month, mday, name, value)

df_climo <- town_tbl %>%
  left_join(climo_quant_long, by = c("town", "buffer", "month", "mday")) %>%
  pivot_wider(names_from = name, values_from = value)

# Smoothing similar to Christine/original fire workflow
rollmean_by_group <- function(df, cols, k = 5) {
  df %>%
    group_by(town, buffer) %>%
    arrange(date, .by_group = TRUE) %>%
    mutate(across(all_of(cols), ~ zoo::rollmean(.x, k = k, align = "center", fill = NA_real_))) %>%
    ungroup()
}

df5day <- rollmean_by_group(town_tbl, cols = intersect(c("bi", "erc", "tmmn", "tmmx", "pr", "eto"), names(town_tbl)), k = 5)
df_fireclimo <- rollmean_by_group(df_climo, cols = intersect(c("bi", "erc"), names(df_climo)), k = 5)

# -----------------------------
# Fire metrics (mirroring original fire workflow)
# -----------------------------
df_fireclimo <- df_fireclimo %>%
  mutate(
    high_bi90 = bi > dailyclimo9625_bi_q0.9,
    high_erc90 = erc > dailyclimo9625_erc_q0.9,
    high_bi95 = bi > dailyclimo9625_bi_q0.95,
    high_erc95 = erc > dailyclimo9625_erc_q0.95,
    sev_bi90 = pmax(bi - dailyclimo9625_bi_q0.9, 0),
    sev_erc90 = pmax(erc - dailyclimo9625_erc_q0.9, 0)
  )

annualfire <- df_fireclimo %>%
  filter(wateryear %in% wateryearselect) %>%
  group_by(town, buffer, wateryear) %>%
  summarise(
    firehighdays_bi90 = sum(high_bi90, na.rm = TRUE),
    firehighdays_erc90 = sum(high_erc90, na.rm = TRUE),
    firehighdays_bi95 = sum(high_bi95, na.rm = TRUE),
    firehighdays_erc95 = sum(high_erc95, na.rm = TRUE),
    fireseverity_bi90 = sum(sev_bi90, na.rm = TRUE),
    fireseverity_erc90 = sum(sev_erc90, na.rm = TRUE),
    firesummerhighdays_bi90 = sum(high_bi90 & season4 == "JJA", na.rm = TRUE),
    firesummerhighdays_erc90 = sum(high_erc90 & season4 == "JJA", na.rm = TRUE),
    firestart_wyday_bi = first_sustained_day(high_bi90, wyday, k = 3),
    firestart_wyday_erc = first_sustained_day(high_erc90, wyday, k = 3),
    firestart_any_bi90 = suppressWarnings(min(wyday[which(high_bi90 %in% TRUE)], na.rm = TRUE)),
    firestart_any_erc90 = suppressWarnings(min(wyday[which(high_erc90 %in% TRUE)], na.rm = TRUE)),
    firestart_k2_erc90 = first_sustained_day(high_erc90, wyday, k = 2),
    firestart_k3_erc90 = first_sustained_day(high_erc90, wyday, k = 3),
    firestart_k3_bi90 = first_sustained_day(high_bi90, wyday, k = 3),
    firemaxrun_bi = max_run(high_bi90),
    firemaxrun_erc = max_run(high_erc90),
    fireseasonlength_bi = {
      idx <- which(high_bi90 %in% TRUE)
      if (length(idx) < 2) NA_integer_ else max(wyday[idx], na.rm = TRUE) - min(wyday[idx], na.rm = TRUE)
    },
    fireseasonlength_erc = {
      idx <- which(high_erc90 %in% TRUE)
      if (length(idx) < 2) NA_integer_ else max(wyday[idx], na.rm = TRUE) - min(wyday[idx], na.rm = TRUE)
    },
    firepeak_bi = max(bi, na.rm = TRUE),
    firepeak_erc = max(erc, na.rm = TRUE),
    firep95_erc = as.numeric(quantile(erc, probs = 0.95, na.rm = TRUE, type = 7)),
    .groups = "drop"
  ) %>%
  mutate(across(starts_with("firestart_any_"), ~ ifelse(is.infinite(.x), NA_real_, .x)))

# -----------------------------
# Heat metrics (standard temperature metrics only)
# -----------------------------
df_temp5 <- df5day %>% mutate(tmean = (tmmn + tmmx) / 2)

dfdaily90climos <- df_climo %>%
  mutate(
    tmaxgt90 = tmmx > dailyclimo9625_tmmx_q0.9,
    tmingt90 = tmmn > dailyclimo9625_tmmn_q0.9
  )

annualtemp <- town_tbl %>%
  filter(wateryear %in% wateryearselect) %>%
  group_by(town, buffer, wateryear) %>%
  summarise(
    meananntemp = mean((tmmn + tmmx) / 2, na.rm = TRUE),
    meanminanntemp = mean(tmmn, na.rm = TRUE),
    meanmaxanntemp = mean(tmmx, na.rm = TRUE),
    maxminjjatemp = max(tmmn[season4 == "JJA"], na.rm = TRUE),
    maxmaxjjatemp = max(tmmx[season4 == "JJA"], na.rm = TRUE),
    meanjjatemp = mean(((tmmn + tmmx) / 2)[season4 == "JJA"], na.rm = TRUE),
    frostdays = sum(tmmn < 0, na.rm = TRUE),
    summerdays = sum(tmmx > 25, na.rm = TRUE),
    .groups = "drop"
  )

annualtempanom <- dfanoms %>%
  filter(wateryear %in% wateryearselect) %>%
  group_by(town, buffer, wateryear) %>%
  summarise(
    meananomanntemp = mean((tmmn + tmmx) / 2, na.rm = TRUE),
    meanminanomanntemp = mean(tmmn, na.rm = TRUE),
    meanmaxanomanntemp = mean(tmmx, na.rm = TRUE),
    maxminanomjjatemp = max(tmmn[season4 == "JJA"], na.rm = TRUE),
    maxmaxanomjjatemp = max(tmmx[season4 == "JJA"], na.rm = TRUE),
    meananomjjatemp = mean(((tmmn + tmmx) / 2)[season4 == "JJA"], na.rm = TRUE),
    .groups = "drop"
  )

annualtemp <- annualtemp %>%
  left_join(
    dfdaily90climos %>%
      filter(wateryear %in% wateryearselect) %>%
      group_by(town, buffer, wateryear) %>%
      summarise(
        tmaxgt90thdays = sum(tmaxgt90, na.rm = TRUE),
        tmingt90thdays = sum(tmingt90, na.rm = TRUE),
        durtmaxgt90th = total_days_in_long_runs(tmaxgt90, min_run = 6),
        .groups = "drop"
      ),
    by = c("town", "buffer", "wateryear")
  ) %>%
  left_join(
    df_temp5 %>%
      filter(wateryear %in% wateryearselect) %>%
      group_by(town, buffer, wateryear) %>%
      summarise(
        gsstartdoy = first_sustained_day(tmean > 10, yday, k = 5),
        gsenddoy   = last_sustained_day(tmean > 10, yday, k = 5),
        .groups = "drop"
      ) %>%
      mutate(gslength = gsenddoy - gsstartdoy),
    by = c("town", "buffer", "wateryear")
  )

# -----------------------------
# Drought metrics
# -----------------------------
df_drought <- town_tbl %>% arrange(town, buffer, date)

calc_dry_dur <- function(x) {
  out <- integer(length(x))
  run <- 0L
  for (i in seq_along(x)) {
    if (is.na(x[i]) || x[i] == 0L) {
      run <- 0L
      out[i] <- 0L
    } else {
      run <- run + 1L
      out[i] <- run
    }
  }
  out
}

df_drought <- df_drought %>%
  group_by(town, buffer) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(drydaydur = calc_dry_dur(drydays)) %>%
  ungroup()

precip_timing <- df_drought %>%
  group_by(town, buffer) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(
    pr20day = zoo::rollsum(pr, 20, align = "left", fill = NA),
    dryday10 = zoo::rollsum(drydays, 10, align = "left", fill = NA),
    pr2day = zoo::rollsum(pr, 2, align = "right", fill = NA)
  ) %>%
  ungroup() %>%
  filter(wateryear %in% wateryearselect)

annualdrought <- df_drought %>%
  filter(wateryear %in% wateryearselect) %>%
  group_by(town, buffer, wateryear) %>%
  summarise(
    annualprecip = sum(pr, na.rm = TRUE),
    annualeto = sum(eto, na.rm = TRUE),
    annualpwd = annualprecip - annualeto,
    coolseasonprecip = sum(pr[season2 == "cool"], na.rm = TRUE),
    coolseasondrydays1mm = sum(pr[season2 == "cool"] < 1, na.rm = TRUE),
    maxdrydur = max(drydaydur, na.rm = TRUE),
    maxcoolseasondrydur = max(drydaydur[season2 == "cool"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    precip_timing %>%
      filter(pr2day >= 10, pr20day >= 5, dryday10 < 10, drydays != 1) %>%
      group_by(town, buffer, wateryear) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      transmute(town, buffer, wateryear, preciponsetwydoy = wyday),
    by = c("town", "buffer", "wateryear")
  ) %>%
  left_join(
    df_drought %>%
      filter(wateryear %in% wateryearselect) %>%
      mutate(wpr = pr * wyday) %>%
      group_by(town, buffer, wateryear) %>%
      summarise(precipCOMwydoy = sum(wpr, na.rm = TRUE) / sum(pr, na.rm = TRUE), .groups = "drop"),
    by = c("town", "buffer", "wateryear")
  )

annualdroughtanoms <- dfanoms %>%
  filter(wateryear %in% wateryearselect) %>%
  group_by(town, buffer, wateryear) %>%
  summarise(
    annualanomprecip = sum(pr, na.rm = TRUE),
    annualanometo = sum(eto, na.rm = TRUE),
    annualanompwd = annualanomprecip - annualanometo,
    coolseasonanomprecip = sum(pr[season2 == "cool"], na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------
# Combine annual metrics
# -----------------------------
annual_all <- annualfire %>%
  full_join(annualtemp, by = c("town", "buffer", "wateryear")) %>%
  full_join(annualtempanom, by = c("town", "buffer", "wateryear")) %>%
  full_join(annualdrought, by = c("town", "buffer", "wateryear")) %>%
  full_join(annualdroughtanoms, by = c("town", "buffer", "wateryear")) %>%
  arrange(town, buffer, wateryear)

annual_all_long <- annual_all %>%
  pivot_longer(cols = -c(town, buffer, wateryear), names_to = "metric", values_to = "value") %>%
  arrange(town, buffer, metric, wateryear)

# -----------------------------
# Question maps
# -----------------------------
question_map_master <- list(
  H1 = c("meananomanntemp", "meanminanomanntemp", "meanmaxanomanntemp", "meananntemp", "meanminanntemp", "meanmaxanntemp"),
  H2 = c("maxminanomjjatemp", "maxmaxanomjjatemp", "meananomjjatemp", "maxminjjatemp", "maxmaxjjatemp", "meanjjatemp"),
  H3 = c("summerdays", "frostdays", "tmaxgt90thdays", "tmingt90thdays"),
  H4 = c("gsstartdoy"),
  H5 = c("gsenddoy"),
  H6 = c("gslength", "durtmaxgt90th"),
  D1 = c("annualprecip", "annualeto", "annualpwd", "annualanomprecip", "annualanometo", "annualanompwd"),
  D2 = c("coolseasonprecip", "coolseasonanomprecip"),
  D3 = c("coolseasondrydays1mm"),
  D4 = c("preciponsetwydoy", "precipCOMwydoy"),
  D5 = c("maxdrydur", "maxcoolseasondrydur"),
  F1 = c("firehighdays_erc90", "firehighdays_bi90", "firehighdays_erc95", "firehighdays_bi95"),
  F2 = c("firesummerhighdays_erc90", "firesummerhighdays_bi90"),
  F3 = c("firepeak_erc", "firepeak_bi", "firep95_erc", "fireseverity_erc90", "fireseverity_bi90"),
  F4 = c("firestart_wyday_erc", "firestart_wyday_bi", "firestart_any_erc90", "firestart_any_bi90", "firestart_k2_erc90", "firestart_k3_erc90", "firestart_k3_bi90"),
  F5 = c("fireseasonlength_erc", "fireseasonlength_bi", "firemaxrun_erc", "firemaxrun_bi")
)

question_labels <- tibble(
  question = names(question_map_master),
  stressor = c(rep("Heat", 6), rep("Drought", 5), rep("Fire", 5)),
  survey_col = c("ID6", "ID7", "ID7a", "ID9", "ID10", "ID11", "ID15", "ID16", "IDI6a", "ID17", "ID18", "ID22", "ID23", "ID24", "ID25", "ID26"),
  question_text = c(
    "Overall, has it gotten hotter or cooler?",
    "During the summer months, has it gotten hotter or cooler?",
    "If summer has changed, how frequently is it hotter?",
    "When does it start getting hotter in the spring?",
    "When does it start getting cooler in the fall?",
    "Do heat conditions last shorter or longer than you remember?",
    "Overall, has it gotten wetter or drier?",
    "During the winter months, has it become wetter or drier?",
    "If winter has changed, how frequently is it drier?",
    "Has there been a change in when rainfall/snowfall starts during the year?",
    "Do drier conditions last shorter or longer than you remember?",
    "Overall, has fire become less common or more common?",
    "During summer, has fire occurred more often or less often?",
    "Have the areas affected by wildfires gotten smaller or larger?",
    "Has there been a change in when the fire season starts during the year?",
    "Do fire conditions last shorter or longer than you remember?"
  )
)

q_metric_map <- tibble(
  question = rep(names(question_map_master), times = lengths(question_map_master)),
  metric = unlist(question_map_master)
) %>%
  left_join(question_labels, by = "question")

# direction is inferred and easy to edit later.
survey_q_map <- question_labels %>%
  mutate(
    direction = c(
      # Heat: leave as normal working assumption
      "normal", "normal", "normal", "reverse", "normal", "normal",
      # Drought: wetter/drier items assumed normal; earlier onset reverse; duration normal
      "normal", "normal", "normal", "reverse", "normal",
      # Fire
      "normal", "normal", "normal", "reverse", "normal"
    )
  )

# -----------------------------
# Trend / anomaly / delta outputs by metric
# -----------------------------
metric_anoms <- annual_all_long %>%
  group_by(town, buffer, metric) %>%
  mutate(
    baseline_mean = mean(value[wateryear %in% baseline_years], na.rm = TRUE),
    anom = value - baseline_mean
  ) %>%
  ungroup()

metric_delta <- annual_all_long %>%
  group_by(town, buffer, metric) %>%
  summarise(
    early_mean = mean(value[wateryear %in% early_years], na.rm = TRUE),
    late_mean = mean(value[wateryear %in% late_years], na.rm = TRUE),
    delta_late_minus_early = late_mean - early_mean,
    .groups = "drop"
  )

metric_arima_diag <- metric_anoms %>%
  group_by(town, buffer, metric) %>%
  summarise(
    years = list(wateryear),
    vals = list(if (use_anomalies) anom else value),
    .groups = "drop"
  ) %>%
  mutate(
    series_complete = map2(years, vals, complete_years_vec),
    diag = map(series_complete, ~ arima_no_drift_qcat_diag(.x, n_sim = 2000, seed = 2026)),
    sen_slope_per_decade = map_dbl(diag, "obs_slope"),
    qcat = map_chr(diag, ~ .x$qcat %||% NA_character_),
    p_upper = map_dbl(diag, ~ if (length(.x$slopes_sim) == 0 || !is.finite(.x$obs_slope)) NA_real_ else mean(.x$slopes_sim >= .x$obs_slope, na.rm = TRUE)),
    p_lower = map_dbl(diag, ~ if (length(.x$slopes_sim) == 0 || !is.finite(.x$obs_slope)) NA_real_ else mean(.x$slopes_sim <= .x$obs_slope, na.rm = TRUE)),
    p_two = 2 * pmin(p_upper, p_lower)
  ) %>%
  select(-years, -vals, -series_complete, -diag) %>%
  left_join(metric_anoms %>% distinct(town, buffer, metric, baseline_mean), by = c("town", "buffer", "metric")) %>%
  left_join(metric_delta, by = c("town", "buffer", "metric")) %>%
  left_join(q_metric_map, by = "metric") %>%
  arrange(stressor, question, metric, town, buffer)

question_outputs <- metric_arima_diag
write_csv(question_outputs, file.path(out_dir, "idaho_question_outputs.csv"))
write_csv(annual_all, file.path(out_dir, "idaho_annual_metrics_wide.csv"))
write_csv(annual_all_long, file.path(out_dir, "idaho_annual_metrics_long.csv"))

# -----------------------------
# Survey data + PDI
# -----------------------------
survey_raw <- read_csv(survey_path, show_col_types = FALSE)

survey_prepped <- survey_raw %>%
  mutate(
    record_id = row_number(),
    town = recode(std_town(Nearest_Town), !!!setNames(names(towns_pretty), std_town(unname(towns_pretty)))),
    tenure_years = to_num(ID2)
  ) %>%
  filter(town %in% towns_keep)

pdi_results <- compute_pdi_tenure(
  survey = survey_prepped,
  annual_long = annual_all_long %>% rename(year = wateryear),
  survey_q_map = survey_q_map,
  q_metric_map = q_metric_map %>% select(question, metric),
  town_col = "town",
  id_col = "record_id",
  tenure_col = "tenure_years",
  end_year = 2025,
  delta_method = "slope",
  k_years = 5,
  min_years = 8
)

write_csv(pdi_results$summary, file.path(out_dir, "idaho_pdi_summary.csv"))
write_csv(pdi_results$person, file.path(out_dir, "idaho_pdi_person.csv"))
write_csv(pdi_results$deltas, file.path(out_dir, "idaho_pdi_deltas.csv"))

# -----------------------------
# Simple plot outputs, question-framed like original fire Rmd
# -----------------------------
plot_metric_timeseries <- function(metric_name, question_code) {
  dat <- annual_all_long %>% filter(metric == metric_name)
  if (nrow(dat) == 0) return(invisible(NULL))
  p <- ggplot(dat, aes(x = wateryear, y = value, group = interaction(town, buffer))) +
    geom_line(alpha = 0.6) +
    facet_grid(buffer ~ town, scales = "free_y", labeller = labeller(town = towns_pretty)) +
    labs(title = paste(question_code, metric_name), x = "Water year", y = metric_name)
  ggsave(file.path(plot_dir, paste0(question_code, "_timeseries_", metric_name, ".png")), p, width = 8, height = 5, dpi = 300)
}

plot_metric_distribution <- function(metric_name, question_code) {
  dat <- annual_all_long %>% filter(metric == metric_name)
  if (nrow(dat) == 0) return(invisible(NULL))
  p <- ggplot(dat, aes(x = value)) +
    geom_histogram(bins = 30) +
    facet_grid(buffer ~ town, scales = "free_x", labeller = labeller(town = towns_pretty)) +
    labs(title = paste(question_code, metric_name, "distribution"), x = metric_name, y = "Count")
  ggsave(file.path(plot_dir, paste0(question_code, "_distribution_", metric_name, ".png")), p, width = 8, height = 5, dpi = 300)
}

walk2(unlist(question_map_master), rep(names(question_map_master), lengths(question_map_master)), plot_metric_timeseries)
walk2(unlist(question_map_master), rep(names(question_map_master), lengths(question_map_master)), plot_metric_distribution)

if (verbose) {
  print(question_outputs, n = 50)
  print(pdi_results$summary, n = 50)
}

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


  person_instr <- deltas %>%
    group_by(.id, .town, tenure_years, question, resp) %>%
    summarise(
      n_metrics = sum(!is.na(delta)),
      instr_delta_raw = if (all(is.na(delta))) NA_real_ else mean(delta, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(question) %>%
    mutate(
      climate_delta_mean_z = {
        if (sd(instr_delta_raw, na.rm = TRUE) > 0) {
          as.numeric(scale(instr_delta_raw))
        } else {
          rep(NA_real_, dplyr::n())
        }
      },
      resp_z = {
        if (sd(resp, na.rm = TRUE) > 0) {
          as.numeric(scale(resp))
        } else {
          rep(NA_real_, dplyr::n())
        }
      }
    ) %>%
    ungroup() %>%
    mutate(
      pdi = ifelse(
        is.na(resp_z) | is.na(climate_delta_mean_z),
        NA_real_,
        resp_z - climate_delta_mean_z
      )
    )
  
  pdi_summary <- person_instr %>%
    group_by(question, .town) %>%
    summarise(
      n = sum(!is.na(pdi)),
      mean_resp = mean(resp, na.rm = TRUE),
      mean_resp_z = mean(resp_z, na.rm = TRUE),
      mean_instr_delta_raw = mean(instr_delta_raw, na.rm = TRUE),
      mean_climate_delta_z = mean(climate_delta_mean_z, na.rm = TRUE),
      mean_pdi = mean(pdi, na.rm = TRUE),
      sd_pdi = sd(pdi, na.rm = TRUE),
      .groups = "drop"
    )
  
  list(person = person_instr, summary = pdi_summary, deltas = deltas)
}
# -----------------------------
# Paths (kept close to original fire workflow)
# -----------------------------
town_dir <- "/mnt/ceph/erichs/git/WWL-climate-extraction/data/ID_data/ID/"
if (!dir.exists(town_dir)) town_dir <- file.path(project_root, "data", "ID_data", "ID")
if (!dir.exists(town_dir)) town_dir <- "/mnt/data"
stopifnot(dir.exists(town_dir))

survey_path <- file.path(project_root, "data/survey", "FW2L_DC.csv")
if (!file.exists(survey_path)) survey_path <- file.path(project_root, "FW2L_DC.csv")
if (!file.exists(survey_path)) survey_path <- "/mnt/data/FW2L_DC.csv"
stopifnot(file.exists(survey_path))

plot_dir <- file.path(project_root, "plots")
out_dir  <- file.path(project_root, "outputs")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Climate file inventory + combined table (mirrors original fire script)
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
    town   = str_remove(stem_nobuf, "_(era5land|gridmet|nclim)$")
  )

bad <- file_index %>% filter(is.na(source) | is.na(town))
if (nrow(bad) > 0) {
  if (verbose) print(bad)
  stop("Some filenames did not match expected pattern Town[_10mibuffer]_(era5land|gridmet|nclim).csv")
}

file_index <- file_index %>%
  mutate(town = std_town(town)) %>%
  filter(town %in% towns_keep)

if (nrow(file_index) == 0) stop("No matching Idaho town climate files were found for the five-town scope.")
if (verbose) print(file_index)

town_data <- file_index %>%
  mutate(obj_name = paste(town, buffer, source, sep = "__")) %>%
  select(obj_name, file) %>%
  tibble::deframe() %>%
  imap(\(f, nm) read_csv(f, show_col_types = FALSE))

inventory <- imap_dfr(town_data, \(df, nm) tibble(object = nm, n = nrow(df), p = ncol(df))) %>%
  arrange(object)
if (verbose) print(inventory, n = Inf)

town_tbl <- file_index %>%
  mutate(data = map(file, \(f) read_csv(f, show_col_types = FALSE))) %>%
  select(town, buffer, source, data) %>%
  unnest(data) %>%
  filter(source == "gridmet")

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

quant_vars <- intersect(c("bi", "erc", "pr", "tmmn", "tmmx"), names(town_tbl))

climo_quant <- town_tbl %>%
  filter(wateryear %in% baseline_years) %>%
  group_by(town, buffer, month, mday) %>%
  summarise(
    across(
      all_of(quant_vars),
      list(
        q0.9 = ~ quantile(.x, probs = 0.9, na.rm = TRUE, type = 7),
        q0.95 = ~ quantile(.x, probs = 0.95, na.rm = TRUE, type = 7)
      ),
      .names = "dailyclimo9625_{.col}_{.fn}"
    ),
    .groups = "drop"
  )

df_climo <- town_tbl %>%
  left_join(climo_quant, by = c("town", "buffer", "month", "mday"))

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
  ) %>% select(-years, -vals, -series_complete, -diag) %>%
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
# Optional ACS demographic context (population-weighted 10-mile buffers)
# Built to mirror the demographic context section in the original fire report,
# but extended so the resulting town-level covariates can be joined to Heat,
# Drought, and Fire PDI summaries.
# -----------------------------
acs_year <- 2022
acs_buffer_miles <- 10
acs_enabled <- all(vapply(
  c("sf", "tigris", "tidycensus", "units"),
  requireNamespace,
  logical(1),
  quietly = TRUE
))

if (acs_enabled) {
  message("ACS demographic context: packages detected. Attempting block-group weighted buffer analysis.")

  acs_result <- tryCatch({
    id_places <- tigris::places(state = "ID", year = acs_year, class = "sf", progress_bar = FALSE) %>%
      dplyr::mutate(town = std_town(NAME)) %>%
      dplyr::filter(town %in% towns_keep)

    if (nrow(id_places) == 0) stop("No Idaho place geometries matched the five target towns.")

    # Use interior representative points so buffering is stable.
    town_pts <- id_places %>%
      dplyr::select(town) %>%
      sf::st_point_on_surface()

    # Pull Idaho block-group ACS with geometry.
    # Poverty is derived from count variables (C17002) after overlap weighting,
    # rather than averaging a precomputed percentage. This is more robust for
    # population-weighted buffer summaries.
    acs_vars <- c(
      total_pop = "B01003_001",
      medinc = "B19013_001",
      pov_denom = "C17002_001",
      pov_below50 = "C17002_002",
      pov_50to99 = "C17002_003",
      edu_denom = "B15003_001",
      bach = "B15003_022",
      mast = "B15003_023",
      prof = "B15003_024",
      doc = "B15003_025",
      m65_66 = "B01001_020",
      m67_69 = "B01001_021",
      m70_74 = "B01001_022",
      m75_79 = "B01001_023",
      m80_84 = "B01001_024",
      m85p = "B01001_025",
      f65_66 = "B01001_044",
      f67_69 = "B01001_045",
      f70_74 = "B01001_046",
      f75_79 = "B01001_047",
      f80_84 = "B01001_048",
      f85p = "B01001_049"
    )

    id_bg <- tidycensus::get_acs(
      geography = "block group",
      state = "ID",
      year = acs_year,
      survey = "acs5",
      variables = acs_vars,
      geometry = TRUE,
      output = "wide",
      cache_table = TRUE
    )

    bg <- id_bg %>%
      dplyr::transmute(
        GEOID,
        total_pop = total_popE,
        medinc = medincE,
        pov_denom = pov_denomE,
        pov_num = pov_below50E + pov_50to99E,
        edu_denom = edu_denomE,
        bachplus_num = bachE + mastE + profE + docE,
        age65_num = m65_66E + m67_69E + m70_74E + m75_79E + m80_84E + m85pE +
          f65_66E + f67_69E + f70_74E + f75_79E + f80_84E + f85pE
      )

    crs_use <- 5070
    town_buf <- town_pts %>%
      sf::st_transform(crs_use) %>%
      sf::st_buffer(units::set_units(acs_buffer_miles, "miles")) %>%
      dplyr::mutate(buffer = "10mibuffer")

    bg_proj <- bg %>% sf::st_transform(crs_use)
    bg_area <- sf::st_area(bg_proj)

    inter <- suppressWarnings(sf::st_intersection(
      town_buf %>% dplyr::select(town, buffer),
      bg_proj %>% dplyr::select(GEOID, total_pop, medinc, pov_denom, pov_num, edu_denom, bachplus_num, age65_num)
    )) %>%
      dplyr::mutate(
        int_area = sf::st_area(geometry)
      ) %>%
      dplyr::left_join(
        bg_proj %>% sf::st_drop_geometry() %>% dplyr::transmute(GEOID, bg_area = as.numeric(bg_area)),
        by = "GEOID"
      ) %>%
      dplyr::mutate(
        overlap_frac = dplyr::if_else(is.finite(bg_area) & bg_area > 0, as.numeric(int_area) / bg_area, NA_real_),
        pop_est = total_pop * overlap_frac,
        # Weight numerator and denominator counts separately, then compute the ratio after aggregation.
        pov_num_w = pov_num * overlap_frac,
        pov_denom_w = pov_denom * overlap_frac,
        bachplus_w = bachplus_num * overlap_frac,
        edu_denom_w = edu_denom * overlap_frac,
        age65_w = age65_num * overlap_frac,
        total_pop_w = total_pop * overlap_frac
      )

    acs_covars <- inter %>%
      sf::st_drop_geometry() %>%
      dplyr::group_by(town, buffer) %>%
      dplyr::summarise(
        n_bg_overlap = dplyr::n(),
        pop_est_in_buffer = sum(pop_est, na.rm = TRUE),
        median_income = stats::weighted.mean(medinc, w = pop_est, na.rm = TRUE),
        poverty_pct = dplyr::if_else(
          is.finite(sum(pov_denom_w, na.rm = TRUE)) & sum(pov_denom_w, na.rm = TRUE) > 0,
          100 * sum(pov_num_w, na.rm = TRUE) / sum(pov_denom_w, na.rm = TRUE),
          NA_real_
        ),
        bachelors_pct = dplyr::if_else(sum(edu_denom_w, na.rm = TRUE) > 0,
                                       100 * sum(bachplus_w, na.rm = TRUE) / sum(edu_denom_w, na.rm = TRUE),
                                       NA_real_),
        age65_pct = dplyr::if_else(sum(total_pop_w, na.rm = TRUE) > 0,
                                   100 * sum(age65_w, na.rm = TRUE) / sum(total_pop_w, na.rm = TRUE),
                                   NA_real_),
        .groups = "drop"
      )

    q_stressor_map <- question_labels %>% dplyr::select(question, stressor)

    pdi_summary_stressor <- pdi_results$summary %>%
      dplyr::left_join(q_stressor_map, by = "question") %>%
      dplyr::group_by(.town, stressor) %>%
      dplyr::summarise(
        n_questions = dplyr::n(),
        mean_pdi = mean(mean_pdi, na.rm = TRUE),
        mean_instr_delta_raw = mean(mean_instr_delta_raw, na.rm = TRUE),
        mean_climate_delta_z = mean(mean_climate_delta_z, na.rm = TRUE),
        .groups = "drop"
      )

    pdi_summary_overall <- pdi_results$summary %>%
      dplyr::group_by(.town) %>%
      dplyr::summarise(
        n_questions = dplyr::n(),
        mean_pdi_all = mean(mean_pdi, na.rm = TRUE),
        mean_instr_delta_raw_all = mean(mean_instr_delta_raw, na.rm = TRUE),
        mean_climate_delta_z_all = mean(mean_climate_delta_z, na.rm = TRUE),
        .groups = "drop"
      )

    pdi_question_acs <- pdi_results$summary %>%
      dplyr::left_join(q_stressor_map, by = "question") %>%
      dplyr::left_join(acs_covars, by = c(".town" = "town"))

    pdi_stressor_acs <- pdi_summary_stressor %>%
      dplyr::left_join(acs_covars, by = c(".town" = "town"))

    pdi_overall_acs <- pdi_summary_overall %>%
      dplyr::left_join(acs_covars, by = c(".town" = "town"))

    list(
      town_points = town_pts,
      buffers = town_buf,
      block_groups = bg_proj,
      overlaps = inter,
      acs_covars = acs_covars,
      pdi_question_acs = pdi_question_acs,
      pdi_stressor_acs = pdi_stressor_acs,
      pdi_overall_acs = pdi_overall_acs
    )
  }, error = function(e) {
    message("ACS demographic context skipped: ", conditionMessage(e))
    NULL
  })

  if (!is.null(acs_result)) {
    acs_covars <- acs_result$acs_covars
    pdi_question_acs <- acs_result$pdi_question_acs
    pdi_stressor_acs <- acs_result$pdi_stressor_acs
    pdi_overall_acs <- acs_result$pdi_overall_acs

    readr::write_csv(acs_covars, file.path(out_dir, "idaho_acs_buffer_covariates.csv"))
    readr::write_csv(pdi_question_acs, file.path(out_dir, "idaho_pdi_question_acs.csv"))
    readr::write_csv(pdi_stressor_acs, file.path(out_dir, "idaho_pdi_stressor_acs.csv"))
    readr::write_csv(pdi_overall_acs, file.path(out_dir, "idaho_pdi_overall_acs.csv"))

    question_display_map <- c(
      H1 = "H1 overall heat",
      H2 = "H2 summer heat",
      H3 = "H3 heat frequency",
      H4 = "H4 spring onset",
      H5 = "H5 fall cooling",
      H6 = "H6 heat duration",
      D1 = "D1 overall wet/dry",
      D2 = "D2 winter wet/dry",
      D3 = "D3 winter dryness",
      D4 = "D4 precip timing",
      D5 = "D5 dry duration",
      F1 = "F1 fire frequency",
      F2 = "F2 summer fire freq.",
      F3 = "F3 fire intensity",
      F4 = "F4 fire start",
      F5 = "F5 fire duration"
    )
    town_display_order <- c("Bovill", "Deary", "Elk River", "Juliaetta", "Kendrick")

    qcat_panel <- question_outputs %>%
      dplyr::filter(buffer == "10mibuffer") %>%
      dplyr::mutate(qcat_num = as.numeric(gsub("Q", "", qcat))) %>%
      dplyr::group_by(question, town) %>%
      dplyr::summarise(qcat_num = round(mean(qcat_num, na.rm = TRUE)), .groups = "drop") %>%
      dplyr::mutate(
        town_display = factor(dplyr::recode(town, !!!towns_pretty), levels = town_display_order),
        question_display = factor(
          question_display_map[as.character(question)],
          levels = rev(unname(question_display_map))
        )
      )

    pdi_panel <- pdi_results$summary %>%
      dplyr::mutate(
        town_display = factor(dplyr::recode(.town, !!!towns_pretty), levels = town_display_order),
        question_display = factor(
          question_display_map[as.character(question)],
          levels = rev(unname(question_display_map))
        )
      )

    p_qcat <- ggplot2::ggplot(
      qcat_panel,
      ggplot2::aes(x = town_display, y = question_display, fill = qcat_num)
    ) +
      ggplot2::geom_tile() +
      ggplot2::geom_text(ggplot2::aes(label = qcat_num), size = 2.5, color = "white") +
      ggplot2::scale_fill_viridis_c(
        name = "qcat",
        limits = c(1, 7),
        breaks = 1:7,
        option = "viridis",
        guide = ggplot2::guide_colorbar(
          title.position = "top",
          barheight = grid::unit(5, "cm"),
          barwidth = grid::unit(0.45, "cm")
        )
      ) +
      ggplot2::coord_equal() +
      ggplot2::labs(title = "ARIMA qcat by question and town", x = NULL, y = NULL) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "plain", size = 11),
        legend.position = "right",
        legend.title = ggplot2::element_text(size = 8),
        legend.text = ggplot2::element_text(size = 8)
      )

    p_pdi <- ggplot2::ggplot(
      pdi_panel,
      ggplot2::aes(x = town_display, y = question_display, fill = mean_pdi)
    ) +
      ggplot2::geom_tile() +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", mean_pdi)), size = 2.3, color = "white") +
      ggplot2::scale_fill_viridis_c(
        name = "Mean PDI",
        option = "viridis",
        guide = ggplot2::guide_colorbar(
          title.position = "top",
          barheight = grid::unit(5, "cm"),
          barwidth = grid::unit(0.45, "cm")
        )
      ) +
      ggplot2::coord_equal() +
      ggplot2::labs(title = "Mean PDI by question and town", x = NULL, y = NULL) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "plain", size = 11),
        legend.position = "right",
        legend.title = ggplot2::element_text(size = 8),
        legend.text = ggplot2::element_text(size = 8)
      )

    png(file.path(plot_dir, "pdi_heatmap_by_stressor.png"), width = 1400, height = 950, res = 150)
    grid::grid.newpage()
    grid::grid.text("Trend distinctiveness and climate-perception alignment", y = 0.98, gp = grid::gpar(fontsize = 16))
    lay <- grid::grid.layout(nrow = 1, ncol = 2, widths = grid::unit(c(1, 1), "null"))
    grid::pushViewport(grid::viewport(y = 0.47, height = 0.86, layout = lay))
    print(p_qcat, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(p_pdi, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
    dev.off()

    p_income <- ggplot2::ggplot(
      pdi_stressor_acs,
      ggplot2::aes(x = median_income, y = mean_pdi, label = .town)
    ) +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::geom_text(nudge_y = 0.05, size = 3, show.legend = FALSE) +
      ggplot2::facet_wrap(~ stressor, scales = "free_y") +
      ggplot2::labs(
        title = "Town mean PDI vs median household income",
        subtitle = "ACS block-group weighted 10-mile buffer context",
        x = "Median household income",
        y = "Mean PDI"
      )
    ggplot2::ggsave(file.path(plot_dir, "pdi_vs_income_by_stressor.png"), p_income, width = 10, height = 4.5, dpi = 300)

    demo_long <- pdi_stressor_acs %>%
      tidyr::pivot_longer(
        cols = c(median_income, poverty_pct, bachelors_pct, age65_pct),
        names_to = "covariate",
        values_to = "covariate_value"
      )

    p_demo <- ggplot2::ggplot(
      demo_long,
      ggplot2::aes(x = covariate_value, y = mean_pdi, label = .town)
    ) +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::geom_text(nudge_y = 0.05, size = 3, show.legend = FALSE) +
      ggplot2::facet_grid(stressor ~ covariate, scales = "free_x") +
      ggplot2::labs(
        title = "Town mean PDI vs demographic context",
        subtitle = "ACS block-group weighted 10-mile buffer",
        x = NULL,
        y = "Mean PDI"
      )
    ggplot2::ggsave(file.path(plot_dir, "pdi_vs_demographics_by_stressor.png"), p_demo, width = 12, height = 8, dpi = 300)
  }
} else {
  message("ACS demographic context skipped: install sf, tigris, tidycensus, and units to enable it.")
}

# -----------------------------
# Simple plot outputs, question-framed like original fire Rmd
# -----------------------------
plot_metric_timeseries <- function(metric_name, question_code) {
  dat <- annual_all_long %>% dplyr::filter(metric == metric_name)
  if (nrow(dat) == 0) return(invisible(NULL))

  trend_df <- dat %>%
    dplyr::group_by(town, buffer) %>%
    dplyr::summarise(
      slope_per_decade = if (sum(is.finite(value) & is.finite(wateryear)) >= 2) {
        unname(stats::coef(stats::lm(value ~ wateryear, data = dplyr::cur_data()))["wateryear"] * 10)
      } else {
        NA_real_
      },
      x = min(wateryear, na.rm = TRUE),
      y = max(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(label = sprintf("Slope = %.2f / decade", slope_per_decade))

  p <- ggplot(dat, aes(x = wateryear, y = value)) +
    geom_line(linewidth = 0.4, color = "black") +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, color = "#2C7FB8") +
    geom_label(
      data = trend_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 2.8,
      label.size = 0.15,
      fill = "white"
    ) +
    facet_grid(buffer ~ town, scales = "free_y", labeller = labeller(town = towns_pretty)) +
    labs(title = paste(question_code, metric_name), x = "Water year", y = metric_name) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "grey70"),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
  ggsave(file.path(plot_dir, paste0(question_code, "_timeseries_", metric_name, ".png")), p, width = 8.5, height = 5.25, dpi = 300)
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

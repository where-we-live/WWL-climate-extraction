# ============================================================
# 0) Libraries
# ============================================================
library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(lubridate)
library(zoo)     # for rollapply
library(ggplot2)

# ============================================================
# 0.5) Verbosity (set FALSE when sourcing from Rmd)
# ============================================================
verbose <- FALSE

# ============================================================
# 1) Directory
# ============================================================
town_dir <- "/mnt/ceph/erichs/git/WWL-climate-extraction/data/ID_data/ID/"
stopifnot(dir.exists(town_dir))

# ============================================================
# 2) Index files + parse filename metadata
# ============================================================
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

if (verbose) print(file_index)

# ============================================================
# 3) Read into list (optional but handy)
# ============================================================
town_data <- file_index %>%
  mutate(obj_name = paste(town, buffer, source, sep = "__")) %>%
  select(obj_name, file) %>%
  tibble::deframe() %>%
  imap(\(f, nm) read_csv(f, show_col_types = FALSE))

inventory <- imap_dfr(town_data, \(df, nm) tibble(object = nm, n = nrow(df), p = ncol(df))) %>%
  arrange(object)
if (verbose) print(inventory, n = Inf)

# ============================================================
# 4) Combined table (this is what you’ll use for metrics)
# ============================================================
town_tbl <- file_index %>%
  mutate(data = map(file, \(f) read_csv(f, show_col_types = FALSE))) %>%
  select(town, buffer, source, data) %>%
  unnest(data)

# ============================================================
# 5) Christine-style time fields from system.index (KEEP THIS ONLY)
# ============================================================
town_tbl <- town_tbl %>%
  mutate(
    date_str = str_extract(`system.index`, "^\\d{8}"),
    date     = as.Date(date_str, format = "%Y%m%d"),
    year     = year(date),
    month    = month(date),
    wateryear = if_else(month >= 10, year + 1L, year),
    doy      = yday(date),
    season4  = case_when(
      month %in% 3:5  ~ "MAM",
      month %in% 6:8  ~ "JJA",
      month %in% 9:11 ~ "SON",
      TRUE            ~ "DJF"
    )
  )

stopifnot(!all(is.na(town_tbl$date)))  # sanity check that dates parsed

library(dplyr)
library(lubridate)

town_tbl <- town_tbl %>%
  mutate(
    yday = yday(date),
    mday = mday(date),
    # Christine-style water-year day (Oct 1 = 1)
    wyday = ifelse(
      year(date) == wateryear & leap_year(year(date)),
      yday + 92,
      ifelse(
        year(date) == wateryear,
        yday + 91,
        ifelse(leap_year(wateryear - 1), yday - 274, yday - 273)
      )
    )
  )

# ============================================================
# 6) Fire metrics from GridMET (bi / erc / fm100 / fm1000)
#    - This follows Christine’s pattern: daily flag -> annual summaries
# ============================================================

# Helper: max consecutive run of TRUEs
max_run <- function(x) {
  x <- as.integer(x)
  if (all(is.na(x))) return(NA_integer_)
  x[is.na(x)] <- 0L
  r <- rle(x)
  if (!any(r$values == 1L)) return(0L)
  max(r$lengths[r$values == 1L])
}

# Helper: first day-of-year where you get k consecutive "high" days
first_sustained_doy <- function(high_vec, doy_vec, k = 3) {
  if (length(high_vec) < k) return(NA_integer_)
  hv <- as.integer(high_vec)
  hv[is.na(hv)] <- 0L
  roll <- zoo::rollapply(hv, width = k, FUN = sum, align = "left", fill = NA)
  idx <- which(roll == k)
  if (length(idx) == 0) return(NA_integer_)
  doy_vec[min(idx)]
}

# Choose thresholds (edit as needed)
thr_bi  <- 50   # Burning Index: example threshold
thr_erc <- 60   # Energy Release Component: example threshold

fire_metrics_gridmet <- town_tbl %>%
  filter(source == "gridmet") %>%
  group_by(town, buffer, source, wateryear) %>%
  summarise(
    n_days = n(),
    
    # --- BI-based "high fire day" ---
    high_bi_days = sum(bi >= thr_bi, na.rm = TRUE),
    high_bi_days_summer = sum(bi >= thr_bi & season4 == "JJA", na.rm = TRUE),
    max_consec_high_bi = max_run(bi >= thr_bi),
    fire_start_doy_bi = first_sustained_doy(bi >= thr_bi, doy, k = 3),
    
    # --- ERC-based "high fire day" (optional) ---
    high_erc_days = sum(erc >= thr_erc, na.rm = TRUE),
    high_erc_days_summer = sum(erc >= thr_erc & season4 == "JJA", na.rm = TRUE),
    max_consec_high_erc = max_run(erc >= thr_erc),
    fire_start_doy_erc = first_sustained_doy(erc >= thr_erc, doy, k = 3),
    
    # --- Fuel moisture summaries (helpful context) ---
    fm100_mean = mean(fm100, na.rm = TRUE),
    fm1000_mean = mean(fm1000, na.rm = TRUE),
    
    .groups = "drop"
  )

print(fire_metrics_gridmet)

# ============================================================
# 7) Quick completeness checks (so you trust what’s being used)
# ============================================================
town_tbl %>%
  group_by(source) %>%
  summarise(
    n = n(),
    bi_nonmissing  = sum(!is.na(bi)),
    erc_nonmissing = sum(!is.na(erc)),
    fm100_nonmissing = sum(!is.na(fm100)),
    .groups = "drop"
  ) %>% print()



####################fire metrics calc

library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(zoo)

# Christine helper (same idea)
quibble <- function(x, q = c(0, 0.9, 0.95), dropNA = TRUE) {
  tibble(x = quantile(x, q, na.rm = dropNA), q = q)
}

# --- 1) Restrict to GridMET + choose WY range like Christine ---
wateryearselect <- 1980:2025

df <- town_tbl %>%
  filter(source == "gridmet", wateryear %in% wateryearselect) %>%
  select(town, buffer, source, `system.index`, date, year, month, wateryear, yday, wyday, mday, season4,
         bi, erc, fm100, fm1000, pr, tmmn, tmmx, vpd, vs)

# --- 2) Optional smoothing like Christine (5-day centered) ---
df5day <- df
df5day[, c("bi","erc")] <- rollmean(df[, c("bi","erc")], 5, align = "center", na.pad = TRUE)

# --- 3) Build daily 1996–2025 climatology quantiles for BI & ERC ---
daily9625_firequants <- df5day %>%
  filter(wateryear %in% 1996:2025) %>%
  group_by(month, mday) %>%
  reframe(across(c(bi, erc), ~ quibble(.x, q = c(0.9, 0.95), dropNA = TRUE))) %>%
  pivot_longer(cols = c(bi, erc), names_to = "var", values_to = "qb")

daily9625_firequants <- daily9625_firequants %>%
  mutate(
    quantile = qb$q,
    value = as.numeric(qb$x),
    name = paste0("dailyclimo9625_", var, "_q", quantile)
  ) %>%
  select(month, mday, name, value)

# attach daily thresholds (q90 and q95) back to all days
df_fireclimo <- df5day %>%
  left_join(daily9625_firequants, by = c("month","mday")) %>%
  pivot_wider(names_from = name, values_from = value)

# --- 4) Define “high fire day” flags (Christine-style: relative to daily climo) ---
df_fireclimo <- df_fireclimo %>%
  mutate(
    high_bi90  = bi  > dailyclimo9625_bi_q0.9,
    high_erc90 = erc > dailyclimo9625_erc_q0.9,
    high_bi95  = bi  > dailyclimo9625_bi_q0.95,
    high_erc95 = erc > dailyclimo9625_erc_q0.95,
    
    # severity = exceedance above threshold (0 if below)
    sev_bi90  = pmax(bi  - dailyclimo9625_bi_q0.9,  0, na.rm = TRUE),
    sev_erc90 = pmax(erc - dailyclimo9625_erc_q0.9, 0, na.rm = TRUE)
  )


first_sustained_wyday <- function(high_vec, wyday_vec, k = 3) {
  hv <- as.integer(high_vec); hv[is.na(hv)] <- 0L
  if (length(hv) < k) return(NA_integer_)
  roll <- zoo::rollapply(hv, width = k, FUN = sum, align = "left", fill = NA)
  idx <- which(roll == k)
  if (length(idx) == 0) NA_integer_ else wyday_vec[min(idx)]
}

# --- 5) Aggregate to WY metrics per town/buffer (this is your “annual metrics” table) ---
annualfire <- df_fireclimo %>%
  group_by(town, buffer, wateryear) %>%
  summarise(
    # Q1: overall frequency
    firehighdays_bi90  = sum(high_bi90,  na.rm = TRUE),
    firehighdays_erc90 = sum(high_erc90, na.rm = TRUE),
    
    # optional: more extreme frequency
    firehighdays_bi95  = sum(high_bi95,  na.rm = TRUE),
    firehighdays_erc95 = sum(high_erc95, na.rm = TRUE),
    
    # optional: severity (proxy for “worse years”)
    fireseverity_bi90  = sum(sev_bi90,  na.rm = TRUE),
    fireseverity_erc90 = sum(sev_erc90, na.rm = TRUE),
    
    # Q2: summer frequency
    firesummerhighdays_bi90  = sum(high_bi90  & season4 == "JJA", na.rm = TRUE),
    firesummerhighdays_erc90 = sum(high_erc90 & season4 == "JJA", na.rm = TRUE),
    
    # Q4: season start timing (sustained run)
    firestart_wyday_bi  = first_sustained_wyday(high_bi90,  wyday, k = 3),
    firestart_wyday_erc = first_sustained_wyday(high_erc90, wyday, k = 3),
    
    # Q5: duration
    firemaxrun_bi  = max_run(high_bi90),
    firemaxrun_erc = max_run(high_erc90),
    
    # optional: season length (first to last high day)
    fireseasonlength_bi = {
      idx <- which(high_bi90 %in% TRUE)
      if (length(idx) < 2) NA_integer_ else (max(wyday[idx]) - min(wyday[idx]))
    },
    fireseasonlength_erc = {
      idx <- which(high_erc90 %in% TRUE)
      if (length(idx) < 2) NA_integer_ else (max(wyday[idx]) - min(wyday[idx]))
    },
    
    # Q3 proxy: “large fire potential” (not area burned)
    firepeak_bi  = max(bi,  na.rm = TRUE),
    firepeak_erc = max(erc, na.rm = TRUE),
    firep95_erc  = as.numeric(quantile(erc, 0.95, na.rm = TRUE)),
    
    .groups = "drop"
  )

print(annualfire)


#######augment to run for k=1  (first high day) or k=3 (first three days that are high)

library(dplyr)
library(zoo)

annualfire2 <- df_fireclimo %>%   # this is the daily table you built right before annualfire
  group_by(town, buffer, wateryear) %>%
  summarise(
    # existing high-day definitions
    high_bi90  = list(high_bi90),
    high_erc90 = list(high_erc90),
    wyday      = list(wyday),
    .groups="drop"
  ) %>%
  mutate(
    firestart_any_bi90  = map2_int(high_bi90,  wyday, ~ { idx <- which(.x %in% TRUE); if(length(idx)==0) NA_integer_ else .y[min(idx)] }),
    firestart_any_erc90 = map2_int(high_erc90, wyday, ~ { idx <- which(.x %in% TRUE); if(length(idx)==0) NA_integer_ else .y[min(idx)] }),
    firestart_k3_bi90   = map2_int(high_bi90,  wyday, ~ first_sustained_wyday(.x, .y, k=3)),
    firestart_k2_erc90  = map2_int(high_erc90, wyday, ~ first_sustained_wyday(.x, .y, k=2)),
    firestart_k3_erc90  = map2_int(high_erc90, wyday, ~ first_sustained_wyday(.x, .y, k=3))
  ) %>%
  select(town, buffer, wateryear, starts_with("firestart"))

annualfire <- annualfire %>%
  left_join(annualfire2, by = c("town","buffer","wateryear"))


###now lets examine

library(dplyr)
library(tidyr)

annualfire_long <- annualfire %>%
  pivot_longer(
    cols = -c(town, buffer, wateryear),
    names_to = "metric",
    values_to = "value"
  ) %>%
  filter(!is.na(value))

# Quick peek
annualfire_long %>% count(metric, sort = TRUE)

##faceted plots

ggplot(
  annualfire_long %>% filter(metric == "firehighdays_erc90"),
  aes(x = wateryear, y = value, group = interaction(town, buffer))
) +
  geom_line(alpha = 0.4) +
  facet_grid(buffer ~ town, scales = "free_y") +
  labs(title = "firehighdays_erc90 by town/buffer", x = "Water year", y = "Days")


###sens slopes

library(trend)   # sens.slope, mk.test
fire_trends <- function(startyears = c(1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)) {
  
  metrics_to_trend <- c(
    "firehighdays_bi90","firehighdays_erc90",
    "firesummerhighdays_bi90","firesummerhighdays_erc90",
    "firestart_wyday_bi","firestart_wyday_erc",
    "firemaxrun_bi","firemaxrun_erc"
  )
  
  dfL <- annualfire_long %>% filter(metric %in% metrics_to_trend)
  
  out <- list()
  
  for (sy in startyears) {
    
    tmp <- dfL %>%
      filter(wateryear >= sy) %>%
      group_by(town, buffer, metric) %>%
      summarise(
        n_years = sum(!is.na(value)),
        sen_slope = {
          v <- na.omit(value)
          if (length(v) < 3) NA_real_ else trend::sens.slope(v)$estimates
        },
        p_value = {
          v <- na.omit(value)
          if (length(v) < 3) NA_real_ else trend::mk.test(v)$p.value
        },
        .groups = "drop"
      ) %>%
      mutate(TimePeriod = paste0(sy, "-2025"))
    
    out[[as.character(sy)]] <- tmp
  }
  
  bind_rows(out)
}

fire_slopes <- fire_trends()
fire_slopes %>% arrange(metric, p_value) %>% print(n = 50)


#distribution check

ggplot(
  annualfire_long %>% filter(metric %in% c("firehighdays_bi90","firehighdays_erc90")),
  aes(x = value)
) +
  geom_histogram(bins = 30) +
  facet_grid(metric ~ buffer, scales = "free_x") +
  labs(title = "Distribution of annual high fire days", x = "Days", y = "Count")


##plot dist check

library(tidyr)
library(ggplot2)

annualfire_long <- annualfire %>%
  pivot_longer(cols = -c(town, buffer, wateryear),
               names_to="metric", values_to="value")

ggplot(annualfire_long %>% filter(metric %in% c("firehighdays_bi90","firehighdays_erc90")),
       aes(x=wateryear, y=value)) +
  geom_line() +
  facet_grid(metric ~ town, scales="free_y") +
  labs(title="Annual high fire days (q90) by town", x="Water year", y="Days")


#buffer facet

ggplot(annualfire_long %>% filter(metric == "firehighdays_erc90"),
       aes(x=wateryear, y=value)) +
  geom_line() +
  facet_grid(buffer ~ town, scales="free_y") +
  labs(title="firehighdays_erc90 by town & buffer", x="Water year", y="Days")





early_years <- 1996:2005
late_years  <- 2016:2025

dfc <- annualfire %>%
  mutate(period = case_when(
    wateryear %in% early_years ~ "1996-2005",
    wateryear %in% late_years  ~ "2016-2025",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(period)) %>%
  select(town, buffer, period, firehighdays_erc90)

ggplot(dfc, aes(x=period, y=firehighdays_erc90)) +
  geom_boxplot() +
  facet_grid(buffer ~ town, scales="free_y") +
  labs(title="ERC90 high fire days: early vs late", x="", y="Days")




# ============================================================
# FIRE QUESTIONS: question-framed plots written to disk
# Output dir:
# /mnt/ceph/erichs/git/WWL-climate-extraction/plots/fire/
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# ---- output directory ----
fire_plot_dir <- "/mnt/ceph/erichs/git/WWL-climate-extraction/plots/fire"
dir.create(fire_plot_dir, recursive = TRUE, showWarnings = FALSE)

# ---- long format like Christine ----
annualfire_long <- annualfire %>%
  pivot_longer(
    cols = -c(town, buffer, wateryear),
    names_to = "metric",
    values_to = "value"
  )

# ---- helpers ----
`%||%` <- function(a, b) if (!is.null(a)) a else b

save_plot <- function(p, filename, width = 14, height = 8, dpi = 300) {
  ggsave(
    filename = file.path(fire_plot_dir, filename),
    plot = p,
    width = width, height = height, dpi = dpi
  )
}

ts_plot <- function(metric_name, title = NULL) {
  dfp <- annualfire_long %>% filter(metric == metric_name, !is.na(value))
  ggplot(dfp, aes(x = wateryear, y = value, group = interaction(town, buffer))) +
    geom_line(alpha = 0.6) +
    facet_grid(buffer ~ town, scales = "free_y") +
    labs(title = title %||% metric_name, x = "Water year", y = metric_name)
}

hist_plot <- function(metric_name, title = NULL) {
  dfp <- annualfire_long %>% filter(metric == metric_name, !is.na(value))
  ggplot(dfp, aes(x = value)) +
    geom_histogram(bins = 30) +
    facet_grid(buffer ~ town, scales = "free_x") +
    labs(title = title %||% paste("Distribution:", metric_name),
         x = metric_name, y = "Count")
}

period_compare_plot <- function(metric_name,
                                early_years = 1996:2005,
                                late_years  = 2016:2025,
                                title = NULL) {
  dfp <- annualfire %>%
    mutate(period = case_when(
      wateryear %in% early_years ~ "early",
      wateryear %in% late_years  ~ "late",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(period)) %>%
    select(town, buffer, period, value = all_of(metric_name))
  
  ggplot(dfp, aes(x = period, y = value)) +
    geom_boxplot() +
    facet_grid(buffer ~ town, scales = "free_y") +
    labs(title = title %||% paste(metric_name, "early vs late"),
         x = "", y = metric_name)
}

# ============================================================
# Q1) Overall frequency
# ============================================================
p <- ts_plot("firehighdays_bi90",  "Q1: Overall frequency (BI q90 high days)")
save_plot(p, "Q1_timeseries_firehighdays_bi90.png")

p <- ts_plot("firehighdays_erc90", "Q1: Overall frequency (ERC q90 high days)")
save_plot(p, "Q1_timeseries_firehighdays_erc90.png")

p <- hist_plot("firehighdays_erc90", "Q1: Distribution of annual ERC q90 high days")
save_plot(p, "Q1_distribution_firehighdays_erc90.png", width = 14, height = 9)

# ============================================================
# Q2) Summer frequency
# ============================================================
p <- ts_plot("firesummerhighdays_bi90",  "Q2: Summer frequency (BI q90 high days, JJA)")
save_plot(p, "Q2_timeseries_firesummerhighdays_bi90.png")

p <- ts_plot("firesummerhighdays_erc90", "Q2: Summer frequency (ERC q90 high days, JJA)")
save_plot(p, "Q2_timeseries_firesummerhighdays_erc90.png")

# ============================================================
# Q3) Magnitude proxy (NOT burned area)
# ============================================================
p <- ts_plot("firepeak_erc", "Q3 proxy: Annual peak ERC (large/impactful fire-weather potential)")
save_plot(p, "Q3proxy_timeseries_firepeak_erc.png")

p <- ts_plot("firep95_erc",  "Q3 proxy: Annual ERC 95th percentile (tail intensity)")
save_plot(p, "Q3proxy_timeseries_firep95_erc.png")

p <- ts_plot("firepeak_bi",  "Q3 proxy: Annual peak BI")
save_plot(p, "Q3proxy_timeseries_firepeak_bi.png")

# ============================================================
# Q4) Timing (start of season)
# ============================================================
p <- ts_plot("firestart_wyday_bi",  "Q4: Fire season start timing (BI q90, sustained run, wyday)")
save_plot(p, "Q4_timeseries_firestart_wyday_bi.png")

p <- ts_plot("firestart_wyday_erc", "Q4: Fire season start timing (ERC q90, sustained run, wyday)")
save_plot(p, "Q4_timeseries_firestart_wyday_erc.png")

p <- period_compare_plot("firestart_wyday_erc",
                         early_years = 1996:2005,
                         late_years  = 2016:2025,
                         title = "Q4: ERC start day (wyday) early vs late")
save_plot(p, "Q4_periodcompare_firestart_wyday_erc.png", width = 14, height = 9)

# ============================================================
# Q5) Duration
# ============================================================
p <- ts_plot("firemaxrun_bi",  "Q5: Max consecutive high-fire run length (BI q90)")
save_plot(p, "Q5_timeseries_firemaxrun_bi.png")

p <- ts_plot("firemaxrun_erc", "Q5: Max consecutive high-fire run length (ERC q90)")
save_plot(p, "Q5_timeseries_firemaxrun_erc.png")

p <- ts_plot("fireseasonlength_erc", "Q5 (optional): Season length (ERC q90: last−first high day)")
save_plot(p, "Q5_timeseries_fireseasonlength_erc.png")

message("Saved fire plots to: ", fire_plot_dir)


# ============================================================
# FIRE QUESTION → METRIC MAP (Option B: expanded)
# ============================================================
fire_map <- list(
  Q1 = c(
    "firehighdays_erc90","firehighdays_bi90",
    "firehighdays_erc95","firehighdays_bi95"
  ),
  Q2 = c(
    "firesummerhighdays_erc90","firesummerhighdays_bi90"
  ),
  Q3 = c(
    "firepeak_erc","firep95_erc","firepeak_bi",
    "fireseverity_erc90","fireseverity_bi90"
  ),
  Q4 = c(
    "firestart_wyday_erc","firestart_wyday_bi",
    "firestart_any_erc90","firestart_any_bi90",
    "firestart_k2_erc90","firestart_k3_erc90","firestart_k3_bi90"
  ),
  Q5 = c(
    "firemaxrun_erc","firemaxrun_bi",
    "fireseasonlength_erc","fireseasonlength_bi"
  )
)

# ============================================================
# FINAL: Christine-style "Question Outputs" for FIRE metrics
# Requires in memory:
#   - annualfire (wide annual metrics table)
#   - fire_map (list Q1..Q5 -> vector of metric names)
#   - forecast, trend loaded
# Also uses helper functions defined here.
# ============================================================

library(dplyr)
library(tidyr)
library(purrr)
library(trend)
library(forecast)
library(stringr)

# ---- choose output directory ----
out_dir <- "/mnt/ceph/erichs/git/WWL-climate-extraction/outputs/fire_questions"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- choose periods (match your earlier choices) ----
baseline_years <- 1996:2025
early_years    <- 1996:2005
late_years     <- 2016:2025

# ---- switch: use anomalies for ARIMA? (doesn't change slope, but can stabilize fits) ----
use_anomalies <- TRUE

# ---- long format ----
annualfire_long <- annualfire %>%
  pivot_longer(cols = -c(town, buffer, wateryear),
               names_to = "metric", values_to = "value") %>%
  arrange(town, buffer, metric, wateryear)

# ---- helpers ----
safe_sen_slope_per_decade <- function(x) {
  x <- as.numeric(x)
  if (sum(!is.na(x)) < 3) return(NA_real_)
  as.numeric(tryCatch(trend::sens.slope(x)$estimates * 10, error = function(e) NA_real_))
}

complete_years_vec <- function(years, vals) {
  yrs <- seq(min(years, na.rm=TRUE), max(years, na.rm=TRUE))
  tibble(wateryear = years, v = vals) %>%
    distinct(wateryear, .keep_all = TRUE) %>%
    right_join(tibble(wateryear = yrs), by = "wateryear") %>%
    arrange(wateryear) %>%
    pull(v)
}

arima_no_drift_qcat_diag <- function(x, n_sim = 2000, seed = 2026) {
  set.seed(seed)
  x <- as.numeric(x)
  x <- x[!is.na(x)]
  if (length(x) < 10) return(list(qcat = NA_character_,
                                  qbreaks = rep(NA_real_, 8),
                                  obs_slope = NA_real_,
                                  slopes_sim = NA_real_))
  
  fit <- tryCatch(
    forecast::auto.arima(x, allowdrift = FALSE, allowmean = TRUE),
    error = function(e) NULL
  )
  if (is.null(fit)) return(list(qcat = NA_character_,
                                qbreaks = rep(NA_real_, 8),
                                obs_slope = NA_real_,
                                slopes_sim = NA_real_))
  
  obs_slope <- safe_sen_slope_per_decade(x)
  
  slopes_sim <- numeric(n_sim)
  for (i in seq_len(n_sim)) {
    simx <- tryCatch(
      as.numeric(simulate(fit, nsim = length(x))),
      error = function(e) {
        r <- tryCatch(na.omit(residuals(fit)), error = function(e2) numeric(0))
        mu <- mean(x, na.rm = TRUE)
        if (length(r) > 5) mu + sample(r, length(x), replace = TRUE)
        else rnorm(length(x), mean = mu, sd = sd(x, na.rm = TRUE))
      }
    )
    slopes_sim[i] <- safe_sen_slope_per_decade(simx)
  }
  
  slopes_sim <- slopes_sim[is.finite(slopes_sim) & !is.na(slopes_sim)]
  if (length(slopes_sim) < 50 || !is.finite(obs_slope)) {
    return(list(qcat = NA_character_, qbreaks = rep(NA_real_, 8),
                obs_slope = obs_slope, slopes_sim = slopes_sim))
  }
  
  qbreaks <- as.numeric(quantile(slopes_sim, probs = seq(0, 1, length.out = 8), na.rm = TRUE, type = 7))
  if (any(diff(qbreaks) <= 0)) qbreaks <- seq(min(slopes_sim), max(slopes_sim), length.out = 8)
  
  labels <- paste0("Q", 1:7)
  qcat <- as.character(cut(obs_slope, breaks = qbreaks, labels = labels,
                           include.lowest = TRUE, right = FALSE))
  if (is.na(qcat) && obs_slope > max(qbreaks)) qcat <- "Q7"
  if (is.na(qcat) && obs_slope < min(qbreaks)) qcat <- "Q1"
  
  list(qcat = qcat, qbreaks = qbreaks, obs_slope = obs_slope, slopes_sim = slopes_sim)
}

# ---- compute anomalies + deltas + ARIMA diagnostics for ALL metrics ----
fire_qmap_df <- tibble(
  metric = unlist(fire_map),
  question = rep(names(fire_map), times = lengths(fire_map))
)

fire_anoms <- annualfire_long %>%
  group_by(town, buffer, metric) %>%
  mutate(
    baseline_mean = mean(value[wateryear %in% baseline_years], na.rm = TRUE),
    anom = value - baseline_mean
  ) %>%
  ungroup()

fire_delta <- annualfire_long %>%
  group_by(town, buffer, metric) %>%
  summarise(
    early_mean = mean(value[wateryear %in% early_years], na.rm = TRUE),
    late_mean  = mean(value[wateryear %in% late_years],  na.rm = TRUE),
    delta_late_minus_early = late_mean - early_mean,
    .groups = "drop"
  )

fire_arima_diag <- fire_anoms %>%
  group_by(town, buffer, metric) %>%
  summarise(
    years = list(wateryear),
    vals  = list(if (use_anomalies) anom else value),
    .groups="drop"
  ) %>%
  mutate(
    series_complete = map2(years, vals, ~ complete_years_vec(.x, .y)),
    diag = map(series_complete, ~ arima_no_drift_qcat_diag(.x, n_sim = 2000, seed = 2026)),
    sen_slope_per_decade = map_dbl(diag, "obs_slope"),
    qcat = map_chr(diag, "qcat"),
    p_upper = map_dbl(diag, ~ mean(.x$slopes_sim >= .x$obs_slope, na.rm=TRUE)),
    p_lower = map_dbl(diag, ~ mean(.x$slopes_sim <= .x$obs_slope, na.rm=TRUE)),
    p_two   = 2 * pmin(p_upper, p_lower)
  ) %>%
  select(town, buffer, metric, sen_slope_per_decade, qcat, p_upper, p_lower, p_two) %>%
  left_join(
    fire_anoms %>% distinct(town, buffer, metric, baseline_mean),
    by = c("town","buffer","metric")
  ) %>%
  left_join(fire_delta, by = c("town","buffer","metric")) %>%
  left_join(fire_qmap_df, by = "metric") %>%
  relocate(question, .before = metric)

# ---- write ONE combined file ----
all_path <- file.path(out_dir, "FIRE_ALL_questions_ARIMA_outputs.csv")
write.csv(fire_arima_diag, all_path, row.names = FALSE)

# ---- write one file per question (Christine-style buckets) ----
for (q in names(fire_map)) {
  dfq <- fire_arima_diag %>% filter(question == q)
  
  # Nice stable file naming
  q_path <- file.path(out_dir, paste0("FIRE_", q, "_ARIMA_outputs.csv"))
  write.csv(dfq, q_path, row.names = FALSE)
}

# ---- optional: write a compact summary like your table view ----
summary_path <- file.path(out_dir, "FIRE_metric_summary_by_Qcat.csv")
fire_arima_diag %>%
  count(question, metric, qcat) %>%
  arrange(question, metric, qcat) %>%
  write.csv(summary_path, row.names = FALSE)

message("Wrote FIRE question outputs to: ", out_dir)
message("Combined file: ", all_path)
message("Summary file: ", summary_path)




##CHECK the CSV files are correct

all_path <- file.path(out_dir, "FIRE_ALL_questions_ARIMA_outputs.csv")
df <- read.csv(all_path)

str(df)
summary(df[, c("sen_slope_per_decade","p_upper","p_lower","p_two")])

# how many missing ARIMA categories / slopes?
with(df, c(
  n = nrow(df),
  qcat_na = sum(is.na(qcat)),
  slope_na = sum(is.na(sen_slope_per_decade)),
  p_two_na = sum(is.na(p_two))
))



# ============================================================
# Tenure-based (person-specific) instrumentation Delta + PDI
# ============================================================

library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

# ---- Helper: safe numeric
to_num <- function(x) suppressWarnings(as.numeric(x))

# ---- Helper: compute a person-specific delta for one metric time series
# delta_method:
#   - "slope"  : OLS slope per decade over the person's lived window
#   - "diff"   : mean(last_k_years) - mean(first_k_years) within the lived window
compute_delta_one_series <- function(df_one, end_year, tenure_years,
                                     min_years = 8,
                                     delta_method = c("slope","diff"),
                                     k_years = 5) {
  delta_method <- match.arg(delta_method)
  
  if (is.na(tenure_years)) return(NA_real_)
  
  # Cap tenure to available record length
  df_one <- df_one %>% filter(!is.na(value))
  if (nrow(df_one) == 0) return(NA_real_)
  
  max_year <- min(end_year, max(df_one$year, na.rm = TRUE))
  # lived window start
  start_year <- max_year - floor(tenure_years) + 1
  
  df_w <- df_one %>%
    filter(year >= start_year, year <= max_year) %>%
    arrange(year)
  
  if (nrow(df_w) < min_years) return(NA_real_)
  
  if (delta_method == "slope") {
    # OLS slope (per year), convert to per decade
    fit <- try(stats::lm(value ~ year, data = df_w), silent = TRUE)
    if (inherits(fit, "try-error")) return(NA_real_)
    slope_per_year <- unname(stats::coef(fit)[["year"]])
    return(10 * slope_per_year)
  }
  
  # delta_method == "diff"
  k <- min(k_years, floor(nrow(df_w) / 2))
  if (k < 2) return(NA_real_)
  
  first_mean <- mean(df_w$value[1:k], na.rm = TRUE)
  last_mean  <- mean(df_w$value[(nrow(df_w)-k+1):nrow(df_w)], na.rm = TRUE)
  return(last_mean - first_mean)
}

# ---- Main: compute person-specific deltas + PDI
# survey: respondent-level survey df with town, tenure, and ID22–ID26
# annual_long: town-year-metric long table (must have columns: town, year, metric, value)
# survey_q_map: maps Q1–Q5 to survey columns (ID22–ID26)
# q_metric_map: maps Q1–Q5 to metric names in annual_long
compute_pdi_tenure <- function(survey,
                               annual_long,
                               survey_q_map,
                               q_metric_map,
                               town_col = "town",
                               id_col   = "record_id",
                               tenure_col = "tenure_years",
                               end_year = 2025,
                               delta_method = c("slope","diff"),
                               k_years = 5,
                               min_years = 8) {
  
  delta_method <- match.arg(delta_method)
  
  # 1) Long perception table (person × question)
  perception_long <- survey %>%
    mutate(
      .town = .data[[town_col]],
      .id   = .data[[id_col]],
      tenure_years = to_num(.data[[tenure_col]])
    ) %>%
    select(.id, .town, tenure_years, all_of(survey_q_map$survey_col)) %>%
    pivot_longer(
      cols = all_of(survey_q_map$survey_col),
      names_to = "survey_item",
      values_to = "resp"
    ) %>%
    left_join(survey_q_map, by = c("survey_item" = "survey_col")) %>%
    mutate(
      resp = to_num(resp),
      # reverse-code only if needed; 1–7 Likert => 8 - resp
      resp = if_else(direction == "reverse" & !is.na(resp), 8 - resp, resp)
    ) %>%
    select(.id, .town, tenure_years, question, resp)
  
  # 2) Person × question × metric rows (expand by q_metric_map)
  person_q_metric <- perception_long %>%
    left_join(q_metric_map, by = "question") %>%
    # q_metric_map should have a column "metric"
    filter(!is.na(metric))
  
  # 3) Join annual series for each town + metric, compute personalized delta
  # Make sure annual_long has the expected columns
  stopifnot(all(c("town","year","metric","value") %in% names(annual_long)))
  
  # Nest by person×question×metric so we can compute delta
  deltas <- person_q_metric %>%
    left_join(annual_long, by = c(".town" = "town", "metric" = "metric")) %>%
    group_by(.id, .town, tenure_years, question, resp, metric) %>%
    summarise(
      delta = compute_delta_one_series(
        df_one = tibble(year = year, value = value),
        end_year = end_year,
        tenure_years = unique(tenure_years),
        min_years = min_years,
        delta_method = delta_method,
        k_years = k_years
      ),
      .groups = "drop"
    )
  
  # 4) Aggregate metrics within question to a single instrumentation delta per person×question
  # (simple mean across linked metrics)
  instr_person_q <- deltas %>%
    group_by(.id, .town, tenure_years, question, resp) %>%
    summarise(
      instr_delta = mean(delta, na.rm = TRUE),
      n_metrics_used = sum(!is.na(delta)),
      .groups = "drop"
    )
  
  # 5) Compute PDI: standardize within each question (across people)
  # PDI = z(perception) - z(instrumentation)
  out <- instr_person_q %>%
    group_by(question) %>%
    mutate(
      z_resp  = as.numeric(scale(resp)),
      z_delta = as.numeric(scale(instr_delta)),
      PDI = z_resp - z_delta
    ) %>%
    ungroup()
  
  out
}
#END










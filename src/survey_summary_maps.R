
# ---- Survey Summary & Maps (Climate Perception) -----------------------------
# Author: (auto-generated template)
# Description: Summarize survey questions by location/demographics and make plots/maps.
# Requires: tidyverse, janitor, sf, tigris, tmap (or ggplot2), scales
# -----------------------------------------------------------------------------

# INSTALL (uncomment if needed)
# install.packages(c("tidyverse", "janitor", "sf", "tigris", "tmap", "scales"))

library(tidyverse)
library(janitor)
library(sf)
library(tigris)
library(tmap)
library(scales)

options(tigris_use_cache = TRUE)

# ---- 0) Paths ---------------------------------------------------------------
# Change this to where your CSV lives:
csv_path <- "../data/W2L_V1_REMOVED OUTLIERS.csv"
out_dir  <- "outputs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- 1) Load data (handle encoding) ----------------------------------------
# readr handles Latin-1 gracefully; adjust if needed.
dat <- readr::read_csv(csv_path, locale = readr::locale(encoding = "Latin1"))

# Quick peek
glimpse(dat)

# ---- 2) Identify key columns (robustly by patterns) ------------------------
# Using regex patterns so we don't rely on exact punctuation/smart quotes.
# Location
col_town   <- tidyselect::vars_select(names(dat), matches("(?i)\\bTown_Name\\b")) |> as.character()
col_zip    <- tidyselect::vars_select(names(dat), matches("(?i)\\bZip_code\\b")) |> as.character()
col_near   <- tidyselect::vars_select(names(dat), matches("(?i)\\bNearest_Town\\b")) |> as.character()

# Demographics
col_age    <- tidyselect::vars_select(names(dat), matches("^\\s*30\\.")) |> as.character()
col_edu    <- tidyselect::vars_select(names(dat), matches("^\\s*31\\.")) |> as.character()
col_occ    <- tidyselect::vars_select(names(dat), matches("^\\s*32\\.")) |> as.character()

# Perception questions (examples)
col_q6     <- tidyselect::vars_select(names(dat), matches("^\\s*6\\."))  |> as.character()
col_q7     <- tidyselect::vars_select(names(dat), matches("^\\s*7\\."))  |> as.character()
col_q7a    <- tidyselect::vars_select(names(dat), matches("^\\s*7a\\.")) |> as.character()
col_q8     <- tidyselect::vars_select(names(dat), matches("^\\s*8\\."))  |> as.character()
col_q9     <- tidyselect::vars_select(names(dat), matches("^\\s*9\\."))  |> as.character()

# If any are missing, warn:
message("Detected columns: ",
        "\n  Town: ", ifelse(length(col_town)>0, col_town, "NOT FOUND"),
        "\n  Zip: ",  ifelse(length(col_zip)>0,  col_zip, "NOT FOUND"),
        "\n  Nearest_Town: ", ifelse(length(col_near)>0, col_near, "NOT FOUND"),
        "\n  Age(30.): ", ifelse(length(col_age)>0, col_age, "NOT FOUND"),
        "\n  Education(31.): ", ifelse(length(col_edu)>0, col_edu, "NOT FOUND"),
        "\n  Occupation(32.): ", ifelse(length(col_occ)>0, col_occ, "NOT FOUND"),
        "\n  Q6(^6.): ", ifelse(length(col_q6)>0, col_q6, "NOT FOUND"),
        "\n  Q7(^7.): ", ifelse(length(col_q7)>0, col_q7, "NOT FOUND"),
        "\n  Q7a(^7a.): ", ifelse(length(col_q7a)>0, col_q7a, "NOT FOUND"),
        "\n  Q8(^8.): ", ifelse(length(col_q8)>0, col_q8, "NOT FOUND"),
        "\n  Q9(^9.): ", ifelse(length(col_q9)>0, col_q9, "NOT FOUND")
)

# ---- 3) Clean minimal set of columns ---------------------------------------
dat_min <- dat |>
  mutate(
    town   = if (length(col_town)) .data[[col_town]] else NA,
    zip    = if (length(col_zip))  .data[[col_zip]]  else NA,
    near_town = if (length(col_near)) .data[[col_near]] else NA,
    age    = if (length(col_age))  .data[[col_age]]  else NA,
    edu    = if (length(col_edu))  .data[[col_edu]]  else NA,
    occ    = if (length(col_occ))  .data[[col_occ]]  else NA,
    q6     = if (length(col_q6))   .data[[col_q6]]   else NA,
    q7     = if (length(col_q7))   .data[[col_q7]]   else NA,
    q7a    = if (length(col_q7a))  .data[[col_q7a]]  else NA,
    q8     = if (length(col_q8))   .data[[col_q8]]   else NA,
    q9     = if (length(col_q9))   .data[[col_q9]]   else NA
  ) |>
  mutate(
    across(everything(), \(x) if (is.character(x)) trimws(x) else x),
    zip = stringr::str_extract(as.character(zip), "\\d{5}") # normalize zip
  )

# ---- 4) Basic summaries -----------------------------------------------------
# Counts by town / nearest_town / zip
counts_by_loc <- dat_min |>
  summarise(
    n = n(),
    .by = c(town, near_town, zip)
  ) |>
  arrange(desc(n))

readr::write_csv(counts_by_loc, file.path(out_dir, "counts_by_location.csv"))

# Age distribution
age_tab <- dat_min |>
  count(age, sort = TRUE)
readr::write_csv(age_tab, file.path(out_dir, "age_distribution.csv"))

# Education distribution
edu_tab <- dat_min |>
  count(edu, sort = TRUE)
readr::write_csv(edu_tab, file.path(out_dir, "education_distribution.csv"))

# Q6/Q7 recode (example categories)
# Adjust these patterns if your survey uses different wording.
recode_temp_dir <- function(x) {
  dplyr::case_when(
    str_detect(str_to_lower(x), "hot")   ~ "Hotter",
    str_detect(str_to_lower(x), "cool")  ~ "Cooler",
    str_detect(str_to_lower(x), "same|no change|about the same") ~ "No change",
    TRUE ~ NA_character_
  )
}

dat_min <- dat_min |>
  mutate(
    q6_dir = recode_temp_dir(q6),
    q7_dir = recode_temp_dir(q7)
  )

# Share "Hotter" by town / zip / age group / education
share_fun <- \(df, grp) {
  df |>
    filter(!is.na(q6_dir)) |>
    summarise(
      n = dplyr::n(),
      hotter = sum(q6_dir == "Hotter", na.rm = TRUE),
      share_hotter = hotter / n,
      .by = {{grp}}
    ) |>
    arrange(desc(share_hotter))
}

share_by_town <- share_fun(dat_min, town)
share_by_zip  <- share_fun(dat_min, zip)
share_by_age  <- share_fun(dat_min, age)
share_by_edu  <- share_fun(dat_min, edu)

readr::write_csv(share_by_town, file.path(out_dir, "share_hotter_by_town.csv"))
readr::write_csv(share_by_zip,  file.path(out_dir, "share_hotter_by_zip.csv"))
readr::write_csv(share_by_age,  file.path(out_dir, "share_hotter_by_age.csv"))
readr::write_csv(share_by_edu,  file.path(out_dir, "share_hotter_by_edu.csv"))

# ---- 5) Plots ---------------------------------------------------------------
# Bar: share hotter by town (top 15)
p_town <- share_by_town |>
  slice_max(order_by = n, n = 15) |>
  ggplot(aes(x = reorder(town, share_hotter), y = share_hotter)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Town", y = "Share 'Hotter' (Q6)",
       title = "Perception that it has gotten hotter (by town)") +
  theme_minimal()

ggsave(file.path(out_dir, "plot_share_hotter_by_town.png"), p_town, width = 8, height = 6, dpi = 150)

# Bar: share hotter by education (top categories)
p_edu <- share_by_edu |>
  slice_max(order_by = n, n = 12) |>
  ggplot(aes(x = reorder(edu, share_hotter), y = share_hotter)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Education", y = "Share 'Hotter' (Q6)",
       title = "Perception that it has gotten hotter (by education)") +
  theme_minimal()

ggsave(file.path(out_dir, "plot_share_hotter_by_education.png"), p_edu, width = 8, height = 6, dpi = 150)

# ---- 6) Map by ZIP (ZCTA) ---------------------------------------------------
# Only if zip column exists:
if (length(col_zip)) {
  # Summarize by ZIP
  zsum <- share_by_zip |>
    filter(!is.na(zip)) |>
    mutate(zip = stringr::str_pad(zip, 5, pad = "0"))

  # Get ZCTAs (2020)
  zctas <- tigris::zctas(year = 2020, cb = TRUE) |>
    st_transform(3857) # web mercator for quick plotting

  zctas_join <- zctas |>
    mutate(ZCTA5CE10 = if ("ZCTA5CE10" %in% names(.)) ZCTA5CE10 else GEOID10,
           ZCTA5CE20 = if ("ZCTA5CE20" %in% names(.)) ZCTA5CE20 else GEOID20) |>
    mutate(zip = dplyr::coalesce(ZCTA5CE20, ZCTA5CE10)) |>
    select(zip, geometry) |>
    left_join(zsum, by = "zip")

  # Quick choropleth with tmap
  tm <- tmap_mode("plot")
  map_zip <- tm_shape(zctas_join) +
    tm_polygons(col = "share_hotter",
                palette = "YlOrRd",
                title = "Share 'Hotter' (Q6)",
                textNA = "No data",
                legend.format = list(fun = function(x) percent(x, accuracy = 1))) +
    tm_layout(title = "Perceived warming by ZIP (Q6)")

  tmap_save(map_zip, filename = file.path(out_dir, "map_share_hotter_by_zip.png"), width = 10, height = 7, dpi = 150)
}

# ---- 7) Multi-question long format (example) --------------------------------
# If you want to summarize multiple perception items at once, pivot longer:
question_cols <- names(dat)[stringr::str_detect(names(dat), "^(\\s*6\\.|\\s*7\\.|\\s*8\\.|\\s*9\\.|eff_mitig_|climate|warming|heat|hot)") ]

dat_long <- dat |>
  select(all_of(c(col_town, col_zip, col_near, col_age, col_edu, question_cols))) |>
  pivot_longer(cols = all_of(question_cols), names_to = "question", values_to = "response")

# Example: percent "Hotter" per question by education
sum_long <- dat_long |>
  mutate(
    edu = if (length(col_edu)) .data[[col_edu]] else NA_character_,
    town = if (length(col_town)) .data[[col_town]] else NA_character_
  ) |>
  group_by(question, edu) |>
  summarise(
    n = n(),
    pct_hotter = mean(str_detect(str_to_lower(response), "hot"), na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(question, desc(pct_hotter))

readr::write_csv(sum_long, file.path(out_dir, "multi_question_by_education.csv"))

# Plot: top 10 questions with highest "hotter" share for a chosen education group
topq <- sum_long |>
  filter(!is.na(edu)) |>
  group_by(question) |>
  summarise(mean_hotter = mean(pct_hotter, na.rm = TRUE), .groups = "drop") |>
  slice_max(order_by = mean_hotter, n = 10) |>
  pull(question)

p_long <- sum_long |>
  filter(question %in% topq) |>
  ggplot(aes(x = edu, y = pct_hotter)) +
  geom_col() +
  facet_wrap(~ question, scales = "free_y") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Education", y = "Share 'Hotter'", title = "Top questions by 'Hotter' share (by education)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_dir, "plot_multi_question_by_education.png"), p_long, width = 12, height = 8, dpi = 150)

# ---- 8) Save a clean CSV subset for downstream use --------------------------
dat_export <- dat_min |>
  select(town, near_town, zip, age, edu, occ, q6, q7, q7a, q8, q9, q6_dir, q7_dir)

readr::write_csv(dat_export, file.path(out_dir, "clean_subset.csv"))

message("Done. Outputs written to: ", normalizePath(out_dir))

# ---- Libraries ----
library(readr)
library(dplyr)
library(stringr)
library(purrr)

# ---- 1) Point to your folder ----
fire_dir <- "/mnt/ceph/erichs/owncloud/W2L/Fire"   # <-- change if yours differs

# sanity check
stopifnot(dir.exists(fire_dir))

# ---- 2) Find CSVs ----
csv_files <- list.files(
  path = fire_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)

if (length(csv_files) == 0) stop("No CSV files found in: ", fire_dir)

print(csv_files)

# ---- 3) Read all CSVs into a named list ----
# names will be file basenames without .csv
fire_data <- csv_files |>
  set_names(nm = tools::file_path_sans_ext(basename(csv_files))) |>
  map(\(f) read_csv(f, show_col_types = FALSE))

# Quick peek: row/col counts
fire_data |>
  imap(\(df, nm) tibble(dataset = nm, n = nrow(df), p = ncol(df))) |>
  bind_rows() |>
  arrange(dataset) |>
  print(n = Inf)

# ---- 4) Convenience objects (if present) ----
# These names match what you're showing in the folder screenshot
if ("CEMS_FWI_change_counties_ID" %in% names(fire_data)) {
  fwi_counties_id <- fire_data[["CEMS_FWI_change_counties_ID"]]
}
if ("CEMS_FWI_change_counties_N" %in% names(fire_data)) {
  fwi_counties_nv <- fire_data[["CEMS_FWI_change_counties_N"]]
}
if ("CEMS_FWI_change_huc10_ID" %in% names(fire_data)) {
  fwi_huc10_id <- fire_data[["CEMS_FWI_change_huc10_ID"]]
}
if ("CEMS_FWI_change_huc10_N" %in% names(fire_data)) {
  fwi_huc10_nv <- fire_data[["CEMS_FWI_change_huc10_N"]]
}
if ("CEMS_FWI_change_huc12_ID" %in% names(fire_data)) {
  fwi_huc12_id <- fire_data[["CEMS_FWI_change_huc12_ID"]]
}
if ("CEMS_FWI_change_huc12_N" %in% names(fire_data)) {
  fwi_huc12_nv <- fire_data[["CEMS_FWI_change_huc12_N"]]
}

# ---- 5) Inspect one dataset's columns quickly ----
# Example:
if (exists("fwi_counties_id")) {
  glimpse(fwi_counties_id)
  names(fwi_counties_id)
}
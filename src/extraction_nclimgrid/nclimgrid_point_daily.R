library(httr)
library(ncdf4)
library(dplyr)
library(lubridate)
library(purrr)

# Get daily NClimGrid at a single point, for a date range
nclimgrid_point_daily <- function(lat, lon,
                                  start_date, end_date,
                                  var = c("tmax", "tmin", "tavg", "prcp"),
                                  accept = c("netCDF", "netCDF4")) {
  var    <- match.arg(var)
  accept <- match.arg(accept)
  
  start_date <- as.Date(start_date)
  end_date   <- as.Date(end_date)
  if (end_date < start_date) stop("end_date must be >= start_date")
  
  # monthly sequence covering [start_date, end_date]
  month_seq <- seq(floor_date(start_date, "month"),
                   floor_date(end_date, "month"),
                   by = "1 month")
  
  out_list <- vector("list", length(month_seq))
  
  for (i in seq_along(month_seq)) {
    ym   <- month_seq[i]
    yr   <- year(ym)
    ym_s <- format(ym, "%Y%m")
    
    # per-month time window, clipped to global [start_date, end_date]
    t_start <- max(start_date, ym)
    t_end   <- min(end_date, ym + months(1) - days(1))
    if (t_start > t_end) next
    
    base_url <- paste0(
      "https://www.ncei.noaa.gov/thredds/ncss/grid/nclimgrid-daily/",
      yr, "/ncdd-", ym_s, "-grd-scaled.nc"
    )
    
    query <- list(
      var        = var,
      latitude   = lat,
      longitude  = lon,
      time_start = paste0(t_start, "T00:00:00Z"),
      time_end   = paste0(t_end,   "T00:00:00Z"),
      timeStride = 1,
      accept     = accept
    )
    
    url <- modify_url(base_url, query = query)
    message("Requesting: ", url)
    
    tmp  <- tempfile(fileext = ".nc")
    resp <- GET(url, write_disk(tmp, overwrite = TRUE))
    
    if (http_error(resp)) {
      warning("HTTP error for ", ym_s, ": ", status_code(resp))
      next
    }
    
    nc <- nc_open(tmp)
    
    # --- time handling (case-insensitive "Day"/"days"/"Hour"/"hours") ---
    time_vals  <- ncvar_get(nc, "time")
    time_units <- ncatt_get(nc, "time", "units")$value
    # e.g. "Day since 1951-01-01 00:00:00"
    
    parts <- strsplit(time_units, " since ")[[1]]
    if (length(parts) != 2) {
      nc_close(nc)
      stop("Unrecognized time units: ", time_units)
    }
    
    unit_str <- tolower(trimws(parts[1]))  # "day" / "days" / "hour(s)"
    origin   <- as.POSIXct(parts[2], tz = "UTC")
    
    if (grepl("^day", unit_str)) {
      time_posix <- origin + ddays(time_vals)
    } else if (grepl("^hour", unit_str)) {
      time_posix <- origin + dhours(time_vals)
    } else {
      nc_close(nc)
      stop("Unhandled time unit: ", parts[1], " in '", time_units, "'")
    }
    # --------------------------------------------------------------------
    
    vals <- ncvar_get(nc, var)
    nc_close(nc)
    
    vals <- as.vector(vals)
    
    df_month <- tibble(
      date  = as.Date(time_posix),
      value = vals
    )
    
    out_list[[i]] <- df_month
  }
  
  out <- bind_rows(out_list) %>%
    filter(date >= start_date, date <= end_date) %>%
    arrange(date)
  
  attr(out, "lat") <- lat
  attr(out, "lon") <- lon
  attr(out, "var") <- var
  
  out
}



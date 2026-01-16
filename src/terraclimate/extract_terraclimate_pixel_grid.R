extract_terraclimate_pixel_grid <- function(
    variables,
    start_year = 1958,
    lon_range = c(-124.848974, -66.93457),
    lat_range = c(24.396308, 49.384358),
    out_dir = "./data/pixel_data_1958_present",
    return = c("rast", "none"),        # "rast" returns SpatRaster(s); "none" just writes files
    write_files = TRUE,
    chunk = c("year", "month"),        # how to write output
    format = c("wide", "long"),        # wide = 1 row per cell per chunk; long = exploded rows
    na_rm = FALSE,
    overwrite = FALSE,
    combine_rasters = FALSE            # if return="rast" and multiple vars, stack them
) {
  # Packages
  suppressPackageStartupMessages({
    library(ncdf4)
    library(terra)
    library(dplyr)
    library(tidyr)
  })
  
  return <- match.arg(return)
  chunk  <- match.arg(chunk)
  format <- match.arg(format)
  
  if (!dir.exists(out_dir) && write_files) dir.create(out_dir, recursive = TRUE)
  
  if (start_year < 1958) stop("start_year must be >= 1958 for TerraClimate aggregated file.")
  
  results_rast <- list()
  
  for (var in variables) {
    message("Processing variable: ", var)
    
    baseurlagg <- paste0(
      "http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_terraclimate_",
      var, "_1958_CurrentYear_GLOBE.nc"
    )
    
    nc <- nc_open(baseurlagg)
    on.exit(try(nc_close(nc), silent = TRUE), add = TRUE)
    
    lon <- ncvar_get(nc, "lon")
    lat <- ncvar_get(nc, "lat")
    
    lon_idx <- which(lon >= min(lon_range) & lon <= max(lon_range))
    lat_idx <- which(lat >= min(lat_range) & lat <= max(lat_range))
    
    if (length(lon_idx) == 0 || length(lat_idx) == 0) {
      stop("No lon/lat indices found for the provided bounding box.")
    }
    
    start_month_index <- (start_year - 1958) * 12 + 1
    start <- c(lon_idx[1], lat_idx[1], start_month_index)
    count <- c(length(lon_idx), length(lat_idx), -1)
    
    # data dims: [lon, lat, time]
    data <- ncvar_get(nc, varid = var, start = start, count = count)
    nc_close(nc)
    
    lon_sub <- lon[lon_idx]
    lat_sub <- lat[lat_idx]
    
    nlon  <- length(lon_sub)
    nlat  <- length(lat_sub)
    ntime <- dim(data)[3]
    
    # Convert to [row(lat), col(lon), time] for terra::rast(array)
    arr <- aperm(data, c(2, 1, 3))  # lat, lon, time
    
    # terra expects row 1 = "top" (max latitude). If lat is ascending (south->north), flip rows.
    if (lat_sub[1] < lat_sub[length(lat_sub)]) {
      arr <- arr[rev(seq_len(nlat)), , ]
    }
    
    r <- terra::rast(
      arr,
      extent = terra::ext(min(lon_sub), max(lon_sub), min(lat_sub), max(lat_sub)),
      crs = "EPSG:4326"
    )
    
    # Create time labels
    t_index <- 0:(ntime - 1)
    year_vec  <- start_year + (t_index %/% 12)
    month_vec <- (t_index %% 12) + 1
    
    names(r) <- sprintf("%s_%04d_%02d", var, year_vec, month_vec)
    
    # Optionally write files in chunks (recommended for pixel-level)
    if (write_files) {
      if (chunk == "year") {
        years <- sort(unique(year_vec))
        for (yy in years) {
          idx <- which(year_vec == yy)
          r_yy <- r[[idx]]
          
          # Wide: one row per cell, 12 month columns (or however many months exist in that year)
          df <- terra::as.data.frame(r_yy, xy = TRUE, cells = TRUE, na.rm = na_rm)
          
          if (format == "wide") {
            # rename layer columns to var_MM (so each year file has 12 consistent month columns)
            mm <- month_vec[idx]
            new_names <- sprintf("%s_%02d", var, mm)
            names(df)[-(1:3)] <- new_names
            df <- df %>%
              mutate(year = yy) %>%
              select(cell, x, y, year, everything())
            
          } else {
            # Long: one row per cell per month (BIG)
            df <- df %>%
              pivot_longer(
                cols = -(cell:x:y),
                names_to = "band",
                values_to = var
              ) %>%
              tidyr::extract(
                band,
                into = c("v", "year", "month"),
                regex = "^(.+)_([0-9]{4})_([0-9]{2})$",
                remove = TRUE,
                convert = TRUE
              ) %>%
              select(cell, x, y, year, month, all_of(var))
          }
          
          out_file <- file.path(out_dir, sprintf("%s_pixels_%d.csv", var, yy))
          if (!file.exists(out_file) || overwrite) {
            write.csv(df, out_file, row.names = FALSE)
          }
        }
      }
      
      if (chunk == "month") {
        for (i in seq_len(ntime)) {
          yy <- year_vec[i]
          mm <- month_vec[i]
          
          df <- terra::as.data.frame(r[[i]], xy = TRUE, cells = TRUE, na.rm = na_rm)
          names(df)[4] <- var
          
          df <- df %>% mutate(year = yy, month = mm) %>%
            select(cell, x, y, year, month, all_of(var))
          
          out_file <- file.path(out_dir, sprintf("%s_pixels_%d_%02d.csv", var, yy, mm))
          if (!file.exists(out_file) || overwrite) {
            write.csv(df, out_file, row.names = FALSE)
          }
        }
      }
    }
    
    results_rast[[var]] <- r
  }
  
  # Return options
  if (return == "none") return(invisible(TRUE))
  
  if (length(results_rast) == 1) return(results_rast[[1]])
  
  if (combine_rasters) {
    return(do.call(c, results_rast))  # stack all variables into one SpatRaster
  } else {
    return(results_rast)              # named list of SpatRaster, one per variable
  }
}


climate_vars <- c("aet","def","PDSI","pet","ppt","q","soil",
                  "srad","swe","tmax","tmin","vap","vpd","ws")

rasters <- extract_terraclimate_pixel_grid(
  variables    = climate_vars,
  start_year   = 1979,
  write_files  = FALSE,
  return       = "rast",
  combine_rasters = FALSE
)

# Example: rasters[["tmax"]]

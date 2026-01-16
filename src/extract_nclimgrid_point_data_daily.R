# NClimGrid analogue of extract_climate_point_data_daily()
# var_type should be one of: "tmax", "tmin", "tavg", "prcp"
extract_nclimgrid_point_data_daily <- function(var_type,
                                               start_date = "1951-01-01",
                                               end_date   = Sys.Date()) {
  
  library(dplyr)
  library(lubridate)

  # coordinates: same file you already use for gridMET
  df_coordinates <- read.csv("./data/coordinates/city_coordinates.csv")
  # expected columns: lon, lat, label
  
  all_climate_data <- data.frame()
  
  for (i in 1:nrow(df_coordinates)) {
    target_lon   <- df_coordinates$lon[i]
    target_lat   <- df_coordinates$lat[i]
    coord_label  <- df_coordinates$label[i]
    
    # get daily series for this point from NClimGrid
    point_df <- nclimgrid_point_daily(
      lat        = target_lat,
      lon        = target_lon,
      start_date = start_date,
      end_date   = end_date,
      var        = var_type
    )
    
    # structure to match your gridMET output
    climate_df <- point_df %>%
      mutate(
        coord = coord_label,
        day   = day(date),
        month = month(date),
        year  = year(date)
      )
    
    # rename "value" -> var_type (e.g., "tmax", "tmin", etc.)
    colnames(climate_df)[which(names(climate_df) == "value")] <- var_type
    
    all_climate_data <- bind_rows(all_climate_data, climate_df)
  }
  
  # Write CSV; change filename here if you want it 100% identical to gridMET
  out_path <- paste0("./data/nclimgrid_point_data/", var_type, "_nclimgrid.csv")
  write.csv(all_climate_data, file = out_path, row.names = FALSE)
  
  message("Wrote: ", out_path)
  
  invisible(all_climate_data)
}

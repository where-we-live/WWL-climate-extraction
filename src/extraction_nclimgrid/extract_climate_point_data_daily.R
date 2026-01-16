#Title: extract_climate_data.R
#Author: Erich Seamon
#Date: 08.12.2025
#Description: this is function which extracts climate data from the UI THREDDS climate server
#the inputs to the function are as follows: 
#
#df_coordinates - these are the point coordinates for the location to extract, in a data frame,
#with lat, lon, and a label.
#var_type - is the climate variable, which includes:

#--------------------------------

extract_climate_point_data_daily <- function(var_type) {
  
  # Load necessary libraries
  library(ncdf4)
  library(dplyr)
  library(lubridate)
  
  df_coordinates <- read.csv("./data/coordinates/city_coordinates.csv")
  
  # Example coordinates
  # df_coordinates <- data.frame(
  #   lon = c(-116.394, -116.647, -116.556, -116.77),
  #   lat = c(46.859, 46.614, 46.799, 46.737),
  #   label = c("Bovill", "Kendrick", "Deary", "Troy")
  # )
  
  # coordinates <- subset(df_coordinates, label == cityname)
  # coordinates <- coordinates[,c(1:2)]
  
  opendap_url <- paste("http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_met_", var_type, "_1979_CurrentYear_CONUS.nc#fillmismatch", sep="")
  
  # Initialize an empty data frame to store data for all coordinates
  all_climate_data <- data.frame()
  
  # Loop over each set of coordinates
  for (i in 1:nrow(df_coordinates)) {
    target_lon <- df_coordinates$lon[i]
    target_lat <- df_coordinates$lat[i]
    coord_label <- df_coordinates$label[i]  # Custom label for each point
    
    # Open the NetCDF file
    nc_data <- nc_open(opendap_url)
    
    # Identify variable name automatically
    var_name <- names(nc_data$var)[1]
    
    # Extract longitude, latitude, and time
    lon <- ncvar_get(nc_data, "lon")
    lat <- ncvar_get(nc_data, "lat")
    time <- ncvar_get(nc_data, "day")
    
    # Find the closest indices for target longitude and latitude
    lon_ind <- which.min(abs(lon - target_lon))
    lat_ind <- which.min(abs(lat - target_lat))
    
    # Set up start and count for data extraction
    start <- c(lon_ind, lat_ind, 1)
    count <- c(1, 1, -1)  # -1 retrieves all time steps
    
    # Extract data
    var_array <- ncvar_get(nc_data, var_name, start = start, count = count)
    nc_close(nc_data)  # Close the NetCDF file
    
    # Convert time to Date format
    start_date <- as.Date("1900-01-01")
    dates <- start_date + time
    
    # Create a data frame for the current coordinate point
    climate_df <- data.frame(date = dates, value = as.numeric(var_array), coord = coord_label)
    
    # Convert temperature from Kelvin to Celsius if var_type is temperature
    if (var_name == "daily_minimum_temperature"|var_name=="daily_maximum_temperature") {
      climate_df$value <- climate_df$value - 273.15
    }
    
    # Rename the value column based on the var_name
    colnames(climate_df)[which(names(climate_df) == "value")] <- var_type
    
    # Add day, month, and year columns
    climate_df <- climate_df %>%
      mutate(day = day(date), month = month(date), year = year(date))
    
    # Append the data to the main data frame
    all_climate_data <- bind_rows(all_climate_data, climate_df)
  }
  
  # Return the extracted climate data
  #return(all_climate_data)
  
  
  write.csv(all_climate_data, file= paste("./data/point_data/", var_type, ".csv", sep=""), row.names = FALSE)
            
}





# install.packages(c("sf", "tigris", "terra", "dplyr", "ggplot2"))
library(sf)
library(tigris)
library(terra)
library(dplyr)
library(ggplot2)

options(tigris_use_cache = TRUE)

# 1) Read an NClimGrid NetCDF raster (any daily/monthly file)
#   CHANGE THIS to the path where you have an NClimGrid file:
nclim_file <- "/mnt/ceph/erichs/git/WWL-climate-extraction/data/example_netcdf/ncdd-202401-grd-scaled.nc"

r <- terra::rast(nclim_file)  # terra raster

# 2) Get Idaho boundary and transform into the NClimGrid CRS
idaho <- states(cb = TRUE, year = 2023) %>%
  filter(NAME == "Idaho") %>%
  st_transform(crs(r))   # match NClimGrid projection

idaho_v <- vect(idaho)   # convert sf -> SpatVector for terra


# 3) Crop & mask the NClimGrid raster to Idaho
r_id <- crop(r, idaho_v)
r_id <- mask(r_id, idaho_v)

# 4) Turn the raster cells into polygons (each cell → 1 polygon)
grid_idaho <- as.polygons(r_id, dissolve = FALSE)
grid_idaho_sf <- st_as_sf(grid_idaho)

ggplot() +
  geom_sf(data = grid_idaho_sf, fill = NA, linewidth = 0.1) +
  geom_sf(data = idaho, fill = NA, color = "red", linewidth = 0.4) +
  ggtitle("Idaho with NClimGrid Cell Boundaries") +
  theme_minimal()




#### Preprocess netcdf of ERA5 altitudes to dataframe ####

#### Libraries ####

library(terra)
library(sf)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
# for loading xarray based .nc data
library(ncdf4)
library(ncdf4.helpers)
library(parallel)
library(ggcorrplot)
library(data.table)

source("src/functions.R")


#### Load Data ####

# TODO: Extend to all files
file <- "data/era5_elevation/ASurf_WFDE5_CRU_v2.1.nc"

# load weather station data (with wind data)
# met_eir <- readr::read_csv("data/met_eireann/final/met_eir.csv.gz")
met_eir <- readr::read_csv("data/met_eireann/final/met_eir_wind.csv.gz")
# take locations only to match to grid points
met_eir_sf <- met_eir %>% 
  distinct(name, county, province, lon, lat) %>% 
  sf::st_as_sf(coords = c("lon", "lat"))
st_crs(met_eir_sf) <- "WGS84"

# load netcdf, convert to dataframe
data <- load_netcdf(file) %>% 
  rename("alt" = ASurf) %>% 
  filter(!is.na(alt))


#### Triangulate to weather stations ####

data_sf <- data %>% 
  # only keep longitudes and latitudes for matching
  distinct(lon = longitude, lat = latitude) %>% 
  st_as_sf(coords = c("lon", "lat"))
st_crs(data_sf) <- "WGS84"

# match locations to closest points in grid
idx <- st_nearest_feature(met_eir_sf, data_sf)

# take these points, attach names of closest stations
data_sf <- data_sf %>% 
  cbind(st_coordinates(data_sf)) %>% 
  rename(longitude = X, latitude = Y) %>% 
  slice(idx) %>% 
  st_drop_geometry() %>% 
  cbind(st_drop_geometry(met_eir_sf))
  
# add labels to data
data_lbl <- left_join(data, data_sf, relationship = "many-to-many") %>% 
  # remove points not corresponding to sites
  filter(!is.na(name)) %>% 
  # remove approximate location, keep only date and name for each wind speed
  select(-c(longitude, latitude))

# join into other data
met_eir_alt <- left_join(met_eir, data_lbl)


#### Save ####

readr::write_csv(met_eir_alt, "data/met_eireann/final/met_eir_wind_alt.csv.gz")

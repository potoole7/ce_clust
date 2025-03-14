##### Preprocess netcdf of ERA5 wind speeds into dataframe ####

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

# ERA5 files downloaded from https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=overview
# 10m u and v wind speed and direction, area: [56, -11, 51, -4]
# file <- "data/download.nc"
files <- list.files("data", pattern = ".nc", full.names = TRUE)

# data <- load_netcdf(file)
# gc()
# # save for easier future access
# data.table::fwrite(data, file = "data/era5_windspeed.csv.gz")

# load from file
# data <- data.table::fread("data/era5_windspeed.csv.gz")
# gc()

# load weather station data
met_eir <- readr::read_csv("data/met_eireann/final/met_eir.csv.gz")
# add rain data for weather stations in six counties as well!
ceda <- readr::read_csv("data/ceda_open/ceda_rain.csv.gz")
met_eir <- bind_rows(met_eir, ceda)
# take locations only to match to grid points
met_eir_sf <- met_eir %>% 
  distinct(name, county, province, lon, lat) %>% 
  sf::st_as_sf(coords = c("lon", "lat"))
st_crs(met_eir_sf) <- "WGS84"

# function to find value corresponding to maximum density estimate for vector
max_density <- \(vec) {
  if (length(vec) == 1) return(vec)
  # return maximum kernel density estimate
  return(with(density(vec), x[which.max(y)]))
}

# function to perform preprocessing on single dataset
preprocess_netcdf <- \(file) {
  
  # load netcdf, convert to dataframe
  print(paste("Loading", file))
  data <- load_netcdf(file)
  gc()
  
  # take daily maximum windspeeds
  data <- data %>% 
    mutate(
      # calculate wind speed and direction from u and v components
      wind_speed = sqrt(u10 ** 2 + v10 ** 2),
      # wind_dir   = atan2(u10, v10) * (180 / pi), # convert to degrees
      wind_dir = (180 + (180 / pi) * atan2(u10, v10)) %% 360,
      # pull day from time
      date       = as_date(time)
    ) %>% 
    # take daily maximum windspeeds and most frequent wind direction
    group_by(longitude, latitude, date) %>% 
    summarise(
      wind_speed = max(wind_speed, na.rm = TRUE), 
      # wind_dir   = mean(wind_dir, na.rm = TRUE), 
      wind_dir     = max_density(wind_dir),
      .groups = "drop"
    ) %>% 
    distinct(longitude, latitude, date, .keep_all = TRUE)
  
  return(data) 
}

# don't run in parallel, uses too much memory
data <- bind_rows(lapply(files, preprocess_netcdf))

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
met_eir_wind <- left_join(met_eir, data_lbl) %>% 
  relocate(wind_speed, .after = "rain")


#### Save ####

# met_eir_wind_old <- readr::read_csv("data/met_eireann/final/met_eir_wind.csv.gz")
# met_eir_wind <- bind_rows(
# bind_rows(
#   met_eir_wind, 
#   select(met_eir_wind_old, -matches("wind_dir"))
# ) |> 
#   arrange(name, date)

readr::write_csv(met_eir_wind, "data/met_eireann/final/met_eir_wind.csv.gz")

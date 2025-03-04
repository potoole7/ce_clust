#### Download Met Eireann Rainfall Data ####

#### Libs ####

library(stringr)
library(readr)
library(dplyr)
library(lubridate)
library(sf)
library(ggplot2)

sf_use_s2(FALSE)


#### Metadata ####

# set max and min dates over which to look at data
min_date <- as_date("1980-01-01")
max_date <- as_date("2022-12-31")

# Only take stations with start date before 1980

# Have we already downloaded the data?
downloaded <- TRUE 
shp_downloaded <- TRUE

# link to shapefiles 
# (following https://www.cianwhite.org/notebooks/mapping_in_r/)
roi_shp_link <- "https://www.cso.ie/en/media/csoie/census/census2011boundaryfiles/26counties/26_Counties.zip"
roi_shp_name <- basename(roi_shp_link)
ni_shp_link <- "https://osni-spatialni.opendata.arcgis.com/datasets/spatialni::osni-open-data-largescale-boundaries-county-boundaries-.zip?outSR=%7B%22latestWkid%22%3A29902%2C%22wkid%22%3A29900%7D"

# download zip files
link <- "https://cli.fusio.net/cli/climate_data/webdata/dly"
nums <- 1000:4000
urls <- paste0(
  link, nums, ".zip"
)

file_path <- "data/met_eireann/files/"
file_name <- basename(urls)
files <- paste0(file_path, file_name)

#### Download and unzip ####

read_url <- function(url, file) {
  tryCatch(
    {
      download.file(url, file, mode = "wb")
    },
    error = function(cond) {
      message(paste("URL does not seem to exist:", url))
      # message("Here's the original error message:")
      # message(conditionMessage(cond))
      NA
    }
  )
}

if (downloaded == FALSE) {
  lapply(seq_along(urls), \(i) read_url(urls[i], files[i]))
  # files successfully downloaded
  files <- list.files(file_path, full.names = TRUE, pattern = ".zip")
  
  # unzip files
  lapply(files, \(x) {
    system(paste("unzip -o", x))
  })
  
  system("mv dly* data/met_eireann/files")
}

roi_shp_file <- paste0(file_path, roi_shp_name)
if (shp_downloaded == FALSE) {
  # Download ROI shapefile from CSO website
  read_url(roi_shp_link, roi_shp_file)
  system(paste("unzip -o", roi_shp_file))
  system(paste0("mv 26\\ Counties ", file_path, "26_counties"))
}

roi_shp_loc <- tolower(str_remove(roi_shp_file, ".zip"))
ni_shp_loc <- "data/met_eireann/files/OSNI_Open_Data_-_Largescale_Boundaries_-_County_Boundaries/"


#### Join into single dataframe ####

# parse csvs
csvs <- list.files(file_path, full.names = TRUE, pattern = ".csv")

# Pull station locations
stations <- vapply(csvs, \(x) {
  str_to_title(str_remove(read_lines(x, n_max = 1), pattern = "Station Name: "))
}, character(1))

# pull lons & lats
locations <- lapply(csvs, \(x) {
  l3 <- read_lines(x, skip = 2, n_max = 1)
  vapply(stringr::str_split(l3, ",")[[1]], parse_number, numeric(1))
})

# finally, pull rainfall data
precip <- lapply(csvs, \(x) {
  # if colnames don't have rain, skip
  if (grepl("rain", read_lines(x, skip = 9, n_max = 1))) {
    readr::read_csv(x, skip = 8, show_col_types = FALSE)
  } else {
    return(NA)
  }
})

# find stations with data other than rain & ignore
non_na <- which(vapply(precip, length, numeric(1)) != 1)

stations <- stations[non_na]
locations <- locations[non_na]
precip <- precip[non_na]

# also remove no rows
has_rows <- which(vapply(precip, nrow, numeric(1)) != 0)
stations <- stations[has_rows]
locations <- locations[has_rows]
precip <- precip[has_rows]

# join station name and location information with 
met_eir <- lapply(seq_along(precip), \(i) {
  precip[[i]]$name <- stations[i]
  precip[[i]]$lat <- locations[[i]][[1]]
  precip[[i]]$lon <- locations[[i]][[2]]
  precip[[i]]$ind <- as.character(precip[[i]]$ind) # avoid bug
  return(precip[[i]])
}) %>% 
  # join
  bind_rows()

#### Parse ####

# parse date from e.g. 01-jan-1943 to date type
# convert month name to numeric by first converting to factor
met_eir$month <- substr(met_eir$date, 4, 6)
months <- unique(met_eir$month)
met_eir$month <- as.character(as.numeric(
  factor(met_eir$month, levels = months)
))

met_eir <- met_eir %>% 
  mutate(
    # pad with 0s
    month = ifelse(nchar(month) == 1, paste0(0, month), month),
    # replace month name with numeric
    date = paste(substr(date, 0, 2), month, substr(date, 8, 11), sep = "/"),
    # convert to date type
    date = as_date(date, format = "%d/%m/%Y")
  ) %>% 
  select(-c(ind, month))

# only take stations with max date = 2022-12-31 (many with end dates in 2023)
met_eir <- filter(met_eir, date < max_date)

# Only take stations with start date before 1980
station_first_dates <- met_eir %>% 
  group_by(name) %>% 
  slice(1) %>% 
  ungroup() %>% 
  distinct(name, date) %>% 
  # only take (74) stations with data before min date
  filter(date <= min_date) 

met_eir <- semi_join(met_eir, station_first_dates, by = "name") %>% 
  filter(date >= min_date)

# ensure there is no missing data
all_levels <- tidyr::crossing(
  "name" = unique(met_eir$name), 
  "date" = unique(met_eir$date)
)

# 80k rows with missing data! Investigate more closely
missing_data <- anti_join(all_levels, met_eir, by = c("name", "date"))
# Can't just remove stations with missing data (67/74!), ignore for now


#### Convert to sf & plot ####

# Load countries in ROI and NI
# TODO: Add provinces
roi_areas <- st_read(paste0(
  roi_shp_loc, "/Census2011_Admin_Counties_generalised20m.shp"
)) %>% 
  janitor::clean_names() %>% 
  select(nuts2name, countyname) %>% 
  mutate(countyname = str_remove(string = countyname, pattern = " County"))

ni_areas <- st_read(paste0(
  ni_shp_loc, "OSNI_Open_Data_-_Largescale_Boundaries_-_County_Boundaries.shp"
)) %>% 
  janitor::clean_names() %>% 
  select(countyname = county_name) %>% 
  mutate(
    nuts2name  = "Northern Ireland", 
    countyname = str_to_title(countyname), 
    # remove silent letters
    countyname = ifelse(countyname == "Londonderry", "Derry", countyname)
  )

areas <- bind_rows(roi_areas, ni_areas)

# transform to more common WGS 84 crs
irl <- rnaturalearth::ne_countries(country = "Ireland", returnclass = "sf")
areas <- st_transform(areas, crs = st_crs(irl))

# join in province
county <- c(
  # Leinster
  "Carlow", "Dublin", "Kildare", "Kilkenny", "Laois", "Longford", "Louth", 
  "Meath", "Offaly", "Westmeath", "Wexford", "Wicklow", 
  # Munster
  "Clare", "Cork",  "Kerry", "Limerick", "Tipperary", "Waterford", 
  # Connacht
  "Galway", "Leitrim", "Mayo", "Roscommon", "Sligo", 
  # Ulster
  "Cavan", "Donegal", "Monaghan", unique(ni_areas$countyname)
)
provinces <- data.frame(
  "countyname" = county, 
  "province" = c(
    rep("Leinster",times = 12), 
    rep("Munster", times = 6),
    rep("Connacht", times = 5), 
    rep("Ulster", times = 9))
)
areas <- left_join(areas, provinces, by = "countyname")

# plot to inspect
ggplot(areas) +
  geom_sf(aes(fill = province), colour = "black", show.legend = FALSE)

# change data to shapefile
met_eir_sf <- st_as_sf(met_eir, coords = c("lon", "lat"))
st_crs(met_eir_sf) <- st_crs(areas)

# add county and province to sites
met_eir_sf <- st_join(
  met_eir_sf,
  select(areas, county = countyname, province)
)

# save object (return geometry to columns)
coords <- st_coordinates(met_eir_sf) %>% 
  as.data.frame() %>% 
  rename(lon = X, lat = Y)
met_eir_save <- st_drop_geometry(met_eir_sf) %>% 
  bind_cols(coords)

# For plot, just take 2023 sum, and extremes (90th percentile only)
met_eir_sf_plt <- met_eir_sf %>% 
  filter(date >= "2022-01-01") %>% 
  group_by(name) %>% 
  summarise(rain = sum(rain, na.rm = TRUE), .groups = "drop") %>% 
  # take 95% percentile (i.e. "extreme" rain only)
  mutate(percentile = quantile(rain, 0.90)) %>% 
  filter(rain >= percentile)

# set to have the same crs
# st_crs(met_eir_sf) <- st_crs(irl)

ggplot(areas) +
  geom_sf(colour = "black", fill = NA) + 
  geom_sf(data = met_eir_sf_plt, aes(colour = rain), size = 2)
# extreme weather seems to happen in Munster, makes sense given fear of 
# South South-Westerly wind at home!
# Some scattered points in Monaghan and Donegal which are interesting, but
# predictably nothing occurring in Leinster
# Interesting none in Connacht, would have thought it'd be battered!

#### Final preprocessing ####

#### Save ####

readr::write_csv(met_eir_save, "data/met_eireann/final/met_eir.csv.gz")
sf::write_sf(areas, "data/met_eireann/final/irl_shapefile.geojson")

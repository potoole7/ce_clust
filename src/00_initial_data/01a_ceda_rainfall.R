#### Preprocess CEDA MIDAS-Open rainfall data for 6 counties ####

#### libs ####

library(sf)
library(parallel)

sf_use_s2(FALSE)

#### Metadata ####

data_loc <- "data/ceda_open"

ncores <- detectCores() - 1

#### Load Data ####

# areas for plotting
areas <- sf::read_sf("data/met_eireann/final/irl_shapefile.geojson")

# pull counties and files for each county
counties <- list.dirs(data_loc, full.names = FALSE, recursive = FALSE)
# data downloaded from https://data.ceda.ac.uk/badc/ukmo-midas-open/data/uk-daily-rain-obs/dataset-version-202407
files <- unlist(lapply(
  counties, \(x) list.files(paste0(data_loc, "/", x), full.names = TRUE)
))

# data <- mclapply(counties, \(county_spec) {
county_spec <- "tyrone"
file <- "data/ceda_open/tyrone/midas-open_uk-daily-rain-obs_dv-202407_tyrone_01551_edenfel_qcv-1_1990.csv"
data <- mclapply(counties, \(county_spec) {
  # files for a county
  files <- list.files(
    paste0(data_loc, "/", county_spec), pattern = ".csv", full.names = TRUE
  )
  
  # site names for each
  sites <- readr::read_tsv(paste0(
    data_loc, "/", county_spec, "/excel_list_station_details.xls"
  ), show_col_types = FALSE) |> 
    janitor::clean_names() |> 
    select(name, latitude, longitude) |> 
    mutate(
      name = tolower(stringr::str_remove(name, ":")), 
      county = stringr::str_to_title(county_spec)
    )
  
  lapply(files, \(file) {
    total_rows <- length(readr::read_lines(file, progress = FALSE))
    loc <- sites[
      vapply(sites$name, \(x) {
        grepl(x, stringr::str_replace_all(file, "-", " "), fixed = TRUE)
      }, FUN.VALUE = logical(1)),
    ]
    # if failing, download new list of locations
    # https://archive.ceda.ac.uk/tools/midas_stations 
    stopifnot("No unique match for site" = nrow(loc) == 1)
    readr::read_csv(
      file, skip = 62, n_max = total_rows - 64, show_col_types = FALSE
    ) |> 
      janitor::clean_names() |> 
      select(date = ob_date, rain = prcp_amt) |> 
      mutate(date = lubridate::as_date(date)) |> 
      cbind(loc) |> 
      relocate(name, county) |> 
      mutate(
        filename = stringr::str_remove(
          stringr::str_remove(
            basename(file), 
            paste0("midas-open_uk-daily-rain-obs_dv-202407_", county_spec, "_")
          ), 
          ".csv"
        )
      )
  })
}, mc.cores = ncores)

# save as it's a pain to reproduce
saveRDS(data, "data/ceda_open/ceda_rain.rds")

data_join <- bind_rows(data) |> 
  select(-filename) |> 
  rename(lon = longitude, lat = latitude) |> 
  mutate(province = "Ulster", name = stringr::str_to_title(name))

# save as dataframe
readr::write_csv(data_join, "data/ceda_open/ceda_rain.csv.gz")

#### Plot ####

# plot to look at sites
ggplot(areas) +
  geom_sf(colour = "black", fill = NA) + 
  geom_sf(
    data = st_to_sf(data_join), 
    aes(colour = name), 
    size = 2
  )
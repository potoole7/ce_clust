#### Download EPA data using RAQSAPI ####

# - Daily summaries of air pollution data from the EPA
# - Specifically looking at PM2.5, PM10, SO2, NO, CO, HC
# - Want to download cities from 2/3 distinct industry/spatial regions in US 
# (~ follow https://www.ondeck.com/resources/biggest-industries-in-the-world):
#   - Rust Belt (Manufacturing): Baltimore, Buffalo, Chicago, Cleveland, Detroit, Pittsburgh, Rochester, and St. Louis
#   - West Coast(Tech): San Francisco, Seattle, Los Angeles, San Diego, Portland, San Jose
#     - Perhaps add Boston and another Tech city here? Although Austin is one..
#   - Texas/South (Construction/oil): Austin, Houston, Dallas, San Antonio, New Orleans, Albuquerque
#   - Maybe also add New England?? Different industry there


#### Libraries ####

# install.packages(pkgs="RAQSAPI", dependencies = TRUE )
library(RAQSAPI)
library(lubridate)
library(dplyr)
library(parallel)


#### Metadata ####

email <- "pot23@bath.ac.uk"
pass <- "taupefox76" # TODO: Have as environment variable!!

# parameter numbers for various pollutants
# TODO: Some parameters not coming through? Remove Carbon Monoxide
params <- c(
  "PM1"   = "81102",
  "PM2.5" = "88101",
  "PM10"  = "81103",
  "SO2"   = "42401",
  "NO"    = "42602",
  "CO"    = "42101",
  "HC"    = "45201" # Benzene
)

# state numeric codes for cities of interest
# TODO: Could change to be a key?
state_codes <- list(
  # Rust Belt
  "MD" = "24", "NY" = "36", "IL" = "17", "OH" = "39", 
  "MI" = "26", "PA" = "42", "NY" = "36", "MO" = "29",
  # West Coast
  "CA" = "06", "WA" = "53", "OR" = "41", "CA" = "06", "CA" = "06", "CA" = "06",
  # Texas/South
  "TX" = "48", "TX" = "48", "TX" = "48", "TX" = "48", "LA" = "22", "NM" = "35"
)

# Temp: Add missing cities
state_codes <- list(
  "CA" = "06", "OR" = "41", "MO" = "29"
)

# countycodes for cities
county_codes <- list(
  # Rust Belt
  "Baltimore"  = "510", "Buffalo" = "029", "Chicago" = "031", "Cleveland" = "035",
  "Detroit"    = "163", "Pittsburgh" = "003", "Rochester" = "055", "St. Louis" = "510",
  # West Coast
  "San Francisco" = "075", "Seattle" = "033", "Los Angeles" = "037", "San Diego" = "073",
  "Portland" = "051", "San Jose" = "085",
  # Texas/South
  "Austin" = "453", "Houston" = "201", "Dallas" = "113", "San Antonio" = "029",
  "New Orleans" = "071", "Albuquerque" = "001"
)

# Temp: just do missing cities
county_codes <- list(
  "Los Angeles" = "037", "Portland" = "051", "St. Louis" = "510"
)

# start and end dates 
end_date <- as_date("2024/01/01")
start_date <- as_date("1980/01/01") # from Minnesota paper

ncores <- min(detectCores() - 1, length(params))
# ncores <- detectCores() - 3


#### Main ####

# signup, if not done already
# aqs_sign_up(email)
# provide credentials
aqs_credentials(username = email, key = pass)

# pull daily summaries for each pollutant and county
data <- mclapply(names(params), \(param) {
# data <- lapply(names(params), \(param) {
  lapply(seq_along(county_codes), \(i) {
    aqs_dailysummary_by_county(
      parameter  = params[[param]],
      bdate      = start_date,
      # bdate      = as_date("2023/01/01"), # test with smaller date range
      edate      = end_date,
      # stateFIPS  = state_code,
      stateFIPS  = state_codes[[i]],
      countycode = county_codes[[i]]
    )
      # dplyr::select(
      #   variable = parameter, 
      #   date     = date_local, 
      #   unit     = units_of_measure, 
      #   mean     = arithmetic_mean, 
      #   state, county, city, 
      #   latitude, longitude # , datum # can ignore datum, negligible diff
      # ) %>% 
      # mutate(date = as_date(date))
  })
}, mc.cores = ncores)
# })
# saveRDS(data, file = "data/epa_air_pollution_data.RDS")

# remove errors
# TODO: Double check this doesn't accidentally remove lots of data?!
data_clean <- lapply(data, \(x) {
  if (inherits(x, "try-error")) return(NULL)
  lapply(x, \(y) {
    if (inherits(y, "try-error")) return(NULL)
    return(y)
  })
})

# join, cleanup and save
data_join <- bind_rows(data_clean) %>% 
  dplyr::select(
    variable = parameter, 
    date     = date_local, 
    unit     = units_of_measure, 
    mean     = arithmetic_mean, 
    state, county, city, 
    latitude, longitude # , datum # can ignore datum, negligible diff
  ) %>% 
  filter(!is.na(latitude), !is.na(longitude)) %>% 
  mutate(
    date = as_date(date),
    variable = case_when(
      variable == "PM10 Total 0-10um STP"    ~ "PM10",
      variable == "PM2.5 - Local Conditions" ~ "PM2.5", 
      TRUE                                   ~ variable
    )
  )

data_join3 <- rbind(data_join, data_join1)
setdiff(cities, data_join3$city)

readr::write_csv(data_join, "data/epa_air_pollution_data.csv")

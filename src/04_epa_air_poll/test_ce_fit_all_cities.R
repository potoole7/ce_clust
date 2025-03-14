#### Fit CE model to Los Angles and San Francisco NO2 and PM10 ####

# As a simplification, take max obs across sites

# TODO Fit for all variables (done)
# TODO Plot a, b values on map
# TODO Fit for each variable separately and compare (data availability will differ)
# TODO Quantify uncertainty in CE model fit (block bootstrap most likely)


#### Libraries ####

library(dplyr, quietly = TRUE)
library(tidyr)
library(sf)
library(ggplot2)
# library(evc)
devtools::load_all("../evc")
library(ggplot2)
library(USAboundaries)
library(gt)

source("src/functions.R")

theme <- evc::evc_theme()

sf::sf_use_s2(FALSE)


#### Functions ####

# Get days with full record of chosen variables for each city
full_record <- \(data) {
  data %>% 
    pivot_wider(names_from = variable, values_from = mean) %>% 
    janitor::clean_names() %>% 
    filter(if_all(any_of(c("no2", "pm10", "so2")), ~ !is.na(.))) %>% 
    return()
}

#### Metadata ####

# cities and variables over which to test CE model
# TODO Plot # days w/ full records for 2/3 vars for each city 
# cities <- c("Los Angeles", "San Francisco")
cities <- c(
  # Baltimore: SO2 missing post 2000 (Investigate on online EPA platform!)
  "Baltimore",
  # Buffalo: No, PM10 missing 2000-~2018, good record for NO2 and SO2 though
  "Buffalo",
  "Chicago", "Cleveland", "Detroit", "Pittsburgh",
  # - Rochester: No, only good record for So2
  "Rochester",
  "St. Louis", "San Francisco",
  # Seattle: NO2 has two strange gaps, try for now!
  "Seattle", "Los Angeles",
  # San Diego some gaps for PM10 (~2017-2019), SO2 (post 2012)
  "San Diego",
  # Portland: Yes, but gap around 2000 for SO2
  "Portland",
  # San Jose: Yes, but SO2 only for pre 2010
  "San Jose",
  # Austin: >= 2 gaps for all variables! (or maybe include after?)
  "Austin",
  "Houston",
  # Dallas: small PM10 and S02 gaps
  "Dallas",
  # San Antonio: SO2 has many gaps, NO2 has gap at start
  "San Antonio",
  # New Orleans: Large gaps and no SO2!
  # Albuquerque: Almost no SO2!
  "New Orleans", "Albuquerque"
)
# TODO: Look into CO2 and O3 as in Vettori paper
variables <- c(
  "NO2" = "Nitrogen dioxide (NO2)", "PM10" = "PM10", "SO2" = "Sulfur dioxide"
)
# Don't look at SO2, too many gaps (for now anyway)
# variables <- c("NO2" = "Nitrogen dioxide (NO2)", "PM10" = "PM10")


#### Load Data ####

data <- readr::read_csv("data/epa_air_pollution_data.csv") %>% 
# data <- readr::read_csv("data/epa_air_pollution_data_2.csv") %>% 
  filter(city %in% cities, variable %in% variables) %>% 
  mutate(
    variable = case_when(
      grepl("NO2", variable)       ~ "NO2", 
      variable == "Sulfur dioxide" ~ "SO2", 
      TRUE ~ variable
    )
  )

# Continental US shapefile for plotting
# areas <- us_states() %>% 
#   # areas <- USAboundaries::us_counties() %>% 
#   filter(!state_abbr %in% c("HI", "AK", "AS", "GU", "MP", "PR", "VI"))


#### Deal with multiple sites ####

# TEMP SOLUTION
# take max obs for each city across all sites (i.e. ignoring site)
data_city <- data %>% 
  group_by(city, variable, date) %>%
  summarise(
    mean = max(mean, na.rm = TRUE),
    # just take first lon and lat vals, won't be that different in same city
    longitude = first(longitude), 
    latitude  = first(latitude),
    .groups = "drop"
  )


#### Exploratory plots ####

# plot ACF and PACF
# TODO Facet per variable and acf/pacf, rather than list
ts_plot <- \(x, city, var, type = c("ACF", "PACF")) {
  # use either acf or pacf function depending on choice
  fun <- acf
  if (type == "PACF") {
    fun <- pacf
  }
  # filter data for specific city
  x_spec <- x %>%
    # filter(name == city)
    filter(city == !!city, variable == var)
  if (nrow(x_spec) == 0) return(NA)
  # ts_spec <- fun(x_spec[[tolower(var)]], plot = FALSE)
  ts_spec <- fun(x_spec$mean, plot = FALSE)
  with(ts_spec, data.frame(lag, acf)) %>% 
    ggplot(aes(lag, acf)) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = lag, yend = 0)) + 
    theme + 
    labs(y = type, x = "Lag", title = paste0(type, ", ", city, ", ", var))
}

# plot for each city and variable
ts_plots <- lapply(cities, function(city) {
  lapply(names(variables), function(var) {
    # print(paste0(city, " - ",  var))
    list(
      "acf"  = ts_plot(data_city, city, var, "ACF"),
      "pacf" = ts_plot(data_city, city, var, "PACF")
    )
  })
})

# TODO Save
# pdf("plots/acf_pacf_plots.pdf")
# ts_plots
# dev.off()


# plot # days w/ full records for combinations of variables
count_record <- \(x, vars = c("NO2", "SO2", "PM10"), name = "n Full") {
  x %>% 
    filter(variable %in% vars) %>% 
    full_record() %>% 
    dplyr::count(city, name = name) %>% 
    return()
}

n_records <- count_record(data_city) %>% 
  full_join(
    count_record(data_city, vars = c("NO2", "PM10"), name = "n NO2 & PM10")
  ) %>% 
  full_join(
    count_record(data_city, vars = c("NO2", "SO2"), name = "n NO2 & SO2")
  ) %>%
  full_join(
    count_record(data_city, vars = c("PM10", "SO2"), name = "n PM10 & SO2")
  )

# simple table to show records
gt(n_records)

# how many cities have > 300 records for each variable
n_records %>% 
  pivot_longer(contains("n")) %>% 
  group_by(name) %>% 
  summarise(n = sum(value >= 300, na.rm = TRUE), .groups = "drop")
# most for NO and SO2, unsuprisingly

# barplot
p <- n_records %>% 
  pivot_longer(contains("n")) %>% 
  ggplot(aes(x = city, y = value, fill = name)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ name, scales = "free_y") + 
  theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p

# conclusion:  

#### CE model setup ####

# change to wide, remove cases where both vars not observed, format for evc
# TODO: Change evc to be able to account for missing data better?
# Could have data in lists for each city and variable, would allow us to 
# estimate marginal distributions using different variables
data_wide <- data_city %>% 
  # Temp: Only look at NO2 and SO2, most cities available
  # filter(variable %in% c("NO2", "SO2")) %>% 
  pivot_wider(names_from = variable, values_from = mean) %>% 
  janitor::clean_names() %>% 
  # filter(if_all(any_of(c("no2", "so2")), ~ !is.na(.))) %>% 
  filter(if_all(any_of(c("no2", "pm10", "so2")), ~ !is.na(.))) %>% 
  rename(name = city) %>% 
  dplyr::select(-date)

# remove cities with < 200 observations
keep_cities <- data_wide %>%
  count(name) %>%
  filter(n >= 400) %>%
  pull(name)
data_wide <- filter(data_wide, name %in% keep_cities)


#### Fit CE ####

fit <- fit_ce(
  data        = data_wide, 
  # vars        = c("no2", "pm10", "so2"),
  vars        = names(data_wide)[4:ncol(data_wide)],
  # marg_prob   = 0.9,
  marg_prob   = list(
    f         = list("response ~ name", "~ name")
    # f = list("response ~ s(lon, lat, k = 50)", "~ s(lon, lat, k = 40)")
  ), 
  cond_prob   = 0.7, # TODO Implement bootstrap plots to decide this
  # formula for evgam fit to excesses over threshold
  # f = list(excess ~ s(longitude, latitude, k = 5), ~s(longitude, latitude, k = 5)), 
  f           = list(excess ~ name, ~name), 
  fit_no_keef = FALSE, # don't fit without Keef constraints
  output_all  = TRUE
)

# pull dependence parameters for each city
dependence <- fit$dependence


#### Plot results ####

# TODO Look to Ireland analysis for what plots to do!

# TODO Plot marginal estimates on map

# TODO Plot dependence parameters on map

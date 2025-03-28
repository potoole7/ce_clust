#### Exploration of US air pollution data ####

# TODO Plot chi-squared on map for different variables

#### Libs ####

# library(tigris)
library(USAboundaries)
library(sf)
library(evc)
# devtools::load_all("../evc")
library(dplyr, quietly = TRUE)
library(tidyr)
library(ggplot2)
library(lubridate)
# devtools::load_all("texmex") # forked texmex which allows for fixed b vals
library(gridExtra)
library(patchwork)
library(evgam)
library(geosphere)
library(latex2exp)
library(ggpattern)

source("src/functions.R")

theme <- evc::evc_theme()

sf::sf_use_s2(FALSE)
# options(tigris_class = "sf")


#### Metadata ####

params <- c("") # Params to keep (see Vettori paper)

cities <- c(
  "Baltimore", "Buffalo", "Chicago", "Cleveland", "Detroit", 
  "Pittsburgh", "Rochester", "St. Louis", "San Francisco", "Seattle", 
  "Los Angeles", "San Diego", "Portland", "San Jose", "Austin", 
  "Houston", "Dallas", "San Antonio", "New Orleans", "Albuquerque"
)

#### Load Data ####

data <- readr::read_csv("data/epa_air_pollution_data.csv") %>% 
# data <- readr::read_csv("data/epa_air_pollution_data_2.csv") %>% 
  filter(
    # !variable %in% c("Carbon monoxide", "PM10-2.5 STP"),
    city %in% cities
  ) %>% 
  mutate(variable = ifelse(variable == "PM10-2.5 STP", "PM10", variable)) %>%
  # remove duplicates
  # TODO: Investigate further why these exist
  distinct(variable, date, mean, state, city, latitude, longitude)
  
# pull shapefile for continental US
areas <- us_states() %>% 
# areas <- USAboundaries::us_counties() %>% 
  filter(!state_abbr %in% c("HI", "AK", "AS", "GU", "MP", "PR", "VI"))

# # convert to sf 
# data_sf <- data %>% 
#   filter(!is.na(longitude)) %>%  # some missing values for Chicago
#   sf::st_as_sf(
#     coords = c("longitude", "latitude"), 
#     crs    = st_crs(areas)
#   )


#### Initial Plotting ####

# tally number of sites for each city
# TODO: What do we do about different sites for different variables?
data %>% 
  group_by(city, variable) %>% 
  distinct(longitude, latitude) %>% 
  summarise(n_sites = n(), .groups = "drop") %>% 
  arrange(desc(n_sites))

# plot time series for each county and variable
data %>%
  ggplot(aes(x = date, y = mean, colour = variable)) +
  geom_line() +
  facet_wrap(~city, scales = "free_y") +
  theme +
  labs(
    title = paste0("Air Pollution Time Series"),
    x = "Date",
    y = "Mean"
  )

# look at air pollution map for e.g. PM10
ggplot(areas) + 
  geom_sf(colour = "black", fill = NA) + 
  geom_sf(
    data = data %>% 
      sf::st_as_sf(
        coords = c("longitude", "latitude"), 
        crs    = st_crs(areas)
      ) %>% 
      filter(
        variable == "PM10", 
        date == as_date("2001-09-22"), 
        !city == "Not in a city"
      ),
    aes(colour = city), 
    size = 3
  ) + 
  theme


#### Investigate multiple sites for each city ####

# add unique identifiers for each site within each city
sites <- data %>% 
  distinct(city, longitude, latitude) %>% 
  group_by(city) %>% 
  mutate(site = row_number()) %>% 
  ungroup()

data <- left_join(data, sites)

# For each city, plot different sites for different variables
ts_plts <- lapply(cities, \(x) {
  data %>%
    filter(city == x) %>%
    # ggplot(aes(x = date, y = mean, colour = variable)) +
    ggplot(aes(x = date, y = mean, colour = factor(site))) +
    geom_line(alpha = 0.5) +
    guides(colour = "none") + 
    facet_wrap(~variable, scales = "free_y") +
    theme +
    labs(
      title = x,
      x = "Date",
      y = "Mean", 
      colour = "Site"
    )
})

# save all plots in pdf
pdf("plots/01_ts_sites.pdf", width = 8, height = 6)
# pdf("plots/01_ts_sites_2.pdf", width = 8, height = 6)
ts_plts
dev.off()

# Tabulate percentage of data available for each site, plot as bar chart
# TODO Fix this, consistently incorrect
data_tbl <- data %>% 
  left_join(
    # find # unique dates for each city (across all vars, not individually)
    data %>% 
      group_by(city) %>% 
      summarise(
        n_obs = as.numeric(max(data$date) - min(data$date)), 
        .groups = "drop"
      )
  ) %>% 
  left_join(
    # find how many of these dates have valid obs for a given/each site
    data %>% 
      filter(!is.na(mean)) %>% 
      # remove duplicates
      # TODO: Look into why these still persist!
      distinct(site, city, variable, date) %>%
      # group_by(city, variable, site) %>% 
      # summarise(
      #   n_non_na = as.numeric(max(data$date) - min(data$date)), 
      #   .groups = "drop"
      # ) %>% 
      count(city, variable, site, name = "n_non_na") %>% 
      identity()
  ) %>% 
  mutate(completeness = n_non_na / n_obs)

data_plt <- data_tbl %>% 
  distinct(site, city, variable, completeness)
  
bar_plts <- lapply(cities, \(x) {
  data_plt %>%
    filter(city == x) %>% 
    ggplot(aes(x = site, y = completeness, fill = factor(site))) +
    facet_wrap(~variable) +
    scale_y_continuous(limits = c(0, 1)) + 
    # geom_bar(stat = "identity") +
    geom_col() + 
    guides(fill = "none") +
    theme +
    labs(
      title = x,
      x = "Site",
      y = "% Completeness",
      fill = "Site"
    )
})

pdf("plots/02_bar_sites.pdf", width = 8, height = 6)
bar_plts
dev.off()

# also do for every city (ignoring site)
# TODO: Functionalise above, repeated twice!

data_tbl_city <- data %>% 
  # left_join(
  #   # find number of unique dates 
  #   data %>% 
  #     group_by(city, variable) %>% 
  #     summarise(n_obs = n_distinct(date), .groups = "drop")
  # ) %>% 
  # left_join(
  #   # find how many of these dates have observations for a given site
  #   data %>% 
  #     filter(!is.na(mean)) %>% 
  #     group_by(city, variable) %>% 
  #     summarise(n_non_na = n_distinct(date), .groups = "drop")
  # ) %>% 
  left_join(
    # find # unique dates for each city (across all vars, not individually)
    data %>% 
      group_by(city, variable) %>% 
      summarise(
        n_obs = as.numeric(max(data$date) - min(data$date)), 
        .groups = "drop"
      )
  ) %>% 
  left_join(
    # find how many of these dates have valid obs for a given/each site
    data %>% 
      filter(!is.na(mean)) %>% 
      # remove duplicates
      # TODO: Look into why these still persist!
      distinct(city, variable, date) %>%
      count(city, variable, name = "n_non_na") %>% 
      identity()
  ) %>%  
  mutate(completeness = n_non_na / n_obs)

# data_tbl %>% 
#   select(variable, city, completeness) %>% 
#   group_by(variable, city) %>% 
#   filter(completeness == max(completeness)) %>% 
#   slice(1) %>% 
#   ungroup() %>% 
#   pivot_wider(names_from = variable, values_from = completeness) %>% 
#   View()

data_plt2 <- data_tbl_city %>% 
  distinct(city, variable, completeness)
  
p <- data_plt2 %>% 
  # ggplot(aes(x = city, y = completeness, fill = factor(city))) +
  # facet_wrap(~variable) +
  # make variable names shorter for plotting
  # TODO: Functionalise
  mutate(variable = case_when(
    grepl("NO2", variable) ~ "NO2",
    grepl("Benzene", variable) ~ "HC",
    grepl("Sulfur dioxide", variable) ~ "SO2",
    TRUE ~ variable
  )) %>% 
  ggplot(aes(x = variable, y = completeness, fill = factor(variable))) +
  facet_wrap(~city) +
  scale_y_continuous(limits = c(0, 1)) + 
  geom_col() + 
  guides(fill = "none") +
  theme +
  labs(
    x = "City",
    y = "% Completeness",
    fill = "Site"
  ) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

ggsave("plots/03_bar_cities.pdf", p, width = 8, height = 6)


#### OLD OLD OLD ####

# How many rows left when only dates common to all counties are kept?
# TODO: Fix
common_dates <- data %>%
  select(variable, date, county, mean) %>%
  arrange(desc(date)) %>%
  group_by(date, variable, county) %>%
  slice(1) %>% # just want one obs per county (multiple sites per county)
  ungroup() %>%
  pivot_wider(names_from = county, values_from = mean) %>%
  filter(across(everything(), ~ !is.na(.))) %>%
  pull(date)


# ~94%, not too bad!!
nrow(filter(data, date %in% common_dates)) / nrow(data)


#### Investigate multiple sites across one city ####

# look at specific city and state
# areas_spec <- counties(state = "IL") %>% 
#   filter(NAME == "Cook")
areas_spec <- counties(state = "CA") %>% 
  filter(NAME == "Los Angeles")
data_spec <- data %>% 
  filter(city == "Los Angeles", variable == "PM10")

# ggplot(areas_spec) + 
# crop to small area around points
data_spec_sf <- data_spec %>% 
  sf::st_as_sf(
    coords = c("longitude", "latitude"), 
    crs    = st_crs(areas)
  )
ggplot(st_crop(areas_spec, st_bbox(data_spec_sf))) +
  geom_sf(colour = "black", fill = NA) + 
  geom_sf(
    data = data_spec %>% 
      sf::st_as_sf(
        coords = c("longitude", "latitude"), 
        crs    = st_crs(areas)
      ) %>% 
      filter(date == as_date("2001-09-22")),
    aes(colour = city), 
    size = 3
  ) + 
  theme


#### Complete Cases ####

# data_spec %>% 
#   filter(variable == "PM10") %>% 
#   dplyr::select(date, mean) %>% 
#   pivot_wider(names_from = date, values_from = mean) %>% 
#   filter(across(everything(), ~ !is.na(.))) %>% 
#   dplyr::select(-variable) %>% 
#   pivot_longer(everything(), names_to = "date", values_to = "mean") %>% 
#   ggplot(aes(x = date, y = mean)) +
#   geom_line() +
#   theme +
#   labs(
#     title = paste0("Air Pollution Time Series"),
#     x = "Date",
#     y = "Mean"
#   )  

sites <- data_spec %>% 
  distinct(longitude, latitude) %>% 
  mutate(site = row_number())

data_spec %>% 
  left_join(sites) %>% 
  ggplot(aes(x = date, y = mean, colour = factor(site))) +
  geom_line()
  

#### Investigate multiple sites across all cities ####

#### Seasonality ####

#### Stationarity ####


#### Fit CE model ####

# TODO: Need to change to wide??
mod_fit <- fit_ce(
  data = data %>% 
    rename(name = city) %>% 
    filter(variable %in% c("PM2.5", "Sulfur dioxide")) %>% 
    pivot_wider(names_from = variable, values_from = mean) %>% 
    filter(if_all(c("PM2.5", "Sulfur dioxide"), ~ !is.na(.))),
  vars = c("PM2.5", "Sulfur dioxide"),
  # arguments to `quantile_thresh`
  # marg_prob = list(
  #   f = list(
  #     "response ~ s(longitude, latitude, k = 3)", 
  #     "~ s(longitude, latitude, k = 3)"
  #   )
  # ), 
  # marg_prob = list(
  #   f = list("response ~ name", "~ name")
  # ),
  marg_prob = 0.9,
  cond_prob = 0.7, # TODO: What conditional probability should I use??
  # formula for evgam fit to excesses over threshold
  # f = list(excess ~ s(longitude, latitude, k = 3), ~s(lon, lat, k = 3)), 
  # f = list(excess ~ name, ~ name), 
  f = NULL, # fit with ~ismev::gpd.fit~ for each location
  fit_no_keef = FALSE, # don't fit without Keef constraints
  output_all = TRUE
)

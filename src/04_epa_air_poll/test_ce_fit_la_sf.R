#### Fit CE model to Los Angles and San Francisco NO2 and PM10 ####

# As a simplification, take max obs across sites
# TODO Plot ACF and PACF
# TODO Block bootstrap for uncertainty

#### Libraries ####

library(dplyr, quietly = TRUE)
library(tidyr)
library(sf)
library(ggplot2)
# library(evc)
devtools::load_all("../evc")
library(ggplot2)
library(USAboundaries)

source("src/functions.R")

theme <- evc::evc_theme()

sf::sf_use_s2(FALSE)

#### Metadata  ####

# cities and variables over which to test CE model
cities <- c("Los Angeles", "San Francisco")
variables <- c(
  "NO2" = "Nitrogen dioxide (NO2)", "PM10" = "PM10", "SO2" = "Sulfur dioxide"
)


#### Load Data ####

data <- readr::read_csv("data/epa_air_pollution_data.csv") %>% 
# data <- readr::read_csv("data/epa_air_pollution_data_2.csv") %>% 
  filter(
    city %in% cities, variable %in% variables
  ) %>% 
  mutate(
    variable = case_when(
      grepl("NO2", variable)       ~ "NO2", 
      variable == "Sulfur dioxide" ~ "SO2", 
      TRUE ~ variable
    )
  )
 

#### Preprocess ####

# take max obs for each city across all sites (i.e. ignoring site)
data_city <- data %>% 
  group_by(city, variable, date) %>%
  summarise(
    mean = max(mean, na.rm = TRUE), .groups = "drop", 
    # just take first lon and lat vals, won't be that different in same city
    longitude = first(longitude), 
    latitude  = first(latitude)
  )

# change to wide, remove cases where both vars not observed, format for evc
# TODO: Change evc to be able to account for missing data better?
# Could have data in lists for each city and variable, would allow us to 
# estimate marginal distributions using different variables
data_wide <- data_city %>% 
  pivot_wider(names_from = variable, values_from = mean) %>% 
  janitor::clean_names() %>% 
  # filter(!is.na(no2), !is.na(pm10)) %>% 
  filter(!is.na(no2), !is.na(pm10), !is.na(so2)) %>% 
  rename(name = city) %>% 
  select(-date)

# 5,532 records, not bad considering!
  
#### Plotting ####

# plot ACF and PACF
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
    list(
      "acf"  = ts_plot(data_city, city, var, "ACF"),
      "pacf" = ts_plot(data_city, city, var, "PACF")
    )
  })
})

#### Fit CE ####

fit <- fit_ce(
  data = data_wide, 
  # vars = c("no2", "pm10", "so2"),
  vars <- names(data_wide)[4:ncol(data_wide)],
  # marg_prob = 0.95, 
  marg_prob = list(
    # f = list("response ~ s(lon, lat, k = 50)", "~ s(lon, lat, k = 40)")
    f = list("response ~ name", "pm10 ~ name")
  ),
  cond_prob = 0.7, # TODO: Implement boostrap plots to decide this
  # formula for evgam fit to excesses over threshold
  # f = list(excess ~ s(lon, lat, k = 40), ~s(lon, lat, k = 40)), 
  f = list(excess ~ name, ~name), 
  fit_no_keef = FALSE, # don't fit without Keef constraints
  output_all = TRUE
)

fit$dependence$`Los Angeles`
fit$dependence$`San Francisco`

# Seem to all have very high a vals, and low positive b vals, indicating 
# strong positive extremal dependence between air pollution variables

# Odd one out is SO2 for LA, doesn't have extremal dependence with PM10

#### Block boostrapping uncertainty ####
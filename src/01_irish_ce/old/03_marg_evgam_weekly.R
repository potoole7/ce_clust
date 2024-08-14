#### Marginal analysis with evgam and INLA, variable threshold ####

# following https://empslocal.ex.ac.uk/people/staff/by223/software/gpd.R
# In paper https://www.tandfonline.com/doi/full/10.1080/01621459.2018.1529596

# Fit varying threshold for evgam  (done)
# Plot exceedances over varying threshold as diagnostic (done)
# Model logscale (done)
# TODO: Test with weekly data, as in KL diverenge clustering paper?
# TODO: Model shape? As in https://www.tandfonline.com/doi/full/10.1080/01621459.2018.1529596 
# TODO: Calculate return levels?
# TODO: Add altitude to plots, to see if there is any link between the two?

# TODO: Investigate 2 rows for Castleisland (Coom) in data_gpd, why not 1?
# TODO: Investigate large wind scale in Wexford, why is this the case? 
# - Do exploratory analysis on wind speeds (and rainfall while you're at it)

#### Libs ####

library(sf)
library(dplyr, quietly = TRUE)
library(tidyr)
library(ggplot2)
library(lubridate)
library(texmex)
library(gridExtra)
library(patchwork)
library(evgam)

source("src/functions.R")

theme <- theme_bw() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 13, face = "bold"),
    strip.background = element_rect(fill = NA, colour = "black"),
    plot.tag = element_text(size = 16, face = "bold"),
    panel.background = element_rect(fill = NA, colour = "black")
  )

sf::sf_use_s2(FALSE)


#### Metadata ####

min_max_dates <- as_date(c("1990-01-01", "2020-12-31"))
all_dates <- seq.Date(min_max_dates[1], min_max_dates[2], by = "day")


#### Load data ####

data <- readr::read_csv(
  "data/met_eireann/final/met_eir_wind.csv.gz",
  show_col_types = FALSE
) %>%
  filter(date %in% all_dates, !is.na(wind_speed))

# load shapefile
areas <- sf::read_sf("data/met_eireann/final/irl_shapefile.geojson")


#### Convert to weekly ####

# Take weekly sum of precipitation and average of wind speed, to account for 
# lag in storm (as in Vignotto, Engelke 2021 study of GB + Ireland)
# TODO: Take weely averages of *daily* wind speed maxima!!
data <- data %>%
  mutate(week = floor_date(date, "week")) %>%
  group_by(week, name, county, province, lon, lat) %>%
  summarise(
    rain = sum(rain, na.rm = TRUE),
    wind_speed = mean(wind_speed, na.rm = TRUE), 
    .groups = "drop"
  ) %>% 
  rename(date = week) # to agree with below code

#### Exploratory Plots ####

# TODO: Plot chi for each location (texmex has function for this)

# TODO: Plot max weekly rainfall

# TODO: Plot max mean wind speed


#### Threshold selection ####

# function to fit varying threshold using quantile regression
# https://empslocal.ex.ac.uk/people/staff/by223/software/gpd.R
quantile_thresh <- function(data, response, tau = .97) {
  fm_ald <- list(
    # response ~ s(lon, lat, k = 50), # location
    formula(paste(response, " ~ s(lon, lat, k = 50)")), # location
    ~ s(lon, lat, k = 40)                               # logscale
  )
  
  # fit the quantile regression model using 97th percentile
  ald_fit <- evgam(fm_ald, data, family = "ald", ald.args = list(tau = tau))
  
  # add threshold to data and filter
  data %>% 
    mutate(
      thresh = predict(ald_fit)$location, 
      # excess = wind_speed - thresh
      excess = !!sym(response) - thresh
    ) %>% 
    filter(excess > 0) %>% 
    return()
}

# Threshold rainfall and wind speed with varying threshold per location
data_rain <- quantile_thresh(data, "rain")
data_ws <- quantile_thresh(data, "wind_speed")

# plot exceedances for each (~25-75/931-1462)
data_thresh_plt <- data %>% 
  distinct(name, lon, lat) %>% 
  left_join(
    count(data_rain, name) %>% 
      rename(n_rain = n)
  ) %>% 
  left_join(
    count(data_rain, name) %>% 
      rename(n_wind = n)
  ) %>% 
  # mutate(across(contains("n_"), \(x) ifelse(is.na(x), 0, x))) %>% 
  pivot_longer(n_rain:n_wind, names_to = "var", values_to = "n") %>% 
  st_as_sf(coords = c("lon", "lat"))
st_crs(data_thresh_plt) <- "WGS84"

ggplot(areas) +
  geom_sf(colour = "black", fill = NA) +
  geom_sf(
    data = data_thresh_plt,
    # aes(fill = n, size = n),
    aes(fill = n),
    size = 6,
    stroke = 1,
    pch = 21
  ) +
  # geom_sf_text(
  #   data = data_thresh_plt,
  #   aes(label = n), size = 6
  # ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~var) +
  # scale_size_continuous(
  #   # breaks = scales,
  #   breaks = c(0, 50, 100, 150),
  #   labels = as.character(c(0, 50, 100, 150)),
  #   guide = "legend"
  # ) +
  # scale_fill_gradientn(
  #   colours = RColorBrewer::brewer.pal(name = "Blues", n = 7),
  #   na.value = "purple",
  #   breaks = c(0, 50, 100, 150),
  #   labels = as.character(c(0, 50, 100, 150)),
  #   guide = "legend"
  # ) +
  labs(fill = "", size = "") +
  guides(fill = guide_legend(), size = guide_legend()) +
  theme +
  theme(legend.position = "right") +
  NULL


#### Fit evgam to rain #### 

# Fit evgam model with spline for spatial location, create predictions
fit_evgam <- \(data, pred_data) {
  # model formula
  f <- list(excess ~ s(lon, lat, k = 12), ~ 1) # TODO: smooth for year?
  # f <- list(excess ~ s(lon, lat, k = 12) + s(date, k = 12), ~ 1) # TODO: smooth for year?
  # fit evgam model
  m <- evgam(f, data = data, family = "gpd")

  # create predictions
  predictions <- predict(m, pred_data, type = "response")
  
  # return model fit and predictions
  return(list(
    "m"           = m,
    "predictions" = predictions
  ))
}

evgam_rain <- fit_evgam(data_rain, data_rain)

# TODO: Functionalise 
data_rain_map <- data_rain %>%
  mutate(scale = evgam_rain$predictions$scale) %>%
  st_as_sf(coords = c("lon", "lat"), remove = FALSE)
st_crs(data_rain_map) <- "WGS84"

# plot scale parameter over space
plot_scale <- \(areas, data_map, scales) {
  ggplot(areas) +
    geom_sf(colour = "black", fill = NA) +
    geom_sf(
      data = data_map,
      aes(fill = scale, size = scale),
      # size = 6,
      stroke = 1,
      pch = 21
    ) +
    # scale_colour_viridis_b() +
    labs(fill = "Scale") +
    scale_fill_gradientn(
      colours = RColorBrewer::brewer.pal(name = "Blues", n = 7)
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_size_continuous(
      breaks = scales,
      labels = as.character(scales),
      guide = "legend"
    ) +
    scale_fill_gradientn(
      colours = RColorBrewer::brewer.pal(name = "Blues", n = 7),
      breaks = scales,
      labels = as.character(scales),
      guide = "legend"
    ) +
    labs(fill = "", size = "") + 
    guides(fill = guide_legend(), size = guide_legend()) +
    theme +
    theme(legend.position = "right") + 
    NULL
}

# (p_scale_rain <- plot_scale(areas, data_rain_map, seq(5, 15, by = 2.5)) + 
(p_scale_rain <- plot_scale(areas, data_rain_map, seq(10, 30, by = 5)) + 
    ggtitle(paste0(
      "GPD scale parameter for rain (shape = ", 
      round(evgam_rain$predictions$shape, 2), 
      ")"
    )) + 
    theme(plot.title = element_text(hjust = 0.5)))

ggsave(
  "plots/gpd_rain_ire_weekly.png", p_scale_rain, width = 6.3, height = 6, units = "in"
)

#### Windspeeds ####

evgam_ws <- fit_evgam(data_ws, data_ws)

data_ws_map <- data_ws %>%
  mutate(scale = evgam_ws$predictions$scale) %>%
  st_as_sf(coords = c("lon", "lat"), remove = FALSE)
st_crs(data_ws_map) <- "WGS84"

# (p_scale_ws <- plot_scale(areas, data_ws_map, seq(0.6, 2, by = 0.1)) + 
(p_scale_ws <- plot_scale(areas, data_ws_map, seq(0.6, 1, by = 0.1)) + 
    ggtitle(paste0(
      "GPD scale parameter for windspeed (shape = ", 
      round(evgam_ws$predictions$shape, 2), 
      ")"
    )) + 
    theme(plot.title = element_text(hjust = 0.5)))

ggsave(
  "plots/gpd_ws_ire_weekly.png", p_scale_ws, width = 6.3, height = 6, units = "in"
)


#### Diagnostic plots (QQ, Probability) ####

pp_qq_plot <- \(evgam_fit, data_map) {
  # calculate exponential margins
  y_tilde <- (1 / evgam_fit$predictions$shape) *
    log(1 + evgam_fit$predictions$shape * (
            data_map$excess / evgam_fit$predictions$scale
          )
    )

  # calculate probability plot points
  # Plot estimated distribution function vs i / n + 1
  prob_plot_df <- data.frame(
    "x" = seq_len(nrow(data_map)) / (nrow(data_map) + 1),
    "y" = sort(1 - exp(-y_tilde))
  )
  # save for Kristina's code
  # readr::write_csv(prob_plot_df, "data/prob_plot_df.csv.gz")
  
  pp <- prob_plot_df %>%
    ggplot(aes(x = x, y = y)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    theme +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Empirical Probability", y = "Model Probability")
  
  
  # QQ plot
  qq_plot_df <- data.frame(
    "x" = -log(1 - (seq_len(nrow(data_map)) / (nrow(data_map) + 1))), # check this
    "y" = sort(y_tilde)
  )
  
  qq_p <- qq_plot_df %>%
    ggplot(aes(x = x, y = y)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    theme +
    scale_x_continuous(limits = c(0, 6)) +
    scale_y_continuous(limits = c(0, 6)) +
    labs(x = "Empirical Quantile", y = "Model Quantile")
  
  return(pp + qq_p)
}

pp_qq_plot(evgam_rain, data_rain_map)
pp_qq_plot(evgam_ws, data_ws_map)
# seem grand!


#### Return levels ####

# TODO: Calculate

#### Fit CE model ####

# Create "dummy" migpd object to fill in with evgam values
dat_mat <- data %>% 
  filter(name == data$name[[1]]) %>% 
  select(rain, wind_speed) %>% 
  as.matrix()
names(dat_mat) <- c("rain", "wind_speed")

# doesn't allow covariates in migpd, but does in evm
temp <- migpd(dat_mat, mqu = 0.95, penalty = "none")
# m <- evm(y = rain, data = data, qu = 0.95, penalty = "none", famuly = "gpd")
# m1 <- update(m, phi = ~lon + lat)

# fill in with GPD vals from evgam fits
# first, make predictions, unique to each location
pred_rain <- fit_evgam(data_rain, data)
pred_ws <- fit_evgam(data_ws, data)

# Pull scale and shape for rain and wind speed
x <- which(!duplicated(data$name))
data_gpd <- data[x, ] %>% 
  distinct(name, county, province, lon, lat) %>% 
  bind_cols(
    rename(pred_rain$predictions[x, ], scale_rain = scale, shape_rain = shape),
    rename(pred_ws$predictions[x, ],   scale_ws = scale,   shape_ws = shape)
  ) %>% 
  # also add number of exceedances over threshold for each
  left_join(
    # data %>% 
    #   filter(rain > thresh_rain) %>% 
    data_rain %>% 
      count(name, thresh) %>% 
      rename(n_rain = n, thresh_rain = thresh)
  ) %>% 
  left_join(
   # data %>% 
   #   filter(wind_speed > thresh_wind) %>% 
    data_ws %>% 
      count(name, thresh) %>% 
      rename(n_wind = n, thresh_wind = thresh)
  ) %>% 
  # fill in NAs (indicating no exceedances) with 0
  mutate(across(contains("n_"), \(x) ifelse(is.na(x), 0, x)))
 
# TODO: Look into this, double check everything is okay!
# create marginal `migpd` objects from `evgam` objects for each site
marginal <- lapply(seq_len(nrow(data_gpd)), \(i) {
  # initialise
  spec_marg <- temp
  # replace data 
  spec_marg$data <- data %>% 
    filter(name == data_gpd$name[i]) %>% 
    select(rain, wind_speed) %>% 
    as.matrix()
  names(spec_marg$data) <- c("rain", "wind_speed")
  # replace thresholds
  # spec_marg$models$rain$threshold <- thresh_rain
  # spec_marg$models$wind_speed$threshold <- thresh_wind
  spec_marg$models$rain$threshold <- data_gpd$thresh_rain[[i]]
  spec_marg$models$wind_speed$threshold <- data_gpd$thresh_wind[[i]]
  # replace coefficients
  spec_marg$models$rain$coefficients[1:2] <- c(
    # data_gpd$scale_rain[i], 
    log(data_gpd$scale_rain[i]),
    data_gpd$shape_rain[i]
  )
  spec_marg$models$wind_speed$coefficients[1:2] <- c(
    # data_gpd$scale_ws[i], 
    log(data_gpd$scale_ws[i]), 
    data_gpd$shape_ws[i]
  )
  
  return(spec_marg)
})
names(marginal) <- paste0(data_gpd$name, " - ", data_gpd$county)

# fit CE dependence model for each site
dependence <- lapply(marginal, \(x) {
  # fit for rain and wind speed
  vars <- c("rain", "wind_speed")
  ret <- lapply(vars, \(col) {
    mexDependence(x, which = col, dqu = 0.95)
  })
  names(ret) <- vars
  return(ret)
})

# pull a and b values for each site
a_b_vals <- lapply(dependence, \(x) {
  list(
    summary(x[[1]])$dependence[1:2],
    summary(x[[2]])$dependence[1:2]
  )
})
names(a_b_vals) <- data_gpd$name

# convert to dataframe
a_b_df <- bind_rows(lapply(seq_along(a_b_vals), \(i) {
  data.frame(
    "name" = names(a_b_vals)[i],
    "var"  = c("rain", "rain", "windspeed", "windspeed"),
    "parameter" = c("a", "b", "a", "b"),
    "value" = c(
      a_b_vals[[i]][[1]][1], a_b_vals[[i]][[1]][2],
      a_b_vals[[i]][[2]][1], a_b_vals[[i]][[2]][2]
    )
    # "a"    = c(a_b_vals[[i]][[1]][1], a_b_vals[[i]][[2]][1]),
    # "b"    = c(a_b_vals[[i]][[1]][2], a_b_vals[[i]][[2]][2])
  )
})) %>% 
  # TODO: Investigate NAs here
  filter(!is.na(name))

# join to previous data (not assigned?)
data_gpd %>% 
  left_join(
    a_b_df %>% 
      mutate(
        var = ifelse(var == "windspeed", "ws", var),
        col = paste(parameter, var, sep = "_")
      ) %>% 
      select(name, col, value) %>% 
      tidyr::pivot_wider(id_cols = name, names_from = col, values_from = value)
  )

# join in area information, change to sf
a_b_sf <- left_join(
  a_b_df, 
  distinct(data, name, lon, lat)
) %>% 
  st_as_sf(coords = c("lon", "lat"))
st_crs(a_b_sf) <- "WGS84"

# check that b values aren't exceptionally negative
sort(a_b_sf[a_b_sf$parameter == "b", ]$value)[1:10]

names <- c("a", "b")
breaks <- list(
  "a" = seq(-1, 1, by = 0.25),
  # "b" = c(seq(-22, 0, by = 2), 0.75)
  # "b" = c(seq(-6, 0, by = 0.5), 0.75)
  "b" = c(seq(-1, 1, by = 0.25))
)
p_lst <- lapply(seq_along(names), \(i) {
  ggplot() +
    geom_sf(data = areas, colour = "black") +
    geom_sf(
      data = a_b_sf %>% 
        filter(parameter == names[i]) %>% 
        # temp
        filter(!(parameter == "b" & value < -5)) %>% 
        identity(),
      # aes(fill = value, size = value),
      aes(fill = value, size = value),
      colour = "black",
      stroke = 1,
      pch = 21
    ) +
    facet_wrap(var ~ parameter, ncol = 2) + 
    scale_size_continuous(
      # breaks = breaks[[i]],
      # labels = as.character(breaks[[i]]), 
      guide = "legend"
    ) +
    # scale_fill_gradientn(
    #   # colours = RColorBrewer::brewer.pal(name = "Blues", n = 7),
    #   # breaks = seq(0, 0.5, by = 0.1),
    #   # labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5"),
    #   breaks = breaks[[i]],
    #   labels = as.character(breaks[[i]]),
    #   guide = "legend"
    # ) +
    scale_fill_gradient2(
      low = "blue3",
      high = "red3",
      na.value = "grey",
      # breaks = breaks[[i]], 
      # labels = as.character(breaks[[i]]),
      guide = "legend"
    ) + 
    guides(fill = guide_legend(), size = guide_legend()) +
    labs(fill = "", size = "") +
    theme
})

# p_lst[[1]] / p_lst[[2]]

# save
# ggsave("plots/a_ce_ire.png", p_lst[[1]], width = 6.3, height = 6, units = "in")
# ggsave("plots/b_ce_ire.png", p_lst[[2]], width = 6.3, height = 6, units = "in")
ggsave("plots/a_b_ce_ire_weekly.png", p_lst[[1]] / p_lst[[2]], width = 16, height = 12, units = "in")

#### Diagnosing problems in mexDependence ####

# look at points which are strange when Laplace transformed
# plot(dependence[[1]][[1]]$margins$data[which(dependence[[1]][[1]]$margins$transformed[, 1] > 2 & dependence[[1]][[1]]$margins$transformed[, 2] > -2), ])
# Why are all extreme points given same rain value???

#### Diagnostics for CE model ####

# easier to look at in a file

# generate plots
# Add county name to location, helps interpretation (done)
plots <- lapply(seq_along(dependence), \(i) { # each location
  lapply(seq_along(dependence[[i]]), \(j) { # rainfall and wind speed
    # generate plots, adding title to first one
    p <- ggplot(dependence[[i]][[j]], plot. = FALSE)
    p[[1]] <- p[[1]] + 
      ggtitle(paste(names(dependence)[i], names(dependence[[i]])[j], sep = " - "))
    return(p)
  })
})

# print plots
# lapply(plots, \(x) {
#   lapply(x, \(y) do.call("grid.arrange", y))
# })

# diagnostic plots for windspeed for every location (save to view)
# pdf("plots.pdf", onefile = TRUE)
# for (i in seq(length(dependence))) {
#   do.call("grid.arrange", ggplot(dependence[[i]][[2]]))  
# }
# dev.off()
pdf("plots/diag_plots_rain_weekly.pdf", onefile = TRUE) 
for (i in seq(length(plots))) {
  do.call("grid.arrange", plots[[i]][[1]])
}
dev.off()

pdf("plots/diag_plots_windspeed_weekly.pdf", onefile = TRUE) 
for (i in seq(length(plots))) {
  do.call("grid.arrange", plots[[i]][[2]])
}
dev.off()

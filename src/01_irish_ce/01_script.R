#### Conditional extremes modelling with evgam used for marginals ####

# following https://empslocal.ex.ac.uk/people/staff/by223/software/gpd.R
# In paper https://www.tandfonline.com/doi/full/10.1080/01621459.2018.1529596

# Fit varying threshold for evgam  (done)
# Plot exceedances over varying threshold as diagnostic (done)
# Model logscale (done)
# Use weekly data, as in KL diverenge clustering paper? (done)
# Try more smoothing for scale parameter, might lead to better results (done)
# - Didn't really..
# Model shape? As in https://www.tandfonline.com/doi/full/10.1080/01621459.2018.1529596 (done)
# - Seems typical to model shape as stationary over time, but not space
# TODO: Try different marginal thresholds for wind to improve QQ/prob plots
# TODO: Add uncertainty to QQ/probability plots

# TODO: Calculate errors around a, b values with bootstrap
# - May want to have any values deemed insig. different to 0 to be 0?

# Compare to "out-of-the-box" estimates from texmex (done)

# TODO: Calculate return levels?
# plot a vs alt etc, to see if there is any link to some predictors? (done)

# Outliers:
# TODO: Investigate 2 rows for Castleisland (Coom) in data_gpd, why not 1?
# TODO: Investigate large wind scale in Wexford, why is this the case? 
# - Do exploratory analysis on wind speeds (and rainfall while you're at it)

#### Libs ####

library(sf)
library(dplyr, quietly = TRUE)
library(tidyr)
library(ggplot2)
library(lubridate)
# library(texmex)
devtools::load_all("texmex") # forked texmex which allows for fixed b vals
library(gridExtra)
library(patchwork)
library(evgam)
library(geosphere)

source("src/functions.R")

theme <- theme_bw() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 16, hjust = 0.5),
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
  # "data/met_eireann/final/met_eir_wind.csv.gz",
  "data/met_eireann/final/met_eir_wind_alt.csv.gz",
  show_col_types = FALSE
) %>%
  filter(date %in% all_dates, !is.na(wind_speed))

# load shapefile
areas <- sf::read_sf("data/met_eireann/final/irl_shapefile.geojson")


#### Take Winter data only ####

data <- data %>% 
  mutate(month = as.numeric(substr(date, 6, 7))) %>% 
  filter(month %in% c(1:3, 10:12)) %>% 
  select(-month)

#### Convert to weekly ####

# function to find value corresponding to maximum density estimate for vector
max_density <- \(vec) {
  if (length(vec) == 1) return(vec)
  # return maximum kernel density estimate
  return(with(density(vec), x[which.max(y)]))
}

# Take weekly sum of precipitation and average of wind speed, to account for 
# lag in storm (as in Vignotto, Engelke 2021 study of GB + Ireland)
# TODO: Take weely averages of *daily* wind speed maxima!!
data <- data %>%
  mutate(week = floor_date(date, "week")) %>%
  group_by(week, name, county, province, lon, lat, alt) %>%
  summarise(
    rain       = sum(rain, na.rm = TRUE),
    wind_speed = mean(wind_speed, na.rm = TRUE),
    wind_dir   = max_density(wind_dir),
    .groups    = "drop"
  ) %>% 
  rename(date = week) # to agree with below code


#### Calculate distance to coast for each site ####

# TODO: replace calls to st_as_sf and st_crs with function
st_to_sf <- function(dat, coords = c("lon", "lat"), crs = "WGS84", ...) {
  out <- dat %>% 
    st_as_sf(coords = coords, ...)
  st_crs(out) <- crs
  return(out)
}

dist2coast <- \(dat, areas) {
  coastline <- st_union(areas)
  coastline <- st_cast(coastline, "MULTILINESTRING")
  
  # locations to find distance to coast for
  locs <- data %>% 
    distinct(name, lon, lat) %>% 
    st_to_sf()

  distances <- as.numeric(st_distance(locs, coastline))
  
  data %>% 
    left_join(data.frame("name" = locs$name, "dist2coast" = distances)) %>% 
    return()
}

data <- dist2coast(data, areas)


#### Data Exploration ####

# Appears to be rain outliers for many sites, investigate dates
# TODO: Add difference between second largest and largest rain reading, for context
# dat_outlier <- data %>% 
#   group_by(name, county) %>% 
#   filter(rain == max(rain)) %>% 
#   ungroup()
#   
# data %>% 
#   distinct(name, county, date, rain) %>% 
#   filter(rain != 0) %>% 
#   arrange(name, rain) %>% 
#   group_by(name) %>% 
#   mutate(rain_diff = lag(rain))

# Plot for one site (Wexford?)
# TODO: add threshold for Wexford?
data %>% 
  filter(name == "Wexford (Newtown W.w.)") %>% 
  ggplot(aes(x = rain, y = wind_speed)) + 
  geom_point() + 
  geom_vline(xintercept = 67.1)

# 2.7% of data kept, in line with quantile chosen
nrow(with(data, data[name == "Wexford (Newtown W.w.)" & rain > 67.1, ])) /
nrow(with(data, data[name == "Wexford (Newtown W.w.)", ]))

# as a test, remove rain over 100 for Wexford
# data <- data %>% 
#   filter(!(name == "Wexford (Newtown W.w.)" & rain > 100))
  

# TODO: Plot chi for each location (texmex has function for this)
# Follow Vignotto, Engelke code
# chiO3 <- chi(winter[, c("O3", "NO")])
# ggplot(chiO3, main =c("Chi" = "Chi: O3 and NO",
#                      "ChiBar" = "Chi-bar: O3 and NO"))

# TODO: Plot max weekly rainfall

# TODO: Plot max mean wind speed

#### Dealing with preferential sampling ####

## Weighting by Poisson point process density
# Assuming 'pp' is a point pattern object created from your data
# pp <- spatstat.geom::ppp(
#   x = data$lon, 
#   y = data$lat, 
#   window = spatstat.geom::owin(range(data$lon), range(data$lat))
# )
# 
# # Calculate density
# density <- density(pp)
# 
# # Assign weights inversely proportional to density
# weights <- 1 / density[pp]

## Kringing
# data2 <- data %>% 
#   mutate(rain = rain + rnorm(n(), 0, 1e-6))
# # Define the gstat object
# g <- gstat::gstat(formula = rain ~ date, locations = ~lon + lat, data = data2)
# 
# # Perform kriging
# kriging_result <- predict(g, newdata = data2)

#### Threshold selection ####

# function to fit varying threshold using quantile regression
# https://empslocal.ex.ac.uk/people/staff/by223/software/gpd.R
# TODO: Vary tau and see how that effects results (QQ plots, etc)
# quantile_thresh <- function(data, response, tau = .95) {
quantile_thresh <- function(data, response, tau = .95, jitter = TRUE) {
  fm_ald <- list(
    # response ~ s(lon, lat, k = 50), # location
    formula(paste(response, " ~ s(lon, lat, k = 50)")), # location
    ~ s(lon, lat, k = 40)                               # logscale
  )
  
  # jitter, if specified, to remove 0s when calculating quantiles
  if (jitter == TRUE) {
    data <- data %>%
      mutate(across(all_of(response), ~ . + rnorm(n(), 0, 1e-6)))
  }
  
  # fit the quantile regression model at tau'th percentile
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
  st_to_sf()

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
  # f <- list(excess ~ s(lon, lat, k = 12), ~ 1) # TODO: smooth for year?
  # f <- list(excess ~ s(lon, lat, k = 12) + s(date, k = 12), ~ 1) 
  # f <- list(excess ~ s(lon, lat, k = 50), ~ 1) # increase smoothing
  f <- list(
    excess ~ s(lon, lat, k = 40), # increase smoothing on scale parameter
    ~ s(lon, lat, k = 40) # shape parameter
  ) 
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
  mutate(
    scale = evgam_rain$predictions$scale,
    shape = evgam_rain$predictions$shape
  ) %>%
  st_to_sf(remove = FALSE)

# plot parameter values over space
# plot_scale <- \(areas, data_map, scales) {
# TODO: Add 
plot_param <- \(areas, data_map, scales, col = "scale") {
  col_sym <- sym(col)
  p <- ggplot(areas) +
    geom_sf(colour = "black", fill = NA) +
    geom_sf(
      data = data_map,
      # aes(fill = scale, size = scale),
      aes(fill = !!col_sym, size = !!col_sym),
      stroke = 1,
      pch = 21
    ) +
    # labs(fill = "Scale") +
    labs(fill = toupper(col)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_size_continuous(
      breaks = scales,
      labels = as.character(scales),
      guide = "legend"
    ) +
    labs(fill = "", size = "") + 
    guides(fill = guide_legend(), size = guide_legend()) +
    theme +
    theme(legend.position = "right")
  
  # add colour depending on range of values (same sign vs centred around 0)
  vals <- data_map[[col]]
  if (min(vals) < 0 && max(vals) > 0) {
    p <- p + 
      scale_fill_gradient2(
      low = "red3",
      high = "blue3",
      na.value = "grey",
      breaks = scales,
      labels = as.character(scales),
      guide = "legend"
    )
  } else {
    p <- p + 
      scale_fill_gradientn(
        colours = RColorBrewer::brewer.pal(name = "Blues", n = 7),
        breaks = scales,
        labels = as.character(scales),
        guide = "legend"
      )  
  }
  return(p)
}

# (p_scale_rain <- plot_param(areas, data_rain_map, seq(10, 30, by = 5)) + 
(p_scale_rain <- plot_param(areas, data_rain_map, seq(10, 40, by = 5)) + 
    ggtitle("GPD scale parameter for rain"))

# (p_shape_rain <- plot_param(areas, data_rain_map, seq(-0.05, 0.15, by = 0.04), "shape") + 
# (p_shape_rain <- plot_param(areas, data_rain_map, seq(-0.1, 0.15, by = 0.04), "shape") + 
(p_shape_rain <- plot_param(
  areas, data_rain_map, sort(c(0, c(0, seq(-0.1, 0.02, by = 0.03)))), "shape"
) + 
    ggtitle("GPD shape parameter for rain"))

ggsave(
  "plots/scale_rain_ire_weekly_winter.png", p_scale_rain, width = 6.3, height = 6, units = "in"
)
ggsave(
  "plots/shape_rain_ire_weekly_winter.png", p_shape_rain, width = 6.3, height = 6, units = "in"
)


#### Windspeeds ####

evgam_ws <- fit_evgam(data_ws, data_ws)

data_ws_map <- data_ws %>%
  # mutate(scale = evgam_ws$predictions$scale) %>%
  mutate(
    scale = evgam_ws$predictions$scale,
    shape = evgam_ws$predictions$shape
  ) %>%
  st_to_sf(remove = FALSE)

# (p_scale_ws <- plot_scale(areas, data_ws_map, seq(0.6, 2, by = 0.1)) + 
# (p_scale_ws <- plot_param(areas, data_ws_map, seq(0.6, 2, by = 0.2)) + 
#     ggtitle(paste0(
#       "GPD scale parameter for windspeed (shape = ", 
#       round(evgam_ws$predictions$shape, 2), 
#       ")"
#     )) + 
#     theme(plot.title = element_text(hjust = 0.5)))

(p_scale_ws <- plot_param(areas, data_ws_map, seq(0.6, 2, by = 0.2)) + 
    ggtitle("GPD scale parameter for wind speed"))

# (p_shape_ws <- plot_param(areas, data_ws_map, seq(-0.35, -0.1, by = 0.05), "shape") + 
(p_shape_ws <- plot_param(areas, data_ws_map, seq(-0.35, -0.1, by = 0.03), "shape") + 
    ggtitle("GPD shape parameter for wind speed"))

ggsave(
  "plots/scale_ws_ire_weekly_winter.png", p_scale_ws, width = 6.3, height = 6, units = "in"
)
ggsave(
  "plots/shape_ws_ire_weekly_winter.png", p_shape_ws, width = 6.3, height = 6, units = "in"
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

# TODO: Add uncertainty bounds! (See ITT2)
pp_qq_plot(evgam_rain, data_rain_map)
pp_qq_plot(evgam_ws, data_ws_map)


#### Return levels ####

# TODO: Calculate

#### Fit marginal model ####

# Create "dummy" migpd object to fill in with evgam values
dat_mat <- data %>% 
  filter(name == data$name[[1]]) %>% 
  select(rain, wind_speed) %>% 
  as.matrix()
names(dat_mat) <- c("rain", "wind_speed")

temp <- migpd(dat_mat, mqu = 0.95, penalty = "none")
# m <- evm(y = rain, data = data, qu = 0.95, penalty = "none", famuly = "gpd")
# m1 <- update(m, phi = ~lon + lat)

# fill in with GPD vals from evgam fits
# first, make predictions, unique to each location
pred_rain <- fit_evgam(data_rain, data)
pred_ws <- fit_evgam(data_ws, data)

# Pull scale and shape for rain and wind speed
# repeated for Dromod (Ruskey)
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
      mutate(thresh = round(thresh, 3)) %>% 
      count(name, thresh) %>% 
      rename(n_rain = n, thresh_rain = thresh)
  ) %>% 
  left_join(
   # data %>% 
   #   filter(wind_speed > thresh_wind) %>% 
    data_ws %>% 
      mutate(thresh = round(thresh, 3)) %>% 
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

#### Compare to "out-of-box" texmex predictions ####

# locs <- unique(data_gpd$name)
# texmex_marginal <- lapply(locs, \(loc) {
#   print(loc)
#   vals <- data %>% 
#     filter(name == loc) %>% 
#     select(rain, wind_speed)
#   threshs <- lapply(list(data_rain, data_ws), \(x) {
#     x %>% 
#       filter(name == loc) %>% 
#       pull(thresh) %>% 
#       first()
#   })
#   return(list(
#     "rain" = evm(y = vals[, 1, drop = TRUE], th = threshs[[1]]),
#     "wind_speed" = evm(y = vals[, 2, drop = TRUE], th = threshs[[2]])
#   ))
# })
# 
# # plot to compare to other values
# rain_coef <- lapply(lapply(texmex_marginal, `[[`, 1), \(x) x$coefficients)
# ws_coef <- lapply(lapply(texmex_marginal, `[[`, 2), \(x) x$coefficients)
# data_gpd_texmex <- data_gpd %>% 
#   mutate(
#     scale_rain_texmex = vapply(rain_coef, \(x) exp(x[[1]]), numeric(1)),
#     shape_rain_texmex = vapply(rain_coef, \(x) x[[2]], numeric(1)),
#     scale_ws_texmex = vapply(ws_coef, \(x) exp(x[[1]]), numeric(1)),
#     shape_ws_texmex = vapply(ws_coef, \(x) x[[2]], numeric(1))
#   )
# 
# comp_plots <- lapply(
#   list("scale_rain", "shape_rain", "scale_ws", "shape_ws"), 
#   \(col) {
#     col_texmex <- paste0(col, "_texmex")
#     x <- data_gpd_texmex %>% 
#       select(name, lon, lat, !!sym(col), !!sym(col_texmex)) %>% 
#       pivot_longer(cols = c(!!sym(col), !!sym(col_texmex)), names_to = "var") %>% 
#       st_to_sf()
#     
#     ggplot() +
#       geom_sf(data = areas, colour = "black") +
#       geom_sf(
#         data = x,
#         aes(fill = value, size = value),
#         colour = "black",
#         stroke = 1,
#         pch = 21
#       ) +
#       facet_wrap(~var) + 
#       theme + 
#       scale_fill_gradientn(
#         colours = RColorBrewer::brewer.pal(name = "Blues", n = 7),
#         na.value = "purple",
#         guide = "legend"
#       )
# })
# 
# pdf("plots/compare_evgam_texmex.pdf")
# do.call("grid.arrange", c(comp_plots, list(ncol = 2, nrow = 2)))
# dev.off()


#### Bootstrap sensitivity analysis on CE threshold selection ####

# plot(mexRangeFit(
#   marginal[[1]], "rain", quantiles = seq(0.6, 0.95, by = 0.05), trace = 11
# ))

# TODO: Repeat for fixed b, easily done through functionalising below!

# function to make bootstrap plots of CE parameters for each location and thresh 
create_bootstrap_plot <- \(
  marginal, quantiles = seq(0.6, 0.9, by = 0.05), fixed_b = FALSE
) {
  # very slow!
  vars <- c("rain", "wind_speed")
  bootstrap_thresh_plots <- lapply(seq_along(marginal), \(i) {
    print(names(marginal)[i])
    lapply(vars, \(col) {
      tryCatch({
        plots <- ggplot(
          mexRangeFit(
            marginal[[i]], 
            col, 
            # quantiles = seq(0.6, 0.9, by = 0.05), 
            quantiles = quantiles,
            fixed_b = fixed_b,
            trace = 11
        ), plot. = FALSE) # don't plot initially
        plots <- lapply(plots, \(x) {
          x + 
            ggtitle(paste(names(marginal)[i], "for", col)) + 
            theme
        })
        return(plots)
      }, error = function(cond) { 
        message(paste0("failed for ", names(marginal)[i], " cond. on ", col, "\n"))
        message(conditionMessage(cond)) 
        return(NA)
      })
    })
  })
  names(bootstrap_thresh_plots) <- names(marginal)  
  return(bootstrap_thresh_plots)
}

# function to save plots of bootstrapped CE parameter values
save_bootstrap_plot <- \(bootstrap_thresh_plots, file) {
  # pdf("plots/bootstrap_weekly_winter.pdf")
  pdf(file)
  lapply(bootstrap_thresh_plots, \(x) { # locations
    lapply(x, \(y) { # rain and wind
      if (length(y) == 1 && is.na(y)) return(NA)
      y[[2]] <- y[[2]] + ggtitle("")
      do.call("grid.arrange", c(y, list(ncol = 1)))
    })
  })
  dev.off()
}

# commented out as takes a long time to produce!
# boostrap_thresh_plots <- create_bootstrap_plot(marginal)
# # save object so it can be loaded again later, saves running above!
# saveRDS(bootstrap_thresh_plots, "bootrstrap_thresh_plots.RDS")
# save_bootstrap_plot(
#   boostrap_thresh_plots,
#   "plots/bootstrap_weekly_winter.pdf"
# )

# threshold of 0.7 seems appropriate for many?
# TODO: Inspect these plots and look for any irregularities

# bootstrap_thresh_plots_fixed_b <- create_bootstrap_plot(marginal, fixed_b = TRUE)
# saveRDS(bootstrap_thresh_plots_fixed_b, "bootrstrap_thresh_plots_fixed_b.RDS")
# save_bootstrap_plot(
#   bootstrap_thresh_plots_fixed_b, 
#   "plots/bootstrap_weekly_winter_fixed_b.pdf"
# )


#### Fit CE model ####

# fit CE dependence model for each site
# Save profile likelihood plots as well
# dependence <- lapply(marginal, \(x) {
# pdf("plots/prof_lik_weekly_winter.pdf", onefile = TRUE)
vars <- c("rain", "wind_speed")
dependence <- lapply(seq_along(marginal), \(i) {
  # print(i)
  # fit for rain and wind speed
  ret <- lapply(vars, \(col) {
    # TODO: Can you supply threshold yourself? or match to thresholding value?
    # TODO: Plot prifle likelihood of dependence model parameters?
    mexDependence(
      # x, 
      marginal[[i]],
      which = col, 
      dqu = 0.7, # seems a reasonable choice of threshold given bootstrap plots
      nOptim = 2, 
      # PlotLikDo = TRUE, 
      # PlotLikTitle = paste0("for ", col, " - ", names(marginal)[i])
      fixed_b = TRUE
    )
  })
  names(ret) <- vars
  return(ret)
})
names(dependence) <- names(marginal)
# dev.off()


#### Bootstrap parameter estimates and calculate quantiles ####

# bootmex()
# mymex <- mex(winter , mqu = .7, dqu = .7, which = "NO")
# myboot <- bootmex(mymex)
# myboot
# plot(myboot,plots="gpd")
# plot(myboot,plots="dependence")

# fails for i <- 3!
# final_bootstraps <- lapply(dependence, \(spec_loc) {
# final_bootstraps <- lapply(seq_along(dependence), \(i) {
#   print(names(dependence)[[i]])
#   spec_loc <- dependence[[i]]
#   lapply(spec_loc, \(spec_loc_var) {
#     boots <- suppressWarnings(suppressMessages(bootmex(spec_loc_var)))
#     # pulled from texmex:::plot.bootmex
#     d2 <- dim(boots$boot[[1]][[2]]) # "dependence"
#     margins <- boots$margins
#     x <- boots$boot
#     co <- unlist(lapply(x, function(z, wh) z[[wh]], wh = which))
#     co <- array(co, dim = c(d2[1], d2[2], length(co)/prod(d2)))
#     lco <- list(length = prod(d2))
#     for (i in 1:d2[2]) {
#       for (j in 1:d2[1]) {
#         lco[[j + d2[1] * (i - 1)]] <- co[j, i, ]
#       }
#     }
#     
#     # a and b values from bootstrap
#     vals <- list(
#       "a" = lco[[1]], 
#       "b" = lco[[2]]
#     )
#     # print(lapply(vals, head))
#     # plot(a_vals, b_vals, main = "Test", xlab = "a", ylab = "b")
#     # test if 0 in 95% CI (i.e. signs are different for 0.025, 0.975 quantiles)
#     tests <- lapply(vals, \(x) {
#       val <- prod(quantile(x, c(0.0025, 0.975), na.rm = TRUE))
#       return(list(val, val < 0))
#     })
#     return(list("values" = vals, "tests" = tests))
#   })
# })

# boots <- bootmex(dependence[[1]][[1]])
# # pulled from texmex:::plot.bootmex
# d2 <- dim(boots$boot[[1]][[2]]) # "dependence"
# margins <- boots$margins
# x <- boots$boot
# co <- unlist(lapply(x, function(z, wh) z[[wh]], wh = which))
# co <- array(co, dim = c(d2[1], d2[2], length(co)/prod(d2)))
# lco <- list(length = prod(d2))
# for (i in 1:d2[2]) {
#   for (j in 1:d2[1]) {
#     lco[[j + d2[1] * (i - 1)]] <- co[j, i, ]
#   }
# }
# 
# vals <- list(
#   "a" = lco[[1]], 
#   "b" = lco[[2]]
# )
# # plot(a_vals, b_vals, main = "Test", xlab = "a", ylab = "b")
# # test if 0 in 95% CI (i.e. signs are different for 0.025, 0.975 quantiles)
# tests <- lapply(vals, \(x) {
#   val <- prod(quantile(x, c(0.0025, 0.975)))
#   return(list(val, val < 0))
# })
# 
# prod(quantile(a_vals, c(0.0025, 0.975)))
# prod(quantile(b_vals, c(0.0025, 0.975)))


#### Plots ####

# pull a and b values for each site
ab_vals <- lapply(dependence, \(x) {
  list(
    summary(x[[1]])$dependence[1:2],
    summary(x[[2]])$dependence[1:2]
  )
})
names(ab_vals) <- data_gpd$name

# convert to dataframe
ab_df <- bind_rows(lapply(seq_along(ab_vals), \(i) {
  data.frame(
    "name" = names(ab_vals)[i],
    "var"  = c("rain", "rain", "windspeed", "windspeed"),
    "parameter" = c("a", "b", "a", "b"),
    "value" = c(
      ab_vals[[i]][[1]][1], ab_vals[[i]][[1]][2],
      ab_vals[[i]][[2]][1], ab_vals[[i]][[2]][2]
    )
    # "a"    = c(ab_vals[[i]][[1]][1], ab_vals[[i]][[2]][1]),
    # "b"    = c(ab_vals[[i]][[1]][2], ab_vals[[i]][[2]][2])
  )
})) %>% 
  # TODO: Investigate NAs here
  filter(!is.na(name))


# join a and b values to elevation and distance to coast
ab_df <- ab_df %>% 
  # left_join(distinct(data, name, county, dist2coast, alt, lat, lon))
  left_join(data %>% 
    group_by(name, county, dist2coast, alt, lat, lon) %>% 
    summarise(wind_dir = max_density(wind_dir), .groups = "drop")
  )

# plot a versus other variables 
# TODO: Add wind direction here, could be interesting!
# TODO: Colour points by county? Or coastal vs inland?
# TODO: Could also colour points for wind speed ~1 as red and ignore in lm 
# (or add additional fitted line for including these points)
pred_p <- ab_df %>% 
  filter(parameter == "a") %>% 
  filter(value < 0.99) %>% 
  pivot_longer(
    # c("alt", "dist2coast"), 
    # c("alt", "dist2coast", "lat"),
    c("alt", "dist2coast", "lat", "lon"), 
    names_to = "predictor", 
    values_to = "pred_val"
  ) %>% 
  ggplot(aes(x = value, y = pred_val)) +
  # geom_point(aes(colour = county)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  ggpubr::stat_cor(cor.coef.name = "rho", method = "pearson") + 
  # facet_grid(var ~ predictor, scales = "free") + 
  facet_wrap(var ~ predictor, scales = "free_y", nrow = 2) + 
  labs(x = "a", y = "value") + 
  theme
pred_p

# TODO: make notes
ggsave("plots/a_val_corrs_weekly_winter.png", pred_p, width = 16, height = 12, units = "in")

ab_df_wide <- ab_df %>% 
  filter(parameter == "a") %>% 
  select(name, county, var, value) %>% 
  pivot_wider(names_from = var, values_from = value)

# just plot a for rain vs a for wind!
# TODO: Lots can be done here
# Add county centroids to plot (done)
# TODO: Add fitted line with and without windspeed high a values 
ab_df_wide %>% 
  # much more highly correlated when removing large wind speeds values for a!
  filter(windspeed < 0.99) %>% 
  # take county wide mean centroid
  group_by(county) %>% 
  # summarise(
  mutate(
    across(c(rain, windspeed), mean, .names = "{.col}_mean") # , .groups = "drop"
  ) %>% 
  ungroup() %>% 
  ggplot(aes(x = rain, y = windspeed)) + 
  # county centroid points
  geom_point(
    aes(x = rain_mean, y = windspeed_mean, colour = county), 
    size = 6, 
    pch = 10, 
    alpha = 0.9
  ) + 
  # site points
  geom_point(aes(colour = county), size = 3, alpha = 0.9) + 
  geom_smooth(colour = "blue", method = "lm") + 
  # doesn't make much difference
  # geom_smooth(aes(x = rain_mean, y = windspeed_mean), colour = "red", method = "lm") + 
  ggpubr::stat_cor(cor.coef.name = "rho", method = "pearson") + 
  scale_x_continuous(limits = c(0, 0.65), expand = c(0, 0)) +
  # scale_y_continuous(limits = c(0, 1), expand = c(0, 0.1)) + 
  scale_y_continuous(limits = c(0, 0.8), expand = c(0, 0)) +
  # scale_y_continuous(limits = c(0, 0.6), expand = c(0, 0)) + # for counties
  theme

# Wind rose and polar plots for plotting a vs wind direction
# TODO: Functionalise to use with wind as well

# Bin wind direction
ab_w_dir <- ab_df_wide %>% 
  # TODO: How to add wind direction??
  left_join(distinct(data, name, wind_dir))

ab_w_dir$direction_bin <- cut(
  ab_w_dir$wind_dir, breaks = seq(0, 360, by = 15), include.lowest = TRUE
)

# function to create wind rose and polar plots
circ_plot <- \(data, response) {
  # Calculate average a for each direction bin
  # summary_a <- aggregate(rain ~ direction_bin, data = ab_w_dir, FUN = mean)
  summary_a <- aggregate(
    formula(paste(response, "~ direction_bin")), 
    data = data, 
    FUN = mean
  )
  
  # Convert direction_bin to numeric values (bin midpoint for better plotting)
  midpoints <- seq(7.5, 352.5, by = 15)  # Midpoints of the bins
  summary_a$midpoint <- midpoints[summary_a$direction_bin]
  
  p1 <- ggplot(
    summary_a, 
    aes(x = factor(midpoint), y = .data[[response]])# , fill = factor(direction_bin))
  ) +
    geom_bar(stat = "identity", color = "white") +
    coord_polar(start = -pi/2) +
    # coord_polar() + 
    scale_x_discrete(
      # name = "Wind Direction (degrees)", 
      name = "test",
      labels = function(x) sprintf("%.0s", x)
    ) +
    labs(y = "a", title = paste0("Wind Rose Plot for ", response)) +
    theme +
    theme(axis.text.x = element_text(angle = 90))
  
  # polar plot
  p2 <- ab_w_dir %>% 
    group_by(name, .data[[response]]) %>% 
    summarise(wind_dir = max_density(wind_dir), .groups = "drop") %>% 
    ggplot(aes(x = wind_dir, y = .data[[response]])) + 
    geom_point(color = "blue", alpha = 0.8) +
    coord_polar(start = -pi/2) +
    labs(
      x = "Wind Direction (degrees)", 
      y = "a", 
      title = paste0("Polar Plot of a for ", response)
    ) +
    theme()
  
  return(list(p1, p2))
}

circ_plot(ab_w_dir, "rain")
circ_plot(ab_w_dir, "windspeed")

# join to previous data (not assigned?)
# data_gpd %>% 
#   left_join(
#     ab_df %>% 
#       mutate(
#         var = ifelse(var == "windspeed", "ws", var),
#         col = paste(parameter, var, sep = "_")
#       ) %>% 
#       select(name, col, value) %>% 
#       tidyr::pivot_wider(id_cols = name, names_from = col, values_from = value)
#   )

# change to sf
ab_sf <- st_to_sf(ab_sf)
# check that b values aren't exceptionally negative (<-1)
sort(ab_sf[ab_sf$parameter == "b", ]$value)[1:10]

names <- c("a", "b")
# TODO: Find out how to force include breaks! See circ survey figure
breaks <- list(
  # "a" = seq(-1, 1, by = 0.25),
  # "a" = sort(c(seq(-1, 1, by = 0.2), -0.5)),
  "a" = seq(-1, 1, by = 0.2),
  # "b" = c(seq(-22, 0, by = 2), 0.75)
  # "b" = c(seq(-6, 0, by = 0.5), 0.75)
  # "b" = c(seq(-1, 1, by = 0.25))
  "b" = seq(-2, 1, by = 0.1)
)
p_lst <- lapply(seq_along(names), \(i) {
  ggplot() +
    geom_sf(data = areas, colour = "black") +
    geom_sf(
      data = ab_sf %>% 
        filter(parameter == names[i]) %>% 
        # temp
        # filter(!(parameter == "b" & value < -5)) %>% 
        identity(),
      # aes(fill = value, size = value),
      aes(fill = value, size = value),
      colour = "black",
      stroke = 1,
      pch = 21
    ) +
    facet_wrap(var ~ parameter, ncol = 2) + 
    scale_size_continuous(
      breaks = breaks[[i]],
      labels = as.character(breaks[[i]]), 
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
      breaks = breaks[[i]], 
      labels = as.character(breaks[[i]]),
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
# ggsave("plots/a_b_ce_ire_weekly_winter.png", p_lst[[1]] / p_lst[[2]], width = 16, height = 12, units = "in")

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

# 53 looks odd!

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
# pdf("plots/diag_plots_rain_weekly_winter.pdf", onefile = TRUE) 
# for (i in seq(length(plots))) {
#   do.call("grid.arrange", plots[[i]][[1]])
# }
# dev.off()
# 
# pdf("plots/diag_plots_windspeed_weekly_winter.pdf", onefile = TRUE) 
# for (i in seq(length(plots))) {
#   do.call("grid.arrange", plots[[i]][[2]])
# }
# dev.off()

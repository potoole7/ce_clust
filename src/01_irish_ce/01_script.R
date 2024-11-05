#### Conditional extremes modelling with evgam used for marginals ####

# following https://empslocal.ex.ac.uk/people/staff/by223/software/gpd.R
# In paper https://www.tandfonline.com/doi/full/10.1080/01621459.2018.1529596

# TODO: make master plotting function for Ireland, lots of code is repeated!

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
library(latex2exp)
library(ggpattern)

source("src/functions.R")
source("src/01_irish_ce/functions.R")

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

# function to remove " - " (county name) from names
rm_cnty <- \(x) {
  vapply(stringr::str_split(x, " - "), `[[`, 1, FUN.VALUE = character(1))
}


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

# Remove locations with invalid CE fits
# loc_nas <- c(
#   "Cloonacool (Lough Easkey)", "Costelloe Fishery", "Crolly (Filter Works)", 
#   "Kilcar (Cronasillagh)", "Kinsale (Vocational School)", "Sneem", 
#   "Wexford (Newtown W.w.)"
# )
# data <- filter(data, !name %in% loc_nas)


#### Take Winter data only, remove weeks with no rain ####

data <- data %>% 
  mutate(month = as.numeric(substr(date, 6, 7))) %>% 
  filter(month %in% c(1:3, 10:12)) %>% 
  dplyr::select(-month) %>% 
  # remove weeks with no rainfall, as in Vignotto 2021 (required/important?)
  filter(rain != 0) %>% 
  identity()


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

data <- dist2coast(data, areas)
# save
readr::write_csv(data, "data/met_eireann/final/met_eir_preprocess.csv.gz")

#### Motivating example plots ####

# locations with highest and lowest rain, to specifically plot
highest_lowest_rain <- data %>% 
  group_by(name) %>% 
  summarise(rain = mean(rain, na.rm = TRUE), .groups = "drop") %>% 
  arrange(rain) %>% 
  slice(c(1, n())) %>% 
  pull(name)

data_plot <- data %>% 
  # mutate(indicator = ifelse(name %in% highest_lowest_rain, name, NA)) %>% 
  mutate(indicator = ifelse(name %in% highest_lowest_rain, name, "other")) %>% 
  arrange(desc(indicator)) %>% 
  st_to_sf()

# First, plot the location of each site
p1 <- ggplot(areas) +
  geom_sf(colour = "black", fill = NA) +
  # points other than two sites
  geom_sf(
    data = filter(data_plot, indicator == "other"),
    colour = "black",
    size = 3,
    alpha = 0.9,
    show.legend = FALSE
  ) +
  # points for two sites with highest and lowest rain
  geom_sf(
    data = filter(data_plot, indicator != "other"),
    # aes(colour = indicator, size = indicator),
    colour = "black",
    alpha = 0.9,
    # show.legend = TRUE
    show.legend = FALSE
  ) +
  scale_colour_manual(values = c(ggsci::pal_nejm()(2))) + 
  scale_size_manual(values = c(4.5, 4.5)) + 
  labs(colour = "", size = "") + 
  theme + 
  # remove axis text
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.key = element_blank()
  )

cols <- ggsci::pal_nejm()(4)
cols[2] <- "black"
p21 <- data_plot %>% 
  filter(name == highest_lowest_rain[2], rain > 0) %>% 
  group_by(name) %>% 
  mutate(across(c(rain, wind_speed), ~ quantile(.x, 0.95, na.rm = TRUE), .names = "quant_{.col}")) %>% 
  ungroup() %>% 
  mutate(col = case_when(
    rain > quant_rain & wind_speed > quant_wind_speed ~ "Both",
    rain > quant_rain & wind_speed <= quant_wind_speed ~ "Rain",
    rain <= quant_rain & wind_speed > quant_wind_speed ~ "Wind",
    TRUE ~ "Neither"
  )) %>% 
  ggplot(aes(x = rain, y = wind_speed)) + 
  # geom_point(aes(colour = name), size = 1.5, alpha = 0.9) + 
  geom_point(aes(colour = col), size = 1.5, alpha = 0.9) + 
  geom_vline(aes(xintercept = quant_rain), linetype = "dashed") + 
  geom_hline(aes(yintercept = quant_wind_speed), linetype = "dashed") + 
  # facet_wrap(~ name, scales = "free_x") + 
  scale_colour_manual(values = cols) + 
  labs(
    # x = "Weekly total precipitation (mm)", 
    x = "precipitation (mm)", 
    y = "wind speed (m/s)", # TODO: What is the unit of ws?
    colour = ""
  ) + 
  guides(colour = "none") + 
  theme + 
  # remove facet labels, colour will do
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(), 
    legend.key = element_blank()
  )  

# Second, plot wind speeds against rain for sites with the highest and lowest rainfall
# TODO: Add 95% quantile lines for both (?)
p2 <- data_plot %>% 
  filter(name %in% highest_lowest_rain, rain > 0) %>% 
  group_by(name) %>% 
  mutate(across(c(rain, wind_speed), ~ quantile(.x, 0.95, na.rm = TRUE), .names = "quant_{.col}")) %>% 
  ungroup() %>% 
  ggplot(aes(x = rain, y = wind_speed)) + 
  # geom_point(aes(colour = name), size = 1.5, alpha = 0.9) + 
  geom_point(aes(colour = name), size = 1.5, alpha = 0.9) + 
  # geom_vline(aes(xintercept = quant_rain)) + 
  # geom_hline(aes(yintercept = quant_wind_speed)) + 
  facet_wrap(~ name, scales = "free_x") + 
  scale_colour_manual(values = c(ggsci::pal_nejm()(2))) + 
  labs(
    # x = "Weekly total precipitation (mm)", 
    x = "precipitation (mm)", 
    y = "wind speed (m/s)", # TODO: What is the unit of ws?
    colour = ""
  ) + 
  theme + 
  # remove facet labels, colour will do
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(), 
    legend.key = element_blank()
  )

# join plots
# TODO: Change size of first plot to be larger!
p_sec_2 <- p1 + 
  (p2 + guides(colour = "none", size = "none")) + 
  # have common legends
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
  
ggsave("latex/plots/01_ire_plot.png", p1, width = 7, height = 6, units = "in")
ggsave("latex/plots/021_mv_extreme.png", p21, width = 7, height = 6, units = "in")
ggsave("latex/plots/02_mot_ex_plot.png", p_sec_2, width = 7, height = 6, units = "in")


#### Data Exploration ####

# How many sites are in the data? How many rows? 

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

data %>% 
  filter(name == data$name[[1]]) %>% 
  ggplot(aes(x = rain, y = wind_speed)) + 
  geom_point() + 
  geom_vline(xintercept = 67.1)


# 2.7% of data kept, in line with quantile chosen
nrow(with(data, data[name == "Wexford (Newtown W.w.)" & rain > 67.1, ])) /
nrow(with(data, data[name == "Wexford (Newtown W.w.)", ]))

# as a test, remove rain over 100 for Wexford
# data <- data %>% 
#   filter(!(name == "Wexford (Newtown W.w.)" & rain > 90))
  

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

# Threshold rainfall and wind speed with varying threshold per location
data_rain <- quantile_thresh(data, "rain")
data_ws <- quantile_thresh(data, "wind_speed")

# data to plot u across space
data_thresh <- data_rain %>% 
  distinct(name, lon, lat, thresh_rain = thresh) %>% 
  left_join(
    data_ws %>% 
    distinct(name, lon, lat, thresh_ws = thresh)
  ) %>% 
  pivot_longer(contains("thresh"), names_to = "var") %>% 
  st_to_sf()

# function for plotting
plot_thresh <- \(x, areas, scales = NULL) {
  ggplot(areas) + 
    geom_sf(colour = "black", fill = NA) +
    geom_sf(
      data = x,
      # aes(fill = n, size = n),
      aes(fill = value, size = value),
      # size = 6,
      stroke = 1,
      pch = 21
    ) +
    # labs(fill = toupper(col)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    # scale_fill_gradientn(
    #   colours = RColorBrewer::brewer.pal(name = "Blues", n = 7),
    #   breaks = scales,
    #   labels = as.character(scales),
    #   guide = "legend"
    # ) +
    viridis::scale_fill_viridis(
      option = "A",
      breaks = scales,
      labels = as.character(scales),
      guide = "legend"
    ) + 
    scale_size_continuous(
      range = c(3, 8),
      breaks = scales,
      labels = as.character(scales),
      guide = "legend"
    ) +
    labs(fill = "", size = "") + 
    guides(fill = guide_legend(), size = guide_legend()) +
    theme +
    theme(
      legend.position = "right",
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.key = element_blank()
    )
}

p_thresh_rain <- plot_thresh(
  filter(data_thresh, var == "thresh_rain"), areas, seq(50, 150, by = 25)
) +
  labs(fill = "u for precipitation", size = "u for precipitation")

p_thresh_ws <- plot_thresh(
  # filter(data_thresh, var == "thresh_ws"), areas, seq(9, 15, by = 2)
  filter(data_thresh, var == "thresh_ws"), areas, seq(8, 11, by = 0.5)
) +
  labs(fill = "u for wind speed", size = "u for wind speed")

p_thresh <- p_thresh_rain | p_thresh_ws
ggsave("latex/plots/031_thresh_plots.png", p_thresh, width = 6.3, height = 6, units = "in")

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
  # labs(fill = "", size = "") +
  # guides(fill = guide_legend(), size = guide_legend()) +
  theme +
  theme(legend.position = "right") +
  NULL


#### Fit evgam to rain #### 

# Fit evgam model with spline for spatial location, create predictions
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
plot_param <- \(data_map, areas, scales, col = "scale") {
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
      range = c(2, 7),
      breaks = scales,
      labels = as.character(scales),
      guide = "legend"
    ) +
    labs(fill = "", size = "") + 
    guides(fill = guide_legend(), size = guide_legend()) +
    theme +
    theme(
      legend.position = "right",
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.key = element_blank()
    )
  
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
      # scale_fill_gradientn(
      #   colours = RColorBrewer::brewer.pal(name = "Blues", n = 7),
      #   breaks = scales,
      #   labels = as.character(scales),
      #   guide = "legend"
      # )  
      viridis::scale_fill_viridis(
        breaks = scales,
        labels = as.character(scales),
        guide = "legend" 
      )
    # p <- p + 
    #   scale_fill_gradient2(
    #   low = "white",
    #   high = "blue3",
    #   na.value = "grey",
    #   breaks = scales,
    #   labels = as.character(scales),
    #   guide = "legend"
    # )
  }
  return(p)
}

# (p_scale_rain <- plot_param(areas, data_rain_map, seq(10, 30, by = 5)) + 
# (p_scale_rain <- plot_param(areas, data_rain_map, seq(10, 40, by = 5)) + 
# lab_scale_rain <- TeX(r"($\sigma_{rain}$)")
lab_scale_rain <- "scale"
p_scale_rain <- plot_param(data_rain_map, areas, seq(10, 45, by = 5)) + 
    # ggtitle("GPD scale parameter for rain"))
    # ggtitle("scale"))
    labs(fill = lab_scale_rain, size = lab_scale_rain) + 
    NULL

# (p_shape_rain <- plot_param(areas, data_rain_map, seq(-0.05, 0.15, by = 0.04), "shape") + 
# (p_shape_rain <- plot_param(areas, data_rain_map, seq(-0.1, 0.15, by = 0.04), "shape") + 
# lab_shape_rain <- TeX(r"($\xi_{rain}$)")
lab_shape_rain <- "shape"
p_shape_rain <- plot_param(
  # areas, data_rain_map, sort(c(0, c(0, seq(-0.1, 0.02, by = 0.03)))), "shape"
  data_rain_map, areas, round(seq(-0.15, 0.15, by = 0.05), 2), "shape"
) + 
    labs(fill = lab_shape_rain, size = lab_shape_rain) + 
    NULL
    # theme(legend.title = element_text(size = 1))

(p_mod_rain <- p_scale_rain + p_shape_rain)

ggsave(
  "latex/plots/032_gpd_rain.png", p_mod_rain, width = 6.3, height = 6, units = "in"
)

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

lab_scale_ws <- TeX(r"($\sigma_{wind}$)")
lab_scale_ws <- "scale"
# p_scale_ws <- plot_param(data_ws_map, areas, seq(0.5, 1, by = 0.1)) + 
p_scale_ws <- plot_param(data_ws_map, areas, seq(0.5, 0.8, by = 0.05)) + 
    # ggtitle("GPD scale parameter for wind speed")
    labs(fill = lab_scale_ws, size = lab_scale_ws) + 
    # theme(legend.title = element_text(size = 15))
    NULL
 
# (p_shape_ws <- plot_param(areas, data_ws_map, seq(-0.35, -0.1, by = 0.05), "shape") + 
# give same colour scheme as rain
lab_shape_ws <- TeX(r"($\xi_{wind}$)")
lab_shape_ws <- "shape"
p_shape_ws <- plot_param(data_ws_map, areas, seq(-0.3, -0.1, by = 0.05), "shape") + 
    # ggtitle("GPD shape parameter for wind speed") + 
    labs(fill = lab_shape_ws, size = lab_shape_ws) + 
    # theme(legend.title = element_text(size = 15)) + 
    scale_fill_gradient2(
      low = "red3",
      high = "blue3",
      na.value = "grey",
      breaks = seq(-0.3, -0.1, by = 0.05),
      labels = as.character(seq(-0.3, -0.1, by = 0.05)),
      guide = "legend"
    )

p_mod_ws <- p_scale_ws + p_shape_ws

ggsave(
  "latex/plots/033_gpd_ws.png", p_mod_ws, width = 6.3, height = 6, units = "in"
)

ggsave(
  "plots/scale_ws_ire_weekly_winter.png", p_scale_ws, width = 6.3, height = 6, units = "in"
)
ggsave(
  "plots/shape_ws_ire_weekly_winter.png", p_shape_ws, width = 6.3, height = 6, units = "in"
)

#### Finite end points for wind speed ####

scales <- seq(10, 25, by = 1)
p_endpoint <- data_ws_map %>% 
  mutate(end_point = thresh - (scale / shape)) %>% 
  plot_param(areas, scales, "end_point") +
  scale_size_continuous(
    range = c(1.5, 10),
    breaks = scales,
    labels = as.character(scales),
    guide = "legend"
  ) + 
  labs(fill = "endpoint (m/s)", size = "endpoint (m/s)")

ggsave(
  "latex/plots/034_ws_endpoint.png", p_endpoint, width = 6.3, height = 6, units = "in"
)


#### Diagnostic plots (QQ, Probability) ####

# save objects
sf::write_sf(data_rain_map, "data_rain_map.geojson")
sf::write_sf(data_ws_map, "data_ws_map.geojson")

pp_qq_plot <- \(evgam_fit, data_map, xlab = "Empirical") {
  
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
    labs(x = xlab, y = "Model")
  
  
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
    labs(x = xlab, y = "")
  
  return(pp + qq_p)
}

# TODO: Add uncertainty bounds! (See ITT2)
pp_qq_rain <- pp_qq_plot(evgam_rain, data_rain_map)
pp_qq_ws <- pp_qq_plot(evgam_ws, data_ws_map)

pp_qq <- pp_qq_rain / pp_qq_ws

ggsave(
  "latex/plots/035_pp_qq.png", pp_qq, width = 6.3, height = 6, units = "in"
)

ggsave(
  "plots/pp_qq_rain_winter.png", pp_qq_rain, width = 6.3, height = 6, units = "in"
)
ggsave(
  "plots/pp_qq_wind_winter.png", pp_qq_ws, width = 6.3, height = 6, units = "in"
)

#### Uncertainty in xi ####



#### Return levels ####

# TODO: Calculate

#### Fit marginal model ####

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
marginal <- gen_marg_migpd(data_gpd, data)
names(marginal) <- paste0(data_gpd$name, " - ", data_gpd$county)


#### Compute Xi statistics for all ####

# calculate chi for one location
# chi_test <- data %>% 
#   # filter(name == data$name[1]) %>% 
#   filter(name == highest_lowest_rain[1]) %>% 
#   dplyr::select(rain, wind_speed) %>% 
#   chi() 
#   
# chi_test %>%   
#   ggplot(main = c("ChiBar" = "", "Chi" = ""))

# for each location, pull out 95th quantile for chibar
# Colour chi grey if not valid
# Plot on map
names <- unique(data$name)
chi_95_df <- bind_rows(lapply(names, \(x) {
  print(paste0(
    round(100 * which(names == x) / length(names), 2), "% completed"
  ))
  chi <- data %>% 
    filter(name == x) %>% 
    dplyr::select(rain, wind_speed) %>% 
    chi()
  
  # whether to show chi or not, based on whether chibar upper extend crosses 1
  show_chi <- !prod(tail(chi$chibar[, 3]) < 1)
  
  loc <- which.min(abs(chi$quantile - 0.95))
  return(data.frame(
    "name" = x,
    "chi" = chi$chi[loc, 2, drop = TRUE], 
    "chibar" = chi$chibar[loc, 2, drop = TRUE], 
    "show_chi" = show_chi
  ))
}))
rownames(chi_95_df) <- NULL

# join in area statistics
chi_95_sf <- chi_95_df %>% 
  pivot_longer(c("chi", "chibar"), names_to = "var") %>% 
  # always show chibar plot
  mutate(show_chi = ifelse(value == "chibar", TRUE, show_chi)) %>% 
  left_join(
  distinct(data, name, lon, lat)
) %>% 
  st_to_sf()

# plot for chibar (haven't time to functionalise!)
scales_chibar <- seq(-0.1, 0.6, by = 0.1)
# lab_chibar <- TeX(r"($\bar{\chi}(0.95)$)")
# lab_chibar <-  expression(bar(chi)(u))
lab_chibar <-  "chibar(0.95)"
point_ranges <- c(2, 6)
chibar_p <- ggplot(areas) +
    geom_sf(colour = "black", fill = NA) +
    geom_sf(
      data = filter(chi_95_sf, var == "chibar"),
      aes(fill = value, size = value),
      stroke = 1,
      pch = 21
    ) +
    # facet_wrap(~var) + 
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    # viridis::scale_fill_viridis(
    #     option = "F",
    #     breaks = scales_chibar,
    #     labels = as.character(scales_chibar),
    #     guide = "legend"
    # ) + 
    scale_fill_gradientn(
      # colours = RColorBrewer::brewer.pal(name = "Blues", n = 7),
      # colours = wesanderson::wes_palette("Zissou1", 8, type = "continuous"), 
      colours = rev(heat.colors(7)),
      breaks = scales_chibar,
      labels = as.character(scales_chibar),
      guide = "legend"
    ) + 
    scale_size_continuous(
      range = point_ranges,
      breaks = scales_chibar,
      labels = as.character(scales_chibar),
      guide = "legend"
    ) +
    labs(fill = lab_chibar, size = lab_chibar) + 
    guides(fill = guide_legend(), size = guide_legend()) +
    theme +
    theme(
      legend.position = "right",
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.key = element_blank()
    )

# save to quickly test plot without killing R
# sf::write_sf(chi_95_sf, "chi_95_sf.geojson")

scales_chi <- seq(-0.1, 0.6, by = 0.1)
# lab_chi <- TeX(r"($\chi$)")
# lab_chi <- TeX(r"($\chi(0.95)$)")
# lab_chi <- expression(chi(u))
lab_chi <- "chi(0.95)"
chi_p <- ggplot(areas) +
    geom_sf(colour = "black", fill = NA) +
    geom_sf(
    # geom_sf_pattern(
      data = filter(chi_95_sf, var == "chi") %>% 
        mutate(show_chi = factor(ifelse(show_chi == TRUE, "yes", "no"))),
      aes(fill = value, size = value),
      # aes(fill = value, size = value, pattern = show_chi) # ,
      stroke = 1,
      pch = 21
    ) +
    # facet_wrap(~var) + 
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    # viridis::scale_fill_viridis(
    #     option = "F",
    #     breaks = scales_chibar,
    #     labels = as.character(scales_chibar),
    #     guide = "legend"
    # ) + 
    scale_fill_gradientn(
      colours = RColorBrewer::brewer.pal(name = "Blues", n = 7),
      # colours = wesanderson::wes_palette("Zissou1", 8, type = "continuous"), 
      # colours = rev(heat.colors(7)),
      breaks = scales_chi,
      labels = as.character(scales_chi),
      guide = "legend"
    ) + 
    scale_size_continuous(
      range = point_ranges,
      breaks = scales_chi,
      labels = as.character(scales_chi),
      guide = "legend"
    ) +
    # scale_pattern_manual(values = c("stripe", "none")) + 
    # scale_pattern_manual(values = c("crosshatch", "none")) + 
    # maximise plot within frame
    coord_sf(expand = FALSE) + 
    labs(fill = lab_chi, size = lab_chi) + 
    guides(fill = guide_legend(), size = guide_legend(), pattern = "none") +
    theme +
    theme(
      legend.position = "right",
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.key = element_blank()
    )

(chi_plots <- chibar_p + chi_p)

ggsave("latex/plots/041_chi_plots.png", chi_plots, width = 6.3, height = 6, units = "in")


#### Compare to "out-of-box" texmex predictions ####

# locs <- unique(data_gpd$name)
# texmex_marginal <- lapply(locs, \(loc) {
#   print(loc)
#   vals <- data %>% 
#     filter(name == loc) %>% 
#     dplyr::select(rain, wind_speed)
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
#       dplyr::select(name, lon, lat, !!sym(col), !!sym(col_texmex)) %>% 
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
# TODO: Repeat after removing weeks with no rain observations!

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
# bootstrap_thresh_plots <- create_bootstrap_plot(marginal)
# # save object so it can be loaded again later, saves running above!
# saveRDS(bootstrap_thresh_plots, "bootrstrap_thresh_plots.RDS")
# save_bootstrap_plot(
#   bootstrap_thresh_plots,
#   "plots/bootstrap_weekly_winter.pdf"
# )

# reload bootstrap plots
bootstrap_thresh_plots <- readRDS("bootstrap_thresh_plots.RDS")
bootstrap_thresh_plots <- bootstrap_thresh_plots[
  !rm_cnty(names(bootstrap_thresh_plots)) %in% loc_nas # remove problem sites
]

# plot for one area
spec_plot <- bootstrap_thresh_plots$`Dublin (Ringsend) - Dublin`
rain_plots <- (spec_plot[[2]][[1]] + 
  labs(y = "a rain | wind speed", title = "")) / 
(spec_plot[[2]][[2]] + 
  labs(y = "b rain | wind speed", title = ""))

wind_plots <- (spec_plot[[1]][[1]] + 
  labs(y = "a wind speed | rain", title = "")) / 
(spec_plot[[1]][[2]] + 
  labs(y = "b wind speed | rain", title = ""))

bootstrap_thresh_plot <- (rain_plots | wind_plots) + plot_layout(ncol = 2, nrow = 1)
  
ggsave("latex/plots/042_bootstrap_thresh.png", bootstrap_thresh_plot, width = 6.3, height = 6, units = "in")

# for fixed b
# bootstrap_thresh_plots_fixed_b <- create_bootstrap_plot(marginal, fixed_b = TRUE)
# saveRDS(bootstrap_thresh_plots_fixed_b, "bootrstrap_thresh_plots_fixed_b.RDS")
# save_bootstrap_plot(
#   bootstrap_thresh_plots_fixed_b, 
#   "plots/bootstrap_weekly_winter_fixed_b.pdf"
# )

#### Fit CE model ####

# fit CE dependence model for each site
dependence <- fit_texmex_dep(
  marginal, 
  mex_dep_args = list(
    start = c(0.01, 0.01), 
    dqu = 0.7,
    fixed_b = FALSE,
    PlotLikDo = FALSE
  ) 
)

# pull names of locations with NA for alpha 
# col sums should equal 6 if dependence fit is valid for rain and wind
loc_nas <- names(marginal)[
  colSums(data.frame(
    # find length of both objects
    vapply(dependence, \(x) vapply(x, length, numeric(1)), numeric(2))
  )) < 6 # test it is the correct size
]

# remove locations with NAs
dependence <- dependence[!names(dependence) %in% loc_nas]

# separate county name from names
if (length(loc_nas) > 0) {
  loc_nas <- vapply(
    stringr::str_split(loc_nas, " - "), `[[`, 1, FUN.VALUE = character(1)
  )
}
data_gpd <- data_gpd[!data_gpd$name %in% loc_nas, ]

# save object
saveRDS(dependence, file = "data/texmex_mexdep_obj.RDS")
# saveRDS(dependence, file = "data/texmex_mexdep_obj_fixed_b.RDS")

#### Untransformed and transformed observations ####

#### Profile likelihood for one location ####

# Fit with unfixed b for Cloone Lake for report plots
mex_cloone <- mexDependence(
  # x, 
  marginal$`Cloone Lake (Caragh River Area) - Kerry`,
  which = "rain", 
  dqu = 0.7, 
  nOptim = 2, 
  fixed_b = FALSE, 
  PlotLikDo = TRUE, 
  # PlotLikTitle = "Profile Likelihood, rain | wind speed"
  PlotLikTitle = ""
)

# mod <- mexDependence(
#   # x, 
#   marginal[[7]],
#   which = "wind_speed", 
#   dqu = 0.7, # seems a reasonable choice of threshold given bootstrap plots
#   nOptim = 2, 
#   fixed_b = fixed_b, 
#   start = start_vals,
#   PlotLikDo = TRUE, 
#   PlotLikTitle = "wind speed | rain"
# )

#### Bootstrap parameter estimates and calculate quantiles ####

# spec_dep_var <- dependence[[i]][[1]]

fixed_b <- FALSE # doesn't work for fixed b
# final_bootstraps <- lapply(seq_along(dependence), \(i) {
#   print(names(dependence)[[i]])
#   spec_dep <- dependence[[i]]
#   print(i)
#   lapply(spec_dep, \(spec_dep_var) {
#     # boots <- suppressWarnings(suppressMessages(
#     #   bootmex(spec_loc_var)
#     # ))
#     boots <- suppressWarnings(suppressMessages(bootmex(
#       spec_dep_var, 
#       nPass = 3,
#       fixed_b = fixed_b # fix b value
#     )))
#     
#     # Calculate credible interval on bootsrapped samples of a and b
#     # pulled from texmex:::plot.bootmex
#     d2 <- dim(boots$boot[[1]][[2]]) # "dependence"
#     margins <- boots$margins
#     x <- boots$boot
#     which <- 2 # dependence rather than marginal parameters
#     
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
#       quantiles <- quantile(x, c(0.025, 0.975), na.rm = TRUE)
#       val <- prod(quantiles)
#       return(list("quantiles" = quantiles, "prod" = val, "test" = val < 0))
#     })
#     names(tests) <- c("a", "b")
#     return(list("boots" = boots, "values" = vals, "tests" = tests))
#   })
# })
# names(final_bootstraps) <- names(dependence)

# save object, takes agggges to make
# saveRDS(final_bootstraps, "data/bootstraps.RDS")
# final_bootstraps <- readRDS("data/bootstraps.RDS")

final_bootstraps <- readRDS("data/bootstraps.RDS")
final_bootstraps <- final_bootstraps[
  !rm_cnty(names(final_bootstraps)) %in% loc_nas # remove problem sites
]

# plot for single location
bootstrap_cloone <- with(
  final_bootstraps$`Cloone Lake (Caragh River Area) - Kerry`,
  data.frame(
    "a_vals" = rain$values$a,
    "b_vals" = rain$value$b,
    "a"      = mex_cloone$dependence$coefficients[1],
    "b"      = mex_cloone$dependence$coefficients[2]
  )
) %>% 
  mutate(
    a_lower = quantile(a_vals, probs = 0.025),
    a_upper = quantile(a_vals, probs = 0.975),
    b_lower = quantile(b_vals, probs = 0.025),
    b_upper = quantile(b_vals, probs = 0.975)
    # a_radius = quantile(a_vals, probs = 0.975) - a
  )

# bootstrap_cloone %>% 
#   ggplot() + 
#   geom_point(aes(x = a_vals, y = b_vals)) +
#   # geom_errorbar(
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   # geom_pointrange(
#   #   data = slice(bootstrap_cloone, 1),
#   #   aes(x = a, y = b, xmin = a_lower, xmax = a_upper, ymin = b_lower, ymax = b_upper), colour = "red"
#   # )
#   # geom_crossbar(
#   #   data = slice(bootstrap_cloone, 1),
#   #   aes(x = a, y = b, xmin = a_lower, xmax = a_upper, ymin = b_lower, ymax = b_upper)
#   # ) + 
#   geom_errorbar(
#     data = slice(bootstrap_cloone, 1),
#     aes(x = a, y = b, ymin = b_lower, ymax = b_upper),
#     width = 0.2
#   ) + 
#   geom_errorbar(
#     data = slice(bootstrap_cloone, 1),
#     aes(x = a, y = b, xmin = a_lower, xmax = a_upper), 
#     width = 0.1
#   ) + 
#   geom_point(aes(x = a, y = b), colour = "red", size = 5) + 
#   # geom_errorbarh(
#   # geom_crossbar(
#   # # geom_errorbar(
#   #   data = slice(bootstrap_cloone, 1),
#   #   aes(x = a, y = b, xmin = a_lower, xmax = a_upper), colour = "red"
#   # ) + 
#   labs(x = "a", y = "b") + 
#   # geom_point(aes(x = a, y = b), colour = "red", size = 4) + 
#   theme

bootstrap_density_p <- ggplot(bootstrap_cloone, aes(a_vals, b_vals)) +
  geom_density_2d(aes(colour = ..level..), linewidth = 2, alpha = 0.8) + 
  scale_color_viridis_c() + 
  geom_point(size = 2) +
  geom_point(aes(x = a, y = b), colour = "red", size = 4) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "a", y = "b") + 
  guides(colour = "none") + 
  theme

ggsave("latex/plots/044_bootstrap_density.png", bootstrap_density_p, width = 6.3, height = 6, units = "in")

# Number of locations where alpha (and beta) doesn't cross 0
# Could be coded better!!
tests <- lapply(final_bootstraps, \(loc) { # loop through locations
 lapply(loc, \(var) { # loop through variables (rain vs ws)
   list("a" = var$tests$a$test , "b" = var$tests$b$test)
 }) 
})

# how many alphas pass 0 for rain?
sum(vapply(tests, \(x) x$rain$a, logical(1)))
sum(vapply(tests, \(x) x$rain$b, logical(1)))
# alphas for wind speed
sum(vapply(tests, \(x) x$wind_speed$a, logical(1)))
# betas for wind speed
sum(vapply(tests, \(x) x$wind_speed$b, logical(1)))

  
#### Plots ####

# pull a and b values for each site
ab_vals <- lapply(dependence, \(x) {
# ab_vals <- lapply(seq_along(dependence), \(i) {
  list(
    summary(x[[1]])$dependence[1:2],
    # summary(dependence[[i]][[1]])$dependence[1:2],
    summary(x[[2]])$dependence[1:2]
    # summary(dependence[[i]][[2]])$dependence[1:2]
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

# save for clustering
readr::write_csv(ab_df, "data/ab_vals.csv")
# readr::write_csv(ab_df, "data/ab_vals_fixed_b.csv")

# plot a versus other variables 
# TODO: Add wind direction here, could be interesting!
# TODO: Colour points by county? Or coastal vs inland?
# TODO: Could also colour points for wind speed ~1 as red and ignore in lm 
# (or add additional fitted line for including these points)
pred_dat <- ab_df %>% 
  filter(parameter == "a") %>% 
  filter(value < 0.99) %>% 
  dplyr::select(-c(dist2coast, alt, wind_dir)) %>% 
  pivot_longer(
    # c("alt", "dist2coast", "lat", "lon"), 
    c("lat", "lon"), 
    names_to = "predictor", 
    values_to = "pred_val"
  ) %>% 
  mutate(
    var = ifelse(var == "rain", "Precipitation", "Wind speed"),
    predictor = ifelse(predictor == "lat", "Latitude", "Longitude")
  )
  
pred_p <- (pred_dat %>% 
  filter(predictor == "Longitude") %>% 
  ggplot(aes(x = value, y = pred_val)) +
  # geom_point(aes(colour = county)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  # ggpubr::stat_cor(
  #   cor.coef.name = "rho", 
  #   method = "pearson", 
  #   digits = 2, 
  #   p.accuracy = 0.001
  # ) + 
  # facet_wrap(var ~ predictor, scales = "free_y", nrow = 2) + 
  facet_wrap(~ var, scales = "free_y", nrow = 2) + 
  labs(x = "a", y = "Longitude") + 
  theme) |
(pred_dat %>% 
   filter(predictor == "Latitude") %>% 
   ggplot(aes(x = value, y = pred_val)) +
   # geom_point(aes(colour = county)) + 
   geom_point() + 
   geom_smooth(method = "lm") + 
   # ggpubr::stat_cor(
   #   cor.coef.name = "rho", 
   #   method = "pearson", 
   #   # position = "jitter", 
   #   label.x.npc = 0,
   #   label.y.npc = 1,
   #   digits = 2, 
   #   p.accuracy = 0.001
   # ) + 
   # facet_wrap(var ~ predictor, scales = "free_y", nrow = 2) + 
   facet_wrap(~ var, scales = "free_y", nrow = 2) + 
   labs(x = "a", y = "Longitude") + 
   theme)

# TODO: make notes
# ggsave("plots/a_val_corrs_weekly_winter.png", pred_p, width = 16, height = 12, units = "in")
ggsave("latex/plots/047_alpha_vs_lon_lat.png", pred_p, width = 16, height = 12, units = "in")

ab_df_wide <- ab_df %>% 
  filter(parameter == "a") %>% 
  dplyr::select(name, county, var, value) %>% 
  pivot_wider(names_from = var, values_from = value)

# just plot a for rain vs a for wind!
# TODO: Lots can be done here
# Add county centroids to plot (done)
# TODO: Add fitted line with and without windspeed high a values 
p_compare <- ab_df_wide %>% 
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
  # geom_point(
  #   aes(x = rain_mean, y = windspeed_mean, colour = county), 
  #   size = 6, 
  #   pch = 10, 
  #   alpha = 0.9
  # ) + 
  # site points
  # geom_point(aes(colour = county), size = 3, alpha = 0.9) + 
  geom_point(size = 3, alpha = 0.9) + 
  geom_smooth(colour = "blue", method = "lm") + 
  # doesn't make much difference
  # geom_smooth(aes(x = rain_mean, y = windspeed_mean), colour = "red", method = "lm") + 
  ggpubr::stat_cor(cor.coef.name = "rho", method = "pearson") + 
  scale_x_continuous(limits = c(0, 0.65), expand = c(0, 0)) +
  # scale_y_continuous(limits = c(0, 1), expand = c(0, 0.1)) + 
  scale_y_continuous(limits = c(0, 0.8), expand = c(0, 0)) +
  # scale_y_continuous(limits = c(0, 0.6), expand = c(0, 0)) + # for counties
  labs(x = "a (precipitation)", y = "a (wind speed)") + 
  theme

# ggsave("plots/p_compare_winter.png", p_compare, width = 6.3, height = 6, units = "in")
ggsave("latex/plots/046_alpha_rain_vs_ws.png", p_compare, width = 6.3, height = 6, units = "in")

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

circ_p_rain <- circ_plot(ab_w_dir, "rain")
circ_p_ws <- circ_plot(ab_w_dir, "windspeed")

ggsave("plots/wind_rose_ws_winter.png", circ_p_rain[[1]], width = 6.3, height = 6, units = "in")
ggsave("plots/wind_rose_rain_winter.png", circ_p_ws[[1]], width = 6.3, height = 6, units = "in")
ggsave("plots/polar_rain_winter.png", circ_p_rain[[2]], width = 6.3, height = 6, units = "in")
ggsave("plots/polar_ws_winter.png", circ_p_ws[[2]], width = 6.3, height = 6, units = "in")

# join to previous data (not assigned?)
# data_gpd %>% 
#   left_join(
#     ab_df %>% 
#       mutate(
#         var = ifelse(var == "windspeed", "ws", var),
#         col = paste(parameter, var, sep = "_")
#       ) %>% 
#       dplyr::select(name, col, value) %>% 
#       tidyr::pivot_wider(id_cols = name, names_from = col, values_from = value)
#   )

# change to sf
ab_sf <- st_to_sf(ab_df)
# check that b values aren't exceptionally negative (<-1)
# also check that all values are equal to fixed value for b, if applying
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
    # facet_wrap(var ~ parameter, ncol = 2) + 
    facet_wrap(~ var, ncol = 2) + 
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
p_lst[[1]]

# save
# ggsave("plots/a_ce_ire_winter.png", p_lst[[1]], width = 6.3, height = 6, units = "in")
ggsave("latex/plots/045_alpha_map.png", p_lst[[1]], width = 6.3, height = 6, units = "in")
# ggsave("plots/b_ce_ire.png", p_lst[[2]], width = 6.3, height = 6, units = "in")
# ggsave("plots/a_b_ce_ire_weekly_winter.png", p_lst[[1]] / p_lst[[2]], width = 16, height = 12, units = "in")

#### Diagnostics for CE model ####

# easier to look at in a file

# generate plots
# Add county name to location, helps interpretation (done)
# TODO: Investigate where a is NA in dependence!
plots <- lapply(seq_along(dependence), \(i) { # each location
  lapply(seq_along(dependence[[i]]), \(j) { # rainfall and wind speed
    # print(paste0(i, "-", j))
    if (is.na(dependence[[i]][[j]]$dependence$coefficients[1])) {
      return(NA)
    }
    
    if (j == 1) {
      lab_x1 <- "F(precipitation)"
      lab_y1 <- "Z wind speed | precipitation"
      lab_x2 <- "Precipitation"
      lab_y2 <- "Wind speed"
    } else {
      lab_x1 <- "F(wind speed)"
      lab_y1 <- "Z precipitation | wind speed"
      lab_x2 <- "Wind speed"
      lab_y2 <- "Precipitation" 
    }
    
    # generate plots, adding title to first one
    p <- ggplot(dependence[[i]][[j]], plot. = FALSE)
    p[[1]] <- p[[1]] + 
      # ggtitle(
      #   paste(names(dependence)[i], names(dependence[[i]])[j], sep = " - ")
      # )
      # labs(x = "F(precipitation)", y = "Z wind speed | precipatation") + 
      labs(x = lab_x1, y = lab_y1) + 
      theme
    # remove second plot, doesn't add much!
    p[[2]] <- p[[3]] + 
      # labs(x = "Precipitation", y = "Wind speed") + 
      labs(x = lab_x2, y = lab_y2) + 
      theme
    p[[3]] <- NULL
    return(p)
  })
})
names(plots) <- names(dependence)

# i <- which(grepl("Cloone Lake", names(plots)))
i <- which(grepl("Cloone Lake", names(plots)))

p_diag <- (plots[[i]][[1]][[1]] / plots[[i]][[1]][[2]]) | 
  (plots[[i]][[2]][[1]] / plots[[i]][[2]][[2]]) 

ggsave("latex/plots/048_diag.png", p_diag, width = 6.3, height = 6, units = "in")
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
pdf("plots/diag_plots_rain_weekly_winter.pdf", onefile = TRUE) 
for (i in seq(length(plots))) {
  if (length(plots[[i]][[1]]) == 1 && is.na(plots[[i]][[1]])) next
  do.call("grid.arrange", plots[[i]][[1]])
}
dev.off()

pdf("plots/diag_plots_windspeed_weekly_winter.pdf", onefile = TRUE) 
for (i in seq(length(plots))) {
  # print(i)
  if (length(plots[[i]][[2]]) == 1 && is.na(plots[[i]][[2]])) next
  do.call("grid.arrange", plots[[i]][[2]])
}
dev.off()

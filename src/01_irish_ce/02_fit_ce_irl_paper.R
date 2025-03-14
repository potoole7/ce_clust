#### Conditional extremes modelling of Irish meterological data ####

# Revising original script for use with Irish data
# Rehashing script from TFR for paper on JS clustering

# TODO Remove sites where dependence model doesn't fit
# TODO Do marginal model fitting, checking and plotting
# TODO Do dependence model fitting, checking and plotting
# TODO Save all plots!
# TODO Add useful functions to package

#### Libs ####

# TODO Need for all of these packages?
library(sf)
# library(evc)
devtools::load_all("../evc")
library(dplyr, quietly = TRUE)
library(tidyr)
library(ggplot2)
library(lubridate)
# library(texmex)
devtools::load_all("texmex") # forked texmex which allows for fixed b vals
library(gridExtra)
library(patchwork)
library(geosphere)
library(latex2exp)
library(ggpattern)
library(parallel)

source("src/functions.R")

theme <- evc::evc_theme()

sf::sf_use_s2(FALSE)

# function to remove " - " (county name) from names
rm_cnty <- \(x) {
  vapply(stringr::str_split(x, " - "), `[[`, 1, FUN.VALUE = character(1))
}


#### Metadata ####

# dqu <- 0.7 # dependence threshold
dqu <- 0.9 # dependence threshold
min_max_dates <- as_date(c("1990-01-01", "2020-12-31"))
all_dates <- seq.Date(min_max_dates[1], min_max_dates[2], by = "day")
fixed_xi <- FALSE # whether to fix shape parameter in GPD margins

ncores <- detectCores() - 1


#### Load Data ####

data <- readr::read_csv(
  "data/met_eireann/final/met_eir_wind.csv.gz",
  # "data/met_eireann/final/met_eir_wind_alt.csv.gz",
  show_col_types = FALSE
) %>%
  filter(date %in% all_dates, !is.na(wind_speed))

# load shapefile
areas <- sf::read_sf("data/met_eireann/final/irl_shapefile.geojson")

# remove problematic locations
data <- data |>
  filter(!name %in% c(
    # remove Castlereagh (Down), outlier in QQ plot for wind speed GPD fit
    "Castlereagh",
    # remove sites where model fails to fit for Keef constraints
    "Crolly (Filter Works)", "Kilcar (Cronasillagh)", "Wexford (Newtown W.w.)"
  ))

# pull just site names, counties and provinces
county_key_df <- data |>
  distinct(name, county, province)


#### Preprocess Data ####

# function to find value corresponding to maximum density estimate for vector
max_density <- \(vec) {
  if (length(vec) == 1) {
    return(vec)
  }
  # return maximum kernel density estimate
  return(with(density(vec), x[which.max(y)]))
}

# take Winter data only
data_winter <- data %>%
  mutate(month = as.numeric(substr(date, 6, 7))) %>%
  filter(month %in% c(1:3, 10:12)) %>%
  dplyr::select(-month)

# Take weekly sum of precipitation and average of wind speed, to account for
# lag in storm (as in Vignotto, Engelke 2021 study of GB + Ireland)
data_week <- data_winter %>%
  mutate(week = floor_date(date, "week")) %>%
  group_by(week, name, county, province, lon, lat) %>%
  # group_by(week, name, county, province, lon, lat, alt) %>%
  summarise(
    rain       = sum(rain, na.rm = TRUE),
    wind_speed = mean(wind_speed, na.rm = TRUE),
    wind_dir   = max_density(wind_dir),
    .groups    = "drop"
  ) %>%
  rename(date = week) |>
  # remove weeks with no rainfall, as in Vignotto 2021 (required/important?)
  filter(rain != 0)

# investigate how many weeks might be missing for each location
locs <- data_week |>
  arrange(desc(province), county) |>
  distinct(name) |>
  pull()
plots <- lapply(locs, \(x) {
  data_week |>
    filter(name == x) |>
    pivot_longer(c(rain, wind_speed), names_to = "var") |>
    ggplot(aes(x = date, y = value)) +
    geom_point() +
    facet_wrap(~var, scales = "free_y") +
    # scale_x_date(date_breaks = "6 months", date_labels = "%Y-%m") +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    ggtitle(paste(
      x,
      data_week$county[data_week$name == x][1],
      data_week$province[data_week$name == x][1],
      sep = " - "
    )) +
    theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
})

# TODO Why are the years 2000 up to end of 2002 missing?!?!
# Investigate wind speed data for this! No download6.rc!!
pdf("plots/irl_explore_missing.pdf", width = 8, height = 8)
plots
dev.off()

data_dates <- unique(data_week$date)
all_dates <- seq.Date(min(data_dates), max(data_dates), by = "week")
(missing_dates <- as_date(setdiff(all_dates, data_dates)))


#### Data exploration ####

# locations with highest and lowest rain, to specifically plot
highest_lowest_rain <- data_week %>%
  group_by(name) %>%
  summarise(rain = mean(rain, na.rm = TRUE), .groups = "drop") %>%
  arrange(rain) %>%
  slice(c(1, n())) %>%
  pull(name)

# plot wind speeds against rain for sites with the highest and lowest rainfall
p_wind_rain <- data_week %>%
  filter(name %in% highest_lowest_rain, rain > 0) %>%
  # add 95% quantile lines for both
  group_by(name) %>%
  mutate(across(
    c(rain, wind_speed), ~ quantile(.x, 0.95, na.rm = TRUE),
    .names = "quant_{.col}"
  )) %>%
  ungroup() %>%
  ggplot(aes(x = rain, y = wind_speed)) +
  geom_point(aes(colour = name), size = 1.5, alpha = 0.8) +
  geom_vline(aes(xintercept = quant_rain), linetype = "dashed") +
  geom_hline(aes(yintercept = quant_wind_speed), linetype = "dashed") +
  facet_wrap(~name, scales = "free_x") +
  scale_colour_manual(values = c(ggsci::pal_nejm()(2))) +
  labs(
    x      = "precipitation (mm)",
    y      = "wind speed (m/s)",
    colour = ""
  ) +
  theme +
  # remove facet labels, colour will do
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.key = element_blank()
  )


#### Chi exploration ####

# Plot chi and chi-squared for one and all locations

# calculate chi for one location
chi_spec <- data_week %>%
  filter(name == highest_lowest_rain[1]) %>%
  dplyr::select(rain, wind_speed) %>%
  texmex::chi()

(chi_plot_spec <- chi_spec %>%
  # use texmex plotting method for chi and chibar
  ggplot(main = c("ChiBar" = "", "Chi" = ""), plot. = FALSE) |>
  lapply(\(x) x + theme) |>
  wrap_plots() +
  # add centred title through patchwork
  patchwork::plot_annotation(
    title = paste0("Tail Dependence for ", highest_lowest_rain[1]),
    theme = theme
  ))

# for each location, pull out 95th quantile for chibar
# Colour chi grey if not valid
# Plot on map
# TODO Functionalise for US air pollution!
names <- unique(data_week$name)
chi_95_df <- bind_rows(mclapply(names, \(x) {
  chi <- data_week %>%
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
}, mc.cores = ncores))
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

# function to plot chi and chi bar on map
scales <- seq(-0.1, 0.6, by = 0.1)
chi_map_plot <- \(
  chi_95_sf,
  var = c("chi", "chibar"),
  scales = seq(-0.1, 0.6, by = 0.1),
  point_ranges = c(2, 6)
) {
  lab <- ifelse(var == "chi", expression(chi(u)), expression(bar(chi)(u)))
  plot_data <- filter(chi_95_sf, var == !!var)
  if (var == "chi") {
    plot_data <- plot_data %>%
      mutate(show_chi = factor(ifelse(show_chi == TRUE, "yes", "no")))
  }

  p <- ggplot(areas) +
    geom_sf(colour = "black", fill = NA) +
    geom_sf(
      # geom_sf_pattern(
      data = plot_data,
      aes(fill = value, size = value),
      # aes(fill = value, size = value, pattern = show_chi),
      pch = 21,
      stroke = 1
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    # add colour scheme afterwards
    # scale_fill_gradientn(
    #   colours = RColorBrewer::brewer.pal(name = "Blues", n = 7),
    #   breaks = scales_chi,
    #   labels = as.character(scales_chi),
    #   guide = "legend"
    # ) +
    scale_size_continuous(
      range  = point_ranges,
      breaks = scales,
      labels = as.character(scales),
      guide  = "legend"
    ) +
    # scale_pattern_manual(values = c("stripe", "none")) +
    # maximise plot within frame
    coord_sf(expand = FALSE) +
    labs(fill = lab, size = lab) +
    guides(fill = guide_legend(), size = guide_legend(), pattern = "none") +
    theme +
    theme(
      legend.position = "right",
      axis.text       = element_blank(),
      axis.ticks      = element_blank(),
      legend.key      = element_blank()
    )

  return(p)
}

# plot chibar and chi
chibar_p <- chi_map_plot(chi_95_sf, "chibar") +
  scale_fill_gradientn(
    colours = rev(heat.colors(7)),
    breaks = scales,
    labels = as.character(scales),
    guide = "legend"
  )
chi_p <- chi_map_plot(chi_95_sf, "chi") +
  scale_fill_gradientn(
    colours = RColorBrewer::brewer.pal(name = "Blues", n = 7),
    breaks = scales,
    labels = as.character(scales),
    guide = "legend"
  )

# combine plots
(chi_plots <- chibar_p + chi_p)

# ggsave("latex/plots/041_chi_plots.png", chi_plots, width = 6.3, height = 6, units = "in")


#### Thresholding ####

# check max and min obs for each site, and approx. how many after thresholding
data_week |>
  group_by(name, county) |>
  summarise(n = n(), .groups = "drop") |>
  filter(n %in% c(min(n), max(n))) |>
  mutate(n_95 = floor(n - 0.95 * n), n_90 = floor(n - 0.90 * n))
# max obs = 37 @ 95% quantile, 74 @ 90% quantile
# min obs = 9 @ 95% quantile, 18 @ 90% quantile, maybe go with 90?


# explore difference between thresholding via evgam and ald and qgam
f_evgam <- list(excess ~ s(lon, lat), ~ s(lon, lat))
if (fixed_xi) {
  f_evgam[[2]] <- ~1
}
shared_args <- list(
  data = data_week,
  vars = c("rain", "wind_speed"),
  # arguments to `quantile_thresh`
  marg_prob = list(
    f = list("response ~ s(lon, lat)", "~ s(lon, lat)"),
    # tau = 0.9
    tau = 0.95
  ),
  # formula for evgam fit to excesses over threshold
  f = f_evgam,
  ncores = max(1, parallel::detectCores() - 2),
  marg_only = TRUE
)
mod_fit_ald <- do.call(
  fit_ce, c(shared_args, list(thresh_fun = evgam_ald_thresh))
)
mod_fit_qgam <- do.call(
  fit_ce, c(shared_args, list(thresh_fun = qgam_thresh))
)

# plot thresholds; decide which seems more reasonable

# first, pull thresholded data together
prep_thresh_data <- \(data_rain, data_ws) {
  data_rain %>%
    distinct(name, lon, lat, thresh_rain = thresh) %>%
    left_join(
      data_ws %>%
        distinct(name, lon, lat, thresh_ws = thresh)
    ) %>%
    pivot_longer(contains("thresh"), names_to = "var")
}
data_thresh_ald <- st_to_sf(prep_thresh_data(
  mod_fit_ald$data_thresh$rain, mod_fit_ald$data_thresh$wind_speed
))
data_thresh_qgam <- st_to_sf(prep_thresh_data(
  mod_fit_qgam$data_thresh$rain, mod_fit_qgam$data_thresh$wind_speed
))

# plot barplots of difference in exceedances
data_thresh_ald |>
  rename(thresh_ald = value) |>
  left_join(
    data_thresh_qgam |>
      rename(thresh_qgam = value) |>
      st_drop_geometry()
  ) |>
  left_join(data_week |> distinct(name, province) |> st_drop_geometry()) |>
  # order plot by province
  arrange(province) |>
  mutate(name = factor(name, levels = unique(name))) |>
  ggplot(aes(x = name, y = thresh_ald - thresh_qgam, fill = province)) +
  geom_bar(stat = "identity") +
  facet_wrap(~var, scales = "free_y") +
  theme +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # remove x-axis labels, too bunched
  theme(axis.text.x = element_blank())
# qgam threshold consistently higher than ald threshold, but not by much


# Plot thresholds on map for both

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

join_plot_thresh <- \(data_thresh, title) {
  p_thresh_rain <- plot_thresh(
    filter(data_thresh, var == "thresh_rain"), areas, seq(50, 150, by = 25)
  ) +
    labs(fill = "u for precipitation", size = "u for precipitation")

  p_thresh_ws <- plot_thresh(
    filter(data_thresh, var == "thresh_ws"), areas, seq(8, 11, by = 0.5)
  ) +
    labs(fill = "u for wind speed", size = "u for wind speed")

  p_thresh_rain + p_thresh_ws +
    patchwork::plot_annotation(
      title = title,
      theme = theme(plot.title = element_text(hjust = 0.5))
    )
}

join_plot_thresh(data_thresh_ald, "ALD")
join_plot_thresh(data_thresh_qgam, "qgam")
# ggsave("latex/plots/031_thresh_plots.png", p_thresh, width = 6.3, height = 6, units = "in")

# plot exceedances for each
plot_exceedances <- \(mod_fit) {
  data_thresh_plt <- data_week %>%
    distinct(name, lon, lat) %>%
    left_join(
      count(mod_fit$data_thresh$rain, name) %>%
        rename(n_rain = n)
    ) %>%
    left_join(
      count(mod_fit$data_thresh$wind_speed, name) %>%
        rename(n_wind = n)
    ) %>%
    pivot_longer(n_rain:n_wind, names_to = "var", values_to = "n") %>%
    st_to_sf()

  ggplot(areas) +
    geom_sf(colour = "black", fill = NA) +
    geom_sf(
      data = data_thresh_plt,
      aes(fill = n),
      size = 6,
      stroke = 1,
      pch = 21
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    facet_wrap(~var) +
    theme +
    theme(legend.position = "right")
}

plot_exceedances(mod_fit_ald)
plot_exceedances(mod_fit_qgam)

# conclusion: Go with qgam
mod_fit <- mod_fit_qgam

# Extract thresholded data and parameter estimates for plotting
data_rain <- mod_fit$data_thresh$rain
data_ws <- mod_fit$data_thresh$wind_speed

data_rain_map <- data_rain %>%
  left_join(mod_fit$evgam_fit$rain$predictions) %>%
  st_to_sf(remove = FALSE) # add location
data_ws_map <- data_ws %>%
  left_join(mod_fit$evgam_fit$wind_speed$predictions) %>%
  st_to_sf(remove = FALSE) # add location

#### Marginal Checking ####

# Probability and Quantile-Quantile plots for checking model fit
# code by Kristina Bratkova!
# TODO Add to package, useful function!
pp_qq_plot <- function(data_map, xlab = "Empirical", debug_q = FALSE) {
  n <- nrow(data_map)
  # empirical probabilities and quantiles
  probs <- seq_len(n) / (n + 1)
  quantiles <- -log(1 - probs)

  # Transform to exponential margins
  # TODO Only requires mod_fit, can get excess from thresholds
  y_tilde <- (1 / data_map$shape) *
    log(1 + data_map$shape * (data_map$excess / data_map$scale))

  # sorted model probabilities and quantiles
  model_probs <- 1 - exp(-y_tilde)
  ord_p <- order(model_probs)
  model_probs <- model_probs[ord_p]
  ord_q <- order(y_tilde)
  # model_quantiles <- sort(y_tilde)
  model_quantiles <- y_tilde[ord_q]

  # Compute confidence bounds using qbeta method
  Ulow <- sapply(seq_len(n), function(i) qbeta(0.025, i, n + 1 - i))
  Uup <- sapply(seq_len(n), function(i) qbeta(0.975, i, n + 1 - i))

  # dataframes for plotting
  prob_plot_df <- data.frame(
    name  = data_map$name[ord_p],
    x     = probs,
    y     = model_probs,
    lower = Ulow,
    upper = Uup
  )
  qq_plot_df <- data.frame(
    name  = data_map$name[ord],
    x     = quantiles,
    y     = model_quantiles,
    lower = -log(1 - Ulow),
    upper = -log(1 - Uup)
  )

  # Probability Plot with confidence intervals
  pp <- ggplot(prob_plot_df, aes(x = x, y = y)) +
    geom_abline(intercept = 0, slope = 1, colour = "#C11432", size = 1.5) +
    geom_point() +
    geom_line(aes(y = upper), linetype = "dashed", colour = "#C11432") +
    geom_line(aes(y = lower), linetype = "dashed", colour = "#C11432") +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#C11432", alpha = 0.2) +
    theme_minimal() +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = xlab, y = "Model")

  # QQ Plot with confidence intervals

  max_x <- max(qq_plot_df$x)
  max_y <- max(c(qq_plot_df$y, qq_plot_df$upper))

  # TEMP: return names of locations with Empirical > 6
  if (debug_q == TRUE) {
    qq_plot_df |>
      filter(x > 6) |>
      pull(name)
  }

  qq_p <- ggplot(qq_plot_df, aes(x = x, y = y)) +
    geom_segment(
      aes(x = 0, y = 0, xend = max(x), yend = max(x)),
      colour = "#C11432",
      size = 1.5
    ) +
    geom_point() +
    geom_line(aes(y = upper), linetype = "dashed", colour = "#C11432") +
    geom_line(aes(y = lower), linetype = "dashed", colour = "#C11432") +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#C11432", alpha = 0.2) +
    theme_minimal() +
    scale_x_continuous(
      limits = c(0, max_x),
      expand = c(0, 0.1)
    ) +
    scale_y_continuous(
      limits = c(0, max_y),
      expand = c(0, 0.1)
    ) +
    labs(x = xlab, y = "")

  return(list(pp, qq_p))
}

# plot for wind speed and rain
pp_qq_rain <- (
  pp_qq_plot(data_rain_map) |>
    wrap_plots()
) +
  patchwork::plot_annotation(title = "Rain", theme = theme)

# TODO Label points with weird values! Look at these locations individually
pp_qq_ws <- (
  pp_qq_plot(data_ws_map) |>
    wrap_plots()
) +
  patchwork::plot_annotation(title = "Wind Speed", theme = theme)

plot_lab <- ifelse(fixed_xi, "fixed_xi", "vary_xi")
ggsave(
  paste0("plots/ire/pp_qq_rain_", plot_lab, ".png"),
  pp_qq_rain,
  width = 6.3,
  height = 6,
  units = "in"
)
ggsave(
  paste0("plots/ire/pp_qq_ws_", plot_lab, ".png"),
  pp_qq_ws,
  width = 6.3,
  height = 6,
  units = "in"
)
# conclusion: allowing shape parameter to vary seems to lead to better fit!



#### Marginal plotting ####

# TODO Plot shape and scale parameters for each location
# TODO: Add to evc?
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
      viridis::scale_fill_viridis(
        breaks = scales,
        labels = as.character(scales),
        guide = "legend"
      )
  }
  return(p)
}

# for rain
lab_scale_rain <- "scale"
p_scale_rain <- plot_param(data_rain_map, areas, seq(10, 45, by = 5)) +
  labs(fill = lab_scale_rain, size = lab_scale_rain)

lab_shape_rain <- "shape"
p_shape_rain <- plot_param(
  data_rain_map, areas, round(seq(-0.15, 0.15, by = 0.05), 2), "shape"
) +
  labs(fill = lab_shape_rain, size = lab_shape_rain)

ggsave(
  "plots/ire/scale_rain_ire.png", p_scale_rain,
  width = 6.3, height = 6, units = "in"
)
ggsave(
  "plots/ire/shape_rain_ire.png", p_shape_rain,
  width = 6.3, height = 6, units = "in"
)

# for wind speed
lab_scale_ws <- "scale"
p_scale_ws <- plot_param(data_ws_map, areas, seq(0.5, 0.8, by = 0.05)) +
  labs(fill = lab_scale_ws, size = lab_scale_ws) +
  NULL

lab_shape_ws <- "shape"
p_shape_ws <- plot_param(data_ws_map, areas, seq(-0.3, -0.1, by = 0.05), "shape") +
  labs(fill = lab_shape_ws, size = lab_shape_ws) +
  scale_fill_gradient2(
    low = "red3",
    high = "blue3",
    na.value = "grey",
    breaks = seq(-0.3, -0.1, by = 0.05),
    labels = as.character(seq(-0.3, -0.1, by = 0.05)),
    guide = "legend"
  )

ggsave(
  "plots/ire/scale_ws_ire.png", p_scale_ws,
  width = 6.3, height = 6, units = "in"
)
ggsave(
  "plots/ire/shape_ws_ire.png", p_shape_ws,
  width = 6.3, height = 6, units = "in"
)

# Plot endpoints where shape parameter is negative
plot_endpoint <- \(data_map, scales, range) {
  data_map %>%
    mutate(end_point = ifelse(shape < 0, thresh - (scale / shape), Inf)) %>%
    plot_param(areas, scales, "end_point") +
    scale_size_continuous(
      range = range,
      breaks = scales,
      labels = as.character(scales),
      guide = "legend"
    ) +
    labs(fill = "endpoint (m/s)", size = "endpoint (m/s)")
}

p_end_rain <- plot_endpoint(data_rain_map, c(500, 1000, 2000, 7000), c(1.5, 10))
p_end_ws <- plot_endpoint(data_ws_map, seq(10, 25, by = 1), c(1.5, 10))

ggsave(
  "plots/ire/endpoint_rain.png", p_end_rain,
  width = 6.3, height = 6, units = "in"
)
ggsave(
  "plots/ire/endpoint_ws.png", p_end_ws,
  width = 6.3, height = 6, units = "in"
)

#### Dependence Modelling ####

dep_args <- shared_args
# add marg fits as data (will skip marginal fitting)
dep_args$data <- mod_fit
dep_args$marg_only <- FALSE
dep_args$ncores <- 1 # for debugging

dep_fit <- do.call(
  fit_ce,
  c(dep_args, list(
    cond_prob   = 0.9,
    fit_no_keef = TRUE,
    # fit_no_keef = FALSE, # causes a lot of models to fail!
    output_all  = TRUE
  ))
)
# fails for Crolly (Filter Works), Kilcar (Cronasillagh) (both Donegal) &
# Wexford (Newtown W.w.)

# TODO Find and remove failed model fits
fail_locs <- sapply(dep_fit$dependence, \(x) any(is.na(unlist(x))))

# save
saveRDS(dep_fit, "data/ire_dep_fit.rds")
# saveRDS(dep_fit, "data/ire_dep_fit.rds")
dep_fit <- readRDS("data/ire_dep_fit.rds")

#### Check residuals ####

# TODO Extend to work with > 2 variables
plot_resid <- \(spec_resid, spec_trans) {
  vars <- names(spec_resid)
  # TODO do for one var, then loop
  # var <- vars[[1]]
  plots <- lapply(vars, \(var) {
    Z <- spec_resid[[var]] # TODO Will also have to loop through these variables!
    if (all(is.na(Z))) {
      return(NA)
    }
    yex <- spec_trans[, var]

    n <- length(Z)
    p <- seq(dqu, 1 - (1 / n), length = n)

    data.frame(p, "resid" = Z) |>
      ggplot(aes(x = p, y = Z)) +
      geom_point(alpha = 0.7) +
      geom_smooth() +
      theme +
      labs(
        x = paste0("F(", var, ")"),
        y = paste0("Z ", colnames(Z), " | ", var)
      )
  })
  names(plots) <- vars
  return(plots)
}

# df to join site names with counties (easier to know where they are)
residuals <- dep_fit$residual
names_df <- data.frame("name" = names(residuals)) |>
  left_join(county_key_df)

resid_plots <- lapply(seq_along(residuals), \(i) {
  plots <- plot_resid(dep_fit$residual[[i]], dep_fit$transformed[[i]])
  lapply(plots, \(x){
    if (all(is.na(x))) {
      return(NA)
    }
    x # +
    # ggtitle(paste(
    #   names_df$name[i],
    #   names_df$county[i],
    #   sep = " - "
    # ))
  })
})
names(resid_plots) <- names_df$name

# pdf("plots/ire/residuals.pdf", width = 8, height = 8)
# resid_plots
# dev.off()

# Plot quantiles of conditional extremes model
# TODO Tidy this code up a lot!!
# cond_var <- "rain"
quantiles <- seq(0.1, by = 0.2, len = 5)
plot_ce_quantiles <- \(
  dep_fit, spec_loc, cond_var = "rain", quantiles = seq(0.1, by = 0.2, len = 5)
) {
  # TODO Functionalise this somehow??
  # take out data for one location
  # TODO Only keep data required for plots
  dep_fit_spec <- list(
    "marginal" = dep_fit$marginal[[spec_loc]],
    "data_thresh" = lapply(dep_fit$data_thresh, \(x) {
      filter(x, name == spec_loc)
    }),
    "original" = dep_fit$original[[spec_loc]],
    "residual" = dep_fit$residual[[spec_loc]],
    "dependence" = dep_fit$dependence[[spec_loc]],
    "transformed" = dep_fit$transformed[[spec_loc]]
  )

  # non-conditioning/other variables
  vars <- names(dep_fit_spec$data_thresh)
  vars <- vars[vars != cond_var]
  n <- nrow(dep_fit_spec$residual[[cond_var]])

  # marginal and dependence parameters for conditioning variables
  # TODO Remove anything unneeded
  marg <- dep_fit_spec$marginal[[cond_var]]
  dep <- dep_fit_spec$dependence[[cond_var]]
  # thresholds (mth, dth) and quantiles (mqu, dqu) for marginal and dependence
  mqu <- dep_fit$arg_vals$marg_prob
  mth <- dep_fit_spec$marginal[[cond_var]]$thresh
  dqu <- dep_fit$arg_vals$cond_prob
  # dth <- dep_fit_spec$dependence[[cond_var]]["dth", ] # Laplace dth!
  dth <- quantile(dep_fit_spec$original[[cond_var]], dqu)
  # Determine x-axis values to estimate CE quantiles at along conditioned var
  xmax <- max(dep_fit_spec$original[[cond_var]])
  # for negative shape GPD will have upper limit
  if (marg$xi < 0) {
    upper <- marg$thresh - (marg$sigma / marg$xi)
  }
  dif <- xmax - dth
  xlim <- c(dth - 0.1 * dif, dth + 1.5 * dif)

  # Determine upper limit of x-axis
  # if (xi < 0 && xlim[2] > upperEnd) {
  #   xlim[2] <- upperEnd
  #   plim <- 1
  # } else {
  #   plim <- pgpd(xlim[2], sigma = sig, xi = xi, u = marThr)
  # }
  if (marg$xi < 0 && xlim[2] > upper) {
    xlim[2] <- upper
    plim <- 1
  } else {
    plim <- evd::pgpd(xlim[[2]], mth, marg$sigma, marg$xi)
  }
  # CDF probabilities to plot at
  p <- seq(dqu, 1 - 1 / n, length = n)
  # take out largest point to avoid Inf in CDF transform
  len <- 501
  plotp <- seq(dqu, plim, len = len)[-len]
  # transform to original scale; these will be x-values in plot
  # plotx <- revTransform(
  #   plotp,
  #   qu = dep$dqu,
  #   data = mar$data[, dep$which],
  #   th = depThr,
  #   sigma = sig,
  #   xi = xi
  # )
  # TODO Gives different results to revTransform, why????
  # plotx <- inv_semi_par_cdf(
  #   F_hat = matrix(plotp),
  #   dat   = dep_fit_spec$original[cond_var],
  #   gpd   = dep_fit_spec$marginal[cond_var] # TODO Need conditional thresh, not marginal thresh!
  # )
  plotx <- revTransform(
    plotp,
    data  = dep_fit_spec$original[[cond_var]],
    qu    = dqu,
    th    = dth,
    sigma = marg$sigma,
    xi    = marg$xi
  )

  # convert probs to Laplace scale (these are values to calculate CE line at)
  # xq <- dep$margins$p2q(plotp)
  xq <- laplace_trans(matrix(plotp))

  # plot on original margins
  # TODO Loop over conditioning variables
  # TODO Extend to > 2 variables
  base_plot <- dep_fit_spec$original |>
    select(all_of(c(cond_var, vars))) |>
    setNames(c("x", "y")) |>
    ggplot(aes(x, y)) +
    geom_point() +
    theme +
    # add vertical line at threshold
    # TODO Loop over variables (i.e. x vs y and y vs x)
    # geom_vline(xintercept = dep_fit_spec$marginal$rain$thresh) +
    geom_vline(xintercept = dth) +
    labs(x = cond_var, y = vars)


  # pull dependence coefficients and quantiles of residuals
  # co <- coef(dep)[, i]
  co <- dep_fit_spec$dependence[[vars]]
  # zq <- quantile(dep$Z[, i], quantiles)
  # zq <- quantile(dep_fit_spec$residual[[cond_var]], quantiles)
  zq <- quantile(dep_fit_spec$residual[[vars]], quantiles)

  # calculates regression lines from quantiles of residuals
  # yq <- sapply(zq, function(z, co, xq) {
  #   co["a"] * xq + co["c"] - co["d"] * log(xq) + xq^co["b"] * z
  # }, xq, co = co)
  yq <- sapply(zq, \(z, xq) {
    # dep$co$a * xq + dep$co$c - dep$co$d * log(xq) + xq^dep$co$b * z
    (co["a", ] * xq) + ((xq^co["b", ]) * z)
  }, xq = xq)

  # transform to original scale
  # ploty <- apply(dep$margins$q2p(yq), 2, revTransform,
  #     data = trns,
  #     qu = qu, th = th, sigma = sigma, xi = xi
  #   )
  # ploty <- inv_semi_par_cdf(
  #   F_hat = yq,
  #   dat   = dep_fit_spec$original[vars],
  #   gpd   = dep_fit_spec$marginal[vars]
  # )
  dth_spec <- quantile(dep_fit_spec$original[[vars]], dqu)
  # ploty <- apply(yq, 2, \(y) {
  ploty <- apply(inv_laplace_trans(yq), 2, \(y) {
    # print(y[1])
    # TODO again, not working, why???
    # inv_semi_par_cdf(
    #   F_hat = matrix(y),
    #   dat   = dep_fit_spec$original[vars],
    #   gpd   = dep_fit_spec$marginal[vars]
    # )
    revTransform(
      y,
      data  = dep_fit_spec$original[[vars]],
      qu    = dqu,
      th    = dth_spec,
      sigma = dep_fit_spec$marginal[[vars]]$sigma,
      xi    = dep_fit_spec$marginal[[vars]]$xi
    )
  })

  # plot CE quantiles recursively
  add_line <- \(p, ploty) {
    if (length(ploty) == 0) {
      return(p)
    } else {
      add_line(
        p +
          geom_line(
            data     = data.frame(x = plotx, y = ploty[, 1]),
            mapping  = aes(x = plotx, y = ploty[, 1]),
            linetype = 2,
            col      = "blue"
          ),
        ploty[, -1, drop = FALSE]
      )
    }
  }
  p <- add_line(base_plot, ploty)

  return(p)
}

quantile_plots <- lapply(names(dep_fit$residual), \(x) {
  print(x)
  list(
    "rain"       = plot_ce_quantiles(dep_fit, x, "rain"),
    "wind_speed" = plot_ce_quantiles(dep_fit, x, "wind_speed")
  )
})

# join resid and quantile plots
diagnostic_plots <- lapply(seq_along(resid_plots), \(i) {
  wrap_plots(resid_plots[[i]]) /
    (quantile_plots[[i]]$rain + quantile_plots[[i]]$wind_speed) +
    # ggtitle(paste(names_df$name[i], names_df$county[i], sep = " - "))
    patchwork::plot_annotation(
      title = paste(names_df$name[i], names_df$county[i], sep = " - "),
      theme = theme
    )
})
names(diagnostic_plots) <- names_df$name

# save
pdf("plots/ire/dep_diagnostic_plots.pdf", width = 8, height = 8)
diagnostic_plots
dev.off()

#### Plot dependence parameters ####

ab_vals <- lapply(dep_fit$dependence, \(x) {
  # ab_vals <- lapply(seq_along(dependence), \(i) {
  list(
    # summary(x[[1]])$dependence[1:2],
    x[[1]][1:2],
    # x[[1]][c("a", "b"), ],
    # summary(x[[2]])$dependence[1:2]
    x[[2]][1:2]
    # x[[2]][c("a", "b"), ]
  )
})
# names(ab_vals) <- data_gpd$name

# convert to dataframe
ab_df <- bind_rows(lapply(seq_along(ab_vals), \(i) {
  data.frame(
    "name" = names(ab_vals)[i],
    "var" = c("rain", "rain", "windspeed", "windspeed"),
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
  left_join(
    data %>%
      # group_by(name, county, dist2coast, alt, lat, lon) %>%
      # group_by(name, county, alt, lat, lon) %>%
      group_by(name, county, lat, lon) %>%
      summarise(wind_dir = max_density(wind_dir), .groups = "drop")
  )

# save for clustering
readr::write_csv(ab_df, "data/ab_vals.csv")
# readr::write_csv(ab_df, "data/ab_vals_fixed_b.csv")

p_compare <- ab_df |>
  pivot_wider(names_from = var, values_from = value) |>
  ggplot(aes(x = rain, y = windspeed)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_smooth(colour = "blue", method = "lm") +
  ggpubr::stat_cor(cor.coef.name = "rho", method = "pearson") +
  facet_wrap(~parameter) +
  theme

# ab_df_wide <- ab_df %>%
#   filter(parameter == "a") %>%
#   dplyr::select(name, county, var, value) %>%
#   # dplyr::select(name, county, var, parameter, value) %>%
#   pivot_wider(names_from = var, values_from = value)
#   # pivot_wider(names_from = c(var, parameter), values_from = value)
#
# # just plot a for rain vs a for wind!
# # TODO: Lots can be done here
# # Add county centroids to plot (done)
# # TODO: Add fitted line with and without windspeed high a values
# p_compare <- ab_df_wide %>%
#   # take county wide mean centroid
#   group_by(county) %>%
#   # summarise(
#   mutate(
#     # across(c(rain, windspeed), mean, .names = "{.col}_mean") # , .groups = "drop"
#     across(c(contains("_a") | contains("_b")), mean, .names = "{.col}_mean") # , .groups = "drop"
#   ) %>%
#   ggplot(aes(x = rain, y = windspeed)) +
#   # county centroid points
#   # geom_point(
#   #   aes(x = rain_mean, y = windspeed_mean, colour = county),
#   #   size = 6,
#   #   pch = 10,
#   #   alpha = 0.9
#   # ) +
#   # site points
#   # geom_point(aes(colour = county), size = 3, alpha = 0.9) +
#   geom_point(size = 3, alpha = 0.9) +
#   geom_smooth(colour = "blue", method = "lm") +
#   # doesn't make much difference
#   # geom_smooth(aes(x = rain_mean, y = windspeed_mean), colour = "red", method = "lm") +
#   ggpubr::stat_cor(cor.coef.name = "rho", method = "pearson") +
#   scale_x_continuous(limits = c(0, 0.65), expand = c(0, 0)) +
#   # scale_y_continuous(limits = c(0, 1), expand = c(0, 0.1)) +
#   scale_y_continuous(limits = c(0, 0.8), expand = c(0, 0)) +
#   # scale_y_continuous(limits = c(0, 0.6), expand = c(0, 0)) + # for counties
#   labs(x = "a (precipitation)", y = "a (wind speed)") +
#   theme
#
# # ggsave("plots/p_compare_winter.png", p_compare, width = 6.3, height = 6, units = "in")
# ggsave("latex/plots/046_alpha_rain_vs_ws.png", p_compare, width = 6.3, height = 6, units = "in")

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
    facet_wrap(~var, ncol = 2) +
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
p_lst[[2]]

#### Apply Jensen-Shannon Divergence method to simulated data ####

#### Libs ####

library(sf)
library(dplyr, quietly = TRUE)
library(tidyr)
library(ggplot2)
devtools::load_all("texmex")

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

#### Metadata ####

vars <- c("rain", "wind_speed")

#### Load Data ####

data <- readRDS("data/sim_dat_3_clust.RDS")

# convert to dataframe
data_df <- bind_rows(lapply(seq_along(data), \(i) {
  n <- length(data[[i]])
  data.frame(
      "rain"       = data[[i]][1:(n / 2)],
      "wind_speed" = data[[i]][((n / 2) + 1):n]
  ) %>% 
    mutate(name = paste0("location_", i))
}))


#### Calculate Conditional Extremes parameters ####


# First, calculate threshold (90th quantile across all locations)
thresh <- apply(data_df[, 1:2], 2, quantile, 0.9)

# for each variable, calculate excess over threshold
data_thresh <- lapply(vars, \(x) {
  data_df %>% 
    select(matches(x), name) %>% 
    mutate(
      thresh = thresh[x],
      excess = !!sym(x) - thresh
    ) %>% 
    filter(excess > 0)
})

# Now fit evgam model for each marginal
# TODO: Is just fitting different models by loc appropriate?
evgam_fit <- lapply(data_thresh, \(x) {
  fit_evgam(
    data = x, 
    pred_data = x,
    # model scale and shape for each location
    f = list(excess ~ name, ~ name) 
  )
})

# Join scale and shape estimates into data
# TODO: Functionalise to work with > 2 variables
data_gpd <- distinct(data_df, name) %>% 
  bind_cols(
    rename(distinct(evgam_fit[[1]]$predictions), scale_rain = scale, shape_rain = shape),
    rename(distinct(evgam_fit[[2]]$predictions), scale_ws = scale, shape_ws = shape),
  )
  
# Now convert marginals to migpd (i.e. texmex format)
marginal <- gen_marg_migpd(data_gpd, data_df)
names(marginal) <- paste0(data_gpd$name, " - ", data_gpd$county)

# Calculate dependence from marginals
dependence <- fit_texmex_dep(
  marginal, 
  mex_dep_args = list(
    start = c(0.01, 0.01), 
    dqu = 0.7,
    fixed_b = FALSE,
    PlotLikDo = FALSE
  ) 
)

#### Calculate divergence upon which to cluster ####

# Perform PAM and k-means clustering, as an exploratory analysis

# pull parameter values for each location
params <- lapply(dependence, pull_params)

# pull Laplace threshold values for each location
thresh <- lapply(dependence, pull_thresh_trans)

# take maximum Laplace thresholds; want to geenra
# thresh_max <- apply(bind_rows(thresh), 2, max)
thresh_max <- lapply(bind_rows(thresh), max)

# list of locs -> list of len 2 of variables, each containing all locs
params <- purrr::transpose(params)

dist_mats <- lapply(seq_along(params), \(i) {
 proxy::dist(
   params[[i]], 
   method = js_div, 
   thresh_max = thresh_max[[i]], 
   data_max = 2 * thresh_max[[i]], 
   n_dat = 10
 )
})

# scree plots
lapply(dist_mats, scree_plot) # looks to be 3 clusters for both
lapply(dist_mats, scree_plot, fun = kmeans) 

# plot clustering for both rain and wind speed
ire_clust_plots <- lapply(dist_mats, \(x) {
  # for rain or wind speed, plot clustering based on PAM and k-means
  lapply(c(pam, kmeans), \(fun) {
    plt_clust(pts, fun(x, 3))
  })
})

# cluster adjacent sites only
dist_mats_adj <- lapply(dist_mats, \(x) {
  ret <- as.matrix(x)
  ret[adj_mat == 0] <- 1e9
  return(ret)
})

lapply(dist_mats_adj, scree_plot)
lapply(dist_mats_adj, scree_plot, fun = kmeans)

ire_clust_adj_plots <- lapply(dist_mats_adj, \(x) {
  # for rain or wind speed, plot clustering based on PAM and k-means
  lapply(c(pam, kmeans), \(fun) {
    plt_clust(pts, fun(x, 3))
  })
})

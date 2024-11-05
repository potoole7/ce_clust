#### Apply Jensen-Shannon Divergence method to simulated data ####

#### Libs ####

library(sf)
library(dplyr, quietly = TRUE)
library(tidyr)
library(ggplot2)
library(cluster)
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
# TODO: Problem, last threshold different to others, investigate
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

# scree plots look to suggest 3 clusters for both
lapply(dist_mats, scree_plot) 
lapply(dist_mats, scree_plot, fun = kmeans) 

# cluster for rain and wind speed using both k-means and PAM
js_clust <- lapply(dist_mats, \(x) {
  lapply(c(pam, kmeans), \(fun) {
    fun(x, 3)
  })
})

# Rand Index for clustering, works well for everything except k-means for wind
clust_rand <- lapply(js_clust, \(x) {
  lapply(x, \(y) {
    if (inherits(x, "pam")) {
      clust_sol <- y$clustering
    } else {
      clust_sol <- as.vector(y$cluster)
    }
    mclust::adjustedRandIndex(
      clust_sol,
      c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3)
    )
  })
})
unlist(clust_rand)

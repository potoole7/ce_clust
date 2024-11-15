#### Apply Jensen-Shannon Divergence method to simulated data ####

#### Libs ####

library(dplyr, quietly = TRUE)
library(tidyr)
library(ggplot2)
library(cluster)
library(evgam)
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


#### Metadata ####

vars <- c("rain", "wind_speed")
prob <- 0.9
cluster_mem <- sort(rep(c(1, 2, 3), 4))

#### Load Data ####

data_mix <- readRDS("data/sim_dat_3_clust.RDS")

# convert to dataframe
data_df <- bind_rows(lapply(seq_along(data_mix), \(i) {
  n <- length(data_mix[[i]])
  data.frame(
      "rain"       = data_mix[[i]][1:(n / 2)],
      "wind_speed" = data_mix[[i]][((n / 2) + 1):n]
  ) %>% 
    mutate(name = paste0("location_", i))
}))


#### Calculate Conditional Extremes parameters ####


# First, calculate threshold (90th quantile across all locations)
thresh <- apply(data_df[, 1:2], 2, quantile, prob)

# for each variable, calculate excess over threshold
data_thresh <- lapply(vars, \(x) {
  data_df %>% 
    dplyr::select(matches(x), name) %>% 
    mutate(
      thresh = thresh[x],
      excess = !!sym(x) - thresh
    ) %>% 
    filter(excess > 0)
})

# Now fit evgam model for each marginal
# TODO: Is just fitting different models by loc appropriate? (Yes I think!)
evgam_fit <- lapply(data_thresh, \(x) {
  fit_evgam(
    data = x, 
    pred_data = x,
    # model scale and shape for each location
    # f = list(excess ~ name, ~ name) 
    f = list(excess ~ name, ~ 1) # keep shape constant
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
    dqu = prob,
    fixed_b = FALSE,
    PlotLikDo = FALSE
  ), 
  fit_no_keef = TRUE # TODO: Had to unconstrain for some locations!
)

# check that all dependence models have run successfully
sapply(dependence, \(x) lapply(x, length))


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
   # data_max = 2 * thresh_max[[i]], 
   data_max = 5 * thresh_max[[i]],
   n_dat = 10
 )
})

# Sum distance matrices for different variables
dist_mat <- do.call(`+`, dist_mats)

# scree plots look to suggest 3 clusters for both
lapply(c(pam, kmeans), \(fun) scree_plot(dist_mat, fun = fun))

# cluster for rain and wind speed using both k-means and PAM
js_clust <- lapply(c(pam, kmeans), \(fun) fun(dist_mat, 3))

# Rand Index for clustering, works well for everything! 
clust_rand <- lapply(js_clust, \(x) {
  if (inherits(x, "pam")) {
    clust_sol <- x$clustering
  } else {
    clust_sol <- as.vector(x$cluster)
  }
  mclust::adjustedRandIndex(clust_sol, cluster_mem)
})
unlist(clust_rand)

# can also use function based on above code!
js_clust(dependence, cluster_mem = cluster_mem)

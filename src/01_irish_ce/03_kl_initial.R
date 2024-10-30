#### Investigatory clustering on a & b vals for CE model ####

# Note: Doesn't work for fixed b when keeping b values, so remove

# TODO: Investigate differences in length for reglines
# TODO: Get data for 6 counties

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

# for clustering
library(cluster)
library(proxy)

# for divergence/distance measures and test statistics
library(philentropy)
library(transport)
library(twosamples)

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


#### Load Data ####

# load original dataset
data <- readr::read_csv("data/met_eireann/final/met_eir_preprocess.csv.gz")

# load shapefile
areas <- sf::read_sf("data/met_eireann/final/irl_shapefile.geojson")

# Fitted conditional extremes models for each site
dependence <- readRDS("data/texmex_mexdep_obj.RDS")
# dependence <- readRDS("data/texmex_mexdep_obj_fixed_b.RDS")

# a and b values
ab_df <- readr::read_csv("data/ab_vals.csv")
# ab_df <- readr::read_csv("data/ab_vals_fixed_b.csv")

# extract point location of each station for plotting on map
pts <- ab_df %>% 
  distinct(lon, lat) %>% 
  st_to_sf()

# remove locations with NAs for a
na_locs <- ab_df %>% 
  filter(is.na(value)) %>% 
  pull(name) %>% 
  unique()

if (length(na_locs) > 0) {
  ab_df <- filter(ab_df, !name %in% na_locs)
  dependence <- dependence[
    !grepl(na_locs, names(dependence), fixed = TRUE)
  ]
}

# load adjacency matrix
adj_mat <- readRDS("data/ire_adj_mat.RDS")

# remove locations with NAs from data
data <- data %>% 
  semi_join(ab_df, by = "name")

# convert ab values to wide format
ab_df_wide <- ab_df %>% 
  mutate(col = paste0(var, "_", parameter)) %>% 
  select(name, county, col, value) %>% 
  pivot_wider(names_from = col, values_from = value)


#### Functions and testing ####

# pull a, b, mu, sigma for a given location
# dep - List of mex objects (i.e. CE models for each var) for single location
pull_params <- \(dep) {
  # fail if not list of mex objects for single location
  stopifnot(is.list(dep))
  stopifnot(all(vapply(dep, class, character(1)) == "mex"))
  # loop through conditioning variables for single location
  return(lapply(dep, \(x) {
    # pull parameter vals
    ret <- as.vector(x$dependence$coefficients)
    names(ret) <- rownames(x$dependence$coefficients)
    # remove c and d if 0 (i.e. Laplace margins, rather than Gumbel)
    if (ret["c"] == 0 && ret["d"] == 0) {
      ret <- ret[!names(ret) %in% c("c", "d")]
    } 
    return(ret)
  }))
}

# pull parameter values for each location
params <- lapply(dependence, pull_params)

# pull thresholds
# dep - List of mex objects (i.e. CE models for each var) for single location
pull_thresh_trans <- \(dep) {
  # fail if not list of mex objects for single location
  stopifnot(is.list(dep))
  stopifnot(all(vapply(dep, class, character(1)) == "mex"))
  # return quantile of transformed data (already calculated in texmex)
  # return(vapply(dep, \(x) x$dependence$dth, numeric(1)))
  return(lapply(dep, \(x) x$dependence$dth))
}

# pull Laplace threshold values for each location
thresh <- lapply(dependence, pull_thresh_trans)

# take maximum Laplace thresholds; want to geenra
thresh_max <- apply(bind_rows(thresh), 2, max)

# How best to split up KL divergence problem? Best to have:
# function to calculate KL divergence for single data point, 
# Function to calculate Jensen-Shannon from this
# Master function that applies this to all data, calculating mean and sd

# function to calculate KL divergence for single data point
# derivation at https://statproofbook.github.io/P/norm-kl.html
kl_gauss <- \(mu1, mu2, var1, var2) {
  return(1 / 2 * (
    (mu2 - mu1)^2 / var2 + (var1 / var2) - log(var1 / var2) - 1
  ))
}

# function to calculate Jensen-Shannon divergence metric from KL divergence
# This metric is symmetric, as required for clustering 
# https://tinyurl.com/526rwy9f
js_gauss <- \(mu1, mu2, var1, var2) {
  # calculate mean, variance for mixture distribution M = (P + Q)/2
  # Sum of normals as in https://online.stat.psu.edu/stat414/lesson/26/26.1
  # (N(u1, v1) + N(u2, v2)) / 2 = N((u1 + u2) / 2, (v1 + v2) / 2^2)
  mu_m <- (mu1 + mu2) / 2 
  var_m <- (var1 + var2) / 4
  
  # calculate JS(P||Q) = ((KL(P||M)) + KL(Q||M))/2
  return(
    (kl_gauss(mu1, mu_m, var1, var_m) + kl_gauss(mu2, mu_m, var1, var_m)) / 2
  )
}

# Function to calculate Jensen-Shannon divergence for each data point
params_x <- params[[1]]$rain
params_y <- params[[2]]$rain
data <- seq(
  thresh_max[[1]], 2 * thresh_max[[1]], length = 10 
)
# js_div <- \(params_x, params_y, data) {
js_div <- \(params_x, params_y, thresh_max, data_max = 2 * thresh_max, n_dat) {
  
  # create data sequence from specified arguments
  data <- seq(thresh_max, data_max, length = n_dat)
  
  # funs to calculate mu and sd for normal dist as in 5.2 of Heff & Tawn '04
  mu_fun <- \(x, data) {
    return(x[["a"]] * data + x[["m"]] * (data ^ x[["b"]]))
  }
  var_fun <- \(x, data) {
    # take square now to avoid in kl_single formula
    return((x[["s"]] * (data ^ x[["b"]]))^2)
  }
  
  # calculate mu and sigma for each data point
  mus <- lapply(list(params_x, params_y), mu_fun, data = data)
  vars <- lapply(list(params_x, params_y), var_fun, data = data)
  
  # Calculate Jensen-Shannon divergence for each data point
  # TODO: How best to summarise across all data points? Sum? Average?
  # return(mapply(
  #   js_gauss,
  #   mu1 = mus[[1]], mu2 = mus[[2]], var1 = vars[[1]], var2 = vars[[2]]
  # ))
  return(sum(mapply(
    js_gauss,
    mu1 = mus[[1]], mu2 = mus[[2]], var1 = vars[[1]], var2 = vars[[2]]
  )))
}


#### Clustering ####

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

ire_clust_plots <- lapply(dist_mats, \(x) {
  plt_clust(pts, pam(x, k = 3, diss = TRUE))
})

# cluster adjacent sites only
dist_mats_adj <- lapply(dist_mats, \(x) {
  ret <- as.matrix(x)
  ret[adj_mat == 0] <- 1e9
  return(ret)
})

lapply(dist_mats_adj, scree_plot)

ire_clust_adj_plots <- lapply(dist_mats_adj, \(x) {
  plt_clust(pts, pam(x, k = 3, diss = TRUE))
})

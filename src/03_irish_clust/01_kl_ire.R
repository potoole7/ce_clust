#### Conditional Extremes clustering for Irish data ####

# Note: Doesn't work for fixed b when keeping b values, so remove

# TODO: Investigate differences in length for reglines
# TODO: Get data for 6 counties

#### Libs ####

library(sf)
library(evc)
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
# find rough mean and variance of bulk rain and wind data across Ireland
data %>% 
  filter(rain != 0, wind_speed != 0) %>% 
  filter(rain < quantile(rain, 0.9), wind_speed < quantile(wind_speed, 0.9)) %>% 
  reframe(across(c("rain", "wind_speed"), \(x) c(mean(x), var(x))))

# load shapefile
areas <- sf::read_sf("data/met_eireann/final/irl_shapefile.geojson")

# Fitted conditional extremes models for each site
dependence <- readRDS("data/texmex_mexdep_obj.RDS")
# sort alphabetically by name
dependence <- dependence[sort(names(dependence))]
# dependence <- readRDS("data/texmex_mexdep_obj_fixed_b.RDS")

# a and b values
ab_df <- readr::read_csv("data/ab_vals.csv") %>% 
  arrange(name)
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
  dplyr::select(name, county, col, value) %>% 
  pivot_wider(names_from = col, values_from = value)


#### Clustering ####

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

# Sum distance matrices for different variables
dist_mat <- do.call(`+`, dist_mats)

# scree plots
scree_plot(dist_mat)
scree_plot(dist_mat, fun = kmeans)

# plot clustering for both rain and wind speed
# ire_clust_plots <- lapply(dist_mats, \(x) {
#   # for rain or wind speed, plot clustering based on PAM and k-means
#   lapply(c(pam, kmeans), \(fun) {
#     plt_clust_map(pts, fun(x, 3))
#   })
# })
# TODO: Any way to plot how likely points are to be in other clusters?
# Silhoutte plot??
ire_clust_plots <- lapply(c(pam, kmeans), \(fun) {
    plt_clust_map(pts, fun(dist_mat, 3))
})

pam_clust2 <- kl_sim_eval(
  dependence, kl_prob = 0.9, k = 2, dist_mat = dist_mat
)
plt_clust_map(pts, areas, pam_clust2)

# cluster adjacent sites only
#dist_mats_adj <- lapply(dist_mats, \(x) {
#  ret <- as.matrix(x)
#  ret[adj_mat == 0] <- 1e9
#  return(ret)
#})
dist_mat_adj <- as.matrix(dist_mat)
dist_mat_adj[adj_mat == 0] <- 1e9

scree_plot(dist_mat_adj)
scree_plot(dist_mat_adj, fun = kmeans)

ire_clust_adj_plots <- lapply(c(pam, kmeans), \(fun) {
    plt_clust_map(pts, fun(dist_mat_adj, 3))
})



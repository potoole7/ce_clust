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

#### Support functions ####

# compute the total within-cluster sum of distances
within_cluster_sum <- function(k, distance_matrix) {
  kmedoids_result <- pam(distance_matrix, k = k)
  return(kmedoids_result$objective[1])  # Total cost (sum of distances)
}

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

#### Functions ####

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
  return(vapply(dep, \(x) x$dependence$dth, numeric(1)))
}

# pull Laplace threshold values for each location
thresh <- lapply(dependence, pull_thresh_trans)

# take maximum Laplace thresholds; want to geenra
thresh_max <- apply(bind_rows(thresh), 2, max)

# Calculate 
to_name <- \(params, thresh_sup, n = 1000) {
  rnorm(n, )
}


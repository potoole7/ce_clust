#### Investigatory clustering on a & b vals for CE model ####

# Note: Doesn't work for fixed b when keeping b values, so remove

# TODO: Investigate differences in length for reglines

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

#### Functions ####

# compute the total within-cluster sum of distances
within_cluster_sum <- function(k, distance_matrix) {
  kmedoids_result <- pam(distance_matrix, k = k, diss = TRUE)
  return(kmedoids_result$objective[1])  # Total cost (sum of distances)
}

# Define the KL divergence function
# kl_divergence <- function(p, q) {
#   p <- p / sum(p)
#   q <- q / sum(q)
#   sum(p * log(p / q), na.rm = TRUE)
# }


#### Load Data ####

# load shapefile
areas <- sf::read_sf("data/met_eireann/final/irl_shapefile.geojson")

# Fitted conditional extremes models for each site
dependence <- readRDS("data/texmex_mexdep_obj_fixed_b.RDS")

# a and b values
# TODO: Save for fixed and varying b in 01_script
# Investigate NAs for Crolly (Filter Works) in Donegal
ab_df <- readr::read_csv("data/ab_vals_fixed_b.csv")

# remove locations with NAs for a
na_locs <- ab_df %>% 
  filter(is.na(value)) %>% 
  pull(name) %>% 
  unique()
if (length(na_locs) > 0) {
  ab_df <- filter(ab_df, !name %in% na_locs)
}

# convert ab values to wide format
ab_df_wide <- ab_df %>% 
  mutate(col = paste0(var, "_", parameter)) %>% 
  select(name, county, col, value) %>% 
  pivot_wider(names_from = col, values_from = value)

#### Calculate Voronoi cells and adjacency matrix ####

# extract point location of each station
pts <- ab_df %>% 
  select(lon, lat) %>% 
  distinct() %>% 
  st_to_sf()

# Calculate voronoi partition for sites
vor <- pts %>% 
  st_union() %>%
  st_voronoi(envelope = st_as_sfc(st_bbox(pts))) %>%
  st_collection_extract(type = "POLYGON") %>% 
  st_sf() %>%  # convert from geometry set to simple feature collection
  identity()

# order voronoi cells to match points
vor <- vor[order(st_nearest_feature(st_centroid(vor), pts)), ]

# check that voronoi cells have been produced correctly
plot(st_geometry(areas)) 
plot(vor, add = TRUE)
plot(pts, col = "blue", pch = 16, add = TRUE)
# test that voronoi cells and points ordered correctly
plot(vor[10, ], col = "blue", fill = NA, lwd = 5, add = TRUE)
plot(pts[10, ], col = "red", pch = 16, add = TRUE)

# calculate adjacency matrix from voronoi cells for stations
adj_mat <- spdep::nb2mat(spdep::poly2nb(vor), style = "B", zero.policy = TRUE) 


#### Simple k-mediod cluster on a and b values ####

# TODO: Add function for euclidean clustering (to repeat for varying b?)

# matrix of a and b values for rain and wind for each site
ab_mat <- ab_df_wide
# if all b values are "fixed", only cluster on a values
if (all(ab_df_wide$rain_b == ab_df_wide$rain_b[1])) {
  ab_mat <- select(ab_mat, contains("_a"))
} else {
  ab_mat <- select(ab_mat, contains("_a"), contains("_b"))
}
ab_mat <- as.matrix(ab_mat)

# dissimilarity matrix, by euclidean distance (simple, prob incorrect approach)
dist_mat_euc <- dist(ab_mat, method = "Euclidean")

# fit PAM clustering algorithm for k = 1:10
total_within_ss <- vapply(
  1:10, within_cluster_sum, dist_mat_euc, FUN.VALUE = numeric(1)
)

# Plot the scree plot
plot(
  1:10, total_within_ss, type = "b", pch = 19, 
  xlab = "Number of Clusters (k)", 
  ylab = "Total Within-Cluster Sum of Distances",
  main = "Scree Plot for K-medoids Clustering"
)
# Seems like there's a very slight elbow at k = 4
# Actually just do 2 to keep interpretation simple
pam_euclid_clust <- pam(dist_mat_euc, k = 3, diss = TRUE)

# plot 
#  Have different size/shape for cluster mediods
plt_clust <- \(pts, pam_clust) {
  pts_plt <- cbind(pts, data.frame("clust" = pam_clust$clustering)) %>% 
    mutate(row = row_number())
  pts_plt$mediod <- ifelse(pts_plt$row %in% pam_clust$medoids, TRUE, FALSE)

  ggplot(areas) + 
    geom_sf(colour = "black", fill = NA) + 
    geom_sf(
      data = pts_plt, 
      aes(colour = factor(clust), shape = mediod, size = as.numeric(mediod)), 
      alpha = 0.8
    ) + 
    scale_shape_discrete(breaks = c(1, 15)) + 
    scale_size_continuous(range = c(3, 6)) +  
    guides(shape = "none", size = "none") + 
    theme
}

plt_clust(pts, pam_euclid_clust)


#### Cluster adjacent sites only ####

# apply adjacency matrix to dissimilarity matrix 
dist_euc_adj <- as.matrix(dist_mat_euc)
dist_euc_adj[adj_mat == 0] <- 1e9

# fit PAM clustering algorithm for k = 1:10
total_within_ss <- vapply(
  1:10, within_cluster_sum, dist_euc_adj, FUN.VALUE = numeric(1)
)

# Plot the scree plot
plot(
  1:10, total_within_ss, type = "b", pch = 19, 
  xlab = "Number of Clusters (k)", 
  ylab = "Total Within-Cluster Sum of Distances",
  main = "Scree Plot for K-medoids Clustering"
)

# Seems like there's a very slight elbow at k = 2
pam_euc_adj_clust <- pam(dist_euc_adj, k = 8, diss = TRUE)

plt_clust(pts, pam_euc_adj_clust)

#### Cluster on regression line values with divergence "metrics" ####

# pull regression line for each location
reglines <- lapply(dependence, \(x) {
  lapply(x, \(y) as.vector(y$dependence$regline))
})

# TODO: Why are reglines different lengths for rain and wind?
# Only different for some locations! also none available for Crolly (Filter Works)
lens <- sapply(reglines, \(x) c(length(x$rain), length(x$wind_speed)))

# TODO: temp solution, for now, only keep for locations 
reglines <- reglines[lens[1, ] == lens[2, ] & lens[1, ] != 0]

# calculate divergences under various metrics:
# function to calculate KL divergence for shifted and normalised reg eqns
calc_kl <- \(x, y, unit = "log", ...) {
  shift_norm <- \(vec, epsilon = 1e-5) {
    # ensure strict positivity by shifting
    out <- vec
    if (any(vec < 0)) {
      out <- out + abs(min(out)) + epsilon
    }
    # normalise
    return(out / sum(out))
  }
  dat <- rbind(shift_norm(x), shift_norm(y))
  # KL(dat, unit = unit, ...)
  distance(dat, method = "kullback-leibler", ...)
}

# output various divergences
divergences <- bind_rows(lapply(reglines, \(x) {
  # list(
  data.frame(
    "kl"          = calc_kl(x[[1]], x[[2]], mute.message = TRUE)[[1]],
    "wasserstein" = wasserstein1d(x[[1]], x[[2]]),
    "ks"          = ks_stat(
      x[[1]] + rnorm(length(x[[1]]), 0, 1e-7), 
      x[[2]] + rnorm(length(x[[2]]), 0, 1e-7)
    ), 
    "cvm"         = cvm_stat(x[[1]], x[[2]]),
    "ad"          = ad_stat(x[[1]], x[[2]])
  )
}))

# perform clustering on each
clust_solns <- apply(divergences, 2, \(x) {
  pam(x, k = 3) # [c("clustering", "medoids")]
})

# plot for each
plots <- lapply(clust_solns, \(x) {
  plt_clust(pts, x)
})



#### Fixed b vs varying b ####
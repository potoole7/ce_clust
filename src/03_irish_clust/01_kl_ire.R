#### Investigatory and Vignotto 2021 clustering on Irish data ####

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

# load shapefile
areas <- sf::read_sf("data/met_eireann/final/irl_shapefile.geojson")

# a and b values from CE model (for simple clustering)
ab_df <- readr::read_csv("data/ab_vals.csv") %>%
  arrange(name)
# ab_df <- readr::read_csv("data/ab_vals_fixed_b.csv")

# remove locations with NAs for a
na_locs <- ab_df %>%
  filter(is.na(value)) %>%
  pull(name) %>%
  unique()

if (length(na_locs) > 0) {
  ab_df <- filter(ab_df, !name %in% na_locs)
}

# remove locations with NAs from data
data <- data %>%
  semi_join(ab_df, by = "name")

# convert ab values to wide format
ab_df_wide <- ab_df %>%
  mutate(col = paste0(var, "_", parameter)) %>%
  select(name, county, col, value) %>%
  pivot_wider(names_from = col, values_from = value)


#### Calculate Voronoi cells and adjacency matrix ####

# extract point location of each station
pts <- ab_df %>%
  distinct(lon, lat) %>%
  st_to_sf()

# calculate adjacency matrix from Voronoi cells
adj_mat <- calc_adj_mat(pts, plot = TRUE)
# save
saveRDS(adj_mat, "data/ire_adj_mat.RDS")


#### Simple k-mediod cluster on a and b values ####

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

# Elbow and silhouette plotting to determine number of clusters
scree_plot(dist_mat_euc) # indicates 2 has slight elbow
sil_boxplot(dist_mat_euc) # looks like 2 as well 

# plot
plt_clust_map(pts, areas, pam(dist_mat_euc, k = 2))
plt_clust_map(pts, areas, pam(dist_mat_euc, k = 3))


#### Cluster adjacent sites only ####

# apply adjacency matrix to dissimilarity matrix
dist_euc_adj <- as.matrix(dist_mat_euc)
dist_euc_adj[adj_mat == 0] <- 1e9

# scree plot
scree_plot(dist_euc_adj)
sil_boxplot(dist_euc_adj) # seems to choose 6

# TODO: Fix, ordering (or adjacency) seems wrong
plt_clust_map(pts, areas, pam(dist_euc_adj, k = 2))
plt_clust_map(pts, areas, pam(dist_euc_adj, k = 6))


#### Vignotto 2021 KL divergence clustering ####

# Split data into lists for each location as required by emp_kl_div
data_lst <- data %>%
  # split data by location
  group_split(name) %>%
  purrr::map(\(x) x %>%
               dplyr::select(rain, wind_speed) %>%
               stack() %>%
               dplyr::select(1) %>%
               as.vector() %>%
               `[[`(1))

# scree plot, looks like k = 3
prob <- 0.88 # fails for prob <- 0.9
# TODO: Investigate location?
kl_mat <- kl_sim_eval(data_lst, prob, k = NULL)$dist_mat

scree_plot(kl_mat, k = 1:5)
sil_boxplot(kl_mat, k = 2:5)

# plot results of clustering
plt_clust_map(pts, areas, kl_sim_eval(data_lst, prob, k = 2, kl_mat))
plt_clust_map(pts, areas, kl_sim_eval(data_lst, prob, k = 3, kl_mat))

plot(silhouette(kl_sim_eval(data_lst, prob, k = 2, kl_mat)))
plot(silhouette(kl_sim_eval(data_lst, prob, k = 3, kl_mat)))

# imposing adjacency
kl_mat_adj <- as.matrix(kl_mat)
kl_mat_adj[adj_mat == 0] <- 1e9

scree_plot(kl_mat_adj, k = 1:8)
sil_boxplot(kl_mat_adj, k = 1:8) # suggests k = 6

# plot results of clustering
# TODO: Fix, points don't look adjacent at all!
plt_clust_map(pts, areas, kl_sim_eval(data_lst, prob, k = 2, kl_mat_adj))
plt_clust_map(pts, areas, kl_sim_eval(data_lst, prob, k = 3, kl_mat_adj))
plt_clust_map(pts, areas, kl_sim_eval(data_lst, prob, k = 6, kl_mat_adj))

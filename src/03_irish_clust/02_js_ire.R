#### Conditional Extremes clustering for Irish data ####

# Note: Doesn't work for fixed b when keeping b values, so remove

# TODO: Investigate differences in length for reglines
# TODO: Get data for 6 counties

#### Libs ####

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
library(evgam)
library(geosphere)
library(circular)

# for clustering
library(cluster)
library(proxy)

# for divergence/distance measures and test statistics
library(philentropy)
library(transport)
library(twosamples)

source("src/functions.R")

sf::sf_use_s2(FALSE)


#### Load Data ####

# load original dataset
data <- readr::read_csv("data/met_eireann/final/met_eir_preprocess.csv.gz")

# pull county for locations
locs <- distinct(data, name, county)

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

# first check names for dependence are the same as in ab_df
unique(ab_df$name) == names(dependence) 
# extract point location of each station for plotting on map
pts <- ab_df %>% 
  distinct(name, lon, lat) %>% 
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

# scree plot (and distance matrix)
set.seed(123)
dist_mat <- js_clust(dependence, scree_k = 1:10)$dist_mat 
# no clear elbow unfortunately, perhaps could say at k = 4?

# look at boxplot of silhoutte widths for different values of k
sil_boxplot(dist_mat, k = 2:6) # looks like k = 2 is best by some distance

# Silhouette plot for different values of k (definitely best for k = 2)
plot(silhouette(pam(dist_mat, k = 2)))
plot(silhouette(pam(dist_mat, k = 3)))
plot(silhouette(pam(dist_mat, k = 4)))

# plot PAM clustering solution (for k = 2 and k = 3)
# TODO: Fix showing cluster centroid
pam_clust2 <- js_clust(dependence, k = 2, dist_mat = dist_mat)
plt_clust_map(pts, areas, pam_clust2)
pam_clust3 <- js_clust(dependence, k = 3, dist_mat = dist_mat)
plt_clust_map(pts, areas, pam_clust3)

# save clustering
saveRDS(pam_clust2, "data/pam_clust2.RDS")

# Plot silhouette coefficients on map
plt_sil_map(
  pts, 
  areas, 
  sil_obj = data.frame(silhouette(pam_clust2)) %>% 
    arrange(rownames(.)),
  medoids = rownames(pam_clust2$medoids)
)

plt_sil_map(
  pts, 
  areas, 
  sil_obj = data.frame(silhouette(pam_clust3)) %>% 
    arrange(rownames(.)),
  medoids = rownames(pam_clust3$medoids)
)

# plot vs covariates (alt, dist2coast, wind_dir)
covar_plt(pts, pam_clust2, data, areas)
covar_plt(pts, pam_clust3, data, areas)


#### Adjacency stipulation clustering ####

# repeat for adjacent sites only
dist_mat_adj <- as.matrix(dist_mat)
dist_mat_adj[adj_mat == 0] <- 1e9

sil_adj <- sil_boxplot(dist_mat_adj, k = 2:10)
sil_adj$plot + 
  scale_y_continuous(limits = c(0, 0.25)) 
sil_adj$sil %>% 
  group_by(k) %>% 
  summarise(sil_width = mean(sil_width)) %>% 
  arrange(desc(sil_width)) %>% 
  slice(1) %>% 
  pull(k)
# suggests k = 6! Valid to use for adjacent distance matrix?

pam_clust_adj2 <- js_clust(dependence, k = 2, dist_mat = dist_mat_adj)
plt_clust_map(pts, areas, pam_clust_adj2)
pam_clust_adj3 <- js_clust(dependence, k = 3, dist_mat = dist_mat_adj)
plt_clust_map(pts, areas, pam_clust_adj3)
pam_clust_adj6 <- js_clust(dependence, k = 6, dist_mat = dist_mat_adj)
plt_clust_map(pts, areas, pam_clust_adj6)

plot(silhouette(pam_clust_adj2))
plot(silhouette(pam_clust_adj3))
plot(silhouette(pam_clust_adj6))

plt_sil_map(
  pts, 
  areas, 
  data.frame(silhouette(pam_clust_adj2)) %>% 
    arrange(rownames(.)),
  medoids = rownames(pam_clust_adj2$medoids)
)

plt_sil_map(
  pts, 
  areas, 
  data.frame(silhouette(pam_clust_adj3)) %>% 
    arrange(rownames(.)),
  medoids = rownames(pam_clust_adj3$medoids)
)

plt_sil_map(
  pts, 
  areas, 
  data.frame(silhouette(pam_clust_adj6)) %>% 
    arrange(rownames(.)),
  medoids = rownames(pam_clust_adj6$medoids)
)

# plot vs covariates (alt, dist2coast, wind_dir)
covar_plt(pts, pam_clust_adj2, data, areas)
covar_plt(pts, pam_clust_adj3, data, areas)
covar_plt(pts, pam_clust_adj6, data, areas)

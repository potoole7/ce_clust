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
  kmedoids_result <- pam(distance_matrix, k = k)
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
  dependence <- dependence[
    !grepl(na_locs, names(dependence), fixed = TRUE)
  ]
}

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
calc_adj_mat <- \(pts, cut_vor = TRUE, plot = FALSE) {
  # Calculate voronoi partition for sites
  vor <- pts %>% 
    st_union() %>%
    st_voronoi(envelope = st_as_sfc(st_bbox(pts))) %>%
    # st_voronoi(envelope = st_as_sfc(st_bbox(areas))) %>%
    st_collection_extract(type = "POLYGON") %>% 
    st_sf() %>%  # convert from geometry set to simple feature collection
    identity()
  
  # order voronoi cells to match points
  vor <- vor[order(st_nearest_feature(st_centroid(vor), pts)), ]
  
  # cutoff voronoi cells from ocean, if desired (stops far away neighbours)
  # TODO: Cut via coastline
  if (cut_vor == TRUE) {
    # slightly smaller than areas bbox
    vor <- st_crop(
      vor, c("xmin" = -10, "ymin" = 51.6, "xmax" = -6, "ymax" = 55.2)
    )
  }
   
  # check that voronoi cells have been produced correctly
  if (plot == TRUE) {
    plot(st_geometry(areas)) 
    plot(vor, add = TRUE)
    plot(pts, col = "blue", pch = 16, add = TRUE)
    # test that voronoi cells and points ordered correctly
    # plot(vor[10, ], col = "blue", fill = NA, lwd = 5, add = TRUE)
    # plot(pts[10, ], col = "red", pch = 16, add = TRUE)
  }
  
  
  # calculate adjacency matrix from voronoi cells for stations
  return(spdep::nb2mat(spdep::poly2nb(vor), style = "B", zero.policy = TRUE))
}

adj_mat <- calc_adj_mat(pts, plot = TRUE)


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

# fit PAM for k = 1:10
scree_plot <- \(dist_mat) {
  total_within_ss <- vapply(
    1:10, within_cluster_sum, dist_mat, FUN.VALUE = numeric(1)
  )
  
  # scree plot
  plot(
    1:10, total_within_ss, type = "b", pch = 19 # , 
    # xlab = "Number of Clusters (k)", 
    # ylab = "Total Within-Cluster Sum of Distances",
    # main = "Scree Plot for K-medoids Clustering"
  )
}
scree_plot(dist_mat_euc)

# Seems like there's a very slight elbow at k = 4
pam_euclid_clust <- pam(dist_mat_euc, k = 4)

# plot 
#  Have different size/shape for cluster mediods
plt_clust <- \(pts, pam_clust) {
  medoids <- pam_clust$medoids
  if (inherits(medoids, "matrix")) medoids <- as.numeric(rownames(medoids))
  # browser()
  pts_plt <- cbind(pts, data.frame("clust" = pam_clust$clustering)) %>% 
    mutate(row = row_number())
  pts_plt$mediod <- ifelse(
    pts_plt$row %in% medoids, TRUE, FALSE
  )

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

# scree plot
scree_plot(dist_euc_adj)

# Seems like there's a very slight elbow at k = 7 
# (When largest cluster breaks up)
pam_euc_adj_clust <- pam(dist_euc_adj, k = 7)

plt_clust(pts, pam_euc_adj_clust)


#### Vignotto 2017 KL divergence clustering on regression lines ####

#### Pull through regression lines ####

# pull regression line for each location
reglines <- lapply(dependence, \(x) {
  lapply(x, \(y) as.vector(y$dependence$regline))
})

# TODO: Why are reglines different lengths for rain and wind?
# Only different for some locations! also none available for Crolly (Filter Works)
lens <- sapply(reglines, \(x) c(length(x$rain), length(x$wind_speed)))

# TODO: temp solution, for now, only keep for locations (41/58 locations)
bool <- (lens[1, ] == lens[2, ]) & (lens[1, ] != 0)
reglines <- reglines[bool]
pts_full <- pts[bool, ]

# split names from county
names_cnty <- stringr::str_split(names(reglines), pattern = " - ")
reglines_df <- bind_rows(lapply(seq_along(reglines), \(i) {
  data.frame(
    "name"       = names_cnty[[i]][[1]],
    "county"     = names_cnty[[i]][[2]],
    "rain"       = reglines[[i]]$rain, 
    "wind_speed" = reglines[[i]]$wind_speed)
}))

### Convert from Laplace to Pareto scale ####

# Convert from Laplace to observation


# Pareto scale: 1/(1 - ecdf(x))

#### Test calculating KL divergence between two locations ####


# function to split data into adjoint sets of extremal behaviour
split_sets <- \(x, thresh_rain, thresh_ws) {
  rain <- filter(x, rain > thresh_rain & wind_speed < thresh_ws)
  ws   <- filter(x, rain < thresh_rain & wind_speed > thresh_ws)
  both <- filter(x, rain > thresh_rain & wind_speed > thresh_ws)
  return(list(rain = rain, ws = ws, both = both))
}

# plot against each other for one location
test <- reglines_df %>% 
  filter(name == reglines_df$name[[1]])
# plot against each other for one location
plot(test$rain, test$wind_speed)
# add vertical and horizontal lines at 70th quantile
# TODO: Won't that just have the same number of points in rain and ws extreme sets?
thresh_rain <- quantile(test$rain, 0.7)
thresh_ws <- quantile(test$wind_speed, 0.7)
graphics::abline(v = thresh_rain, lty = 2)
graphics::abline(h = thresh_ws, lty = 2)

# split into 3 sets (extreme rain only, extreme ws only, both)
A <- split_sets(test, thresh_rain, thresh_ws)
points(A$rain$rain, A$rain$wind_speed, col = "red")
points(A$ws$rain, A$ws$wind_speed, col = "blue")
points(A$both$rain, A$both$wind_speed, col = "green")

# calculate empirical prop of data belonging to each of above sets
calc_emp_props <- \(A) {
  nrows <- vapply(A, nrow, numeric(1))
  return(nrows / sum(nrows))
}

# data for another location
test2 <- reglines_df %>% 
  filter(name == unique(reglines_df$name)[[2]])
A2 <- split_sets(
  test2, 
  quantile(test2$rain, 0.7), 
  quantile(test2$wind_speed, 0.7)
)

# calculate KL divergence between this site and another 
kl_divergence <- \(A1, A2) {
  # calculate empirical proportion of thresholded points in each set
  emp_prop1 <- calc_emp_props(A1)
  emp_prop2 <- calc_emp_props(A2)
  
  # KL divergence based on set proportions
  kl <- (1 / 2) * sum(vapply(seq_along(emp_prop1), \(i) {
    (emp_prop1[i] - emp_prop2[i]) * 
      (log(emp_prop1[i] / emp_prop2[i]))
  }, numeric(1)))
  
  return(kl)
}
kl_divergence(A, A2)


# Create wrapper function around each intermediate function above
kl_diss_fun <- \(df1, df2, prob = 0.7) {
  # browser()
  df_lst <- as.list(environment())[1:2] # pull dfs into list
  # TODO: Missing conversion step (Laplace -> Obs -> Pareto scale)
  # threshold and split into disjoint sets
  # TODO: Change thresholding to risk function as in Vignotto 2021
  A_lst <- lapply(df_lst, \(x) {
    split_sets(x, quantile(x$rain, prob), quantile(x$wind_speed, prob))
  })
  # calculate kl divergence between df1 and df2
  return(kl_divergence(A_lst[[1]], A_lst[[2]]))
}
kl_diss_fun(test, test2)

#### Calculate dissimilarity matrix based on KL divergence ####

locs <- unique(reglines_df$name)
# loop through each location
kl_div_lst <- lapply(locs, \(spec_loc) {
  # calculate KL divergence between spec_loc and other locs (including itself)
  kl_divs <- vapply(locs, \(other_loc) {
    kl_diss_fun(
      filter(reglines_df, name == spec_loc),
      filter(reglines_df, name == other_loc)
    ) 
  }, numeric(1))
})

kl_div_mat <- as.matrix(bind_rows(kl_div_lst))


#### Cluster with PAM and inspect results ####

# scree plot
scree_plot(kl_div_mat)

# elbow at k = 4
pam_kl <- pam(kl_div_mat, k = 4)
plt_clust(pts_full, pam_kl)

# do the same with adjacency matrix
adj_mat_full <- calc_adj_mat(pts_full)

kl_div_mat_adj <- kl_div_mat
kl_div_mat_adj[adj_mat_full == 0] <- 1e9

scree_plot(kl_div_mat_adj)

# Seems like there's a very slight elbow at k = 7 
# (When largest cluster breaks up)
pam_kl_adj <- pam(kl_div_mat_adj, k = 4)
plt_clust(pts_full, pam_kl_adj)


#### Cluster on regression line values with divergence "metrics" ####

# Old!

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
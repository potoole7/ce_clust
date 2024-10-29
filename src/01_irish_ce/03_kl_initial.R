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
# TODO: Easier to do as dataframe than in lists?
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
# thresh_max <- as.list(apply(bind_rows(thresh), 2, max))

# Function to sample from normal dist (likelihood for single conditional) for 
# Y_{-i}|Y_{i} for data exceeding global threshold
# TODO: Check normal distribution sampled from correctly, particularly logic
# params - CE parameters for one location
# thresh - Thresholds for one location
# thresh_max - Global maximum threshold
# data_max - Maximum data point to go up to, defaults as twice thresh_max
# data_len - How many data points between thresh_max and data_max should we 
# sample from a normal distribution for?
# n - Number of normal samples to generate for each data point
gen_ce_norm <- \(
  params, thresh, thresh_max, data_max = 2 * thresh_max, data_len = 10, n = 1000
) {
  # print(data_max)
  # print(data_len)
  # parameters
  # TODO: Can pull in a more programmatic way
  a <- params[["a"]] # alpha 
  b <- params[["b"]] # beta
  m <- params[["m"]] # mu nuisance parameter
  s <- params[["s"]] # sigma nuisance
  
  # Take sequence of length data_len from global max thresh to specified max
  # TODO: Should we sample from this using `sample()`?? Wouldn't be the same 
  # across different locations then, unless we did outside of function and/or 
  # had same seed
  data <- seq(thresh_max, data_max, length = data_len)
  
  # generate normal samples for each data point
  return(vapply(data, \(x) {
    # Likelihood for single conditional (5.2 H & T 2004)
    rnorm(
      n     = n,
      # mean  = a * x + m * (x ^ b),
      mean  = a * x + m * (x ^ b), #  + 10, # shift to the right -> non-neg
      sd    = s * (x ^ b)
    )
  }, FUN.VALUE = numeric(n)))
}

set.seed(123)
gen_ce_norm(
  params = params[[1]]$rain, 
  thresh = thresh[[1]][["rain"]], 
  thresh_max = thresh_max[["rain"]]
)

set.seed(123)
x <- gen_ce_norm(
  params = params[[1]]$rain, 
  thresh = thresh[[1]][["rain"]], 
  thresh_max = thresh_max[["rain"]]
)

set.seed(123)
y <- gen_ce_norm(
  params = params[[2]]$rain, 
  thresh = thresh[[2]][["rain"]], 
  thresh_max = thresh_max[["rain"]]
)

# Compare normal samples using KL divergence
# TODO: Could functionalise much of this??
# x - Matrix 
kl_compare <- \(x, y) {
  
  # convert x and y from matrix to vector form, if required
  x_fun <- x
  y_fun <- y 
  if (is.matrix(x)) {
    x_fun <- as.vector(x)
  }
  if (is.matrix(y)) {
    y_fun <- as.vector(y)
  }
  stopifnot(is.vector(x_fun) && is.vector(y_fun))
  
  # shift upwards if required, won't change shape of distribution
  # TODO: However, does appear to effect value of KL divergence!
  if (any(x_fun < 0) || any(y_fun < 0)) {
    epsilon <- 1E-8
    min_val <- min(c(x_fun, y_fun)) - epsilon
    message(paste0("Shifting by ", min_val, " so samples are non-negative"))
    # min_val <- -5 + epsilon
    # min_val <- -10 + epsilon
    x_fun <- x_fun - min_val
    y_fun <- y_fun - min_val
  }
  
  # divide by sum to get proper distributions
  x_fun <- x_fun / sum(x_fun)
  y_fun <- y_fun / sum(y_fun)
  
  plot(
    density(x_fun, bw = 3.807e-6), 
    xaxt = "n", 
    yaxt = "n", 
    xlim = c(-0.00005, 0.00025), 
    ylim = c(0, 40000)
  )
  # axis(1, at = seq(0, 0.00020, by = 0.00005))
  axis(1, at = seq(0, 0.00020, by = 0.0001))
  axis(2, at = seq(0, 40000, by = 10000))
  lines(density(y_fun), col = "red")
   
  # calculate JS divergence between x and y (symmetric and a proper metric)
  # 
  return(philentropy::distance(
    rbind(x_fun, y_fun), 
    # method = "kullback-leibler", 
    method = "jensen-shannon", 
    unit = "log", 
    mute.message = TRUE
  ))
}

#### Clustering ####

# for each location, generate normal samples
norm_samples_loc <- lapply(seq_along(params), \(i) {
  lapply(seq_along(params[[i]]), \(j) {
    gen_ce_norm(
      params     = params[[i]][[j]], 
      thresh     = thresh[[i]][[j]],
      thresh_max = thresh_max[[j]]
    )
  })
})

# pull for rain and wind speed
# norm_samples_rain <- lapply(norm_samples, \(x) as.vector(x[[1]]))
# norm_samples_ws <- lapply(norm_samples, \(x) as.vector(x[[2]]))

# list of locs -> list of len 2 of variables, each containing all locs
norm_samples_var <- purrr::transpose(norm_samples_loc)

# calculate distance matrices
# dist_rain <- dist(x = norm_samples_rain, method = kl_compare)
# dist_ws <- dist(x = norm_samples_rain, method = kl_compare)
dist_mats <- lapply(norm_samples_var, dist, method = kl_compare)

# scree plots
lapply(dist_mats, scree_plot)

# cluster based on both
# TODO: Not working at all! Need to fix
plt_clust(pts, pam(dist_mats[[2]], k = 2, diss = TRUE)) 

# cluster adjacent sites only
dist_mats_adj <- lapply(dist_mats, \(x) {
  ret <- as.matrix(x)
  ret[adj_mat == 0] <- 1e9
  return(ret)
})

lapply(dist_mats_adj, scree_plot)

plt_clust(pts, pam(dist_mats_adj[[1]], k = 3, diss = FALSE)) 

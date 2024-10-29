#### Generate  simulated dataset ####

#### libs ####
library(copula)
library(evd)
library(cluster)

# source functions
source("src/functions.R")
source("src/01_irish_ce/functions.R")


#### Gaussian mixture with t-copula ####


#### Gaussian distribution w/ resampling ####

# Want to create a mixture of Gaussian distributions and t-copula
# - For each location, need data for two variables, e.g. rain & wind speed
# - Gaussian data will be uncorrelated and make up bulk data
# - t-copulas will provide extreme data, correlation between variables will be 
# prespecified and used to define known cluster structure
# - "Normal" cluster methods must not perform well, while extreme methods 
# hopefully will


#### Gaussian distribution w/ resampling (OLD OLD OLD) ####

# Want to create a mixture of Gaussian distributions and t-copula
# - Gaussians will be uncorrelated and make up bulk data
# - t-copulas will provide extreme data and be correlated in known clusters

# metadata
seed_number <- 123
n <- 10000 # number of samples to take
n_locs <- 12 # number of "locations"
# cluster_mem <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3) # known cluster membership
cluster_mem <- c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3)
n_vars <- 2 # number of variables per location
# mix_p <- c("mvn" = 0.7, "gpd" = 0.3) # mixture percentages
mix_p <- c("mvn" = 0.5, "gpd" = 0.5) # mixture percentages
# Normal distribution parameters (using standard normal with 0 correlation)
mu <- rep(0, n_vars)
sigma <- diag(1, n_vars)
# sigma <- diag(0.5, n_vars)

# Generate MVN data
# data_mvn <- mvrnorm(n, mu, sigma)
# list of locs, each with mvn distribution with cols for each var
set.seed(seed_number)
data_mvn <- lapply(seq_len(n_locs), \(i) {
  MASS::mvrnorm(n, mu, sigma)
})

vapply(data_mvn, \(x) round(cor(x)[1, 2], 3), numeric(1))

# t-copula correlation matrix
# rho <- matrix(0, nrow = n_locs, ncol = n_locs) # initialise
# initialise
rho <- list(
  "cluster_1" = sigma,
  "cluster_2" = sigma,
  "cluster_3" = sigma
)

# separate into known clusters by changing correlation
# rho[1:3, 1:3] <- 0.8  
# rho[4:6, 4:6] <- 0.8  
# rho[7:10, 7:10] <- 0.8  
# diag(rho) <- 1
rho[[1]][1, 2] <- rho[[1]][2, 1] <- 0.95
rho[[2]][1, 2] <- rho[[2]][2, 1] <- 0.5
rho[[3]][1, 2] <- rho[[3]][2, 1] <- 0.05

# create t-copula with above correlation structure
# t_cop <- tCopula(
#   param = rho[lower.tri(rho)], 
#   dim = n_locs, 
#   df = 1, # low df <=> heavier tails
#   dispstr = "un"
# )
set.seed(seed_number)
t_cop <- lapply(seq_len(n_locs), \(i) {
  
  # switch rho such that locs 1-3, locs 4-6 and locs 7-10 form three clusters
  rho_spec <- rho[[1]]
  # if (i %in% 4:6) {
  if (i %in% 5:8) {
    rho_spec <- rho[[2]]
  # } else if (i %in% 7:10) {
  } else if (i %in% 9:12) {
    rho_spec <- rho[[3]]
  }
  # print(rho_spec)
  
  # create t-copula
  tCopula(
    param = rho_spec[lower.tri(rho_spec)], 
    dim = n_vars, 
    df = 1, # low df <=> heavier tails
    dispstr = "un"
  )
})

# Generate uniform samples from copula
# u <- rCopula(n, t_cop)  
# # Transform to GPD marginals 
# # TODO: Transform (all? some?) to Laplace margins
# data_gpd <- qgpd(u, loc = 0, scale = 1, shape = 0.2)

set.seed(seed_number)
data_gpd <- lapply(t_cop, \(x) {
  # Generate uniform samples from copula
  u <- rCopula(n, x)  
  # Transform to GPD marginals 
  # TODO: Transform (all? some?) to Laplace margins
  # TODO: Change location parameter? I think this is where I'm going wrong
  # evd::qgpd(u, loc = 0, scale = 1, shape = 0.2)
  evd::qgpd(u, loc = 3, scale = 1, shape = 0.5)
})

# test: Correlation for GPD margins follows known correlation structure
# proxy::simil(lapply(seq_len(ncol(data_gpd)), function(i) data_gpd[, i]))
vapply(data_gpd, \(x) round(cor(x)[1, 2], 3), numeric(1))

# mixture model
# data_mix <- rbind(
#   # sample from normal distributions
#   apply(data_mvn, 2, sample, size = n * mix_p[[1]]),
#   # sample from GPD margins
#   apply(data_gpd, 2, sample, size = n * mix_p[[2]])
# )
# TODO: Change to Laplace margins?? will happen in Vignotto method anyway
set.seed(seed_number)
data_mix <- lapply(seq_along(data_mvn), \(i) {
  # TODO: Make programmatic
  mix_p[[1]] * data_mvn[[i]] + 
    mix_p[[2]] * data_gpd[[i]]
})

plot(data_mix[[1]])
# should within reason remove most of normal samples
apply(data_mvn[[1]], 2, max)
apply(data_gpd[[1]], 2, quantile, 0.9)
apply(data_mix[[1]], 2, quantile, 0.9)

# TODO: Plot for each location
# check that there is little correlation here
vapply(data_mix, \(x) round(cor(x)[1, 2], 3), numeric(1))

# TODO: cluster using "normal" methods 
# (i.e. dissimilarity matrix based on correlation between variables)

# cluster using Vignotto 2021 method:
# split data into list of locations where variables are stacked
data_lst <- lapply(data_mix, \(x) c(x[, 1], x[, 2]))

# Calculate dissimilarity matrix 
# prob <- 0.8
prob <- 0.9
kl_mat <- proxy::dist(
  data_lst, method = emp_kl_div, print = FALSE, prob = prob
)

scree_plot(kl_mat, 1:5) # obvious choice seems to be 3

# cluster based on dissimilarity
pam_kl_clust <- pam(kl_mat, k = 3)
pam_kl_clust$clustering

#### Generate  simulated dataset ####

#### libs ####
library(copula)
library(evd)
library(cluster)
library(dplyr)
# devtools::install_github("r-forge/lcmix/pkg")
library(lcmix)

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


#### Gaussian/Gamma distribution w/ resampling (OLD OLD OLD) ####

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
# mix_p <- c("mvn" = 0.7, "gpd" = 0.3) # mixture percentages
mix_p <- c("mvn" = 0.3, "gpd" = 0.7) # mixture percentages
# Normal distribution parameters (using standard normal with 0 correlation)
# mu <- rep(0, n_vars)
# sigma <- diag(1, n_vars)
# sigma <- diag(0.5, n_vars)
# Gamma parameters
# gam_mean <- c(23.2, 6.44)
# gam_var <- c(304, 2.17)
# parameter values roughly coming from Ireland dataset
# gam_mean <- 20
# gam_var <- 200
gam_mean <- 2
gam_var <- 10
alpha <- (gam_mean ** 2) / gam_var
beta <- (gam_mean) / gam_var

# Generate MVN data
# data_mvg <- mvrnorm(n, mu, sigma)
# list of locs, each with mvn distribution with cols for each var
set.seed(seed_number)
data_mvg <- lapply(seq_len(n_locs), \(i) {
  rmvgamma(n, alpha, beta, corr = diag(rep(1, n_vars)))
})

par(mfrow = c(2, 2))
# lapply(data_gpd[c(1, 5, 10)], plot)
lapply(data_mvg[c(1, 5, 7, 12)], plot)
par(mfrow = c(1, 1))


vapply(data_mvg, \(x) round(cor(x)[1, 2], 3), numeric(1))

# t-copula correlation matrix
# rho <- matrix(0, nrow = n_locs, ncol = n_locs) # initialise
# initialise
id_mat <- diag(1, n_vars)
rho <- list(
  "cluster_1" = id_mat,
  "cluster_2" = id_mat # ,
  # "cluster_3" = sigma
)

# separate into known clusters by changing correlation
# rho[1:3, 1:3] <- 0.8  
# rho[4:6, 4:6] <- 0.8  
# rho[7:10, 7:10] <- 0.8  
# diag(rho) <- 1
rho[[1]][1, 2] <- rho[[1]][2, 1] <- 0.95
# rho[[2]][1, 2] <- rho[[2]][2, 1] <- 0.65
# rho[[3]][1, 2] <- rho[[3]][2, 1] <- 0.01

rho[[2]][1, 2] <- 0.1

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
  df <- 1
  rho_spec <- rho[[1]]
  # for three clusters
  # if (i %in% 5:8) {
  #   rho_spec <- rho[[2]]
  # } else if (i %in% 9:12) {
  #   rho_spec <- rho[[3]]
  # }
  # for two clusters
  if (i %in% 7:12) {
    df <- 7
    rho_spec <- rho[[2]]
    # return(data_mvg[[i]])
  }
  
  # create t-copula (Gaussian copula for other cluster??)
  tCopula(
    param = rho_spec[lower.tri(rho_spec)], 
    dim = n_vars, 
    df = df, # low df <=> heavier tails
    dispstr = "un"
  )
})

# Generate uniform samples from copula
# u <- rCopula(n, t_cop)  
# # Transform to GPD marginals 
# # TODO: Transform (all? some?) to Laplace margins
# data_gpd <- qgpd(u, loc = 0, scale = 1, shape = 0.2)

set.seed(seed_number)
data_gpd <- lapply(seq_along(t_cop), \(i) {
  # TODO: Most basic of testing, no "extreme" data from t-Copula for 2nd clust
  # if (i %in% 7:12) return(data_mvg[[i]])
  # Generate uniform samples from copula
  u <- rCopula(n, t_cop[[i]])  
  # Transform to GPD marginals 
  # TODO: Transform (all? some?) to Laplace margins
  # TODO: Change location parameter? I think this is where I'm going wrong
  # evd::qgpd(u, loc = 0, scale = 1, shape = 0.2)
  # evd::qgpd(u, loc = 3, scale = 1, shape = 0.05)
  # evd::qgpd(u, loc = 3, scale = 1, shape = 0.01)
  # evd::qgpd(u, loc = 0, scale = 3, shape = -0.5)
  
  # for Gamma (parameters chosen so max is ~ 300, in line with Irish data)
  # evd::qgpd(u, loc = quantile(data_mvg[[i]], 0.99), scale = 3, shape = 0.3)
  # evd::qgpd(u, loc = quantile(data_mvg[[i]], 0.99), scale = 3, shape = -0.05)
  evd::qgpd(
    u, 
    # loc = quantile(data_mvg[[i]], 0.99),
    loc = 0.99 * max(data_mvg[[i]]), 
    scale = 3, 
    # shape = 0.3
    shape = -0.05
  )
})

# test: Correlation for GPD margins follows known correlation structure
# proxy::simil(lapply(seq_len(ncol(data_gpd)), function(i) data_gpd[, i]))
vapply(data_gpd, \(x) round(cor(x)[1, 2], 3), numeric(1))

# Plot var1 vs var2 for first, second and third clusters
# TODO: Appears to still be many co-occurring extremes for cluster 3
par(mfrow = c(2, 2))
# lapply(data_gpd[c(1, 5, 10)], plot)
lapply(data_gpd[c(1, 5, 7, 12)], plot)
par(mfrow = c(1, 1))

# mixture model
# data_mix <- rbind(
#   # sample from normal distributions
#   apply(data_mvg, 2, sample, size = n * mix_p[[1]]),
#   # sample from GPD margins
#   apply(data_gpd, 2, sample, size = n * mix_p[[2]])
# )
set.seed(seed_number)
data_mix <- lapply(seq_along(data_mvg), \(i) {
  mix_p[[1]] * data_mvg[[i]] + 
    mix_p[[2]] * data_gpd[[i]]
})

par(mfrow = c(2, 2))
# lapply(data_mix[c(1, 5, 10)], plot)
lapply(data_mix[c(1, 5, 7, 12)], plot)
par(mfrow = c(1, 1))

# should within reason remove most of normal samples
# apply(data_mvg[[1]], 2, max)
# apply(data_gpd[[1]], 2, quantile, 0.9)
apply(data_mix[[1]], 2, quantile, c(0.9, 0.95, 0.99))

# TODO: Plot for each location
# check that there is little correlation here
vapply(data_mix, \(x) round(cor(x)[1, 2], 3), numeric(1))

# TODO: cluster using "normal" methods 
# (i.e. dissimilarity matrix based on correlation between variables)

# cluster using Vignotto 2021 method:
# split data into list of locations where variables are stacked
# TODO: May be error here!! Plots don't match
data_lst <- lapply(data_mix, as.vector)

# Calculate dissimilarity matrix 
# prob <- 0.8
prob <- 0.95
kl_mat <- proxy::dist(
  data_lst, method = emp_kl_div, print = FALSE, prob = prob
)

# save data
saveRDS(data_lst, file = "data/sim_dat_2_clust.RDS")

scree_plot(kl_mat, 1:5) # choice of 2/3 seems correct

# walk through function
source("src/01_irish_ce/functions.R")
# debugonce(emp_kl_div)
# debugonce(pareto_trans)
# debugonce(partition_max)
emp_kl_div(
  data_lst[[1]], 
  # data_lst[[2]], 
  data_lst[[7]], 
  # prob = prob, 
  prob = 0.95, 
  plot = TRUE, 
  print = TRUE
)

# cluster based on dissimilarity
set.seed(seed_number)
# pam_kl_clust <- pam(kl_mat, k = 3)
pam_kl_clust <- pam(kl_mat, k = 2)
pam_kl_clust$clustering

mclust::adjustedRandIndex(
  pam_kl_clust$clustering, 
  # c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3)
  c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2)
)

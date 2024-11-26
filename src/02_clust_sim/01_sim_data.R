#### Generate simulated dataset ####

#### libs ####
library(copula)
library(evc)
library(evd)
library(cluster)
library(dplyr)
# devtools::install_github("r-forge/lcmix/pkg")
library(lcmix)

# source functions
source("src/functions.R")


# Want to create a mixture of Gamma distributions and t-copula
# - For each location, need data for two variables, e.g. rain & wind speed
# - Gamma data will be uncorrelated and make up bulk data
# - t-copulas will provide extreme data, correlation between variables will be 
# pre-specified and used to define known cluster structure


#### MV Gamma ####

# Want to create a mixture of Gamma distributions and t-copula
# - Gaussians will be uncorrelated and make up bulk data
# - t-copulas will provide extreme data and be correlated in known clusters

# metadata
seed_number <- 123
n <- 10000 # number of samples to take
n_locs <- 12 # number of "locations"
# cluster_mem <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3) # known cluster membership
cluster_mem <- c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3)
n_vars <- 2 # number of variables per location
mix_p <- c("mvn" = 0.3, "gpd" = 0.7) # mixture percentages
# GPD parameters
scale_gpd <- 6 
shape_gpd <- -0.05
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
gam_mean <- 1
gam_var <- 5
alpha <- (gam_mean ** 2) / gam_var
beta <- (gam_mean) / gam_var

# Generate MV Gamma data for two variables at each location
set.seed(seed_number)
data_mvg <- lapply(seq_len(n_locs), \(i) {
  rmvgamma(n, alpha, beta, corr = diag(rep(1, n_vars)))
})

par(mfrow = c(2, 2))
lapply(data_mvg[c(1, 5, 10)], plot)
# lapply(data_mvg[c(1, 5, 7, 12)], plot)
par(mfrow = c(1, 1))

# confirm bulk data for variables at each location have no correlation 
# TODO: Is this realistic? Won't they have some correlation?
vapply(data_mvg, \(x) round(cor(x)[1, 2], 3), numeric(1))

#### t-Copula w/ GPD margins ####

# t-copula correlation matrix
# initialise
id_mat <- diag(1, n_vars)
rho <- list(
  "cluster_1" = id_mat,
  "cluster_2" = id_mat,
  "cluster_3" = id_mat
)

# separate into known clusters by changing correlation
# rho[1:3, 1:3] <- 0.8  
# rho[4:6, 4:6] <- 0.8  
# rho[7:10, 7:10] <- 0.8  
# diag(rho) <- 1
rho[[1]][1, 2] <- rho[[1]][2, 1] <- 0.95
# rho[[2]][1, 2] <- 0.1
rho[[2]][1, 2] <- rho[[2]][2, 1] <- 0.5
rho[[3]][1, 2] <- rho[[3]][2, 1] <- 0.1

# create t-copula with above correlation structure
# t_cop <- tCopula(
#   param = rho[lower.tri(rho)], 
#   dim = n_locs, 
#   df = 1, # low df <=> heavier tails
#   dispstr = "un"
# )
# TODO: DF for t-copula very important, determine how well clustering recovered
set.seed(seed_number)
t_cop <- lapply(seq_len(n_locs), \(i) {
  
  # switch rho such that locs 1-3, locs 4-6 and locs 7-10 form three clusters
  df <- 1
  rho_spec <- rho[[1]]
  # for three clusters
  if (i %in% 5:8) {
    df <- 5
    rho_spec <- rho[[2]]
  } else if (i %in% 9:12) {
    df <- 10
    rho_spec <- rho[[3]]
  }
  # for two clusters
  # if (i %in% 7:12) {
  #   df <- 7
  #   rho_spec <- rho[[2]]
  # }
  
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
  evd::qgpd(
    u, 
    # set threshold
    # loc = quantile(data_mvg[[i]], 0.99),
    loc = max(data_mvg[[i]]),
    scale = scale_gpd, 
    shape = shape_gpd
  )
})

# test: Correlation for GPD margins follows known correlation structure
# proxy::simil(lapply(seq_len(ncol(data_gpd)), function(i) data_gpd[, i]))
vapply(data_gpd, \(x) round(cor(x)[1, 2], 3), numeric(1))

# Plot var1 vs var2 for first, second and third clusters
# TODO: Appears to still be many co-occurring extremes for cluster 3
par(mfrow = c(2, 2))
lapply(data_gpd[c(1, 5, 10)], plot)
# lapply(data_gpd[c(1, 5, 7, 12)], plot)
par(mfrow = c(1, 1))

#### Mixture Model ####

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
lapply(data_mix[c(1, 5, 10)], plot)
# lapply(data_mix[c(1, 5, 7, 12)], plot)
par(mfrow = c(1, 1))

# TEMP: Have data as purely from t-copula with GPD margins
# data_mix <- data_gpd

# check correlation after combining
vapply(data_mix, \(x) round(cor(x)[1, 2], 3), numeric(1))


#### cluster using Vignotto 2021 method ####

# split data into list of locations where variables are stacked
data_lst <- lapply(data_mix, as.vector)

# elbow plot
kl_mat <- kl_sim_eval(data_lst, kl_prob = 0.9)[[1]] # clear elbow at k = 3

# cluster at k = 3 and evaluate performance
(pam_clust <- kl_sim_eval(
  data_lst, 
  kl_prob = 0.9, 
  k = 3, 
  dist_mat = kl_mat, 
  cluster_mem = c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3)
))

# also look at silhouette summary and plot of clustering solution
sil <- silhouette(pam_clust$pam)
summary(sil)
plot(sil)

# boxplot for silhouettes, shows k = 3 is best, as expected
sil <- sil_boxplot(dist_mat, 2:5)

# save data
# saveRDS(data_lst, file = "data/sim_dat_2_clust.RDS")
saveRDS(data_lst, file = "data/sim_dat_3_clust.RDS")

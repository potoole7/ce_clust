#### Simulation sensitivity analysis of Vignotto 2021 ####

#### libs ####

library(copula)
library(evd)
library(cluster)
library(dplyr)
library(tidyr)
library(lcmix)
library(ggplot2)

# source functions
source("src/functions.R")
source("src/01_irish_ce/functions.R")


#### Metadata ####

seed_number <- 123
n <- 10000 # number of samples to take
n_locs <- 12 # number of "locations"
cluster_mem <- sort(rep(1:2, n_locs / 2)) # two clusters for half of locations each
n_vars <- 2 # number of variables per location
mix_p <- c("gauss_cop" = 0.3, "t_cop" = 0.7) # mixture percentages
# Normal distribution parameters (using standard normal with 0 correlation)
mu <- rep(0, n_vars)
sigma <- diag(1, n_vars)
# GPD parameters
scale_gpd <- 1
shape_gpd <- -0.05


#### Simulate data ####

# Function parameters
# correlation parameters for each cluster
cor_norm <- c(0.1, 0.1)
cor_t <- c(0.9, 0.1)
# degrees of freedom for t-Copula
df_t <- c(1, 5)
# quantile for Vignotto 2021 method
prob <- 0.9

# simulate from Gaussian copula with Gaussian margins
# TODO: How does using Gaussian copula differ from MVN distribution?
data_norm <- lapply(seq_len(n_locs), \(i){
  cor <- cor_norm[[1]]
  if (i > floor(n_locs / 2)) {
    cor <- cor_norm[[2]]
  }
  
  cop_norm <- normalCopula(cor, dim = n_vars, dispstr = "un")
  u <- rCopula(n, cop_norm)
  return(qnorm(u, mean = mu, sd = sigma))
})

# check correlation is correct
vapply(data_norm, \(x) round(cor(x)[1, 2], 3), numeric(1))

# simulate from t-Copula with GPD margins
data_gpd <- lapply(seq_len(n_locs), \(i) {
  cor <- cor_t[[1]]
  df <- df_t[[1]]
  if (i > floor(n_locs / 2)) {
    cor <- cor_t[[2]]
  df <- df_t[[2]]
  }
  cop_t <- copula::tCopula(cor, dim = n_vars, df = df_t[[2]], dispstr = "un")
  u <- rCopula(n, cop_t)
  return(evd::qgpd(
    p     = u,
    loc   = max(data_norm[[i]]), 
    scale = scale_gpd, 
    shape = shape_gpd
  ))
})
vapply(data_gpd, \(x) round(cor(x)[1, 2], 3), numeric(1))

# mixture
data_mix <- lapply(seq_len(n_locs), \(i){
  mix_p[[1]] * data_norm[[i]] + mix_p[[2]] * data_gpd[[i]]
})

vapply(data_mix, \(x) round(cor(x)[1, 2], 3), numeric(1))

# check plots
par(mfrow = c(2, 2))
lapply(data_mix[c(1, 5, 7, 12)], plot)
par(mfrow = c(1, 1))

# KL divergence between areas using Vignotto 2021 method
kl_mat <- proxy::dist(
  data_mix, method = emp_kl_div, print = FALSE, prob = prob
)

# clustering solution
pam_kl_clust <- pam(kl_mat, k = 2)
# evaluate quality
mclust::adjustedRandIndex(
  pam_kl_clust$clustering, 
  cluster_mem
)

kl_sim_eval <- \(
  n_locs = 12,     # number of locations
  n = 1e4,         # number of simulations
  cor_norm,        # bulk correlation for two clusters from Gaussian copula
  params_norm,     # normal marginal parameters (same for both)
  cor_t,           # extreme correlation for two clusters from t-copula
  df_t,            # degrees of freedom for t-copula
  params_gpd,      # GPD margin parameters
  mix_p,           # mixture weights
  kl_prob,         # Extremal quantile
  ret_data = FALSE # Whether to return data (may be memory intensive)
) {
  # many arguments must be of length 2 for 2 clusters
  stopifnot(all(vapply(
    list(cor_norm, params_norm, cor_t, df_t, params_gpd, mix_p), 
    \(x) length(x) == 2, logical(1)
  )))
  # prob must be valid probability
  stopifnot(0 <= kl_prob && kl_prob <= 1)
  stopifnot(sum(mix_p) == 1)
  
  # Simulate from Gaussian Copula with normal margins
  data_norm <- lapply(seq_len(n_locs), \(i){
    cor <- cor_norm[[1]]
    if (i > floor(n_locs / 2)) {
      cor <- cor_norm[[2]]
    }
    
    cop_norm <- normalCopula(cor, dim = n_vars, dispstr = "un")
    u <- rCopula(n, cop_norm)
    return(qnorm(u, mean = params_norm[1], sd = params_norm[2]))
  })

  # simulate from t-Copula with GPD margins
  data_gpd <- lapply(seq_len(n_locs), \(i) {
    cor <- cor_t[[1]]
    df <- df_t[[1]]
    if (i > floor(n_locs / 2)) {
      cor <- cor_t[[2]]
    df <- df_t[[2]]
    }
    cop_t <- copula::tCopula(cor, dim = n_vars, df = df_t[[2]], dispstr = "un")
    u <- rCopula(n, cop_t)
    return(evd::qgpd(
      p     = u,
      loc   = max(data_norm[[i]]), 
      scale = scale_gpd, 
      shape = shape_gpd
    ))
  })
  
  # mixture
  data_mix <- lapply(seq_len(n_locs), \(i){
    mix_p[[1]] * data_norm[[i]] + mix_p[[2]] * data_gpd[[i]]
  })
  
  # KL divergence between areas using Vignotto 2021 method
  kl_mat <- proxy::dist(
    data_mix, method = emp_kl_div, print = FALSE, prob = kl_prob
  )
  
  # clustering solution
  pam_kl_clust <- pam(kl_mat, k = 2)
  # evaluate quality
  adj_rand <- mclust::adjustedRandIndex(
    pam_kl_clust$clustering, 
    cluster_mem
  )
  print(paste0("adjusted Rand index" = adj_rand))
  
  # return object
  ret <- list(
    "data"     = list(
      "data_norm" = data_norm, 
      "data_gpd"  = data_gpd,
      "data_mix"  = data_mix
    ), 
    "pam"      = pam_kl_clust,
    "adj_rand" = adj_rand
  )
  
  # returning data each time may be quite memory intensive!
  if (ret_data == FALSE) {
    ret <- ret[names(ret) != "data"]
  }
  return(ret)
}

# test (works for kl_prob from 0.9 to 0.99!)
set.seed(seed_number)
kl_sim_eval(
  n_locs = 12, 
  n = 1e4, 
  cor_norm = c(0.5, 0.5), 
  cor_t = c(0.8, 0.3),
  params_norm = c(0, 1),
  df_t = c(1, 5),
  params_gpd = c(1, -0.05),
  mix_p = c(0.5, 0.5),
  kl_prob = 0.99
)


#### Grid search ####

# create grid 
# TODO: Reduce number of parameters (Currently 128k, can't possibly run so many times)
grid <- tidyr::crossing(
  # Don't need all cor_vals for both (i.e. (0.3, 0.8) is the same as (0.8, 0.3)
  cor_norm1 = seq(0, 1, by = 0.1),
  cor_norm2 = seq(0, 1, by = 0.1),
  # cor_norm1 = 0.1,
  # cor_norm2 = 0.1,
  # cor_norm = seq(0, 1, by = 0.1),
  # t-copula correlation
  cor_t1    = seq(0.1, 0.9, by = 0.2),
  cor_t2    = seq(0.2, 1, by = 0.2),
  # Degrees of freedom for t-copula
  # df_t1     = c(1, 5, 10),
  # df_t2     = c(1, 5, 10),
  df_t1   = 1,
  df_t2   = 1,
  # mixture percentages (must sum to 1)
  mix_p1  = 0.5,
  mix_p2  = 0.5,
  # extremal quantiles
  # kl_prob = c(0.9, 0.95, 0.99)
  kl_prob = 0.9
) %>% 
  filter(
    # use same correlation in both clusters for Gaussian copula
    cor_norm1 == cor_norm2,
    # mixture weights must sum to 1
    mix_p1 + mix_p2 == 1
  )


#### Plot ####




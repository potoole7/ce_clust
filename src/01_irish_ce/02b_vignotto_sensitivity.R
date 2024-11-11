#### Simulation sensitivity analysis of Vignotto 2021 ####

#### libs ####

library(copula)
library(evd)
library(cluster)
library(dplyr)
library(tidyr)
library(lcmix)
library(ggplot2)
library(parallel)

# source functions
source("src/functions.R")
source("src/01_irish_ce/functions.R")

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


#### Metadata ####

seed_number <- 123
n <- 1e3 # number of samples to take
n_locs <- 12 # number of "locations"
cluster_mem <- sort(rep(1:2, n_locs / 2)) # two clusters for half of locations each
n_vars <- 2 # number of variables per location
mix_p <- c("gauss_cop" = 0.5, "t_cop" = 0.5) # mixture percentages
# Normal distribution parameters (using standard normal with 0 correlation)
# mu <- rep(0, n_vars)
# sigma <- diag(1, n_vars)
# GPD parameters
scale_gpd <- 1
shape_gpd <- -0.05

# correlation parameters for each cluster
# cor_gauss <- c(0.1, 0.1)
cor_gauss <- 0.1
cor_t <- c(0.9, 0.1)
# degrees of freedom for t-Copula
# df_t <- c(1, 5)
df_t <- c(3, 3)
# quantile for Vignotto 2021 method
prob <- 0.9

# Number of cores to use for parallel computation
n_cores <- detectCores() - 1 


#### Simulate data and test clustering ####

# simulate from Gaussian copula with GPD margins
# (This will form asymptotically independent data)
# TODO: How does using Gaussian copula differ from MVN distribution?
gauss_cop <- lapply(seq_len(n_locs), \(i){
  # cor <- cor_gauss[[1]]
  # if (i > floor(n_locs / 2)) {
  #   cor <- cor_gauss[[2]]
  # }
  
  # cop_norm <- normalCopula(cor, dim = n_vars, dispstr = "un")
  cop_norm <- normalCopula(cor_gauss, dim = 2, dispstr = "un")
  u <- rCopula(n, cop_norm)
  # return(qnorm(u, mean = mu, sd = sigma))
  evd::qgpd(
    p     = u,
    # loc   = max(gauss_cop[[i]]), 
    loc   = 0,
    scale = scale_gpd, 
    shape = shape_gpd
  )
})

# check correlation is correct
# vapply(gauss_cop, \(x) round(cor(x)[1, 2], 3), numeric(1))
par(mfrow = c(2, 2))
lapply(gauss_cop[c(1, 5, 7, 10)], plot)
par(mfrow = c(1, 1))

# simulate from t-Copula with GPD margins
t_cop <- lapply(seq_len(n_locs), \(i) {
  cor <- cor_t[[1]]
  df <- df_t[[1]]
  if (i > floor(n_locs / 2)) {
    cor <- cor_t[[2]]
    df <- df_t[[2]]
  }
  cop_t <- copula::tCopula(cor, dim = 2, df = df_t[[2]], dispstr = "un")
  u <- rCopula(n, cop_t)
  return(evd::qgpd(
    p     = u,
    # loc   = max(gauss_cop[[i]]), 
    loc   = 0,
    scale = scale_gpd, 
    shape = shape_gpd
  ))
})

vapply(t_cop, \(x) round(cor(x)[1, 2], 3), numeric(1))
par(mfrow = c(2, 2))
lapply(t_cop[c(1, 5, 7, 10)], plot)
par(mfrow = c(1, 1))

# mixture
# data_mix <- lapply(seq_len(n_locs), \(i){
#   mix_p[[1]] * gauss_cop[[i]] + mix_p[[2]] * t_cop[[i]]
# })
data_mix <- lapply(seq_len(n_locs), \(i) {
  x <- nrow(gauss_cop[[i]])
  y <- nrow(t_cop[[i]])
  # sample mix_p * nrow(gauss_cop) rows from gauss_cop (same for t_cop)
  rbind(
    gauss_cop[[i]][sample(seq_len(x), size = mix_p[[1]] * x), ],
    t_cop[[i]][sample(seq_len(y), size = mix_p[[2]] * y), ]
  )
})

vapply(data_mix, \(x) round(cor(x)[1, 2], 3), numeric(1))
par(mfrow = c(2, 2))
lapply(data_mix[c(1, 5, 7, 10)], plot)
par(mfrow = c(1, 1))

# test dissimilarity for one
emp_kl_div(
  data_mix[[1]], data_mix[[10]], prob = 0.95, plot = TRUE
)

# KL divergence between areas using Vignotto 2021 method
kl_mat <- proxy::dist(
  data_mix, method = emp_kl_div, print = FALSE, prob = 0.9
)

# clustering solution
pam_kl_clust <- pam(kl_mat, k = 2)
# evaluate quality
mclust::adjustedRandIndex(
  pam_kl_clust$clustering, 
  cluster_mem
)


#### Functionalise ####

# Function to generate copula data
sim_cop_dat <- \(
  n_locs = 12,     # number of locations
  n = 1e4,         # number of simulations
  cor_gauss,       # bulk correlation for two clusters from Gaussian copula
  # params_norm,     # normal marginal parameters (same for both)
  cor_t,           # extreme correlation for two clusters from t-copula
  df_t,            # degrees of freedom for t-copula
  params_gpd,      # GPD margin parameters
  mix_p            # mixture weights
) {
  # many arguments must be of length 2 for 2 clusters
  stopifnot(all(vapply(
    # list(cor_gauss, params_norm, cor_t, df_t, params_gpd, mix_p), 
    list(cor_gauss, cor_t, df_t, params_gpd, mix_p), 
    \(x) length(x) == 2, logical(1)
  )))
  stopifnot(sum(mix_p) == 1)
  
  # Simulate from Gaussian Copula with GPD margins
  gauss_cop <- lapply(seq_len(n_locs), \(i){
    # pull correlation specified for each cluster
    cor <- cor_gauss[[1]]
    if (i > floor(n_locs / 2)) {
      cor <- cor_gauss[[2]]
    }
    cop_norm <- normalCopula(cor, dim = 2, dispstr = "un")
    # cop_norm <- normalCopula(cor_mat, dim = 2, dispstr = "un")
    u <- rCopula(n, cop_norm)
    # return(qnorm(u, mean = mu, sd = sigma))
    evd::qgpd(
      p     = u,
      # loc   = max(gauss_cop[[i]]), 
      loc   = 0,
      scale = scale_gpd, 
      shape = shape_gpd
    )
  })

  # simulate from t-Copula with GPD margins
  t_cop <- lapply(seq_len(n_locs), \(i) {
    cor <- cor_t[[1]]
    df <- df_t[[1]]
    if (i > floor(n_locs / 2)) {
      cor <- cor_t[[2]]
      df <- df_t[[2]]
    }
    cop_t <- copula::tCopula(cor, dim = 2, df = df_t[[2]], dispstr = "un")
    u <- rCopula(n, cop_t)
    return(evd::qgpd(
      p     = u,
      # loc   = max(gauss_cop[[i]]), 
      loc   = 0,
      scale = scale_gpd, 
      shape = shape_gpd
    ))
  })
  
  # mixture
  # data_mix <- lapply(seq_len(n_locs), \(i){
  #   mix_p[[1]] * gauss_cop[[i]] + mix_p[[2]] * t_cop[[i]]
  # })
  data_mix <- lapply(seq_len(n_locs), \(i) {
    x <- nrow(gauss_cop[[i]])
    y <- nrow(t_cop[[i]])
    # sample mix_p * nrow(gauss_cop) rows from gauss_cop (same for t_cop)
    rbind(
      gauss_cop[[i]][sample(seq_len(x), size = mix_p[[1]] * x), ],
      t_cop[[i]][sample(seq_len(y), size = mix_p[[2]] * y), ]
    )
  })
 return(list(
   "gauss_cop" = gauss_cop, 
   "t_cop"     = t_cop,
   "data_mix"  = data_mix
 ))
}

# TODO: Separate data generation and KL divergence/clustering procedures
# Function to 
kl_sim_eval <- \(
  data_mix, # Mixture data from copulas
  kl_prob, # Extremal quantile
  ...
) {
  
  # prob must be valid probability
  stopifnot(0 <= kl_prob && kl_prob <= 1)
  
  # KL divergence between areas using Vignotto 2021 method
  kl_mat <- proxy::dist(
    data_mix, method = emp_kl_div, print = FALSE, prob = kl_prob, ...
  )
  
  # clustering solution
  pam_kl_clust <- pam(kl_mat, k = 2)
  # evaluate quality
  adj_rand <- mclust::adjustedRandIndex(
    pam_kl_clust$clustering, 
    cluster_mem
  )
  # print(paste0("adjusted Rand index" = adj_rand))
  
  # return object
  return(list(
    "pam"      = pam_kl_clust,
    "adj_rand" = adj_rand
  ))
}

# test 
set.seed(seed_number)
data <- sim_cop_dat(
  n_locs = 12, 
  n = 1e3, 
  # cor_gauss = c(0.5, 0.5), 
  cor_gauss = c(0.1, 0.1), 
  cor_t = c(0.9, 0.1),
  # df_t = c(1, 5),
  df_t = c(3, 3),
  params_gpd = c(1, -0.05),
  mix_p = c(0.5, 0.5)
)
data_mix <- data$data_mix
kl_sim_eval(data_mix, kl_prob = 0.92)


#### Grid search ####

# create grid 
# TODO: Reduce number of parameters (Currently 128k, can't possibly run so many times)
grid <- tidyr::crossing(
  # Don't need all cor_vals for both (i.e. (0.3, 0.8) is the same as (0.8, 0.3)
  cor_gauss1 = seq(0, 1, by = 0.1),
  cor_gauss2 = seq(0, 1, by = 0.1),
  # cor_gauss1 = 0.1,
  # cor_gauss2 = 0.1,
  # cor_gauss = seq(0, 1, by = 0.1),
  # t-copula correlation
  cor_t1    = seq(0.1, 0.9, by = 0.2),
  cor_t2    = seq(0, 1, by = 0.2),
  # Degrees of freedom for t-copula
  # df_t1     = c(1, 5, 10),
  # df_t2     = c(1, 5, 10),
  # df_t1   = 1,
  # df_t2   = 1,
  df_t1   = 3,
  df_t2   = 3,
  # mixture percentages (must sum to 1)
  mix_p1  = 0.5,
  mix_p2  = 0.5,
  # extremal quantiles
  # kl_prob = c(0.9, 0.95, 0.99)
  kl_prob = 0.9
) %>% 
  filter(
    # use same correlation in both clusters for Gaussian copula
    cor_gauss1 == cor_gauss2,
    # mixture weights must sum to 1
    mix_p1 + mix_p2 == 1
  )

# run kl_sim_eval for each row in grid
# results_grid <- bind_rows(lapply(seq_len(nrow(grid)), \(i) {
# TODO: Parallelise
set.seed(seed_number)
results_grid <- bind_rows(mclapply(seq_len(nrow(grid)), \(i) {
  
  print(paste0("Progress: ", round(i / nrow(grid), 3) * 100, "%"))
  
  row <- grid[i, , drop = FALSE]
  data_mix <- with(row, sim_cop_dat(
    n = 1e3,
    cor_gauss  = c(cor_gauss1, cor_gauss2),
    cor_t      = c(cor_t1, cor_t2),
    df_t       = c(df_t1, df_t2),
    params_gpd = c(scale_gpd, shape_gpd),
    mix_p      = c(0.5, 0.5)
  ))$data_mix
  
  kl_clust <- tryCatch({
    kl_sim_eval(data_mix, kl_prob = row$kl_prob)
  # if an error is produced, return a dummy list
  }, error = function(cond) {
    return(list(
      "adj_rand" = NA,
      "pam" = list("clustering" = NA)
    ))
  })
   
  return(cbind(
    row,
    data.frame(
      "adj_rand" = kl_clust$adj_rand, 
      "membership" = paste0(kl_clust$pam$clustering, collapse = "-")
    )
  ))
# }))
}, mc.cores = n_cores))

# save
saveRDS(results_grid, file = "data/vignotto_grid_search_res.RDS")


#### Plot ####

ggplot(results_grid, aes(x = cor_t1, y = adj_rand, colour = factor(cor_t2))) +
  geom_point() +
  geom_line() +
  facet_wrap(~ cor_gauss1, ncol = 3) +
  theme

ggplot(results_grid, aes(x = cor_gauss1, y = adj_rand)) +
  geom_point() +
  geom_line() +
  # facet_wrap(~ cor_gauss1, ncol = 3) +
  facet_grid(cor_t1 ~ cor_t2) + 
  theme
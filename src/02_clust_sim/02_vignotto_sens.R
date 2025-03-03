#### Simulation sensitivity analysis of Vignotto 2021 ####

#### libs ####

library(copula)
library(evc)
library(evd)
library(dplyr)
library(tidyr)
library(lcmix)
library(ggplot2)
library(parallel)

# source functions
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
cor_gauss <- 0.5
cor_t <- c(0.8, 0.2)
# degrees of freedom for t-Copula
# df_t <- c(1, 5)
df_t <- c(5, 5)
# quantile for Vignotto 2021 method
prob <- 0.9

# Number of cores to use for parallel computation
n_cores <- detectCores() - 1 


#### Simulate data and test clustering ####

# test data simulation function
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
  # cor_t1    = seq(0, 1, by = 0.1),
  # cor_t2    = seq(0, 1, by = 0.1),
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
# TODO: Functionalise as used in other scripts
n_times <- 100
results_vec <- lri_vec <- lri_mean_vec <- vector(length = n_times)
set.seed(seed_number)
# results_grid <- bind_rows(lapply(seq_len(nrow(grid)), \(i) {
# results_grid <- bind_rows(mclapply(1:3, \(i) { # test
results_grid <- bind_rows(mclapply(seq_len(nrow(grid)), \(i) {

  print(paste0("Progress: ", round(i / nrow(grid), 3) * 100, "%"))
  system(sprintf(
    'echo "\n%s\n"', 
    paste0(round(i / nrow(grid), 3) * 100, "% completed", collapse = "")
  ))
  
  row <- grid[i, , drop = FALSE]
  for (j in seq_len(n_times)) {
    data_mix <- with(row, sim_cop_dat(
      n          = 1e3,
      cor_gauss  = c(cor_gauss1, cor_gauss2),
      cor_t      = c(cor_t1, cor_t2),
      df_t       = c(df_t1, df_t2),
      params_gpd = c(scale_gpd, shape_gpd),
      mix_p      = c(0.5, 0.5)
    ))$data_mix
  
    kl_clust <- tryCatch({
      kl_sim_eval(
        data_mix, kl_prob = row$kl_prob, k = 2, cluster_mem = cluster_mem
      )
    # if an error is produced, return a dummy list
    }, error = function(cond) {
      return(list("adj_rand" = NA))
    })
    results_vec[[j]] <- kl_clust$adj_rand
    
    # also calculate local rand index
    lri <- local_rand_index(kl_clust$pam$clustering, cluster_mem)
    # concatonate to string to store in vector
    lri_vec[[j]] <- paste0(lri, collapse = "-")
    # also average by cluster
    lri_mean_vec[[j]] <- paste0(vapply(unique(cluster_mem), \(x) {
      mean(lri[cluster_mem == x])
    }, numeric(1)), collapse = "-")
  }
  
  # return(cbind(row, "adj_rand" = results_vec))
  return(cbind(
    row, 
    "adj_rand"        = results_vec, 
    "local_rand"      = lri_vec,
    "mean_local_rand" = lir_mean_vec 
  ))
# })
}, mc.cores = n_cores))

# calculate means and add to results
results_grid_mean <- results_grid %>% 
  group_by(across(1:last_col(1))) %>% 
  summarise(mean_rand = mean(adj_rand, na.rm = TRUE), .groups = "drop")

results_grid <- results_grid %>% 
  left_join(results_grid_mean)

# save
saveRDS(results_grid, file = "data/vignotto_grid_search_res.RDS")
results_grid <- readRDS("data/vignotto_grid_search_res.RDS")


#### Plot ####

p1 <- results_grid %>% 
  ggplot() + 
  geom_point(
    aes(x = cor_gauss1, y = adj_rand), 
    colour = "black", alpha = 0.05, size = 0.8
  ) + 
  geom_point(
    aes(x = cor_gauss1, y = mean_rand), 
    colour = "#BC3C29FF", 
    size = 2
  ) + 
  geom_line(
    aes(x = cor_gauss1, y = mean_rand), 
    colour = "#BC3C29FF", linewidth = 1.2
  ) + 
  facet_grid(cor_t1 ~ cor_t2) + 
  labs(y = "Adjusted Rand Index", x = "Gaussian corr") + 
  theme + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

ggsave(plot = p1, "plots/01a_vignotto_sensitivity.png")

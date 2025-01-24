#### Simulation sensitivity analysis of Jensen-Shannon divergence method ####

# More realistic example:
# - 60 locations (similar to Ireland data)
# - 3 clusters (for 20 locations each)
# - Add perturbation to copula correlation parameters

#### libs ####

library(dplyr, quietly = TRUE)
library(evc)
# devtools::load_all("../evc")
library(tidyr)
library(ggplot2)
library(evd)
library(evgam)
# devtools::load_all("texmex")
library(cluster)
library(lcmix)
library(copula)
library(parallel)

# source functions
source("src/functions.R")


#### Metadata ####

seed_number <- 123
n <- 1e3 # number of samples to take
n_locs <- 60 # number of "locations"
# three clusters 
n_clust <- 3
# # for a third of locations each
# cluster_mem <- sort(rep(seq_len(n_clust), n_locs / n_clust)) 
# 1st cluster has 30, second 20, last 10
cluster_mem <- c(rep(1, 30), rep(2, 20), rep(3, 10))
n_vars <- 2 # number of variables per location
perturb_cor <- TRUE # perturb correlation within cluster for each location
mix_p <- c("gauss_cop" = 0.5, "t_cop" = 0.5) # mixture percentages
cor_gauss <- 0.5 # constant Guassian correlation throughout
# GPD parameters
scale_gpd <- 1
shape_gpd <- -0.05
marg_prob <- 0.9
kl_prob <- 0.9
conf_level <- 0.95 # confidence level for CIs in plot

# Number of cores to use for parallel computation
n_cores <- detectCores() - 1 


#### Grid search ####

# create grid 
set.seed(seed_number)
grid <- tidyr::crossing(
  # Don't need all cor_vals for both (i.e. (0.3, 0.8) is the same as (0.8, 0.3)
  # cor_gauss = seq(0, 1, by = 0.1),
  # keep Gaussian correlation constant
  cor_gauss = cor_gauss,
  # t-copula correlation
  cor_t1    = seq(0, 1, by = 0.3),
  cor_t2    = seq(0.1, 1, by = 0.3),
  cor_t3    = seq(0.2, 1, by = 0.3),
  # Degrees of freedom for t-copula
  df_t      = 3,
  # mixture percentages (must sum to 1)
  mix_p1    = 0.5,
  mix_p2    = 0.5,
  # extremal quantiles
  # kl_prob = c(0.9, 0.95, 0.99)
  kl_prob   = kl_prob
) %>% 
  # mixture weights must sum to 1
  filter(mix_p1 + mix_p2 == 1) %>% 
  # round to 1 decimal place (weird floating point error happening here)
  mutate(across(contains("cor"), \(x) round(x, 1))) %>% 
  identity()

# run kl_sim_eval for each row in grid
# n_times <- 50
n_times <- 100
results_vec <- vector(length = n_times)
set.seed(seed_number)
# # test: seq w/ biggest cor diffs between clusters; should be easy to cluster
# test_seq <- with(
#   grid, 
#   which(cor_t1 == 0 & cor_t2 == 0.4 & cor_t3 == 0.8)
# )
# results_grid <- bind_rows(mclapply(test_seq, \(i) {
results_grid <- bind_rows(mclapply(seq_len(nrow(grid)), \(i) {
   
  system(sprintf(
    'echo "\n%s\n"', 
    paste0(round(i / nrow(grid), 3) * 100, "% completed", collapse = "")
  ))
  
  row <- grid[i, , drop = FALSE]
  for (j in seq_len(n_times)) {
    # generate simulation data for given parameter set
    # source("src/functions.R")
    # debugonce(sim_cop_dat)
    data_mix <- with(row, sim_cop_dat(
      n           = 1e3,
      n_locs      = n_locs,
      n_clust     = n_clust,
      cluster_mem = cluster_mem,
      cor_gauss   = rep(cor_gauss, n_clust),
      cor_t       = c(cor_t1, cor_t2, cor_t3),
      df_t        = rep(df_t, n_clust),
      params_gpd  = c(scale_gpd, shape_gpd),
      mix_p       = c(0.5, 0.5),
      perturb_cor = perturb_cor # perturb cor for different locations
    ))$data_mix
    
    clust_res <- tryCatch({
      # if an error is produced, return a dummy list
      dependence <- fit_ce(
        data_mix, 
        marg_prob   = marg_prob,
        cond_prob   = row$kl_prob, 
        fit_no_keef = TRUE
      )
      js_clust(dependence, k = n_clust, cluster_mem = cluster_mem)
    }, error = function(cond) {
      return(list("adj_rand" = NA))
    })
    results_vec[[j]] <- clust_res$adj_rand
  }
  return(cbind(row, "adj_rand" = results_vec))
}, mc.cores = n_cores))

# remove NAs, give warning
nas <- is.na(results_grid$adj_rand)
if (sum(nas) > 0) {
  message(paste0("Warning: ", sum(nas), " NAs in results, removing"))
  results_grid <- results_grid[!nas, ]
}

# Summarise results to add confidence levels to plot
results_grid_sum <- summarise_sens_res(results_grid, conf_level = conf_level)
# tally rand index occurrences for given parameter set
results_grid_tally <- results_grid %>% 
  # round so you don't have lots of similar values
  # mutate(adj_rand_rnd = round(adj_rand, 2)) %>% 
  # group_by(across(!contains("_rand")), adj_rand_rnd) %>% 
  group_by(across(!contains("_rand")), adj_rand) %>% 
  tally(name = "n_rand") %>% 
  ungroup()
# add both to results
results_grid <- results_grid %>% 
  # for repeated joining
  dplyr::select(-(ends_with("_rand") & !matches("adj_rand"))) %>% 
  left_join(results_grid_sum) %>% 
  left_join(results_grid_tally)

# save 
# TODO: Change to csv file!
saveRDS(
  results_grid, 
  # file = paste0("data/js_grid_search_res_n_loc_", n_locs, "_dqu_", kl_prob, "_marg_0.9.RDS")
  file = paste0(
    "data/js_grid_search_res_n_loc_", n_locs, 
    "_dqu_", kl_prob, 
    "_gauss_cor_", cor_gauss,
    ".RDS"
  )
)

# # load data and redo results_grid_sum and results_grid_tally
# # May be floating point error in some parameter values due to loading from RDS!
results_grid <- readRDS(
  # paste0("data/js_grid_search_res_n_loc_", n_locs, "_dqu_", kl_prob, "_marg_0.9.RDS")
  paste0(
    "data/js_grid_search_res_n_loc_", n_locs, 
    "_dqu_", kl_prob, 
    "_gauss_cor_", cor_gauss,
    ".RDS"
  )
)
results_grid_sum <- summarise_sens_res(results_grid, conf_level = conf_level)
results_grid_tally <- results_grid %>%
  group_by(across(!contains("_rand")), adj_rand) %>%
  tally(name = "n_rand")


#### Plot ####
# remove unlikely scenario of having perfect t-copula correlation
results_grid_plt <- filter(results_grid, cor_t1 < 1, cor_t2 < 1)
res_grid_sum_plt <- filter(results_grid_sum, cor_t1 < 1, cor_t2 < 1)
# don't want too many dots, round and sum up counts
# res_grid_tally_plt <- results_grid_tally %>% 
#   filter(cor_t1 < 1 & cor_t2 < 1) %>% 
#   mutate(adj_rand = round(adj_rand, 1)) %>% 
#   group_by(across(!contains("_rand")), adj_rand) %>% 
#   summarise(n_rand = sum(n_rand), .groups = "drop")

p1 <- results_grid_plt %>% 
  ggplot() + 
  geom_point(
    aes(x = cor_t3, y = adj_rand), 
    colour = "black", alpha = 0.05, size = 1
  ) +   
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) + 
  # Confidence intervals around mean line
  geom_ribbon(
    data = res_grid_sum_plt, 
    aes(x = cor_t3, ymin = lower_rand, ymax = upper_rand), 
    fill = "#C11432",
    alpha = 0.3,
    size = 2
  ) + 
  geom_line(
    data = res_grid_sum_plt, 
    aes(x = cor_t3, y = mean_rand), 
    colour = "#C11432",
    linewidth = 1
  ) + 
  facet_grid(cor_t1 ~ cor_t2) +  
  labs(y = "Adjusted Rand Index", x = "rho_t_3") + 
  evc_theme() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + 
  ggsci::scale_fill_nejm() + 
  ggsci::scale_colour_nejm()
p1


smooth_plt <- \(data, col) {
  data %>% 
    ggplot() + 
    # plot points, and smooth line through them
    geom_point(
      aes(x = .data[[col]], y = adj_rand), 
      colour = "black", alpha = 0.05, size = 1
    ) +   
    geom_smooth(
      aes(x = .data[[col]], y = adj_rand), 
      colour    = "#C11432",
      fill      = "#C11432",
      formula   = y ~ splines::bs(x = x, knots = 3),
      se        = TRUE,
      linewidth = 1
    ) + 
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) + 
    facet_grid(cor_t1 ~ cor_t2) +  
    # labs(y = "Adjusted Rand Index", x = "Gaussian corr") + 
    labs(y = "Adjusted Rand Index", x = "") + 
    evc_theme() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + 
    ggsci::scale_fill_nejm() + 
    ggsci::scale_colour_nejm()
}

p1 <- smooth_plt(results_grid_plt, "cor_gauss")
p2 <- smooth_plt(results_grid_plt, "cor_t3")

data %>% 
  ggplot() + 
  # plot points, and smooth line through them
  geom_point(
    aes(x = .data[[col]], y = adj_rand), 
    colour = "black", alpha = 0.05, size = 1
  ) +   
  geom_smooth(
    aes(x = .data[[col]], y = adj_rand), 
    colour = "#C11432",
    fill = "#C11432",
    formula = y ~ splines::bs(x = x, knots = 3),
    se = TRUE,
    linewidth = 1
  ) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) + 
  facet_grid(cor_t1 ~ cor_t2) +  
  # labs(y = "Adjusted Rand Index", x = "Gaussian corr") + 
  labs(y = "Adjusted Rand Index", x = "") + 
  evc_theme() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + 
  ggsci::scale_fill_nejm() + 
  ggsci::scale_colour_nejm()
  

# ternary plot
res_grid_sum_plt %>%
  dplyr::select(contains("cor_t"), mean_rand) %>% 
  ggtern(
    aes(x = cor_t1, y = cor_t2, z = cor_t3) 
  ) + 
  stat_density_tern(geom="polygon",aes(fill=..level..,alpha=..level..)) + 
  geom_point(aes(colour = mean_rand), size = 2) + 
  viridis::scale_colour_viridis() + 
  tern_limits(T = 1.05, L = 1.05, R = 1.05) + 
  guides(color = "none", fill = "none", alpha = "none")


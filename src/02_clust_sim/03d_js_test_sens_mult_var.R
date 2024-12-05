#### Apply Jensen-Shannon Divergence method to simulated data ####

#### Libs ####

library(dplyr, quietly = TRUE)
# library(evc)
devtools::load_all("../evc")
library(tidyr)
library(ggplot2)
library(cluster)
library(evgam)
library(parallel)
devtools::load_all("texmex")

source("src/functions.R")


#### Metadata ####

# metadata for original testing of method for > 2 vars
seed_number <- 123
# number of "locations"
n_locs <- 12 
# Multivariate (> 2-dimensions)
# vars <- c("rain", "wind_speed")
vars <- c("rain", "wind_speed", "temperature") # TODO: Required?
n_vars <- length(vars)
prob <- 0.9 # TODO: May have to lower for more than two variables?
cluster_mem <- sort(rep(c(1, 2), 6)) # known cluster membership
n <- 1e4 # number of samples to take
mix_p <- c("gauss_cop" = 0.5, "t_cop" = 0.5) # mixture percentages
# Copula correlation parameters
# TODO: May have to not hardcode??
cor_gauss <- c(0.5, 0.5)
cor_t = c(0.9, 0.1)
# t-copula df
df_t <- c(3, 3)
# GPD parameters
scale_gpd <- 1
shape_gpd <- -0.05
# confidence level for summary and plots
conf_level <- 0.95

# Number of cores to use for parallel computation
n_cores <- detectCores() - 1 


#### Simulate Data ####

set.seed(seed_number)
data <- sim_cop_dat(
  n_locs = n_locs, 
  n_vars = n_vars,
  n = n, 
  cor_gauss = cor_gauss,
  cor_t = cor_t,
  df_t = df_t,
  params_gpd = c(scale_gpd, shape_gpd),
  mix_p = c(0.5, 0.5)
)
data_mix <- data$data_mix


#### Calculate Conditional Extremes parameters ####

# TODO: Extend to work with list of locations with stacked vector to split
dependence <- fit_ce(
  data        = data_mix, 
  vars        = vars,
  marg_prob   = prob, 
  cond_prob   = prob, 
  fit_no_keef = TRUE
)

# check that all dependence models have run successfully
sapply(dependence, \(x) lapply(x, length))


#### Calculate divergence upon which to cluster ####

# Perform PAM and k-means clustering, as an exploratory analysis

# first, produce distance matrix and associated elbow plot
js_mat <- js_clust(dependence)$dist_mat # suggests k = 2, as desired

# silhouette plot, also suggests k = 2, works really well here
sil_boxplot(js_mat, k = 2:6)

# cluster and assess performance
js_clust(dependence, k = 2, dist_mat = js_mat, cluster_mem = cluster_mem)


#### Sensitivity Analysis ####

# create grid 
grid <- tidyr::crossing(
  cor_gauss1 = seq(0, 1, by = 0.1),
  cor_gauss2 = seq(0, 1, by = 0.1),
  # t-copula correlation
  cor_t1    = seq(0.1, 0.9, by = 0.2),
  cor_t2    = seq(0, 1, by = 0.2),
  # Degrees of freedom for t-copula
  df_t1   = 3,
  df_t2   = 3,
  # mixture percentages (must sum to 1)
  mix_p1  = 0.5,
  mix_p2  = 0.5,
  # extremal quantiles
  # kl_prob = c(0.9, 0.95, 0.99)
  kl_prob = prob 
) %>% 
  filter(
    # use same correlation in both clusters for Gaussian copula
    cor_gauss1 == cor_gauss2,
    # mixture weights must sum to 1
    mix_p1 + mix_p2 == 1
  )

# run kl_sim_eval for each row in grid
# n_times <- 50
n_times <- 100
results_vec <- vector(length = n_times)
set.seed(seed_number)
# results_grid <- bind_rows(lapply(seq_len(nrow(grid)), \(i) {
results_grid <- bind_rows(mclapply(seq_len(nrow(grid)), \(i) {
  
  print(paste0("Progress: ", round(i / nrow(grid), 3) * 100, "%"))
  system(sprintf(
    'echo "\n%s\n"', 
    paste0(round(i / nrow(grid), 3) * 100, "% completed", collapse = "")
  ))
  
  row <- grid[i, , drop = FALSE]
  for (j in seq_len(n_times)) {
    # generate simulation data for given parameter set
    data_mix <- with(row, sim_cop_dat(
      n_locs     = n_locs, 
      n_vars     = n_vars,
      n          = 1e3,
      cor_gauss  = c(cor_gauss1, cor_gauss2),
      cor_t      = c(cor_t1, cor_t2),
      df_t       = c(df_t1, df_t2),
      params_gpd = c(scale_gpd, shape_gpd),
      mix_p      = c(0.5, 0.5)
    ))$data_mix
    
    # if an error is produced, return a dummy list
    clust_res <- tryCatch({
      # fit CE model
      dependence <- fit_ce(
        data        = data_mix, 
        vars        = c("rain", "ws", "temp"),
        marg_prob   = row$kl_prob,
        cond_prob   = row$kl_prob
      )
      # cluster based on results
      js_clust(dependence, k = 2, cluster_mem = cluster_mem)
    }, error = function(cond) {
      return(list("adj_rand" = NA))
    })
    results_vec[[j]] <- clust_res$adj_rand
  }
  
  # return(cbind(
  #   row,
  #   data.frame(
  #     "adj_rand" = kl_clust$adj_rand, 
  #     "membership" = paste0(kl_clust$pam$clustering, collapse = "-")
  #   )
  # ))
  return(cbind(row, "adj_rand" = results_vec))
  # }))
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
  file = paste0("data/js_sens_res_nvar_", n_vars, "_dqu_", prob, ".RDS")
  # file = paste0("data/js_grid_search_res_dqu_", kl_prob, "_marg_0.9.RDS")
)
# # load data and redo results_grid_sum and results_grid_tally
# # May be floating point error in some parameter values due to loading from RDS!
# results_grid <- readRDS(
#   paste0("data/js_grid_search_res_dqu_", kl_prob, ".RDS")
#   # paste0("data/js_grid_search_res_dqu_", kl_prob, "_marg_0.9.RDS")
# )
# results_grid_sum <- summarise_sens_res(results_grid, conf_level = conf_level)
# results_grid_tally <- results_grid %>%
#   group_by(across(!contains("_rand")), adj_rand) %>%
#   tally(name = "n_rand")

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
# p1 <- res_grid_tally_plt %>% 
  # mutate(adj_rand = round(adj_rand, 2)) %>% 
  # TODO: Temp for fixing plots
  # filter(cor_t1 == 0.7 & cor_t2 == 0.4) %>% 
  ggplot() + 
  geom_point(
    aes(x = cor_gauss1, y = adj_rand), 
    colour = "black", alpha = 0.05, size = 1
  ) +   
  # TODO: Size background points by number of occurrences
  # TODO: Or maybe do alpha by size?
  # geom_boxplot(aes(x = factor(cor_gauss1), y = adj_rand)) + 
  # geom_point(
  #   aes(x = cor_gauss1, y = adj_rand, size = n_rand, stroke = 1.5), 
  #   fill = "white",
  #   colour = "black",
  #   alpha = 0.4,
  #   shape = 21
  # ) +
  # facet_grid(cor_t1 ~ cor_t2) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) + 
  # Confidence intervals around mean line
  geom_ribbon(
    data = res_grid_sum_plt, # %>% filter(cor_t1 == 0.3, cor_t2 == 0.6),
    aes(x = cor_gauss1, ymin = lower_rand, ymax = upper_rand), 
    fill = "#C11432",
    alpha = 0.3,
    size = 2
  ) + 
  geom_line(
    data = res_grid_sum_plt, # %>% filter(cor_t1 == 0.3, cor_t2 == 0.6),
    aes(x = cor_gauss1, y = mean_rand), 
    colour = "#C11432",
    linewidth = 1
  ) + 
  facet_grid(cor_t1 ~ cor_t2) +  
  labs(y = "Adjusted Rand Index", x = "Gaussian corr") + 
  evc_theme() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + 
  ggsci::scale_fill_nejm() + 
  ggsci::scale_colour_nejm()
p1

ggsave(
  plot = p1, 
  paste0("plots/01e_js_sens_3_var_dqu_", prob, ".png"), 
  width = 10, 
  height = 7
  # paste0("plots/01b_js_sensitivity_dqu_", kl_prob, "_marg_0.9.png")
)

#### Compare to two variables ####


results_grid2 <- readRDS("data/js_grid_search_res_dqu_0.9.RDS")

results_grid_join <- bind_rows(
  results_grid %>% 
    distinct() %>% 
    mutate(ind = "3 Variables") %>% 
    relocate(ind), 
  results_grid2 %>% 
    distinct() %>% 
    mutate(ind = "2 Variables") %>% 
    relocate(ind)
) %>% 
  # fix weird floating point problem from loading from .RDS file
  mutate(across(contains("cor"), \(x) round(x, 1)))

p2 <- results_grid_join %>% 
  ggplot() + 
  # points
  # geom_point(aes(x = cor_gauss1, y = mean_rand, colour = ind), size = 2) +   
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) + 
  # Confidence intervals around mean line
  geom_ribbon(
    aes(x = cor_gauss1, ymin = lower_rand, ymax = upper_rand, fill = ind), 
    # fill = "#C11432",
    alpha = 0.3,
    size = 2, 
  ) + 
  geom_line(
    aes(x = cor_gauss1, y = mean_rand, colour = ind), 
    # colour = "#C11432",
    linewidth = 1,
    show.legend = FALSE
  ) + 
  facet_grid(cor_t1 ~ cor_t2) +  
  labs(
    y = "Adjusted Rand Index", 
    x = "Gaussian corr", 
    fill = "# Variables per location"
  ) + 
  evc_theme() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + 
  guides(fill = guide_legend(override.aes = list(alpha = 1))) + 
  ggsci::scale_fill_nejm() + 
  ggsci::scale_colour_nejm()
p2

ggsave(
  plot = p2, 
  paste0("plots/01f_js_sens_3_var_compare_dqu_", prob, ".png"), 
  width = 10, 
  height = 7
  # paste0("plots/01b_js_sensitivity_dqu_", kl_prob, "_marg_0.9.png")
)

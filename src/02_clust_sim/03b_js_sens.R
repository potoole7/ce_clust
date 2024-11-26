#### Simulation sensitivity analysis of Jensen-Shannon divergence method ####

#### libs ####

library(dplyr, quietly = TRUE)
library(evc)
library(tidyr)
library(ggplot2)
library(evd)
library(evgam)
devtools::load_all("texmex")
library(cluster)
library(lcmix)
library(copula)
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
# GPD parameters
scale_gpd <- 1
shape_gpd <- -0.05
marg_prob <- 0.9
kl_prob <- 0.9

# Number of cores to use for parallel computation
n_cores <- detectCores() - 1 

#### Grid search ####

# create grid 
grid <- tidyr::crossing(
  # Don't need all cor_vals for both (i.e. (0.3, 0.8) is the same as (0.8, 0.3)
  cor_gauss1 = seq(0, 1, by = 0.1),
  cor_gauss2 = seq(0, 1, by = 0.1),
  # cor_gauss1 = 0.1,
  # cor_gauss2 = 0.1,
  # t-copula correlation
  cor_t1    = seq(0.1, 0.9, by = 0.2),
  cor_t2    = seq(0, 1, by = 0.2),
  # Degrees of freedom for t-copula
  # df_t1     = c(1, 5, 10),
  # df_t2     = c(1, 5, 10),
  df_t1   = 3,
  df_t2   = 3,
  # mixture percentages (must sum to 1)
  mix_p1  = 0.5,
  mix_p2  = 0.5,
  # extremal quantiles
  # kl_prob = c(0.9, 0.95, 0.99)
  kl_prob = kl_prob
) %>% 
  filter(
    # use same correlation in both clusters for Gaussian copula
    cor_gauss1 == cor_gauss2,
    # mixture weights must sum to 1
    mix_p1 + mix_p2 == 1
  )

# run kl_sim_eval for each row in grid
n_times <- 100
set.seed(seed_number)
# results_grid <- bind_rows(lapply(seq_len(nrow(grid)), \(i) {
results_grid <- bind_rows(mclapply(seq_len(nrow(grid)), \(i) {
  
  print(paste0("Progress: ", round(i / nrow(grid), 3) * 100, "%"))
  system(sprintf(
    'echo "\n%s\n"', 
    paste0(round(i / nrow(grid), 3) * 100, "% completed", collapse = "")
  ))
  
  row <- grid[i, , drop = FALSE]
  results_vec <- vector(length = n_times)
  for (j in seq_len(n_times)) {
    # generate simulation data for given parameter set
    data_mix <- with(row, sim_cop_dat(
      n          = 1e3,
      cor_gauss  = c(cor_gauss1, cor_gauss2),
      cor_t      = c(cor_t1, cor_t2),
      df_t       = c(df_t1, df_t2),
      params_gpd = c(scale_gpd, shape_gpd),
      mix_p      = c(0.5, 0.5)
    ))$data_mix
    
    js_clust <- tryCatch({
      # if an error is produced, return a dummy list
      dependence <- fit_ce(
        data_mix, 
        marg_prob   = marg_prob,
        cond_prob   = row$kl_prob,
        split_data  = FALSE
      )
      js_clust(dependence, k = 2, cluster_mem = cluster_mem)
    }, error = function(cond) {
      return(list("adj_rand" = NA))
    })
    results_vec[[j]] <- js_clust$adj_rand
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

# calculate means and add to results
results_grid_mean <- results_grid %>% 
  group_by(across(1:last_col(1))) %>% 
  summarise(mean_rand = mean(adj_rand, na.rm = TRUE), .groups = "drop")

results_grid <- results_grid %>% 
  left_join(results_grid_mean)

# save 
# TODO: Change to csv file!
saveRDS(
  results_grid, 
  file = paste0("data/js_grid_search_res_dqu_", kl_prob, ".RDS")
  # file = paste0("data/js_grid_search_res_dqu_", kl_prob, "_marg_0.9.RDS")
)
results_grid <- readRDS(
  paste0("data/js_grid_search_res_dqu_", kl_prob, ".RDS")
  # paste0("data/js_grid_search_res_dqu_", kl_prob, "_marg_0.9.RDS")
)


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
p1

# ggsave(plot = p1, "plots/01b_js_sensitivity.png")
ggsave(
  plot = p1, 
  paste0("plots/01b_js_sensitivity_dqu_", kl_prob, ".png")
  # paste0("plots/01b_js_sensitivity_dqu_", kl_prob, "_marg_0.9.png")
)


#### Compare to Vignotto ####

results_grid_vig <- readRDS("data/vignotto_grid_search_res.RDS")

preprocess_fun <- \(x, name) {
 ret <- x %>% 
  dplyr::select(-adj_rand) %>% 
  distinct()
 names(ret)[names(ret) == "mean_rand"] <- name
 return(ret)
}

results_grid_join <- preprocess_fun(results_grid, "rand_js") %>% 
  # kl probs may not match
  dplyr::select(-kl_prob) %>% 
  left_join(preprocess_fun(results_grid_vig, "rand_vig")) %>% 
  pivot_longer(contains("rand_"))

p2 <- results_grid_join %>% 
  mutate(name = ifelse(name == "rand_js", "CE", "Vignotto")) %>% 
  ggplot(aes(x = cor_gauss1, y = value, colour = name)) + 
  geom_point(size = 2) + 
  geom_line(linewidth = 1.2) + 
  facet_grid(cor_t1 ~ cor_t2) + 
  labs(y = "Adjusted Rand Index", x = "Gaussian corr") + 
  theme + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + 
  ggsci::scale_colour_nejm()
p2

ggsave(
  plot = p2, 
  paste0("plots/01c_sensitivity_dqu_", kl_prob, ".png"), 
  # paste0("plots/01c_sensitivity_dqu_", kl_prob, "_marg_0.9.png")
  width = 10, 
  height = 7
)


#### Compare sensitivity of Jensen-Shannon method for different quantiles ####

# pull relative files
# files <- list.files("data", pattern = "_marg_", full.names = TRUE)
# stopifnot(length(files) == 6)
# 
# # pull marginal probabilities for each (not labelled in data)
# marg_probs <- as.numeric(vapply(stringr::str_split(files, "_marg_"), \(x) {
#   stringr::str_remove(x[[2]], ".RDS")
# }, character(1)))
# 
# # pull data, labelling marginal probabilities
# data_js <- bind_rows(lapply(seq_along(files), \(i){
#   mutate(readRDS(files[[i]]), "marg_prob" = marg_probs[i])
# })) %>% 
#   mutate(indicator = paste0("dqu = ", kl_prob, ", mqu = ", marg_prob))
# 
# # plot
# data_js %>% 
#   dplyr::select(-adj_rand) %>% 
#   distinct() %>% 
#   filter(kl_prob != 0.99) %>% 
#   # ggplot(aes(x = cor_gauss1, y = value, colour = indicator)) + 
#   ggplot(aes(x = cor_gauss1, y = mean_rand, colour = indicator)) + 
#   geom_point(size = 2) + 
#   geom_line(linewidth = 1.2) + 
#   facet_grid(cor_t1 ~ cor_t2) + 
#   labs(y = "Adjusted Rand Index", x = "Gaussian corr") + 
#   theme + 
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + 
#   ggsci::scale_colour_nejm()
#   
# tabulate which is best in every case (i.e. every box)
# data_js %>% 
#   dplyr::select(-adj_rand) %>% 
#   distinct() %>% 
#   group_by(cor_t1, cor_t2) %>% 
#   filter(mean_rand == max(mean_rand, na.rm = FALSE)) %>% 
#   ungroup() %>% 
#   dplyr::count(cor_t1, cor_t2, indicator)

#### Simulation sensitivity analysis of Jensen-Shannon silhouette ####

# Checking whether silhouette plot chooses correct number of clusters

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
conf_level <- 0.95 # confidence level for CIs in plot

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
n_times <- 50
results_vec <- vector(length = n_times)
set.seed(seed_number)
results_grid <- bind_rows(mclapply(seq_len(nrow(grid)), \(i) {
# results_grid <- bind_rows(mclapply(c(1, 50, 100, 150, 200, 250, 300), \(i) {
  
  print(paste0("Progress: ", round(i / nrow(grid), 3) * 100, "%"))
  system(sprintf(
    'echo "\n%s\n"', 
    paste0(round(i / nrow(grid), 3) * 100, "% completed", collapse = "")
  ))
  
  row <- grid[i, , drop = FALSE]
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
    
    sil_res <- tryCatch({
      # Fit CE model, if an error is produced, return a dummy list
      dependence <- fit_ce(
        data_mix, 
        marg_prob   = marg_prob,
        cond_prob   = row$kl_prob, 
        fit_no_keef = TRUE
      )
      # find JS distance matrix
      dist <- js_clust(dependence)$dist_mat
      # Calculate silhouette coefficient for 2-10 clusters
      sil <- sil_boxplot(dist, k = 2:10, show_plt = FALSE)$sil
      # Find best number of clusters
      sil %>% 
        group_by(k) %>% 
        summarise(average = mean(sil_width, na.rm = TRUE), .groups = "drop") %>% 
        filter(average == max(average)) %>% 
        pull(k) %>% 
        return()
    }, error = function(cond) {
      return(NA)
    })
    results_vec[[j]] <- sil_res
  }
  
  # return(cbind(
  #   row,
  #   data.frame(
  #     "adj_rand" = kl_clust$adj_rand, 
  #     "membership" = paste0(kl_clust$pam$clustering, collapse = "-")
  #   )
  # ))
  return(cbind(row, "sil_k" = results_vec))
# }))
}, mc.cores = n_cores))

# remove NAs, give warning
nas <- is.na(results_grid$k)
if (sum(nas) > 0) {
  message(paste0("Warning: ", sum(nas), " NAs in results, removing"))
  results_grid <- results_grid[!nas, ]
}

# Save object
# saveRDS(results_grid, file = paste0("data/js_sil_best_k.RDS"))
results_grid <- readRDS(paste0("data/js_sil_best_k.RDS"))


#### Plot ####

# plot occurrence of values of k for different correlation values
# TODO: Add red tp points at 2
p1 <- results_grid %>% 
  group_by(cor_gauss1, cor_t1, cor_t2) %>% 
  dplyr::count(sil_k) %>% 
  mutate(
    col = case_when(
      n == max(n) ~ "Most frequent",
      sil_k == 2  ~  "True value (not found)",
      TRUE        ~  "Else"
    ), 
    col = factor(
      col, levels = c("Most frequent", "True value (not found)", "Else")
    )
  ) %>% 
  ggplot() + 
  geom_point(aes(x = cor_gauss1, y = sil_k, size = n, colour = col)) + 
  scale_y_continuous(
    limits = c(2, 9), 
    breaks = seq(2, 8, by = 2), 
    expand = c(0.15, 0)
  ) +
  scale_size_continuous(range = c(1.5, 6.5)) + 
  facet_grid(cor_t1 ~ cor_t2) + 
  evc_theme() + 
  labs(y = "Adjusted Rand Index", x = "Gaussian corr", colour = "") + 
  scale_colour_manual(values = c(ggsci::pal_nejm()(2), "black")) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
p1

ggsave(
  plot = p1, 
  paste0("plots/04_js_sil_best_k.png"), 
  width = 15, height = 9
  # paste0("plots/01b_js_sensitivity_dqu_", kl_prob, "_marg_0.9.png")
)

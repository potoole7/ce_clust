#### Simulation sensitivity analysis of Jensen-Shannon divergence method ####

# NOTE Unlikely to make final (or draft) paper!

# Testing for different # data points and multiples of maximum threshold values

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
conf_level <- 0.95 # confidence level for CIs

# Number of cores to use for parallel computation
n_cores <- detectCores() - 2

#### Grid search ####

# create grid
grid <- tidyr::crossing(
  cor_gauss1 = seq(0, 1, by = 0.1),
  cor_gauss2 = seq(0, 1, by = 0.1),
  # t-copula correlation
  # cor_t1    = seq(0.1, 0.9, by = 0.2),
  # cor_t2    = seq(0, 1, by = 0.2),
  # TODO: More t-copula correlation values???
  cor_t1 = 0.4,
  cor_t2 = 0.7,
  # Degrees of freedom for t-copula
  df_t1 = 3,
  df_t2 = 3,
  # mixture percentages (must sum to 1)
  mix_p1 = 0.5,
  mix_p2 = 0.5,
  # extremal quantiles
  kl_prob = kl_prob,
  # data point parameterisation
  dat_max_mult = 2:5,
  # n_dat = seq(10, 50, by = 10)
  # n_dat = seq(1, 9, by = 2)
  n_dat = c(1:5, 10, 20, 30)
  # n_dat <- c(seq(1, 9, by = 1), seq(10, 50, by = 10))
) %>%
  filter(
    # use same correlation in both clusters for Gaussian copula
    cor_gauss1 == cor_gauss2,
    # mixture weights must sum to 1
    mix_p1 + mix_p2 == 1
  )

n_times <- 100
results_vec <- vector(length = n_times)
set.seed(seed_number)
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
      n          = 1e3,
      cor_gauss  = c(cor_gauss1, cor_gauss2),
      cor_t      = c(cor_t1, cor_t2),
      df_t       = c(df_t1, df_t2),
      # params_gpd = c(scale_gpd, shape_gpd),
      mix_p      = c(0.5, 0.5),
      qfun       = evd::qgpd,
      qargs      = c("scale" = scale_gpd, "shape" = shape_gpd)
    ))$data_mix

    clust_res <- tryCatch(
      {
        # if an error is produced, return a dummy list
        dependence <- fit_ce(
          data_mix,
          marg_prob   = marg_prob,
          cond_prob   = row$kl_prob,
          fit_no_keef = TRUE
        )
        js_clust(
          dependence,
          k            = 2,
          cluster_mem  = cluster_mem,
          dat_max_mult = row$dat_max_mult,
          n_dat        = row$n_dat
        )
      },
      error = function(cond) {
        return(list("adj_rand" = NA))
      }
    )
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
# saveRDS(results_grid, file = "data/js_grid_search_n_points.RDS")
results_grid <- readRDS(file = "data/js_grid_search_n_points.RDS")


# remove NAs, give warning
# TODO: Lots of NAs, investigate!!!!
nas <- is.na(results_grid$adj_rand)
if (sum(nas) > 0) {
  message(paste0("Warning: ", sum(nas), " NAs in results, removing"))
  results_grid <- results_grid[!nas, ]
}


#### Plot ####

# TODO: Average by dat_max_mult, average by n_dat, plot
p1 <- results_grid %>%
  group_by(cor_gauss1, dat_max_mult) %>%
  # summarise(mean_rand = mean(adj_rand, na.rm = TRUE), .groups = "drop") %>%
  summarise(
    across(
      ends_with("_rand") & !matches("adj_rand"), \(x) mean(x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  ggplot() +
  # geom_point(
  #   aes(x = cor_gauss1, y = mean_rand, colour = factor(dat_max_mult)),
  #   size = 2
  # ) +
  # geom_ribbon(
  #   aes(x = cor_gauss1, ymin = lower_rand, ymax = upper_rand, fill = factor(dat_max_mult)),
  #   alpha = 0.3,
  #   linewidth = 2
  # ) +
  geom_line(
    aes(x = cor_gauss1, y = mean_rand, colour = factor(dat_max_mult)),
    linewidth = 1.2,
    show.legend = FALSE
  ) +
  ggsci::scale_colour_nejm() +
  labs(
    y = "Adjusted Rand Index",
    x = "Gaussian corr",
    fill = "times threshold"
  ) +
  theme +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
p1

ggsave(
  plot = p1,
  paste0("plots/01c_js_sens_times_thresh_dqu_", kl_prob, ".png")
  # paste0("plots/01b_js_sensitivity_dqu_", kl_prob, "_marg_0.9.png")
)


p2 <- results_grid_sum %>%
  # filter(n_dat %in% 1:10) %>%
  group_by(cor_gauss1, n_dat) %>%
  summarise(
    across(
      ends_with("_rand") & !matches("adj_rand"), \(x) mean(x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  ggplot() +
  # geom_point(
  #   aes(x = cor_gauss1, y = mean_rand, colour = factor(n_dat)),
  #   size = 2
  # ) +
  # geom_ribbon(
  #   aes(x = cor_gauss1, ymin = lower_rand, ymax = upper_rand, fill = factor(n_dat)),
  #   alpha = 0.1,
  #   linewidth = 2,
  #   show.legend = FALSE
  # ) +
  geom_line(
    aes(x = cor_gauss1, y = mean_rand, colour = factor(n_dat)),
    linewidth = 1.2,
    # show.legend = FALSE
  ) +
  ggsci::scale_colour_nejm() +
  labs(
    y = "Adjusted Rand Index",
    x = "Gaussian corr",
    fill = "# points"
  ) +
  theme +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
p2

ggsave(
  plot = p2,
  paste0("plots/01d_js_sens_n_dat_dqu_", kl_prob, ".png")
  # paste0("plots/01b_js_sensitivity_dqu_", kl_prob, "_marg_0.9.png")
)

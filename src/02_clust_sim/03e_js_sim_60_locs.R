#### Simulation sensitivity analysis of Jensen-Shannon divergence method ####

# Assess performance of Jensen-Shannon clustering algorithm
# algorithm on simulation dataset consisting of mixture of Gaussian and
# t-copulas with GPD margins, through a grid search of
# reasonable parameter values (particularly copula corr pars)

# More "realistic" example:
# - 60 locations (similar to Ireland data)
# - 3 clusters (for 20 locations each), each with 30, 20, 10 members
# - Add slight perturbation to copula correlation parameters
# - As we have 3 clusters, need to keep Gaussian correlation fixed throughout
#   in order to be able to plot the results against the 3 t-copula corrs, which
#   we do a grid search over.

#### libs ####

library(dplyr, quietly = TRUE)
devtools::load_all("../evc")
library(tidyr)
library(ggplot2)
library(parallel)
# library(ggtern)

source("src/functions.R")


#### Metadata ####

seed_number <- 123
n <- 1e3 # number of samples to take
n_locs <- 60 # number of "locations"
# three clusters
n_clust <- 3
# # for a third of locations each (1st cluster has 30, 2nd 20, 3rd 10)
cluster_mem <- c(rep(1, 30), rep(2, 20), rep(3, 10))
n_vars <- 2 # number of variables per location (e.g. rain and windspeed)
perturb_cor <- TRUE # perturb correlation within cluster for each location
mix_p <- c("gauss_cop" = 0.5, "t_cop" = 0.5) # mixture percentages
# constant Gaussian copula correlation throughout, otherwise too hard to vis
cor_gauss <- 0.5
# GPD parameters
scale_gpd <- 1
shape_gpd <- -0.05
marg_prob <- 0.9
cond_prob <- 0.9
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
  # t-copula correlations (don't have overlapping values)
  cor_t1 = seq(0, 1, by = 0.3),
  cor_t2 = seq(0.1, 1, by = 0.3),
  cor_t3 = seq(0.2, 1, by = 0.3),
  # Degrees of freedom for t-copula
  df_t = 3,
  # mixture percentages (must sum to 1)
  mix_p1 = 0.5,
  mix_p2 = 0.5,
  # extremal quantiles
  cond_prob = cond_prob
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
    data_mix <- with(row, sim_cop_dat(
      n = 1e3,
      n_locs = n_locs,
      n_clust = n_clust,
      cluster_mem = cluster_mem,
      cor_gauss = rep(cor_gauss, n_clust),
      cor_t = c(cor_t1, cor_t2, cor_t3),
      df_t = rep(df_t, n_clust),
      # params_gpd  = c(scale_gpd, shape_gpd),
      mix_p = c(0.5, 0.5),
      perturb_cor = perturb_cor, # perturb cor for different locations
      qfun = evd::qgpd,
      qargs = c("scale" = scale_gpd, "shape" = shape_gpd)
    ))$data_mix

    clust_res <- tryCatch(
      {
        # if an error is produced, return a dummy list
        dependence <- fit_ce(
          data_mix,
          marg_prob   = marg_prob,
          cond_prob   = row$cond_prob,
          fit_no_keef = TRUE
        )
        js_clust(dependence, k = n_clust, cluster_mem = cluster_mem)
      },
      error = function(cond) {
        return(list("adj_rand" = NA))
      }
    )
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
readr::write_csv(
  results_grid,
  file = paste0(
    "data/js_grid_search_res_n_loc_", n_locs,
    "_dqu_", cond_prob,
    "_gauss_cor_", cor_gauss,
    ".csv"
  )
)

# # load data and redo results_grid_sum and results_grid_tally
# # May be floating point error in some parameter values due to loading from RDS!
results_grid <- readr::read_csv(
  paste0(
    "data/js_grid_search_res_n_loc_", n_locs,
    "_dqu_", cond_prob,
    "_gauss_cor_", cor_gauss,
    ".csv"
  )
)
results_grid_sum <- summarise_sens_res(results_grid, conf_level = conf_level)
results_grid_tally <- results_grid %>%
  group_by(across(!contains("_rand")), adj_rand) %>%
  tally(name = "n_rand")
results_grid <- results_grid %>%
  # for repeated joining
  dplyr::select(-(ends_with("_rand") & !matches("adj_rand"))) %>%
  left_join(results_grid_sum) %>%
  left_join(results_grid_tally)


#### Plot ####
# remove unlikely scenario of having perfect t-copula correlation
results_grid_plt <- filter(results_grid, cor_t1 < 1, cor_t2 < 1)
res_grid_sum_plt <- filter(results_grid_sum, cor_t1 < 1, cor_t2 < 1)

p1 <- results_grid_plt %>%
  ggplot() +
  geom_point(
    # aes(x = cor_t3, y = adj_rand),
    # use cor_t1 as x-axis since it's the only column with 4 unique values
    aes(x = cor_t1, y = adj_rand),
    colour = "black", alpha = 0.05, size = 1
  ) +
  # Confidence intervals around mean line
  geom_ribbon(
    data = res_grid_sum_plt,
    # aes(x = cor_t3, ymin = lower_rand, ymax = upper_rand),
    aes(x = cor_t1, ymin = lower_rand, ymax = upper_rand),
    fill = "#C11432",
    alpha = 0.3,
    size = 2
  ) +
  geom_line(
    data = res_grid_sum_plt,
    # aes(x = cor_t3, y = mean_rand),
    aes(x = cor_t1, y = mean_rand),
    colour = "#C11432",
    linewidth = 1
  )

# TODO save
# Can clearly see that clustering is best where correlation parameters in
# t-copula for different clusters are the most different
# Unsurprisingly, cor_tx = (0.1, 0.5, 0.9) has the best ARI, reflecting this
p1 <- common_plot_items(
  p1,
  # xlab = expression(rho[t[3]]),
  xlab = expression(rho[t[1]]),
  # x_sec_lab = expression(rho[t[2]]),
  x_sec_lab = expression(rho[t[3]]),
  # y_sec_lab = expression(rho[t[1]]),
  y_sec_lab = expression(rho[t[2]]),
  facet_form = (cor_t2 ~ cor_t3)
) +
  scale_x_continuous(breaks = unique(results_grid_plt$cor_t1)) +
  NULL

ggsave(
  plot = p1,
  # paste0("latex/plots/01e_js_sens_3_var_dqu_", cond_prob, ".png"),
  paste0("latex/plots/01c_js_sens_sim_60_locs_dqu_", cond_prob, ".png"),
  width = 10,
  height = 7
)


# Above plot, but this time with smoothed LOESS line rather than geom_line
smooth_plt <- \(data, col) {
  data %>%
    ggplot() +
    # plot points, and smooth line through them
    geom_point(
      # aes(x = .data[[col]], y = adj_rand),
      aes(x = {{ col }}, y = adj_rand),
      colour = "black", alpha = 0.05, size = 1
    ) +
    geom_smooth(
      # aes(x = .data[[col]], y = adj_rand),
      aes(x = {{ col }}, y = adj_rand),
      colour = "#C11432",
      fill = "#C11432",
      # formula = y ~ splines::bs(x = x, knots = 3),
      formula = y ~ splines::bs(x = x, knots = length(unique(data[[col]]))),
      se = TRUE,
      linewidth = 1
    ) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    # facet_grid(cor_t1 ~ cor_t2) +
    facet_grid(cor_t2 ~ cor_t3) +
    labs(y = "ARI", x = "") +
    evc_theme() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
    ggsci::scale_fill_nejm() +
    ggsci::scale_colour_nejm()
}

# doesn't work any more as cor_gauss is fixed!
p1_smooth <- smooth_plt(results_grid_plt, cor_gauss)
p2_smooth <- smooth_plt(results_grid_plt, cor_t1)

# ternary plot
# TODO: Fix (do I even bother??)
res_grid_sum_plt %>%
  dplyr::select(contains("cor_t"), median_rand) %>%
  ggtern::ggtern(aes(x = cor_t1, y = cor_t2, z = cor_t3, value = median_rand)) +
  ggtern::stat_interpolate_tern(
    geom = "polygon",
    formula = value ~ x + y,
    method = "gam",
    # Fill contours by interpolated levels of median_rand
    aes(fill = ..level..)
  ) +
  geom_point(aes(colour = median_rand), size = 2) + # Add points for reference
  ggtern::tern_limits(T = 1.05, L = 1.05, R = 1.05) +
  viridis::scale_fill_viridis(option = "magma") +
  viridis::scale_colour_viridis(option = "magma") +
  theme_classic() +
  theme(legend.position = "right")

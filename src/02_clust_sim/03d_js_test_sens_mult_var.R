#### Apply Jensen-Shannon Divergence method in multivariate (>2) case ####

# Assess performance of Jensen-Shannon clustering algorithm
# algorithm on simulation dataset consisting of mixture of Gaussian and
# t-copulas with GPD margins, through a grid search of
# reasonable parameter values (particularly copula corr pars)

# This time, we show how this method works for more than 2 variables, in
# contrast to the Vignotto method which only works in the bivariate case


#### Libs ####

library(dplyr, quietly = TRUE)
devtools::load_all("../evc")
library(tidyr)
library(ggplot2)
library(parallel)

source("src/functions.R")


#### Metadata ####

# metadata for original testing of method for > 2 vars
seed_number <- 123
# number of "locations"
n_locs <- 12
# Multivariate (> 2-dimensions)
vars <- c("rain", "wind_speed", "temperature") # dummy names for 3 vars
n_vars <- length(vars)
cluster_mem <- sort(rep(c(1, 2), 6)) # known cluster membership
n <- 1e4 # number of samples to take
mix_p <- c("gauss_cop" = 0.5, "t_cop" = 0.5) # mixture percentages
df_t <- 3 # degrees of freedom for t-copula
# GPD parameters
scale_gpd <- 1
shape_gpd <- -0.05
marg_prob <- 0.9 # marginal probability for GPD
cond_prob <- 0.9 # quantile for conditional probability
conf_level <- 0.9 # confidence level for CIs in plot

# Number of cores to use for parallel computation
n_cores <- detectCores() - 1


#### Test for > 2 variables ####

# First, we're going to test that the method works for > 2 variables

# Copula parameters for two clusters
cor_gauss <- c(0.5, 0.5)
cor_t <- c(0.9, 0.1)
df_t <- rep(df_t, 2)

# simulate data with three variables at each location
set.seed(seed_number)
data <- sim_cop_dat(
  n_locs     = n_locs,
  n_vars     = n_vars,
  n          = n,
  cor_gauss  = cor_gauss,
  cor_t      = cor_t,
  df_t       = df_t,
  # params_gpd = c(scale_gpd, shape_gpd),
  mix_p      = c(0.5, 0.5),
  qfun       = evd::qgpd,
  qargs      = c("scale" = scale_gpd, "shape" = shape_gpd)
)
data_mix <- data$data_mix

# Fit CE model
dependence <- fit_ce(
  data        = data_mix,
  vars        = vars,
  marg_prob   = marg_prob,
  cond_prob   = cond_prob,
  fit_no_keef = TRUE
)$dependence

# check that all dependence models have run successfully
sapply(dependence, \(x) lapply(x, length))


# Perform PAM and k-means clustering, as an exploratory analysis

# first, produce distance matrix and associated elbow plot
js_mat <- js_clust(dependence)$dist_mat # suggests k = 2, as desired

# cluster and assess performance
js_clust(dependence, k = 2, dist_mat = js_mat, cluster_mem = cluster_mem)
# adj_rand == 1, as desired (perfect clustering!)


#### Sensitivity Analysis ####

# create grid
grid <- tidyr::crossing(
  cor_gauss1 = seq(0, 1, by = 0.1),
  cor_gauss2 = seq(0, 1, by = 0.1),
  # t-copula correlation
  cor_t1 = seq(0.1, 0.9, by = 0.2),
  cor_t2 = seq(0, 1, by = 0.2),
  # Degrees of freedom for t-copula
  df_t1 = 3,
  df_t2 = 3,
  # mixture percentages (must sum to 1)
  mix_p1 = 0.5,
  mix_p2 = 0.5,
  # extremal quantiles
  kl_prob = cond_prob
) %>%
  filter(
    # use same correlation in both clusters for Gaussian copula
    cor_gauss1 == cor_gauss2
  )

# run kl_sim_eval for each row in grid
n_times <- 500
results_vec <- lri_vec <- lri_mean_vec <- vector(length = n_times)
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
      n_locs     = n_locs,
      n_vars     = n_vars,
      n          = 1e3,
      cor_gauss  = c(cor_gauss1, cor_gauss2),
      cor_t      = c(cor_t1, cor_t2),
      df_t       = c(df_t1, df_t2),
      # params_gpd = c(scale_gpd, shape_gpd),
      mix_p      = c(0.5, 0.5),
      qfun       = evd::qgpd,
      qargs      = c("scale" = scale_gpd, "shape" = shape_gpd)
    ))$data_mix

    # if an error is produced, return a dummy list
    clust_res <- tryCatch(
      {
        # fit CE model
        dependence <- fit_ce(
          data        = data_mix,
          vars        = vars,
          marg_prob   = marg_prob,
          cond_prob   = cond_prob,
          fit_no_keef = TRUE
        )$dependence

        # cluster based on results
        js_clust(dependence, k = 2, cluster_mem = cluster_mem)
      },
      error = function(cond) {
        return(list("adj_rand" = NA))
      }
    )
    results_vec[[j]] <- clust_res$adj_rand

    # also calculate local rand index
    lri <- local_rand_index(clust_res$pam$clustering, cluster_mem)
    # concatenate to string to store in vector
    lri_vec[[j]] <- paste0(lri, collapse = "_")
    # also average by cluster
    lri_mean_vec[[j]] <- paste0(vapply(unique(cluster_mem), \(x) {
      mean(lri[cluster_mem == x])
    }, numeric(1)), collapse = "-")
  }

  return(cbind(
    row,
    "adj_rand"        = results_vec,
    "local_rand"      = lri_vec,
    "mean_local_rand" = lri_mean_vec
  ))
  # }))
}, mc.cores = n_cores))

# remove NAs, give warning
nas <- is.na(results_grid$adj_rand)
if (sum(nas) > 0) {
  message(paste0("Warning: ", sum(nas), " NAs in results, removing"))
  results_grid <- results_grid[!nas, ]
}

# TODO Also functionalise? Repeated throughout scripts
# Summarise results to add confidence levels to plot
results_grid_sum <- summarise_sens_res(results_grid, conf_level = conf_level)
# tally rand index occurrences for given parameter set
results_grid_tally <- results_grid %>%
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
  paste0("data/js_sens_res_nvar_", n_vars, "_dqu_", cond_prob, ".csv")
)

# load data and redo results_grid_sum and results_grid_tally
results_grid <- readr::read_csv(
  paste0("data/js_sens_res_nvar_", n_vars, "_dqu_", cond_prob, ".csv")
)
results_grid_sum <- summarise_sens_res(results_grid, conf_level = conf_level)
results_grid_tally <- results_grid %>%
  group_by(across(!contains("_rand")), adj_rand) %>%
  tally(name = "n_rand")
# add both to results, as before
results_grid <- results_grid %>%
  dplyr::select(-(ends_with("_rand") & !matches("adj_rand"))) %>%
  left_join(results_grid_sum) %>%
  left_join(results_grid_tally)


#### Plot ####

# remove unlikely scenario of having perfect t-copula correlation
results_grid_plt <- filter(results_grid, cor_t1 < 1, cor_t2 < 1)
res_grid_sum_plt <- filter(results_grid_sum, cor_t1 < 1, cor_t2 < 1)

p1 <- results_grid_plt |>
  ggplot() +
  geom_point(
    aes(x = cor_gauss1, y = adj_rand),
    colour = "black", alpha = 0.05, size = 1
  ) +
  # Confidence intervals around mean line
  geom_ribbon(
    data = res_grid_sum_plt,
    aes(x = cor_gauss1, ymin = lower_rand, ymax = upper_rand),
    fill = "#C11432",
    alpha = 0.3,
    size = 2
  ) +
  geom_line(
    data = res_grid_sum_plt,
    aes(x = cor_gauss1, y = median_rand),
    colour = "#C11432",
    linewidth = 1
  )
(p1 <- common_plot_items(p1))

ggsave(
  plot = p1,
  # paste0("plots/01e_js_sens_3_var_dqu_", cond_prob, ".png"),
  paste0("latex/plots/sim_01d_js_sens_3_var_dqu_", cond_prob, ".png"),
  width = 10,
  height = 7
)


#### Compare to two variables ####

# TODO Bug here, data is the very same!
# results_grid2 <- readRDS("data/js_grid_search_res_dqu_0.9.RDS")
results_grid2 <- readr::read_csv("data/js_grid_search_res_dqu_0.9.csv")

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
  geom_line(
    aes(x = cor_gauss1, y = median_rand, colour = ind),
    linewidth = 1,
    show.legend = FALSE
  ) +
  # Confidence intervals around mean line
  geom_ribbon(
    aes(x = cor_gauss1, ymin = lower_rand, ymax = upper_rand, fill = ind),
    alpha = 0.3,
    size = 2,
  ) +
  guides(fill = guide_legend(override.aes = list(alpha = 1)))

(p2 <- common_plot_items(p2))

ggsave(
  plot = p2,
  # paste0("plots/01f_js_sens_3_var_compare_dqu_", cond_prob, ".png"),
  paste0("latex/plots/01f_js_sens_3_var_compare_dqu_", cond_prob, ".png"),
  width = 10,
  height = 7
)

#### Apply Jensen-Shannon Divergence method in multivariate (>2) case ####

# Assess performance of Jensen-Shannon clustering algorithm
# algorithm on simulation dataset consisting of mixture of Gaussian and
# t-copulas with GPD margins, through a grid search of
# reasonable parameter values (particularly copula corr pars)

# This time, we show how this method works for more than 2 variables, in
# contrast to the Vignotto method which only works in the bivariate case

# Previously, we used the same correlation parameter across all variable
# comparisons.
# Now, let's vary that! So, for three variables, we specify a 3 x 3 correlation
# matrix for the t-copula:
#  \[
#  R_t^{(s)} =
#    \begin{pmatrix}
#  1 & \rho^{(s)}_{12} & \rho^{(s)}_{13} \\
#  \rho^{(s)}_{12} & 1 & \rho^{(s)}_{23} \\
#  \rho^{(s)}_{13} & \rho^{(s)}_{23} & 1
#  \end{pmatrix}
#  \]
# So we specify rho_{12}, rho_{13} and rho_{23}, rather than a single row

# We keep rho_{Gauss} fixed throughout this experiment


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
n_vars <- 3
n_clust <- 2
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

# parameters for MC estimation of JSGa
n_mc <- 500 # number of MC samples
mc_method <- "laplace_trunc2" # truncate Laplace using empirical dist quant
laplace_cap <- 0.99 # Take 99th quantile


# Number of cores to use for parallel computation
n_cores <- detectCores() - 1


#### Test for 3 variables ####

# First, we're going to test that the method works for > 2 variables

# Copula parameters for two clusters
cor_gauss <- c(0.5, 0.5)
# cor_t <- c(0.9, 0.1)
cor_t <- list(
  c(0.8, 0.55, 0.3), # entries to correlation matrix
  c(0.55, 0.3, 0.8)
)
df_t <- rep(df_t, 2)

# simulate data with three variables at each location
set.seed(seed_number)
# data <- sim_cop_dat(
data <- sim_cop_dat_vec(
  n_locs     = n_locs,
  n_vars     = n_vars,
  n          = n,
  # cor_gauss  = cor_gauss,
  cor_gauss  = c(0.3, 0.3),
  cor_t      = cor_t,
  df_t       = df_t,
  # params_gpd = c(scale_gpd, shape_gpd),
  mix_p      = c(0.5, 0.5),
  qfun       = evd::qgpd,
  qargs      = c("scale" = scale_gpd, "shape" = shape_gpd)
)
data_gauss <- data$gauss_cop
data_t <- data$t_cop
data_mix <- data$data_mix

# plot
# plot(data_gauss[[12]])
# points(data_t[[12]], col = "red", pch = 19, cex = 0.5)

# Fit CE model
ce_fit <- fit_ce(
  data        = data_mix,
  # vars        = vars,
  marg_prob   = marg_prob,
  cond_prob   = cond_prob,
  fit_no_keef = TRUE,
  output_all  = TRUE
)
dependence <- ce_fit$dependence

# check that all dependence models have run successfully
sapply(dependence, \(x) lapply(x, length))

# Perform PAM and k-means clustering, as an exploratory analysis

# first, produce distance matrix and associated elbow plot
js_mat <- js_clust(
  ce_fit,
  n = n_mc, mc_method = mc_method, laplace_cap = laplace_cap
)$dist_mat # suggests k = 2, as desired

# cluster and assess performance
js_clust(
  dependence,
  k = 2, dist_mat = js_mat, cluster_mem = cluster_mem,
  n = n_mc, mc_method = mc_method, laplace_cap = laplace_cap
)
# adj_rand == 1, as desired (perfect clustering!)


#### Sensitivity Analysis ####

# patterns of trivariate extremal dependence within t-copula corr matrix
# TODO May have to soften if too many correlation matrices are removed
# low_val <- 0.1
# low_val <- 0.2
low_val <- 0.3
med_val <- 0.5
# high_val <- 0.8
# high_val <- 0.7
high_val <- 0.6
rho_patterns <- tibble::tribble(
  ~pattern, ~rho12, ~rho13, ~rho23,
  # homogeneous baseline
  "homog_low", low_val, low_val, low_val,
  "homog_mid", med_val, med_val, med_val,
  "homog_high", high_val, high_val, high_val,

  # one-high / one-low (mild heterogeneity)
  "one_high_12", high_val, med_val, low_val,
  "one_high_23", med_val, low_val, high_val,
  "one_high_13", low_val, high_val, med_val,

  # two-high / one-low
  "two_high_12_13", high_val, high_val, med_val,
  "two_high_12_23", high_val, med_val, high_val,
  "two_high_13_23", med_val, high_val, high_val
) |>
  rowwise() %>%
  # check that correlation matrices are valid
  mutate(valid = is_valid_cor_vec(c(rho12, rho13, rho23), n_vars = 3)) %>%
  ungroup() %>%
  filter(valid) |>
  select(-valid)

grid <- crossing(
  rename(
    rho_patterns,
    pattern1 = pattern,
    rho12_1 = rho12,
    rho13_1 = rho13,
    rho23_1 = rho23
  ),
  rename(
    rho_patterns,
    pattern2 = pattern,
    rho12_2 = rho12,
    rho13_2 = rho13,
    rho23_2 = rho23
  ),
  cor_gauss = seq(0, 1, by = 0.2)
) %>%
  # # Drop diagonal: same t-pattern in both clusters
  # filter(pattern1 != pattern2) %>%
  mutate(
    df_t1 = 3,
    df_t2 = 3,
    mix_p1 = 0.5,
    mix_p2 = 0.5,
    kl_prob = cond_prob
  )

# Number of simulation scenarios
nrow(grid)

# run kl_sim_eval for each row in grid
# n_times <- 500
n_times <- 500
# results_vec <- lri_vec <- lri_mean_vec <- vector(length = n_times)
results_vec <- vector(length = n_times)
set.seed(seed_number)
# i <- 11
i <- 1
results_grid <- bind_rows(mclapply(seq_len(nrow(grid)), \(i) {
  print(paste0("Progress: ", round(i / nrow(grid), 3) * 100, "%"))
  system(sprintf(
    'echo "\n%s\n"',
    paste0(round(i / nrow(grid), 3) * 100, "% completed", collapse = "")
  ))

  row <- grid[i, , drop = FALSE]
  for (j in seq_len(n_times)) {
    # generate simulation data for given parameter set
    data_mix <- with(row, sim_cop_dat_vec(
      n_vars = n_vars,
      n_locs = n_locs,
      n = n,
      cor_gauss = c(cor_gauss, cor_gauss),
      # cor_t      = c(cor_t1, cor_t2),
      cor_t = list(
        c(rho12_1, rho13_1, rho23_1),
        c(rho12_2, rho13_2, rho23_2)
      ),
      df_t = c(df_t1, df_t2),
      # params_gpd = c(scale_gpd, shape_gpd),
      mix_p = c(0.5, 0.5),
      qfun = evd::qgpd,
      qargs = c("loc" = 0, "scale" = scale_gpd, "shape" = shape_gpd)
    ))$data_mix

    # transform to Laplace margins
    data_mix_trans <- lapply(data_mix, \(x) {
      laplace_trans(evd::pgpd(x, loc = 0, scale = scale_gpd, shape = shape_gpd))
    })

    clust_res <- tryCatch(
      {
        ce_fit <- lapply(seq_along(data_mix_trans), \(l) {
          o <- ce_optim(
            Y         = data_mix_trans[[l]],
            dqu       = row$kl_prob,
            control   = list(maxit = 1e6),
            constrain = FALSE
          )
        })

        # "true" number of clusters known from simulation design
        js_clust(
          ce_fit,
          trans = data_mix_trans,
          k = n_clust, cluster_mem = cluster_mem,
          n = n_mc, mc_method = mc_method, laplace_cap = laplace_cap
        )
      },
      error = function(cond) {
        return(list("adj_rand" = NA))
      }
    )
    results_vec[[j]] <- clust_res$adj_rand

    # # also calculate local rand index
    # lri <- local_rand_index(clust_res$pam$clustering, cluster_mem)
    # # concatenate to string to store in vector
    # lri_vec[[j]] <- paste0(lri, collapse = "_")
    # # also average by cluster
    # lri_mean_vec[[j]] <- paste0(vapply(unique(cluster_mem), \(x) {
    #   mean(lri[cluster_mem == x])
    # }, numeric(1)), collapse = "-")
  }

  return(cbind(
    row,
    "adj_rand" = results_vec # ,
    # "local_rand"      = lri_vec,
    # "mean_local_rand" = lri_mean_vec
  ))
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
  paste0("data/js_sens_hetero_nvar_", n_vars, "_dqu_", cond_prob, ".csv")
)

# load data and redo results_grid_sum and results_grid_tally
results_grid <- readr::read_csv(
  paste0("data/js_sens_hetero_nvar_", n_vars, "_dqu_", cond_prob, ".csv")
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

results_grid_plt <- results_grid %>%
  mutate(
    facet_lab1 = glue::glue("({rho12_1}, {rho13_1}, {rho23_1})"),
    facet_lab2 = glue::glue("({rho12_2}, {rho13_2}, {rho23_2})")
  )

p1 <- results_grid_plt %>%
  ggplot(aes(x = cor_gauss, y = mean_rand)) +
  geom_ribbon(
    aes(ymin = lower_rand, ymax = upper_rand),
    fill = "#C11432",
    alpha = 0.25
  ) +
  geom_line(
    colour = "#C11432",
    linewidth = 1
  ) +
  facet_grid(facet_lab1 ~ facet_lab2) +
  # scale_x_continuous(
  #   sec.axis = sec_axis(
  #     ~ .,
  #     name = x_sec_lab,
  #     breaks = NULL,
  #     labels = NULL
  #   )
  # ) +
  scale_y_continuous(
    # sec.axis = sec_axis(
    #   ~ .,
    #   name = y_sec_lab,
    #   breaks = NULL,
    #   labels = NULL
    # ),
    breaks = seq(0, 1, by = 0.2),
    # breaks = seq(-0.2, 1, by = 0.2),
    limits = c(-0.2, 1)
  ) +
  labs(
    x = expression(rho["Gauss"]),
    y = "ARI",
    fill = "",
    colour = ""
  ) +
  evc_theme() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5),
    strip.text = element_text(size = 8)
  ) +
  ggsci::scale_fill_nejm() +
  ggsci::scale_colour_nejm()
p1

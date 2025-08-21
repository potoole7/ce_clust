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
# vars <- c("rain", "wind_speed", "temperature") # dummy names for 3 vars
# n_vars <- length(vars)
n_vars <- 2
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

bind_rows(
  bind_rows(
    data.frame(data_gauss[[1]]) |>
      mutate(name = "Gaussian Copula"),
    data.frame(data_t[[1]]) |>
      mutate(name = "t-Copula")
  ) |>
    mutate(cor_gauss = cor_gauss[1], cor_t = cor_t[1]),
  bind_rows(
    data.frame(data_gauss[[12]]) |>
      mutate(name = "Gaussian Copula"),
    data.frame(data_t[[12]]) |>
      mutate(name = "t-Copula")
  ) |>
    mutate(cor_gauss = cor_gauss[1], cor_t = cor_t[2]),
) |>
  group_by(cor_t) |>
  mutate(
    quant_x = quantile(X1, 0.99),
    quant_y = quantile(X2, 0.99)
  ) |>
  ungroup() |>
  ggplot(aes(x = X1, y = X2, colour = name)) +
  geom_point(alpha = 0.4) +
  facet_wrap(~cor_t, scale = "fixed") +
  geom_vline(aes(xintercept = quant_x), linetype = "dashed") +
  geom_hline(aes(yintercept = quant_y), linetype = "dashed") +
  labs(x = "", y = "") +
  guides(colour = "none") +
  evc::evc_theme() +
  scale_x_continuous(
    sec.axis = sec_axis(
      ~.,
      name = expression(rho[t]),
      breaks = NULL,
      labels = NULL
    )
  )

# Fit CE model
# TODO
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
# i <- 11
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
      n_vars     = n_vars,
      n_locs     = n_locs,
      n          = n,
      cor_gauss  = c(cor_gauss1, cor_gauss2),
      cor_t      = c(cor_t1, cor_t2),
      df_t       = c(df_t1, df_t2),
      # params_gpd = c(scale_gpd, shape_gpd),
      mix_p      = c(0.5, 0.5),
      qfun       = evd::qgpd,
      qargs      = c("loc" = 0, "scale" = scale_gpd, "shape" = shape_gpd)
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
  # geom_ribbon(
  #   data = res_grid_sum_plt,
  #   aes(x = cor_gauss1, ymin = lower_rand, ymax = upper_rand),
  #   fill = "#C11432",
  #   alpha = 0.3,
  #   size = 2
  # ) +
  # geom_line(
  #   data = res_grid_sum_plt,
  #   aes(x = cor_gauss1, y = median_rand),
  #   colour = "#C11432",
  #   linewidth = 1
  # ) +
  geom_smooth(
    aes(x = cor_gauss1, y = adj_rand),
    method = "loess",
    colour = "#C11432",
    se = TRUE,
    alpha = 0.7,
    linewidth = 1.2
  ) +
  NULL
(p1 <- common_plot_items(p1))

# TODO Save manually to avoid clipping
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
  # points
  # geom_point(aes(x = cor_gauss1, y = mean_rand, colour = ind), size = 2) +
  # Confidence intervals around mean line
  # geom_ribbon(
  #   aes(x = cor_gauss1, ymin = lower_rand, ymax = upper_rand, fill = ind),
  #   # fill = "#C11432",
  #   alpha = 0.3,
  #   size = 2
  # ) +
  # geom_line(
  #   aes(x = cor_gauss1, y = mean_rand, colour = ind),
  #   # colour = "#C11432",
  #   linewidth = 1
  # ) +
  # # TODO Add colour and alpha to uncertainty
  geom_smooth(
    aes(x = cor_gauss1, y = adj_rand, colour = ind),
    method = "loess",
    # colour = "#C11432",
    # se = TRUE,
    se = FALSE,
    alpha = 0.7,
    linewidth = 1.2,
    key_glyph = draw_key_rect # have rectangle in legend, rather than line
  ) +
  ggsci::scale_colour_nejm() +
  guides(
    colour = guide_legend(override.aes = list(
      fill     = ggsci::pal_nejm()(2),
      colour   = NA, # no border
      size     = 5,
      linetype = 0, # hide line
      alpha    = 1
    ))
  )
common_plot_items(p2)

ggsave(
  plot = p2,
  # paste0("plots/01f_js_sens_3_var_compare_dqu_", cond_prob, ".png"),
  paste0("latex/plots/01f_js_sens_3_var_compare_dqu_", cond_prob, ".png"),
  width = 10,
  height = 7
)

#### Grid search for 5, 8, 10 variables ####

# create grid
grid_mult <- tidyr::crossing(
  cor_gauss1 = seq(0, 1, by = 0.1),
  cor_gauss2 = seq(0, 1, by = 0.1),
  # t-copula correlation
  # cor_t1 = c(0.6),
  cor_t1 = c(0.4),
  # cor_t2 = c(0.7),
  cor_t2 = c(0.5),
  # Degrees of freedom for t-copula
  df_t1 = 3,
  df_t2 = 3,
  # mixture percentages (must sum to 1)
  mix_p1 = 0.5,
  mix_p2 = 0.5,
  # extremal quantiles
  kl_prob = cond_prob,
  # n_vars = c(5, 8, 10) # number of variables to simulate
  n_vars = c(2, 3, 4, 5)
) %>%
  filter(
    # use same correlation in both clusters for Gaussian copula
    cor_gauss1 == cor_gauss2,
    cor_t1 == cor_t2 - 0.1 # only want most difficult case where cors are close
  )

# run kl_sim_eval for each row in grid
# n_times <- 1
# n_times <- 10
n_times <- 200
results_vec <- lri_vec <- lri_mean_vec <- vector(length = n_times)
set.seed(seed_number)

# i <- 11
results_grid_mult <- bind_rows(mclapply(seq_len(nrow(grid_mult)), \(i) {
  print(paste0("Progress: ", round(i / nrow(grid_mult), 3) * 100, "%"))
  system(sprintf(
    'echo "\n%s\n"',
    paste0(round(i / nrow(grid_mult), 3) * 100, "% completed", collapse = "")
  ))

  row <- grid_mult[i, , drop = FALSE]
  for (j in seq_len(n_times)) {
    # generate simulation data for given parameter set
    data_mix <- with(row, sim_cop_dat(
      n_vars     = n_vars,
      n_locs     = n_locs,
      n          = n,
      cor_gauss  = c(cor_gauss1, cor_gauss2),
      cor_t      = c(cor_t1, cor_t2),
      df_t       = c(df_t1, df_t2),
      # params_gpd = c(scale_gpd, shape_gpd),
      mix_p      = c(0.5, 0.5),
      qfun       = evd::qgpd,
      qargs      = c("loc" = 0, "scale" = scale_gpd, "shape" = shape_gpd)
    ))$data_mix

    # transform to Laplace margins
    data_mix_trans <- lapply(data_mix, \(x) {
      laplace_trans(evd::pgpd(x, loc = 0, scale = scale_gpd, shape = shape_gpd))
    })

    clust_res <- tryCatch(
      {
        ce_fit <- lapply(seq_along(data_mix_trans), \(l) {
          # TODO Investigate warnings for i = 1!
          o <- ce_optim(
            Y         = data_mix_trans[[l]],
            dqu       = row$kl_prob,
            control   = list(maxit = 1e6),
            constrain = FALSE
          )
        })

        # "true" number of clusters known from simulation design
        js_clust(ce_fit, k = n_clust, cluster_mem = cluster_mem)
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

# save
# readr::write_csv(
#   results_grid_mult,
#   paste0("data/js_sens_res_nvar_", n_vars, "_dqu_", cond_prob, "mult_var.csv")
# )

results_grid_mult <- readr::read_csv(
  paste0("data/js_sens_res_nvar_", n_vars, "_dqu_", cond_prob, "mult_var.csv")
)

# plot
results_grid_mult_sum <- summarise_sens_res(results_grid_mult, conf_level = conf_level)
results_grid_mult_tally <- results_grid_mult %>%
  group_by(across(!contains("_rand")), adj_rand) %>%
  tally(name = "n_rand")
# add both to results, as before
results_grid_mult <- results_grid_mult %>%
  dplyr::select(-(ends_with("_rand") & !matches("adj_rand"))) %>%
  left_join(results_grid_mult_sum) %>%
  left_join(results_grid_mult_tally)

results_grid_mult_plt <- filter(results_grid_mult, cor_t1 < 1, cor_t2 < 1)
res_grid_sum_plt <- filter(results_grid_mult_sum, cor_t1 < 1, cor_t2 < 1)

p1_mult <- results_grid_mult_plt |>
  ggplot() +
  # geom_point(
  #   aes(
  #     x = cor_gauss1, y = adj_rand, colour = factor(n_vars),
  #     # colour = "black",
  #     alpha = 0.05, size = 1
  #   )
  # ) +
  # Confidence intervals around mean line
  # geom_ribbon(
  #   data = res_grid_sum_plt,
  #   aes(x = cor_gauss1, ymin = lower_rand, ymax = upper_rand, fill = factor(n_vars)),
  #   # fill = "#C11432",
  #   alpha = 0.3,
  #   size = 2
  # ) +
  # geom_line(
  #   data = res_grid_sum_plt,
  #   aes(x = cor_gauss1, y = median_rand, colour = factor(n_vars)),
  #   # colour = "#C11432",
  #   linewidth = 1,
  #   alpha = 0.7
  # ) +
  geom_smooth(
    aes(x = cor_gauss1, y = adj_rand, colour = factor(n_vars), fill = factor(n_vars)),
    method = "loess",
    # colour = "#C11432",
    se = TRUE,
    alpha = 0.4,
    linewidth = 1.2
  ) +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  ggsci::scale_colour_nejm() +
  evc::evc_theme() +
  labs(
    x      = expression(rho["Gauss"]),
    y      = "ARI",
    fill   = "Copula Dimension",
    colour = "Copula Dimension"
  ) +
  scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.01))
  )
p1_mult

ggsave(
  plot = p1_mult,
  # paste0("plots/01e_js_sens_3_var_dqu_", cond_prob, ".png"),
  paste0("latex/plots/sim_01d_js_sens_3_var_dqu_", cond_prob, "mult_var.png"),
  width = 10,
  height = 7
)

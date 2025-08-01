#### Compare new and old clustering methods ####

# Assess performance of Jensen-Shannon clustering algorithm vs Vignotto
# algorithm on simulation dataset consisting of mixture of Gaussian and
# t-copulas with GPD margins, through a grid search of
# reasonable parameter values (particularly copula corr pars)

#### libs ####

library(dplyr, quietly = TRUE)
# library(evc)
devtools::load_all("../evc")
# devtools::load_all("../evc2")
library(tidyr)
library(ggplot2)
library(parallel)

# source functions
source("src/functions.R")


#### Metadata ####

seed_number <- 123
n <- 1e3 # number of samples to take
# n_locs <- 12 # number of "locations"
n_locs <- 30
n_clust <- 2 # number of clusters
cluster_mem <- sort(rep(1:n_clust, n_locs / 2)) # two clusters for half of locations each
n_vars <- 2 # number of variables per location
mix_p <- c("gauss_cop" = 0.5, "t_cop" = 0.5) # mixture percentages
# GPD parameters
scale_gpd <- 1
shape_gpd <- -0.05
df_t <- 3 # degrees of freedom for t-copulas
marg_prob <- 0.9 # marginal probability for GPD
kl_prob <- 0.9 # quantile for conditional probability
# TODO Change confidence intervals to use quantiles rather than normal assumption (??)
conf_level <- 0.9 # confidence level for CIs in plot
n_mc <- 50 # number of Monte Carlo simulations to use

# Number of cores to use for parallel computation
n_cores <- detectCores() - 1


#### Grid search ####

# create grid
grid <- tidyr::crossing(
  # Don't need all cor_vals for both (i.e. (0.3, 0.8) is the same as (0.8, 0.3)
  cor_gauss1 = seq(0, 1, by = 0.1),
  # cor_gauss1 = 0.7,
  cor_gauss2 = seq(0, 1, by = 0.1),
  # cor_gauss2 = 0.7,
  # t-copula correlation
  cor_t1 = seq(0.1, 0.9, by = 0.2),
  # cor_t1 = seq(0.3, 0.7, by = 0.2),
  # cor_t1 = 0.7,
  cor_t2 = seq(0, 0.8, by = 0.2),
  # cor_t2 = seq(0.2, 0.8, by = 0.2),
  # Degrees of freedom for t-copula
  df_t1 = df_t,
  df_t2 = df_t,
  # mixture percentages (must sum to 1)
  mix_p1 = 0.5,
  mix_p2 = 0.5,
  # extremal quantile
  kl_prob = kl_prob
) %>%
  # use same correlation in both clusters for Gaussian copula
  filter(cor_gauss1 == cor_gauss2)

n_cores <- min(c(n_cores, nrow(grid)))

# run kl_sim_eval for each row in grid
# n_times <- 500
n_times <- 200
# n_times <- 100
results_vec_orig <- results_vec_laplace <- results_vec_lap_trunc <- results_vec_lap_trunc2 <-
  results_vec_emp <- results_vec_hybrid <-
  lri_vec <- lri_mean_vec <- vector(length = n_times)
set.seed(seed_number)
# i <- nrow(grid) - 1
# i <- 5
results_grid <- bind_rows(mclapply(seq_len(nrow(grid)), \(i) {
  # results_grid <- bind_rows(lapply(seq_len(nrow(grid)), \(i) {
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

    ce_fit <- lapply(seq_along(data_mix_trans), \(l) {
      o <- ce_optim(
        Y         = data_mix_trans[[l]],
        dqu       = row$kl_prob,
        control   = list(maxit = 1e6),
        constrain = FALSE
      )
    })
    clust_res_orig <- js_clust(
      ce_fit,
      k = n_clust, cluster_mem = cluster_mem, n = 1000, mc_method = "uniform" # use_laplace = FALSE
    )
    clust_res_laplace <- js_clust(
      ce_fit,
      k = n_clust, cluster_mem = cluster_mem, n = 1000, mc_method = "laplace"
    )
    clust_res_lap_trunc <- js_clust(
      ce_fit,
      k = n_clust, cluster_mem = cluster_mem, n = 1000, mc_method = "laplace_trunc",
      laplace_cap = 0.95
    )
    clust_res_lap_trunc2 <- js_clust(
      ce_fit,
      k = n_clust, cluster_mem = cluster_mem, n = 1000, mc_method = "laplace_trunc",
      laplace_cap = 0.975
    )
    clust_res_emp <- js_clust(
      ce_fit,
      trans = data_mix_trans,
      k = n_clust, cluster_mem = cluster_mem, n = 1000, mc_method = "empirical" # use_laplace = TRUE
    )
    clust_res_hybrid <- js_clust(
      ce_fit,
      trans = data_mix_trans,
      k = n_clust, cluster_mem = cluster_mem, n = 1000, mc_method = "hybrid" # use_laplace = TRUE
    )

    # if (clust_res_orig$adj_rand != clust_res_new$adj_rand) {
    #   message("Warning: Different adj_rand values for original and new clustering methods!")
    #   browser()
    #   saveRDS(data_mix_trans, "data_mix_trans_test.RDS")
    #   saveRDS(row, "row_test.RDS")
    # }

    results_vec_orig[[j]] <- clust_res_orig$adj_rand
    results_vec_laplace[[j]] <- clust_res_laplace$adj_rand
    results_vec_lap_trunc[[j]] <- clust_res_lap_trunc$adj_rand
    results_vec_lap_trunc2[[j]] <- clust_res_lap_trunc2$adj_rand
    results_vec_emp[[j]] <- clust_res_emp$adj_rand
    results_vec_hybrid[[j]] <- clust_res_hybrid$adj_rand

    # also calculate local rand index
    # TODO Functionalise? Repeated in other scripts
    # lri <- local_rand_index(cluster_mem, clust_res$pam$clustering)
    # # concatenate to string to store in vector
    # lri_vec[[j]] <- paste0(lri, collapse = "_")
    # # also average by cluster
    # lri_mean_vec[[j]] <- paste0(vapply(unique(cluster_mem), \(x) {
    #   mean(lri[cluster_mem == x])
    # }, numeric(1)), collapse = "-")
  }

  return(cbind(
    row,
    "adj_rand_orig" = results_vec_orig,
    "adj_rand_laplace" = results_vec_laplace,
    "adj_rand_lap_trunc" = results_vec_lap_trunc,
    "adj_rand_lap_trunc2" = results_vec_lap_trunc2,
    "adj_rand_emp" = results_vec_emp,
    "adj_rand_hybrid" = results_vec_hybrid
    # "local_rand"      = lri_vec,
    # "mean_local_rand" = lri_mean_vec
  ))
}, mc.cores = n_cores))
# }))

# remove NAs, give warning
nas <- is.na(results_grid$adj_rand)
sum(nas)
if (sum(nas) > 0) {
  message(paste0("Warning: ", sum(nas), " NAs in results, removing"))
  results_grid <- results_grid[!nas, ]
}

# Summarise results to add confidence levels to plot
# results_grid_sum <- summarise_sens_res(results_grid, conf_level = conf_level)
# # tally rand index occurrences for given parameter set
# results_grid_tally <- results_grid %>%
#   group_by(across(!contains("_rand")), adj_rand) %>%
#   tally(name = "n_rand") %>%
#   ungroup()
# # add both to results
# results_grid <- results_grid %>%
#   # for repeated joining
#   dplyr::select(-(ends_with("_rand") & !matches("adj_rand"))) %>%
#   left_join(results_grid_sum) %>%
#   left_join(results_grid_tally)

# TODO save
# readr::write_csv(results_grid, "js_new.csv.gz")
# readr::write_csv(results_grid, "js_new_emp.csv.gz")
# readr::write_csv(results_grid, "js_new_emp_trunc.csv.gz")
readr::write_csv(results_grid, "js_new_emp_trunc_mult.csv.gz")

# results_grid <- readr::read_csv("js_new_emp_trunc_mult.csv.gz")


#### Plot ####


# which method is highest most frequently?
adj_rand_cols <- results_grid[, grepl("^adj_rand", names(results_grid))]

(comparison_counts <- sort(sapply(names(adj_rand_cols), function(col) {
  row_sums <- apply(adj_rand_cols, 1, function(row) {
    sum(row[col] > row[names(adj_rand_cols) != col])
  })
  sum(row_sums)
}), decreasing = TRUE))


res_plt <- results_grid |>
  # select(-adj_rand) |>
  # pivot_longer(c(adj_rand_orig, adj_rand_new, adj_rand_emp),
  pivot_longer(c(adj_rand_orig, adj_rand_emp, adj_rand_hybrid, adj_rand_laplace, adj_rand_lap_trunc, adj_rand_lap_trunc2),
    names_to = "ind", values_to = "adj_rand"
  ) |>
  mutate(
    ind = case_when(
      ind == "adj_rand_orig" ~ "Original",
      ind == "adj_rand_emp" ~ "Empirical",
      ind == "adj_rand_hybrid" ~ "Hybrid",
      ind == "adj_rand_laplace" ~ "Laplace",
      ind == "adj_rand_lap_trunc" ~ "Laplace Trunc 0.95",
      ind == "adj_rand_lap_trunc2" ~ "Laplace Trunc 0.975",
      TRUE ~ ind
    )
  )

common_plot_items(
  res_plt |>
    ggplot() +
    geom_smooth(
      aes(x = cor_gauss1, y = adj_rand, colour = ind, fill = ind),
      alpha = 0.5,
      method = "loess",
      # colour = "#C11432",
      se = TRUE,
      linewidth = 1.5
    )
)

# would also like to split and plot individually
plot_split <- group_split(res_plt, cor_t1, cor_t2, .keep = TRUE) |>
  lapply(\(x) {
    common_plot_items(
      x |>
        ggplot() +
        geom_smooth(
          aes(x = cor_gauss1, y = adj_rand, colour = ind, fill = ind),
          alpha = 0.5,
          method = "loess",
          # colour = "#C11432",
          se = TRUE,
          linewidth = 1.5
        ) +
        scale_y_continuous(limits = c(0, 1))
    )
  })

# what about plotting the median?
res_plt_med <- res_plt |>
  group_by(cor_gauss1, cor_t1, cor_t2, ind) |>
  summarise(
    median_rand = median(adj_rand, na.rm = TRUE),
    lower_rand = quantile(adj_rand, 0.05, na.rm = TRUE),
    upper_rand = quantile(adj_rand, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

res_plt_med %>%
  group_by(cor_gauss1, cor_t1, cor_t2) %>%
  filter(median_rand == max(median_rand)) %>%
  ungroup() |>
  count(ind, sort = TRUE)

common_plot_items(
  res_plt_med |>
    ggplot() +
    # geom_smooth(
    geom_line(
      aes(x = cor_gauss1, y = median_rand, colour = ind),
      alpha = 0.5,
      # method = "loess",
      # colour = "#C11432",
      # se = TRUE,
      linewidth = 1.5
    ) +
    # geom_ribbon(
    #   aes(x = cor_gauss1, ymin = lower_rand, ymax = upper_rand, fill = ind),
    #   alpha = 0.3
    # ) +
    NULL
)




# remove unlikely scenario of having perfect t-copula correlation
results_grid_plt <- filter(results_grid, cor_t1 < 1, cor_t2 < 1)
res_grid_sum_plt <- filter(results_grid_sum, cor_t1 < 1, cor_t2 < 1)

p1 <- results_grid_plt |>
  ggplot() +
  geom_point(
    aes(x = cor_gauss1, y = adj_rand),
    colour = "black",
    alpha = 0.05,
    # size = 1
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
  # Confidence intervals around mean line
  # geom_ribbon(
  #   data = res_grid_sum_plt, # %>% filter(cor_t1 == 0.3, cor_t2 == 0.6),
  #   aes(x = cor_gauss1, ymin = lower_rand, ymax = upper_rand),
  #   fill = "#C11432",
  #   alpha = 0.3,
  #   size = 2
  # ) +
  # geom_line(
  #   data = res_grid_sum_plt, # %>% filter(cor_t1 == 0.3, cor_t2 == 0.6),
  #   aes(x = cor_gauss1, y = median_rand),
  #   colour = "#C11432",
  #   linewidth = 1
  # ) +
  geom_smooth(
    aes(x = cor_gauss1, y = adj_rand),
    method = "loess",
    colour = "#C11432",
    se = TRUE,
    linewidth = 1
  ) +
  NULL
common_plot_items(p1)


# load original results
set.seed(seed_number)
results_grid_orig <- readr::read_csv("data_pre_mc/js_grid_search_res_dqu_0.9.csv") |>
  mutate(across(contains("cor_"), \(x) round(x, 1))) |>
  filter(cor_t1 == grid$cor_t1, cor_t2 == grid$cor_t2) |>
  # randomly take 200 simulations
  # group_by(cor_gauss1, cor_gauss2) |>
  # slice_sample(n = n_times, replace = FALSE) |>
  # ungroup() |>
  identity()

results_grid_join <- bind_rows(
  results_grid %>%
    # mutate(ind = "CE") %>%
    mutate(ind = "MC") %>%
    relocate(ind) |>
    identity(),
  results_grid_orig %>%
    mutate(ind = "Old") %>%
    relocate(ind) %>%
    dplyr::select(-matches("mean_rand")) %>%
    identity()
) %>%
  # fix weird floating point problem from loading from .RDS file
  mutate(across(contains("cor"), \(x) round(x, 1))) |>
  # remove unlikely scenario of having perfect t-copula correlation
  filter(cor_t1 < 1, cor_t2 < 1) |>
  identity()

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

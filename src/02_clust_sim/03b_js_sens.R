#### Simulation sensitivity analysis of Jensen-Shannon divergence method ####

# Assess performance of Jensen-Shannon clustering algorithm vs Vignotto
# algorithm on simulation dataset consisting of mixture of Gaussian and
# t-copulas with GPD margins, through a grid search of
# reasonable parameter values (particularly copula corr pars)

#### libs ####

library(dplyr, quietly = TRUE)
# library(evc)
devtools::load_all("../evc")
library(tidyr)
library(ggplot2)
library(parallel)

# source functions
source("src/functions.R")


#### Metadata ####

seed_number <- 123
n <- 1e3 # number of samples to take
n_locs <- 12 # number of "locations"
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

# Number of cores to use for parallel computation
n_cores <- detectCores() - 1


#### Grid search ####

# create grid
grid <- tidyr::crossing(
  # Don't need all cor_vals for both (i.e. (0.3, 0.8) is the same as (0.8, 0.3)
  cor_gauss1 = seq(0, 1, by = 0.1),
  cor_gauss2 = seq(0, 1, by = 0.1),
  # t-copula correlation
  cor_t1 = seq(0.1, 0.9, by = 0.2),
  cor_t2 = seq(0, 1, by = 0.2),
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
      n          = n,
      cor_gauss  = c(cor_gauss1, cor_gauss2),
      cor_t      = c(cor_t1, cor_t2),
      df_t       = c(df_t1, df_t2),
      params_gpd = c(scale_gpd, shape_gpd),
      mix_p      = c(0.5, 0.5)
    ))$data_mix

    clust_res <- tryCatch(
      {
        # if an error is produced, return a dummy list
        dependence <- fit_ce(
          data_mix,
          marg_prob   = marg_prob,
          cond_prob   = row$kl_prob,
          fit_no_keef = TRUE
        )$dependence
        # "true" number of clusters known from simulation design
        js_clust(dependence, k = n_clust, cluster_mem = cluster_mem)
      },
      error = function(cond) {
        return(list("adj_rand" = NA))
      }
    )
    results_vec[[j]] <- clust_res$adj_rand

    # also calculate local rand index
    # TODO Functionalise? Repeated in other scripts
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
  paste0("data/js_grid_search_res_dqu_", kl_prob, ".csv")
)

# load data and redo results_grid_sum and results_grid_tally
results_grid <- readr::read_csv(
  paste0("data/js_grid_search_res_dqu_", kl_prob, ".csv")
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

# TODO Add visualisations for Local Rand Index!

# remove unlikely scenario of having perfect t-copula correlation
results_grid_plt <- filter(results_grid, cor_t1 < 1, cor_t2 < 1)
res_grid_sum_plt <- filter(results_grid_sum, cor_t1 < 1, cor_t2 < 1)

# TODO Improvements to plotting labels
# library(latex2exp)
# library(ggtext)

# TODO How to structure confidence intervals? CLI, credible interval? LOESS?
# TODO How to label facets? If with maths symbols? How to match LaTeX font?
# TODO Many elements common across other plots, functionalise
p1 <- results_grid_plt |>
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
    aes(x = cor_gauss1, y = median_rand),
    colour = "#C11432",
    linewidth = 1
  ) +
  # geom_smooth(
  #   aes(x = cor_gauss1, y = adj_rand),
  #   method = "loess",
  #   colour = "#C11432",
  #   se = TRUE,
  #   linewidth = 1
  # ) +
  NULL

ggsave(
  plot = common_plot_items(p1),
  # paste0("plots/01b_js_sensitivity_dqu_", kl_prob, ".png"),
  paste0("latex/plots/sim_01a_js_sensitivity_dqu_", kl_prob, ".png"),
  width = 10,
  height = 7
)

# (try) save using LaTeX font for maths symbols
# tikzDevice::tikz(
#   # file = paste0("plots/01b_js_sensitivity_dqu_", kl_prob, ".tex"),
#   file = paste0("test.tex"),
#   width = 10,
#   height = 7
# )
# p1
# tikzDevice::dev.off()


#### Compare to Vignotto ####

# load Vignotto results
# results_grid_vig <- readRDS("data/vignotto_grid_search_res.RDS")
results_grid_vig <- readr::read_csv("data/vignotto_grid_search_res.csv")

results_grid_join <- bind_rows(
  results_grid %>%
    distinct() %>%
    mutate(ind = "CE") %>%
    relocate(ind),
  results_grid_vig %>%
    mutate(ind = "Vignotto") %>%
    relocate(ind) %>%
    dplyr::select(-matches("mean_rand")) %>%
    summarise_sens_res() %>%
    distinct()
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
  geom_ribbon(
    aes(x = cor_gauss1, ymin = lower_rand, ymax = upper_rand, fill = ind),
    # fill = "#C11432",
    alpha = 0.3,
    size = 2
  ) +
  geom_line(
    aes(x = cor_gauss1, y = mean_rand, colour = ind),
    # colour = "#C11432",
    linewidth = 1
  ) +
  # # TODO Add colour and alpha to uncertainty
  # geom_smooth(
  #   aes(x = cor_gauss1, y = adj_rand, colour = ind),
  #   method = "loess",
  #   colour = "#C11432",
  #   se = TRUE,
  #   linewidth = 1
  # ) +
  # ovverride alpha for fill
  guides(
    colour = "none",
    fill   = guide_legend(override.aes = list(alpha = 1))
  )

ggsave(
  plot = common_plot_items(p2),
  # paste0("plots/01c_sensitivity_dqu_", kl_prob, ".png"),
  paste0("latex/plots/sim_01b_ce_vs_vi_dqu_", kl_prob, ".png"),
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

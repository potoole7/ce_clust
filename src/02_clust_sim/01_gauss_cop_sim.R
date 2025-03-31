#### Simulation study for Gaussian copula data ####

# Assess performance of Jensen-Shannon clustering algorithm on simulation
# dataset generated from Gaussian copulas, for which we have certain
# theoretical guarentees
# Do so for a number of simulations and correlation parameters

# NOTE
# for df -> \infty, t-copula -> Gaussian copula, so alpha_t > alpha_gauss??

# Add start values to grid search (done)
# Add k estimation to grid search and rerun (but still use true k) (done)
# TODO Also save lines from plot, can combine for each simulation as before!
# Combine density plots pre and post clustering (done)
# Also produce boxplots for parameter estimates pre- and post-clustering (done)
# TODO Replace a and b with alpha and beta in plots
# How should I limit y-axes? Should I have them free or not?
# (done, let x-axes be free, y-axes fixed (to show difference in density)

# TODO Fix for varying cond_prob

#### libs ####

library(dplyr, quietly = TRUE)
devtools::load_all("../evc")
library(tidyr)
library(ggplot2)
library(patchwork)
library(parallel)

source("src/functions.R")

cores <- detectCores() - 1


#### Metadata ####

seed_number <- 123
n <- 1e4 # number of samples to take
n_locs <- 12 # number of "locations"
n_vars <- 2 # number of variables per location
mix_p <- c("gauss_cop" = 1, "t_cop" = 0) # mixture percentages (all Gaussian!)
# cond_prob <- 0.9 # quantile for conditional probability
cond_probs <- c(0.9, 0.95, 0.98, 0.99) # quantiles for conditional probability
max_clust <- n_locs - 1 # maximum number of clusters to choose from

# number of clusters
# n_clust <- 2 # number of clusters
n_clusts <- c(2, 3)
# vars <- c("rain", "wind_speed") # dummy names for variables

# fix b to a given value? Run for both
fixed_bs <- c(TRUE, FALSE)

stopifnot("Only supported for 2 or 3 clusters" = n_clusts %in% c(2, 3))

#### Functions ####

# pull parameters from dependence object
get_pars <- \(ce_fit) lapply(ce_fit, \(x) lapply(x, `[[`, "params"))

# function to find K thorugh combination of TWGSS, AIC and LRT
# TODO Allow use of cluster start values, only works with default for now
find_k <- \(ce_fit, data_mix, max_clust, fixed_b = FALSE) {
  # Pull JS distance matrix
  dist <- js_clust(ce_fit, scree_k = 1:max_clust)$dist_mat
  # extract TWGSS
  twgss <- evc::scree_plot(dist, k = 1:max_clust)

  # Calculate AIC and perform Likelihood Ratio test (i.e. model-based criteria)
  # initialise list of dependence models
  dep_clust_lst <- rep(
    get_pars(ce_fit),
    max_clust
  )
  # aic <- lr <- vector(mode = "list", length = max_clust - 1)
  aic <- lr <- vector(mode = "list", length = max_clust)
  # for (k in 2:max_clust) {
  for (k in 1:max_clust) {
    # cluster
    pam_fit <- js_clust(dist_mat = dist, k = k, return_dist = TRUE)
    # refit CE model
    dependence_clust <- fit_optim_clust(
      pam_fit$pam$clustering,
      data_mix,
      n_vars = n_vars,
      cond_prob = cond_prob,
      trans_fun = trans_fun,
      # TODO Fix so that we can use cluster start values
      # Could use correlation in data_mix sample before sampling??
      # start_vals = start_vals
      start_vals = c(0.01, 0.01),
      fixed_b = fixed_b
    )
    # extract parameters (remove residuals)
    dependence_clust <- get_pars(dependence_clust)
    dep_clust_lst[[k]] <- dependence_clust

    # calculate AIC
    aic[[k]] <- ce_aic(dependence_clust)

    # calculate likelihood ratio statistic
    # if (k > 2) {
    if (k > 1) {
      # df = 4 (n pars) * 2 (n_vars) (# more parameters in "full" model)
      df <- 4 * n_vars
      lr[[k]] <- lrt_test(
        dep_clust_lst[[k - 1]], dep_clust_lst[[k]],
        df = df
      )
    }
  }

  aic <- unlist(aic)
  # lr <- c(NA, NA, unlist(lapply(lr, "[[", "p_value")))
  lr <- c(NA, unlist(lapply(lr, "[[", "p_value")))

  # find best k as mode of choice for TWGSS, AIC and LRT
  get_mode <- function(x) {
    uniq_x <- unique(x)
    tab <- table(x)
    return(uniq_x[which.max(tab)])
  }
  # use elbow finding algorithm to choose k from TWGSS and AIC
  k_twgss <- find_elbow(twgss)
  k_aic <- find_elbow(aic)

  # find first LR test statistic w/ p-value < .05
  k_lr <- which(lr < 0.05)[1]

  # if no LR test statistic is significant, choose from elbow plots
  elbow_k <- get_mode(c(k_twgss, k_aic))
  # also do so if LRT chooses lowest k but other methods choose lower k,
  # as our LRT can only be computed for k > 2 (cannot cluster with k == 1)
  if (length(k_lr) == 0 || (k_lr == 2 && elbow_k < 2)) {
    k_lr <- elbow_k
  }

  # if (length(k_lr) == 0) k_lr <- 2
  ks <- c("k_twgss" = k_twgss, "k_aic" = k_aic, "k_lr" = k_lr)
  k <- get_mode(ks)

  # plot to confirm
  p <- data.frame("AIC" = -aic, "TWGSS" = twgss, k = seq_along(twgss)) |>
    pivot_longer(AIC:TWGSS) |>
    ggplot(aes(x = k, y = value, colour = name)) +
    geom_line() +
    geom_vline(xintercept = n_clust, linetype = "dashed") +
    facet_wrap(~name, scales = "free_y") +
    evc::evc_theme() +
    labs(
      title = paste0(
        "Elbow Plots for AIC and TWGSS, true k = ",
        n_clust,
        ", ",
        "LRT choice = ",
        k_lr,
        ", DQU = ",
        cond_prob * 100,
        "%"
      ),
      x = "k",
      y = ""
    ) +
    scale_x_continuous(breaks = 1:max_clust, limits = c(1, max_clust)) +
    guides(colour = "none") +
    ggsci::scale_colour_nejm()

  # return optimal k, k estimates from each method, and plot
  return(list(
    "k"        = k,
    "k_method" = ks,
    "plot"     = p
  ))
}


#### Test Example ####

# metadata for example:
# correlation for two clusters
# cor_gauss <- c(0.3, 0.9)
cond_prob <- 0.9 # conditional probability
# n_clust <- 2 # number of clusters
n_clust <- 3
cor_gauss <- seq(0.3, 0.9, length.out = n_clust) # correlation for clusters)
cluster_mem <- rep(1:n_clust, each = n_locs / n_clust) # cluster membership
fixed_b <- FALSE

# calculate theoretical alpha and beta values
start_vals <- lapply(cor_gauss, \(x) c(sign(x) * x^2, 1 / 2))

# TODO Can the same thing be achieved with qlaplace here and plaplace below??
set.seed(seed_number)
data <- sim_cop_dat(
  n_locs      = n_locs,
  n_vars      = n_vars,
  n           = n,
  n_clust     = n_clust,
  cluster_mem = cluster_mem,
  cor_gauss   = cor_gauss,
  mix_p       = c(1, 0), # only generate from Gaussian copula
  qfun        = qnorm
)
data_mix <- data$data_mix

# store original Gaussian data
data_mix_gauss <- data_mix
# transform normal data to CDF and convert to Laplace margins
trans_fun <- \(x) laplace_trans(pnorm(x))
data_mix <- lapply(data_mix, trans_fun)

# cluster with minimum and max rho
min_max_rho <- c(
  which(cluster_mem == min(cluster_mem))[1],
  which(cluster_mem == max(cluster_mem))[1]
)

# plot for locs in clusters 1 & 2, to show difference in sample correlations
# gaussian
plot(data_mix_gauss[[min_max_rho[1]]], xlab = "", ylab = "")
points(data_mix_gauss[[min_max_rho[[2]]]], col = "blue")
# laplace
plot(data_mix[[min_max_rho[1]]], xlab = "", ylab = "")
points(data_mix[[min_max_rho[2]]], col = "blue")

# look at chi and chibar for both
do.call(`+`, ggplot(texmex::chi(data_mix[[min_max_rho[1]]]), plot. = FALSE)) /
  do.call(`+`, ggplot(texmex::chi(data_mix[[min_max_rho[2]]]), plot. = FALSE))


# look at the correlation between variables in the extreme right tails
sapply(data_mix[c(min_max_rho)], \(x) {
  cor(x[x[, 1] > quantile(x[, 1], cond_prob) & x[, 2] > quantile(x[, 2], cond_prob), ])[[2]]
})

# fit CE to each location
ce_fit <- lapply(seq_along(data_mix), \(i) {
  o <- ce_optim(
    # Y         = x,
    Y = data_mix[[i]],
    dqu = cond_prob,
    control = list(maxit = 1e6),
    constrain = FALSE,
    # start = c(0.01, 0.01)
    start = start_vals[[cluster_mem[i]]],
    fixed_b = fixed_b # Should we fix b?
  )
})
# Returns list of len n_locs, w/ lists of len n_vars, of resids & params

# function to extract a and b values
extract_pars <- \(ce_fit) {
  x <- ce_fit[[1]]
  y <- ce_fit[[1]][[1]]
  bind_rows(lapply(ce_fit, \(x) { # loop through "locations"
    data.frame(
      vapply(x, \(y) { # loop through "variables"
        y$params[c("a", "b"), ]
      }, numeric(n_vars)),
      "parameter" = c("a", "b"),
      row.names = NULL
    )
  }), .id = "location") |>
    pivot_longer(var_1:var_2) |>
    mutate(location = as.numeric(location)) |>
    arrange(location, name, parameter)
}

pars_df <- extract_pars(ce_fit) |>
  # label cluster
  # mutate(cluster = ifelse(
  #   location %in% which(cluster_mem == 1), 1, 2
  # ))
  mutate(cluster = cluster_mem[location])

diff_from_start <- \(pars_df, start_vals, n = n_locs) {
  n_start_vals <- length(unlist(start_vals))

  pars_df |>
    mutate(
      # join start values
      # TODO See if this still works for two clusters!
      # start = unlist(lapply(start_vals, rep, times = n)),
      start = unlist(
        lapply(start_vals, rep, times = nrow(pars_df) / n_start_vals)
      ),
      # calculate difference
      diff = abs(value - start)
    ) |>
    group_by(cluster, parameter) |>
    summarise(
      across(c(value, start, diff), mean, .names = "{.col}_mean"),
      across(c(value, diff), sd, .names = "{.col}_sd"),
      .groups = "drop"
    ) |>
    relocate(
      cluster:parameter,
      starts_with("start"),
      starts_with("value"),
      starts_with("diff")
    )
}

diff_from_start(pars_df, start_vals)
# can clearly see uncertainty gets lower for higher rho values; makes sense
# as for these we have better convergence for bivariate normal to
# asymptotic (i.e. theoretical) a and b values
# In particular, beta -> 1/2 as rho -> 1, as expected

# Pull JS distance (for later clustering)
dist <- js_clust(ce_fit, scree_k = 1:max_clust)$dist_mat

# Find optimal k value using AIC, TWGSS and LRT
# TODO Investigate why k is 7 for LRT when k_true = 3!!
k_est <- find_k(ce_fit, data_mix_gauss, max_clust = 10, fixed_b = fixed_b)
(k_est$k_method) # see individual estimates
k_est$plot # plot looks great!!
k <- k_est$k
k == n_clust

# cluster for k = n_clust
pam_fit <- js_clust(
  dist_mat = dist, k = k, cluster_mem = cluster_mem
)
pam_fit$adj_rand == 1 # perfect clustering

clust_mem <- pam_fit$pam$clustering
# correct clustering if required
# if (clust_mem[1] == 2) clust_mem <- 3 - clust_mem
# TODO Fix this!!
# if (clust_mem[1] != 1) {
#   # redo so that clust mem 1 is 1, 2 is 2, etc
# }

ce_fit_clust <- fit_optim_clust(
  clust_mem = clust_mem,
  # data_mix, # should be data_mix_gaus, since we transform inside function!
  data_mix_gauss,
  n_vars,
  cond_prob,
  trans_fun = trans_fun,
  # start_vals = c(0.01, 0.01),
  start_vals = start_vals,
  fixed_b = fixed_b
)

# Slightly worse, but good estimates!
# TODO Need to fix beta = 1/2 and repeat this, see if it helps
pars_df_clust <- extract_pars(ce_fit_clust) |>
  rename(cluster = location)

# Compare estimates
diff_from_start(pars_df_clust, start_vals, n = 2)
# estimates ever so slightly worse but essentially the same
# however, uncertainty reduced naturally with more samples, as expected
# Interested to see below with higher conditional quantiles if the disparity
# pre- and post-clustering is greater


#### Grid Search ####

# grid of rho_gauss values to test
# grid <- tidyr::crossing(
#   cor_gauss1 = seq(0.1, 0.9, by = 0.2),
#   cor_gauss2 = seq(0, 0.8, by = 0.2),
#   cond_prob  = cond_probs # conditional probabilities to estimate CE model at
# )

# want following cases (rather than "exhaustive" grid search):
# rhos close together (low rho)
# rhos close together (medium rho)
# rhos close together (high rho)
# rhos far apart (one low, one med, one high)
grid <- data.frame(
  cor_gauss1 = c(0.1, 0.4, 0.7, 0.9),
  cor_gauss2 = c(0.2, 0.5, 0.8, 0.5),
  cor_gauss3 = c(0.3, 0.6, 0.9, 0.1)
) |>
  tidyr::crossing(
    "cond_prob" = cond_probs,
    "n_clust"   = n_clusts,
    "fixed_b"   = fixed_bs
  ) |>
  # remove redundant cor_gauss3 values
  mutate(cor_gauss3 = ifelse(n_clust == 2, NA, cor_gauss3))

# number of simulations to run for each
# n_times <- 2
n_times <- 500
# n_times <- 200
# objects for saving parameter values and cluster assignment
results_orig <- results_clust <- clust_assign <- vector(
  mode = "list", length = n_times
)
set.seed(seed_number)
results_grid <- mclapply(seq_len(nrow(grid)), \(i) {
  # results_grid <- lapply(seq_len(nrow(grid)), \(i) {
  print(paste0("Progress: ", round(i / nrow(grid), 3) * 100, "%"))
  system(sprintf(
    'echo "\n%s\n"',
    paste0(round(i / nrow(grid), 3) * 100, "% completed", collapse = "")
  ))

  print(paste0("i = ", i))

  # calculate theoretical alpha and beta values
  # start_vals <- lapply(grid[i, ], \(x) c("a" = sign(x) * x^2, "b" = 1 / 2))
  row <- grid[i, , drop = FALSE]
  if (row$n_clust == 2) {
    row <- select(row, -cor_gauss3)
  }

  # x clusters for 1/x of locations each
  cluster_mem <- rep(seq_len(row$n_clust), each = n_locs / row$n_clust)

  # start_vals <- lapply(grid[i, 1:2], \(x) c("a" = sign(x) * x^2, "b" = 1 / 2))
  start_vals <- lapply(
    row[, grepl("cor_gauss", names(row))], \(x) {
      c("a" = sign(x) * x^2, "b" = 1 / 2)
    }
  )

  for (j in seq_len(n_times)) {
    print(paste0("j = ", j))

    # generate simulation data
    data_mix <- sim_cop_dat(
      n_locs = n_locs,
      n_vars = n_vars,
      n = n,
      n_clust = row$n_clust,
      # cor_gauss = c(cor_gauss1, cor_gauss2),
      cor_gauss = unlist(row[, grepl("cor_gauss", names(row))]),
      mix_p = c(1, 0), # only generate from Gaussian copula
      qfun = qnorm,
      quiet = TRUE
    )$data_mix

    # transform to Laplace margins
    data_mix_gauss <- data_mix
    data_mix <- lapply(data_mix, trans_fun)

    # fit CE
    ce_fit <- lapply(seq_along(data_mix), \(l) {
      o <- ce_optim(
        Y = data_mix[[l]],
        dqu = row$cond_prob,
        control = list(maxit = 1e6),
        constrain = FALSE,
        start = start_vals[[cluster_mem[l]]], # start values for each cluster
        fixed_b = row$fixed_b
      )
    })

    # extract parameter estimates
    results_orig[[j]] <- extract_pars(ce_fit) |>
      # label cluster
      mutate(
        # cluster = ifelse(
        #   location %in% which(cluster_mem == 1), 1, 2
        # ),
        cluster = cluster_mem[location],
        sim = j # label simulation number
      )

    # Pull JS distance (for later clustering)
    dist <- js_clust(ce_fit, scree_k = 1:max_clust)$dist_mat

    # Find optimal k value using AIC, TWGSS and LRT
    # TODO Fix, currently causes an error
    k_estimate <- find_k(
      ce_fit, data_mix_gauss,
      max_clust = 10, fixed_b = fixed_b
    )

    # cluster and refit
    pam_fit <- js_clust(
      dist_mat = dist,
      k = row$n_clust, # NOTE Still use "true" k, easier to plot!!
      return_dist = TRUE
    )

    # correct clustering if required
    clust_mem <- pam_fit$pam$clustering

    # refit CE on clustered data
    ce_fit_clust <- fit_optim_clust(
      # pam_fit$pam$clustering,
      clust_mem,
      # data_mix,
      data_mix_gauss,
      n_vars = n_vars,
      cond_prob = row$cond_prob,
      trans_fun = trans_fun,
      start_vals = start_vals,
      fixed_b = row$fixed_b
    )
    # also save clustering object
    clust_assign[[j]] <- paste0(pam_fit$pam$clustering, collapse = "-")

    # extract parameters
    results_clust[[j]] <- extract_pars(ce_fit_clust) |>
      rename(cluster = location) |>
      mutate(
        # estimates for k (overall and for different methods)
        k_est    = k_estimate$k,
        k_aic    = k_estimate$k_method[["k_aic"]],
        k_twgss  = k_estimate$k_method[["k_twgss"]],
        k_lr     = k_estimate$k_method[["k_lr"]],
        sim      = j, # simulation number
        adj_rand = pam_fit$adj_rand # clustering performance
      )
  }
  # Return parameters pre and post-clustering
  return(list(
    "results_orig" = cbind(grid[i, ], bind_rows(results_orig)),
    "results_clust" = cbind(grid[i, ], bind_rows(results_clust)),
    "clust_assign" = cbind(
      grid[i, ], data.frame("clust" = unlist(clust_assign))
    )
  ))
}, mc.cores = cores)
# })

# save
file <- "data/sim_gauss_cop_grid.RDS"
# saveRDS(results_grid, file)

# load previous results
results_grid <- readRDS(file)


#### Analyse Simulations ###

# How to analyse?
# 1) Check for where clustering wasn't perfect (should be close cor_gauss vals)
# 2) Plot densities for par estimates pre- & post-clustering, with true vals
# How best to facet/plot?? might just need boxplots!!
# TODO Boxplots also!! Christian requested (densities are nice though hehehe)

# pull results
results_orig <- bind_rows(lapply(results_grid, `[[`, "results_orig")) |>
  mutate(cor_gauss3 = ifelse(n_clust == 2, NA, cor_gauss3))
results_clust <- bind_rows(lapply(results_grid, `[[`, "results_clust")) |>
  mutate(cor_gauss3 = ifelse(n_clust == 2, NA, cor_gauss3))
clust_assign <- bind_rows(lapply(results_grid, `[[`, "clust_assign")) |>
  mutate(cor_gauss3 = ifelse(n_clust == 2, NA, cor_gauss3))
# clust_assign <- unlist(lapply(results_grid, `[[`, "clust_assign"))

# check that clusters have been assigned properly
# TODO Check when clustering not assigned perfectly
head(sort(table(clust_assign$clust), decreasing = TRUE))
# all occur where correlations are very close together and <= 0.5, as expected
(clust_misclass <- clust_assign |>
  filter(!clust %in% c(
    "1-1-1-1-1-1-2-2-2-2-2-2",
    "1-1-1-1-2-2-2-2-3-3-3-3"
  )) |>
  count(across(contains("cor_gauss")), fixed_b, n_clust) |>
  arrange(n_clust, fixed_b))

# proportion of perfect clusterings
clust_assign |>
  filter(clust %in% c(
    "1-1-1-1-1-1-2-2-2-2-2-2",
    "1-1-1-1-2-2-2-2-3-3-3-3"
  )) |>
  nrow() %>%
  `/`(nrow(clust_assign))


# function to produce parsible title from grid search parameters
# TODO CHange from here to reflect addition of fixed_b and n_clust
lab_grid <- \(x) {
  x |>
    mutate(
      value_true = case_when(
        parameter == "a" & cluster == 1 ~ sign(cor_gauss1) * cor_gauss1^2,
        parameter == "a" & cluster == 2 ~ sign(cor_gauss2) * cor_gauss2^2,
        parameter == "a" & cluster == 3 ~ sign(cor_gauss3) * cor_gauss3^2,
        parameter == "b" ~ 1 / 2
      ),
      grid = paste0(
        "rho[Gauss[1]]==", round(cor_gauss1, 2),
        "~\",\"~",
        "rho[Gauss[2]]==", round(cor_gauss2, 2)
      )
    ) |>
    mutate(
      grid = ifelse(
        n_clust == 3,
        paste0(
          grid,
          "~\",\"~",
          "rho[Gauss[3]]==", round(cor_gauss3, 2)
        ),
        grid
      )
    )
}
# title_fun <- \(x) parse(text = unique(x$grid)) # parse title
# function to produce density plot/hist of parameter estimates w/ true vals
title_func <- \(y) {
  # parse(text = unique(y$grid))
  parse(text = paste0(
    unique(y$grid),
    ifelse(
      # all(unique(y$fixed_b) == TRUE), paste0("\", ~beta = 1/2"), ""
      all(unique(y$fixed_b) == TRUE), "~\",\"~beta==1/2", ""
    )
  ))
}
# dens_plot_fun <- \(x, title_fun = \(y) parse(text = unique(y$grid))) {
dens_plot_fun <- \(x, title_fun = title_func) {
  # remove bs if fixed, pointless plotting exercise
  if (all(x$fixed_b == TRUE)) {
    x <- x |>
      filter(parameter != "b")
  }

  p <- x |>
    ggplot(aes(x = value, fill = factor(cluster))) +
    geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.5) +
    geom_density(alpha = 0.5) +
    geom_vline(
      aes(xintercept = value_true, colour = factor(cluster)),
      linetype = "dashed",
      linewidth = 2
    ) +
    ggsci::scale_fill_nejm() +
    ggsci::scale_colour_nejm() +
    labs(
      # title = "Original",
      title = title_fun(x),
      # title =
      x = "Value",
      y = "Density"
    ) +
    guides(fill = "none", colour = "none") +
    evc::evc_theme()

  # facet by pre/post-clustering, if available (if not, just by parameter)
  if ("ind" %in% names(x)) {
    p <- p +
      facet_wrap(parameter ~ ind, scales = "free_x")
  } else {
    p <- p +
      facet_wrap(~parameter, scales = "free_y")
  }
  return(p)
}

# plot densities pre-clustering
results_orig_plt <- lab_grid(results_orig) |>
  # group_split(grid)
  # group_split(grid, fixed_b, n_clust, .keep = TRUE)
  group_split(grid, fixed_b, n_clust, cond_prob, .keep = TRUE)
p_orig <- lapply(results_orig_plt, dens_plot_fun)

# plot densities post-clustering
results_clust_plt <- lab_grid(results_clust) |>
  # group_split(grid)
  group_split(grid, fixed_b, n_clust, cond_prob)
p_clust <- lapply(results_clust_plt, dens_plot_fun)

# combine plots
results_combined_plt <- lapply(seq_along(results_orig_plt), \(i) {
  results_orig_plt[[i]] |>
    mutate(ind = "Pre-clustering") |>
    bind_rows(
      results_clust_plt[[i]] |>
        mutate(ind = "Post-clustering")
    ) |>
    mutate(
      ind = factor(ind, levels = c("Pre-clustering", "Post-clustering"))
    )
})
p_combined <- lapply(results_combined_plt, dens_plot_fun)

plot_file <- "plots/01_sim/00_gauss_cop_densities.pdf"
if (fixed_b) plot_file <- gsub(".pdf", "_fixed_b.pdf", plot_file)
pdf(plot_file, width = 8, height = 5)
p_combined
dev.off()

# also want to produce boxplots
boxplot_fun <- \(x, title_fun = title_func) {
  if (all(x$fixed_b == TRUE)) {
    x <- x |>
      filter(parameter != "b")
  }

  p <- x |>
    ggplot(aes(x = factor(cluster), y = value, fill = factor(cluster))) +
    geom_boxplot() +
    # add true values
    geom_hline(
      aes(yintercept = value_true, colour = factor(cluster)),
      linetype = "dashed",
      linewidth = 2
    ) +
    ggsci::scale_fill_nejm() +
    ggsci::scale_colour_nejm() +
    labs(
      title = title_fun(x),
      x = "Clusters",
      y = "Value"
    ) +
    guides(fill = "none", colour = "none") +
    evc::evc_theme() +
    # remove axis text and labels
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank()
    )

  if ("ind" %in% names(x)) {
    p <- p +
      # facet_grid(ind ~ parameter, scales = "free")
      # facet_wrap(ind ~ parameter, scales = "free")
      facet_wrap(parameter ~ ind, scales = "free_x")
  } else {
    p <- p +
      facet_wrap(~parameter, scales = "free")
  }
  return(p)
}

# plot boxplots for pre- and post-clustering results
# TODO Decide whether to keep on same y-scale or not?
p_box_combined <- lapply(results_combined_plt, boxplot_fun)

boxplot_file <- "plots/01_sim/00_gauss_cop_boxplots.pdf"
if (fixed_b) boxplot_file <- gsub(".pdf", "_fixed_b.pdf", boxplot_file)
pdf(boxplot_file, width = 8, height = 5)
p_box_combined
dev.off()

# Also plot boxplots with multiple different configurations on the one plot
results_orig_boxplt <- lab_grid(results_orig) |>
  # group_split(fixed_b, n_clust, cond_prob)
  group_split(fixed_b, n_clust, grid)

results_clust_boxplt <- lab_grid(results_clust) |>
  # group_split(fixed_b, n_clust, cond_prob)
  group_split(fixed_b, n_clust, grid)

# TODO Latex font doesn't match eslewhere, so annoying!
# TODO Add vertical lines for true values, rather than x's
# TODO Functionalise!
p_box_full <- lapply(seq_along(results_orig_boxplt), \(i) {
  results_orig_boxplt[[i]] |>
    mutate(ind = "Pre-clustering") |>
    bind_rows(
      results_clust_boxplt[[i]] |>
        mutate(ind = "Post-clustering")
    ) |>
    mutate(
      ind = factor(ind, levels = c("Pre-clustering", "Post-clustering"))
    )
}) |>
  lapply(\(x) {
    if (all(x$fixed_b == TRUE)) {
      x <- x |>
        filter(parameter != "b")
    }

    x_plot <- x |>
      # mutate(grid = paste0("rho[gauss] == c(", cor_gauss1, ", ", cor_gauss2, ")")) |>
      mutate(
        grid = ifelse(
          n_clust == 3,
          paste0(
            "c(", cor_gauss1, ", ", cor_gauss2, ", ", cor_gauss3, ")"
          ),
          paste0("c(", cor_gauss1, ", ", cor_gauss2, ")")
        ),
        cond_prob = factor(paste0(cond_prob * 100, "%"))
      ) |>
      # only show lower and upper dependence quantile, to highlight differences
      filter(cond_prob %in% c("90%", "99%"))

    p <- x_plot |>
      ggplot(aes(x = cond_prob, y = value, fill = factor(cluster))) +
      geom_boxplot(outliers = FALSE) +
      # add true values (as x's)
      # geom_point(
      #   # aes(x = grid, y = value_true),
      #   aes(x = cond_prob, y = value_true),
      #   colour = "black",
      #   size = 4,
      #   pch = 4,
      #   stroke = 2,
      #   position = position_dodge(width = 0.8) # match position of boxes
      # ) +
      # add true values for alphas
      geom_hline(
        data = filter(x_plot, parameter == "a"),
        aes(yintercept = value_true, colour = factor(cluster)),
        linetype = "dashed",
        linewidth = 1,
        alpha = 0.7
      ) +
      geom_hline(
        data = filter(x_plot, parameter == "b"),
        aes(yintercept = value_true),
        linetype = "dashed",
        colour = "black",
        linewidth = 1,
        alpha = 0.7
      ) +
      # add true value for beta (same, at 1/2)
      # scale_x_discrete(labels = function(x) parse(text = x)) +
      # labs(x = parse(text = "rho[Gauss]")) +
      labs(
        title = parse(text = unique(x$grid)),
        x = "Dependence Quantile",
        y = "Parameter value"
      ) +
      evc::evc_theme() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      guides(fill = "none", colour = "none") +
      ggsci::scale_fill_nejm() +
      ggsci::scale_colour_nejm()

    if (all(x$fixed_b == TRUE)) {
      # p <- p + labs(title = parse(text = "beta==1/2"))
      p <- p +
        labs(title = parse(text = paste0(unique(x$grid), "~\",\"~beta==1/2")))
    }

    if ("ind" %in% names(x)) {
      p <- p +
        # facet_wrap(parameter ~ ind)
        facet_grid(
          rows = vars(parameter),
          cols = vars(ind),
          scales = "free_y",
          labeller = labeller(parameter = c("a" = "alpha", "b" = "beta"), .default = label_parsed)
        )
    } else {
      p <- p +
        facet_wrap(~parameter, scales = "fixed_x")
    }
    p
  })

boxplot_full_file <- "plots/01_sim/00_gauss_cop_boxplots_full.pdf"
pdf(boxplot_full_file, width = 8, height = 5)
p_box_full
dev.off()


# Notes:
# For varying b:
# - Seems like pre- and post-clustering, alpha and beta estimates are very similar,
# but variance is reduced, naturally as sample size increases considerably
# - In general, alpha is quite well estimated, but beta is not as well estimated
# - Beta is better estimated for higher rho_gauss values, while for lower ones
# alpha (i.e. the relationship) is quite small and thus it might be hard to
# estimate beta?
# As highlighted above, what would happen if we fixed beta??
# - Does lots better! But still not great for smaller alpha/independence
# - Particularly does better post- clustering
# - However, density is slightly less "peaked", as model might be exploring
# more of the parameter space for a?

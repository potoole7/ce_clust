#### Choose the best k for PAM clustering algorithm ####

# TODO Change kl_prob as well! Want to see how well algorithm does for diff dqu
# TODO Why doesn't silhouette work? Need different definition

# TODO Save plot(s) as 01_e ...

#### Libs ####

library(dplyr, quietly = TRUE)
# library(evc)
devtools::load_all("../evc")
library(tidyr)
library(ggplot2)
library(evd)
library(evgam)
# devtools::load_all("texmex")
library(cluster)
library(lcmix)
library(copula)
library(parallel)

# source functions
source("src/functions.R")


#### Metadata ####

seed_number <- 123
n_locs <- 60 # number of "locations" (needs got be higher for more clusters)
n_clust <- 11
# cluster_mem <- ceiling(1:n_locs / (n_locs / n_clust))
n_vars <- 2 # number of variables per location
mix_p <- c("gauss_cop" = 0.5, "t_cop" = 0.5) # mixture percentages
# GPD parameters
scale_gpd <- 1
shape_gpd <- -0.05
marg_prob <- 0.9 # marginal threshold quantile
# kl_prob <- 0.9 # conditional threshold quantile
kl_prob <- seq(0.85, 0.975, by = 0.025) # conditional threshold quantiles
# n <- 1e3 # number of samples to take
n_exceed <- 100 # number of exceedances to have for each conditional threshold
n <- round(n_exceed / (1 - kl_prob))
conf_level <- 0.95 # confidence level for CIs in plot
# try values of k from 2 to max_clust
max_clust <- 15

# dataframe with number of sims required for each dqu to have 100 exceedances
exceed_df <- data.frame("kl_prob" = kl_prob, "n" = n)

# Number of cores to use for parallel computation
n_cores <- detectCores() - 1

# Change name of saved files depending on options chosen above
save_file <- paste0("data/js_best_k_clust_", n_clust)
plot_file <- paste0("plots/js_best_k_clust_", n_clust)
if (length(kl_prob) > 1) {
  save_file <- paste0(save_file, "_dqu")
  plot_file <- paste0(plot_file, "_dqu")
}
# if (length(n) > 1) {
#   save_file <- paste0(save_file, "_eq_n_exceed")
#   plot_file <- paste0(plot_file, "_eq_n_exceed")
# }
save_file <- paste0(save_file, ".RDS")
plot_file <- paste0(plot_file, ".pdf")

#### Functions ####
# TODO Move to functions file and/or `evc` package
# TODO Ensure functions are "pure"

# choose optimal k (from http://sherrytowers.com/2013/10/24/k-means-clustering/)
# opt_k <- \(x) {
#   v = -diff(x)
#   nv = length(v)
#   fom = v[1:(nv-1)]/v[2:nv]
#   return(which.max(fom) + 1)
# }

# calculate AIC for a given CE fit
# TODO Add BIC
calc_inf <- \(ll, n_par, n = NULL, type = "AIC") {
  return(switch(type,
    "AIC" = (-2 * ll) + (2 * n_par) # ,
    # "BIC"  = (-2 * ll) + (log(n) * n_par)
  ))
}

# Function to calculate AIC from CE model output
ce_aic <- \(dependence) {
  # pull LLs for each model
  lls <- vapply(dependence, \(x) {
    vapply(x, `[`, "ll", , FUN.VALUE = numeric(1))
  }, numeric(length(dependence[[1]])))
  # number of parameters = 4 * number of clusters * number of variables
  # TODO Will this number of larger for > 2 variables? for 3 it'll be 6!
  # TODO Test for 3 variables before putting in package or anything
  n_par <- 4 * length(dependence) * length(dependence[[1]])
  # n par for different AICs per model
  # n_par <- 4 * length(dependence[[1]])
  # combined_ll
  ll <- sum(lls)
  # ll <- mean(lls)
  # return AIC
  return(calc_inf(ll, n_par))
  # calculate AIC for each model
  # return(vapply(lls, \(x) {
  #   calc_inf(x, n_par)
  # }, numeric(1)))
}

# fun to refit CE model with given clustering
fit_ce_clust <- \(clust_mem, data_mix) {
  clusts <- unique(clust_mem)
  # create new data_mix with cluster membership
  data_mix_clust <- lapply(clusts, \(x) {
    do.call(rbind, data_mix[clust_mem == x])
  })

  # refit model and return
  return(fit_ce(
    data_mix_clust,
    vars = paste0("col_", seq_len(n_vars)),
    # TODO Functionalise these inputs, may want to vary kl_prob!
    marg_prob = marg_prob,
    cond_prob = kl_prob,
    # f           = NULL, # fit models with ismev
    fit_no_keef = TRUE,
    output_all = FALSE
  ))
}


#### Grid search ####

# TODO: Write note to explain grid setup
# setup grid
grid <- bind_rows(lapply(2:n_clust, \(i) {
  # data.frame(t(c(seq(0, 1, length.out = i), rep(NA, times = n_clust - i))))
  data.frame(t(c(seq(0, 0.95, length.out = i), rep(NA, times = n_clust - i))))
})) %>%
  setNames(paste0("cor_t", seq_len(n_clust))) %>%
  mutate(
    # true number of clusters
    k_true = seq_len(nrow(.)) + 1,
    # cor_gauss = 0.1, # 0.3, # keep constant
    # Degrees of freedom for t-copula
    df_t   = 3,
    # extremal quantiles
    # kl_prob = kl_prob,
    # mixture percentages (must sum to 1)
    mix_p  = 0.5
  ) %>%
  # extremal quantiles and sample sizes
  crossing(kl_prob = kl_prob) %>%
  left_join(exceed_df, by = "kl_prob")
# crossing inverts grid for some reason, undo
grid <- grid[nrow(grid):1, ]


# run kl_sim_eval for each row in grid
# TODO Increase n_times to 500
n_times <- 100
# initialise vectors to store results for each repetition of simulations
# elb_vec <- sil_vec <- aic_vec <- vector(length = n_times)
results_vec <- mem_vec <- vector(length = n_times)
elb_vec <- sil_vec <- aic_vec <- vector(length = n_times)
# i <- 3
# i <- 13
set.seed(seed_number)
results_grid <- bind_rows(mclapply(seq_len(nrow(grid)), \(i) {
  # results_grid <- bind_rows(lapply(seq_len(nrow(grid)), \(i) {

  print(paste0("Progress: ", round(i / nrow(grid), 3) * 100, "%"))
  system(sprintf(
    'echo "\n%s\n"',
    paste0(round(i / nrow(grid), 3) * 100, "% completed", collapse = "")
  ))

  # pull specific row
  row <- grid[i, , drop = FALSE]
  print(row)
  # pull correlations in row which are non NA
  cor_t <- row[, grepl("cor_t", names(row))]
  cor_t <- gauss_t <- cor_t[!is.na(cor_t)]
  n_cor <- length(cor_t)
  # cluster membership known to split locations evenly
  cluster_mem <- ceiling(seq_len(n_locs) / (n_locs / n_cor))
  # if (i == 3) {
  #   # cor_t <- c(0, 0.4, 0.6, 0.8)
  #   cor_t <- c(0, 0.4, 0.7, 0.9)
  # }
  for (j in seq_len(n_times)) {
    # generate simulation data for given parameter set
    data_mix <- with(row, sim_cop_dat(
      n           = n,
      n_locs      = n_locs,
      n_clust     = n_cor,
      # cor_gauss   = rep(cor_gauss, n_cor),
      # Important: set cor_gauss same as cor_t, to have more pronounced diffs
      cor_gauss   = cor_t,
      cor_t       = cor_t,
      df_t        = rep(df_t, n_cor),
      params_gpd  = c(scale_gpd, shape_gpd),
      mix_p       = c(0.5, 0.5),
      cluster_mem = cluster_mem
    ))$data_mix

    # see how different they are (works for 4 clusters)
    # plot(data_mix[[1]])
    # plot(data_mix[[46]])
    # lapply(data_mix[c(1, 46)], as.data.frame) %>%
    #   bind_rows(.id = "location") %>%
    #   setNames(c("location", "x", "y")) %>%
    #   ggplot(aes(x = x, y = y, colour = location)) +
    #   geom_density_2d() +
    #   evc::evc_theme()

    # Fit CE model, if an error is produced, return a dummy list and skip
    dependence <- tryCatch(
      {
        fit_ce(
          data_mix,
          vars = paste0("col_", seq_len(n_vars)),
          marg_prob = marg_prob,
          cond_prob = row$kl_prob,
          # f           = NULL, # fit models with ismev
          fit_no_keef = TRUE,
          output_all = FALSE
        )
      },
      error = function(cond) {
        # TODO: Fill in what goes here
        return(1)
      }
    )
    # If "true" k was known, how well would clustering do?
    clust_res <- js_clust(dependence, k = n_cor, cluster_mem = cluster_mem)
    results_vec[[j]] <- clust_res$adj_rand
    # store clustering membership (true values known as cluster_mem)
    mem_vec[[j]] <- NA
    if (clust_res$adj_rand < 1) {
      mem_vec[[j]] <- paste0(clust_res$pam$clustering, collapse = "_")
    }

    # Pull JS distance matrix
    # TODO: Don't return plot
    dist <- js_clust(dependence, scree_k = 1:max_clust)$dist_mat

    # check for three clusters that distances are distinct
    # plot(dist)
    # summary(as.matrix(dist)[, 1][2:20])
    # summary(as.matrix(dist)[, 1][21:40])
    # summary(as.matrix(dist)[, 1][41:60])

    # plot(dist)
    # image(as.matrix(dist))
    # image_lst[[i]] <- image(as.matrix(dist))
    # summary(as.matrix(dist)[, 1][2:15])
    # summary(as.matrix(dist)[, 1][16:30])
    # summary(as.matrix(dist)[, 1][31:45])
    # summary(as.matrix(dist)[, 1][46:60])

    # Calculate TWGSS for elbow
    # TODO: Don't return plot
    twgss <- evc::scree_plot(dist, k = 1:max_clust)
    # TODO Just return pasted twgss
    # TODO Find better method of "algorithmically" choosing k from elbow
    # elb_vec[[j]] <- opt_k(twgss)
    elb_vec[[j]] <- paste0(round(twgss, 2), collapse = "_")

    # NOTE Silhouette is not correctly implemented for this method; ignore
    # # Calculate silhouette coefficient for 2-max_clust clusters
    # sil_boxplot(dist, k = 2:max_clust, show_plt = TRUE)$plot
    # sil <- sil_boxplot(dist, k = 2:max_clust, show_plt = FALSE)$sil
    # # Find silhouette coefficient for different clusterings
    # sil_vec[[j]] <- sil %>%
    #   group_by(k) %>%
    #   summarise(mean = mean(sil_width, na.rm = TRUE), .groups = "drop") %>%
    #   filter(mean == max(mean)) %>%
    #   pull(k)

    # cluster for different values of k, refit CE model and calculate AIC
    aic <- vapply(2:max_clust, \(k) {
      # aic <- lapply(2:max_clust, \(k) {
      # cluster
      pam_fit <- js_clust(
        dependence,
        dist_mat = dist, k = k, return_dist = TRUE
      )
      # refit CE model
      # TODO Investigate no exceedances message from this! May just be wrong
      dependence_clust <- fit_ce_clust(pam_fit$pam$clustering, data_mix)
      # calculate AIC
      ce_aic(dependence_clust)
    }, numeric(1))

    # aic_vec[[j]] <- opt_k(aic)
    aic_vec[[j]] <- paste0(round(aic, 0), collapse = "_")

    # # plot for multiple aic values per k
    # bind_rows(
    #   lapply(aic, \(x) data.frame("aic" = x)),
    #   .id = "k"
    # ) %>%
    #   mutate(k = as.numeric(k) + 1) %>%
    #   group_by(k) %>%
    #   mutate(row = row_number()) %>%
    #   ungroup() %>%
    #   ggplot() +
    #   # geom_boxplot(aes(x = k, y = -aic, group = k)) +
    #   geom_line(aes(x = row, y = -aic, colour = factor(k)), size = 2) +
    #   evc::evc_theme()
    #
    # vapply(aic, sum, numeric(1)) %>%
    #   as.data.frame() %>%
    #   setNames("aic") %>%
    #   mutate(k = 2:max_clust) %>%
    #   ggplot(aes(x = k, y = -aic)) +
    #   geom_point() +
    #   geom_line()

    # plot(-aic, type = "b")
  }

  # return(cbind(row, "sil_k" = sil_vec, "elb_k" = elb_vec, "aic_k" = aic_vec))
  # TODO Add local rand index
  return(cbind(
    row,
    "adj_rand" = results_vec,
    "membership" = mem_vec,
    # "k_true"   = n_cor, # true number of clusters
    "elb" = elb_vec,
    "aic" = aic_vec
  ))
}, mc.cores = n_cores))

# remove NAs, give warning
nas <- is.na(results_grid$adj_rand)
if (sum(nas) > 0) {
  message(paste0("Warning: ", sum(nas), " NAs in results, removing"))
  results_grid <- results_grid[!nas, ]
}

# save
saveRDS(results_grid, file = save_file)
results_grid <- readRDS(save_file)



#### Plotting ####

# For each k val, plot elbow plots on top of each other with some alpha
# Like PPD plots in Bayesian stats

# test plot: plot first simulation
k <- 4
par(mfrow = c(1, 2))
plot(
  as.numeric(stringr::str_split(
    results_grid[results_grid$k_true == k, ]$elb[1], "_"
  )[[1]]),
  type = "b", ylab = paste0("TWGSS (k_true = ", k, ")")
)
plot(
  -as.numeric(stringr::str_split(
    results_grid[results_grid$k_true == k, ]$aic[1], "_"
  )[[1]]),
  type = "b", ylab = paste0("-AIC (k_true = ", k, ")")
)
par(mfrow = c(1, 1))

# plot AIC and TWGSS elbows across all simulations and "true" k values
k_true_vals <- unique(results_grid$k_true)

# fun to pull out elbow and AIC values for a given number of clusters
# TODO Move above
pull_elb_aic <- \(x) {
  bind_rows(lapply(seq_len(nrow(x)), \(i) {
    # elbow/TWGSS (stored in character string seperated by underscores)
    elb <- as.numeric(stringr::str_split(x$elb[i], "_")[[1]])
    # AIC
    aic <- c(NA, as.numeric(stringr::str_split(x$aic[i], "_")[[1]]))
    return(cbind(
      dplyr::select(x[i, , drop = FALSE], -c(elb, aic)),
      "elb" = elb,
      "aic" = aic,
      "k" = seq_along(elb),
      row.names = NULL
    ))
  }))
}

# extract AIC and TWGSS values across all simulations for all true k
# aic_elb_vals_lst <- mclapply(k_true_vals, \(k) {
#   results_grid %>%
#     filter(k_true == k) %>%
#     # simulation number, indicating each unique simulation performed
#     mutate(sim_n = row_number()) %>%
#     pull_elb_aic()
# }, mc.cores = n_cores)
aic_elb_vals_lst <- results_grid %>%
  group_split(k_true, kl_prob, .keep = TRUE) %>%
  mclapply(\(df) {
    df %>%
      # simulation number, indicating each unique simulation performed
      mutate(sim_n = row_number()) %>%
      pull_elb_aic()
  }, mc.cores = n_cores)

# plot for each true value of k
# TODO: Functionalise? Pretty nice plot!
plots <- mclapply(aic_elb_vals_lst, \(x) {
  k_true <- x$k_true[1]
  x %>%
    pivot_longer(cols = c("elb", "aic")) %>%
    mutate(
      value = ifelse(name == "aic", -value, value),
      name = ifelse(name == "elb", "TWGSS", "-AIC")
    ) %>%
    ggplot(aes(x = k, y = value, colour = name, group = sim_n)) +
    geom_line(alpha = 0.2) +
    geom_vline(xintercept = k_true, linetype = "dashed") +
    facet_wrap(~name, scales = "free_y") +
    evc::evc_theme() +
    labs(
      title = paste0(
        "Elbow Plots for AIC and TWGSS, true k = ",
        k_true,
        ", ",
        "DQU = ",
        x$kl_prob[1] * 100,
        "%, ",
        "n_sim = ",
        n_times
      ),
      x = "k",
      y = ""
    ) +
    scale_x_continuous(breaks = 1:max_clust, limits = c(1, max_clust)) +
    guides(colour = "none") +
    ggsci::scale_colour_nejm()
}, mc.cores = n_cores)

pdf(file = plot_file, width = 10, height = 10)
plots
dev.off()

# TODO: Visualise Rand Index (get from 03b_js_sens.R) (needed??)

# TODO: Visualise Local Rand Index (needed??)

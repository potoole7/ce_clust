#### Bootstrap uncertainty in CE estimates pre- and post-clustering ####

# Want to see if our clustering improves the conditional extremes parameter
# uncertainty by bootstrapping the estimates from our simulation example

# NOTE Plot objects for paper are saved below if you just want to manually
# save them to avoid clipping LaTeX

# What will plot be in paper???

#### libs ####

devtools::load_all("../evc")
library(dplyr)
library(boot)
library(parallel)
library(ggplot2)
library(tidyr)

# source simulation and bootstrapping functions
source("src/functions.R")

ncores <- detectCores() - 1


#### Metadata ####

seed_number <- 123
n <- 1e3 # number of samples to take
n_locs <- 12 # number of "locations"
n_clust <- 2 # number of clusters
cluster_mem <- sort(rep(1:n_clust, n_locs / n_clust)) # two clusters for half of locations each
n_vars <- 2 # number of variables per location
mix_p <- c("gauss_cop" = 0.5, "t_cop" = 0.5) # mixture percentages
# t-copula df
df_t <- 3
# Gaussian and t-copula correlation parameters
# Choose nice values where we know the model should work well
cor_gauss_ex <- 0.5
cor_t_ex <- c(0.1, 0.9)
# GPD parameters
scale_gpd <- 1
shape_gpd <- -0.05
# Marginal and conditional thresholds before and after clustering
cond_prob <- 0.9
cond_prob_clust <- 0.9
marg_prob <- 0.98 # marginal quantile to calculate expectation in bootstrap
# R <- 100
R <- 500 # number of bootstrap samples

# parameters for MC estimation of JSGa
n_mc <- 500 # number of MC samples
mc_method <- "laplace_trunc2" # truncate Laplace using empirical dist quant
laplace_cap <- 0.99 # Take 99th quantile

# Number of cores to use for parallel computation
n_cores <- detectCores() - 1


#### Bootstrap alpha and beta parameter values before clustering ####

# generate simulation data where we know model works well
set.seed(seed_number)
data_mix <- sim_cop_dat(
  n = n,
  cor_gauss = c(cor_gauss_ex, cor_gauss_ex),
  cor_t = cor_t_ex,
  df_t = c(df_t, df_t),
  mix_p = c(0.5, 0.5),
  qfun = evd::qgpd,
  qargs = c("loc" = 0, "scale" = scale_gpd, "shape" = shape_gpd)
)$data_mix

# fit CE model to save plots of data pre and post-transformation
ce <- fit_ce(
  data_mix,
  marg_prob = 0.95,
  cond_prob = cond_prob,
  fit_no_keef = TRUE,
  ncores = n_cores,
  output_all = TRUE
)

# plot original vs transformed data
par(mfrow = c(1, 2))
# with(ce$original[[4]], plot(rain, wind_speed, xlab = "", ylab = ""))
with(ce$original[[4]], plot(var1, var2, xlab = "", ylab = ""))
plot(ce$transformed[[1]], xlab = "", ylab = "")
par(mfrow = c(1, 1))

plt_orig <- ce$original[[4]] |>
  select(-name) |>
  mutate(
    ind = "Original",
    col_var = case_when(
      # wind_speed == max(wind_speed) ~ 1, # case with weak dependence
      var2 == max(var2) ~ 1, # case with weak dependence
      # wind_speed > 5 & rain > 6 ~ 2, # case with strong dependence
      var2 > 5 & var1 > 6 ~ 2, # case with strong dependence
      # case with strong negative dependence
      TRUE ~ 4
    ),
    row_n = row_number()
  ) |>
  mutate(col_var = ifelse(
    row_n == which(ce$transformed[[4]][, 1] < -6 & ce$transformed[[4]][, 1] < -5),
    3,
    col_var
  )) |>
  # pivot_longer(rain:wind_speed)
  pivot_longer(var1:var2)

plt_df <- bind_rows(
  plt_orig,
  data.frame(
    # "rain"       = ce$transformed[[4]][, 1],
    "var1" = ce$transformed[[4]][, 1],
    # "wind_speed" = ce$transformed[[4]][, 2],
    "var2" = ce$transformed[[4]][, 2],
    "ind" = "Laplace",
    "row_n" = seq_len(nrow(ce$transformed[[4]]))
  ) |>
    # pivot_longer(rain:wind_speed) |>
    pivot_longer(var1:var2) |>
    mutate(col_var = plt_orig$col_var)
) |>
  pivot_wider(names_from = name, values_from = value) |>
  mutate(ind = factor(ind, levels = c("Original", "Laplace")))

p <- plt_df |>
  # ggplot(aes(x = rain, y = wind_speed)) +
  ggplot(aes(x = var1, y = var2)) +
  # geom_point(aes(colour = factor(col_var), size = factor(col_var)), alpha = 0.7) +
  geom_point() +
  geom_point(
    data = plt_df[plt_df$col_var != 4, ],
    aes(colour = factor(col_var)),
    size = 5
  ) +
  facet_wrap(~ind, scales = "free") +
  evc_theme() +
  labs(x = "", y = "") +
  scale_colour_manual(values = ggsci::pal_nejm()(3)) +
  guides(colour = "none") +
  NULL

# TODO Save again if fixable?
# ggsave("latex/plots/transform.png", p, width = 8, height = 6)

# Function to transform data, fit CE and perform bootstrapping
fit_boot <- \(data_mix, ncores = 1, ret_ce = FALSE) {
  # transform to Laplace margins
  data_mix_trans <- lapply(data_mix, \(x) {
    laplace_trans(evd::pgpd(x, loc = 0, scale = scale_gpd, shape = shape_gpd))
  })

  # fit CE
  ce_fit <- lapply(seq_along(data_mix_trans), \(i) {
    o <- ce_optim(
      Y = data_mix_trans[[i]],
      dqu = cond_prob,
      control = list(maxit = 1e6),
      constrain = FALSE
    )
  })

  # extract dependence parameters
  dep <- lapply(ce_fit, \(x) {
    lapply(x, `[[`, "params")
  })

  # run bootstrap
  set.seed(seed_number)
  bootstrap_res <- boot_ce_dep(
    dep,
    data_mix,
    data_mix_trans,
    marg_pars = list(loc = 0, scale = scale_gpd, shape = shape_gpd),
    cond_prob = cond_prob,
    marg_prob = 0.98,
    R = R,
    trace = 10,
    ncores = ncores
  )
  if (ret_ce) {
    # return(list(bootstrap_res, ce_fit))
    # return(list(bootstrap_res, ce_fit, data_mix_trans))
    return(list(
      "bootstrap_res" = bootstrap_res,
      "ce_fit" = ce_fit,
      "data_mix_trans" = data_mix_trans
    ))
  }
  return(bootstrap_res)
}

bootstrap_run <- fit_boot(data_mix, ncores = ncores, ret_ce = TRUE)
ce_fit <- bootstrap_run$ce_fit
data_mix_trans <- bootstrap_run$data_mix_trans
bootstrap_res <- bootstrap_run$bootstrap_res

# extract the bootstrapped data for each cluster
boot_est <- purrr::transpose(lapply(bootstrap_res, \(x) { # loop through bootstrap samples
  lapply(seq_len(n_clust), \(i) {
    unname(unlist(x[cluster_mem == i]))
  })
}))

# convert to dataframe to plot
boot_est_df <- bind_cols(lapply(boot_est, unlist)) |>
  setNames(c(1, 2)) |>
  pivot_longer(everything(), names_to = "cluster") |>
  mutate(cluster = as.numeric(cluster)) |>
  mutate(cluster = ifelse(cluster == 1, cor_t_ex[1], cor_t_ex[2]))

# TODO Return to add post-clustering bootstrap results
plot_boot <- \(boot_est_df) {
  boot_est_df |>
    mutate(cluster = factor(cluster, levels = cor_t_ex)) |>
    ggplot(aes(x = value)) +
    geom_histogram(aes(y = ..density.., fill = cluster), bins = 30, alpha = 0.7) +
    geom_density(aes(fill = cluster), alpha = 0.5) +
    facet_wrap(
      ~cluster,
      scales = "fixed",
      labeller = labeller(cluster = \(x) {
        return(ifelse(x == 0.1, "rho[t] == 0.1", "rho[t] == 0.9"))
      }, .default = label_parsed)
    ) +
    labs(x = "", y = "") +
    evc::evc_theme() +
    guides(fill = "none") +
    ggsci::scale_fill_nejm()
}

plot_boot(boot_est_df)


#### Bootstrap **after** clustering ####

# cluster for k = 2
k <- n_clust
clust <- js_clust(
  ce_fit,
  trans = data_mix_trans,
  k = k, cluster_mem = cluster_mem,
  n = n_mc, mc_method = mc_method, laplace_cap = laplace_cap,
  return_dist = TRUE
)
# check clustering perfect
clust$adj_rand == 1

# assign data based on clustering results to cluster centroids
data_mix_clust <- lapply(seq_len(k), \(x) {
  do.call(rbind, data_mix[clust$pam$clustering == x])
})

# repeat bootstrapping, this time on clustered observations
bootstrap_res_clust <- fit_boot(data_mix_clust, ncores = ncores)

# extract the bootstrapped data for each cluster
boot_est_clust <- purrr::transpose(lapply(bootstrap_res_clust, \(x) {
  lapply(x, \(y) unname(unlist(y)))
}))

boot_est_clust_df <- bind_cols(lapply(boot_est_clust, unlist)) |>
  setNames(c(1, 2)) |>
  pivot_longer(everything(), names_to = "cluster") |>
  mutate(cluster = as.numeric(cluster)) |>
  mutate(cluster = ifelse(cluster == 1, cor_t_ex[1], cor_t_ex[2]))

# plot just the clustered data
plot_boot(boot_est_clust_df)

p_dat <- bind_rows(
  boot_est_df,
  boot_est_clust_df,
  .id = "ind" # pre vs post-clustering
) |>
  mutate(
    cluster = factor(cluster, levels = cor_t_ex),
    ind = factor(case_when(
      ind == 1 ~ "Pre-clustering",
      TRUE ~ "Post-clustering"
    ), levels = c("Pre-clustering", "Post-clustering"))
  )
# TODO Investigate why histogram for pre-clustering is higher than density?
p <- p_dat |>
  ggplot(aes(x = value)) +
  # geom_histogram(aes(y = ..density.., fill = ind), bins = 30, alpha = 0.5) +
  # geom_histogram(aes(fill = ind), bins = 30, alpha = 0.5) +
  geom_density(aes(y = ..density.., fill = ind), alpha = 0.5) +
  facet_wrap(
    ~cluster,
    scales = "free_y",
    labeller = labeller(cluster = \(x) {
      return(ifelse(x == 0.1, "rho[t] == 0.1", "rho[t] == 0.9"))
    }, .default = label_parsed)
  ) +
  labs(x = "Conditional Expectation", y = "", fill = "") +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  evc::evc_theme()
p

saveRDS(p, "plots/sim_01e_bootstrap.rds") # save plot object
ggsave("latex/plots/sim_01e_bootstrap.png", p, width = 8, height = 5)

# also do boxplot
p_box <- p_dat |>
  ggplot(aes(x = cluster, y = value, fill = ind)) +
  geom_boxplot(width = 1, outliers = FALSE, key_glyph = "rect") +
  labs(x = parse(text = "rho[t]"), y = "Conditional Expectation", fill = "") +
  evc::evc_theme() +
  ggsci::scale_fill_nejm() +
  scale_x_discrete(expand = c(0.25, 0.25))
p_box

saveRDS(p_box, "plots/sim_01e_bootstrap_box.rds")
ggsave("latex/plots/sim_01e_bootstrap_box.png", p_box, width = 8, height = 5)

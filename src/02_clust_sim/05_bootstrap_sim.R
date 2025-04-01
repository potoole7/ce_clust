#### Bootstrap uncertainty in CE estimates pre- and post-clustering ####

# Want to see if our clustering improves the conditional extremes parameter
# uncertainty by bootstrapping the estimates from our simulation example

# What will plot be in paper???

#### libs ####

devtools::load_all("../evc")
library(dplyr)
library(boot)
library(parallel)
library(ggplot2)
library(tidyr)

# source simulation functions
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
    return(list(bootstrap_res, ce_fit))
  }
  return(bootstrap_res)
}

bootstrap_res <- fit_boot(data_mix, ncores = ncores, ret_ce = TRUE)
ce_fit <- bootstrap_res[[2]]
bootstrap_res <- bootstrap_res[[1]]

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
clust <- js_clust(ce_fit, k = k, cluster_mem = cluster_mem)
# check clustering perfect
clust$adj_rand == 1

# assign data based on clustering results to cluster centroids
data_mix_clust <- lapply(seq_len(k), \(x) {
  do.call(rbind, data_mix[clust$pam$clustering == x])
})

# repeat bootstrapping
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
p <- p_dat |>
  ggplot(aes(x = value)) +
  # geom_histogram(aes(y = ..density.., fill = ind), bins = 30, alpha = 0.5) +
  geom_density(aes(fill = ind), alpha = 0.5) +
  facet_wrap(
    ~cluster,
    scales = "free_y",
    labeller = labeller(cluster = \(x) {
      return(ifelse(x == 0.1, "rho[t] == 0.1", "rho[t] == 0.9"))
    }, .default = label_parsed)
  ) +
  labs(x = "Conditional Expectation", y = "", fill = "") +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  evc::evc_theme() +
  ggsci::scale_fill_nejm()
p

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
ggsave("latex/plots/sim_01e_bootstrap_box.png", p_box, width = 8, height = 5)

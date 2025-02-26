#### See how uncertainty changes with bootstrapping ####

# TODO: Need to use `sim` = "balanced"

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
cluster_mem <- sort(rep(1:2, n_locs / 2)) # two clusters for half of locations each
n_vars <- 2 # number of variables per location
mix_p <- c("gauss_cop" = 0.5, "t_cop" = 0.5) # mixture percentages
# GPD parameters
scale_gpd <- 1
shape_gpd <- -0.05
# marg_prob <- 0.9
marg_prob <- 0.95
# cond_prob <- 0.9
cond_prob <- 0.95
conf_level <- 0.95 # confidence level for CIs in plot

# Number of cores to use for parallel computation
n_cores <- detectCores() - 1 

#### Test ####

# generate simulation data where we know model works well
set.seed(seed_number)
data_mix <- sim_cop_dat(
  n          = n,
  # cor_gauss  = c(0.7, 0.7),
  cor_gauss  = c(0.4, 0.4),
  cor_t      = c(0.1, 0.9),
  df_t       = c(3, 3),
  params_gpd = c(scale_gpd, shape_gpd),
  mix_p      = c(0.5, 0.5)
)$data_mix

# fit CE
ce_fit <- fit_ce(
  data_mix, 
  marg_prob   = marg_prob,
  cond_prob   = cond_prob,
  fit_no_keef = TRUE, 
  output_all  = TRUE
)

# run bootstrap
bootstrap_res <- boot_ce(
  ce_fit, 
  R = 100, 
  trace = 10, 
  ncores = ncores
)

# plot margins for both theoretical clusters
plot_marg <- \(x) {
  x |> 
    ggplot(aes(x = value, fill = parameter)) +
    # ggplot(aes(x = value, fill = factor(cluster))) +
    # also add histograms
    geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.7) +
    geom_density(alpha = 0.5) +
    # facet_wrap(~ parameter + vars, scales = "free", nrow = 2) +
    facet_wrap(~ cluster + parameter + vars, scales = "free", nrow = 2) +
    # facet_grid(vars ~ parameter, scales = "free") +
    evc::evc_theme() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
}
bootstrap_res$marginal |> 
  mutate(
    cluster = ifelse(name %in% paste0("location_", 1:6), "cor_t = 0.1", "cor_t = 0.9")
  ) |> 
  plot_marg()
  
# do the same with dependence
# TODO Very interesting that some results for no corr give 1 and -1! Say to C
plot_dep <- \(x) {
  x |> 
    ggplot(aes(x = value, fill = parameter)) +
    geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.7) +
    geom_density(alpha = 0.5) +
    facet_wrap(~ cluster + parameter + vars, scales = "free", nrow = 2) +
    scale_x_continuous(limits = c(-1, 1)) +
    evc::evc_theme() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
}
bootstrap_res$dependence |> 
  mutate(
    cluster = ifelse(name %in% paste0("location_", 1:6), "cor_t = 0.1", "cor_t = 0.9"), 
    vars = paste0(vars, " | ", cond_var)
  ) |> 
  plot_dep()
  

#### Bootstrap on clustered observations ####

# cluster for k = 2
k <- 2
clust <- js_clust(ce_fit$dependence, k = k, cluster_mem = cluster_mem)
# check clustering perfect
clust$adj_rand == 1

# assign data based on clustering results to cluster centroids
data_mix_clust <- lapply(seq_len(k), \(x) {
  do.call(rbind, data_mix[clust$pam$clustering == x])
})

# TODO repeat bootstrapping and parameter extraction process
evc_fit_clust <- fit_ce(
  data_mix_clust, 
  marg_prob   = marg_prob,
  cond_prob   = cond_prob,
  fit_no_keef = TRUE, 
  output_all  = TRUE
)
bootstrap_res_clust <- boot_ce(
  evc_fit_clust, 
  R = 100, 
  trace = 10, 
  ncores = ncores
)

bootstrap_res_clust$marginal |> 
  mutate(
    cluster = ifelse(name == "location_1", "clust_1", "clust_2")
  ) |> 
  plot_marg()

# dramatically better! But median of densities seems pretty similar?!
plt_dep_data_clust <- bootstrap_res_clust$dependence |>
  mutate(
    cluster = ifelse(name == "location_1", "clust_1", "clust_2"), 
    vars = paste0(vars, " | ", cond_var)
  )
plot_dep(plt_dep_data_clust)

# yep, alpha values very similar!!
# Better when using cond_prob = 0.95 rather than 0.9
plt_dep_data_clust |> 
  group_by(cluster, vars, parameter) |> 
  summarise(value = median(value, na.rm = TRUE), .groups = "drop") |> 
  pivot_wider(names_from = cluster, values_from = value) |> 
  arrange(parameter, vars)

par(mfrow = c(1, 2))
plot(evc_fit_clust$transformed$location_1, main = "cor_t = 0.1")
plot(evc_fit_clust$transformed$location_2, main = "cor_t = 0.9")
par(mfrow = c(1, 1))

# TODO plot dependence parameter confidence intervals for both
conf_level <- 0.95
z <- qnorm((1 + conf_level) / 2)


#### Grid search ####

# TODO Do for a number of simulations!! 
# Needed?? Maybe just doing this once is sufficient? Show Christian anyway!
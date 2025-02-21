#### See how uncertainty changes with bootstrapping ####

#### libs ####

devtools::load_all("../evc")
library(dplyr)
library(boot)
library(parallel)

# source simulation functions
source("src/functions.R")

#### Functions ####

# TODO Would be much easier if fit_ce returns an object with methods!!
# TODO Look at texmex package for implementation
# TODO Ensure extreme observations in sample, sometimes not (see texmex)
# TODO Change pipes
boot_ce <- \(data, indices) {
  
  # subset and split by location
  # TODO Anyway to take same subset of indices for each location?
  data_subset <- data[indices, ] %>% 
    group_split(name) %>% 
    # convert back to matrix
    lapply(\(x) as.matrix(x[, -ncol(x)]))
  
  # parameters of interest
  pars <- c("a", "b")
  n_par <- length(pars) # n_pars
  n_par_per_var <- n_par * ncol(data_subset[[1]]) # n_pars * n_vars
  
  # fit CE
  # TODO Where do names "rain" and "wind speed" come from???
  # TODO Add tryCatch for when fit_ce fails
  dependence <- fit_ce(
    data_subset, 
    # marg_prob   = 0.9,
    marg_prob = list(
      f = list("response ~ name", "~ name"), tau = 0.9, jitter = TRUE
    ),
    cond_prob   = 0.9,
    fit_no_keef = TRUE
  )
  
  # extract a and b parameter values
  # TODO Add method in package
  param_vals <- as.vector(vapply(dependence, \(x) {
    as.vector(vapply(x, \(y) {
      y[pars, ]
    }, numeric(n_par)))
  }, numeric(n_par_per_var)))
  
  # label appropriately
  names(param_vals) <- paste(
    # location names
    # names(dependence), 
    rep(names(dependence), each = n_par_per_var),
    # parameter names pasted to variable names
    as.vector(outer(pars, names(dependence[[1]]), paste, sep = "_")), 
    sep = "_"
  )
  return(param_vals)
}


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
marg_prob <- 0.9
kl_prob <- 0.9
conf_level <- 0.95 # confidence level for CIs in plot

# Number of cores to use for parallel computation
n_cores <- detectCores() - 1 

#### Test ####

# generate simulation data where we know model works well
set.seed(seed_number)
data_mix <- sim_cop_dat(
  n          = n,
  cor_gauss  = c(0.7, 0.7),
  cor_t      = c(0.1, 0.9),
  df_t       = c(3, 3),
  params_gpd = c(scale_gpd, shape_gpd),
  mix_p      = c(0.5, 0.5)
)$data_mix

# convert data to dataframe as required by boot (split within boot_ce, by loc)
data_mix_df = bind_rows(lapply(seq_along(data_mix), \(i) {
  ret <- data.frame(data_mix[[i]])
  ret$name <- paste0("location_", i)
  return(ret)
}))

# bootstrap CE fit to this data
# TODO Look into the different options for boot
# R = 100: ~64.7s, R = 500: ~ 349s (6 mins), R = 1000 ~ 680s (11 mins) (linear)
debugonce(boot_ce)
boot_res <- boot(
  data = data_mix_df, 
  statistic = boot_ce, 
  R = 100,
  # R = 1000,
  parallel = "multicore", 
  ncpus = n_cores
)

# Extract parameter estimates and related bca confidence intervals
param_est <- boot_res$t0

# TODO Look into confidence interval for first parameter
# TODO Look into why bca not working (may just be too difficult!)
(ci1 <- boot.ci(boot_res, type = c("norm", "perc", "basic"), index = 1))




# extract confidence intervals
conf_int <- mclapply(seq_along(param_est), function(i) {
  # Get the BCa confidence intervals for the i-th coefficient
  # TODO Parallelise
  # ci <- boot.ci(boot_res, type = "norm", index = i)
  # ci <- boot.ci(boot_res, type = "perc", index = i)
  ci <- boot.ci(boot_res, type = c("norm", "perc"), index = i)
  
  # Return the lower and upper bounds for both confidence intervals
  return(c(ci$normal[4:5], ci$percent[4:5]))
})


# Now, fit CE to this data and cluster (k known to be 2)
dependence <- fit_ce(
  data_mix, 
  marg_prob   = 0.9,
  cond_prob   = 0.9,
  fit_no_keef = TRUE
)
clust <- js_clust(dependence, k = 2, cluster_mem = cluster_mem)
# check clustering perfect
clust$adj_rand == 1

# TODO assign data based on clustering results to cluster centroids
# data_mix_clust...

# TODO repeat bootstrapping and parameter extraction process

# TODO Compare results (does clustering improve parameter estimates?)


#### Grid search ####

# TODO Do for a number of simulations!! (Maybe show Christian results of this first?)
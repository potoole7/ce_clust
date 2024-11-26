#### Test evc package functions to see if they're working ####

# Rough test of package functions

#### Libs ####

library(dplyr)
# library(evc)
devtools::load_all("../evc")

#### Load Data ####

data_mix <- readRDS("data/sim_dat_3_clust.RDS")
cluster_mem <- c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3)

#### Vignotto 2021 ####

kl_sim_eval(
  data_mix, 
  kl_prob     = 0.9, 
  k           = 3, 
  cluster_mem = cluster_mem 
)

#### Fit CE model and cluster ####

dependence <- fit_ce(
  data_mix, 
  # marg_prob   = 0.9,
  marg_prob = list(
    f      = list("response ~ name", "~ name"), 
    tau    = .95, 
    jitter = TRUE
  ),
  cond_prob   = 0.9,
  split_data  = TRUE 
)

js_clust(
  dependence,  
  nclust = 3, 
  cluster_mem = cluster_mem
)

#### Plotting ####

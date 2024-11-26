#### Apply Jensen-Shannon Divergence method to simulated data ####

#### Libs ####

library(dplyr, quietly = TRUE)
library(evc)
library(tidyr)
library(ggplot2)
library(cluster)
library(evgam)
devtools::load_all("texmex")

source("src/functions.R")

theme <- theme_bw() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 13, face = "bold"),
    strip.background = element_rect(fill = NA, colour = "black"),
    plot.tag = element_text(size = 16, face = "bold"),
    panel.background = element_rect(fill = NA, colour = "black")
  )


#### Metadata ####

vars <- c("rain", "wind_speed")
prob <- 0.9
cluster_mem <- sort(rep(c(1, 2, 3), 4))

#### Load Data ####

data_mix <- readRDS("data/sim_dat_3_clust.RDS")

# convert to dataframe
data_df <- bind_rows(lapply(seq_along(data_mix), \(i) {
  n <- length(data_mix[[i]])
  data.frame(
      "rain"       = data_mix[[i]][1:(n / 2)],
      "wind_speed" = data_mix[[i]][((n / 2) + 1):n]
  ) %>% 
    mutate(name = paste0("location_", i))
}))


#### Calculate Conditional Extremes parameters ####

dependence <- fit_ce(
  data_df, 
  marg_prob   = prob, 
  cond_prob   = prob, 
  split_data  = TRUE, 
  fit_no_keef = TRUE
)

# check that all dependence models have run successfully
sapply(dependence, \(x) lapply(x, length))


#### Calculate divergence upon which to cluster ####

# Perform PAM and k-means clustering, as an exploratory analysis

# first, elbow plot
js_mat <- js_clust(dependence)$dist_mat # suggests k = 3

# cluster and assess performance
js_clust(dependence, k = 3, dist_mat = js_mat, cluster_mem = cluster_mem)

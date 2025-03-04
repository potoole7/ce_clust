#### Conditional Extremes clustering for Irish data ####

#### Libs ####

library(sf)
# library(evc)
devtools::load_all("../evc")
library(dplyr, quietly = TRUE)
library(tidyr)
library(ggplot2)
# # devtools::load_all("texmex") # forked texmex which allows for fixed b vals
library(circular)

source("src/functions.R")

sf::sf_use_s2(FALSE)

#### Metadata ####

max_k <- 10 # maximum number of clusters to look at 
# dqu <- 0.7
dqu <- 0.9

#### Load Data ####

# load original dataset
data <- readr::read_csv("data/met_eireann/final/met_eir_preprocess.csv.gz")

# pull county for locations
locs <- distinct(data, name, county)

# load shapefile
areas <- sf::read_sf("data/met_eireann/final/irl_shapefile.geojson")

# Fitted conditional extremes models for each site
dependence <- readRDS("data/texmex_mexdep_obj.RDS")
# sort alphabetically by name
dependence <- dependence[sort(names(dependence))]
# dependence <- readRDS("data/texmex_mexdep_obj_fixed_b.RDS")

# a and b values
ab_df <- readr::read_csv("data/ab_vals.csv") %>% 
  arrange(name)
# ab_df <- readr::read_csv("data/ab_vals_fixed_b.csv")

# first check names for dependence are the same as in ab_df
all(unique(ab_df$name) == names(dependence))

# remove locations with NAs for a
na_locs <- ab_df %>% 
  filter(is.na(value)) %>% 
  pull(name) %>% 
  unique()

if (length(na_locs) > 0) {
  ab_df <- filter(ab_df, !name %in% na_locs)
  dependence <- dependence[
    !names(dependence) %in% na_locs
  ]
}

# extract point location of each station for plotting on map
pts <- ab_df %>% 
  distinct(name, lon, lat) %>% 
  st_to_sf()

# calculate adjacency matrix from Voronoi cells
adj_mat <- calc_adj_mat(pts, plot = TRUE)

# remove locations with NAs from data
data <- data %>% 
  semi_join(ab_df, by = "name")

# convert ab values to wide format
ab_df_wide <- ab_df %>% 
  mutate(col = paste0(var, "_", parameter)) %>% 
  dplyr::select(name, county, col, value) %>% 
  pivot_wider(names_from = col, values_from = value)


#### Choose k ####

# scree plot (and distance matrix)
set.seed(123)
clust_obj <- js_clust(dependence, scree_k = 1:max_k)
dist_mat <- clust_obj$dist_mat
twgss <- clust_obj$total_within_ss
# no clear elbow unfortunately

# AIC
# TODO Put in package
# calculate AIC for a given CE fit
calc_inf <- \(ll, n_par, n = NULL, type = "AIC") {
  return(switch(
    type, 
    "AIC"  = (-2 * ll) + (2 * n_par)
  ))
}

# Function to calculate AIC from CE model output
# TODO Put in package
ce_aic <- \(dependence) {
  # pull LLs for each model
  lls <- vapply(dependence, \(x) {
    vapply(x, `[`, "ll", , FUN.VALUE = numeric(1))
  }, numeric(length(dependence[[1]])))
  # number of parameters = 4 * number of clusters * number of variables
  n_par <- 4 * length(dependence) * length(dependence[[1]])
  # n par for different AICs per model
  # combined_ll
  ll <- sum(lls)
  # return AIC
  return(calc_inf(ll, n_par))
}

# cluster for different values of k, refit CE model and calculate AIC
# TODO specific to this analysis, move to top of script
calc_aic <- \(dist_mat, k_vals) {
  aic <- vapply(k_vals, \(k) {
    # cluster
    pam_fit <- js_clust(
      dependence, dist_mat = dist_mat, k = k # , return_dist = TRUE
    )
    
    # pull medoids and their location information
    medoids_df <- data %>% 
      filter(name %in% pam_fit$medoids) %>% 
      distinct(across(name:lat))
    
    # add data to new clusters
    data_clust <- data %>% 
      left_join(
        # match locations to their cluster centroids
        data.frame(
          "name"   = names(pam_fit$clustering), 
          "medoid" = pam_fit$medoids[pam_fit$clustering])
      ) %>% 
      # replace location information with medoid information
      select(-names(medoids_df)) %>% 
      left_join(medoids_df, by = c("medoid" = "name")) %>% 
      relocate(medoid:lat, .after = date) %>% 
      rename(name = medoid) %>% 
      arrange(name, desc(date))
    
    # evgam formula will change depending on number of clusters
    # default
    f_marg_prob <- list("response ~ name", "~ name")
    f <- list(excess ~ name, ~1)
    if (k > 3) {
        # f = list("response ~ s(lon, lat, k = 50)", "~ s(lon, lat, k = 40)")
      f_marg_prob <- list(
        # "response ~ s(lon, lat, k = 50)", 
        paste0("response ~ s(lon, lat, k = ", k, ")"), # 50)", 
        # "~ s(lon, lat, k = 40)"
        paste0("~ s(lon, lat, k = ", k, ")")
      )
      f <-  list(
        # excess ~ s(lon, lat, k = 40), 
        paste0("excess ~ s(lon, lat, k = ", k, ")"),
        # ~s(lon, lat, k = 40)
        paste0("~s(lon, lat, k = ", k, ")")
      )
    }
    
    # refit CE model
    ce_fit <- fit_ce(
      data = data_clust,
      vars = c("rain", "wind_speed"), 
      # arguments to `quantile_thresh`
      marg_prob = list(
        # f = list("response ~ s(lon, lat, k = 50)", "~ s(lon, lat, k = 40)")
        f = f_marg_prob
      ), 
      cond_prob = dqu, # dqu
      # formula for evgam fit to excesses over threshold
      # f = list(excess ~ s(lon, lat, k = 40), ~s(lon, lat, k = 40)), 
      f = f,
      fit_no_keef = FALSE, # don't fit without Keef constraints
      ncores = max(1, parallel::detectCores() - 2),
      output_all = TRUE
    )
  
    # calculate AIC
    ce_aic(ce_fit$dependence)
  }, numeric(1))
  
  return(aic)
}

aic <- calc_aic(dist_mat, 2:max_k)

# plot aic and twgss together nicely
elb_plt <-\(twgss, aic) {
  data.frame("elb" = twgss, "aic" = c(NA, aic), "k" = 1:max_k) %>% 
    pivot_longer(cols = c("elb", "aic")) %>%
    mutate(
      value = ifelse(name == "aic", -value, value),
      name = ifelse(name == "elb", "TWGSS", "-AIC")
    ) %>% 
    ggplot(aes(x = k, y = value, colour = name)) + #, group = sim_n)) + 
    # geom_line(alpha = 0.2) + 
    geom_line() + 
    # geom_vline(xintercept = k_true, linetype = "dashed") +
    facet_wrap(~name, scales = "free_y") + 
    evc::evc_theme() + 
    labs(
      title = "Elbow Plots for AIC, TWGSS",
      x     = "k", 
      y     = ""
    ) + 
    scale_x_continuous(breaks = 1:max_k, limits = c(1, max_k)) + 
    guides(colour = "none") + 
    ggsci::scale_colour_nejm()
}

elb_plt(twgss, aic)
# AIC seems to suggest k = 6, while TWGSS has no clear elbow 

# Cluster for reasonable k, save results for later refitting of CE
pam_clust <- lapply(2:max_k, \(k) { 
  js_clust(dependence, k = k, dist_mat = dist_mat)
})
saveRDS(pam_clust, "data/pam_clust.RDS")

# save map plots to see which seems most sensible
pam_map_p_lst <- lapply(seq_along(pam_clust), \(k) {
  plt_clust_map(pts, areas, pam_clust[[k]]) + 
    ggtitle(paste0("k" = k + 1))
})

pdf("plots/pam_map_p_lst.pdf")
pam_map_p_lst
dev.off()

# Conclusion: Unsure! k = 3 or 6 probably best, but no clear elbow for TWGSS

#### Clustering ####

# Cluster for "correct" k
pam_clust3 <- js_clust(dependence, k = 3, dist_mat = dist_mat)
pam_clust6 <- js_clust(dependence, k = 6, dist_mat = dist_mat)

# save clustering
saveRDS(pam_clust3, "data/pam_clust3.RDS")
saveRDS(pam_clust6, "data/pam_clust6.RDS")

# plot vs covariates
covar_plt(pts, pam_clust3, data, areas)
covar_plt(pts, pam_clust6, data, areas)

#### Adjacency stipulation clustering ####

# repeat for adjacent sites only
dist_mat_adj <- as.matrix(dist_mat)
dist_mat_adj[adj_mat == 0] <- 1e9
dist_mat_adj <- as.dist(dist_mat_adj)

# calculate TWGSS and AIC
clust_obj_adj <- js_clust(
  dependence, dist_mat = dist_mat_adj, scree_k = 1:max_k
)
twgss_adj <- clust_obj_adj$total_within_ss
aic_adj <- calc_aic(dist_mat_adj, 2:max_k)

elb_plt(twgss_adj, aic_adj) # perhaps k = 8?

# cluster and save data and plots
pam_clust_adj <- lapply(2:max_k, \(k) { 
  js_clust(dependence, k = k, dist_mat = dist_mat_adj)
})
saveRDS(pam_clust_adj, "data/pam_clust_adj.RDS")

pam_map_p_lst_adj <- lapply(seq_along(pam_clust_adj), \(k) {
  plt_clust_map(pts, areas, pam_clust_adj[[k]]) + 
    ggtitle(paste0("k" = k + 1))
})

pdf("plots/pam_map_p_lst_adj.pdf")
pam_map_p_lst_adj
dev.off()

# Plot for optimal k value
pam_clust_adj8 <- js_clust(dependence, k = 8, dist_mat = dist_mat_adj)
plt_clust_map(pts, areas, pam_clust_adj8)

# plot vs covariates (alt, dist2coast, wind_dir)
covar_plt(pts, pam_clust_adj8, data, areas)


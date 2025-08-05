#### Testing truncated Laplace on results of CE model for Ireland ####

# What data do we need?
# - CE fits for each site
# - areas, station points and elevations for plotting

# Look at:
# - Clustering solutions
# - 2-norm, mean, max and condition number for dissimilarity matrices at each

# What order to select different quantiles etc? Proposal:
# - CE bootstrapping decides dependence quantile (0.85)
# - Laplace truncation quantile possibly decided by looking at divergence
# matrices, or clustering results
# - finally, k chosen based on simple scree plot for favoured previous choices

# TODO Add elevation to plots
# TODO Look at each variable separately, as well as together!
# TODO Maybe try for different dependence quantiles also??

#### Libs ####

library(sf)
# library(evc)
devtools::load_all("../evc")
# devtools::load_all("../evc_mc")
library(dplyr, quietly = TRUE)
library(tidyr)
library(ggplot2)
library(lubridate)
library(terra)
library(elevatr)
# library(texmex)
# devtools::load_all("texmex") # forked texmex which allows for fixed b vals
library(gridExtra)
library(patchwork)
library(geosphere)
library(latex2exp)
library(ggpattern)
library(parallel)

source("src/functions.R")

sf::sf_use_s2(FALSE)

theme <- ggplot2::theme_bw() +
  ggplot2::theme(
    legend.position  = "bottom",
    plot.title       = ggplot2::element_text(size = 16, hjust = 0.5),
    axis.text        = ggplot2::element_text(size = 12),
    axis.title       = ggplot2::element_text(size = 14, face = "bold"),
    legend.text      = ggplot2::element_text(size = 12),
    strip.text       = ggplot2::element_text(size = 13, face = "bold"),
    strip.background = ggplot2::element_rect(fill = NA, colour = "black"),
    plot.tag         = ggplot2::element_text(size = 16, face = "bold"),
    panel.background = ggplot2::element_rect(fill = NA, colour = "black")
  )

#### Metadata ####

# number of MC samples to take
n_mc <- 1000
# n_mc <- 200

k <- 3 # number of clusters

# dqu <- c(0.85, 0.88, 0.9)
# dqu <- 0.85 # happiest with this quantile from bootstrapping
# TODO Expand to more Laplace quantiles?
laplace_q <- seq(0.7, 0.95, by = 0.05)
# laplace_q <- c(0.8, 0.85, 0.9)

# grid <- crossing(
#   dqu = dqu,
#   laplace_q = laplace_q
# )


#### Load Data ####

# load Conditional Extremes fits
ce_fit <- readRDS("data/ce_fit.RDS")
dep_fit <- readRDS("data/dep_fit.RDS")

# load shapefile site locations and elevation raster for plotting
areas <- sf::read_sf("data/met_eireann/final/irl_shapefile.geojson")
pts <- sf::read_sf("data/met_eireann/final/irl_points.geojson")
pts <- pts |>
  filter(name %in% names(dep_fit$dependence))
elev_df <- readr::read_csv(
  "data/met_eireann/final/irl_elev.csv",
  show_col_types = FALSE
) |>
  mutate(elev_bin = factor(elev_bin, levels = unique(elev_bin)))


#### Cluster for different truncation points ####

set.seed(123)
# x <- laplace_q[[1]]
# TODO Split by variable!!
clust_obj_lst <- lapply(laplace_q, \(x) {
  # First, calculate divergence matrix and TWGSS for scree plot
  dist_obj <- js_clust(
    dep_fit$dependence,
    k = NULL,
    scree_k = 1:5,
    n = n_mc,
    # mc_method = "uniform"
    mc_method = "laplace_trunc",
    laplace_cap = x,
    return_dist = TRUE
  )

  # next, cluster
  # TODO (for k = 3, presumed to always be best but needs checking!)
  clust_obj <- js_clust(
    dist_mat = dist_obj$dist_mat,
    dep_fit$dependence,
    k = k
  )

  # return distance matrices for each variable, as well as clustering + TWGSS
  return(list(
    "dist_mats" = dist_obj$dist_mats,
    "clust_obj" = clust_obj,
    "twgss"     = dist_obj$total_within_ss
  ))
})

# quickly plot TWGSS elbow for each (k = 3 best for each!)
bind_rows(lapply(seq_along(clust_obj_lst), \(i) data.frame(
  "laplace_q" = laplace_q[i],
  "twgss"     = clust_obj_lst[[i]]$twgss
))) |>
  group_by(laplace_q) %>%
  mutate(row = row_number()) %>%
  ungroup() %>%
  ggplot() +
  geom_line(aes(x = row, y = twgss)) +
  geom_point(aes(x = row, y = twgss)) +
  facet_wrap(~laplace_q, scales = "fixed") +
  theme

#### Compare divergences ####

# extract divergences
div_lst <- lapply(clust_obj_lst, `[[`, "dist_mats")

# transpose from laplace_q -> variable list to vice versa
div_lst <- purrr::transpose(div_lst)

# take first variable as test
div_lst_spec <- div_lst[[1]]

# calculate mean, max, 2-norm and condition number for each
div_mean <- sapply(div_lst_spec, \(x) mean(x))
div_max <- sapply(div_lst_spec, \(x) max(x))
div_norm <- sapply(div_lst_spec, \(x) norm(as.matrix(x), type = "2"))
div_kappa <- sapply(div_lst_spec, \(x) {
  kappa(as.matrix(x), exact = TRUE, triangular = TRUE)
})

# combine into dataframe
div_df <- data.frame(
  "laplace_q" = laplace_q,
  "mean"      = div_mean,
  "max"       = div_max,
  "norm"      = div_norm,
  "kappa"     = div_kappa
)

# plot raw values
div_df |>
  pivot_longer(cols = -laplace_q) |>
  ggplot(aes(x = laplace_q, y = value)) +
  geom_line() +
  geom_point() +
  facet_wrap(~name, scales = "free_y") +
  theme
# comments:

# plot first derivative
div_df |>
  mutate(
    across(c(mean, max, norm, kappa), \(x) {
      c(NA, diff(x) / diff(laplace_q))
    })
  ) |>
  pivot_longer(cols = -laplace_q) |>
  ggplot(aes(x = laplace_q, y = value)) +
  geom_line() +
  geom_point() +
  facet_wrap(~name, scales = "free_y") +
  theme
# comments: First big bump occurs at 0.85, so 0.8 may be a good candidate?
# lines up well with simulation results!


#### Plot clustering solutions on map ####

# plot each clustering solution
plt_lst <- lapply(seq_along(clust_obj_lst), \(i) {
  plt_clust_map(
    pts,
    areas,
    clust_obj_lst[[i]]$clust_obj,
    elev_df = elev_df,
    rm_elev_leg = TRUE
  ) +
    ggtitle(
      label = paste0("Laplace truncation at ", laplace_q[i])
    )
})

# save to inspect
pdf(
  file = "laplace_trunc_test.pdf",
  width = 12, height = 8
)
plt_lst
dev.off()
# conclusion: They all look pretty good except 0.7 and 0.95, others are largely
# similar, 0.75 and 0.9 have nice median sites if that's important!

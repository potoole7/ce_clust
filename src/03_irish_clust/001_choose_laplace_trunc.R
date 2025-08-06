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

# TODO:
# Add elevation to plots (done)
# TODO Comment on condition number and 2-norm of dissimilarity matrices
# Look at each variable (rain and wind) separately, as well as together! (done)
# TODO Allow for different Laplace caps for each variable?? "Best" seems
# different for rain and wind speed, somewhat unsurprisingly
# TODO Maybe try different dependence quantiles also?? Load from other file

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

vars <- c("rain", "wind_speed")

# number of MC samples to take
n_mc <- 1000
# n_mc <- 200

k <- 3 # number of clusters

# dqu <- c(0.85, 0.88, 0.9)
# dqu <- 0.85 # happiest with this quantile from bootstrapping
# TODO Expand to more Laplace quantiles?
# laplace_q <- seq(0.7, 0.95, by = 0.05)
laplace_q <- seq(0.65, 0.95, by = 0.05)
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
x <- laplace_q[[1]]
# TODO Parallelise if memory allows
clust_obj_lst <- lapply(laplace_q, \(x) {
  print(x)
  # First, calculate divergence matrix and TWGSS for scree plot
  dist_obj <- js_clust(
    dep_fit$dependence,
    k = NULL,
    # scree_k = 1:5,
    scree_k = list(1:5, 1:5), # create scree plot for each variable
    n = n_mc,
    # mc_method = "uniform"
    mc_method = "laplace_trunc",
    laplace_cap = x,
    return_dist = TRUE
  )

  # # next, cluster
  # # TODO (for k = 3, presumed to always be best but needs checking!)
  # clust_obj <- js_clust(
  #   dist_mat = dist_obj$dist_mat,
  #   dep_fit$dependence,
  #   k = k
  # )

  # next, cluster across both variables (together, and separately)
  dist_mats <- c(
    list(dist_obj$dist_mat), # combined distance matrix
    dist_obj$dist_mats # separate distance matrices for each variable
  )
  clust_obj_var <- lapply(dist_mats, \(y) {
    js_clust(
      dist_mat = y,
      dep_fit$dependence,
      k = k
    )
  })

  twgss <- c(
    list(scree_plot(dist_mats[[1]], k = 1:5)), # TWGSS for combined diss mat
    dist_obj$total_within_ss
  )

  # return distance matrices for each variable, as well as clustering + TWGSS
  return(list(
    "dist_mats" = dist_mats,
    # "clust_obj" = clust_obj,
    "clust_obj" = clust_obj_var,
    "twgss"     = twgss
  ))
})
saveRDS(clust_obj_lst, "data/laplace_trunc_clust_obj.RDS")

# quickly plot TWGSS elbow for each (k = 3 best for each!)
p_elbow <- bind_rows(lapply(seq_along(clust_obj_lst), \(i) {
  data.frame(
    "laplace_q" = laplace_q[i],
    # "twgss"     = clust_obj_lst[[i]]$twgss
    "twgss_combined" = clust_obj_lst[[i]]$twgss[[1]],
    "twgss_rain" = clust_obj_lst[[i]]$twgss[[2]],
    "twgss_wind_speed" = clust_obj_lst[[i]]$twgss[[3]]
  )
})) |>
  group_by(laplace_q) %>%
  mutate(row = row_number()) %>%
  ungroup() %>%
  pivot_longer(contains("twgss"), values_to = "twgss") |>
  # mutate(name = ifelse(grepl("rain", name), "Rain", "Wind Speed")) |>
  mutate(
    name = case_when(
      name == "twgss_rain" ~ "Rain",
      name == "twgss_wind_speed" ~ "Wind Speed",
      TRUE ~ "Both"
    )
  ) |>
  ggplot() +
  geom_line(aes(x = row, y = twgss, colour = name)) +
  geom_point(aes(x = row, y = twgss, colour = name)) +
  facet_wrap(~laplace_q, scales = "fixed") +
  theme +
  labs(x = "", y = "TWGSS", colour = "") +
  ggsci::scale_colour_nejm()
p_elbow


#### Compare divergences ####

# TODO Also add combined divergence matrices!!

# extract divergences
div_lst <- lapply(clust_obj_lst, `[[`, "dist_mats")

# transpose from laplace_q -> variable list to vice versa
div_lst <- purrr::transpose(div_lst)

# take first variable as test
# div_lst_spec <- div_lst[[1]]
#
# div_mean <- sapply(div_lst_spec, \(x) mean(x))
# div_max <- sapply(div_lst_spec, \(x) max(x))
# div_norm <- sapply(div_lst_spec, \(x) norm(as.matrix(x), type = "2"))
# div_kappa <- sapply(div_lst_spec, \(x) {
#   kappa(as.matrix(x), exact = TRUE, triangular = TRUE)
# })
#
# # combine into dataframe
# div_df <- data.frame(
#   "laplace_q" = laplace_q,
#   "mean"      = div_mean,
#   "max"       = div_max,
#   "norm"      = div_norm,
#   "kappa"     = div_kappa
# )
#
# # plot raw values
# # div_df |>
#   pivot_longer(cols = -laplace_q) |>
#   ggplot(aes(x = laplace_q, y = value)) +
#   geom_line() +
#   geom_point() +
#   facet_wrap(~name, scales = "free_y") +
#   theme
# # comments:
#
# # plot first derivative
# div_df |>
#   mutate(
#     across(c(mean, max, norm, kappa), \(x) {
#       c(NA, diff(x) / diff(laplace_q))
#     })
#   ) |>
#   pivot_longer(cols = -laplace_q) |>
#   ggplot(aes(x = laplace_q, y = value)) +
#   geom_line() +
#   geom_point() +
#   facet_wrap(~name, scales = "free_y") +
#   theme
# # comments: First big bump occurs at 0.85, so 0.8 may be a good candidate?
# # lines up well with simulation results!


# calculate mean, max, 2-norm and condition number for each
div_df <- bind_rows(lapply(seq_along(div_lst), \(i) {
  x <- div_lst[[i]]
  bind_rows(lapply(seq_along(x), \(j) {
    y <- x[[j]]
    data.frame(
      "mean" = mean(y),
      "max" = max(y),
      "norm" = norm(as.matrix(y), type = "2"),
      "kappa" = kappa(as.matrix(y), exact = TRUE, triangular = TRUE),
      "laplace_q" = laplace_q[j],
      "variable" = c("Both", vars)[i]
    )
  }))
}))

# plot raw values
div_df |>
  pivot_longer(cols = 1:4) |>
  # max, mean and 2-norm all tell us much the same thing
  filter(name %in% c("kappa", "norm")) |>
  ggplot(aes(x = laplace_q, y = value, colour = variable)) +
  geom_line() +
  geom_point() +
  # change y-axis to scientific notation
  scale_y_continuous(labels = scales::scientific) +
  facet_wrap(variable ~ name, scales = "free_y", nrow = 3) +
  # facet_grid(variable ~ name, scales = "free") +
  theme
# Comments
# Both:
# Rain:
# Wind speed:
# Conclusion:

# also plot first derivatives
div_df |>
  mutate(
    across(c(mean, max, norm, kappa), \(x) {
      # c(NA, diff(x) / diff(laplace_q))
      c(NA, diff(x))
    })
  ) |>
  pivot_longer(cols = 1:4) |>
  filter(name %in% c("kappa", "norm")) |>
  filter(laplace_q > 0.65) |>
  ggplot(aes(x = laplace_q, y = value, colour = variable)) +
  ggtitle("Derivatives") +
  geom_line() +
  geom_point() +
  scale_y_continuous(labels = scales::scientific) +
  facet_wrap(variable ~ name, scales = "free_y", nrow = 3) +
  theme
# Comments:


#### Plot clustering solutions on map ####

# # plot each clustering solution
# plt_lst <- lapply(seq_along(clust_obj_lst), \(i) {
#   plt_clust_map(
#     pts,
#     areas,
#     clust_obj_lst[[i]]$clust_obj,
#     elev_df = elev_df,
#     rm_elev_leg = TRUE
#   ) +
#     ggtitle(
#       label = paste0("Laplace truncation at ", laplace_q[i])
#     )
# })

plt_lst <- lapply(seq_along(clust_obj_lst), \(i) {
  wrap_plots(list(
    plt_clust_map(
      pts, areas, clust_obj_lst[[i]]$clust_obj[[1]],
      elev_df = elev_df,
      rm_elev_leg = TRUE
    ) +
      ggtitle("Both"),
    plt_clust_map(
      pts, areas, clust_obj_lst[[i]]$clust_obj[[2]],
      elev_df = elev_df,
      rm_elev_leg = TRUE
    ) +
      ggtitle("Rain"),
    plt_clust_map(
      pts, areas, clust_obj_lst[[i]]$clust_obj[[3]],
      elev_df = elev_df,
      rm_elev_leg = TRUE
    ) +
      ggtitle("Wind Speed")
  )) +
    patchwork::plot_annotation(
      title = paste0("Laplace truncation at ", laplace_q[i]),
      theme = evc::evc_theme()
    ) +
    NULL
})

ggsave("laplace_test0.png", plt_lst[[1]], width = 12, height = 12)
ggsave("laplace_test1.png", plt_lst[[2]], width = 12, height = 12)
ggsave("laplace_test2.png", plt_lst[[3]], width = 12, height = 12)
ggsave("laplace_test3.png", plt_lst[[4]], width = 12, height = 12)
ggsave("laplace_test4.png", plt_lst[[5]], width = 12, height = 12)
ggsave("laplace_test5.png", plt_lst[[6]], width = 12, height = 12)
ggsave("laplace_test6.png", plt_lst[[7]], width = 12, height = 12)

# save to inspect (takes ages!!)
# TODO Why isn't this saving properly??
pdf(
  file = "laplace_trunc_test.pdf",
  width = 12, height = 8
)
plt_lst
dev.off()
# conclusion: They all look pretty good except 0.7 and 0.95, others are largely
# similar, 0.75 and 0.9 have nice median sites if that's important!
# So pretty nice to have relatively stable clustering across truncation points!

# comments:

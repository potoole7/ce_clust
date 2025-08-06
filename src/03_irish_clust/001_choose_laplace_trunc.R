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
# Allow for different Laplace caps for each variable?? (done)
# TODO Maybe try different dependence quantiles also?? Load from other file
# TODO Might be better to plot individually for each variable, which might
# help compare Laplace quantiles more easily (could facet?)


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

n_cores <- parallel::detectCores() - 1

# dqu <- c(0.85, 0.88, 0.9)
# dqu <- 0.85 # happiest with this quantile from bootstrapping
# TODO Expand to more Laplace quantiles?
laplace_q <- seq(0.7, 0.95, by = 0.05)
# laplace_q <- seq(0.65, 0.95, by = 0.05)
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
# TODO Expand to more quantiles, not that expensive to compute!
# clust_obj_lst <- lapply(laplace_q, \(x) {
laplace_q_gran <- seq(
  min(laplace_q),
  max(laplace_q),
  by = abs(diff(laplace_q)[1] / 2) # twice as granular
)
laplace_q_gran <- sort(unique(c(laplace_q_gran, 0.975, 0.99, 0.999)))

# Calc divergence, TWGSS and clustering sol for given Lap trunc quantile
clust_obj_fun <- \(laplace_cap) {
  print(laplace_cap)
  # First, calculate divergence matrix and TWGSS for scree plot
  dist_obj <- js_clust(
    dep_fit$dependence,
    k = NULL,
    # scree_k = 1:5,
    scree_k = list(1:5, 1:5), # create scree plot for each variable
    n = n_mc,
    # mc_method = "uniform"
    mc_method = "laplace_trunc",
    laplace_cap = laplace_cap,
    return_dist = TRUE
  )

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

  names(dist_mats) <- names(clust_obj_var) <- names(twgss) <- c(
    "combined", vars
  )

  # return distance matrices for each variable, as well as clustering + TWGSS
  return(list(
    "dist_mats" = dist_mats,
    # "clust_obj" = clust_obj,
    "clust_obj" = clust_obj_var,
    "twgss"     = twgss
  ))
}

clust_obj_lst <- mclapply(laplace_q_gran, \(x) {
  clust_obj_fun(x)
}, mc.cores = n_cores)
names(clust_obj_lst) <- laplace_q_gran
saveRDS(clust_obj_lst, "data/laplace_trunc_clust_obj.RDS")
# clust_obj_lst <- readRDS("data/laplace_trunc_clust_obj.RDS")

# quickly plot TWGSS elbow for each (k = 3 best for each! Nice)
p_elbow <- bind_rows(lapply(seq_along(clust_obj_lst), \(i) {
  data.frame(
    # "laplace_q" = laplace_q[i],
    "laplace_q" = laplace_q_gran[i],
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
      # "laplace_q" = laplace_q[j],
      "laplace_q" = laplace_q_gran[j],
      "variable" = c("Both", vars)[i]
    )
  }))
}))

plot_div <- \(div_df) {
  # plot raw values
  p_raw <- div_df |>
    pivot_longer(cols = 1:4) |>
    # max, mean and 2-norm all tell us much the same thing
    filter(name %in% c("kappa", "norm")) |>
    mutate(name = factor(name, levels = c("norm", "kappa"))) |>
    ggplot(aes(x = laplace_q, y = value, colour = variable)) +
    geom_line() +
    geom_point() +
    # change y-axis to scientific notation
    scale_x_continuous(breaks = laplace_q_gran) +
    scale_y_continuous(labels = scales::scientific) +
    facet_wrap(variable ~ name, scales = "free", nrow = 3) +
    # facet_grid(variable ~ name, scales = "free") +
    theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # also plot first derivatives
  p_der <- div_df |>
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
    scale_x_continuous(breaks = laplace_q_gran) +
    scale_y_continuous(labels = scales::scientific) +
    facet_wrap(variable ~ name, scales = "free_y", nrow = 3) +
    theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  return(list("raw" = p_raw, "derivative" = p_der))
}

plot_div(div_df)
# Comments:

plot_div(filter(div_df, laplace_q > 0.75, laplace_q <= 0.9))

div_df |>
  # filter(laplace_q <= 0.95) |>
  group_by(variable) |>
  filter(kappa == min(kappa)) |>
  ungroup()
# kappa lowest for highest quantiles?? Surprising!


#### Plot clustering solutions on map ####

# function to plot each clustering solution, faceting by variable
facet_map_plot <- \(clust_obj_lst, laplace_q) {
  clust_obj_lst <- clust_obj_lst[names(clust_obj_lst) %in% laplace_q]
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

  return(plt_lst)
}
# plt_lst <- facet_map_plot(clust_obj_lst, laplace_q)
plt_lst <- facet_map_plot(clust_obj_lst, laplace_q_gran)

lapply(seq_along(plt_lst), \(i) {
  print(i)
  print(i)
  if (i < 10) {
    i_str <- paste0("0", i - 1)
  } else {
    i_str <- as.character(i - 1)
  }
  ggsave(
    paste0("laplace_test0", i_str, ".png"),
    plt_lst[[i]],
    width = 12,
    height = 10
  )
})

# save to inspect (takes ages!!)
# TODO Why isn't this saving properly??
# pdf(
#   file = "laplace_trunc_test.pdf",
#   width = 12, height = 8
# )
# plt_lst
# dev.off()

# Comments:
# Rain: Same for < 0.8, 2 points different for 0.85, then 1 diff for 0.9
# => Clustering is stable to choice of Laplace truncation
# Wind: Also same for <0.8, at 0.85 more in West cluster,  for 0.9 more points
# from East go to Centre, but also one to West cluster, so probably best not
# to go that far
# Both: Interestingly the least stable across truncation quantiles!
# is the same for 0.8 and 0.85, at 0.9 stubborn

#### Different Laplace truncation across variables ####

# also look at clustering solutions for wind_speed at 0.85 and different rain levels
laplace_q_wind <- 0.85
# laplace_q_wind <- lapply(laplace_q, \(x) {
laplace_q_wind <- lapply(laplace_q_gran, \(x) {
  c("rain" = x, "wind_speed" = laplace_q_wind)
})

clust_obj_lst_wind <- mclapply(laplace_q_wind, \(x) {
  clust_obj_fun(x)
}, mc.cores = n_cores)
names(clust_obj_lst_wind) <- sapply(laplace_q_wind, `[[`, "rain")

# plot each clustering solution
plt_lst_wind <- facet_map_plot(clust_obj_lst_wind, sapply(laplace_q_wind, `[[`, "rain"))
lapply(seq_along(plt_lst_wind), \(i) {
  print(i)
  if (i < 10) {
    i_str <- paste0("0", i - 1)
  } else {
    i_str <- as.character(i - 1)
  }
  ggsave(
    paste0("laplace_test_wind", i_str, ".png"),
    plt_lst_wind[[i]],
    width = 12,
    height = 10
  )
})
# Conclusions: Not much difference, at 0.9 again there is no red point in
# the South

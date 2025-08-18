#### Testing truncated Laplace on results of CE model for Ireland ####

# As opposed to `001_choose_laplace_trunc.R` (which uses quantiles of exp
# distribution to truncate Laplace), here we test two different empirical
# methods of choosing the Lapace truncation point:
# (i) Take 99th percentile of the empirical (Laplace-trasformed) data
# (ii) Take minima of sitewise maxima

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

# TODO Expand to more DQU values? (just 0.85 for now)
# dqu <- c(0.85, 0.88, 0.9)
emp_q <- c(0.95, 0.99, 0.995) # truncate Laplace at 99th percentile of empirical distribution

# grid <- crossing(
#   dqu = dqu,
#   emp_q = emp_q
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

# Calc divergence, TWGSS and clustering sol for given parameters
clust_obj_fun <- \(
  params
) {
  print(unlist(params))
  stopifnot(is.list(params))
  stopifnot(all(names(params) == c("mc_method", "laplace_cap")))
  # First, calculate divergence matrix and TWGSS for scree plot
  # dist_obj <- js_clust(
  #   dep_fit$dependence,
  #   k = NULL,
  #   # scree_k = 1:5,
  #   scree_k = list(1:5, 1:5), # create scree plot for each variable
  #   n = n_mc,
  #   # mc_method = "uniform"
  #   mc_method = "laplace_trunc",
  #   laplace_cap = laplace_cap,
  #   return_dist = TRUE
  # )

  # First, calculate divergence matrix and TWGSS for scree plot
  common_pars <- list( # parameters not in `params`
    dep_fit$dependence,
    k = NULL,
    scree_k = list(1:5, 1:5), # create scree plot for each variable
    n = n_mc,
    return_dist = TRUE,
    trans = dep_fit$transformed
  )
  dist_obj <- do.call(
    js_clust,
    c(common_pars, params) # add parameters from `params`
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
    # TODO Check if I should be doing for all elements of dist_mats??
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

# different setups to try
params_lst <- list(
  # Cutoff laplace at empirical quantile
  # TODO Loop over emp_q
  list(
    mc_method = "laplace_trunc2",
    laplace_cap = 0.95 # should be shit
  ),
  list(
    mc_method = "laplace_trunc2",
    laplace_cap = 0.99
  ),
  list(
    mc_method = "laplace_trunc2",
    laplace_cap = 0.995 # should be shit
  ),
  # Cutoff laplace at minima of sitewise maxima
  list(
    mc_method = "laplace_trunc3",
    laplace_cap = NULL
  )
)

# Run clustering for each set of parameters
set.seed(123)
# debugonce(js_clust)
# clust_obj_lst <- mclapply(params_lst, \(x) {
clust_obj_lst <- lapply(params_lst, \(x) {
  clust_obj_fun(x)
  # }, mc.cores = n_cores)
})

# give appropriate names
names(clust_obj_lst) <- sapply(
  params_lst, \(x) ifelse(
    is.null(x[["laplace_cap"]]),
    # NA,
    "min_of_sitewise_maxima",
    # x[["laplace_cap"])
    paste0("emp_q_", as.character(x[["laplace_cap"]]))
  )
)


#### Elbow plots ####

# quickly plot TWGSS elbow for each (k = 3 best for each! Nice)
p_elbow <- bind_rows(lapply(seq_along(clust_obj_lst), \(i) {
  data.frame(
    # "laplace_q" = laplace_q[i],
    "method" = params_lst[[i]]$mc_method,
    "emp_q" = ifelse(
      is.null(params_lst[[i]]$laplace_cap),
      NA,
      params_lst[[i]]$laplace_cap
    ),
    "twgss_combined" = clust_obj_lst[[i]]$twgss[[1]],
    "twgss_rain" = clust_obj_lst[[i]]$twgss[[2]],
    "twgss_wind_speed" = clust_obj_lst[[i]]$twgss[[3]]
  )
})) |>
  group_by(emp_q) %>%
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
  facet_wrap(~emp_q, scales = "fixed", nrow = 1) +
  theme +
  labs(x = "", y = "TWGSS", colour = "") +
  ggsci::scale_colour_nejm()
p_elbow
# Conclusion: All suggest 3, but interesting to note that we go further into
# the tails, the TWGSS (and therefore the similarity) appears to have increased

#### Map plot of clustering solutions ####

# function to plot each clustering solution, faceting by variable
# TODO Add titles
facet_map_plot <- \(clust_obj_lst) {
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
        # title = paste0("Laplace truncation at ", laplace_q[i]),
        title = names(clust_obj_lst)[i],
        theme = evc::evc_theme()
      ) +
      NULL
  })

  return(plt_lst)
}
plt_lst <- facet_map_plot(clust_obj_lst)

lapply(seq_along(plt_lst), \(i) {
  print(i)
  if (i < 10) {
    i_str <- paste0("0", i - 1)
  } else {
    i_str <- as.character(i - 1)
  }
  ggsave(
    # paste0("laplace_test0", i_str, ".png"),
    paste0("laplace_plots/emp_laplace_test0", i_str, ".png"),
    plt_lst[[i]],
    width = 12,
    height = 10
  )
})
# Conclusion:
# Emp 0.99 and sitewise maxima are both good
# Problems:
# - Both have West cluster point in Wicklow for wind (can be explained by mountains??)
# - Emp 0.99 has Donegal site has blue for both rain and wind but red together,
#   while in sitewise maxima site is at least red for rain
#   - Fine as we can just say combining variables highlights this difference
# Conclusion: Go with 99th empirical quantile


# best model is for 99th empirical quantile; take results for that
clust_obj <- clust_obj_lst$emp_q_0.99
saveRDS(clust_obj, "data/clust_obj_emp_laplace.RDS")

#### Compare Conditional Extremes output before and after clustering ####

# Refit conditional extremes model on clusters rather than individual locations
# compare results

#### Libs ####

library(sf)
library(evc)
# devtools::load_all("../evc")
library(dplyr, quietly = TRUE)
library(tidyr)
library(ggplot2)
library(lubridate)
# library(texmex)
devtools::load_all("texmex") # forked texmex which allows for fixed b vals
library(gridExtra)

source("src/functions.R")

sf::sf_use_s2(FALSE)


#### Load Data ####

# load original dataset
data <- readr::read_csv("data/met_eireann/final/met_eir_preprocess.csv.gz")

# pull county for locations
locs <- distinct(data, name, county)

# load shapefile
areas <- sf::read_sf("data/met_eireann/final/irl_shapefile.geojson")

# Fitted conditional extremes models for each site (i.e. orig unclustered fit)
dep_orig <- readRDS("data/texmex_data/texmex_mexdep_obj.RDS")
# # sort alphabetically by name
dep_orig <- dep_orig[sort(names(dep_orig))]

# also load marginals
marg_orig <- readRDS("data/texmex_data/texmex_marginal_obj.RDS")

# load clustering solution
clust <- readRDS("data/texmex_data/pam_clust2.RDS")

#### Setup ###

# change clustering solution into dataframe
clust_df <- clust$clustering %>% 
  as.data.frame() %>% 
  setNames("cluster") %>% 
  mutate(name = names(clust$clustering)) %>% 
  as_tibble() # remove rownames

# join to data, remove any locations with no cluster
data_clust_full <- data %>% 
  left_join(clust_df) %>% 
  filter(!is.na(cluster))

# group locations in the same cluster together
data_clust <- data_clust_full %>% 
  group_by(date, name = paste0("cluster_", cluster)) %>% 
  summarise(
    rain       = sum(rain, na.rm = TRUE),
    wind_speed = mean(wind_speed, na.rm = TRUE),
    .groups    = "drop"
  )


#### Fit CE w/ cluster specific parameters ####

mod_fit <- fit_ce(
  data_clust,
  cond_prob   = 0.7,  # TODO: Is this the best conditional prob? Bootstrap!
  # marg_prob   = 0.9,
  # cluster specific thresholds
  marg_prob   = list(
    f   = list("response ~ name", "~ name"), 
    tau = 0.9
  ),
  f           = list(excess ~ name, ~name), 
  fit_no_keef = FALSE, 
  output_all = TRUE
)
dep <- mod_fit$dependence
marg <- mod_fit$marginal


#### Fit CE w/ spatially varying cluster specific parameters ####

mod_fit_stat <- data_clust_full %>% 
  group_split(cluster, .keep = TRUE) %>% 
  lapply(\(x) {
    fit_ce(
      x,
      cond_prob   = 0.7,  # TODO: Is this the best conditional prob? Bootstrap!
      marg_prob   = list(
        f   = list("response ~ s(lon, lat, k = 15)", "~ s(lon, lat, k = 15)"),
        tau = 0.85
      ), 
      f           = list(excess ~ s(lon, lat, k = 15), ~s(lon, lat, k = 15)), 
      fit_no_keef = FALSE, 
      output_all = TRUE
    )
  })

dep_stat <- c(mod_fit_stat[[1]]$dependence, mod_fit_stat[[2]]$dependence)
marg_stat <- c(mod_fit_stat[[1]]$marginal, mod_fit_stat[[2]]$marginal)


#### Plot ####

# TODO: a, b map plot for orig, cluster, and spatially varying cluster

# First, extract ab values for each clustering solution
extract_ab <- \(dep, df = NULL) {
  ab_vals <- lapply(dep, \(x) {
    list(
      x[[1]]$dependence$coefficients[1:2], 
      x[[2]]$dependence$coefficients[1:2]
    )
  })
  ab_df <- bind_rows(lapply(seq_along(ab_vals), \(i) {
    data.frame(
      "name" = names(ab_vals)[i],
      "var"  = c("rain", "rain", "windspeed", "windspeed"),
      "parameter" = c("a", "b", "a", "b"),
      "value" = c(
        ab_vals[[i]][[1]][1], ab_vals[[i]][[1]][2],
        ab_vals[[i]][[2]][1], ab_vals[[i]][[2]][2]
      )
    )
  })) %>% 
    filter(!is.na(name))
  return(ab_df)
}

# join in plotting info and plot
ab_map_plot <- \(ab_df) {
  ab_sf <- dplyr::distinct(data_clust_full, name, lon, lat, cluster) %>% 
    left_join(ab_df) %>% 
    st_to_sf()
  
  scales <- seq(-1, 1, by = 0.2)
  # plot
  ggplot(areas) +
    geom_sf(colour = "black", fill = NA) +
    geom_sf(
      data = ab_sf,
      aes(fill = value, shape = factor(cluster), size = value), # , colour = value, , )
      # size = 6,
    ) +
    # facet_wrap(~var) + 
    facet_grid(var ~ parameter) + 
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    # default circles and triangles as shape
    scale_shape_manual(values = c(21, 24)) +
    scale_fill_gradient2(
      low = "blue3",
      high = "red3",
      na.value = "grey",
      breaks = scales, 
      labels = as.character(scales),
      guide = "legend"
    ) + 
    # scale_colour_gradient2(
    #   low = "blue3",
    #   high = "red3",
    #   na.value = "grey",
    #   breaks = scales, 
    #   labels = as.character(scales),
    #   guide = "legend"
    # ) + 
    scale_size_continuous(
      range = c(2, 7),
      breaks = scales,
      labels = as.character(scales),
      guide = "legend"
    ) +
    labs(fill = "", size = "") + 
    guides(
      fill = guide_legend(), 
      # size = guide_legend(), 
      shape = guide_legend(), 
      colour = "none"
    ) +
    evc_theme() + 
    theme(
      legend.position = "right",
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.key = element_blank()
    )
}

# original estimates
ab_orig <- extract_ab(dep_orig)
(p1 <- ab_map_plot(ab_orig))

p2 <- extract_ab(dep) %>% 
  mutate(cluster = ifelse(name == "cluster_1", 1, 2)) %>% 
  select(-name) %>% 
  full_join(
    distinct(data_clust_full, name, lon, lat, cluster), 
    relationship = "many-to-many"
  ) %>% 
  ab_map_plot()
p2

p3 <- extract_ab(dep_stat) %>% 
  full_join(
    distinct(data_clust_full, name, lon, lat, cluster), 
    # relationship = "many-to-many"
  ) %>% 
  ab_map_plot()
p3

# TODO: Plot difference between p1 and p2 and p3
# p2 - p1
extract_ab(dep) %>% 
  mutate(cluster = ifelse(name == "cluster_1", 1, 2)) %>%
  select(-name) %>% 
  full_join(distinct(data_clust_full, name, lon, lat, cluster)) %>% 
  left_join(rename(ab_orig, value_orig = value)) %>%
  mutate(
    value = value - value_orig
  ) %>%
  ab_map_plot()

# p3 - p1
extract_ab(dep_stat) %>% 
  full_join(distinct(data_clust_full, name, lon, lat, cluster)) %>% 
  left_join(rename(ab_orig, value_orig = value)) %>%
  mutate(
    value = value - value_orig
  ) %>%
  ab_map_plot()

# TODO: Something on uncertainty in marginal and/or dependence parameters??
# Can I create profile likelihood plots from dependence objects??
prof_lik <- mexDependence(
  marg_orig$`Aughnasheelan (Miskawn)`,
  which = "rain", 
  dqu = 0.7, 
  nOptim = 2, 
  PlotLikDo = TRUE, 
  # PlotLikTitle = "Profile Likelihood, rain | wind speed"
  PlotLikTitle = ""
)

prof_lik <- mexDependence(
  marg$cluster_1,
  which = "rain", 
  dqu = 0.7, 
  nOptim = 2, 
  PlotLikDo = TRUE, 
  # PlotLikTitle = "Profile Likelihood, rain | wind speed"
  PlotLikTitle = ""
)

prof_lik <- mexDependence(
  marg_stat$`Aughnasheelan (Miskawn)`,
  which = "rain", 
  dqu = 0.7, 
  nOptim = 2, 
  PlotLikDo = TRUE, 
  # PlotLikTitle = "Profile Likelihood, rain | wind speed"
  PlotLikTitle = ""
)

ggplot(marg_orig$`Aughnasheelan (Miskawn)`)
ggplot(marg$cluster_1)
ggplot(marg_stat$`Aughnasheelan (Miskawn)`)

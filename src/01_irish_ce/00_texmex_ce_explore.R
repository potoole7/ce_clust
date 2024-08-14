#### Investigating issues with texmex CE fitting to evgam output ####

#### Libs ####

library(sf)
library(dplyr, quietly = TRUE)
library(tidyr)
library(ggplot2)
library(lubridate)
library(texmex)
library(gridExtra)
library(patchwork)
library(evgam)

source("src/functions.R")

theme <- theme_bw() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 13, face = "bold"),
    strip.background = element_rect(fill = NA, colour = "black"),
    plot.tag = element_text(size = 16, face = "bold"),
    panel.background = element_rect(fill = NA, colour = "black")
  )

sf::sf_use_s2(FALSE)

#### Load data ####

# Load from `src/01_irish_ce/02_marginals_evgam_var_thresh.R`
a_b_sf <- sf::read_sf("a_b_sf.geojson")
marginal <- readRDS("marginal.RDS")
dependence <- readRDS("dependence.RDS")

#### Investigate large sigma in summary ####

# find name of location with minimum b value
min_b_name <- a_b_sf %>% 
  filter(value == min(value)) %>% 
  pull(name)

# find dependence and marginal for this
spec_marginal <- marginal[[min_b_name]]
spec_dependence <- dependence[[min_b_name]]

# confirm b value of ~-8 for rain (note: no longer that low!)
# summary(spec_dependence[[1]])

# Investigate very large (589) sigma value for marginal summary, which doesn't
# match coefficient?!
summary(spec_marginal)

# from texmex:::coef.migpd, need to have log of scale parameter in 
# coefficients, rather than actual value!! (done)

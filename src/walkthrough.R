#### Walkthrough of mod and clust procedure for Gaussian copula data ####

#### Libs ####

library(dplyr, quietly = TRUE)
devtools::load_all("../evc")
# devtools::load_all("../evc_mc/")
library(tidyr)
library(ggplot2)
library(patchwork)
library(parallel)

source("src/functions.R")


#### Metadata ####

# TODO Comment
seed_number <- 123
n <- 1e4
n_locs <- 4
n_vars <- 2
mix_p <- c("gauss_cop" = 1, "t_cop" = 0)
# mix_p <- c("gauss_cop" = 0.5, "t_cop" = 0.5)
cond_prob <- 0.95
n_clust <- 2
cluster_mem <- rep(1:n_clust, each = n_locs / n_clust) # cluster membership
cor_gauss <- c(0.1, 0.9)
# cor_gauss <- c(0.3, 0.3)
cor_t <- NULL
# cor_t <- c(0.1, 0.9)
df_t <- NULL
# df_t <- c(3, 3)


#### Simulate and Transform Data ####

set.seed(seed_number)
# TODO should I use t-copula/mixture instead??
data <- sim_cop_dat(
  n_locs      = n_locs,
  n_vars      = n_vars,
  n           = n,
  n_clust     = n_clust,
  cluster_mem = cluster_mem,
  cor_gauss   = cor_gauss,
  cor_t       = cor_t,
  df_t        = df_t,
  mix_p       = mix_p,
  qfun        = qnorm # standard normal margins
)
data_mix_gauss <- data$data_mix

# plot margins for all locations
plot_marg <- \(x, add_quant = FALSE) {
  p <- bind_rows(lapply(x, as.data.frame), .id = "loc") |>
    setNames(c("loc", "Y1", "Y2")) |>
    mutate(
      gauss_corr = ifelse(loc %in% c(1, 2), 0.9, 0.1),
      quant_y1   = quantile(Y1, cond_prob),
      quant_y2   = quantile(Y2, cond_prob)
    ) |>
    ggplot(aes(x = Y1, y = Y2, colour = factor(gauss_corr))) +
    geom_point() +
    facet_wrap(~loc) +
    evc::evc_theme() +
    labs(
      x = expression(Y[1]),
      y = expression(Y[2]),
      colour = expression(rho[Gauss])
    ) +
    guides(colour = guide_legend(override.aes = list(size = 7))) + # , shape = 15))) +
    # remove facet labels
    theme(strip.text = element_blank())

  if (add_quant == TRUE) {
    p <- p +
      geom_vline(aes(xintercept = quant_y1), linetype = "dashed") +
      geom_hline(aes(yintercept = quant_y2), linetype = "dashed")
  }
  return(p)
}

# plot for original Gaussian margins
(p1 <- plot_marg(data_mix_gauss))

# transform normal data to CDF and convert to Laplace margins
trans_fun <- \(x) laplace_trans(pnorm(x))
data_mix_trans <- lapply(data_mix_gauss, trans_fun)

(p2 <- plot_marg(data_mix_trans, add_quant = TRUE))


#### Fit CE model ####

# we know what asymptotic CE parameter values will be for Gaussian copula
start_vals <- lapply(cor_gauss, \(x) c(sign(x) * x^2, 1 / 2))

# fit CE to each location
ce_fit <- lapply(seq_along(data_mix_trans), \(i) {
  o <- ce_optim(
    Y = data_mix_trans[[i]],
    dqu = cond_prob,
    control = list(maxit = 1e6),
    constrain = FALSE,
    start = start_vals[[cluster_mem[i]]]
  )
})
# Returns list of len n_locs, w/ lists of len n_vars, of resids & params

# Tabulate nicely (done in LaTeX)
extract_pars <- \(ce_fit) {
  # # debugging
  # x <- ce_fit[[1]]
  # y <- ce_fit[[1]][[1]]
  bind_rows(lapply(ce_fit, \(x) { # loop through "locations"
    data.frame(
      sapply(x, \(y) { # loop through "variables"
        y$params[c("a", "b", "m", "s", "dth"), ]
        # }, numeric(n_vars)),
      }),
      "parameter" = c("a", "b", "m", "s", "dth"),
      row.names = NULL
    )
  }), .id = "location") |>
    rename(Y1 = var_1, Y2 = var_2) |>
    # have parameter as columns, with Y1 and Y2 as rows
    pivot_longer(
      cols = c("Y1", "Y2"),
      names_to = "cond_variable",
      values_to = "value"
    ) |>
    pivot_wider(
      names_from = parameter,
      values_from = value
    )
}

pars_df <- extract_pars(ce_fit)


#### Cluster ####

# sample from Laplace upper tail (== exponential sample above thresh)
yi <- pars_df |>
  group_split(cond_variable) |>
  lapply(\(x) {
    thresh_max <- max(x$dth)
    print(thresh_max)
    # samp <- thresh_max + rexp(200, rate = 1)
    samp <- thresh_max + rexp(5, rate = 1) # only do 5 for this example
    return(samp)
  })
vars <- unique(pars_df$cond_variable)
names(yi) <- vars

yi[[2]] |> round(3)

# as an example, estimate MVNs when conditioning on Y2 for sites 1 and 3
locs <- c(1, 3)
mvn_ex <- bind_rows(lapply(seq_along(locs), \(i) {
  pars_df |>
    # filter(location == 1, cond_variable == "Y2") |>
    filter(location == locs[i], cond_variable == "Y2") |>
    mutate(
      mu = paste(round(a * yi[[i]] + (yi[[i]]^b) * m, 2), collapse = ", "),
      sigma = paste(round((yi[[i]]^b) * s, 2), collapse = ", ")
    )
}))

paste0(
  "N(",
  mvn_ex$mu,
  ", ",
  "diag(",
  mvn_ex$sigma,
  ")",
  ")"
)

# Perform clustering
profvis::profvis({
  clust_sol <- js_clust(
    ce_fit,
    k = 2, cluster_mem = cluster_mem, return_dist = TRUE, n = 500
  )
})

# show dissimilarity matrices
lapply(
  lapply(c(clust_sol$dists, list(clust_sol$dist)), as.matrix),
  round, 3
)

# Clustering solution
list("sol" = clust_sol$pam$clustering, "adj_rand" = clust_sol$adj_rand)

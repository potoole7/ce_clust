#### Derive own MLE for conditional extremes model ####

# Derive custom CE model, compare to texmex estimates

#### Libs ####

devtools::load_all("../evc")
# library(evc)
# devtools::load_all("texmex")
library(ismev)
library(copula)

source("src/functions.R")


### Metadata ####

n <- 1e3
cor <- 0.6
dim <- 4
df <- 3
xi <- 0.1
sigma <- 1
v <- 10
dqu <- 0.7


#### Simulate data ####

# sim data from 1 "location" as mixture of Gaussian + t-cop w/ GPD margins
set.seed(123)
dat <- sim_cop_dat(
  n_locs     = 1, 
  n_vars     = dim, 
  n          = n, 
  n_clust    = 1,
  cor_gauss  = cor, 
  cor_t      = cor,
  df_t       = df, 
  params_gpd = c(xi, sigma)
)$data_mix[[1]]


#### Fit marginal model ####

#f it GPD using ismev, rather than texmex
set.seed(123)
gpd <- lapply(seq_len(ncol(dat)), \(i) {
  ismev::gpd.fit(dat[, i], threshold = quantile(dat[, i], 0.9))
})

# compare to texmex
set.seed(123)
mex_marg <- texmex::migpd(dat, mth = c(gpd[[1]]$threshold, gpd[[2]]$threshold))

# quite different! May just be due to uncertainty in xi, scale though
gpd_orig <- gpd
lapply(seq_along(gpd_orig), \(i) {
  rbind(
    c(log(gpd_orig[[i]]$mle[1]), gpd_orig[[i]]$mle[1]), 
    mex_marg$models[[i]]$par
  )
})

# Use estimates from texmex to ensure dependence model is correctly estimated
gpd <- lapply(gpd, \(x) {
  list("sigma" = x$mle[1], "xi" = x$mle[2], "thresh" = x$threshold)
})
gpd <- lapply(seq_along(gpd), \(i) {
  return(list(
    "sigma"  = exp(mex_marg$models[[i]]$par[1]),
    "xi"     = mex_marg$models[[i]]$par[2],
    "thresh" = mex_marg$models[[i]]$threshold
  ))
})


#### Convert to CDF ####

F_hat <- semi_par(dat, gpd)

# check to see if CDF created correctly
sapply(seq_len(dim), \(i) {
  plot(dat[, i], F_hat[, i])
})


#### Marginal transformation to Laplace ####

# transform to Laplace margins
Y <- laplace_trans(F_hat)

# check shape
sapply(seq_len(dim), \(i) {
  plot(dat[, i], Y[, i])
})

# check this matches texmex 
Y_mex <- texmex:::mexTransform(
  mex_marg, 
  margins = list(
    p2q = function (p) ifelse(p < 0.5, log(2 * p), -log(2 * (1 - p))), 
    q2p = function (q) ifelse(q < 0, exp(q) / 2, 1 - 0.5 * exp(-q))
  )
)$transformed

lapply(seq_len(dim), \(i) {
  plot(Y[, i], Y_mex[, i], xlim = c(0, 5))
  norm(Y[, i, drop = FALSE] - Y_mex[, i, drop = FALSE])
})
norm(Y - Y_mex)


#### Fit CE ####

o <- ce_optim(
  Y = Y, 
  dqu = dqu,
  start = c("a" = 0.01, "b" = 0.01), 
  control = list(maxit = 1e6), 
  constrain = TRUE,
  v = 10, 
  aLow = -1
)

# compare to results from mexDependence
mexfit <- texmex::mexDependence(mex_marg, which = 1, dqu = dqu) 
# compare marginals
mexfit$dependence$coefficients[c(1, 2, 5, 6), ]
o[[1]]

norm(mexfit$dependence$coefficients[c(1, 2, 5, 6), ] - o[[1]])

# TODO: Calculate quantiles of CE model (needed?)

#### Include in fit_ce, compare to texmex & fun above ####

evc_fit <- fit_ce(
  dplyr::mutate(as.data.frame(dat), name = "loc_1"), # TODO: Add error if name not present
  vars = paste0("V", seq_len(dim)),
  # marg_prob = list(
  #   tau = 0.95, 
  #   jitter = TRUE
  # ),
  marg_prob = 0.9,
  f = NULL,
  cond_prob = dqu, 
  output_all = TRUE
)

# check margins are the same
vapply(seq_along(evc_fit$marginal$loc_1), \(i) {
  norm(as.matrix(unlist(evc_fit$marginal$loc_1[[i]]) - unlist(gpd[[i]])))
}, numeric(1))

# check dependence parameter estimates are the same
norm(evc_fit$dependence$loc_1$V1 - o[[1]])

#### multiple locations w/ `evgam` ####

set.seed(123)
dat3 <- sim_cop_dat(
  n_locs     = 3, 
  n_vars     = dim, 
  n          = n, 
  n_clust    = 1,
  cor_gauss  = cor, 
  cor_t      = cor,
  df_t       = df, 
  params_gpd = c(xi, sigma)
)$data_mix

# devtools::load_all("../evc")
# debugonce(fit_ce)
evc_fit <- fit_ce(
  data = dat3,
  vars = paste0("V", seq_len(dim)),
  marg_prob = list(
    f = list("response ~ name", "~ name"),
    tau = 0.95, 
    jitter = TRUE
  ),
  # marg_prob = 0.9,
  # f = NULL,
  cond_prob = dqu
)
 
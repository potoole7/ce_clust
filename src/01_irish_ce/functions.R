#### Functions ####

#### Closed form normal-normal KL/Jensen-Shannon Divegence ####

# pull a, b, mu, sigma for a given location
# dep - List of mex objects (i.e. CE models for each var) for single location
pull_params <- \(dep) {
  # fail if not list of mex objects for single location
  stopifnot(is.list(dep))
  # stopifnot(all(vapply(dep, class, character(1)) == "mex"))
  # loop through conditioning variables for single location
  return(lapply(dep, \(x) {
    # pull parameter vals
    ret <- as.vector(x$dependence$coefficients)
    names(ret) <- rownames(x$dependence$coefficients)
    # remove c and d if 0 (i.e. Laplace margins, rather than Gumbel)
    if (ret["c"] == 0 && ret["d"] == 0) {
      ret <- ret[!names(ret) %in% c("c", "d")]
    } 
    return(ret)
  }))
}

# pull thresholds
# dep - List of mex objects (i.e. CE models for each var) for single location
pull_thresh_trans <- \(dep) {
  # fail if not list of mex objects for single location
  stopifnot(is.list(dep))
  stopifnot(all(vapply(dep, class, character(1)) == "mex"))
  # return quantile of transformed data (already calculated in texmex)
  return(lapply(dep, \(x) x$dependence$dth))
}

# function to calculate KL divergence for single data point
# derivation at https://statproofbook.github.io/P/norm-kl.html
kl_gauss <- \(mu1, mu2, var1, var2) {
  return(1 / 2 * (
    (mu2 - mu1)^2 / var2 + (var1 / var2) - log(var1 / var2) - 1
  ))
}

# function to calculate Jensen-Shannon divergence metric from KL divergence
# This metric is symmetric, as required for clustering 
# https://tinyurl.com/526rwy9f
js_gauss <- \(mu1, mu2, var1, var2) {
  # calculate mean, variance for mixture distribution M = (P + Q)/2
  # Sum of normals as in https://online.stat.psu.edu/stat414/lesson/26/26.1
  # (N(u1, v1) + N(u2, v2)) / 2 = N((u1 + u2) / 2, (v1 + v2) / 2^2)
  mu_m <- (mu1 + mu2) / 2 
  var_m <- (var1 + var2) / 4
  
  # calculate JS(P||Q) = ((KL(P||M)) + KL(Q||M))/2
  return(
    (kl_gauss(mu1, mu_m, var1, var_m) + kl_gauss(mu2, mu_m, var1, var_m)) / 2
  )
}

# Function to calculate Jensen-Shannon divergence for each data point
js_div <- \(params_x, params_y, thresh_max, data_max = 2 * thresh_max, n_dat) {
  
  # test that input vectors have correct conditional extremes parameters
  stopifnot(is.vector(params_x) && is.vector(params_y))
  stopifnot(
    all(c(names(params_x), names(params_y)) == rep(c("a", "b", "m", "s"), 2))
  )
  
  # create data sequence from specified arguments
  data <- seq(thresh_max, data_max, length = n_dat)
  
  # funs to calculate mu and sd for normal dist as in 5.2 of Heff & Tawn '04
  mu_fun <- \(x, data) {
    return(x[["a"]] * data + x[["m"]] * (data ^ x[["b"]]))
  }
  var_fun <- \(x, data) {
    # take square now to avoid in kl_single formula
    return((x[["s"]] * (data ^ x[["b"]]))^2)
  }
  
  # calculate mu and sigma for each data point
  mus <- lapply(list(params_x, params_y), mu_fun, data = data)
  vars <- lapply(list(params_x, params_y), var_fun, data = data)
  
  # Calculate Jensen-Shannon divergence for each data point
  # TODO: How best to summarise across all data points? Sum? Average?
  return(sum(mapply(
    js_gauss,
    mu1 = mus[[1]], mu2 = mus[[2]], var1 = vars[[1]], var2 = vars[[2]]
  )))
}


#### Simulation functions ####

# Function to generate copula data
sim_cop_dat <- \(
  n_locs = 12,     # number of locations
  n = 1e4,         # number of simulations
  cor_gauss,       # bulk correlation for two clusters from Gaussian copula
  # params_norm,     # normal marginal parameters (same for both)
  cor_t,           # extreme correlation for two clusters from t-copula
  df_t,            # degrees of freedom for t-copula
  params_gpd,      # GPD margin parameters
  mix_p            # mixture weights
) {
  # many arguments must be of length 2 for 2 clusters
  stopifnot(all(vapply(
    # list(cor_gauss, params_norm, cor_t, df_t, params_gpd, mix_p), 
    list(cor_gauss, cor_t, df_t, params_gpd, mix_p), 
    \(x) length(x) == 2, logical(1)
  )))
  stopifnot(sum(mix_p) == 1)
  
  # Simulate from Gaussian Copula with GPD margins
  gauss_cop <- lapply(seq_len(n_locs), \(i){
    # pull correlation specified for each cluster
    cor <- cor_gauss[[1]]
    if (i > floor(n_locs / 2)) {
      cor <- cor_gauss[[2]]
    }
    cop_norm <- normalCopula(cor, dim = 2, dispstr = "un")
    # cop_norm <- normalCopula(cor_mat, dim = 2, dispstr = "un")
    u <- rCopula(n, cop_norm)
    # return(qnorm(u, mean = mu, sd = sigma))
    evd::qgpd(
      p     = u,
      # loc   = max(gauss_cop[[i]]), 
      loc   = 0,
      scale = scale_gpd, 
      shape = shape_gpd
    )
  })

  # simulate from t-Copula with GPD margins
  t_cop <- lapply(seq_len(n_locs), \(i) {
    cor <- cor_t[[1]]
    df <- df_t[[1]]
    if (i > floor(n_locs / 2)) {
      cor <- cor_t[[2]]
      df <- df_t[[2]]
    }
    cop_t <- copula::tCopula(cor, dim = 2, df = df_t[[2]], dispstr = "un")
    u <- rCopula(n, cop_t)
    return(evd::qgpd(
      p     = u,
      # loc   = max(gauss_cop[[i]]), 
      loc   = 0,
      scale = scale_gpd, 
      shape = shape_gpd
    ))
  })
  
  # mixture
  # data_mix <- lapply(seq_len(n_locs), \(i){
  #   mix_p[[1]] * gauss_cop[[i]] + mix_p[[2]] * t_cop[[i]]
  # })
  data_mix <- lapply(seq_len(n_locs), \(i) {
    x <- nrow(gauss_cop[[i]])
    y <- nrow(t_cop[[i]])
    # sample mix_p * nrow(gauss_cop) rows from gauss_cop (same for t_cop)
    rbind(
      gauss_cop[[i]][sample(seq_len(x), size = mix_p[[1]] * x), ],
      t_cop[[i]][sample(seq_len(y), size = mix_p[[2]] * y), ]
    )
  })
 return(list(
   "gauss_cop" = gauss_cop, 
   "t_cop"     = t_cop,
   "data_mix"  = data_mix
 ))
}

#### Vignotto 2021 ####

# Function to calculate KL divergence for 2 cluster simulation dataset
kl_sim_eval <- \(
  data_mix, # Mixture data from copulas
  kl_prob,  # Extremal quantile
  k = 2,    # number of clusters
  cluster_mem = NULL, # known cluster membership, for ARI
  ...
) {
  
  # prob must be valid probability
  stopifnot(0 <= kl_prob && kl_prob <= 1)
  
  # KL divergence between areas using Vignotto 2021 method
  kl_mat <- proxy::dist(
    data_mix, method = emp_kl_div, print = FALSE, prob = kl_prob, ...
  )
  
  # clustering solution
  pam_kl_clust <- pam(kl_mat, k = k)
  ret <- list("pam" = pam_kl_clust)
  # evaluate quality
  if (!is.null(cluster_mem)) {
    adj_rand <- mclust::adjustedRandIndex(
      pam_kl_clust$clustering, 
      cluster_mem
    )  
    # print(paste0("adjusted Rand index" = adj_rand))
    ret <- c(ret, list("adj_rand" = adj_rand))
  }
  return(ret)
}


# function to convert to Pareto scale
# divide by N + 1, not N (so do i / n + 1)
# x: Vector representing e.g. rain
pareto_trans <- \(x) {
  stopifnot(is.vector(x))
  
  # order and sort data
  x_ord <- order(x)
  x_sort <- x[x_ord]
  
  # calculate ECDF
  n <- length(x)
  ecdf_vals <- (seq_len(n)) / (n + 1)
  # convert back to original order
  ecdf_vals_x_ord <- numeric(n)
  ecdf_vals_x_ord[x_ord] <- ecdf_vals
  
  # pareto transform 
  return(1 / (1 - ecdf_vals_x_ord))
}

# function to compute bivariate risk function
# x: list of vectors representing different 
# fun: fun used to calculate "risk", defaults to max as easier to partition
risk_fun <- \(x, fun = max) {
  # test list
  stopifnot(is.list(x)) 
  # test equal length
  stopifnot(length(unique(vapply(x, length, numeric(1)))) == 1) 
  
  # TODO: Improve/speedup (could have x as a dataframe!)
  risk_vec <- vector(length = length(x[[1]]))
  for (i in seq_along(x[[1]])) {  
    risk_vec[[i]] <- fun(
      # pull ith entry in each vector in the list x
      vapply(x, `[[`, i, FUN.VALUE = numeric(1))
    )
  }
  return(risk_vec)
}

# Partition into 3 subsets, calculate empirical proportion in each
# TODO: Only works for max risk fun, unsure how to do for sum...
# TODO: Fix, obviously wrong since it doesn't identify times where both 
# are high!
partition_max <- \(x, prob, plot = FALSE) {
  
  # x must be a list of length 2 (rain and wind speed)
  stopifnot(is.list(x) && length(x) == 2) 
  
  # convert to dataframe, easier to subset (must follow this order!)
  df <- data.frame(
    rain       = x[[1]], 
    wind_speed = x[[2]], # TODO: Make x a dataframe already
    # TODO: Fix this, wrong somehow!!!
    R          = risk_fun(x, max)
  ) 
  
  # calculate quantile of risk function
  qu <- quantile(df$R, prob = prob)[[1]]
  
  # partition into 3 subsets
  df <- df %>% 
    mutate(extreme = case_when(
      rain > qu & wind_speed > qu ~ "both",
      rain > qu                    ~ "rain",
      wind_speed > qu              ~ "wind_speed",
      TRUE                          ~ NA
    ))
  
  if (plot) {
    df %>% 
      ggplot() +
      geom_point(aes(x = rain, y = wind_speed, colour = extreme)) + 
      geom_vline(xintercept = qu, linetype = "dashed") + 
      geom_hline(yintercept = qu, linetype = "dashed") + 
      ggsci::scale_colour_nejm() + 
      labs(y = "wind speed", title = paste0(
        "Hazard subsets, q_u = ", round(qu, 3), " (prob = ", prob, ")"
      )) + 
      theme
  }
  
  # return list with points falling into each category
  return(list(
      "rain"       = sum(df$extreme == "rain", na.rm = TRUE),
      "wind_speed" = sum(df$extreme == "wind_speed", na.rm = TRUE),
      "both"       = sum(df$extreme == "both", na.rm = TRUE)
    )
  )
}

# Function to calculate proportion 
calc_prop <- \(x) {
  stopifnot(is.list(x) & names(x) == c("rain", "wind_speed", "both"))
  # denom <- length(unlist(x))
  denom <- sum(unlist(x))
  ret <- lapply(x, \(y) y / denom)
  names(ret) <- names(x)
  return(ret)
}

# function to calculate KL divergence between any two locations
emp_kl_div <- \(x, y, convert_pareto = TRUE, prob = 0.9, print = TRUE, plot = FALSE) {
  # stopifnot(length(x) == length(y))
  
  # split x and y in half (rain vs wind speed)
  # remove NAs from x and y from padding matrix
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  n <- length(x)
  m <- length(y) # need for m???
  # TODO: Make more general for multivariate (not bivariate) case
  df_lst <- list(
  # x_df <- data.frame(
    "x" = data.frame(
      "rain"       = x[1:(n / 2)],
      "wind_speed" = x[((n / 2) + 1):n]
    ),
    "y" = data.frame(
    "rain"       = y[1:(m / 2)],
    "wind_speed" = y[((m / 2) + 1):m]
    )
  )
  
  if (print) {
    n_both <- lapply(df_lst, \(z) {
      nrow(filter(
        z, 
        rain > quantile(rain, prob), 
        wind_speed > quantile(wind_speed, prob)
      ))
    })
    print(paste0("Number both extreme for x before rescaling: ", n_both[[1]]))
    print(paste0("Number both extreme for y before rescaling: ", n_both[[2]]))
  }
   
  # convert to Pareto scale
  df_lst_par <- df_lst
  if (convert_pareto) {
    df_lst_par <- lapply(df_lst, \(z) {
      data.frame(apply(z, 2, pareto_trans))
    })
    # plot on pareto margins
    if (plot) {
      par(mfrow = c(1, 2)) # TODO: Make more general for MV case
      lapply(df_lst_par, plot)
      par(mfrow = c(1, 1))
    }
  }
  
  # partition into 3 subsets
  df_part <- lapply(df_lst_par, \(z) {
    partition_max(list(z[, 1], z[, 2]), prob = prob)
  })
  
  # fail if there are any 0s
  stopifnot(
    "No set can contain zeros, reduce `prob`" = all(unlist(df_part) != 0)
  )
  
  # calculate proportions of partitions
  df_part_prop <- lapply(df_part, \(z) {
    denom <- sum(unlist(z)) # denominator is # extreme obs
    lapply(z, `/`, denom) # find prop of extreme obs for each disjoint set
  })
  
  # calculate proportions of partitions
  x_part <- df_part_prop[[1]]
  y_part <- df_part_prop[[2]]
  if (print) {
    print(paste0("both for x: ", round(x_part$both, 3)))
    print(paste0("both for y: ", round(y_part$both, 3)))
  }
  sum_vals <- vector(length = length(x_part))
  for (i in seq_along(sum_vals)) { 
    sum_vals[i] <- (x_part[[i]] - y_part[[i]]) * 
                      log(x_part[[i]] / y_part[[i]])
  }
  # sum_vals[is.nan(sum_vals) | is.infinite(sum_vals)] <- 0
  if (any(is.nan(sum_vals) | is.infinite(sum_vals), na.rm = TRUE)) return(NA)
  return((1  / 2) * sum(sum_vals))
}



#### Conditional Extremes functions ####

# Function to calculate and cluster on the Jensen-Shannon Divergence for 
# the conditional extremes model
# TODO: functionalise n_dat
js_clust <- \(
  dependence,
  nclust = 3,
  cluster_mem = NULL
) {
  
  # TODO: Add stopifnot clause for class of dependence
   
  # pull parameter values for each location
  params <- lapply(dependence, pull_params)
  
  # pull Laplace threshold values for each location
  # TODO: Problem, last threshold different to others, investigate
  thresh <- lapply(dependence, pull_thresh_trans)
  
  # take maximum Laplace thresholds; want to geenra
  # thresh_max <- apply(bind_rows(thresh), 2, max)
  thresh_max <- lapply(bind_rows(thresh), max)
  
  # list of locs -> list of len 2 of variables, each containing all locs
  params <- purrr::transpose(params)
  
  dist_mats <- lapply(seq_along(params), \(i) {
    proxy::dist(
      params[[i]], 
      method = js_div, 
      thresh_max = thresh_max[[i]], 
      data_max = 2 * thresh_max[[i]], 
      n_dat = 10
    )
  })
  
  # sum distance matrices over different variables together
  dist_mat <- do.call(`+`, dist_mats)
  # scree plots look to suggest 3 clusters for both
  # lapply(dist_mats, scree_plot) 
  # lapply(dist_mats, scree_plot, fun = kmeans) 
  
  # cluster for rain and wind speed using both k-means and PAM
  # js_clust <- lapply(dist_mats, \(x) {
  #   lapply(c(pam, kmeans), \(fun) {
  #     fun(x, 3)
  #   })
  # })
  # cluster for rain and wind speed using PAM
  pam_js_clust <- pam(dist_mat, k = nclust)
  ret <- list("pam" = pam_js_clust)
  # evaluate quality
  if (!is.null(cluster_mem)) {
    adj_rand <- mclust::adjustedRandIndex(
      pam_js_clust$clustering, 
      cluster_mem
    )
    # print(paste0("adjusted Rand index" = adj_rand))
    ret <- c(ret, list("adj_rand" = adj_rand))
  }
  
  return(ret)
}

# Function to fit CE model
fit_ce <- \(
  data_mix, 
  vars = c("rain", "wind_speed"), 
  marg_prob = 0.9,
  cond_prob = 0.9,
  f = list(excess ~ name, ~ 1), # keep shape constant for now
  split_data = TRUE
) {
  
  # convert to dataframe
  # TODO: Tidy up, quite verbose
  if (split_data) {
    data_df <- bind_rows(lapply(seq_along(data_mix), \(i) {
      n <- length(data_mix[[i]])
      data.frame(
        "rain"       = data_mix[[i]][1:(n / 2)],
        "wind_speed" = data_mix[[i]][((n / 2) + 1):n]
      ) %>% 
        mutate(name = paste0("location_", i))
    }))
  } else {
    data_df <- bind_rows(lapply(seq_along(data_mix), \(i) {
      as.data.frame(data_mix[[i]]) %>% 
        mutate(name = paste0("location_", i))
    }))
    names(data_df)[1:2] <- c("rain", "wind_speed")
  }
  
  # First, calculate threshold (90th quantile across all locations)
  thresh <- apply(data_df[, 1:2], 2, quantile, marg_prob)
  
  # for each variable, calculate excess over threshold
  data_thresh <- lapply(vars, \(x) {
    data_df %>% 
      dplyr::select(matches(x), name) %>% 
      mutate(
        thresh = thresh[x],
        excess = !!sym(x) - thresh
      ) %>% 
      filter(excess > 0)
  })
  
  # Now fit evgam model for each marginal
  # TODO: Is just fitting different models by loc appropriate? (Yes I think!)
  evgam_fit <- lapply(data_thresh, \(x) {
    fit_evgam(
      data = x, 
      pred_data = x,
      # model scale and shape for each location
      # f = list(excess ~ name, ~ name) 
      f = f 
    )
  })
  
  # Join scale and shape estimates into data
  # TODO: Functionalise to work with > 2 variables
  data_gpd <- distinct(data_df, name) %>% 
    bind_cols(
      rename(
        distinct(evgam_fit[[1]]$predictions), 
        scale_rain = scale, 
        shape_rain = shape),
      rename(
        distinct(evgam_fit[[2]]$predictions), 
        scale_ws = scale, 
        shape_ws = shape),
    )
  
  # Now convert marginals to migpd (i.e. texmex format)
  marginal <- gen_marg_migpd(data_gpd, data_df)
  names(marginal) <- paste0(data_gpd$name, " - ", data_gpd$county)
  
  # Calculate dependence from marginals
  return(fit_texmex_dep(
    marginal, 
    mex_dep_args = list(dqu = cond_prob), 
    fit_no_keef = TRUE 
  ))
}

# function to fit varying threshold using quantile regression
# https://empslocal.ex.ac.uk/people/staff/by223/software/gpd.R
# TODO: Vary tau and see how that effects results (QQ plots, etc)
# quantile_thresh <- function(data, response, tau = .95) {
quantile_thresh <- function(data, response, tau = .95, jitter = TRUE) {
  fm_ald <- list(
    # response ~ s(lon, lat, k = 50), # location
    formula(paste(response, " ~ s(lon, lat, k = 50)")), # location
    ~ s(lon, lat, k = 40)                               # logscale
  )
  
  # jitter, if specified, to remove 0s when calculating quantiles
  if (jitter == TRUE) {
    data <- data %>%
      mutate(across(all_of(response), ~ . + rnorm(n(), 0, 1e-6)))
  }
  
  # fit the quantile regression model at tau'th percentile
  ald_fit <- evgam(fm_ald, data, family = "ald", ald.args = list(tau = tau))
  
  # add threshold to data and filter
  data %>% 
    mutate(
      thresh = predict(ald_fit)$location, 
      # excess = wind_speed - thresh
      excess = !!sym(response) - thresh
    ) %>% 
    filter(excess > 0) %>% 
    return()
}

# Fit evgam model with spline for spatial location, create predictions
fit_evgam <- \(
  data, 
  pred_data,
  # formula used in evgam, fitting to both scale and shape parameters
  f = list(
    excess ~ s(lon, lat, k = 40), # increase smoothing on scale parameter
    ~ s(lon, lat, k = 40) # shape parameter
  ) 
) {
  # model formula
  # f <- list(
  #   excess ~ s(lon, lat, k = 40), # increase smoothing on scale parameter
  #   ~ s(lon, lat, k = 40) # shape parameter
  # ) 
  # fit evgam model
  m <- evgam::evgam(f, data = data, family = "gpd")

  # create predictions
  predictions <- evgam:::predict.evgam(m, pred_data, type = "response")
  
  # return model fit and predictions
  return(list(
    "m"           = m,
    "predictions" = predictions
  ))
}

# create marginal `migpd` objects from `evgam` objects for each site
# - data_gpd: scale and shape parameters for each location
# - data: Data for each location
gen_marg_migpd <- \(data_gpd, data, mqu = 0.95) {
  # Create "dummy" migpd object to fill in with evgam values
  dat_mat <- data %>% 
    filter(name == data$name[[1]]) %>% 
    dplyr::select(rain, wind_speed) %>% 
    as.matrix()
  names(dat_mat) <- c("rain", "wind_speed")
  
  temp <- texmex::migpd(dat_mat, mqu = mqu, penalty = "none")
  # m <- evm(y = rain, data = data, qu = 0.95, penalty = "none", famuly = "gpd")
  # m1 <- update(m, phi = ~lon + lat)
  
  marginal <- lapply(seq_len(nrow(data_gpd)), \(i) {
    # initialise
    # browser()
    spec_marg <- temp
    # replace data 
    spec_marg$data <- data %>% 
      filter(name == data_gpd$name[i]) %>% 
      dplyr::select(rain, wind_speed) %>% 
      as.matrix()
    names(spec_marg$data) <- c("rain", "wind_speed")
    # replace thresholds
    # spec_marg$models$rain$threshold <- thresh_rain
    # spec_marg$models$wind_speed$threshold <- thresh_wind
    spec_marg$models$rain$threshold <- data_gpd$thresh_rain[[i]]
    spec_marg$models$wind_speed$threshold <- data_gpd$thresh_wind[[i]]
    # replace coefficients
    spec_marg$models$rain$coefficients[1:2] <- c(
      # data_gpd$scale_rain[i], 
      log(data_gpd$scale_rain[i]),
      data_gpd$shape_rain[i]
    )
    spec_marg$models$wind_speed$coefficients[1:2] <- c(
      # data_gpd$scale_ws[i], 
      log(data_gpd$scale_ws[i]), 
      data_gpd$shape_ws[i]
    )
    
    return(spec_marg)
  })
  return(marginal)
}

# Fit CE dependence model for each site
fit_texmex_dep <- \(
  marginal, 
  vars = c("rain", "wind_speed"), 
  mex_dep_args = list(
    start = c(0.01, 0.01), 
    dqu = 0.7,
    fixed_b = FALSE,
    PlotLikDo = FALSE
  ),
  fit_no_keef = FALSE
) {
  dependence <- lapply(seq_along(marginal), \(i) {
    # fit for rain and wind speed
    ret <- lapply(vars, \(col) {
      # TODO: Can you supply threshold yourself? or match to thresholding value?
      # TODO: Plot profile likelihood of dependence model parameters?
      mod <- do.call(
        texmex::mexDependence, 
        args = c(list(x = marginal[[i]], which = col), mex_dep_args)
      )
      
      # if model didn't optimise with Keef 2013 constrains, return NA
      ll <- mod$dependence$loglik
      if (is.na(ll) || abs(mod$dependence$loglik) > 1e9) {
        message("Model not fitting properly under Keef constraints")
        if (fit_no_keef) {
          mod <- do.call(
            texmex::mexDependence, 
            args = c(
              list(x = marginal[[i]], which = col, constrain = FALSE), 
              mex_dep_args
            )
          )
        } else {
          return(NA)
        }
      }
      return(mod)
    })
    names(ret) <- vars
    return(ret)
  })
  names(dependence) <- names(marginal)
  return(dependence)
}


#### Utility functions ####

# calculate adjacency matrix from Voronoi cells
calc_adj_mat <- \(pts, cut_vor = TRUE, plot = FALSE) {
  # Calculate voronoi partition for sites
  vor <- pts %>% 
    st_union() %>%
    st_voronoi(envelope = st_as_sfc(st_bbox(pts))) %>%
    # st_voronoi(envelope = st_as_sfc(st_bbox(areas))) %>%
    st_collection_extract(type = "POLYGON") %>% 
    st_sf() %>%  # convert from geometry set to simple feature collection
    identity()
  
  # order voronoi cells to match points
  vor <- vor[order(st_nearest_feature(st_centroid(vor), pts)), ]
  
  # cutoff voronoi cells from ocean, if desired (stops far away neighbours)
  # TODO: Cut via coastline
  if (cut_vor == TRUE) {
    # slightly smaller than areas bbox
    vor <- st_crop(
      vor, c("xmin" = -10, "ymin" = 51.6, "xmax" = -6, "ymax" = 55.2)
    )
  }
   
  # check that voronoi cells have been produced correctly
  if (plot == TRUE) {
    plot(st_geometry(areas)) 
    plot(vor, add = TRUE)
    plot(pts, col = "blue", pch = 16, add = TRUE)
    # test that voronoi cells and points ordered correctly
    # plot(vor[10, ], col = "blue", fill = NA, lwd = 5, add = TRUE)
    # plot(pts[10, ], col = "red", pch = 16, add = TRUE)
  }
  
  
  # calculate adjacency matrix from voronoi cells for stations
  return(spdep::nb2mat(spdep::poly2nb(vor), style = "B", zero.policy = TRUE))
}

# compute the total within-cluster sum of distances
# TODO: Create methods for the below functions to differ for PAM vs k-means
within_cluster_sum <- function(k, distance_matrix, fun = cluster::pam, ...) {
  clust_res <- fun(distance_matrix, k, ...)
  if (inherits(clust_res, "kmeans")) {
    return(clust_res$tot.withinss)
  } else if (inherits(clust_res, "pam")) {
    return(clust_res$objective[1])
  } else {
    stop("Clustering class not currently supported")
  }
}

# compute and produce scree plot
scree_plot <- \(dist_mat, k = 1:10, fun = cluster::pam, ...) {
  total_within_ss <- vapply(
    k, within_cluster_sum, dist_mat, fun = fun, ..., FUN.VALUE = numeric(1)
  )
  
  # scree plot
  plot(
    k, total_within_ss, type = "b", pch = 19 # , 
    # xlab = "Number of Clusters (k)", 
    # ylab = "Total Within-Cluster Sum of Distances",
    # main = "Scree Plot for K-medoids Clustering"
  )
  return(total_within_ss)
}

# plot clustering solution on map
plt_clust <- \(pts, clust_obj) {
  
  if (inherits(clust_obj, "kmeans")) {
    clust_element <- "cluster"  
    medoids <- NA
  } else if (inherits(clust_obj, "pam")) {
    clust_element <- "clustering"  
    medoids <- clust_obj$medoids
    if (inherits(medoids, "matrix")) medoids <- as.numeric(rownames(medoids))
  } else {
    stop("Clustering class not currently supported")
  }
  
  pts_plt <- cbind(pts, data.frame("clust" = clust_obj[[clust_element]])) %>% 
    mutate(
      row = row_number()
    ) %>% 
    mutate(mediod = ifelse(row %in% medoids, TRUE, FALSE))
  
  ggplot(areas) + 
    geom_sf(colour = "black", fill = NA) + 
    geom_sf(
      data = pts_plt, 
      aes(colour = factor(clust), shape = mediod, size = as.numeric(mediod)), 
      alpha = 0.8
    ) + 
    scale_shape_discrete(breaks = c(1, 15)) + 
    scale_size_continuous(range = c(3, 6)) +  
    guides(shape = "none", size = "none") + 
    theme
}

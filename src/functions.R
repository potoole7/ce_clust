#### Choosing K via AIC and TWGSS elbows ####

# function to perform likelihood ratio test for two CE model fits
lrt_test <- \(results_k, results_k1, df) {
  # Extract total log-likelihoods by summing all "ll" rows across all locations and variables
  sum_ll <- \(results) {
    # Sum "ll" across locations and vars
    sum(sapply(results, function(loc) {
      sum(sapply(loc, function(var_matrix) var_matrix["ll", ]))
    }))
  }

  ll_k <- sum_ll(results_k)
  ll_k1 <- sum_ll(results_k1)

  # Compute LRT statistic
  LRT_statistic <- -2 * (ll_k - ll_k1)

  # Compute p-value from chi-square distribution
  p_value <- 1 - pchisq(LRT_statistic, df = df)

  # Return results as a list
  return(list(LRT_statistic = LRT_statistic, df = df, p_value = p_value))
}


# function to find elbow in TWGSS plot
# From https://ieeexplore.ieee.org/document/5961514
find_elbow <- function(twgss) {
  n <- length(twgss)
  # Coordinates of the first and last points
  p1 <- c(1, twgss[1])
  p2 <- c(n, twgss[n])

  # Compute the distances from each point to the line (p1, p2)
  distances <- sapply(1:n, function(i) {
    p <- c(i, twgss[i])
    # Distance from point p to the line through p1 and p2
    abs((p2[2] - p1[2]) * p[1] - (p2[1] - p1[1]) * p[2] +
      p2[1] * p1[2] - p2[2] * p1[1]) /
      sqrt((p2[2] - p1[2])^2 + (p2[1] - p1[1])^2)
  })

  best_k <- which.max(distances)
  return(best_k)
}

# calculate AIC for a given CE fit
# TODO Add BIC
calc_inf <- \(ll, n_par, n = NULL, type = "AIC") {
  return(switch(type,
    "AIC" = (-2 * ll) + (2 * n_par) # ,
    # "BIC"  = (-2 * ll) + (log(n) * n_par)
  ))
}

# Function to calculate AIC from CE model output
ce_aic <- \(dep) {
  # pull LLs for each model
  lls <- vapply(dep, \(x) {
    vapply(x, `[`, "ll", , FUN.VALUE = numeric(1))
  }, numeric(length(dep[[1]])))
  # number of parameters = 4 * number of clusters * number of variables
  # TODO Will this number of larger for > 2 variables? for 3 it'll be 6!
  # TODO Test for 3 variables before putting in package or anything
  n_par <- 4 * length(dep) * length(dep[[1]])
  # n par for different AICs per model
  # n_par <- 4 * length(dep[[1]])
  # combined_ll
  ll <- sum(lls)
  # ll <- mean(lls)
  # return AIC
  return(calc_inf(ll, n_par))
  # calculate AIC for each model
  # return(vapply(lls, \(x) {
  #   calc_inf(x, n_par)
  # }, numeric(1)))
}

# fun to refit CE model with given clustering
fit_ce_clust <- \(clust_mem, data_mix, n_vars = 2, marg_prob, cond_prob) {
  clusts <- unique(clust_mem)
  # create new data_mix with cluster membership
  data_mix_clust <- lapply(clusts, \(x) {
    do.call(rbind, data_mix[clust_mem == x])
  })

  # refit model and return
  return(fit_ce(
    data_mix_clust,
    vars = paste0("col_", seq_len(n_vars)),
    # TODO Functionalise these inputs, may want to vary kl_prob!
    marg_prob = marg_prob,
    cond_prob = cond_prob,
    # f           = NULL, # fit models with ismev
    fit_no_keef = TRUE,
    output_all = FALSE
  ))
}

# fun to refit just CE model (no marginals)
fit_optim_clust <- \(
  clust_mem,
  data_mix,
  n_vars = 2,
  cond_prob,
  trans_fun = NULL,
  start_vals = c(0.01, 0.01),
  fixed_b = FALSE
) {
  clusts <- unique(clust_mem)

  # create new data_mix with cluster membership
  data_mix_clust <- lapply(clusts, \(x) {
    do.call(rbind, data_mix[clust_mem == x])
  })

  # transform from one margin to another, if desired
  if (!is.null(trans_fun)) {
    data_mix_clust <- lapply(data_mix_clust, trans_fun)
  }

  # refit model and return
  # return(fit_ce(
  #   data_mix_clust,
  #   vars = paste0("col_", seq_len(n_vars)),
  #   # TODO Functionalise these inputs, may want to vary kl_prob!
  #   marg_prob = marg_prob,
  #   cond_prob = cond_prob,
  #   # f           = NULL, # fit models with ismev
  #   fit_no_keef = TRUE,
  #   output_all = FALSE
  # ))
  return(lapply(seq_along(data_mix_clust), \(i) {
    if (is.list(start_vals)) {
      start_vals_spec <- start_vals[[i]]
    } else {
      start_vals_spec <- start_vals
    }

    o <- ce_optim(
      Y = data_mix_clust[[i]],
      dqu = cond_prob,
      control = list(maxit = 1e6),
      constrain = FALSE,
      # TODO Should I change this to what initial guess should be from theory??
      # start     = c(0.01, 0.01)
      start = start_vals_spec,
      fixed_b = fixed_b
    )
  }))
}


#### Simulation plots  ####

# Function to add common plot items to each plot of adjusted rand index
common_plot_items <- \(
  p,
  xlab = expression(rho["Gauss"]),
  x_sec_lab = expression(rho[t[1]]),
  y_sec_lab = expression(rho[t[2]]),
  facet_form = (cor_t1 ~ cor_t2)
) {
  p +
    # facet_grid(cor_t1 ~ cor_t2) +
    facet_grid(facet_form) +
    scale_x_continuous(
      sec.axis = sec_axis(
        ~.,
        # name = expression(rho[t[1]]),
        name   = x_sec_lab,
        breaks = NULL,
        labels = NULL
      )
    ) +
    scale_y_continuous(
      sec.axis = sec_axis(
        ~.,
        # name = expression(rho[t[2]]),
        name   = y_sec_lab,
        breaks = NULL,
        labels = NULL
      ),
      breaks = seq(0, 1, by = 0.2)
    ) +
    labs(
      # x = expression(rho["Gauss"]),
      # x = TeX("$\\rho_{t_1}$"),
      x = xlab,
      # y = "Adjusted Rand Index",
      y = "ARI",
      fill = "",
      colour = ""
    ) +
    evc_theme() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 0.5),
      # axis.title = element_markdown(size = 12, family = "serif")
    ) +
    ggsci::scale_fill_nejm() +
    ggsci::scale_colour_nejm()
}


#### Bootstrapping ####

# function to perform bootstrapping for CE model
# TODO Move to package
# TODO Have work for > 2 variables
boot_ce <- \(fit, R = 100, trace = 10, ncores = 1) {
  # check if fit has correct structure
  stopifnot(
    "fit must be an object returned by fit_ce with output_all = TRUE" =
      all(c("marginal", "arg_vals") %in% names(fit))
  )

  # Temp: only works for ismev fit currently
  # stopifnot("Currently only supported for ")

  # pull data
  arg_vals <- fit$arg_vals # original arguments to fit_ce object
  marginal <- fit$marginal # marginal GPDs
  orig <- fit$original # original data
  transformed <- fit$transformed # Laplace transformed data
  dependence <- fit$dependence # dependence parameters

  # extract marginal and dependence thresholds
  # TODO Implement for multiple locations
  # thresh_dep <- lapply(dependence, \(x) {
  #   vapply(x, \(y) {
  #     y[rownames(y) == "dth", ]
  #   }, numeric(length(marginal[[1]]) - 1))[1, ]
  # })
  thresh_dep <- lapply(dependence, \(x) {
    res <- vapply(x, \(y) {
      y[rownames(y) == "dth", , drop = FALSE]
    }, numeric(length(x) - 1))
    if (!is.matrix(res)) {
      res <- matrix(res, nrow = 1)
      colnames(res) <- names(x)
    }
    res[1, , drop = FALSE]
  })

  # TODO Temp: Using same thresh for each location, need to expand!
  thresh_marg <- lapply(marginal, \(x) {
    vapply(x, `[[`, "thresh", FUN.VALUE = numeric(1))
  })
  thresh_marg <- thresh_marg[[1]]

  # Parallel setup
  # TODO Implement below, not using now
  apply_fun <- ifelse(ncores == 1, lapply, parallel::mclapply)
  ext_args <- NULL
  if (ncores > 1) {
    ext_args <- list(mc.cores = ncores)
  }
  loop_fun <- \(...) {
    do.call(apply_fun, c(list(...), ext_args))
  }

  # Pull start values
  # TODO Make `coef` method for dependence, would be much easier than this loop!
  start <- lapply(dependence, \(x) {
    lapply(x, \(y) {
      # scale back towards zero in case point est on edge of original parameter
      # space and falls off edge of constrained space for bootstrap sample
      y[c("a", "b"), , drop = FALSE] * 0.75
    })
  })

  # TODO Temp: remove
  # marg_loc  <- marginal[[1]]
  # trans_loc <- transformed[[1]]
  # dep_loc   <- dependence[[1]]
  # orig_loc <- orig[[1]]
  # cond_var <- "O3"
  # n_pass <- 3
  # Function to prepare bootstrapped data for each location and conditioned var
  prep_boot_loc_dat <- \(
    marg_loc, trans_loc, dep_loc, orig_loc, thresh_dep, n_pass = 3
    # marg_loc, trans_loc, dep_loc, orig_loc, cond_var, thresh_dep, n_pass = 3
  ) {
    # pull bootstrap sample
    # (must be different for each loc as nrows may differ)
    indices <- sample(seq_len(nrow(trans_loc)), replace = TRUE)
    trans_loc_boot <- trans_loc[indices, ]

    # TODO test this to see if it ever fails?
    stopifnot(names(thresh_dep) == colnames(trans_loc))

    # Reorder bootstrap sample to have the same order as original data
    test <- FALSE
    # which <- which(names(marg_loc) %in% cond_var)
    while (test == FALSE) {
      # for (j in 1:(dim(g)[[2]])){ # loop across columns
      for (j in seq_along(marg_loc)) {
        # replace ordered Yi with ordered sample from standard Laplace CDF
        u <- matrix(runif(nrow(trans_loc_boot)))
        trans_loc_boot[
          order(trans_loc_boot[, j]), j
        ] <- sort(laplace_trans(u)[, 1])
      }
      # need exceedance in cond. var, and also in other vars for these rows
      # TODO I fit CE for every variable, not a specific one; how to handle?
      # TODO Is this correct???
      # if (sum(trans_loc_boot[, which] > thresh_dep[which]) > 1 &&
      #     all(trans_loc_boot[trans_loc_boot[, which] > thresh_dep[which], which] > 0)) {
      #   test <- TRUE
      # }
      if (any(colSums(sweep(trans_loc_boot, 2, thresh_dep, FUN = ">")) > 0)) {
        test <- TRUE
      }
    }

    # convert variables to original scale
    # TODO Check if this is correct! Look at plots for texmex version
    orig_loc_boot <- inv_semi_par_cdf(
      # inverse Laplace transform to CDF
      inv_laplace_trans(trans_loc_boot),
      # original data for where semiparametric thresholding occurs
      select(orig_loc, -matches("name")),
      marg_loc # marginal GPD parameters
    )
    colnames(orig_loc_boot) <- names(marg_loc)

    # test for no marg exceedances over sampled points, if so resample w/ nPass
    max_vals <- apply(orig_loc_boot, 2, max, na.rm = TRUE)
    marg_thresh <- vapply(marg_loc, `[[`, "thresh", FUN.VALUE = numeric(1))
    if (!all(max_vals > marg_thresh)) {
      return(list(NA))
    }
    return(orig_loc_boot)
  }

  # perform bootstrapping
  # TODO Allow parallel computation
  # boot_fits <- lapply(seq_len(R), \(i) {
  boot_fits <- loop_fun(seq_len(R), \(i) {
    if (i %% trace == 0) {
      system(sprintf("echo %s", paste(i, "replicates done")))
    }
    dat_boot <- lapply(seq_along(marginal), \(j) { # loop through locations
      # prepare bootstrapped data for location j
      dat_spec <- prep_boot_loc_dat(
        marginal[[j]],
        transformed[[j]],
        dependence[[j]],
        orig[[j]],
        thresh_dep[[j]]
      )
      # check if NA, if so then rerun
      if (all(is.na(dat_spec)) && n_pass > 1) {
        for (i in seq_len(n_pass - 1)) {
          dat_spec <- prep_boot_loc_dat(
            marginal[[j]], transformed[[j]], dependence[[j]]
          )
          if (!all(is.na(dat_spec))) {
            break
          }
        }
      } else if (all(is.na(dat_spec))) {
        stop(paste0(
          "Failed to generate bootstrapped data after ", n_pass, " attempts"
        ))
      }
      # label data with locations
      as_tibble(dat_spec) |>
        mutate(name = names(marginal)[j])
    })
    names(dat_boot) <- names(marginal)

    # refit CE model using same parameters as before, force output_all = TRUE
    fit_boot <- do.call(
      fit_ce,
      c(
        # TODO Investigate why this won't work without bind_rows
        # TODO Add fixed values for dependence and marginal thresholds
        list(
          data = bind_rows(dat_boot),
          start = start,
          # start = c("a" = 0.01, "b" = 0.01),
          # TODO Implement arg_val for different locations in fit_ce
          marg_val = thresh_marg,
          marg_prob = NULL
        ),
        arg_vals[names(arg_vals) != "marg_prob"],
        output_all = TRUE
      )
    )
    return(list(
      # extract a and b, all that's important
      "dependence" = lapply(fit_boot$dependence, \(x) {
        lapply(x, \(y) {
          y[c("a", "b"), ]
        })
      }),
      # extract xi and sigma
      "marginal" = lapply(fit_boot$marginal, \(x) {
        lapply(x, \(y) unlist(y[!names(y) == "thresh"]))
      })
    ))
  })

  # extract a and b parameters for dependence and xi and sigma for marginal
  # into nice structure:

  # for marginal, want df with locs and vars and sigma and xi estimates for each

  # pull margins for each bootstrap sample
  marg_out <- lapply(boot_fits, `[[`, "marginal") |>
    # Extract location and parameters from each list element
    # for each bootstrap margin
    purrr::map_df(~ {
      # loop through locations
      purrr::map_dfr(.x, \(loc_data) {
        # Create a tidy data frame for each location's parameters, keeping location as it is
        loc_data |>
          # loop through variables, extract parameters and name accordingly
          purrr::map_dfr(~ tibble(
            parameter = names(.x),
            value = .x
          ), .id = "vars")
      }, .id = "name")
    })

  # TODO Describe dep_out, reduce Chat GPT comments
  # TODO Fix for > 1 location
  # extract dependence from each bootstrap sample
  dep_out <- lapply(boot_fits, `[[`, "dependence") |>
    # lapply(\(boot_sample) {
    loop_fun(\(boot_sample) {
      # boot_sample <- lapply(boot_fits, `[[`, "dependence")[[1]]
      # again, loop through locations for a given bootstrap sample
      purrr::map_dfr(names(boot_sample), function(loc_name) {
        # loc <- boot_sample[[1]]
        loc <- boot_sample[[loc_name]]

        # loop over variables
        purrr::map_dfr(names(loc), function(cond_var) {
          # Get the matrix for the current conditioning variable
          var_mat <- loc[[cond_var]]

          # be careful if var_mat is not a matrix
          if (is.matrix(var_mat)) {
            var_names <- colnames(var_mat)
          } else {
            var_names <- names(loc)[!names(loc) == cond_var]
            var_mat <- as.matrix(var_mat)
          }

          # Create a tidy data frame for this location's matrix
          tibble(
            parameter = rep(c("a", "b"), each = ncol(var_mat)),
            # vars      = rep(colnames(var_mat), times = 2),
            vars      = rep(var_names, times = 2),
            value     = c(var_mat["a", ], var_mat["b", ]),
            cond_var  = rep(cond_var, ncol(var_mat) * 2),
            name      = rep(loc_name, ncol(var_mat) * 2)
          )
        })
      })
    })

  # Combine all the lists of data frames into one final data frame
  dep_out <- bind_rows(dep_out)

  return(list("marginal" = marg_out, "dependence" = dep_out))
}





#### Sensitivity analysis ####

# calculate local rand index
# (from https://github.com/Xuanjie-Shao/NonStaExtDep/blob/main/Functions/Utils.R)
local_rand_index <- \(P0, P1) {
  # for each partition P1, true partition P0
  if (length(P0) == length(P1)) {
    D <- length(P0)
  } else {
    return(0)
  }
  beta.all <- vector(length = D)
  for (i in seq_len(D)) {
    betai <- num_ss <- num_dd <- 0
    for (j in seq_len(D)[-i]) {
      if (P0[i] == P0[j] && P1[i] == P1[j]) {
        num_ss <- num_ss + 1
      } else if (P0[i] != P0[j] && P1[i] != P1[j]) {
        num_dd <- num_dd + 1
      }
    }
    betai <- num_ss + num_dd
    beta.all[i] <- betai
  }
  return(beta.all / D)
}

# TODO: Add shared plotting functions

# function to calculate confidence interval for sensitivity analysis results
# Calculated as credible interval, rather than from CLT and normality assumption
summarise_sens_res <- \(results_grid, conf_level = 0.95) {
  # z-score corresponding to confidence level
  # z <- qnorm((1 + conf_level) / 2)

  results_grid %>%
    # remove any previous summaries as they mess up grouping
    dplyr::select(-(ends_with("_rand") & !matches("adj_rand"))) %>%
    relocate(adj_rand, .after = everything()) %>% # ensure is last col
    group_by(across(1:last_col(1))) %>%
    summarise(
      mean_rand    = mean(adj_rand, na.rm = TRUE),
      # need median to give central estimate for confidence interval!
      median_rand  = median(adj_rand, na.rm = TRUE),
      # # (z score corresponding to confidence level) * standard error
      # marg_rand  = z * (sd(adj_rand, na.rm = TRUE) / sqrt(n())),
      # # CIs
      # lower_rand = mean_rand - marg_rand,
      # upper_rand = mean_rand + marg_rand,
      lower_rand   = quantile(adj_rand, (1 - conf_level) / 2, na.rm = TRUE),
      upper_rand   = quantile(adj_rand, 1 - (1 - conf_level) / 2, na.rm = TRUE),
      .groups      = "drop"
    ) %>%
    # dplyr::select(-marg_rand) %>% # not required for plotting
    return()
}


#### Plotting functions ####

# plot clustering vs dist2coast and altitude
covar_plt <- \(pts, clust_obj, data_full, areas) {
  plt_dat <- pts %>%
    st_drop_geometry() %>%
    mutate(clust = clust_obj$clustering) %>%
    left_join(
      data_full %>%
        # average wind direction across locations
        group_by(name) %>%
        mutate(
          wind_dir = as.numeric(mean(circular(
            wind_dir,
            type = "angles", units = "degrees", modulo = "2pi"
          )))
        ) %>%
        ungroup() %>%
        select(name, dist2coast, wind_dir, alt) %>%
        distinct()
    ) %>%
    # pivot_longer(dist2coast:alt, names_to = "var", values_to = "value") %>%
    identity()

  # histogram of different variables for different clusters
  # TODO: add uncertainty bounds based on variability
  p1 <- plt_dat %>%
    group_by(clust) %>%
    # summarise(value = mean(value), .groups = "drop") %>%
    summarise(
      across(c("dist2coast", "alt"), mean),
      wind_dir = as.numeric(mean(circular(
        wind_dir,
        type = "angles", units = "degrees", modulo = "2pi"
      ))),
      .groups = "drop"
    ) %>%
    pivot_longer(dist2coast:wind_dir, names_to = "var", values_to = "value") %>%
    ggplot() +
    geom_bar(
      aes(x = factor(clust), y = value, fill = factor(clust)),
      stat = "identity"
    ) +
    facet_wrap(~var, scales = "free") +
    labs(x = "Cluster", fill = "Cluster") +
    ggsci::scale_fill_nejm() +
    evc_theme()

  # plot density of different variables for both clusters
  p2 <- plt_dat %>%
    pivot_longer(dist2coast:alt, names_to = "var", values_to = "value") %>%
    ggplot() +
    geom_density(
      aes(x = value, fill = factor(clust)),
      alpha = 0.5
    ) +
    facet_wrap(~var, scales = "free") +
    evc_theme() +
    labs(fill = "Cluster") +
    guides(fill = guide_legend(override.aes = list(alpha = 1))) +
    ggsci::scale_fill_nejm()

  return(list("hist" = p1, "density" = p2))
}


#### Simulation functions ####

# Function to generate copula data
sim_cop_dat <- \(
  n_locs = 12, # number of locations
  n_vars = 2, # number of variables at each location
  n = 1e4, # number of simulations for each variable
  n_clust = 2, # number of clusters
  cluster_mem = NULL, # desired membership, if NULL then evenly split
  cor_gauss = NULL, # bulk correlation for n_clust clusters from Gaussian copula
  # params_norm,         # normal marginal parameters (same for both)
  cor_t = NULL, # extreme correlation for n_clust clusters from t-copula
  df_t = NULL, # degrees of freedom for t-copula
  # params_gpd = NULL, # GPD margin parameters
  mix_p = c(0.5, 0.5), # mixture weights
  perturb_cor = FALSE, # perturb correlation for each location within clusters
  perturb_val = 0.05, # value to perturb by if desired
  quiet = FALSE, # do return messages
  qfun = evd:::qgpd,
  qargs = NULL
) {
  # many arguments must be of length n_clust
  # stopifnot(all(vapply(
  #   # list(cor_gauss, cor_t, df_t, params_gpd, mix_p),
  #   list(cor_gauss, cor_t, df_t),
  #   \(x) length(x) == n_clust, logical(1)
  # )))

  # checks for mixture percentages
  stopifnot(is.numeric(mix_p) && length(mix_p) == 2 && sum(mix_p) == 1)

  # if (!quiet && is.null(params_gpd)) {
  #   message("No GPD parameters supplied, assuming Laplace margins")
  # }

  # Simulate from Gaussian Copula with GPD margins
  gen_gauss_cop <- \() {
    gauss_cop <- lapply(seq_len(n_locs), \(i) {
      # pull correlation specified for each cluster
      # if membership unspecified, assign equal membership to each cluster
      if (is.null(cluster_mem)) {
        group <- ceiling(i / (n_locs / n_clust))
      } else {
        group <- cluster_mem[i]
      }
      # Assign the corresponding value to the result
      cor <- cor_gauss[group]

      # set correlation matrix/vector with correct dimensions
      cor <- rep(cor, n_vars * (n_vars - 1) / 2)

      # create (Gaussian) copula object
      cop_norm <- copula::normalCopula(cor, dim = n_vars, dispstr = "un")
      # simulate uniform draws from copula
      u <- copula::rCopula(n, cop_norm)
      # transform to GPD margins, if GPD parameters specified
      # if (!is.null(params_gpd)) {
      #   evd::qgpd(
      #     p     = u,
      #     loc   = 0,
      #     scale = params_gpd[1],
      #     shape = params_gpd[2]
      #   )
      #   # Otherwise, assume Laplace margins (i.e. we skip GPD and go to CE)
      # } else {
      #   # TODO Change to be called qlaplace
      #   # TODO This could just be an rlaplace/rmvlaplace function!
      #   evc:::inv_laplace_trans(u)
      # }
      # use supplied quantile function; defaults to default parameters
      if (is.null(qargs)) {
        return(do.call(qfun, list(u)))
      } else {
        return(do.call(qfun, list(u, qargs)))
      }
    })
    return(gauss_cop)
  }

  if (mix_p[[1]] == 0) {
    if (!quiet) {
      message("Data from t-copula only")
    }
    gauss_cop <- NULL
  } else {
    gauss_cop <- gen_gauss_cop()
  }

  # simulate from t-Copula with GPD margins
  # TODO Functionalise repeated code, lierally the same code bar supplying df!
  gen_t_cop <- \() {
    t_cop <- lapply(seq_len(n_locs), \(i) {
      if (is.null(cluster_mem)) {
        group <- ceiling(i / (n_locs / n_clust))
      } else {
        group <- cluster_mem[i]
      }
      # Assign the corresponding value to the result
      cor <- cor_t[group]
      df <- df_t[group]

      # optionally perturb correlation for each location within clusters
      if (perturb_cor) {
        cor <- cor + runif(1, -perturb_val, perturb_val)
        cor <- pmax(pmin(cor, 1), 0) # ensure within [0, 1]
      }

      cor <- rep(cor, n_vars * (n_vars - 1) / 2)
      cop_t <- copula::tCopula(cor, dim = n_vars, df = df, dispstr = "un")
      u <- copula::rCopula(n, cop_t)
      # if (!is.null(params_gpd)) {
      #   evd::qgpd(
      #     p     = u,
      #     loc   = 0,
      #     scale = params_gpd[1],
      #     shape = params_gpd[2]
      #   )
      #   # Otherwise, assume Laplace margins (i.e. we skip GPD and go to CE)
      # } else {
      #   message("No GPD parameters supplied, assuming Laplace margins")
      #   evc:::inv_laplace_trans(u)
      # }
      if (is.null(qargs)) {
        return(do.call(qfun, list(u)))
      } else {
        return(do.call(qfun, list(u, qargs)))
      }
    })
    return(t_cop)
  }


  if (mix_p[[2]] == 0) {
    if (!quiet) {
      message("Data from Gaussian copula only")
    }
    t_cop <- NULL
  } else {
    t_cop <- gen_t_cop()
  }

  # Return margins of t-copula if Gaussian copula is NULL, and vice versa
  if (is.null(gauss_cop)) {
    data_mix <- t_cop
  } else if (is.null(t_cop)) {
    data_mix <- gauss_cop
    # Return mixture of Gaussian and t-copulas if desired
  } else {
    data_mix <- lapply(seq_len(n_locs), \(i) {
      x <- nrow(gauss_cop[[i]])
      y <- nrow(t_cop[[i]])
      # sample mix_p * nrow(gauss_cop) rows from gauss_cop (same for t_cop)
      rbind(
        gauss_cop[[i]][sample(seq_len(x), size = mix_p[[1]] * x), ],
        t_cop[[i]][sample(seq_len(y), size = mix_p[[2]] * y), ]
      )
    })
  }
  return(list(
    "gauss_cop" = gauss_cop,
    "t_cop"     = t_cop,
    "data_mix"  = data_mix
  ))
}


#### handy sf functions ####

# convert df to sf
st_to_sf <- function(dat, coords = c("lon", "lat"), crs = "WGS84", ...) {
  out <- dat %>%
    st_as_sf(coords = coords, ...)
  st_crs(out) <- crs
  return(out)
}

# calculate distance to coast
dist2coast <- \(dat, areas, convert_to_coastline = TRUE) {
  coastline <- areas
  # convert areas to coastline, if required
  if (convert_to_coastline) {
    coastline <- st_union(coastline)
    coastline <- st_cast(coastline, "MULTILINESTRING")
  }

  # locations to find distance to coast for
  locs <- data %>%
    distinct(name, lon, lat) %>%
    st_to_sf()

  distances <- as.numeric(st_distance(locs, coastline))

  data %>%
    left_join(data.frame("name" = locs$name, "dist2coast" = distances)) %>%
    return()
}


#### Load netcdf data ####

# load netcdf and return as data.table object
load_netcdf <- function(file) {
  # open connection to file
  nc_connection <- nc_open(file)
  # load all data from netcdf file
  data <- lapply(nc_connection$var, function(x) {
    ncvar_get(nc_connection, x)
  })
  gc()
  names(data) <- names(nc_connection$var)

  data <- data[!names(data) %in% c("number", "expver")]


  # data[[1]][, , 1] %>%
  #   as.matrix() %>%
  #   as.data.frame() %>%
  #   mutate(row = row_number()) %>%
  #   pivot_longer(contains("V")) %>%
  #   mutate(name = as.numeric(stringr::str_remove(name, "V"))) %>%
  #   ggplot() +
  #   geom_raster(aes(x = row, y = name, fill = value))

  # create dummy data.table to store data in
  lon_vals <- nc_connection$dim$longitude$vals
  lat_vals <- nc_connection$dim$latitude$vals
  if (is.null(lon_vals)) {
    lon_vals <- nc_connection$dim$lon$vals
    lat_vals <- nc_connection$dim$lat$vals
  }
  if ("valid_time" %in% names(nc_connection$dim)) {
    nc_connection$dim$time <- nc_connection$dim$valid_time
  }
  time_vals <- nc_connection$dim$time$vals
  # ret <- data.table(
  #   "longitude" = rep(lon_vals, each = length(lat_vals) * length(time_vals)),
  #   "latitude"  = rep(
  #     rep(lat_vals, each = length(time_vals)), times = length(lon_vals)
  #   ),
  #   "time"      = rep(time_vals, times = length(lon_vals) * length(lat_vals)),
  #   "u10"       = NA,
  #   "v10"       = NA,
  #   "msl"       = NA
  # )
  # gc()

  # only add time vals if available
  vals <- list(lon_vals, lat_vals)
  names <- c("longitude", "latitude")
  if (!is.null(time_vals)) {
    vals <- c(vals, list(time_vals))
    names <- c(names, "time")
  }
  ret <- lapply(seq_along(data), \(i) {
    dimnames(data[[i]]) <- vals
    melt(
      data[[i]],
      varnames = names,
      value.name = names(data)[i]
    )
  })
  rm(data)
  gc()

  # (specific to wind speed) Join list variables
  if ("u10" %in% as.vector(vapply(ret, names, character(ncol(ret[[1]]))))) {
    ret <- ret[[1]] %>%
      bind_cols(select(ret[[2]], "v10")) %>%
      # bind_cols(select(ret[[3]], "msl")) %>%
      # filter(time == 1) %>%
      # filter(time == 920424) %>%
      # ggplot() +
      # geom_raster(aes(x = longitude, y = latitude, fill = u10))
      identity()
  } else {
    ret <- ret[[1]]
  }

  # add data to this data.table
  # for (i in seq_len(length(data))) {
  # ret[, (3 + i)] <- as.vector(matrix(
  # data[[i]], nrow = dim(data[[i]])[3], byrow = TRUE
  # ))
  # data[[i]] <- NA
  # gc()
  # }
  # rm(data); gc()

  # preprocess time
  # time since this date
  if (!is.null(time_vals)) {
    units <- nc_connection$dim$time$units
    reference_date <- as.POSIXct(
      # substr(nc_connection$dim$time$units, 13, 31), tz = "UTC"
      gsub(".*?(\\d{4}-\\d{2}-\\d{2}).*", "\\1", units),
      tz = "UTC"
    )
    if (grepl("seconds", units)) {
      ret$time <- reference_date + ret$time
    } else if (grepl("hours", units)) {
      ret$time <- reference_date + ret$time * 3600
    }
  }

  ret <- setDT(ret)

  # close connection
  nc_close(nc_connection)

  return(ret)
}

# create matrix, as required by extremis
create_matrix <- function(dat, col) {
  ret <- dat %>%
    select(name, date, !!col) %>%
    pivot_wider(names_from = name, values_from = !!col, id_cols = date) %>%
    drop_na()

  dates <- unique(ret$date)

  ret <- ret %>%
    select(-date) %>%
    as.matrix()
  rownames(ret) <- as.character(dates)
  return(ret)
}

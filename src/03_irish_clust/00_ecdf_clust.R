#### Cluster Conditional Extremes using ECDF only (no marginal GPD model) ####

# Don't fit marginal model, just use all ECDF (and PIT to convert to Laplace)!
# See if CE model fit is improved for variety of scenarios:
# - Increased dependence threshold for wind ()
# - Fix alpha to be positive
# - Fix beta (to be positive)???
# - Bootstrapped starting values for alpha, beta ()
# Also (i) clustering and (ii) refitting to see if bootstrap CIs improve!

# TODO Split this file into multiple files/scripts, too much going on here!

# TODO Could use bootstrapped estimates as start values?
# TODO Can we also produce profile likelihood plots??

# TODO Fix adjacency stipulation!!! (Can't really do with k-medoids)

# For meeting:
# TODO Facet clustering solution for various DQU
# Plot bootstrapped uncertainty before/after clustering (done)

# TODO Final diagnostics for chosen models!
# TODO May have to choose number of clusters after Laplace truncation point?

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

# parameters for MC estimation of JSGa
n_mc <- 1000 # number of MC samples
mc_method <- "laplace_trunc2" # truncate Laplace using empirical dist quant
laplace_cap <- 0.99 # Take 99th quantile (for choice, see `002_emp_laplace_trunc.R)``


#### Misc Functions ####

# get raster of elevations from areas file
get_elev_rast <- \(areas, z = 9, bins = NULL, labels = NULL) {
  # convert to projected CRS
  areas_proj <- st_transform(areas, IrishGrid = 29902)

  # pull elevation data from Open Elevation API (DEM = Digital Elevation Model)
  dem <- terra::rast(get_elev_raster(
    locations = areas_proj, z = z, clip = "locations"
  ))

  # Mask DEM to the exact MULTIPOLYGON footprint
  dem_masked <- terra::mask(dem, vect(areas_proj))

  # Convert raster to a data.frame for ggplot
  df_dem <- as.data.frame(dem_masked, xy = TRUE, na.rm = TRUE)
  names(df_dem) <- c("x", "y", "elevation")
  # filter out negative elevations
  df_dem <- filter(df_dem, elevation >= 0)

  if (!is.null(bins)) {
    df_dem <- df_dem %>%
      mutate(
        elev_bin = factor(cut(
          elevation,
          # breaks = c(seq(0, 500, by = 100), Inf),
          breaks = bins,
          # labels = c(
          #   "200–300 m", "300–400 m", "400–500 m", "> 500 m"
          # ),
          labels = labels,
          right  = FALSE
        ), levels = labels)
      )
  }

  return(df_dem)
}

# pull parameters from dependence object
get_pars <- \(ce_fit) lapply(ce_fit, \(x) lapply(x, `[[`, "params"))

# function to remove " - " (county name) from names
rm_cnty <- \(x) {
  vapply(stringr::str_split(x, " - "), `[[`, 1, FUN.VALUE = character(1))
}

# function to convert to ECDF
calc_ecdf <- \(x) {
  apply(x, 2, \(dat_spec) {
    dat_spec_ord <- order(dat_spec)
    dat_spec_sort <- dat_spec[dat_spec_ord]

    # calculate ECDF
    m <- length(dat_spec)
    ecdf_vals <- (seq_len(m)) / (m + 1)
    # convert back to original order
    ecdf_dat_ord <- numeric(m)
    ecdf_dat_ord[dat_spec_ord] <- ecdf_vals
    ecdf_dat_ord
  })
}

# perform PIT on ECDF to get data on Laplace margins
trans_fun <- \(x) {
  laplace_trans(calc_ecdf(x))
}

# pull dependence parameters and residuals out separately
pull_element <- \(x, element) {
  if (element %in% names(x)) {
    return(x[[element]])
  } else {
    lapply(x, pull_element, element)
  }
}

# pull dependence paramaters into df
dep_to_df <- \(dep) {
  # pull a and b values for each location and variable
  ab_vals <- lapply(dep, \(x) {
    list(x[[1]][1:2], x[[2]][1:2])
  })

  # convert to dataframe
  bind_rows(lapply(seq_along(ab_vals), \(i) {
    data.frame(
      "name" = names(ab_vals)[i],
      "vars" = c("rain", "rain", "wind_speed", "wind_speed"),
      "parameter" = c("a", "b", "a", "b"),
      "value" = c(
        ab_vals[[i]][[1]][1], ab_vals[[i]][[1]][2],
        ab_vals[[i]][[2]][1], ab_vals[[i]][[2]][2]
      )
    )
  })) |>
    # TODO: Investigate NAs here
    filter(!is.na(name))
}

# function to calculate adjacency matrix (slightly different to evc version)
calc_adj_mat2 <- \(pts, cut_vor = TRUE, plot = FALSE, areas = NULL) {
  # pts_proj <- sf::st_transform(pts, 2157)

  # # Calculate voronoi partition for sites
  # vor <- pts_proj |>
  #   sf::st_union() |>
  #   sf::st_voronoi(envelope = sf::st_as_sfc(sf::st_bbox(pts_proj))) |>
  #   sf::st_collection_extract(type = "POLYGON") |>
  #   sf::st_sf() |> # convert from geometry set to simple feature collection
  #   identity()

  # # order voronoi cells to match points
  # # vor <- vor[order(sf::st_nearest_feature(sf::st_centroid(vor), pts_proj)), ]
  # order_index <- order(sf::st_nearest_feature(sf::st_centroid(vor), pts_proj))
  # vor <- vor[order_index, ]

  # # Use the same ordering to get names from pts
  # pts_names <- pts_proj$name[order_index]

  # # cutoff voronoi cells from ocean, if desired (stops far away neighbours)
  # # TODO: Generalise more! Functionalise box below
  # if (cut_vor == TRUE) {
  #   # slightly smaller than areas bbox
  #   # vor <- sf::st_crop(
  #   #   vor, c("xmin" = 95000, ymin = 177000, xmax = 1110000, ymax = 129000)
  #   # )
  # }

  # # check that Voronoi cells have been produced correctly
  # if (plot == TRUE && !is.null(areas)) {
  #   plot(sf::st_geometry(areas))
  #   # plot(vor, add = TRUE)
  #   plot(sf::st_crop(
  #     sf::st_transform(
  #       vor,
  #       "WGS84"
  #     ),
  #     c("xmin" = -10, "ymin" = 51.6, "xmax" = -5.7, "ymax" = 55.2)
  #   ), add = TRUE)
  #   plot(pts, col = "blue", pch = 16, add = TRUE)
  # }

  # # calculate adjacency matrix from voronoi cells for stations
  # adj_mat <- spdep::nb2mat(spdep::poly2nb(vor), style = "B", zero.policy = TRUE)

  # rownames(adj_mat) <- colnames(adj_mat) <- pts_names


  # Compute Voronoi tessellation
  coords <- st_coordinates(pts)
  voronoi <- deldir::deldir(coords[, 1], coords[, 2])

  # Extract neighbours (sharing Voronoi edges)
  nb <- deldir::tile.list(voronoi) # list of Voronoi tiles
  neighbours <- spdep::nb2listw(
    spdep::cell2nb(length(nb), 1),
    style = "B"
  )$neighbours

  # Create adjacency matrix from neighbours
  adj_mat <- spdep::nb2mat(neighbours, style = "B", zero.policy = TRUE)
  # set diagonal to 1
  diag(adj_mat) <- 1

  return(adj_mat)
}

#### Plotting functions ####

# check residuals for single location
plot_resid <- \(spec_resid) {
  vars <- names(spec_resid)
  plots <- lapply(vars, \(var) {
    Z <- spec_resid[[var]]

    # For plotting, tidy up variable names
    cond_var <- stringr::str_replace_all(var, "_", " ")
    lhs_var <- stringr::str_replace_all(colnames(Z), "_", " ")

    if (all(is.na(Z))) {
      return(NA)
    }

    n <- length(Z)

    dqu_spec <- dqu
    if (length(dqu) > 1) {
      dqu_spec <- dqu[vars == var]
    }
    p <- seq(dqu_spec, 1 - (1 / n), length = n)

    data.frame(p, "resid" = Z) |>
      ggplot(aes(x = p, y = Z)) +
      geom_point(alpha = 0.7) +
      geom_smooth() +
      theme +
      labs(
        x = paste0("F(", cond_var, ")"),
        # y = paste0("Z ", colnames(Z), " | ", cond_var)
        y = paste0("Z ", lhs_var, " | ", cond_var)
      )
  })
  names(plots) <- vars
  return(plots)
}

# produce residual plots for all locations
resid_plot <- \(fit) {
  # df to join site names with counties (easier to know where they are)
  residuals <- fit$residual
  names_df <- data.frame("name" = names(residuals)) |>
    left_join(county_key_df)

  resid_plots <- lapply(seq_along(residuals), \(i) {
    plots <- plot_resid(fit$residual[[i]])
    # TODO Add parameter estimates
    plots_spec <- lapply(seq_along(plots), \(j){
      print(j)
      if (all(is.na(plots[[j]]))) {
        return(NA)
      }
      pars <- fit$dependence[[i]][[j]][c("a", "b"), ]
      plots[[j]] +
        # ggtitle(paste(
        #   names_df$name[i],
        #   names_df$county[i],
        #   sep = " - "
        # ))
        ggtitle(
          paste0("a = ", round(pars[1], 3), ", b = ", round(pars[2], 3))
        )
    })
    # join plots
    wrap_plots(plots_spec) +
      plot_annotation(
        title = paste0(
          names_df$name[i],
          " (",
          names_df$county[i],
          ")"
        ),
        theme = theme
      )
  })
  names(resid_plots) <- names_df$name

  return(resid_plots)
}

# plot quantiles of conditional expectation at single location (for single var)
plot_ce_quantiles <- \(
  dep_fit, spec_loc, cond_var = "rain", quantiles = seq(0.1, by = 0.2, len = 5)
) {
  # TODO Functionalise this somehow??
  # take out data for one location
  dep_fit_spec <- list(
    "original" = dep_fit$original[[spec_loc]],
    "residual" = dep_fit$residual[[spec_loc]],
    "dependence" = dep_fit$dependence[[spec_loc]],
    "transformed" = dep_fit$transformed[[spec_loc]]
  )

  # non-conditioning/other variables
  vars <- names(dep_fit_spec$residual)
  vars <- vars[vars != cond_var]
  n <- nrow(dep_fit_spec$residual[[cond_var]])

  # dependence parameters for conditioning variables
  dep <- dep_fit_spec$dependence[[cond_var]]
  # dependence thresholds and quantiles (on laplace scale)
  dqu <- dep_fit$arg_vals$cond_prob
  if (length(dqu) > 1) {
    dqu <- dqu[[cond_var]]
  }
  # dth <- quantile(dep_fit_spec$original[[cond_var]], dqu)
  dth <- dep_fit_spec$dependence[[cond_var]]["dth", ]
  # Determine x-axis values to estimate CE quantiles at along conditioned var
  # TODO Change to Laplace scale
  # xmax <- max(dep_fit_spec$original[[cond_var]])
  xmax <- max(dep_fit_spec$transformed[, cond_var])
  dif <- xmax - dth
  xlim <- c(dth - 0.1 * dif, dth + 1.5 * dif)

  # Determine upper limit of x-axis
  # if (marg$xi < 0 && xlim[2] > upper) {
  #   xlim[2] <- upper
  #   plim <- 1
  # } else {
  #   plim <- evd::pgpd(xlim[[2]], mth, marg$sigma, marg$xi)
  # }
  plim <- 1
  # CDF probabilities to plot at
  p <- seq(dqu, 1 - 1 / n, length = n)
  # take out largest point to avoid Inf in CDF transform
  len <- 501
  plotp <- seq(dqu, plim, len = len)[-len]
  # transform to original scale; these will be x-values in plot
  # TODO Gives different results to revTransform, why????
  # plotx <- inv_semi_par_cdf(
  #   F_hat = matrix(plotp),
  #   dat   = dep_fit_spec$original[cond_var],
  #   gpd   = dep_fit_spec$marginal[cond_var] # TODO Need conditional thresh, not marginal thresh!
  # )
  # plotx <- revTransform(
  #   plotp,
  #   data  = dep_fit_spec$original[[cond_var]],
  #   qu    = dqu,
  #   th    = dth,
  #   sigma = marg$sigma,
  #   xi    = marg$xi
  # )
  plotx <- as.vector(laplace_trans(plotp))

  # convert probs to Laplace scale (these are values to calculate CE line at)
  # xq <- dep$margins$p2q(plotp)
  xq <- laplace_trans(plotp)

  # plot on original margins
  # TODO Loop over conditioning variables
  # TODO Extend to > 2 variables
  # base_plot <- dep_fit_spec$original |>
  #   select(all_of(c(cond_var, vars))) |>
  labels <- lapply(c(cond_var, vars), stringr::str_replace_all, "_", " ")
  base_plot <- data.frame(dep_fit_spec$transformed) |>
    setNames(c("x", "y")) |>
    ggplot(aes(x, y)) +
    geom_point() +
    theme +
    # add vertical line at threshold
    geom_vline(xintercept = dth) +
    # labs(x = cond_var, y = vars)
    labs(x = labels[[1]], y = labels[[2]])


  # pull dependence coefficients and quantiles of residuals
  # co <- coef(dep)[, i]
  co <- dep_fit_spec$dependence[[vars]]
  # zq <- quantile(dep$Z[, i], quantiles)
  # zq <- quantile(dep_fit_spec$residual[[cond_var]], quantiles)
  zq <- quantile(dep_fit_spec$residual[[vars]], quantiles)

  # calculates regression lines from quantiles of residuals
  # yq <- sapply(zq, function(z, co, xq) {
  #   co["a"] * xq + co["c"] - co["d"] * log(xq) + xq^co["b"] * z
  # }, xq, co = co)
  yq <- sapply(zq, \(z, xq) {
    # dep$co$a * xq + dep$co$c - dep$co$d * log(xq) + xq^dep$co$b * z
    (co["a", ] * xq) + ((xq^co["b", ]) * z)
  }, xq = xq)

  # # transform to original scale
  # # ploty <- apply(dep$margins$q2p(yq), 2, revTransform,
  # #     data = trns,
  # #     qu = qu, th = th, sigma = sigma, xi = xi
  # #   )
  # # ploty <- inv_semi_par_cdf(
  # #   F_hat = yq,
  # #   dat   = dep_fit_spec$original[vars],
  # #   gpd   = dep_fit_spec$marginal[vars]
  # # )
  # # dth_spec <- quantile(dep_fit_spec$original[[vars]], dqu)
  # # ploty <- apply(yq, 2, \(y) {
  # ploty <- apply(inv_laplace_trans(yq), 2, \(y) {
  #   # print(y[1])
  #   # TODO again, not working, why???
  #   # inv_semi_par_cdf(
  #   #   F_hat = matrix(y),
  #   #   dat   = dep_fit_spec$original[vars],
  #   #   gpd   = dep_fit_spec$marginal[vars]
  #   # )
  #   revTransform(
  #     y,
  #     data  = dep_fit_spec$original[[vars]],
  #     qu    = dqu,
  #     th    = dth_spec,
  #     sigma = dep_fit_spec$marginal[[vars]]$sigma,
  #     xi    = dep_fit_spec$marginal[[vars]]$xi
  #   )
  # })
  ploty <- yq
  dth_spec <- dth

  # plot CE quantiles recursively
  add_line <- \(p, ploty) {
    if (length(ploty) == 0) {
      return(p)
    } else {
      add_line(
        p +
          geom_line(
            data     = data.frame(x = plotx, y = ploty[, 1]),
            mapping  = aes(x = plotx, y = ploty[, 1]),
            linetype = 2,
            col      = "blue"
          ),
        ploty[, -1, drop = FALSE]
      )
    }
  }
  p <- add_line(base_plot, ploty)

  return(p)
}

# produce quantile plots in all locations
quant_plot <- \(dep_fit) {
  lapply(names(dep_fit$residual), \(x) {
    tryCatch(
      {
        list(
          "rain"       = plot_ce_quantiles(dep_fit, x, "rain"),
          "wind_speed" = plot_ce_quantiles(dep_fit, x, "wind_speed")
        )
      },
      error = function(e) {
        message("Error in plotting quantiles for ", x, ": ", e$message)
        return(NA)
      }
    )
  })
}

# join resid and quantile plots
diag_plot <- \(resid_plots, quantile_plots, names_df = NULL) {
  lapply(seq_along(resid_plots), \(i) {
    tryCatch(
      {
        p <- wrap_plots(resid_plots[[i]]) /
          (quantile_plots[[i]]$rain + quantile_plots[[i]]$wind_speed)
        if (!is.null(names_df)) {
          p <- p +
            patchwork::plot_annotation(
              title = paste(names_df$name[i], names_df$county[i], sep = " - "),
              theme = theme
            )
        }
        return(p)
      },
      error = function(e) {
        message("Error in joining plots for ", names(resid_plots)[i], ": ", e$message)
        return(NA)
      }
    )
  })
}

# function for plotting alpha, beta estimates
map_plot <- \(ab_df, data, n_breaks = 8, range = c(1, 6), elev_df = NULL) {
  ab_sf <- ab_df |>
    left_join(distinct(data, name, lon, lat)) |>
    st_to_sf()

  # check that b values aren't exceptionally negative (<-1)
  # also check that all values are equal to fixed value for b, if applying
  # sort(ab_sf[ab_sf$parameter == "b", ]$value)[1:10]

  # plot parameter values for each parameter and variable
  names <- c("a", "b")
  p_lst <- lapply(seq_along(names), \(i) {
    # browser()
    p <- ggplot()
    # add elevation raster
    if (!is.null(elev_df)) {
      # plot bins if available, if not plot directly (less distinct colours)
      if ("elev_bin" %in% names(elev_df)) {
        p <- p +
          geom_tile(
            data = elev_df,
            aes(x = x, y = y, fill = elev_bin),
            width = diff(range(elev_df$x)) / length(unique(elev_df$x)),
            height = diff(range(elev_df$y)) / length(unique(elev_df$y)),
            show.legend = FALSE
          ) +
          labs(x = "", y = "") +
          # scale_fill_viridis_d(
          #   name = "Elevation",
          #   option = "D",
          #   direction = 1
          # )
          scico::scale_fill_scico_d(
            palette = "bilbao",
            name = "Elevation\n(m)",
            direction = -1
          )
      } else {
        p <- p +
          geom_tile(
            data = elev_df,
            aes(x = x, y = y, fill = elevation),
            show.legend = FALSE
          ) +
          labs(x = "", y = "") +
          scale_fill_viridis(
            name = "Elevation\n(m)",
            option = "D",
            limits = c(0, 1500)
          )
      }
      p <- p + geom_sf(data = areas, fill = NA, colour = "black")
    } else {
      p <- p + geom_sf(data = areas, fill = NA, colour = "black")
    }
    p +
      # geom_sf(data = areas, fill = NA, colour = "white") +
      coord_sf(expand = FALSE) + # remove padding around plot
      ggnewscale::new_scale_fill() +
      geom_sf(
        data = ab_sf %>%
          filter(parameter == names[i]) %>%
          mutate(vars = paste0(parameter, " - ", vars)),
        aes(fill = value, size = value),
        colour = "black",
        stroke = 1,
        pch = 21
      ) +
      scale_size_continuous(
        breaks = scales::extended_breaks(n = n_breaks),
        range = range,
        guide = "legend"
      ) +
      scale_fill_gradient2(
        low = "blue3",
        high = "red3",
        na.value = "grey",
        breaks = scales::extended_breaks(n = n_breaks),
        guide = "legend"
      ) +
      guides(fill = guide_legend(), size = guide_legend()) +
      labs(fill = "", size = "") +
      facet_wrap(
        ~vars,
        ncol = 2,
        labeller = as_labeller(c(
          # "a - rain"       = "alpha ~ ' - ' ~ 'Rain | Wind Speed'",
          # "a - wind_speed" = "alpha ~ ' - ' ~ 'Wind Speed | Rain'",
          # "a - wind_speed" = "alpha ~ ' - ' ~ 'Wind Speed | Rain'",
          # "b - rain"       = "beta ~ ' - ' ~ 'Rain | Wind Speed'",
          # "b - wind_speed" = "beta ~ ' - ' ~ 'Wind Speed | Rain'"
          "a - rain"       = "alpha ~ ' - ' ~ 'Precipitation | Wind Speed'",
          "a - wind_speed" = "alpha ~ ' - ' ~ 'Wind Speed | Precipitation'",
          "a - wind_speed" = "alpha ~ ' - ' ~ 'Wind Speed | Precipitation'",
          "b - rain"       = "beta ~ ' - ' ~ 'Precipitation | Wind Speed'",
          "b - wind_speed" = "beta ~ ' - ' ~ 'Wind Speed | Precipitation'"
        ), default = label_parsed)
      ) +
      theme
  })
  return(p_lst)
}


#### Bootstrapping functions ####

# function to perform bootstrapping for CE model for ECDF
boot_ce_ecdf <- \(
  dep,
  orig,
  transformed,
  # marg_pars = c(loc = 0, scale = 1, shape = -0.05),
  cond_prob = 0.9,
  # marg_prob = 0.98, # marginal quantile to calculate expectation after boots
  R = 100,
  trace = 10,
  ncores = 1,
  fixed_b = FALSE,
  ...
) {
  # pull data
  # TODO Should pull fixed_b from this!
  arg_vals <- list("cond_prob" = cond_prob) # TODO Add any others here
  dependence <- dep # dependence parameters

  # calculate marginal threshold
  # TODO Will have to change for real data, as can't use qgpd
  # marg_val <- do.call(evd::qgpd, c(list(marg_prob), marg_pars))

  # extract marginal and dependence thresholds
  # TODO Implement for multiple locations
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

  # Parallel setup
  apply_fun <- ifelse(ncores == 1, lapply, parallel::mclapply)
  ext_args <- NULL
  if (ncores > 1) {
    ext_args <- list(mc.cores = ncores)
  }
  loop_fun <- \(...) {
    do.call(apply_fun, c(list(...), ext_args))
  }

  # Pull start values
  start <- lapply(dependence, \(x) {
    lapply(x, \(y) {
      # scale back towards zero in case point est on edge of original parameter
      # space and falls off edge of constrained space for bootstrap sample
      if (fixed_b) {
        y[c("a", "b"), , drop = FALSE] * c(0.75, 1)
      } else {
        y[c("a", "b"), , drop = FALSE] * 0.75
      }
    })
  })

  # Function to prepare bootstrapped data for each location and conditioned var
  prep_boot_loc_dat <- \(
    dep_loc, trans_loc, orig_loc, thresh_dep, n_pass = 3
  ) {
    # pull bootstrap sample
    # (must be different for each loc as nrows may differ)
    indices <- sample(seq_len(nrow(trans_loc)), replace = TRUE)
    trans_loc_boot <- trans_loc[indices, ]

    # Reorder bootstrap sample to have the same order as original data
    test <- FALSE
    # which <- which(names(marg_loc) %in% cond_var)
    while (test == FALSE) {
      for (j in seq_along(dep_loc)) {
        # replace ordered Yi with ordered sample from standard Laplace CDF
        u <- matrix(runif(nrow(trans_loc_boot)))
        trans_loc_boot[
          order(trans_loc_boot[, j]), j
        ] <- sort(laplace_trans(u)[, 1])
      }
      # need exceedance in cond. var, and also in other vars for these rows
      if (any(colSums(sweep(trans_loc_boot, 2, thresh_dep, FUN = ">")) > 0)) {
        test <- TRUE
      }
    }

    # convert variables to original scale
    # orig_loc_boot <- inv_semi_par_cdf(
    #   # inverse Laplace transform to CDF
    #   inv_laplace_trans(trans_loc_boot),
    #   # original data for where semiparametric thresholding occurs
    #   orig_loc,
    #   # TODO Have to change if parameters are different
    #   lapply(seq_len(ncol(orig_loc)), \(i) {
    #     setNames(marg_pars[c("scale", "shape", "loc")], c("sigma", "xi", "thresh"))
    #   })
    # )

    # convert from Laplace scale back to ECDF (no need to get original data)
    ecdf_loc_boot <- inv_laplace_trans(trans_loc_boot)

    # test for no marg exceedances over sampled points, if so resample w/ nPass
    # max_vals <- apply(orig_loc_boot, 2, max, na.rm = TRUE)
    # marg_thresh <- marg_val # TODO Will have to change if vars don't have same margins
    # if (!all(max_vals > marg_thresh)) {
    #   return(list(NA))
    # }
    return(ecdf_loc_boot)
  }

  # perform bootstrapping
  boot_fits <- loop_fun(seq_len(R), \(i) {
    if (i %% trace == 0) {
      system(sprintf("echo %s", paste(i, "replicates done")))
    }
    dat_ecdf_boot <- lapply(seq_along(dependence), \(j) { # loop through locations
      # prepare bootstrapped data for location j
      # TODO Change argument order
      dat_spec <- prep_boot_loc_dat(
        dependence[[j]],
        transformed[[j]],
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
      return(dat_spec)
    })

    # transform from GPD to Laplace margins
    # dat_boot_trans <- lapply(dat_boot, \(x) {
    #   laplace_trans(do.call(evd::pgpd, c(list(x), marg_pars)))
    # })
    dat_boot_trans <- lapply(dat_ecdf_boot, laplace_trans)

    # refit CE model using same dependence quantile
    # TODO Could change to mapply?
    fit_boot <- lapply(seq_along(dat_boot_trans), \(j) {
      # if (i == 1 && j == 2) {
      # browser()
      # }
      o <- ce_optim(
        Y = dat_boot_trans[[j]],
        dqu = cond_prob,
        start = start[[j]],
        control = list(maxit = 1e6),
        fixed_b = fixed_b,
        ...
      )
    })
    # fit_boot <- mapply(\(dat, start_val) {
    #   o <- ce_optim(
    #     # Y = dat_boot_trans[[j]],
    #     Y = dat,
    #     dqu = cond_prob,
    #     # start = start[[j]],
    #     start = start_val,
    #     control = list(maxit = 1e6),
    #     ...
    #   )
    # }, dat_boot_trans, start)

    # extract just parameters
    fit_boot_pars <- lapply(fit_boot, \(x) lapply(x, `[[`, "params"))
    names(fit_boot_pars) <- names(dep)

    # Calculate conditional expectation at marg_val
    # return(lapply(fit_boot_pars, \(x) {
    #   lapply(x, \(y) {
    #     y["a", ] * marg_val + (marg_val^(y["b", ])) * y["m", ]
    #   })
    # }))
    # dep_out <- lapply(boot_fits, `[[`, "dependence") |>
    dep_out <- purrr::map_dfr(names(fit_boot_pars), function(loc_name) {
      # loc <- boot_sample[[1]]
      loc <- fit_boot_pars[[loc_name]]

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

    # Combine all the lists of data frames into one final data frame
    dep_out <- bind_rows(dep_out)
  })

  return(boot_fits)
}

# function to plot bootstrapped parameter estimates for different quantiles
plot_boot_quant <- \(
  data_lst,
  quantiles = seq(0.5, 0.9, by = 0.1),
  constrain = TRUE,
  R = 10,
  ncores = 1,
  county_key_df = NULL,
  start = c(0.01, 0.01),
  ...
) {
  # transform data
  Y_lst <- lapply(data_lst, trans_fun)

  # fit model at each threshold
  ce_fit_quant <- mclapply(quantiles, \(q) {
    ret <- lapply(Y_lst, \(y) {
      o <- ce_optim(
        Y = y,
        dqu = q,
        constrain = constrain,
        start = start,
        nruns = 3,
        ...
      )
    })
    locs <- names(ret)
    ret <- lapply(c("resid", "params"), \(x) {
      setNames(pull_element(ret, x), locs)
    })
    names(ret) <- c("residual", "dependence")
    return(ret)
  }, mc.cores = ncores)

  # extract parameter estimates and label by quantile (for plotting)
  ab_df <- bind_rows(lapply(seq_along(ce_fit_quant), \(i) {
    dep_to_df(ce_fit_quant[[i]]$dependence) |>
      mutate(quantile = quantiles[i])
  }))

  if (!is.null(county_key_df)) {
    ab_df <- ab_df |>
      left_join(county_key_df, by = c("name"))
  }

  # now for each model fit, bootstrap
  # TODO Again, need to make these parameters
  boot_fit_quant <- mclapply(seq_along(ce_fit_quant), \(i) {
    boot_ce_ecdf(
      dep = ce_fit_quant[[i]]$dependence,
      orig = data_lst,
      transformed = Y_lst,
      cond_prob = quantiles[i],
      # marg_prob = 0.98, # marginal quantile to calculate expectation after boots
      # R = 500,
      R = R,
      trace = R + 1, # don't show
      constrain = constrain,
      ...
    )
  }, mc.cores = ncores)

  # join all together
  boot_ab_df <- bind_rows(lapply(seq_along(boot_fit_quant), \(i) {
    ret <- bind_rows(boot_fit_quant[[i]], .id = "run") |>
      mutate(quantile = quantiles[i])
  }))

  loc_df <- ab_df |>
    distinct(across(matches("name") | matches("county")))
  plots <- lapply(seq_along(locs), \(i) {
    title <- ifelse(
      is.null(county_key_df),
      loc_df$name[[i]],
      paste0(loc_df$name[i], "-", loc_df$county[i])
    )

    ab_df |>
      filter(name == loc_df$name[[i]]) |>
      ggplot(aes(x = quantile, y = value)) +
      geom_point(
        data = boot_ab_df |>
          filter(name == loc_df$name[[i]]) |>
          mutate(vars = ifelse(vars == "windspeed", "wind_speed", vars)),
        colour = "black",
        size = 2,
        alpha = 0.7
      ) +
      geom_line() +
      geom_point(size = 5, colour = "orange", alpha = 0.7) +
      facet_wrap(~ parameter + vars, scales = "free") +
      theme +
      scale_x_continuous(breaks = quantiles) +
      labs(
        x     = "Quantile",
        y     = "Parameter estimate",
        title = title
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  names(plots) <- loc_df$name

  return(plots)
}


#### Metadata ####

# marginal target threshold (wind speed requires higher threshold)
# mqu <- c("rain" = 0.9, "wind_speed" = 0.98)
# dqu <- 0.9 # dependence threshold
# dqu2 <- c(0.9, 0.95)
# dqu2 <- c(0.9, 0.95) # increased threshold for wind speed
min_max_dates <- as_date(c("1990-01-01", "2020-12-31"))
all_dates <- seq.Date(min_max_dates[1], min_max_dates[2], by = "day")
fixed_xi <- TRUE # whether to fix shape parameter in GPD margins
seed_number <- 123
max_k <- 10

ncores <- detectCores() - 1


#### Load Data ####

data <- readr::read_csv(
  "data/met_eireann/final/met_eir_wind.csv.gz",
  # "data/met_eireann/final/met_eir_wind_alt.csv.gz",
  show_col_types = FALSE
) %>%
  filter(date %in% all_dates, !is.na(wind_speed))

data <- data |>
  # Remove sites with negative alpha values
  # filter(!name %in% c(
  #   "Belfast Newforge", "Belmont", "Carmoney", "Castlereagh",
  #   "Dungonnell Filters No 2", "Marble Arch Caves", "Pomeroy Primary School",
  #   "Woodburn North"
  # ))
  # remove six counties :(
  filter(!county %in% c(
    "Down", "Antrim", "Fermanagh", "Tyrone", "Armagh", "Derry"
  ))

# load shapefile
areas <- sf::read_sf("data/met_eireann/final/irl_shapefile.geojson")

# pull elevation for areas
# TODO May have to add another bin for 600m, since over that is a mountain!
# elev_df <- get_elev_rast(
#   areas,
#   z = 9,
#   bins = c(seq(0, 600, by = 100), Inf),
#   labels = c(
#     "0–100 m", "100–200 m",
#     "200–300 m", "300–400 m", "400–500 m", "500-600 m", "> 600 m"
#   )
# )
# readr::write_csv(
#   elev_df,
#   "data/met_eireann/final/irl_elev.csv"
# )
elev_df <- readr::read_csv(
  "data/met_eireann/final/irl_elev.csv",
  show_col_types = FALSE
)

# TODO Change height bins?? > 600m is a mountain, so would make sense!
if (!is.factor(elev_df$elev_bin)) {
  elev_df <- elev_df |>
    mutate(
      elev_bin = stringr::str_remove(elev_bin, " m"),
      elev_bin = factor(
        elev_bin,
        # levels = c(
        #   "0–100 m", "100–200 m",
        #   "200–300 m", "300–400 m", "400–500 m", "500-600 m", "> 600 m"
        # )
        levels = c(
          "0–100", "100–200",
          "200–300", "300–400", "400–500", "500-600", "> 600"
        )
      )
    )
}

# pull just site names, counties and provinces
county_key_df <- data |>
  distinct(name, county, province)
readr::write_csv(
  county_key_df,
  "data/met_eireann/final/ire_county_key.csv"
)

# extract point location of each station for plotting on map
pts <- data %>%
  distinct(name, lon, lat) %>%
  st_to_sf()
sf::write_sf(pts, "data/met_eireann/final/irl_points.geojson")


# calculate adjacency matrix from Voronoi cells,
# to restrict clustering to adjacent sites
# TODO Check if correct
# adj_mat <- calc_adj_mat2(pts)


#### Preprocess Data ####

# take Winter data only
data_winter <- data %>%
  mutate(month = as.numeric(substr(date, 6, 7))) %>%
  filter(month %in% c(1:3, 10:12)) %>%
  dplyr::select(-month)

# Take weekly sum of precipitation and average of wind speed, to account for
# lag in storm (as in Vignotto, Engelke 2021 study of GB + Ireland)
data_week <- data_winter %>%
  mutate(week = floor_date(date, "week")) %>%
  group_by(week, name, county, province, lon, lat) %>%
  summarise(
    rain       = sum(rain, na.rm = TRUE),
    wind_speed = mean(wind_speed, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  rename(date = week) |>
  # remove weeks with no rainfall, as in Vignotto 2021 (required/important?)
  filter(rain != 0)

# for many problem locs, there doesn't seem to be much relationship between
# whether extremes happen together or separately?
# data_week |>
#   filter(
#     # name == "Belfast Newforge"
#     # name == "Belmont"
#     # name == "Dungonnell Filters No 2"
#     # name == "Lough Navar Forest"
#     # name == "Marble Arch Caves"
#     name == "Woodburn North"
#   ) |>
#   mutate(across(
#     c(rain, wind_speed), ~ quantile(.x, 0.90, na.rm = TRUE),
#     .names = "quant_{.col}"
#   )) %>%
#   ggplot(aes(x = rain, y = wind_speed)) +
#   geom_point(aes(), size = 1.5, alpha = 0.8) +
#   geom_vline(aes(xintercept = quant_rain), linetype = "dashed") +
#   geom_hline(aes(yintercept = quant_wind_speed), linetype = "dashed") +
#   evc::evc_theme()


#### Exploration of sites ####

# Identify sites with highest and lowest rainfall
highest_lowest_rain <- data_week %>%
  group_by(name) %>%
  summarise(rain = mean(rain, na.rm = TRUE), .groups = "drop") %>%
  arrange(rain) %>%
  slice(c(1, n())) %>%
  pull(name)

# highest_lowest_rain <- c("Malahide Castle", "Dublin (Ringsend)")
highest_lowest_rain <- c("Dublin (Ringsend)", "Kilcar (Cronasillagh)")


# plot sites on left and locations with highest and lowest rainfall on right
data_plot <- data_week %>%
  # mutate(indicator = ifelse(name %in% highest_lowest_rain, name, NA)) %>%
  mutate(indicator = ifelse(name %in% highest_lowest_rain, name, "other")) %>%
  arrange(desc(indicator)) %>%
  st_to_sf()

# elevation plot (can superimpose others on top)
# TODO What to do with longitude/latitude?
p_terrain <- ggplot() +
  geom_tile(
    data = elev_df,
    aes(x = x, y = y, fill = elev_bin),
    width = diff(range(elev_df$x)) / length(unique(elev_df$x)),
    height = diff(range(elev_df$y)) / length(unique(elev_df$y))
  ) +
  scico::scale_fill_scico_d(
    palette = "bilbao",
    name = "Elevation (m)",
    direction = -1
  ) +
  # scale_fill_manual(
  #   values = c(
  #     "#FFFFFF", "#C5C2B2", "#B19E68",
  #     "#A6785B", "#9B5352", "#6D1F23"
  #   )
  # ) +
  labs(x = "", y = "", fill = "Elevation (m)") +
  theme +
  # TODO position legend in bottom right of plot?
  # theme(
  #   legend.direction = "vertical",
  #   legend.position = c(0.75, 0.1),
  #   legend.title = element_text(size = 12, face = "bold"),
  #   legend.text = element_text(size = 12)
  # ) +
  ggnewscale::new_scale_fill() +
  geom_sf(data = areas, colour = "black", fill = NA) +
  # remove padding around plot
  coord_sf(expand = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  NULL

# ggsave(filename = "test.png", p_terrain)

# plot the location of each site
p1 <- p_terrain +
  # points other than two sites
  geom_sf(
    data = filter(data_plot, indicator == "other"),
    colour = "black",
    # size = 3,
    size = 4,
    # change shape to hollow circle to better see elevation
    shape = 21,
    stroke = 0.9,
    alpha = 0.9,
    show.legend = FALSE
  ) +
  # points for two sites with highest and lowest rain
  geom_sf(
    data = filter(data_plot, indicator != "other"),
    aes(colour = indicator, size = indicator),
    # have thicker border to make it more visible
    # shape = 21,
    # stroke = 2,
    alpha = 0.9,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = c(ggsci::pal_nejm()(2))) +
  # scale_size_manual(values = c(4.5, 4.5)) +
  scale_size_manual(values = c(5, 5)) +
  labs(colour = "", size = "") +
  coord_sf(expand = FALSE) + # remove padding around plot
  # theme +
  evc_theme(legend.position = NULL, nejm_pal = FALSE) +
  # remove axis text
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank() # ,
    # legend.key = element_blank()
  )

# ggsave(filename = "test_map_plot.png", plot = p1, width = 10, height = 8)

cols <- ggsci::pal_nejm()(4)
cols[2] <- "black"
p21 <- data_plot %>%
  # filter(name == highest_lowest_rain[2], rain > 0) %>%
  filter(name %in% highest_lowest_rain, rain > 0) |>
  group_by(name) %>%
  mutate(across(c(rain, wind_speed), ~ quantile(.x, 0.95, na.rm = TRUE), .names = "quant_{.col}")) %>%
  ungroup() %>%
  mutate(col = case_when(
    rain > quant_rain & wind_speed > quant_wind_speed ~ "Both",
    rain > quant_rain & wind_speed <= quant_wind_speed ~ "Rain",
    rain <= quant_rain & wind_speed > quant_wind_speed ~ "Wind",
    TRUE ~ "Neither"
  )) %>%
  ggplot(aes(x = rain, y = wind_speed)) +
  # geom_point(aes(colour = name), size = 1.5, alpha = 0.9) +
  geom_point(aes(colour = col), size = 1.5, alpha = 0.9) +
  # TODO Investigate why these aren't showing up???
  geom_vline(aes(xintercept = quant_rain), linetype = "dashed") +
  geom_hline(aes(yintercept = quant_wind_speed), linetype = "dashed") +
  # facet_wrap(~ name, scales = "free_x") +
  scale_colour_manual(values = cols) +
  labs(
    # x = "Weekly total precipitation (mm)",
    x = "precipitation (mm)",
    y = "wind speed (m/s)",
    colour = ""
  ) +
  guides(colour = "none") +
  facet_wrap(~name, scales = "fixed") +
  # theme +
  # remove facet labels, colour will do
  # theme(
  #   strip.background = element_blank(),
  #   strip.text.x = element_blank(),
  #   legend.key = element_blank()
  # ) +
  NULL

# Second, plot wind speeds against rain for sites with the highest and lowest rainfall
# TODO: Add 95% quantile lines for both (?)
p2 <- data_plot %>%
  filter(name %in% highest_lowest_rain, rain > 0) %>%
  group_by(name) %>%
  # mutate(across(c(rain, wind_speed), ~ quantile(.x, 0.95, na.rm = TRUE), .names = "quant_{.col}")) %>%
  mutate(
    quant_rain = quantile(rain, 0.95, na.rm = TRUE),
    quant_wind_speed = quantile(wind_speed, 0.95, na.rm = TRUE)
  ) |>
  ungroup() %>%
  # ggplot(aes(x = rain, y = wind_speed)) +
  ggplot(aes(x = wind_speed, y = rain)) +
  # geom_point(aes(colour = name), size = 1.5, alpha = 0.9) +
  geom_point(aes(colour = name), size = 1.5, alpha = 0.9) +
  geom_vline(aes(xintercept = quant_wind_speed), linetype = "dashed") +
  geom_hline(aes(yintercept = quant_rain), linetype = "dashed") +
  facet_wrap(~name, scales = "fixed") +
  scale_colour_manual(values = c(ggsci::pal_nejm()(2))) +
  scale_x_continuous(
    # limits = c(0, 12),
    limits = c(2, 12),
    # breaks = seq(2.5, 12.5, by = 2.5),
    # breaks = seq(0, 12, by = 2),
    breaks = seq(2, 12, by = 2),
    expand = c(0, 0.2)
  ) +
  scale_y_continuous(
    limits = c(0, 310),
    expand = c(0.02, 0)
  ) +
  labs(
    # x = "Weekly total precipitation (mm)",
    # x = "precipitation (mm)",
    y = "precipitation (mm)",
    # y = "wind speed (m/s)", # TODO: What is the unit of ws?
    x = "wind speed (m/s)", # TODO: What is the unit of ws?
    colour = ""
  ) +
  theme +
  # remove facet labels, colour will do
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.key = element_blank()
  )
# p2

# join plots
# TODO: Change size of first plot to be larger!
p_sec_2 <- p1 +
  (p2 + guides(colour = "none", size = "none")) +
  # have common legends
  # plot_layout(guides = "collect") &
  # # theme(legend.position = "bottom") +
  # theme(legend.position = "left") +
  NULL

ggsave("plot_1_test.png", p_sec_2, width = 12, height = 6, units = "in")
# ggsave("latex/plots/02_mot_ex_plot.png", p_sec_2, width = 10, height = 6, units = "in")
# ggsave("latex/plots/02_mot_ex_plot_elev.png", p_sec_2, width = 12, height = 6, units = "in")


#### Chi exploration ####

# Plot chi and chi-squared for one location
loc <- "Costelloe Fishery"
chi_spec <- data_week %>%
  # filter(name == highest_lowest_rain[2]) %>%
  filter(name == loc) |>
  dplyr::select(rain, wind_speed) %>%
  texmex::chi()

(chi_plot_spec <- chi_spec %>%
  # use texmex plotting method for chi and chibar
  ggplot(main = c("ChiBar" = "", "Chi" = ""), plot. = FALSE) |>
  lapply(\(x) x + theme) %>%
  # wrap_plots() +
  `[[`(2) +
  # add centred title through patchwork
  patchwork::plot_annotation(
    title = paste0("Tail Dependence, ", loc),
    theme = theme
  ))
# saveRDS(chi_plot_spec, file = "temp.RDS")
ggsave("latex/plots/chi_plot_spec.png", chi_plot_spec, width = 6, height = 6, units = "in")

# for each location, pull out 95th quantile for chibar
# Colour chi grey if not valid
# Plot on map
# TODO Functionalise for US air pollution!
names <- unique(data_week$name)
chi_95_df <- bind_rows(lapply(names, \(x) {
  chi <- data_week %>%
    filter(name == x) %>%
    dplyr::select(rain, wind_speed) %>%
    texmex::chi()

  # whether to show chi or not, based on whether chibar upper extend crosses 1
  show_chi <- !prod(tail(chi$chibar[, 3]) < 1)

  loc <- which.min(abs(chi$quantile - 0.95))
  return(data.frame(
    "name" = x,
    "chi" = chi$chi[loc, 2, drop = TRUE],
    "chibar" = chi$chibar[loc, 2, drop = TRUE],
    "show_chi" = show_chi
  ))
  # }, mc.cores = ncores))
}))
rownames(chi_95_df) <- NULL

# join in area statistics
chi_95_sf <- chi_95_df %>%
  pivot_longer(c("chi", "chibar"), names_to = "var") %>%
  # always show chibar plot
  mutate(show_chi = ifelse(value == "chibar", TRUE, show_chi)) %>%
  left_join(
    distinct(data, name, lon, lat)
  ) %>%
  st_to_sf()

# function to plot chi and chi bar on map
scales <- seq(-0.1, 0.6, by = 0.1)
chi_map_plot <- \(
  chi_95_sf,
  var = c("chi", "chibar"),
  scales = seq(-0.1, 0.6, by = 0.1),
  point_ranges = c(2, 6),
  rm_axis = TRUE
) {
  # lab <- ifelse(var == "chi", expression(chi(u)), expression(bar(chi)(u)))
  lab <- ifelse(var == "chi", expression(chi(0.95)), expression(bar(chi)(0.95)))
  plot_data <- filter(chi_95_sf, var == !!var)
  if (var == "chi") {
    plot_data <- plot_data %>%
      mutate(show_chi = factor(ifelse(show_chi == TRUE, "yes", "no")))
  }

  p <- ggplot(areas) +
    geom_sf(colour = "black", fill = NA) +
    geom_sf(
      # geom_sf_pattern(
      data = plot_data,
      aes(fill = value, size = value),
      # aes(fill = value, size = value, pattern = show_chi),
      pch = 21,
      stroke = 1
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    # add colour scheme afterwards
    # scale_fill_gradientn(
    #   colours = RColorBrewer::brewer.pal(name = "Blues", n = 7),
    #   breaks = scales_chi,
    #   labels = as.character(scales_chi),
    #   guide = "legend"
    # ) +
    scale_size_continuous(
      range  = point_ranges,
      breaks = scales,
      labels = as.character(scales),
      guide  = "legend"
    ) +
    # scale_pattern_manual(values = c("stripe", "none")) +
    # maximise plot within frame
    coord_sf(expand = FALSE) +
    labs(fill = lab, size = lab) +
    guides(fill = guide_legend(), size = guide_legend(), pattern = "none") +
    theme +
    theme(legend.position = "right", legend.key = element_blank())

  # remove axis text and ticks if required
  if (rm_axis == TRUE) {
    p <- p +
      theme(
        axis.text  = element_blank(),
        axis.ticks = element_blank()
      )
  }

  return(p)
}

# plot chibar and chi
chibar_p <- chi_map_plot(chi_95_sf, "chibar", rm_axis = FALSE) +
  scale_fill_gradientn(
    colours = rev(heat.colors(7)),
    breaks = scales,
    labels = as.character(scales),
    guide = "legend"
  )
# chi_p <- chi_map_plot(chi_95_sf, "chi") +
chi_p <- chi_map_plot(chi_95_sf, "chi", rm_axis = FALSE) +
  scale_fill_gradientn(
    colours = RColorBrewer::brewer.pal(name = "Blues", n = 7),
    breaks = scales,
    labels = as.character(scales),
    guide = "legend"
  )

# combine plots
(chi_plots <- chibar_p + chi_p)

# save for just chi (presentation) and both
ggsave("latex/plots/041_chi_plots.png", chi_plots, width = 6.3, height = 6, units = "in")
ggsave("latex/plots/chi_map_ire.png", chi_p, width = 6, height = 6, units = "in")

# Also combine chi plot with elevation plot
# p_terrain_chi <- p_terrain + (chi_p + theme(legend.position = "bottom"))
p_terrain_chi <- (p_terrain + theme(legend.position = "right")) + chi_p
# p_terrain_chi

ggsave("latex/plots/043_chi_elevation.png", p_terrain_chi, width = 8, height = 6)


#### ECDF -> Laplace Transformation ####

# split data by location, convert to list of matrices
data_lst <- data_week |>
  group_split(name, .keep = TRUE) # keep name to label list afterwards
locs <- purrr::map_chr(data_lst, ~ as.character(.x$name[1]))
names(data_lst) <- locs
data_lst <- lapply(data_lst, \(x) {
  matrix(
    c(x$rain, x$wind_speed),
    ncol = 2,
    dimnames = list(NULL, c("rain", "wind_speed"))
  )
})

# convert from ECDF to Laplace margins
Y_lst <- lapply(data_lst, trans_fun)


#### Select dependence quantile via bootstrapping ####

# plot params at various quantiles for each loc, w/ bootstrapped samples
quant_plots <- plot_boot_quant(
  data_lst,
  quantiles = sort(c(seq(0.5, 0.9, by = 0.1), 0.85, 0.88, 0.95, 0.98)),
  constrain = TRUE,
  R = 30,
  ncores = ncores,
  # ncores = 1,
  county_key_df = county_key_df
)
saveRDS(quant_plots, "plots/tests/boot_quantiles.RDS")

pdf("plots/tests/boot_quantiles.pdf", width = 14, height = 8)
quant_plots
dev.off()

# notes:
# - Seem to get worse after around 90th quantile (or at least for some locations)
# - (somewhat) frequent pattern whereby params -> 0 as quantile -> 1
#   Suppose that's from an increasing lack of data?
# - Examples where beta is very low seem to be from problems with start values
#   or convergence issues with optimisation
# - Examples where a & b estimates aren't inside bootstrap sample, particularly
#   for higher quantile levels (>= 90) (Woodburn North Antrim, Lough Navar Forest Fermanagh)

# Individual problem locs:
# rain:
# - Belfast Newforge (Antrim) - a estimates outside more realistic bootstrapped ones,
#  (and negative), strangely except for 0.9 quantile
#  Far more obs extreme for just rain than for rain + wind speeds
# - Belmont Antrim - a estimates very strange < 0.95, maybe just remove
# - Costelloe Fishery: wide bootstrap CI @ 90%
#   Same with Derryhillagh (outside boots), Glenamaddy Galway,
# - Dungonnell Filters No 2 Antrim: a estimates way below bootstrapped ones
# - Same with Pomeroy Primary School (Tyrone)
# - Marbe Arch Caves Fermanagh has strange a (also windspeed), may just remove
# - Lough Navar Fermanagh also weird though! Bootstrapped estimates often
#   negative
# wind speed:
# - Carmoney Derry a values start positive and end negative
#   Same with Cashel, Castleisland
# - Dungonell Filters No 2 also also very large bootstrap CIs


#### Fit CE ####

# decision: 0.85 for both? Seems to avoid instabilities observed for large q
dqu <- 0.85
# dqu <- c("rain" = 0.88, "wind_speed" = 0.85)
# dqu <- 0.9

# fit CE model
ce_fit <- lapply(Y_lst, \(x) {
  o <- ce_optim(
    Y = x,
    # dqu = 0.7,
    dqu = dqu,
    # dqu = 0.95,
    cond_var = c("rain", "wind_speed"),
    control = list(maxit = 1e6),
    # constrain = FALSE,
    constrain = TRUE,
    start = c(0.01, 0.01),
    fixed_b = FALSE, # TODO limit a, b to be positive
    nruns = 3
  )
})
saveRDS(ce_fit, "data/ce_fit.RDS")
ce_fit <- readRDS(file = "data/ce_fit.RDS")

# pull dependence parameters and residuals out separately
dep_fit <- lapply(c("resid", "params"), \(x) {
  setNames(pull_element(ce_fit, x), locs)
})
names(dep_fit) <- c("residual", "dependence")

# add transformed data etc for plotting functions
dep_fit$transformed <- Y_lst
dep_fit$original <- data_week
dep_fit$arg_vals <- list("cond_prob" = dqu)
# saveRDS(dep_fit, file = "data/dep_fit.RDS")
dep_fit <- readRDS("data/dep_fit.RDS")

# pull dependence paramaters into df
ab_df <- dep_to_df(dep_fit$dependence)

# OLD OLD OLD (before removing NI sites)
# why so many locations with low alpha values?
ab_df |>
  filter(parameter == "a" & value < -0.2)
# problem locations from before have low alpha, might be a good idea to
# limit to positive values
# All locations here are from NI, may indicate problem with data? Remove!

ce_fit_final <- ce_fit
dep_fit_final <- dep_fit


#### CE Diagnostics ####

# produce residual plots for all locations
resid_plots <- resid_plot(dep_fit)

# produce quantile plots in all locations
quantile_plots <- quant_plot(dep_fit)
names(quantile_plots) <- names(resid_plots)

# join resid and quantile plots
diag_plots <- diag_plot(resid_plots, quantile_plots, arrange(county_key_df, name))
names(diag_plots) <- names(resid_plots)

pdf(paste0("plots/tests/residuals_dqu_", dqu, ".pdf"), width = 8, height = 8)
diag_plots
dev.off()

# plot results on map
map_plots <- map_plot(
  ab_df,
  data_week,
  n_breaks = 8,
  elev_df = elev_df
)

# ggsave(filename = "test.png", map_plots[[1]] / map_plots[[2]], width = 12, height = 12)
gc()
saveRDS(map_plots, "latex/plots/map_plots.RDS")
map_plot_save <- map_plots[[1]] / map_plots[[2]]
# saveRDS(map_plot_save, "test.RDS")
# map_plot_save

ggsave(
  filename = "latex/plots/ire_ce_new.png",
  map_plot_save,
  width = 12, height = 12
)

# choose one (nicer) diagnostic plot to include in presentation
p_spec <- diag_plots$`Costelloe Fishery` +
  plot_annotation(
    title = "Costelloe Fishery",
    theme = theme
  )
ggsave(p_spec,
  filename = "latex/plots/041_costelloe_fishery.png",
  width = 6, height = 6, units = "in"
)


# Comments:
# W/o NI counties, individual plots actually look really good!
# Is there a plot for every location???
# However, some very small beta values observed for 9/72 locations!
# Why might that be??
ab_df |>
  filter(value < -0.2) |>
  left_join(county_key_df, by = "name")


#### Limit alpha to be positive ####

# ce_fit_alpha <- lapply(Y_lst, \(x) {
#   o <- ce_optim(
#     Y = x,
#     dqu = dqu,
#     cond_var = c("rain", "wind_speed"),
#     control = list(maxit = 1e6),
#     # constrain = FALSE,
#     constrain = TRUE,
#     start = c(0.01, 0.01),
#     aLow = 0,
#     nruns = 3
#   )
# })
#
# dep_fit_alpha <- lapply(c("resid", "params"), \(x) {
#   setNames(pull_element(ce_fit_alpha, x), locs)
# })
# names(dep_fit_alpha) <- c("residual", "dependence")
#
# # add transformed data etc for plotting functions
# dep_fit_alpha$transformed <- Y_lst
# dep_fit_alpha$original <- data_week
# dep_fit_alpha$arg_vals <- list("cond_prob" = dqu)
#
# ab_df_alpha <- dep_to_df(dep_fit_alpha$dependence)
#
# # TODO Investigate differences
# # TODO Investigate failed plots
# resid_plots_alpha <- resid_plot(dep_fit_alpha)
# quantile_plots_alpha <- quant_plot(dep_fit_alpha)
# diag_plots_alpha <- diag_plot(
#   resid_plots_alpha,
#   quantile_plots_alpha,
#   arrange(county_key_df, name)
# )
# names(diag_plots_alpha) <- names(resid_plots_alpha)
#
#
# pdf(paste0(
#   "plots/tests/residuals_dqu_",
#   dqu,
#   "_alpha.pdf"
# ), width = 8, height = 8)
# # resid_plots_alpha
# diag_plots_alpha
# dev.off()
#
# map_plots_alpha <- map_plot(
#   ab_df_alpha,
#   data_week,
#   n_breaks = 8
# )
#
# # investigate locations which previously had very low alpha values
# prob_locs <- ab_df |>
#   filter(parameter == "a" & value < -0.2) |>
#   distinct(name) |>
#   pull()
# ab_df_alpha |>
#   filter(name %in% prob_locs) |>
#   left_join(
#     ab_df |>
#       filter(name %in% prob_locs) |>
#       rename(value_prev = value)
#   )
# # b values seem slightly higher, and a values shrink to basically 0
#
#
# # find when different to previous fit; this is where plots will differ
# ab_df_alpha |>
#   left_join(
#     rename(ab_df, "value_prev" = value)
#   ) |>
#   mutate(diff = abs(value - value_prev)) |>
#   filter(diff > 0.001) |>
#   arrange(desc(diff))
# # pretty much same as above, biggest changes are just from having alpha >= 0
#
# # conclusion: fix alpha > 0 as it will reduce impact of strange location on
# # clustering
#
# ce_fit_final <- ce_fit_alpha
# dep_fit_final <- dep_fit_alpha


#### Bootstrap ####

# boot_fit <- boot_ce_ecdf(
#   dep = dep_fit_final$dependence,
#   orig = data_lst,
#   transformed = Y_lst,
#   # cond_prob = 0.9,
#   cond_prob = dqu,
#   # marg_prob = 0.98, # marginal quantile to calculate expectation after boots
#   R = 500,
#   trace = 10,
#   ncores = ncores,
#   constrain = TRUE
# )
# saveRDS(boot_fit, "data/ecdf_boot_fit.RDS")
boot_fit <- readRDS("data/ecdf_boot_fit.RDS")

boot_df <- bind_rows(boot_fit, .id = "run") |>
  group_by(name, parameter, vars) |>
  mutate(mean = mean(value, na.rm = TRUE)) |>
  ungroup()

# plot for each location
boot_df |>
  # filter(name == "Armagh") |>
  filter(name == "Malahide Castle") |>
  ggplot(aes(x = vars, y = value, fill = parameter)) +
  geom_boxplot() +
  # facet_wrap(~ parameter + vars) +
  facet_wrap(~ parameter + vars, scales = "free") +
  theme

boot_df |>
  filter(name == "Malahide Castle") |>
  ggplot(aes(x = value, fill = parameter)) +
  geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.7) +
  geom_density(alpha = 0.5) +
  # geom_vline(
  #   aes(xintercept = mean),
  #   linetype = "dashed",
  #   color = "black"
  # ) +
  facet_wrap(~ parameter + vars, scales = "free") +
  # facet_grid(
  #   rows = vars(vars),
  #   cols = vars(parameter),
  #   scales = "free_x"
  # ) +
  theme +
  guides(fill = "none")

# also plot a vs b for each, add the original model estimate and the
# multivariate (geometric) median

# calculate multivariate median
boot_med_df <- boot_df |>
  select(-c(mean, cond_var)) |>
  pivot_wider(names_from = parameter, values_from = value) |>
  group_split(name, vars, .keep = TRUE) |>
  lapply(\(x) {
    ret <- Gmedian::Gmedian(x[, c("a", "b")])
    data.frame("a" = ret[1], "b" = ret[2], "vars" = x$vars[1], name = x$name[1])
  }) |>
  bind_rows()

# loc <- "Crolly (Filter Works)"
loc <- "Malahide Castle"
boot_df |>
  select(-c(mean, cond_var)) |>
  filter(name == loc) |>
  pivot_wider(names_from = parameter, values_from = value) |>
  ggplot() +
  geom_point(aes(x = a, y = b, colour = vars)) +
  # plot model estimate
  geom_point(
    data = ab_df |>
      filter(name == loc) |>
      pivot_wider(names_from = parameter, values_from = value),
    aes(x = a, y = b),
    # large X
    size = 5,
    stroke = 2,
    colour = "black",
    shape = 4
  ) +
  # add geometric median
  geom_point(
    data = boot_med_df |>
      filter(name == loc),
    aes(x = a, y = b),
    # large +
    size = 5,
    stroke = 2,
    colour = "black",
    shape = 3
  ) +
  facet_wrap(~vars, scales = "free") +
  theme


#### Clustering: choose k ####

# extract JS distance matrix
set.seed(123)
clust_obj <- js_clust(
  dep_fit_final$dependence,
  trans = dep_fit_final$transformed,
  scree_k = 1:max_k,
  n = n_mc,
  mc_method = mc_method,
  laplace_cap = laplace_cap
  # mc_method = "uniform"
)
dist_mat <- clust_obj$dist_mat

# look at 2-norm for each location
sort(apply(as.matrix(dist_mat), 2, mean))
sort(apply(as.matrix(dist_mat), 2, norm, type = "2"))
# interesting that Malahide and Ringsend, where rain was lowest, are
# notable outliers (especially Malahide)

# TODO function from Gaussian copula simulations, also needs to be packaged up!
# Function to find optimal k value via TWGSS and AIC
find_k <- \(
  ce_fit,
  data_mix, # list of data at different locations
  max_clust,
  cond_prob = 0.9,
  fixed_b = FALSE,
  adj_mat = NULL,
  spec_vars = NULL,
  ...
) {
  # Pull JS distance matrix
  dist <- js_clust(
    ce_fit,
    scree_k = 1:max_clust,
    spec_vars = spec_vars,
    n = n_mc
  )$dist_mat

  # apply adjacency if specified
  if (!is.null(adj_mat)) {
    dist_adj <- as.matrix(dist)
    dist_adj[adj_mat == 0] <- 1e9 # arbitrarily large number
    dist <- as.dist(dist_adj)
  }

  # extract TWGSS
  twgss <- evc::scree_plot(dist, k = 1:max_clust)

  # Calculate AIC and perform Likelihood Ratio test (i.e. model-based criteria)
  # initialise list of dependence models
  dep_clust_lst <- rep(
    get_pars(ce_fit),
    max_clust
  )
  aic <- lr <- vector(mode = "list", length = max_clust)
  for (k in 1:max_clust) {
    # cluster
    # debugonce(js_clust)
    pam_fit <- js_clust(
      dist_mat = dist,
      k = k,
      return_dist = TRUE,
      spec_vars = spec_vars,
      n = n_mc
    )
    # refit CE model
    dependence_clust <- fit_optim_clust(
      pam_fit$pam$clustering,
      # data_mix,
      data_lst,
      n_vars = 2,
      cond_prob = cond_prob,
      trans_fun = trans_fun,
      start_vals = c(0.01, 0.01),
      fixed_b = fixed_b,
      ...
    )
    # extract parameters (remove residuals)
    dependence_clust <- get_pars(dependence_clust)
    dep_clust_lst[[k]] <- dependence_clust

    # calculate AIC
    aic[[k]] <- ce_aic(dependence_clust)

    # calculate likelihood ratio statistic
    # if (k > 2) {
    # if (k > 1) {
    #   # df = 4 (n pars) * 2 (n_vars) (# more parameters in "full" model)
    #   df <- 4 * n_vars
    #   lr[[k]] <- lrt_test(
    #     dep_clust_lst[[k - 1]], dep_clust_lst[[k]],
    #     df = df
    #   )
    # }
  }

  aic <- unlist(aic)
  # lr <- c(NA, NA, unlist(lapply(lr, "[[", "p_value")))
  # lr <- c(NA, unlist(lapply(lr, "[[", "p_value")))

  # find best k as mode of choice for TWGSS, AIC and LRT
  get_mode <- function(x) {
    uniq_x <- unique(x)
    tab <- table(x)
    return(uniq_x[which.max(tab)])
  }
  # use elbow finding algorithm to choose k from TWGSS and AIC
  k_twgss <- find_elbow(twgss)
  k_aic <- find_elbow(aic)

  # find first LR test statistic w/ p-value < .05
  # k_lr <- which(lr < 0.05)[1]

  # if no LR test statistic is significant, choose from elbow plots
  elbow_k <- get_mode(c(k_twgss, k_aic))
  # also do so if LRT chooses lowest k but other methods choose lower k,
  # as our LRT can only be computed for k > 2 (cannot cluster with k == 1)
  # if (length(k_lr) == 0 || (k_lr == 2 && elbow_k < 2)) {
  #   k_lr <- elbow_k
  # }

  # ks <- c("k_twgss" = k_twgss, "k_aic" = k_aic, "k_lr" = k_lr)
  ks <- c("k_twgss" = k_twgss, "k_aic" = k_aic)
  k <- get_mode(ks)

  # plot to confirm
  p <- data.frame("AIC" = -aic, "TWGSS" = twgss, k = seq_along(twgss)) |>
    pivot_longer(AIC:TWGSS) |>
    ggplot(aes(x = k, y = value, colour = name)) +
    geom_line() +
    # geom_vline(xintercept = n_clust, linetype = "dashed") +
    facet_wrap(~name, scales = "free_y") +
    theme +
    labs(
      title = paste0(
        # "Elbow Plots for AIC and TWGSS, true k = ",
        # n_clust,
        "Elbow Plots for AIC and TWGSS",
        # ", ",
        # "LRT choice = ",
        # k_lr,
        ", DQU = ",
        cond_prob * 100,
        "%"
      ),
      x = "k",
      y = ""
    ) +
    scale_x_continuous(breaks = 1:max_clust, limits = c(1, max_clust)) +
    guides(colour = "none") +
    ggsci::scale_colour_nejm()

  # return optimal k, k estimates from each method, and plot
  return(list(
    "k"        = k,
    "k_method" = ks,
    "plot"     = p
  ))
}

# TODO: Remove elevation legend and add back in after
k_obj <- find_k(
  ce_fit_final,
  data_lst,
  max_clust = 10,
  cond_prob = dqu,
  fixed_b = FALSE,
  cond_var = c("rain", "wind_speed")
)
# saveRDS(k_obj, "k_obj.RDS")
# readRDS(k_obj, "k_obj.RDS")

# take k from TWGSS, looks like a pretty clear elbow at 3!
k_obj$plot
k_obj$k_method
k <- k_obj$k_method["k_twgss"]


#### Cluster and plot ####

# cluster for k = 3!
if (!exists("dist_mat")) {
  set.seed(123)
  clust_obj <- js_clust(
    dep_fit_final$dependence,
    scree_k = 1:max_k,
    n = n_mc
  )
  dist_mat <- clust_obj$dist_mat
}
k <- 3
pam_fit <- js_clust(dist_mat = dist_mat, k = k, n = n_mc)
# plot on map
# ggsave(filename = "test.png", plt_clust_map(pts, areas, pam_fit))
# TODO
plt_clust_map(
  pts, areas, pam_fit,
  plot_medoids = FALSE,
  elev_df = elev_df,
  rm_elev_leg = FALSE
)

# also look at k = 2, k = 4
# most sites in Western cluster, but Kilcar is noticably an outlier
p_k2 <- plt_clust_map(pts, areas, js_clust(dist_mat = dist_mat, k = 2, n = n_mc), rm_elev_leg = TRUE)
# Only Malahide in 4th cluster, interestingly!
p_k4 <- plt_clust_map(pts, areas, js_clust(dist_mat = dist_mat, k = 4, n = n_mc), rm_elev_leg = TRUE)

# join and save for supplementary materials
p_k_alt <- wrap_plots(list(
  p_k2 +
    NULL,
  p_k4 +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
))
ggsave(
  filename = "latex/plots/cluster_dqu_k_2_4.png",
  p_k_alt,
  width = 12,
  height = 8
)


#### Sensitivity analysis to choice of DQU ####

# sensible quantiles
# quantiles <- c(0.85, 0.88, 0.9, 0.92, 0.95)
quantiles <- c(0.85, 0.88, 0.9)

# function to fit CE model and cluster for different DQU
quant_fit <- \(q, k_spec = 3, ...) {
  print(paste0("Fitting DQU = ", q * 100, "%"))
  ce_fit_spec <- lapply(Y_lst, \(x) {
    o <- ce_optim(
      Y = x,
      dqu = q,
      cond_var = c("rain", "wind_speed"),
      control = list(maxit = 1e6),
      constrain = TRUE,
      # aLow = 0,
      start = c(0.01, 0.01),
      nruns = 3
    )
  })
  locs <- names(ce_fit_spec)

  dep_fit_spec <- lapply(c("resid", "params"), \(x) {
    setNames(pull_element(ce_fit_spec, x), locs)
  })
  names(dep_fit_spec) <- c("residual", "dependence")

  # cluster for k = 3
  # pam_fit_spec <- js_clust(
  #   dep_fit_spec$dependence,
  #   trans = Y_lst,
  #   # k = 3,
  #   k = k_spec,
  #   scree_k = list(1:5, 1:5),
  #   n = n_mc,
  #   mc_method = mc_method,
  #   laplace_cap = laplace_cap,
  #   return_dist = TRUE,
  #   ...
  # )

  # TODO Investigate why dist_obj$total_within_ss is only of length two??
  dist_obj <- js_clust(
    dep_fit_spec$dependence,
    trans = Y_lst,
    # k = 3,
    # k = k_spec,
    k = NULL,
    scree_k = list(1:5, 1:5), # TODO Why length 2? 2 variables?
    n = n_mc,
    mc_method = mc_method,
    laplace_cap = laplace_cap,
    return_dist = TRUE,
    ...
  )

  # pull out distance matrices for
  dist_mats <- list(
    dist_obj$dist_mat,
    dist_obj$dist_mats
  )

  twgss <- c(
    list(scree_plot(dist_mats[[1]], k = 1:5)), # TWGSS for combined diss mat
    dist_obj$total_within_ss
  )


  # only keep PAM object
  pam_fit_spec <- js_clust(
    dist_mat = dist_obj$dist_mat,
    dep_fit$dependence,
    k = k
  )

  #

  # plot
  map_plot_spec <- plt_clust_map(
    pts, areas, pam_fit_spec,
    plot_medoids = FALSE,
    # plot_medoids = TRUE,
    elev_df = elev_df
  ) +
    # ggtitle(paste0("DQU = ", q * 100, "%")) +
    ggtitle(paste0(q * 100, "%")) +
    guides(colour = "none", fill = "none")

  return(list(
    "ce_fit" = ce_fit_spec,
    "dep_fit" = dep_fit_spec,
    "pam_fit" = pam_fit_spec,
    "map_plot" = map_plot_spec
  ))
}

fit_quant <- lapply(quantiles, quant_fit)
saveRDS(fit_quant, file = "fit_quant.RDS")
fit_quant <- readRDS("fit_quant.RDS")

# TODO Add elevation to these plots!!
# p_sens <- wrap_plots(lapply(fit_quant, `[[`, "map_plot")[1:3])
p_sens <- wrap_plots(lapply(seq_along(fit_quant), \(i) {
  p_ret <- fit_quant[[i]]$map_plot
  # remove latitude for all but first plot
  if (i > 1) {
    p_ret <- p_ret +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  return(p_ret)
}))

# seems to be pretty stable for 85-90!

# save plot
pdf("plots/tests/cluster_dqu.pdf", width = 8, height = 8)
lapply(fit_quant, `[[`, "map_plot")
dev.off()

# save
ggsave(
  plot = p_sens,
  filename = "latex/plots/cluster_dqu.png",
  width = 12, height = 8
)

# also add medoids
add_medoid <- \(fit_quant) {
  lapply(seq_along(fit_quant), \(i) {
    q <- quantiles[[i]]
    plt_clust_map(
      pts,
      areas,
      # pam_fit_spec,
      fit_quant[[i]]$pam_fit,
      # plot_medoids = FALSE
      plot_medoids = TRUE
    ) +
      ggtitle(paste0("DQU = ", q * 100, "%")) +
      guides(colour = "none", fill = "none")
  })
}
plt_quant_med <- add_medoid(fit_quant)

pdf("plots/tests/cluster_dqu_medoid.pdf", width = 8, height = 8)
plt_quant_med
dev.off()


#### Cluster each variable separately ####

vars <- c("rain", "wind_speed")
# quantiles <- 0.85
fit_quant_sep <- lapply(vars, \(x) {
  print(paste0("Fitting for ", x))
  lapply(quantiles, quant_fit, spec_vars = x)
})
saveRDS(fit_quant_sep, file = "fit_quant_sep.RDS")
fit_quant_sep <- readRDS("fit_quant_sep.RDS")

# combined for rain and wind speed seperately
(p_rain <- wrap_plots(lapply(fit_quant_sep[[1]], `[[`, "map_plot")[1:3]) +
  patchwork::plot_annotation(title = "Rain", theme = theme))
# For rain, DQU 85 seems to have best results, somewhat less stable thereafter
(p_wind <- wrap_plots(lapply(fit_quant_sep[[2]], `[[`, "map_plot")[1:3]) +
  patchwork::plot_annotation(title = "Wind Speed", theme = theme))
# For wind speed, same, but median is a bit weird!

# plot for all
print_quant_plt <- \(fit_quant, fit_quant_sep, quantiles) {
  lapply(quantiles, \(q_spec) {
    x <- which(quantiles == q_spec)
    p <- wrap_plots(list(
      fit_quant[[x]]$map_plot +
        ggtitle("Both"),
      fit_quant_sep[[1]][[x]]$map_plot +
        # ggtitle("Rain | Wind Speed") +
        ggtitle("Precipitation | Wind Speed") +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()),
      fit_quant_sep[[2]][[x]]$map_plot +
        # ggtitle("Wind Speed | Rain") +
        ggtitle("Wind Speed | Precipitation") +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    )) +
      patchwork::plot_annotation(
        paste0("DQU = ", q_spec * 100, "%"),
        theme = evc::evc_theme()
      ) +
      NULL
  })
}

pdf("plots/tests/cluster_dqu_sep.pdf", width = 10, height = 10)
print_quant_plt(fit_quant, fit_quant_sep, quantiles)
dev.off()

# add medoids
fit_quant_sep_med <- lapply(fit_quant_sep, add_medoid)
pdf("plots/tests/cluster_dqu_sep_medoid.pdf", width = 10, height = 10)
print_quant_plt(
  lapply(plt_quant_med, \(x) list("map_plot" = x)),
  lapply(fit_quant_sep_med, \(x) lapply(x, \(y) list("map_plot" = y))),
  quantiles
)
dev.off()

# combine all three for DQU = 0.85, our chosen DQU
x <- which(quantiles == dqu)
p <- wrap_plots(list(
  fit_quant[[x]]$map_plot +
    ggtitle("Both"),
  fit_quant_sep[[1]][[x]]$map_plot +
    # ggtitle("Wind Speed | Rain") +
    ggtitle("Wind Speed | Precipitation") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()),
  fit_quant_sep[[2]][[x]]$map_plot +
    # ggtitle("Rain | Wind Speed") +
    ggtitle("Precipitation | Wind Speed") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
)) +
  # patchwork::plot_annotation("DQU = 85%", theme = evc::evc_theme()) +
  NULL

# TODO Improve plot!
ggsave(
  # filename = "latex/plots/cluster_dqu_sep.png",
  # filename = "plots/tests/cluster_dqu_sep_new.png",
  filename = "latex/plots/cluster_dqu_sep_new.png",
  p,
  width = 12,
  height = 8
)
# TODO Create 9 map plots for each likely DQU


##### Adjacency ####

# TODO Fix, adjacency matrix is not correct!

# dist_mat_adj <- as.matrix(dist_mat)
# # ensure order of columns and rows is the same as adjacency matrix
# ord <- match(rownames(adj_mat), colnames(dist_mat_adj))
# dist_mat_adj <- dist_mat_adj[ord, ord]
#
# # double check
# all(rownames(dist_mat_adj) == rownames(adj_mat))
# all(colnames(dist_mat_adj) == colnames(adj_mat))
# all(rownames(dist_mat_adj) == colnames(dist_mat_adj))
#
# dist_mat_adj[adj_mat == 0] <- 1e9 # arbitrarily large number
# # diagonal entries need to be 0
# diag(dist_mat_adj) <- 0

# Compute distance matrix using adjacency stipulation
dist_mat_adj <- as.matrix(dist_mat)
dist_mat_adj[adj_mat == 0] <- 10^40
dist_mat_adj <- as.dist(dist_mat_adj)

# apply adjacency matrix
k_obj_adj <- find_k(
  ce_fit,
  data_lst,
  max_clust = 10,
  cond_prob = dqu,
  fixed_b = FALSE,
  cond_var = c("rain", "wind_speed"),
  adj_mat = adj_mat
)
k_obj_adj$plot
k_obj_adj$k_method
# no clear elbow anywhere

k_adj <- 3
pam_fit_adj <- js_clust(dist_mat = dist_mat_adj, k = k_adj, n = n_mc)

# plot on map
# TODO Check that point ordering is definitely correct!
# ggsave("test_adj.png", plt_clust_map(pts, areas, pam_fit_adj))
plt_clust_map(pts, areas, pam_fit_adj)

# also look at k = 2, k = 4
plt_clust_map(pts, areas, js_clust(dist_mat = dist_mat_adj, k = 2, n = n_mc))
plt_clust_map(pts, areas, js_clust(dist_mat = dist_mat_adj, k = 4, n = n_mc))


#### Refit CE, plot and bootstrap uncertainty ####

# TODO Decide how best to "combine" information across cluster;
# should days be added together (i.e. sum of rain & max wind), or
# should we just use as more data??

# reassign points to their cluster medoids
clust_assign_df <- data.frame(
  "name" = names(pam_fit$clustering), "clust" = pam_fit$clustering
)
data_week_clust <- data_week |>
  left_join(clust_assign_df) |>
  mutate(medoid = ifelse(name %in% pam_fit$medoids, TRUE, FALSE))

med_loc <- data_week_clust |>
  filter(medoid == TRUE) |>
  distinct(clust, lon, lat)

data_week_clust <- data_week_clust |>
  select(-c(lon, lat)) |>
  left_join(med_loc) |>
  mutate(name = paste0("cluster_", clust)) |>
  group_by(name, date, lon, lat) |>
  summarise(
    rain = sum(rain, na.rm = TRUE),
    wind_speed = max(wind_speed, na.rm = TRUE),
    .groups = "drop"
  )

data_lst_clust <- data_week_clust |>
  group_split(name, .keep = TRUE) # keep name to label list afterwards
locs <- purrr::map_chr(data_lst_clust, ~ as.character(.x$name[1]))
names(data_lst_clust) <- locs
data_lst_clust <- lapply(data_lst_clust, \(x) {
  matrix(
    c(x$rain, x$wind_speed),
    ncol = 2,
    dimnames = list(NULL, c("rain", "wind_speed"))
  )
})

# Bootstrap choice of DQU
quant_plots_clust <- plot_boot_quant(
  data_lst_clust,
  quantiles = sort(c(seq(0.5, 0.9, by = 0.1), 0.85, 0.88, 0.92)),
  constrain = TRUE,
  # R = 10,
  R = 30,
  # ncores = ncores,
  ncores = 1,
  # start = c(0.01, 0.1),
  # fixed_b = TRUE,
  aLow = 0
  # county_key_df = county_key_df
)

# # seems that things are stable up to around this dependence quantile (??)
# # dqu_clust <- 0.85
# Maybe just stick with previous dependence quantile, to show reduction in
# uncertainty
dqu_clust <- dqu
dqu_clust <- 0.88

# Refit conditional extremes
# TODO Check if constraining is a good idea!
Y_lst_clust <- lapply(data_lst_clust, trans_fun)
ce_fit_clust <- lapply(Y_lst_clust, \(x) {
  o <- ce_optim(
    Y = x,
    dqu = dqu_clust,
    cond_var = c("rain", "wind_speed"),
    control = list(maxit = 1e6),
    constrain = TRUE,
    # start = c(0.01, 0.1),
    # fixed_b = TRUE,
    nruns = 3
    # aLow = 0
  )
})

# Diagnostic plots
dep_fit_clust <- lapply(c("resid", "params"), \(x) {
  setNames(pull_element(ce_fit_clust, x), locs)
})
names(dep_fit_clust) <- c("residual", "dependence")
dep_fit_clust$transformed <- Y_lst_clust
dep_fit_clust$original <- data_week_clust
dep_fit_clust$arg_vals <- list("cond_prob" = dqu)

# produce residual plots for all locations
resid_plots <- resid_plot(dep_fit_clust)

# produce quantile plots in all locations
quantile_plots <- quant_plot(dep_fit_clust)

# join resid and quantile plots
diag_plots <- diag_plot(resid_plots, quantile_plots)
names(diag_plots) <- names(resid_plots)

# save first plot to show
ggsave(
  filename = "latex/plots/diag_plots_postclust_new.png",
  diag_plots[[1]],
  width = 8,
  height = 8
)

# Plot alpha, beta
ab_df_clust <- dep_to_df(dep_fit_clust$dependence)

map_plots_clust <- map_plot(
  ab_df_clust,
  data_week_clust,
  n_breaks = 6,
  # range = c(5, 10),
  range = c(4, 10),
  elev_df = elev_df
)

ggsave(
  # map_plots_clust[[1]] / map_plots_clust[[2]],
  (map_plots_clust[[1]] + guides(fill = "none", size = "none")) /
    map_plots_clust[[2]],
  filename = "latex/plots/ab_map_post_clust_new.png",
  width = 12, height = 12
)

# TODO Don't plot on map!
ab_df_clust |>
  # filter(parameter == "a") |>
  # ggplot(aes(x = value, fill = vars)) +
  ggplot(aes(y = value, x = name)) +
  geom_dotplot(aes(fill = vars), binaxis = "y", stackdir = "center", binwidth = 0.05) +
  # facet_wrap(~ parameter + vars, scales = "free") +
  facet_wrap(~ parameter + vars, scales = "free") +
  labs(x = "") +
  theme +
  guides(fill = "none")

# bootstrap conditional extremes
boot_fit_clust <- boot_ce_ecdf(
  dep = dep_fit_clust$dependence,
  orig = data_lst_clust,
  transformed = Y_lst_clust,
  cond_prob = dqu_clust,
  # R = 100,
  R = 500,
  trace = 10,
  ncores = ncores,
  constrain = TRUE,
  # fixed_b = TRUE,
  aLow = 0
)

# TODO Plot (hopefully reduced!) boostrap uncertainty
boot_df_clust <- bind_rows(boot_fit_clust, .id = "run") |>
  group_by(name, parameter, vars) |>
  mutate(mean = mean(value, na.rm = TRUE)) |>
  ungroup()

boot_df_clust |>
  # filter(parameter == "a") |>
  ggplot(aes(x = vars, y = value, fill = vars)) +
  geom_boxplot(aes(fill = vars), alpha = 0.7) +
  # facet_wrap(~ parameter + vars) +
  # facet_wrap(~ parameter + vars, scales = "free") +
  facet_wrap(~ parameter + name, scales = "free") +
  theme +
  # guides(fill = "none") +
  NULL

boot_df_clust |>
  # filter(parameter == "a") |>
  ggplot(aes(x = value, fill = vars)) +
  geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.7) +
  geom_density(alpha = 0.5) +
  # facet_wrap(~ parameter + vars, scales = "free") +
  facet_wrap(~ parameter + name, scales = "free") +
  theme +
  # guides(fill = "none") +
  NULL

# Plot idea:
# "dumbell" plot of bootstrapped a and b values for each location,
# with cluster-wise uncertainty at bottom

plot_df <- boot_df |>
  left_join(clust_assign_df) |>
  bind_rows(boot_df_clust) |>
  mutate(
    name_clust = stringr::str_remove(name, "cluster_"),
    name = ifelse(
      grepl("cluster_", name),
      stringr::str_replace(stringr::str_to_title(name), "_", " "),
      name
    ),
    clust = ifelse(is.na(clust), name_clust, clust)
  ) |>
  # for each location, estimate 95% CI for a and b
  group_by(name, clust, parameter, vars) |>
  summarise(
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975),
    mean = mean(value),
    .groups = "drop"
  ) |>
  mutate(diff = upper - lower)

# pull mean uncertainty by cluster
plot_df <- plot_df |>
  bind_rows(
    plot_df |>
      group_by(clust, parameter, vars) |>
      summarise(
        diff = mean(diff),
        .groups = "drop"
      ) |>
      mutate(name = paste0("Mean (Cluster ", clust, ")"))
  ) |>
  # use Markdown bold for cluster rows
  mutate(
    name = ifelse(
      grepl("Cluster", name), glue("<b>{name}</b>"), name
    )
  )

# TODO Change colour depending on whether it's above/below clust diff
cols <- c("black", rev(ggsci::pal_nejm()(2)))
boot_comp_plots <- plot_df |>
  group_split(clust, parameter, vars) |>
  lapply(\(x) {
    # pull names for factor levels
    names_spec <- unique(x$name)
    (clust_name <- names_spec[grepl("Cluster", names_spec)][1])
    mean_name <- names_spec[grepl("Mean", names_spec)]
    names_spec <- names_spec[!names_spec %in% c(clust_name, mean_name)]

    # pull size of CI for clustered fit, want to add to plot
    clust_diff <- x |>
      filter(name == clust_name) |>
      pull(diff)

    x_plot <- x |>
      mutate(
        name = factor(
          name,
          levels = c(clust_name, mean_name, sort(names_spec, decreasing = TRUE))
        ),
        lower_ci = factor(case_when(
          diff == clust_diff ~ 1,
          diff > clust_diff ~ 2,
          TRUE ~ 3
        )),
        clust = paste0("Cluster ", clust, ", ", parameter, " (", vars, ")"),
      )

    lower <- sum(x_plot[!grepl("Cluster", x_plot$name), "lower_ci"] == 2)
    x_plot$clust <- paste0(
      x_plot$clust,
      " (",
      lower,
      "/",
      nrow(x_plot) - 2,
      " lower 95% CI)"
    )

    cols_spec <- cols
    if (lower == 0) {
      cols_spec <- cols_spec[-2]
    }

    # browser()
    # print(paste0("For ", x_plot$clust[1], ", ", sum(x_plot$lower_ci == 2), "/", nrow(x_plot) - 1))

    x_plot |>
      ggplot() +
      geom_segment(
        aes(
          x      = name,
          # y      = lower,
          y      = 0,
          # yend   = upper
          yend   = diff,
          colour = lower_ci
        )
      ) +
      # set colours to colour palette specified above
      scale_colour_manual(values = cols_spec) +
      geom_hline(
        data = x_plot |>
          filter(
            name == clust_name,
            vars == x$vars[1],
            parameter == x$parameter[1]
          ),
        aes(yintercept = diff),
        linetype = "dashed"
      ) +
      facet_wrap(~clust) +
      coord_flip(clip = "off", expand = TRUE) + # TODO Do I want to expand??
      theme +
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = ggtext::element_markdown()
      ) +
      guides(colour = "none")
  })

# pdf("plots/tests/boot_compare.pdf", width = 8, height = 8)
pdf(
  paste0("plots/tests/boot_compare_dqu_", dqu_clust, ".pdf"),
  width = 8,
  height = 8
)
boot_comp_plots
dev.off()

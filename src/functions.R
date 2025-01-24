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
  for (i in seq_len(D)){
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
  return(beta.all/D)
}

# TODO: Add shared plotting functions

# function to calculate confidence interval for sensitivity analysis results
summarise_sens_res <- \(results_grid, conf_level = 0.95) {
  # z-score corresponding to confidence level
  z <- qnorm((1 + conf_level) / 2)
  
  results_grid %>% 
    # remove any previous summaries as they mess up grouping
    dplyr::select(-(ends_with("_rand") & !matches("adj_rand"))) %>% 
    relocate(adj_rand, .after = everything()) %>% # ensure is last col
    group_by(across(1:last_col(1))) %>% 
    summarise(
      mean_rand  = mean(adj_rand, na.rm = TRUE), 
      # (z score corresponding to confidence level) * standard error
      marg_rand  = z * (sd(adj_rand, na.rm = TRUE) / sqrt(n())),
      # CIs
      lower_rand = mean_rand - marg_rand,
      upper_rand = mean_rand + marg_rand,
      .groups    = "drop"
    ) %>% 
    dplyr::select(-marg_rand) %>%  # not required for plotting
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
            wind_dir, type = "angles", units = "degrees", modulo = "2pi"
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
        wind_dir, type = "angles", units = "degrees", modulo = "2pi"
      ))),
      .groups = "drop"
    ) %>% 
    pivot_longer(dist2coast:wind_dir, names_to = "var", values_to = "value") %>% 
    ggplot() +
    geom_bar(
      aes(x = factor(clust), y = value, fill = factor(clust)),
      stat = "identity"
    ) +
    facet_wrap(~ var, scales = "free") + 
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
    facet_wrap(~ var, scales = "free") + 
    evc_theme() + 
    labs(fill = "Cluster") + 
    guides(fill = guide_legend(override.aes = list(alpha = 1))) + 
    ggsci::scale_fill_nejm()
  
  return(list("hist" = p1, "density" = p2)) 
}


#### Simulation functions ####

# Function to generate copula data
# TODO: Extend to simulate more than two variables
# TODO: Extend to simulate more than two clusters
sim_cop_dat <- \(
  n_locs      = 12,    # number of locations
  n_vars      = 2,     # number of variables at each location
  n           = 1e4,   # number of simulations for each variable
  n_clust     = 2,     # number of clusters
  cluster_mem = NULL,  # desired membership, if NULL then evenly split
  cor_gauss,           # bulk correlation for n_clust clusters from Gaussian copula
  # params_norm,         # normal marginal parameters (same for both)
  cor_t,               # extreme correlation for n_clust clusters from t-copula
  df_t,                # degrees of freedom for t-copula
  params_gpd,          # GPD margin parameters
  mix_p = c(0.5, 0.5), # mixture weights
  perturb_cor = FALSE, # perturb correlation for each location within clusters
  perturb_val = 0.05   # value to perturb by if desired
) {
  # many arguments must be of length n_clust
  # stopifnot(all(vapply(
  #   # list(cor_gauss, cor_t, df_t, params_gpd, mix_p), 
  #   list(cor_gauss, cor_t, df_t), 
  #   \(x) length(x) == n_clust, logical(1)
  # )))
  # stopifnot(sum(mix_p) == 1)
  
  # Simulate from Gaussian Copula with GPD margins
  gauss_cop <- lapply(seq_len(n_locs), \(i){
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
    # transform to GPD margins
    evd::qgpd(
      p     = u,
      loc   = 0,
      scale = params_gpd[1], 
      shape = params_gpd[2] 
    )
  })
  
  # simulate from t-Copula with GPD margins
  # TODO: Functinoalise repeated behaviour here from above?
  t_cop <- lapply(seq_len(n_locs), \(i) {
    if (is.null(cluster_mem)) {
      group <- ceiling(i / (n_locs / n_clust))
    } else {
      group <- cluster_mem[i]
    }
    # Assign the corresponding value to the result
    cor <- cor_t[group]
    df  <- df_t[group]
    
    # optionally perturb correlation for each location within clusters
    if (perturb_cor) {
      cor <- cor + runif(1, -perturb_val, perturb_val)
      cor <- pmax(pmin(cor, 1), 0) # ensure within [0, 1]
    }
    
    cor <- rep(cor, n_vars * (n_vars - 1) / 2)
    cop_t <- copula::tCopula(cor, dim = n_vars, df = df, dispstr = "un")
    u <- copula::rCopula(n, cop_t)
    return(evd::qgpd(
      p     = u,
      loc   = 0,
      scale = params_gpd[1], 
      shape = params_gpd[2] 
    ))
  })
  
  # mixture
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
    lon_vals <-  nc_connection$dim$lon$vals
    lat_vals <- nc_connection$dim$lat$vals
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
  rm(data); gc()

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
    reference_date <- as.POSIXct(
      substr(nc_connection$dim$time$units, 13, 31), tz = "UTC"
    )
    ret$time <- reference_date + ret$time * 3600
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

#### Sensitivity analysis ####

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
    cop_norm <- copula::normalCopula(cor, dim = 2, dispstr = "un")
    # cop_norm <- normalCopula(cor_mat, dim = 2, dispstr = "un")
    u <- copula::rCopula(n, cop_norm)
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
    u <- copula::rCopula(n, cop_t)
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

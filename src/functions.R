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

  ret <- lapply(seq_along(data), \(i) {
    dimnames(data[[i]]) <- list(lon_vals, lat_vals, time_vals)
    melt(
      data[[i]],
      varnames = c("longitude", "latitude", "time"),
      value.name = names(data)[i]
    )
  })
  rm(data); gc()

  ret <- ret[[1]] %>% 
    bind_cols(select(ret[[2]], "v10")) %>%
    # bind_cols(select(ret[[3]], "msl")) %>%
    # filter(time == 1) %>% 
    # filter(time == 920424) %>% 
    # ggplot() + 
    # geom_raster(aes(x = longitude, y = latitude, fill = u10))
    identity()

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
  reference_date <- as.POSIXct(
    substr(nc_connection$dim$time$units, 13, 31), tz = "UTC"
  )
  ret$time <- reference_date + ret$time * 3600

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

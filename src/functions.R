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

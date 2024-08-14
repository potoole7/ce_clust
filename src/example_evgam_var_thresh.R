#### evgam spatial GPD fit ####

# following https://empslocal.ex.ac.uk/people/staff/by223/software/gpd.R
# In paper https://www.tandfonline.com/doi/full/10.1080/01621459.2018.1529596

library(chron)
library(fields)
library(evgam)
library(maps)

#### Variable location parameter for GPD ####

# data(USgusts)
load("data/USgusts.rda")
USgusts <- subset(USgusts, month %in% c("08", "09", "10")) # Autumn only
# USgusts <- USgusts %>% 
#   filter(
#     month %in%c("08")
#   )

# here we use a thin-plate spline to capture spatial variation in the
# quantile, i.e. the asymmetric Laplace distribution's location
# and scale parameters
fmla.ald <- list(
  gust ~ s(lon, lat, k = 50), # location
  ~ s(lon, lat, k = 40)       # logscale
)

# fit the quantile regression model using 97th percentile
# might take a minute or two
# fit.ald <- fevgam(fmla.ald, USgusts, family="ald", tau=.97)
fit.ald <- evgam(fmla.ald, USgusts, family="ald", ald.args = list(tau=.97))

# this can help identify redundant smooths or model misspecification
summary(fit.ald)
# plot smooths
plot(fit.ald)

#### Fit evgam GPD to location and shape parameters

# create a data frame with just threshold exceedances
USgusts$threshold <- predict(fit.ald)$location
USgusts$excess <- USgusts$gust - USgusts$threshold
USgusts.gpd <- subset(USgusts, excess > 0)

fmla.gpd <- list(
  excess ~ s(lon, lat, k=50), # scale
  ~ s(lon, lat, k=40) # shape
)
# fit.gpd <- fevgam(fmla.gpd, USgusts.gpd, family = "gpd")
fit.gpd <- evgam(fmla.gpd, USgusts.gpd, family = "gpd")
USgusts.gpd <- cbind(USgusts.gpd, predict(fit.gpd, type = "response"))

summary(fit.gpd)
plot(fit.gpd)


#### Return levels ####

# map of US
state <- map("state")

# Create square grid over US
lon.plot <- seq(min(state$x, na.rm =TRUE), max(state$x, na.rm = TRUE), l = 100)
lat.plot <- seq(min(state$y, na.rm =TRUE), max(state$y, na.rm = TRUE), l = 100)
df.plot <- expand.grid(lon = lon.plot, lat = lat.plot)

# Subset so as to only be over the continental US
over.land <- !is.na(map.where("state", df.plot$lon, df.plot$lat))
df.plot <- subset(df.plot, over.land)

# return period * no. days
return.period <- 100 * 92

#### 100-year return levels without extremal index (not advisable) ####

df.plot <- cbind(df.plot, predict(fit.gpd, df.plot, type="response"))
df.plot$threshold <- predict(fit.ald, df.plot, type="response")$location
# TODO: Check in Coles book
df.plot$rl1 <- df.plot$threshold + 
  df.plot$scale * (((1 - fit.ald$tau) * return.period) ^ df.plot$shape - 1) / 
  df.plot$shape

plot.mat1 <- matrix(NA, length(lon.plot), length(lat.plot))
plot.mat1[over.land] <- df.plot$rl1
image.plot(
  lon.plot, lat.plot, plot.mat1, xlab = "longitude", ylab = "latitude"
)
map("state", add = TRUE)


#### 100-year return levels with extremal index ####

# the following does a probability integral transform to convert to
# censored unit Frechet margins
USgusts$exi <- USgusts$excess > 0
USgusts$frech <- -1/log(fit.ald$tau)
USgusts$frech[
  USgusts$exi
] <- -1/log(1 - (1 - fit.ald$tau) * 
       (1 - pgpd(USgusts.gpd$excess, 0, USgusts.gpd$scale, USgusts.gpd$shape)))

# this creates an integer time variable 
# which is needed for calculating consecutive days' maxima
USgusts$chron <- as.numeric(
  chron(paste(USgusts$month, USgusts$day, USgusts$year, sep="/"))
)

# split data frame by weather station
# so that running max can be calculated for 
# each station
USgusts.exi <- split(USgusts, USgusts$WMO_id)

# replace Frechet variable with two-observation running max
USgusts.exi <- do.call(
  rbind, 
  lapply(USgusts.exi, dfrunmax, cons = "chron", ynm = "frech")
)

# estimate extremal index with 
# thin-plate spline for spatial variation
# and probit link
# smooth parameter estimate tends not to converge too well
# when there's only one smoothing parameter
# fit.exi <- fevgam(frech ~ s(lon, lat), USgusts.exi, family="exi", exilink="probit")
fit.exi <- evgam(
  frech ~ s(lon, lat), 
  USgusts.exi, 
  family = "exi", # extremal index estimation 
  # exilink = "probit"
  # exi.args = list(link = "probit") # needs ID?
  exi.args = list(
    id = "exi", # TODO: Wrong! Need to look into
    link = "probit"
  ) 
)

# plot estimated spline for extremal index (excludes constant)
plot(fit.exi)

# df.plot$exi <- predict(fit.exi, df.plot, type="response")$location
df.plot$exi <- predict(fit.exi, df.plot, type="response")$exi
df.plot$rl2 <- df.plot$threshold + 
  df.plot$scale * 
    (((1 - fit.ald$tau) * return.period * df.plot$exi)^df.plot$shape - 1) / 
      df.plot$shape

plot.mat2 <- matrix(NA, length(lon.plot), length(lat.plot))
plot.mat2[over.land] <- df.plot$rl2
image.plot(
  lon.plot, lat.plot, plot.mat2, xlab = "longitude", ylab = "latitude"
)
map("state", add = TRUE)

# compare 100-year return level estimates with and without extremal index
image.plot(plot.mat2 / plot.mat1)
map("state", add = TRUE)

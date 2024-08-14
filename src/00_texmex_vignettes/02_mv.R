#### texmex vignettes: Multivariate Conditional Extremes ####

# See Heffernan & Tawn 2004 for CE explanation
# Notes: 
# - Implements both Gumbel and Laplacian PITs
# - Has parameter constraints from Keef 2013 paper
# - Has really nice summaries for CE model!
# - Notice how fast fitting is! Granted, only 532 obs, but still!



#### Libs ####
library(texmex) # inc. rain, portporie, winter data
library(gridExtra)
palette(c("black", "purple", "cyan", "orange"))
set.seed(123)

#### Exploration ####

head(winter)
# columns: O3, N02, NO, SO2, PM10

# Dependence of different variables
GGally::ggpairs(winter)

# Extremal dependence
chiO3 <- chi(winter[, c("O3", "NO")]) # calculate measure of mv dependence
ggplot( # another ggplot method
  chiO3, 
  main = c("Chi"="Chi: O3 and NO", "ChiBar"="Chi-bar: O3 and NO")
)

# CI for chi-bar includes 0 as quantile -> 1 => asymp. independent

#### Fit model ####

# fit marginal and multivariate dependence models
mex.O3 <- mex( # wrapper around migpd and mexDependence
  winter, 
  mqu =.7, # marginal quantile for thresholding
  penalty = "none", 
  which = "O3" # condition on 03
)
mex.NO2 <- mex(winter, mqu =.7, penalty = "none", which = "NO2")
mex.NO <- mex(winter, mqu =.7, penalty = "none", which = "NO")
mex.SO2 <- mex(winter, mqu =.7, penalty = "none", which ="SO2")
mex.PM10 <- mex(winter, mqu =.7, penalty = "none", which = "PM10")

# first block gives marginal model results:
# - Threshold for each variable
# - Prob of exceeding threshold
# - GPD parameters
# - Upper end point
# Second block gives a, b, c, d, m and s for Depedence model
# - a and b give asymptotic dependence structure for pairs of variables
# - c and d not required for Laplace margins
# - what are m and s? Something to do with Gaussian assumption for residuals
# - m and s are the mean and sd for Z, ideally mean should be close to 0
# for Gaussian assumption
summary(mex.O3)

# shows centred + scaled residuals (Z), should be independent
ggplot(mex.O3)

# fit marginal GPDs
marg <- migpd(winter, mqu = 0.7, penalty = "none")
mrf <- mexRangeFit(marg, "NO", trace=11)
ggplot(mrf)
# Plots suggest should use different starting values for a, b!

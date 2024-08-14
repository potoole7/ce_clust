#### Test "hack" of texmex with ismev on Winter data ####

#### libs ####

library(dplyr, quietly = TRUE)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(ismev)
library(texmex)

# texmex_mqu <- 0.7 # from vignette
texmex_mqu <- 0.6
# results almost completely equal when texmex_mqu == ismev_mqu == 0.7
# ismev_mqu <- texmex_mqu
# ismev_mqu <- 0.8 # different so as to change results 
# ismev_mqu <- 0.65 # lower than texmex_mqu
ismev_mqu <- 0.8


#### texmex fit ####

# Fit "out-of-box" texmex marginal GPDs
texmex_marginal <- migpd(winter, mqu = texmex_mqu, penalty = "none")

# Fit corresponding CE models 
texmex_o3 <- mexDependence(texmex_marginal, which = "O3")
texmex_no2 <- mexDependence(texmex_marginal, which = "NO2")
texmex_no <- mexDependence(texmex_marginal, which = "NO")

# look at plot for O3
ggplot(texmex_o3)

#### "Hacked" fit using `ismev` ####

# fit ismev model for each
ismev_fits <- lapply(seq_len(ncol(winter)), \(i) {
  # fit GPD using `ismev`
  spec_gpd <- gpd.fit(
    xdat = winter[, i], 
    # Use slightly lower threshold to try get similar results
    # threshold = texmex_marginal$models[[i]]$threshold[[1]] - 3
    threshold = quantile(winter[, i], ismev_mqu)
  )
  # keep coefficients
  spec_gpd$mle[[1]] <- log(spec_gpd$mle[[1]])
  names(spec_gpd$mle) <- c("phi: ", "xi: ")
  return(spec_gpd$mle) 
})

# initialise marginal models fit using ismev
ismev_marginals <- texmex_marginal
ismev_marginals$mqu[] <- ismev_mqu
# replace values from texmex with those from ismev using slightly lower thresh
ismev_marginals$models <- lapply(seq_along(texmex_marginal$models), \(i) {
 x <- texmex_marginal$models[[i]] 
 x$threshold <- quantile(winter[, i], ismev_mqu)
 x$coefficients <- ismev_fits[[i]]
 return(x)
})
names(ismev_marginals$models) <- names(texmex_marginal$models)

ismev_o3 <- mexDependence(ismev_marginals, which = "O3")

#### Comparison ####

# compare a and b values
texmex_o3$dependence$coefficients
# ismev_o3$dependence$coefficients

# for same value (@0.7 and 0.8), parameter values are practically the same

# for texmex_mqu = 0.7, ismev_mqu = 0.8, sign of b for NO2, PM10 and SO2 and 
# a for NO changes!
# Model seems to be sensitive to change in mqu here, but may also just be because of 
# low parameter values and weak signal in this case
# Parameter value changes seem to make sense
# Could also be because of starting values for a and b?

# for ismev_mqu = 0.65, b for NO2, NO and PM10 changes signs
# However, b is very close to 0 for all anyway, so little difference in results, 
# and may just due to sensitity of analysis to threshold choice

# for ismev_mqu = 0.7, texmex_mqu = 0.8:
# a sign different for NO2, b different for NO2, NO, PM10
# again, value of b appears to be very sensitive to different thresholds

# compare diagnostic plots
ggplot(texmex_o3)
ggplot(ismev_o3)
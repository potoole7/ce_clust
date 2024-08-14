#### texmex vignettes: Marginal GPD ####

# Notes:
# GPD parameterised in terms of phi = log(sigma)
# plots have ggplot methods, interesting! Palette choice easy as well

#### Libs ####

library(texmex) # inc. rain, portporie, winter data
library(gridExtra)
palette(c("black", "purple", "cyan", "orange"))
set.seed(123)


#### Threshold Selection ####

# stability plot
grfRain <- gpdRangeFit(rain, umax = 35)
# mean residual life plot
mrlRain <- mrl(rain)

# Uses ggplot.gpdRangeFit(grfRain) method, how is this implementated?!?
g1 <- ggplot(grfRain) 
g2 <- ggplot(mrlRain)

grid.arrange(
  g1[[1]] + ggtitle("Stability plot, scale parameter"),
  g1[[2]] + ggtitle("Stability plot, shape parameter"), 
  g2 + ggtitle("MRL plot"), 
  ncol = 2
)

# selects threshold of ~20


#### Fit GPD ####

# Fit Extreme Value Model with default GPD family
(rain.fit <- evm(rain, th = 20))
ggplot(rain.fit)


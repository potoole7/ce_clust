#### Test bootstrapping vs `texmex` ####

# Testing implementation of boostrapping in evc vs texmex for texmex example

#### Libs ####

devtools::load_all("texmex")
devtools::load_all("../evc")
library(dplyr, warn.conflicts = FALSE)
library(ggplot2)
library(parallel)

ncores <- detectCores() - 1

source("src/functions.R")

#### Load Data ####

data <- texmex::winter

#### texmex ####

# fit mex object
tex_mex <- mex(winter , mqu = .7, dqu = .7, which = "NO")
# bootstrap
tex_boot <- bootmex(tex_mex)
# print and plot
tex_boot
par(mfrow = c(2, 5))
plot(tex_boot, plots = "gpd")
par(mfrow = c(2, 2))
plot(tex_boot, plots = "dependence")
par(mfrow = c(1, 1))

#### evc ####

# fit CE model using evc
# TODO Look at NO2, may be a bug there, estimates differ from texmex, rest fine
evc_fit <- fit_ce(
  data       = mutate(data, name = "location_1"), 
  vars       = names(data),
  f          = NULL, # fit with ismev
  marg_prob  = 0.7,
  cond_prob  = 0.7, 
  output_all = TRUE
)


#### bootstrap evc fit ####

# run bootstrap
bootstrap_res <- boot_ce(evc_fit, R = 100, trace = 10, ncores = ncores)

# ### Comparison ####

# density plot of xi and sigma values
bootstrap_res$marginal |> 
  left_join(
    as_tibble(summary(tex_boot$simpleMar)$co[c("sigma", "xi"), ]) |> 
    mutate(parameter = c("sigma", "xi")) |> 
    tidyr::pivot_longer(
      O3:PM10, names_to = "vars", values_to = "texmex_value"
    )
  ) |> 
  ggplot(aes(x = value, fill = parameter)) +
  # also add histograms
  geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.7) +
  geom_density(alpha = 0.5) +
  geom_vline(aes(xintercept = texmex_value), linetype = "dashed") +
  facet_wrap(~ parameter + vars, scales = "free", nrow = 2) +
  # facet_grid(vars ~ parameter, scales = "free") +
  evc::evc_theme() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
  
# density plot of a and b for NO | O3 (also add bootmex estimates)
bootstrap_res$dependence |> 
  filter(cond_var == "NO") |> 
  left_join(
    as_tibble(tex_boot$simpleDep[c(1, 2), ]) |> 
      mutate(parameter = c("a", "b")) |> 
      tidyr::pivot_longer(
        O3:PM10, names_to = "vars", values_to = "texmex_value"
      )
  ) |> 
  ggplot(aes(x = value, fill = parameter)) + 
  geom_histogram(aes(y = ..density..), bins = 20, alpha = 0.8) +
  geom_density(alpha = 0.4) + 
  geom_vline(aes(xintercept = texmex_value), linetype = "dashed") +
  facet_wrap(~ parameter + vars, scales = "free", nrow = 2) + 
  evc::evc_theme() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))


#### Test for > 1 location ####

# perturb data slightly and add as second "location"
data_mult <- bind_rows(
  mutate(data, name = "location_1"),
  mutate(
    data, 
    across(everything(), ~ jitter(.x, amount = 0.1)), # perturb slightly
    name = "location_2"
  )
)

# fit CE
# start values for each location and variable
start_mat <- matrix(c(0.1, 0.1), nrow = 2, ncol = ncol(data) - 1)
rownames(start_mat) <- c("a", "b")
locs <- unique(data_mult$name)
start <- lapply(locs, \(loc) {
  vars <- names(data)
  ret <- lapply(vars, \(var) {
    # start value for each dependence model for 
    ret <- start_mat
    colnames(ret) <- vars[vars != var]
    ret
  })
  names(ret) <- vars
  return(ret)
})
names(start) <- locs
start

evc_fit_mult <- fit_ce(
  data       = data_mult,
  vars       = names(data),
  start      = start,
  marg_prob  = 0.7,
  cond_prob  = 0.7, 
  output_all = TRUE
)

# bootstrap
debugonce(boot_ce)
bootstrap_res_mult <- boot_ce(
  evc_fit_mult, R = 100, trace = 10, ncores = ncores
)

bootstrap_res_mult$marginal |> 
  # filter(name == "location_1") |> 
  ggplot(aes(x = value, fill = parameter)) +
  # also add histograms
  geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.7) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ name + parameter + vars, scales = "free") + # , nrow = 2) +
  # facet_wrap(~ parameter + vars, scales = "free", nrow = 2) + 
  evc::evc_theme() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

bootstrap_res_mult$dependence |> 
  # filter(name == "location_1") |> 
  filter(cond_var == "NO") |>
  ggplot(aes(x = value, fill = parameter)) +
  # also add histograms
  geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.7) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ name + parameter + vars, scales = "free") + # , nrow = 2) +
  # facet_wrap(~ parameter + vars, scales = "free", nrow = 2) + 
  evc::evc_theme() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
 
 

#### Test with varying thresholds and varying marginal parameters ####
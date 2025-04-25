#### Testing different marginal thresholds ####

# - need to run soe of 02_fit_ce_irl_paper.R first!

#### Metadata ####

mqu_rain <- 0.9
mqu_ws <- c(0.9, 0.95, 0.96, 0.97, 0.98)
fixed_xi_vec <- c(FALSE, TRUE)

meta_df <- tidyr::crossing(
  mqu_rain = mqu_rain,
  mqu_ws = mqu_ws,
  fixed_xi = fixed_xi_vec
)

#### test ####

# should probably do plots separate to model fits but too lazy!
plot_fits <- lapply(seq_len(nrow(meta_df)), \(i) {
  system(sprintf(
    'echo "\n%s\n"',
    paste0(round(i / nrow(meta_df), 3) * 100, "% completed", collapse = "")
  ))

  # pull out metadata
  mqu <- unlist(meta_df[i, c("mqu_rain", "mqu_ws")])
  fixed_xi <- meta_df$fixed_xi[i]

  # fit marginal model with specified arguments
  f_evgam <- list(excess ~ s(lon, lat), ~ s(lon, lat))
  if (fixed_xi == TRUE) {
    f_evgam[[2]] <- ~1 # model xi with intercept (i.e. fixed value) only
  }
  shared_args <- list(
    data = data_week,
    vars = c("rain", "wind_speed"),
    # arguments to `quantile_thresh`
    # marg_prob = c(0.9, 0.95),
    marg_prob = list(
      f = list("response ~ s(lon, lat)", "~ s(lon, lat)"),
      tau = mqu
    ),
    # formula for evgam fit to excesses over threshold
    f = f_evgam,
    ncores = max(1, parallel::detectCores() - 2),
    marg_only = TRUE
  )
  # fit model
  mod_fit_qgam <- do.call(
    fit_ce, c(shared_args, list(thresh_fun = qgam_thresh))
  )
  mod_fit <- mod_fit_qgam

  # pull out estimates of marginals
  data_rain <- mod_fit$data_thresh$rain
  data_ws <- mod_fit$data_thresh$wind_speed

  data_rain_map <- data_rain %>%
    left_join(mod_fit$evgam_fit$rain$predictions) %>%
    st_to_sf(remove = FALSE) # add location
  data_ws_map <- data_ws %>%
    left_join(mod_fit$evgam_fit$wind_speed$predictions) %>%
    st_to_sf(remove = FALSE) # add location

  # generate PP-QQ plots, label appropriately and return
  set.seed(seed_number)
  pp_qq_rain <- (
    pp_qq_plot(data_rain_map) |>
      lapply(\(x) x + labs(x = "")) |>
      wrap_plots()
  ) +
    patchwork::plot_annotation(title = "Rain", theme = theme)

  pp_qq_ws <- (
    pp_qq_plot(data_ws_map) |>
      wrap_plots()
  ) +
    patchwork::plot_annotation(title = "Wind Speed", theme = theme)

  # p <- pp_qq_rain / pp_qq_ws +
  #   patchwork::plot_annotation(
  #     title = paste0(
  #       "mqu_rain: ", mqu[1],
  #       ", mqu_ws: ", mqu[2],
  #       ", fixed_xi: ", fixed_xi
  #     ),
  #     theme = theme
  #   )

  # return(p)
  # return(list("p" = p, "mod_fit" = mod_fit))
  return(list(
    "p" = list(
      pp_qq_rain,
      pp_qq_ws,
      title = paste0(
        "mqu_rain: ", mqu[1],
        ", mqu_ws: ", mqu[2],
        ", fixed_xi: ", fixed_xi
      )
    ),
    "mod_fit" = mod_fit
  ))
})

plots <- lapply(plot_fits, `[[`, "p")

# only keep rain plots for first two plots, then only look at wind
plots <- lapply(seq_along(plots), \(i) {
  if (i <= 2) {
    p <- plots[[i]][[1]] / plots[[i]][[2]] +
      patchwork::plot_annotation(
        title = plots[[i]][[3]],
        theme = theme
      )
  } else {
    p <- plots[[i]][[2]] +
      patchwork::plot_annotation(title = plots[[i]][[3]])
  }
  return(p)
})

# save plots
pdf(
  file = "plots/tests/pp_qq_marg_vary_mqu_fixed_xi.pdf",
  width = 8,
  height = 8
)
plots
dev.off()

# investigate particularly poor location for wind speed
# fit_spec <- plot_fits[[6]]$mod_fit # mqu_ws = 0.96, fixed_xi = TRUE
fit_spec <- mod_fit

data_ws <- fit_spec$data_thresh$wind_speed

data_ws_map <- data_ws %>%
  left_join(fit_spec$evgam_fit$wind_speed$predictions) %>%
  st_to_sf(remove = FALSE)

# plot QQ
p <- pp_qq_plot(data_ws_map)[[2]]

# pull name of location where empirical value is very high!
# locations where empirical value is quite high
(locs <- pp_qq_plot(data_ws_map)[[2]]$data |>
  arrange(desc(x)) |>
  head(5))

# specific location where it's bad
# loc <- "Castlereagh"
# loc <- "Marble Arch Caves"
# loc <- "Ringsend"
# loc <- first(locs$name)
loc <- "Crolly (Filter Works)"

# look at QQ plot without problem location (looks really good!)
pp_qq_plot(filter(data_ws_map, name != loc))[[1]]
pp_qq_plot(filter(data_ws_map, name != loc))[[2]]

# investigate QQ plot for just this location
p_loc <- pp_qq_plot(filter(data_ws_map, name == loc))
wrap_plots(p_loc)
# seems to be one crazy point really! Can also investigate this
p_loc[[2]]$data |>
  arrange(desc(x)) |>
  head(3)

# plot rain, windspeed and rain vs wind speed for Ringsend
# data_winter |>
data_week |>
  filter(name == loc) |>
  # filter(date > "2002-01-01" & date < "2003-01-01") |>
  pivot_longer(c(rain, wind_speed), names_to = "var") |>
  ggplot(aes(x = date, y = value)) +
  geom_point() +
  facet_wrap(~var, scales = "free_y") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  # scale_x_date(date_breaks = "1 month", date_labels = "%Y %m") +
  ggtitle(paste(
    loc,
    data_week$county[data_week$name == loc][1],
    data_week$province[data_week$name == loc][1],
    sep = " - "
  )) +
  theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

data_week %>%
  # add 95% quantile lines for both
  filter(name == loc) |>
  mutate(across(
    c(rain, wind_speed), ~ quantile(.x, 0.95, na.rm = TRUE),
    .names = "quant_{.col}"
  )) %>%
  ggplot(aes(x = rain, y = wind_speed)) +
  geom_point(aes(colour = name), size = 1.5, alpha = 0.8) +
  geom_vline(aes(xintercept = quant_rain), linetype = "dashed") +
  geom_hline(aes(yintercept = quant_wind_speed), linetype = "dashed") +
  facet_wrap(~name, scales = "free_x") +
  scale_colour_manual(values = c(ggsci::pal_nejm()(2))) +
  labs(
    x      = "precipitation (mm)",
    y      = "wind speed (m/s)",
    colour = ""
  ) +
  theme +
  # remove facet labels, colour will do
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.key = element_blank()
  )

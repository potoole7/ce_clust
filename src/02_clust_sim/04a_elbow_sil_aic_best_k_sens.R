#### Simulation sensitivity analysis of Jensen-Shannon silhouette ####

# Testing silhouette, elbow and AIC for choosing number of clusters

#### libs ####

library(dplyr, quietly = TRUE)
library(evc)
library(tidyr)
library(ggplot2)
library(evd)
library(evgam)
devtools::load_all("texmex")
library(cluster)
library(lcmix)
library(copula)
library(parallel)

# source functions
source("src/functions.R")


#### Metadata ####

seed_number <- 123
n <- 1e3 # number of samples to take
n_locs <- 60 # number of "locations" (needs got be higher for more clusters)
n_clust <- 11
# cluster_mem <- ceiling(1:n_locs / (n_locs / n_clust))
n_vars <- 2 # number of variables per location
mix_p <- c("gauss_cop" = 0.5, "t_cop" = 0.5) # mixture percentages
# GPD parameters
scale_gpd <- 1
shape_gpd <- -0.05
marg_prob <- 0.9
kl_prob <- 0.9
conf_level <- 0.95 # confidence level for CIs in plot

# Number of cores to use for parallel computation
n_cores <- detectCores() - 1 


#### Grid search ####

grid <- bind_rows(lapply(2:n_clust, \(i) {
  data.frame(t(c(seq(0, 1, length.out = i), rep(NA, times = n_clust - i))))
})) %>% 
  setNames(paste0("cor_t", seq_len(n_clust))) %>% 
  mutate(
    cor_gauss = 0.1, # 0.3, # keep constant
    # Degrees of freedom for t-copula
    df_t   = 3,
    # mixture percentages (must sum to 1)
    mix_p  = 0.5,
    # extremal quantiles
    kl_prob = kl_prob
  )

# run kl_sim_eval for each row in grid
n_times <- 5
elb_vec <- sil_vec <- aic_vec <- vector(length = n_times)
set.seed(seed_number)
results_grid <- bind_rows(mclapply(seq_len(nrow(grid)), \(i) {
  
  print(paste0("Progress: ", round(i / nrow(grid), 3) * 100, "%"))
  system(sprintf(
    'echo "\n%s\n"', 
    paste0(round(i / nrow(grid), 3) * 100, "% completed", collapse = "")
  ))
  
  # pull specific row
  row <- grid[i, , drop = FALSE]
  # pull correlations in row which are non NA
  cor_t <- row[, grepl("cor_t", names(row))]
  cor_t <- cor_t[!is.na(cor_t)]
  n_cor <- length(cor_t)
  print(cor_t)
  for (j in seq_len(n_times)) {
    
    # generate simulation data for given parameter set
    data_mix <- with(row, sim_cop_dat(
      n           = 1e3, 
      n_locs      = n_locs, 
      n_clust     = n_cor,
      cor_gauss   = rep(cor_gauss, n_cor),
      cor_t       = cor_t,
      df_t        = rep(df_t, n_cor),
      params_gpd  = c(scale_gpd, shape_gpd),
      mix_p       = c(0.5, 0.5), 
      cluster_mem = ceiling(seq_len(n_locs) / (n_locs / n_cor))
    ))$data_mix
    
    # Fit CE model, if an error is produced, return a dummy list and skip
    dependence <- tryCatch({
      fit_ce(
        data_mix, 
        vars = paste0("col_", seq_len(n_vars)), 
        marg_prob   = marg_prob,
        cond_prob   = row$kl_prob, 
        fit_no_keef = TRUE, 
        output_all = FALSE
      )
    }, error = function(cond) {
      # TODO: Fill in what goes here
      return(1)
    })
    # # gives correct answer!
    # js_clust(
    #   dependence,
    #   k = 3,
    #   cluster_mem =  ceiling(seq_len(n_locs) / (n_locs / n_cor))
    # )

    # Pull JS distance matrix
    dist <- js_clust(dependence, scree_k = 2:15)$dist_mat
    # check for three clusters that distances are distinct
    # plot(dist)
    # summary(as.matrix(dist)[, 1][2:20])
    # summary(as.matrix(dist)[, 1][21:40])
    # summary(as.matrix(dist)[, 1][41:60])

    plot(dist)
    image(as.matrix(dist))
    summary(as.matrix(dist)[, 1][2:15])
    summary(as.matrix(dist)[, 1][16:30])
    summary(as.matrix(dist)[, 1][31:45])
    summary(as.matrix(dist)[, 1][46:60])

    # Calculate TWGSS for elbow
    twgss <- evc::scree_plot(dist, k = 1:15)
    # choose optimal k (from http://sherrytowers.com/2013/10/24/k-means-clustering/)
    opt_k <- \(x) {
      v = -diff(x)
      nv = length(v)
      fom = v[1:(nv-1)]/v[2:nv]
      return(which.max(fom) + 1)
    }
    elb_vec[[j]] <- opt_k(twgss)
    
    # Calculate silhouette coefficient for 2-15 clusters
    sil_boxplot(dist, k = 2:15, show_plt = TRUE)$plot
    sil <- sil_boxplot(dist, k = 2:15, show_plt = FALSE)$sil
    # Find silhouette coefficient for different clusterings
    sil_vec[[j]] <- sil %>% 
      group_by(k) %>% 
      summarise(mean = mean(sil_width, na.rm = TRUE), .groups = "drop") %>% 
      filter(mean == max(mean)) %>% 
      pull(k)
    
    # Calculate AIC for different solutions
    pam_aic <- \(fit) {
      # -2 * fit$objective[[2]] + (2 * length(fit$medoids))
      -2 * fit$objective[[2]] + (4 * length(fit$medoids))
    }
    aic <- vapply(1:15, \(k) {
      pam_fit <- pam(dist, k = k)
      pam_aic(pam_fit)
    }, numeric(1))
    plot(-aic, type = "b")
    aic_vec[[j]] <- opt_k(aic)
    
    # Adjusted AIC 
    # pam_aicc <- \(fit) {
    #   n <- length(fit$clustering)
    #   k <- length(fit$medoids)
    #   -2 * fit$objective[[2]] + (2 * k) + (2 * k * (k + 1)) / (n - k - 1)
    # }
    # aicc <-  vapply(1:15, \(k) {
    #   pam_fit <- pam(dist, k = k)
    #   pam_aicc(pam_fit)
    # }, numeric(1))
    
    # BIC
    # pam_bic <- \(fit) {
    #   n <- length(fit$clustering)
    #   k <- length(fit$medoids)
    #   -2 * fit$objective[[2]] + (k * log(n))
    # }
    # bic <- vapply(1:15, \(k) {
    #   pam_fit <- pam(dist, k = k)
    #   pam_bic(pam_fit)
    # }, numeric(1))
    # plot(-bic, type = "b")
    # opt_k(bic)
  }
  
  return(cbind(row, "sil_k" = sil_vec, "elb_k" = elb_vec, "aic_k" = aic_vec))
}, mc.cores = n_cores))

# remove NAs, give warning
nas <- is.na(results_grid$k)
if (sum(nas) > 0) {
  message(paste0("Warning: ", sum(nas), " NAs in results, removing"))
  results_grid <- results_grid[!nas, ]
}

# Save object
# saveRDS(results_grid, file = paste0("data/js_sil_best_k.RDS"))
# saveRDS(results_grid, file = paste0("data/js_sil_best_k_clust_5.RDS"))
saveRDS(results_grid, file = paste0("data/js_sil_best_k_clust_11.RDS"))
# results_grid <- readRDS(paste0("data/js_sil_best_k_clust_5.RDS"))
# results_grid <- readRDS(paste0("data/js_sil_best_k_clust_5.RDS"))
# results_grid <- readRDS(paste0("data/js_sil_best_k_clust_11.RDS"))

table(results_grid$sil_k)
table(results_grid$elb_k)
table(results_grid$aic_k)



#### Plot ####

# plot occurrence of values of k for different correlation values
# TODO: Add red tp points at 2
p1 <- results_grid %>% 
  group_by(cor_gauss1, cor_t1, cor_t2) %>% 
  dplyr::count(sil_k) %>% 
  mutate(
    col = case_when(
      n == max(n) ~ "Most frequent",
      sil_k == 2  ~  "True value (not found)",
      TRUE        ~  "Else"
    ), 
    col = factor(
      col, levels = c("Most frequent", "True value (not found)", "Else")
    )
  ) %>% 
  ggplot() + 
  geom_point(aes(x = cor_gauss1, y = sil_k, size = n, colour = col)) + 
  scale_y_continuous(
    limits = c(2, 9), 
    breaks = seq(2, 8, by = 2), 
    expand = c(0.15, 0)
  ) +
  scale_size_continuous(range = c(1.5, 6.5)) + 
  facet_grid(cor_t1 ~ cor_t2) + 
  evc_theme() + 
  labs(y = "Adjusted Rand Index", x = "Gaussian corr", colour = "") + 
  scale_colour_manual(values = c(ggsci::pal_nejm()(2), "black")) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
p1

ggsave(
  plot = p1, 
  paste0("plots/04_js_sil_best_k.png"), 
  width = 15, height = 9
  # paste0("plots/01b_js_sensitivity_dqu_", kl_prob, "_marg_0.9.png")
)

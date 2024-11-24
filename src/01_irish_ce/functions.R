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
    cop_norm <- normalCopula(cor, dim = 2, dispstr = "un")
    # cop_norm <- normalCopula(cor_mat, dim = 2, dispstr = "un")
    u <- rCopula(n, cop_norm)
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
    u <- rCopula(n, cop_t)
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

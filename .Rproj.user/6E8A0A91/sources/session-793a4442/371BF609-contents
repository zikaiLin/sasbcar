#' Fit Model using Rcpp
#'
#' @param X A matrix of predictors.
#' @param y A vector of response.
#' @param partition A vector of partition, defaults to NULL.
#' @param radius_candidats A vector of candidate radius parameters, defaults to c(1,2).
#' @param radius_update_thinning An integer specifying thinning interval for radius updates, defaults to 50.
#' @param burnin An integer specifying the number of burn-in iterations, defaults to 2000.
#' @param thinning An integer specifying the thinning interval for MCMC sampling, defaults to 1.
#' @param mcmc_sample An integer specifying the number of MCMC samples to draw, defaults to 2000.
#' @param prior_max A numeric value specifying the maximum prior value of probability map, defaults to 1.
#' @param prior_min A numeric value specifying the minimum prior value of probability map, defaults to 0.
#' @param threshold A numeric value specifying the threshold for inclusion probabilities, defaults to 0.6.
#' @param verbose An integer specifying the verbosity at MCMC sampling steps, defaults to 200 steps.
#' @param kernel_poly_degree An integer specifying the polynomial degree for the squared exponential kernel, defaults to 20.
#' @param kernel_concentration A numeric value specifying the concentration parameter for the squared exponential kernel, defaults to 0.01.
#' @param kernel_smoothness A numeric value specifying the smoothness parameter for the squared exponential kernel, defaults to 50.
#' @param seed An integer specifying the random seed for reproducibility, defaults to 1234.
#' @importFrom BayesGPfit GP.eigen.funcs.fast
#' @importFrom stats rnorm rgamma rbinom pnorm var
#' @export



fit_model <- function(X,
                      y,
                      partition = NULL, 
                      radius_candidats = c(1,2),
                      radius_update_thinning = 50, 
                      burnin = 2000,
                      thinning = 1, 
                      mcmc_sample = 2000, 
                      prior_max = 1,
                      prior_min = 0,
                      threshold = 0.6,
                      verbose = 200,
                      kernel_poly_degree = 20,
                      kernel_concentration = 0.01, 
                      kernel_smoothness = 50, 
                      seed = 1234){

  

  # Get the non-zero columns
  non_zero_column = as.integer(apply(X, 2, function(x) (sum(x)>0)))
  exclude_column = 1- non_zero_column

  # Initialized coefficients from prior
  theta_init = matrix(stats::rnorm((kernel_poly_degree+1)*length(non_zero_column),0,1)
                            ,(kernel_poly_degree+1),length(non_zero_column))

  alpha_init = stats::rnorm(1,0,1)
  z_init = rep(alpha_init,nrow(X))
  delta_init = stats::rbinom(ncol(X),1,0.5) # random initial indicator
  a_init = 1/stats::rgamma(1,0.5,1) # random initial a

  for (j in 1:length(non_zero_column)) {
    temp <- as.matrix(X[,j])
    B <- BayesGPfit::GP.eigen.funcs.fast(temp, 
                                        poly_degree = kernel_poly_degree,
                                        a = kernel_concentration,
                                        b = kernel_smoothness)
    z_init <- z_init+ delta_init[j]* B %*% theta_init[,j]
  }
  
  # Variance partition 
  if(is.null(partition)){
    prior_selection_prob = apply(X,2,stats::var)
    partition = matrix(prior_selection_prob, 28, 28)
    partition[prior_selection_prob>0.05] = 1
    partition[prior_selection_prob<=0.05] = 2
  }

  # Update the partition matrix probability map
  pj_mat = update_pj_mat(matrix(delta_init, 28, 28), 
                         partition,
                         radius_partition = rep(1, length(unique(partition))))
  # Soften the partition probability
  pj_mat[pj_mat == 1] <- 0.95
  pj_mat[pj_mat == 0] <- 0.05
  
  model = Basis_Expansion(
    y = y,
    X = X,
    partition = partition,
    initial_alpha = alpha_init,
    initial_a = a_init,
    initial_delta = delta_init,
    initial_z = z_init ,
    initial_theta = theta_init,
    poly_degree = kernel_poly_degree,
    poly_a = kernel_concentration,
    poly_b = kernel_smoothness,
    oc_theta = TRUE,
    oc_alpha = TRUE,
    oc_delta = TRUE,
    oc_z = TRUE,
    oc_sigma = TRUE,
    oc_a = TRUE,
    burnin = burnin,
    thinning = thinning,
    mcmc_sample = mcmc_sample,
    threshold = 0.6,
    initial_sigma_sq = 1,
    initial_prob_mat = pj_mat,
    radius_candidats = radius_candidats,
    verbose = verbose,
    radius_update_thinning = radius_update_thinning,
    excluded_vox = exclude_column,
    prior_max = prior_max,
    prior_min = prior_min
  )
  return(model)
}
#' Predict using a trained model
#'
#' This function takes a trained model and new data as input and returns predictions.
#'
#' @param model The trained model object.
#' @param new_data The new data for prediction.
#' @param kernel_poly_degree An integer specifying the polynomial degree for the squared exponential kernel, defaults to 20.
#' @param kernel_concentration A numeric value specifying the concentration parameter for the squared exponential kernel, defaults to 0.01.
#' @param kernel_smoothness A numeric value specifying the smoothness parameter for the squared exponential kernel, defaults to 50.
#' @return A vector of predictions.
#' @importFrom BayesGPfit GP.eigen.funcs.fast
#' @importFrom stats rbinom pnorm



predict_model <- function(model, new_data,
                          kernel_poly_degree = 20,
                          kernel_concentration = 0.01, 
                          kernel_smoothness = 50) {
  # Extract posterior mean values from the model
  delta_post = model$post_mean$delta
  theta_post = model$post_mean$theta
  non_zero_column <- non_zeros_voxels(new_data)
  predict_training_p_model_1 = rep(0, nrow(new_data))
  
  for (j in 1:ncol(new_data)) {
    if (non_zero_column[j] == 1) {
      temp <- as.matrix(new_data[, j])
      B <- BayesGPfit::GP.eigen.funcs.fast(temp, 
                                           poly_degree = kernel_poly_degree, 
                                           a = kernel_concentration, 
                                           b = kernel_smoothness)
      predict_training_p_model_1 <-  predict_training_p_model_1 + delta_post[j] * B %*% theta_post[, , j]
    }
  }
  
  # Generate predictions
  predictions <- stats::rbinom(nrow(new_data), 1, stats::pnorm(predict_training_p_model_1))
  
  # Return the predictions
  return(predictions)
}

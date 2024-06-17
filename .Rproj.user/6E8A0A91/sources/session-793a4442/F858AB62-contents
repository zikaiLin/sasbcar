#' This is a utility function that adds two numbers
#' @param training The training data
#' @return A vector of integers indicating the non-zero voxels
non_zeros_voxels <- function(training) {
    return(as.integer(apply(training, 2, function(x) (sum(x)>0))))
}

#' This is a utility function that subtracts two numbers
#' @param training The training data, matrix
#' @param thresholding The thresholding values, prior selection probability greater than thresholding is 1, else 2
#' @return The difference between the two numbers
#' @importFrom stats var


partition_var <- function(training, thresholding=0.05) {
    width = ncol(training)
    height = nrow(training)
    prior_selection_prob = apply(training,2,stats::var)
    partition = matrix(prior_selection_prob, nrow=height, ncol=width)
    # Thresholding
    partition[prior_selection_prob>thresholding] = 1
    partition[prior_selection_prob<=thresholding] = 2
    return(partition)
}

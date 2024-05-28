#' Calculate scale by the population standard deviation, without centering
#'
#' This function allows you to scale vectors of a matrix by their population standard deviation without centering; we assume our sample is the population.
#' @param X numeric matrix
#' @param center Logical: should X be centered? Defaults to FALSE.
#' @keywords internal
scale_varp <- function(X, center = FALSE){

  # calculate centering values (if needed)
  if (center){
    center_vals <- apply(X, 2, function(j) mean(j, na.rm=TRUE))
    X <- apply(X, 2, function(j) j - mean(j, na.rm=TRUE))
  } else {
    center_vals <- rep(0, ncol(X))
  }

  # calculate scaling values
  scale_vals <- apply(X, 2, function(j) sqrt(mean(j^2, na.rm = TRUE)))

  # scale X matrix
  scaled_X <- scale(X, center = center_vals, scale = scale_vals)

  return(list(scale_vals = scale_vals,
              center_vals = center_vals,
              scaled_X = scaled_X))

}

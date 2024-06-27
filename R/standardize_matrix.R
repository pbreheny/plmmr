#' A helper function to standardize matrices 
#'
#' @param X a matrix 
#' @param penalty_factor a vector of 0s and 1s, passed from `plmm()` 
#'
#' @returns a list with the standardized matrix and its details 
#' @keywords internal
#'
standardize_matrix <- function(X, penalty_factor){
  std_X <- ncvreg::std(X)
  std_X_details <- list(center = attr(std_X, 'center'),
                        scale = attr(std_X, 'scale'),
                        ns = attr(std_X, 'nonsingular'))
  if (!is.null(penalty_factor)) {
    # if user supplied penalty factor, make sure to adjust for columns that 
    # may have been constant features 
    penalty_factor <- penalty_factor[std_X_details$ns]
  }
  
  return(list(std_X = std_X,
              std_X_details = std_X_details,
              penalty_factor = penalty_factor))
}

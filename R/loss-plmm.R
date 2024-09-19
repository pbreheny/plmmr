#' Loss method for "plmm" class
#'
#' @param y Observed outcomes (response) vector
#' @param yhat Predicted outcomes (response) vector
#'
#' @returns A numeric vector of the squared-error loss values for the given
#' observed and predicted outcomes
#'
#'
#' @export
#'
#' @examples
#' admix_design <- create_design(X = admix$X, outcome_col = admix$y)
#' fit <- plmm(design = admix_design, K = relatedness_mat(admix$X))
#' yhat <- predict(object = fit, newX = admix$X, type = 'lp', lambda = 0.05)
#' head(plmm_loss(yhat = yhat, y = admix$y))
#'
plmm_loss <- function(y, yhat) {
  val <- (y - yhat)^2
  return(val)
}

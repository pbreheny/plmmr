#' Coef method for "cv.plmm" class
#'
#' @param object An object of class "cv.plmm."
#' @param lambda A numeric vector of lambda values.
#' @param which Vector of lambda indices for which coefficients to return. Defaults to lambda index with minimum CVE.
#' @param ... Additional arguments.
#' @export

coef.cv.plmm <- function(object, lambda, which = object$min, ...){
  coef.plmm(object$fit, lambda = lambda, which = which, ...)
}
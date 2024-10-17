#' Coef method for "cv_plmm" class
#'
#' @param object An object of class "cv_plmm."
#' @param lambda A numeric vector of lambda values.
#' @param which Vector of lambda indices for which coefficients to return. Defaults to lambda index with minimum CVE.
#' @param ... Additional arguments (not used).
#'
#' @rdname coef.cv_plmm
#'
#' @returns Returns a named numeric vector. Values are the coefficients of the
#' model at the specified value of either `lambda` or `which`. Names are the
#' values of `lambda`.
#'
#' @export
#'
#' @examples
#' admix_design <- create_design(X = admix$X, y = admix$y)
#' cv_fit <- cv_plmm(design = admix_design, return_fit = TRUE)
#' head(coef(cv_fit))
#'
coef.cv_plmm <- function(object, lambda, which = object$min, ...){
  coef(object$fit, lambda = lambda, which = which, ...)
}

#' A summary function for `cv_plmm` objects
#'
#' @param object A `cv_plmm` object
#' @param ...  Not used
#'
#' @return The return value is an object with S3 class `summary.cv_plmm`. The class has its own print method and contains the following list elements:
#' * `lambda_min`: The lambda value at the minimum cross validation error
#' * `lambda_1se`: The maximum lambda value within 1 standard error of the minimum cross validation error
#' * `penalty`: The penalty applied to the fitted model
#' * `nvars`: The number of variables selected at lambda_min
#' * `cve`: The cross validation error at all folds
#' * `min`: The minimum cross validation error
#' * `fit`: The `plmm` fit used in the cross validation
#'
#' @rdname summary.cv_plmm
#'
#' @export
#'
#' @examples
#' admix_design <- create_design(X = admix$X, y = admix$y)
#' cv_fit <- cv_plmm(design = admix_design)
#' summary(cv_fit)
summary.cv_plmm <- function(object, ...) {
  structure(
    list(
      lambda_min = object$lambda_min,
      lambda_1se = object$lambda.1se,
      penalty = object$fit$penalty,
      nvars = predict(object$fit, type = "nvars", lambda = object$lambda_min),
      cve = object$cve,
      min = object$min,
      fit = object$fit
    ),
    class = "summary.cv_plmm"
  )
}

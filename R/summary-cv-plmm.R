#' A summary function for cv_plmm objects
#'
#' @param object A \code{cv_plmm} object
#' @param lambda The regularization parameter value at which inference should be reported. Can choose a numeric value, 'min', or '1se'. Defaults to 'min.'
#' @param ...  Not used
#' @return The return value is an object with S3 class `summary.cv_plmm`. The class has its own print method and contains the following list elements:
#' * `lambda_min`: The lambda value at the minimum cross validation error
#' * `lambda.1se`: The maximum lambda value within 1 standard error of the minimum cross validation error
#' * `penalty`: The penalty applied to the fitted model
#' * `nvars`: The number of non-zero coefficients at the selected lambda value
#' * `cve`: The cross validation error at all folds
#' * `min`: The minimum cross validation error
#' * `fit`: The \code{plmm} fit used in the cross validation
#'
#' if `returnBiasDetails = TRUE`, two more items are returned:
#' * `bias`: The mean bias of the cross validation
#' * `loss`: The loss at each value of `lambda`
#'
#'
#' @rdname summary.cv_plmm
#'
#' @export
#'
#' @examples
#' admix_design <- create_design(X = admix$X, outcome_col = admix$y)
#' cv_fit <- cv_plmm(design = admix_design)
#' summary(cv_fit)

summary.cv_plmm <- function(object, lambda = "min", ...){

  # determine the number of nonzero coefficients at specified lambda value

  if(lambda == "min"){
    # nvars (tells number of non-zero coefficients)
    nvars <- predict(object$fit, type="nvars", lambda=object$lambda_min)
  }

  if(lambda == "1se"){
    # nvars (tells number of non-zero coefficients)
    nvars <- predict(object$fit, type="nvars", lambda=object$lambda.1se)
  }

  if(is.numeric(lambda)){
    if(!(lambda %in% object$fit$lambda)) stop("The user-specified lambda is not one of the lambda values used to fit the model.")

    nvars <- predict(object$fit, type="nvars", lambda = lambda)
  }

  # TODO: think about what else should go here
  out <- structure(list(lambda_min = object$lambda_min,
                        lambda.1se = object$lambda.1se,
                        penalty = object$fit$penalty,
                        nvars = nvars,
                        cve = object$cve,
                        min = object$min,
                        fit = object$fit),
                   class = "summary.cv_plmm")
  if("Bias" %in% names(object)){
    out$bias <- object$Bias
    out$loss <- object$Loss
  }

  return(out)

}

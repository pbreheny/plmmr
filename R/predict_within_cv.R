#' Predict method to use in cross-validation (within \code{cvf})
#'
#' @param fit  A list with the components returned by `plmm_fit`.
#' @param testX A design matrix used for computing predicted values (i.e, the test data).
#' @param type A character argument indicating what type of prediction should be returned. Passed from `cvf()`,
#'             Options are "lp," "coefficients," "vars," "nvars," and "blup." See details.
#' @param fbm Logical: is trainX an FBM object? If so, this function expects that testX is also an FBM. The two X matrices must be stored the same way.
#' @param Sigma_11 Variance-covariance matrix of the training data. Extracted from `estimated_Sigma` that is generated using all observations. Required if \code{type == 'blup'}.
#' @param Sigma_21 Covariance matrix between the training and the testing data. Extracted from `estimated_Sigma` that is generated using all observations. Required if \code{type == 'blup'}.
#'
#' @returns A numeric vector of predicted values
#'
#' @details
#' Define beta-hat as the coefficients estimated at the value of lambda that minimizes cross-validation error (CVE). Then options for `type` are as follows:
#'  * 'lp' (default): uses the linear predictor (i.e., product of test data and estimated coefficients) to predict test values of the outcome.
#'    Note that this approach does not incorporate the correlation structure of the data.
#'
#'  * 'blup' (acronym for Best Linear Unbiased Predictor): adds to the 'lp' a value that represents the estimated random effect.
#'  This addition is a way of incorporating the estimated correlation structure of data into our prediction of the outcome.
#'
#' Note: the main difference between this function and the `predict.plmm()` method is that
#' here in CV, the standardized testing data (std_test_X), Sigma_11, and Sigma_21 are calculated in `cvf()` instead of the function defined here.
#'
#' @keywords internal
#'
predict_within_cv <- function(fit,
                              testX,
                              type,
                              fbm = FALSE,
                              Sigma_11 = NULL,
                              Sigma_21 = NULL) {

  # make sure X is in the correct format...
  # case 1: testX is filebacked
  fbm_flag <- inherits(testX,"big.matrix")

  # format dim. names
  if(is.null(dim(fit$beta_vals))) {
    # case 1: fit$beta_vals is a vector
    names(fit$beta_vals) <- lam_names(fit$lambda)
  } else {
    # case 2: fit$beta_vals is a matrix
    colnames(fit$beta_vals) <- lam_names(fit$lambda)
  }

  # calculate the estimated mean values for test data
  a <- fit$beta_vals[1,]
  b <- fit$beta_vals[-1,,drop=FALSE]
  Xb <- sweep(testX %*% b, 2, a, "+")

  # for linear predictor, return mean values
  if (type=="lp") {
    return(drop(Xb))
  }

  # for blup, will incorporate the estimated variance
  if (type == "blup") {
    # TODO: to find the inverse of Sigma_11 using Woodbury's formula? think on this...
    resid_train <- (drop(fit$y) - fit$std_Xbeta)
    ranef <- Sigma_21 %*% chol2inv(chol(Sigma_11)) %*% resid_train
    blup <- Xb + ranef
    return(blup)
  }

}

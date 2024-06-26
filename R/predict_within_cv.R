#' Predict method to use in cross-validation (within \code{cvf})
#'
#' @param fit  A list with the components returned by `plmm_fit`.
#' @param trainX The standardized design matrix of training data, *pre-rotation*.
#' @param testX A design matrix used for computing predicted values (i.e, the test data).
#' @param og_scale_beta testX is on the scale of the original data, so we need the beta_vals that are untransformed to match that scale.
#'       See `plmm_fit()` and `untransform()` for details.
#' @param type A character argument indicating what type of prediction should be returned. Passed from `cvf()`,
#'             Options are "lp," "coefficients," "vars," "nvars," and "blup." See details.
#' @param fbm Logical: is trainX an FBM object? If so, this function expects that testX is also an FBM. The two X matrices must be stored the same way.
#' @param idx Vector of indices of the penalty parameter \code{lambda} at which predictions are required. By default, all indices are returned.
#' @param V11 Variance-covariance matrix of the training data. Extracted from `estimated_V` that is generated using all observations. Required if \code{type == 'blup'}.
#' @param V21 Covariance matrix between the training and the testing data. Extracted from `estimated_V` that is generated using all observations. Required if \code{type == 'blup'}.
#' @param ... Additional optional arguments
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
#' here in CV, predictions are made on the *standardized* scale (i.e., both
#' the trainX and testX data come from std_X). The `predict.plmm()` method
#' makes predictions on the scale of X (the original scale)
#'
#' @keywords internal

predict_within_cv <- function(fit,
                              trainX,
                              testX,
                              og_scale_beta,
                              type,
                              fbm = FALSE,
                              idx=1:length(fit$lambda),
                              V11 = NULL,
                              V21 = NULL, ...) {

  # make sure X is in the correct format...
  # case 1: testX is filebacked
  fbm_flag <- inherits(testX,"big.matrix")

  # get beta values (for nonsingular features) from fit
  std_scale_beta <- fit$std_scale_beta[,idx,drop = FALSE]

  # format dim. names
  if(is.null(dim(std_scale_beta))) {
    # case 1: std_scale_beta is a vector
    names(std_scale_beta) <- lam_names(fit$lambda)
  } else {
    # case 2: std_scale_beta is a matrix
    colnames(std_scale_beta) <- lam_names(fit$lambda)
  }

  # calculate the estimated mean values for test data
  a <- og_scale_beta[1,]
  b <- og_scale_beta[-1,,drop=FALSE]
  Xb <- sweep(testX %*% b, 2, a, "+") # testX is on the original (pre-standardization) scale

  # for linear predictor, return mean values
  if (type=="lp") return(drop(Xb))

  # for blup, will incorporate the estimated variance
  if (type == "blup") {
    # covariance comes from selected rows and columns from estimated_V that is generated in the overall fit (V11, V21)
    # test1 <- V21 %*% chol2inv(chol(V11)) # true
    # TODO: to find the inverse of V11 using Woodbury's formula? think on this...
    a_train <- std_scale_beta[1,]
    b_train <- std_scale_beta[-1,,drop=FALSE]
    Xb_train <- sweep(trainX %*% b_train, 2, a_train, "+")
    resid_train <- (drop(fit$centered_y) - Xb_train)
    ranef <- V21 %*% (chol2inv(chol(V11)) %*% resid_train)
    blup <- Xb + ranef

    return(blup)
  }

}

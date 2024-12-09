#' Predict method to use in cross-validation (within \code{cvf})
#'
#' @param fit  A list with the components returned by `plmm_fit`.
#' @param trainX The training data, *pre-standardization* and *pre-rotation*
#' @param trainY The training outcome, *not centered*. Only needed if `type = 'blup'`
#' @param testX A design matrix used for computing predicted values (i.e, the test data).
#' @param std_X_details A list with 3 items:
#'  * 'center': the centering values for the columns of `X`
#'  * 'scale': the scaling values for the non-singular columns of `X`
#'  * 'ns': indices of nonsingular columns in `std_X`. Note: this is the vector we really need here!
#' @param type A character argument indicating what type of prediction should be returned. Passed from `cvf()`,
#'             Options are "lp," "coefficients," "vars," "nvars," and "blup." See details.
#' @param fbm Logical: is trainX an FBM object? If so, this function expects that testX is also an FBM. The two X matrices must be stored the same way.
#' @param Sigma_11 Variance-covariance matrix of the training data. Extracted from `estimated_Sigma` that is generated using all observations. Required if \code{type == 'blup'}.
#' @param Sigma_21 Covariance matrix between the training and the testing data. Extracted from `estimated_Sigma` that is generated using all observations. Required if \code{type == 'blup'}.
#' @param ... Additional optional arguments
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
#' here in CV, predictions are made on the *standardized* scale (i.e., both
#' the trainX and testX data come from std_X). The `predict.plmm()` method
#' makes predictions on the scale of X (the original scale)
#'
#' @keywords internal
predict_within_cv <- function(fit,
                              trainX,
                              trainY = NULL,
                              testX,
                              std_X_details,
                              type,
                              fbm = FALSE,
                              Sigma_11 = NULL,
                              Sigma_21 = NULL, ...) {

  # make sure X is in the correct format...
  # case 1: testX is filebacked
  fbm_flag <- inherits(testX,"big.matrix")

  train_scale_beta <- fit$std_scale_beta

  # format dim. names
  if(is.null(dim(fit$std_scale_beta))) {
    # case 1: fit$std_scale_beta is a vector
    names(fit$std_scale_beta) <- lam_names(fit$lambda)
  } else {
    # case 2: fit$std_scale_beta is a matrix
    colnames(fit$std_scale_beta) <- lam_names(fit$lambda)
  }

  # adjust the dimension of the estimated coefficients
  train_scale_b_og_dim <- adjust_beta_dimension(std_scale_beta = fit$std_scale_beta,
                                                p = ncol(testX),
                                                std_X_details = std_X_details,
                                                fbm_flag = fbm_flag)
  # calculate the estimated mean values for test data (on the *standardized* scale)
  a <- train_scale_b_og_dim[1,]
  b <- train_scale_b_og_dim[-1,,drop=FALSE]
  Xb <- sweep(testX %*% b, 2, a, "+")

  # for linear predictor, return mean values
  if (type=="lp") {
    return(drop(Xb))
  }

  # for blup, will incorporate the estimated variance
  if (type == "blup") {
    # covariance comes from selected rows and columns from estimated_Sigma that
    #   is generated in the overall fit (Sigma_11, Sigma_21)
    # TODO: to find the inverse of Sigma_11 using Woodbury's formula? think on this...
    aa <- fit$std_scale_beta[1,]
    bb <- fit$std_scale_beta[-1,]
    Xb_train <- sweep(trainX %*% bb, 2, aa, "+")
    resid_train <- (drop(trainY) - Xb_train)
    ranef <- Sigma_21 %*% chol2inv(chol(Sigma_11)) %*% resid_train
    blup <- Xb + ranef
    return(blup)
  }

}

#' Predict method to use in cross-validation (within `cvf()`)
#'
#' @param fit  A list with the components returned by `plmm_fit`.
#' @param testX A design matrix used for computing predicted values (i.e, the test data).
#' @param type A character argument indicating what type of prediction should be returned. Passed from `cvf()`.
#'             Options are "lp," "coefficients," "vars," "nvars," and "blup." See details.
#' @param fbm Logical: is `trainX` a filebacked `big.matrix` object? If so, this function expects that `testX` is also an FBM. The two X matrices must be stored the same way.
#' @param Sigma_21 Covariance matrix between the training and the testing data. Required if `type == 'blup'`.
#'
#' @return A numeric vector of predicted values
#'
#' @details
#'
#'  * `lp` (linear predictor): uses the product of `testX` and the beta coefficients of `fit` to predict new values of the outcome. This does not incorporate the correlation structure of the data.
#'
#'  * `blup` (acronym for Best Linear Unbiased Predictor): adds to the `lp`` a value that represents the estimated random effect. This addition is a way of incorporating
#'     the estimated correlation structure of data into our prediction of the outcome.
#'
#'  * `coefficients`: returns the estimated beta-hat
#'
#'  * `vars`: returns the _indices_ of variables (e.g., SNPs) with nonzero coefficients at each value of lambda. EXCLUDES intercept.
#'
#'  * `nvars`: returns the _number_ of variables (e.g., SNPs) with nonzero coefficients at each value of lambda. EXCLUDES intercept.
#'
#' Note: the main difference between this function and the `predict.plmm()` method is that
#' here in CV, the standardized testing data (`std_test_X`), `Sigma_11`, and `Sigma_21` are calculated in `cvf()` instead of the function defined here.
#'
#' @keywords internal
#'
predict_within_cv <- function(fit,
                              testX,
                              type,
                              fbm = FALSE,
                              Sigma_21 = NULL) {
  # format dim. names
  if (is.null(dim(fit$beta_vals))) {
    # case 1: fit$beta_vals is a vector
    names(fit$beta_vals) <- lam_names(fit$lambda)
  } else {
    # case 2: fit$beta_vals is a matrix
    colnames(fit$beta_vals) <- lam_names(fit$lambda)
  }

  # calculate the estimated mean values for test data
  a <- fit$beta_vals[1, ]
  b <- fit$beta_vals[-1, , drop = FALSE]
  Xb <- sweep(testX %*% b, 2, a, "+")

  # for linear predictor, return mean values
  if (type == "lp") {
    Xb
  } else if (type == "blup") {
    # for blup, will incorporate the estimated variance
    compute_blup(fit, Xb, Sigma_21, idx = seq_along(fit$lambda))
  }

}

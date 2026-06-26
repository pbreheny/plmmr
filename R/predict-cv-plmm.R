#' Predict method for `cv_plmm` class
#'
#' @param object    An object of class `cv_plmm`.
#' @param newX      Matrix of values at which predictions are to be made (not used for `type` = "coefficients", "vars", or "nvars").
#'                  This can be either a filebacked `big.matrix` or a matrix object. Note: Columns of this argument must be named!
#' @param type      A character argument indicating what type of prediction should be
#'                  returned. Options are "lp," "coefficients," "vars," "nvars," and "blup." See details.
#' @param X         Optional: if `type = 'blup'` and the model was fit in-memory, the design matrix used to fit the model represented in `object` must be supplied.
#'                  When supplied, this design matrix will be standardized using the center/scale values in `object$std_X_details`, so please **do not** standardize this matrix before supplying here.
#'                  **Note**: If the model was fit file-backed, then the filepath to the `.bk` file with this standardized design matrix is returned as `std_X` in the fit supplied to `object`.
#' @param lambda    A numeric vector of regularization parameter `lambda` values
#'                  at which predictions are requested.
#' @param idx       Vector of indices of regularization parameter `lambda` at which
#'                  predictions are requested.  By default, this is the lambda index which minimizes the cross-validation error.
#' @param ...       Additional optional arguments
#'
#' @return Depends on the `type` - see Details
#'
#' @details
#' Define beta-hat as the coefficients estimated at the value of lambda that minimizes cross-validation error (CVE). Then options for `type` are as follows:
#'
#'  * `lp` (linear predictor): uses the product of `newX` and the beta coefficients of `object` to predict new values of the outcome. This does not incorporate the correlation structure of the data.
#'
#'  * `blup` (acronym for Best Linear Unbiased Predictor): adds to the `lp` a value that represents the estimated random effect. This addition is a way of incorporating
#'  the estimated correlation structure of data into our prediction of the outcome.
#'
#'  * `coefficients`: returns the estimated beta-hat
#'
#'  * `vars`: returns the _indices_ of variables (e.g., SNPs) with nonzero coefficients at each value of lambda. EXCLUDES intercept.
#'
#'  * `nvars`: returns the _number_ of variables (e.g., SNPs) with nonzero coefficients at each value of lambda. EXCLUDES intercept.
#'
#' @rdname predict.cv_plmm
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' train_idx <- sample(1:nrow(admix$X), 100)
#' # Note: ^ shuffling is important here! Keeps test and train groups comparable.
#' train <- list(X = admix$X[train_idx,], y = admix$y[train_idx])
#' train_design <- create_design(X = train$X, y = train$y)
#'
#' test <- list(X = admix$X[-train_idx,], y = admix$y[-train_idx])
#' fit <- cv_plmm(design = train_design)
#'
#' pred1 <- predict(object = fit, newX = test$X, X = train$X) # Minimum CVE lambda
#' pred2 <- predict(object = fit, newX = test$X, X = train$X, idx = fit$min1se) # 1 SE lambda
predict.cv_plmm <- function(
  object,
  newX,
  type = c("blup", "coefficients", "vars", "nvars", "lp"),
  X,
  lambda,
  idx = object$min,
  ...
) {
  predict(object$fit, newX = newX, type = type, X = X, lambda = lambda, idx = idx, ...)
}

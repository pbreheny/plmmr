#' Predict method for plmm class
#'
#' @param object    An object of class \code{plmm}.
#' @param newX      Matrix of values at which predictions are to be made (not used for
#'                  `type="coefficients"` or for some of the `type` settings in `predict`).
#'                  This can be either a FBM object or a 'matrix' object.
#'                  **Note**: Columns of this argument must be named!
#' @param type      A character argument indicating what type of prediction should be
#'                  returned. Options are "lp," "coefficients," "vars," "nvars," and "blup." See details.
#' @param X         Optional: if \code{type = 'blup'} and the model was fit in-memory, the design matrix used to fit the model represented in \code{object} must be supplied.
#'                  This design matrix will be standardized using the center/scale values in \code{object$std_X_details}, so please **do not** standardize this matrix before supplying here.
#'                  **Note**: If the model was fit file-backed, then the filepath to the .bk file with this matrix is returned as 'std_X' in the fit supplied to 'object'.
#' @param lambda    A numeric vector of regularization parameter \code{lambda} values
#'                  at which predictions are requested.
#' @param idx       Vector of indices of the penalty parameter \code{lambda} at which
#'                  predictions are required. By default, all indices are returned.
#' @param ...       Additional optional arguments
#'
#' @details
#' Define beta-hat as the coefficients estimated at the value of lambda that minimizes cross-validation error (CVE). Then options for `type` are as follows:
#'
#'  * 'response' (default): uses the product of newX and beta-hat to predict new values of the outcome. This does not incorporate the correlation structure of the data.
#'  For the stats folks out there, this is simply the linear predictor.
#'
#'  * 'blup' (acronym for Best Linear Unbiased Predictor): adds to the 'response' a value that represents the esetimated random effect. This addition is a way of incorporating
#'  the estimated correlation structure of data into our prediction of the outcome.
#'
#'  * 'coefficients': returns the estimated beta-hat
#'
#'  * 'vars': returns the _indices_ of variables (e.g., SNPs) with nonzero coefficients at each value of lambda. EXCLUDES intercept.
#'
#'  * 'nvars': returns the _number_ of variables (e.g., SNPs) with nonzero coefficients at each value of lambda. EXCLUDES intercept.
#'
#'
#' @rdname predict.plmm
#'
#' @returns Depends on the `type` - see Details
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
#' fit <- plmm(design = train_design)
#'
#' # make predictions for all lambda values
#'  pred1 <- predict(object = fit, newX = test$X, type = "lp")
#'  pred2 <- predict(object = fit, newX = test$X, type = "blup", X = train$X)
#'
#' # look at mean squared prediction error
#' mspe <- apply(pred1, 2, function(c){crossprod(test$y - c)/length(c)})
#' min(mspe)
#'
#' mspe_blup <- apply(pred2, 2, function(c){crossprod(test$y - c)/length(c)})
#' min(mspe_blup) # BLUP is better
#'
#' # compare the MSPE of our model to a null model, for reference
#' # null model = intercept only -> y_hat is always mean(y)
#' crossprod(mean(test$y) - test$y)/length(test$y)
#'
#'
predict.plmm <- function(object,
                         newX,
                         type=c("blup", "coefficients", "vars", "nvars", "lp"),
                         X = NULL,
                         lambda,
                         idx=1:length(object$lambda),
                         ...) {

  # if predictions are to be made, make sure X is in the correct format...
  if (!missing(newX)){
    # case 1: newX is a big.matrix
    if (inherits(newX, "big.matrix")){
      fbm_flag <- TRUE
    } else {
      # case 2: X is in-memory
      fbm_flag <- FALSE
    }

  }
  # check format for beta values
  ifelse(inherits(object$beta_vals, "Matrix"),
         sparse_flag <- TRUE,
         sparse_flag <- FALSE)

  # prepare other arguments
  type <- match.arg(type)
  beta_vals <- coef(object, lambda, which=idx, drop=FALSE) # includes intercept
  p <- nrow(beta_vals)-1
  n <- nrow(object$linear_predictors)

  # addressing each type:

  if (type=="coefficients") return(beta_vals)

  if (type=="nvars") return(apply(beta_vals[-1, , drop=FALSE]!=0, 2, sum)) # don't count intercept

  if (type=="vars") return(drop(apply(beta_vals[-1, , drop=FALSE]!=0, 2, FUN=which))) # don't count intercept

  # calculate linear predictors with new data
  a <- beta_vals[1,]
  b <- beta_vals[-1,,drop=FALSE]
  Xb <- sweep(newX %*% b, 2, a, "+")

  if (type=="lp") {
    return(drop(Xb))
  }

  # Note: the BLUP is constructed on the scale of the standardized model-fitting data;
  #   This is the scale at which object$K was constructed
  if (type == "blup"){
    # check dimensions -- must have same number of features in test & train data
    if (object$p != ncol(newX)){stop("the X from the model fit and newX do not have the same number of features - please make these align\n")}

   # check for singularity -- this keeps us from scaling by a 0 value
    singular <- setdiff(seq(1:length(object$std_X_details$center)),
                        object$std_X_details$ns)
    if (length(singular) >= 1) object$std_X_details$scale[singular] <- 1

    # use center/scale values from the X in the model fit to standardize both X and newX
    std_X <- scale(X,
                   center = object$std_X_details$center,
                   scale = object$std_X_details$scale)
    if (fbm_flag) {
      std_test_info <- .Call("big_std",
                             newX@address,
                             as.integer(count_cores()),
                             object$std_X_details$center,
                             object$std_X_details$scale,
                             PACKAGE = "plmmr")
      newX@address <- std_test_info$std_X # now, newX is standardized
    } else {
      std_newX <- scale(newX,
                        center = object$std_X_details$center,
                        scale = object$std_X_details$scale)
    }

    if (fbm_flag) {
      Sigma_11 <- construct_variance(fit = object)
      const <- (object$eta/ncol(newX))
      XXt <- bigalgebra::dgemm(TRANSA = 'N',
                               TRANSB = 'T',
                               A = newX,
                               B = object$std_X)
      Sigma_21 <- const*XXt
      Sigma_21 <- Sigma_21[,] # convert to in-memory matrix
      # TODO: need to address the case where newX has columns that were
      # removed from object$std_X
    } else {
      Sigma_11 <- construct_variance(fit = object)
      Sigma_21 <- object$eta * (1/p)*tcrossprod(std_newX, std_X)
    }
    resid_old <- drop(object$y) - object$std_Xbeta
    ranef <- Sigma_21 %*% (chol2inv(chol(Sigma_11)) %*% resid_old)
    blup <- drop(Xb + ranef)

    return(blup)
  }

}

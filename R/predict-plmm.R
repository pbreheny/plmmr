#' Predict method for plmm class
#'
#' @param object    An object of class \code{plmm}.
#' @param newX      Matrix of values at which predictions are to be made (not used for
#'                  `type="coefficients"` or for some of the `type` settings in `predict`).
#'                  This can be either a FBM object or a 'matrix' object.
#'                  **Note**: Columns of this argument must be named!
#' @param type      A character argument indicating what type of prediction should be
#'                  returned. Options are "lp," "coefficients," "vars," "nvars," and "blup." See details.
#' @param lambda    A numeric vector of regularization parameter \code{lambda} values
#'                  at which predictions are requested.
#' @param idx       Vector of indices of the penalty parameter \code{lambda} at which
#'                  predictions are required. By default, all indices are returned.
#' @param X         Original design matrix (not including intercept column)
#'                  from object. Required if \code{type == 'blup'} and object is too large to be
#'                  returned in `plmm` object. Again, **columns must be named!**
#' @param y         Original continuous outcome vector from object.
#'                  Required if \code{type == 'blup'}.
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
#'  * 'vars': returns the _indicies_ of variables (e.g., SNPs) with nonzero coefficients at each value of lambda. EXCLUDES intercept.
#'
#'  * 'nvars': returns the _number_ of variables (e.g., SNPs) with nonzero coefficients at each value of lambda. EXCLUDES intercept.
#'
#'
#' @rdname predict.plmm
#'
#' @returns Depends on the `type` - see Details
#'
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' train_idx <- sample(1:nrow(admix$X), 100)
#' # Note: ^ shuffling is important here! Keeps test and train groups comparable.
#' train <- list(X = admix$X[train_idx,], y = admix$y[train_idx])
#' test <- list(X = admix$X[-train_idx,], y = admix$y[-train_idx])
#' fit <- plmm(X = train$X, y = train$y)
#'
#' # make predictions for all lambda values
#'  pred1 <- predict(object = fit, newX = test$X, type = "blup", X = train$X, y = train$y)
#'
#' # look at mean squared prediction error
#' mspe <- apply(pred1, 2, function(c){crossprod(test$y - c)/length(c)})
#' min(mspe)
#'
#' # compare the MSPE of our model to a null model, for reference
#' # null model = intercept only -> y_hat is always mean(y)
#' crossprod(mean(test$y) - test$y)/length(test$y)
#'
#'
predict.plmm <- function(object,
                         newX,
                         type=c("lp", "coefficients", "vars", "nvars", "blup"),
                         lambda,
                         idx=1:length(object$lambda),
                         X,
                         y,
                         ...) {

  # object type checks
  if (!missing(X)){
    if (!identical(class(X), class(newX))) {
      stop("\nFor now, the classes of X and newX must match (we plan to extend/enhance this
           further in the future). bigstatsr::as_FBM() and/or bigmemory::as.big.matrix()
           if you need to convert the type of one of your matrices.")
    }
  }


  # if predictions are to be made, make sure X is in the correct format...
  if (!missing(newX)){
    # case 1: newX is an FBM
    if (inherits(newX,"FBM")){
      fbm_flag <- TRUE
      # convert to big.matrix (FBM cannot multiply with dgCMatrix type of beta_vals)
      newX <- fbm2bm(newX)
    } else if (inherits(newX, "big.matrix")){
      # case 2: newX is a big.matrix
      fbm_flag <- TRUE
    } else {
      # case 3: X is in-memory
      fbm_flag <- FALSE
    }

  }

  # check format for beta_vals
  ifelse(inherits(object$beta_vals, "Matrix"),
         sparse_flag <- TRUE,
         sparse_flag <- FALSE)

  # prepare other arguments
  type <- match.arg(type)
  beta_vals <- coef(object, lambda, which=idx, drop=FALSE) # includes intercept
  p <- object$p
  n <- object$n

  # addressing each type:

  if (type=="coefficients") return(beta_vals)

  if (type=="nvars") return(apply(beta_vals[-1, , drop=FALSE]!=0, 2, sum)) # don't count intercept

  if (type=="vars") return(drop(apply(beta_vals[-1, , drop=FALSE]!=0, 2, FUN=which))) # don't count intercept

    a <- beta_vals[1,]
    b <- beta_vals[-1,,drop=FALSE]
    Xb <- sweep(newX %*% b, 2, a, "+")

  if (type=="lp") return(drop(Xb))

  if (type == "blup"){ # assuming eta of X and newX are the same
    if (fbm_flag) stop("\nBLUP prediction outside of cross-validation is not yet implemented for filebacked data. This will be available soon.")
    if (missing(X)) stop("The design matrix is required for BLUP calculation. Please supply the no-intercept design matrix to the X argument.")
    if (missing(y) & is.null(object$y)) stop("The vector of outcomes is required for BLUP calculation. Please either supply it to the y argument, or set returnX=TRUE in the plmm function.")
    # check dimensions -- must have same number of features
    if (ncol(X) != ncol(newX)){stop("\nX and newX do not have the same number of features - please make these align")}
    if (!is.null(object$y)) y <- object$y

    V11 <- construct_variance(fit = object)
    V21 <- object$eta * (1/p) * tcrossprod(newX,X) # same as V21_check
    Xb_old <- sweep(X %*% b, 2, a, "+")
    resid_old <- drop(y) - Xb_old

    ranef <- V21 %*% (chol2inv(chol(V11)) %*% resid_old)
    blup <- drop(Xb + ranef)

    return(blup)
  }

}

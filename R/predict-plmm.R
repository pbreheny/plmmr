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
#'  pred2 <- predict(object = fit, newX = test$X, type = "blup")
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
                         type=c("lp", "coefficients", "vars", "nvars", "blup"),
                         lambda,
                         idx=1:length(object$lambda),
                         ...) {

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

  if (type=="lp") {
    a <- beta_vals[1,]
    b <- beta_vals[-1,,drop=FALSE]
    Xb <- sweep(newX %*% b, 2, a, "+")
    return(drop(Xb))
  }

  if (type == "blup"){
    if (fbm_flag) stop("\nBLUP prediction outside of cross-validation is not yet implemented for filebacked data. This will be available soon.")
    # check dimensions -- must have same number of features in test & train data
    if (length(object$std_X_details$center) != ncol(newX)){stop("\nX and newX do not have the same number of features - please make these align")}

    # below, we subset the columns to include only the nonsingular features from
    #   the training data -- we don't have estimated beta coefs. for these features!
    singular <- setdiff(seq(1:length(object$std_X_details$center)),
                        object$std_X_details$ns)
    object$std_X_details$scale[singular] <- 1
    # object$std_X_details$center[singular] <- 0
    std_newX <- scale(newX,
                      center = object$std_X_details$center,
                      scale = object$std_X_details$scale)

    train_scale_beta_og_dim <- adjust_beta_dimension(object$std_scale_beta,
                                                     p = ncol(newX),
                                                     std_X_details = object$std_X_details,
                                                     fbm_flag = fbm_flag)
    a <- train_scale_beta_og_dim[1,]
    b <- train_scale_beta_og_dim[-1,,drop=FALSE]
    Xb <- sweep(std_newX %*% b, 2, a, "+")

    Sigma_11 <- construct_variance(fit = object)
    Sigma_21 <- object$eta * (1/p)*tcrossprod(std_newX[,object$std_X_details$ns],
                                              object$std_X)
    Xb_old <- sweep(object$std_X %*% object$std_scale_beta[-1,], 2,
                    object$std_scale_beta[1,], "+")
    resid_old <- drop(object$y) - Xb_old

    ranef <- Sigma_21 %*% (chol2inv(chol(Sigma_11)) %*% resid_old)
    blup <- drop(Xb + ranef)

    return(blup)
  }

}

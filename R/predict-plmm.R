#' Predict method for plmm class
#'
#' @param object An object of class \code{plmm}.
#' @param newX Design matrix used for computing predicted values if requested. This 
#' can be either a FBM object or a matrix object. **If you supply an FBM object here, 
#' this function assumes that this matrix has been done.** If you need to do standardization,
#' see `big_std()` and `process_plink()`.
#'  **Columns must be named!**
#' @param type A character argument indicating what type of prediction should be 
#' returned. Options are "lp," "coefficients," "vars," "nvars," and "blup." See details. 
#' @param lambda A numeric vector of regularization parameter \code{lambda} values 
#' at which predictions are requested.
#' @param idx Vector of indices of the penalty parameter \code{lambda} at which 
#' predictions are required. By default, all indices are returned.
#' @param X Original design matrix (not including intercept column) 
#' from object. Required if \code{type == 'blup'} and object is too large to be 
#' returned in `plmm` object. Again, **columns must be named!**
#' @param y Original continuous outcome vector from object. 
#' Required if \code{type == 'blup'}. 
#' @param K An optional list or matrix as returned by `choose_K()`. 
#' @param ... Additional optional arguments
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
#' @export
#'
#' @examples 
#' set.seed(123)
#' train_idx <- sample(1:nrow(admix$X), 100) 
#' # Note: ^ shuffling is important here! Keeps test and train groups comparable.
#' train <- list(X = admix$X[train_idx,], y = admix$y[train_idx])
#' test <- list(X = admix$X[-train_idx,], y = admix$y[-train_idx])
#' fit <- plmm(X = train$X, y = train$y, K = relatedness_mat(train$X))
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
#'  \dontrun{
#'   # file-backed example 
#'   plmm(X = "~/tmp_files/penncath_lite", penalty = "lasso", trace = T, returnX = FALSE) -> foo
#'   pen <- bigsnpr::snp_attach("~/tmp_files/penncath_lite.rds")
#'   y_hat <- predict(foo, newX = pen$genotypes)
#'   y <- pen$fam$affection
#'   # notice: many of the rows in pen$genotypes correspond to samples with a missing phenotype
#'   # to assess the quality of the fit 
#'   names(y) <- pen$fam$family.ID
#'   y_idx <- which(as.character(pen$std_X_rownames) %in% names(y))
#'   crossprod(y - y_hat[y_idx,50])/length(y) # estimate of MSPE
#'  }
#'  
#'  
predict.plmm <- function(object,
                         newX,
                         type=c("lp", "coefficients", "vars", "nvars", "blup"),
                         lambda,
                         idx=1:length(object$lambda),
                         X,
                         y,
                         K = NULL,
                         ...) {
  
  # object type checks
  if (!missing(X)){
    if (!identical(class(X), class(newX))) {
      stop("\nFor now, the classes of X and newX must match (we plan to extend/enhance this 
           further in the future). See plmmr::process_plink() and bigstatsr::as_FBM()
           if you need to convert the type of one of your matrices.")
    }
  }
  
  if (!missing(newX)){
    ifelse(grepl("FBM", class(newX)),
           fbm_flag <- TRUE,
           fbm_flag <- FALSE)
    
  }
  
  ifelse(grepl("Matrix", class(object$beta_vals)),
         sparse_flag <- TRUE,
         sparse_flag <- FALSE)
  
  
  # prepare other arguments 
  type <- match.arg(type)
  beta_vals <- coef.plmm(object, lambda, which=idx, drop=FALSE) # includes intercept 
  p <- object$p 
  n <- object$n 

  # addressing each type: 
  
  if (type=="coefficients") return(beta_vals)

  if (type=="nvars") return(apply(beta_vals[-1, , drop=FALSE]!=0, 2, sum)) # don't count intercept
  
  if (type=="vars") return(drop(apply(beta_vals[-1, , drop=FALSE]!=0, 2, FUN=which))) # don't count intercept
  
  if (fbm_flag){
    # add column of 1s for intercept 
    newX_with_intcpt <- bigstatsr::FBM(init = 1,
                                nrow = newX$nrow,
                                ncol = newX$ncol + 1) 
    # fill in other columns with values of newX
    bigstatsr::big_apply(newX,
                         a.FUN = function(X, ind, res){
                           res[,ind+1] <- X[,ind]
                         },
                         a.combine = cbind,
                         res = newX_with_intcpt)
    # convert to big.matrix (FBM cannot multiply with dgCMatrix type of beta_vals)
    bm_newX <- fbm2bm(newX_with_intcpt)
    # calculate linear predictor 
    Xb <- bm_newX %*% beta_vals
  } else {
    a <- beta_vals[1,]
    b <- beta_vals[-1,,drop=FALSE]
    Xb <- sweep(newX %*% b, 2, a, "+")
  }
  
  if (type=="lp") return(drop(Xb))
  
  if (type == "blup"){ # assuming eta of X and newX are the same 
    if (fbm_flag) stop("\nBLUP prediction outside of cross-validation is not yet implemented for filebacked data. This will be available soon.")
    if (missing(X)) stop("The design matrix is required for BLUP calculation. Please supply the no-intercept design matrix to the X argument.") 
    if (missing(y) & is.null(object$y)) stop("The vector of outcomes is required for BLUP calculation. Please either supply it to the y argument, or set returnX=TRUE in the plmm function.")
    
    if (!is.null(object$y)) y <- object$y
   #  ns <- ifelse(object$ns_idx)
    if (is.null(K)){
      # TODO: consider whether using the standardized scale for prediction would be better....
    # standardize (this won't affect predicted y)
    # std_X <- ncvreg::std(X)
    # std_newX <- ncvreg::std(newX)
    # std_X_details <- list(ns = attr(std_X, 'nonsingular'),
    #                       scale = attr(std_X, 'scale'),
    #                       center = attr(std_X, 'center'))
      V11 <- object$eta*relatedness_mat(X, fbm = fbm_flag, std = FALSE) + (1 - object$eta)*diag(object$n)
    } else {
      V11 <- object$eta*K + (1 - object$eta)*diag(object$n)
    }
    
    # check dimensions -- must have same number of features
    # even from the same dataset, this can be an issue with features 'becoming' constant...
    p1 <- ncol(X)+1
    p2 <- ncol(newX)+1
    p3 <- nrow(beta_vals)
    if (!identical(p1, p2) | !identical(p1, p3)) {
      shared_p <- intersect(colnames(X), colnames(newX)) # this is why we need column names
      X <- X[,shared_p]
      newX <- newX[,shared_p]
      beta_vals <- beta_vals[c("(Intercept)",shared_p),]
    }
    
    # cannot use U, S when computing V21, because nv=0. V is needed to restore X 
    V21 <- object$eta * (1/p) * tcrossprod(newX,X) # same as V21_check 
    Xb_old <- sweep(X %*% b, 2, a, "+")
    resid_old <- drop(y) - Xb_old
    
    ranef <- V21 %*% (chol2inv(chol(V11)) %*% resid_old)
    blup <- drop(Xb + ranef)

    return(blup)
  }

}

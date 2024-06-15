#' PLMM fit: a function that fits a PLMM using the values returned by plmm_prep()
#'
#' @param prep A list as returned from \code{plmm_prep}
#' @param std_X_details A list with components 'center' (values used to center X), 'scale' (values used to scale X), and 'ns' (indices for nonsignular columns of X)
#' @param eta_star The ratio of variances (passed from plmm())
#' @param penalty.factor A multiplicative factor for the penalty applied to each coefficient. If supplied, penalty.factor must be a numeric vector of length equal to the number of columns of X. The purpose of penalty.factor is to apply differential penalization if some coefficients are thought to be more likely than others to be in the model. In particular, penalty.factor can be 0, in which case the coefficient is always in the model without shrinkage.
#' @param fbm_flag Logical: is std_X an FBM object? Passed from `plmm()`.
#' @param penalty The penalty to be applied to the model. Either "MCP" (the default), "SCAD", or "lasso".
#' @param gamma The tuning parameter of the MCP/SCAD penalty (see details). Default is 3 for MCP and 3.7 for SCAD.
#' @param alpha Tuning parameter for the Mnet estimator which controls the relative contributions from the MCP/SCAD penalty and the ridge, or L2 penalty. alpha=1 is equivalent to MCP/SCAD penalty, while alpha=0 would be equivalent to ridge regression. However, alpha=0 is not supported; alpha may be arbitrarily small, but not exactly 0.
#' @param lambda.min The smallest value for lambda, as a fraction of lambda.max. Default is .001 if the number of observations is larger than the number of covariates and .05 otherwise.
#' @param nlambda Length of the sequence of lambda. Default is 100.
#' @param lambda A user-specified sequence of lambda values. By default, a sequence of values of length nlambda is computed, equally spaced on the log scale.
#' @param eps Convergence threshold. The algorithm iterates until the RMSD for the change in linear predictors for each coefficient is less than eps. Default is \code{1e-4}.
#' @param max.iter Maximum number of iterations (total across entire path). Default is 10000.
#' @param convex (future idea; not yet incorporated) convex Calculate index for which objective function ceases to be locally convex? Default is TRUE.
#' @param dfmax (future idea; not yet incorporated) Upper bound for the number of nonzero coefficients. Default is no upper bound. However, for large data sets, computational burden may be heavy for models with a large number of nonzero coefficients.
#' @param init Initial values for coefficients. Default is 0 for all columns of X.
#' @param warn Return warning messages for failures to converge and model saturation? Default is TRUE.
#' @param returnX Return the standardized design matrix along with the fit? By default, this option is turned on if X is under 100 MB, but turned off for larger matrices to preserve memory.
#' @param ... Additional arguments that can be passed to `biglasso::biglasso_simple_path()`
#' @returns  A list with these components:
#'   * n: # of rows in X
#'   * p: # of columns in X (including constant features)
#'   * y: outcome on original scale
#'   * std_X_details: list with 'center' and 'scale' vectors, same as `plmm_prep()`
#'   * s: eigenvalues of K
#'   * U: eigenvectors of K
#'   * rot_X: X on the rotated (i.e., transformed) scale.
#'    Note that the dimensions of `rot_X` are likely to be different than those of X.
#'   * rot_y: y on the rotated scale
#'   * stdrot_X: X on the rotated scale once it has been re-standardized.
#'   * lambda: vector of tuning parameter values
#'   * stdrot_scale_beta: the coefficients estimated on the scale of `stdrot_X`
#'   * std_scale_beta: the coefficients estimated on the scale of `std_X`
#'   * linear.predictors: the product of `stdrot_X` and `b`
#'    (linear predictors on the transformed and restandardized scale)
#'   * eta: a number (double) between 0 and 1 representing the estimated
#'    proportion of the variance in the outcome attributable to population/correlation
#'    structure.
#'   * iter: numeric vector with the number of iterations needed in model fitting
#'    for each value of `lambda`
#'   * converged: vector of logical values indicating whether the model fitting
#'    converged at each value of `lambda`
#'   * loss: vector with the numeric values of the loss at each value of `lambda`
#'    (calculated on the ~rotated~ scale)
#'   * penalty: character string indicating the penalty with which the model was
#'    fit (e.g., 'MCP')
#'   * penalty.factor: vector of indicators corresponding to each predictor,
#'    where 1 = predictor was penalized.
#'   * gamma: numeric value indicating the tuning parameter used for the SCAD or
#'    lasso penalties was used. Not relevant for lasso models.
#'   * alpha: numeric value indicating the elastic net tuning parameter.
#'   * ns: the indices for the nonsingular values of X
#'   * feature_names: formatted column names of the design matrix
#'   * nlambda: number of lambda values used in model fitting
#'   * eps: tolerance ('epsilon') used for model fitting
#'   * max.iter: max. number of iterations per model fit
#'   * warn: logical - should warnings be given if model fit does not converge?
#'   * init: initial values for model fitting
#'   * trace: logical - should messages be printed to the console while models are fit?
#'
#' @keywords internal
#'

plmm_fit <- function(prep,
                     std_X_details,
                     eta_star,
                     penalty.factor,
                     fbm_flag,
                     penalty,
                     gamma = 3,
                     alpha = 1,
                     lambda.min,
                     nlambda = 100,
                     lambda,
                     eps = 1e-04,
                     max.iter = 10000,
                     convex = TRUE,
                     dfmax = prep$p + 1,
                     init = NULL,
                     warn = TRUE,
                     returnX = TRUE,
                     ...){

  # error checking ------------------------------------------------------------
  if (penalty=="MCP" && gamma <= 1) stop("gamma must be greater than 1 for the MC penalty", call.=FALSE)
  if (penalty=="SCAD" && gamma <= 2) stop("gamma must be greater than 2 for the SCAD penalty", call.=FALSE)
  if (nlambda < 2) stop("nlambda must be at least 2", call.=FALSE)
  if (alpha <= 0) stop("alpha must be greater than 0; choose a small positive number instead", call.=FALSE)

  if(prep$trace){cat("\nBeginning rotation ('preconditioning').")}

  # rotate data ----------------------------------------------------------------
  if(!fbm_flag) {
    w <- (prep$eta * prep$s + (1 - prep$eta))^(-1/2)
    wUt <- sweep(x = t(prep$U), MARGIN = 1, STATS = w, FUN = "*")
    rot_X <- wUt %*% prep$std_X
    rot_y <- wUt %*% prep$y

    # re-standardize rot_X
    stdrot_X <- ncvreg::std(rot_X)
    stdrot_X_details <- list(center = attr(stdrot_X, "center"),
                            scale = attr(stdrot_X, "scale"))
  } else {
    rot_res <- rotate_filebacked(prep)
    stdrot_X <- rot_res$stdrot_X
    rot_y <- rot_res$rot_y
    stdrot_X_details <- list(center = rot_res$stdrot_X_center,
                             scale = rot_res$stdrot_X_scale)
  }


  if (prep$trace)(cat("\nRotation (preconditiong) finished at ",
                      format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

  # set up lambda -------------------------------------------------------
  if (prep$trace) cat("\nSetting up lambda/preparing for model fitting.")
  if (missing(lambda)) {
    lambda <- setup_lambda(X = stdrot_X,
                           y = rot_y,
                           alpha = alpha,
                           nlambda = nlambda,
                           lambda.min = lambda.min,
                           penalty.factor = penalty.factor)
    user.lambda <- FALSE
  } else {
    # make sure (if user-supplied sequence) is in DESCENDING order
    if (length(lambda) > 1 ) {
      if (max(diff(lambda)) > 0) stop("\nUser-supplied lambda sequence must be in descending (largest -> smallest) order")
    }
    nlambda <- length(lambda)
    user.lambda <- TRUE
  }

  # placeholders for results ---------------------------------
  # setting up 'init' as below is not needed when this is called from plmm or cv_plmm
  # as those user-facing functions set up 'init' as a vector of zeroes
  if (is.null(init)){
    init <- rep(0, ncol(stdrot_X))
  }

  if(!fbm_flag){
    r <- drop(rot_y - stdrot_X %*% init)
    # cat("\nhead(residuals, 10):", head(r, 10))
    # cat("\ntail(residuals, 10):", tail(r, 10))
    linear.predictors <- matrix(NA, nrow = nrow(stdrot_X), ncol=nlambda)
    stdrot_scale_beta <- matrix(NA, nrow=ncol(stdrot_X), ncol=nlambda)
  } else {
    r <- rot_y - stdrot_X%*%as.matrix(init) # again, using bigalgebra method here
    # cat("\nhead(residuals, 10):", head(r, 10))
    # cat("\ntail(residuals, 10):", tail(r, 10))
  }

  iter <- integer(nlambda)
  converged <- logical(nlambda)
  loss <- numeric(nlambda)
  # population var is 1 since we re-standardized (pass this arg to fit function)
  xtx <- rep(1, ncol(stdrot_X))

  # main attraction -----------------------------------------------------------
  if (prep$trace) cat("Beginning model fitting.\n")
  if (!fbm_flag) {
    # set up progress bar -- this can take a while
    if(prep$trace){pb <- utils::txtProgressBar(min = 0, max = nlambda, style = 3)}
    for (ll in 1:nlambda){
      lam <- lambda[ll]
      res <- ncvreg::ncvfit(X = stdrot_X,
                            y = rot_y,
                            init = init,
                            r = r,
                            xtx = xtx,
                            penalty = penalty,
                            gamma = gamma,
                            alpha = alpha,
                            lambda= lam,
                            eps = eps,
                            max.iter = max.iter,
                            penalty.factor = penalty.factor,
                            warn = warn)
      stdrot_scale_beta[, ll] <- init <- res$beta
      linear.predictors[,ll] <- stdrot_X%*%(res$beta)
      iter[ll] <- res$iter
      converged[ll] <- ifelse(res$iter < max.iter, TRUE, FALSE)
      loss[ll] <- res$loss
      r <- res$resid
      if(prep$trace){utils::setTxtProgressBar(pb, ll)}
    }

    # reverse the POST-ROTATION standardization on estimated betas
    std_scale_beta <- matrix(0,
                               nrow = nrow(stdrot_scale_beta) + 1,
                               ncol = ncol(stdrot_scale_beta))
    bb <-  stdrot_scale_beta/stdrot_X_details$scale
    std_scale_beta[-1,] <- bb
    std_scale_beta[1,] <- mean(prep$y) - crossprod(stdrot_X_details$center, bb)

  } else {
    res <- biglasso::biglasso_path(
      X = stdrot_X,
      y = rot_y,
      r = r,
      init = init,
      xtx = xtx,
      penalty = penalty,
      lambda = lambda,
      alpha = alpha,
      gamma = gamma,
      eps = eps,
      max.iter = max.iter,
      penalty.factor = penalty.factor,
      ...)

    stdrot_scale_beta <- res$beta
    linear.predictors <- stdrot_X %*% stdrot_scale_beta

    iter <- res$iter
    converged <- ifelse(iter < max.iter, TRUE, FALSE)
    loss <- res$loss
    r <- res$resid

    # reverse the POST-ROTATION standardization on estimated betas
    # NB: the intercept of a PLMM is always the mean of y. We prove this in our methods work.
    std_scale_beta <- Matrix::sparseMatrix(i = rep(1,ncol(stdrot_scale_beta)),
                                             j = 1:ncol(stdrot_scale_beta),
                                             x = mean(prep$y),
                                             dims = c(nrow(stdrot_scale_beta) + 1,
                                                      ncol = ncol(stdrot_scale_beta)))
    bb <-  stdrot_scale_beta/stdrot_X_details$scale
    std_scale_beta[-1,] <- bb
    # browser()
    std_scale_beta[1,] <- mean(prep$y) - crossprod(stdrot_X_details$center, bb)

  }
  if (prep$trace) {
    cat("Model fitting finished at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '\n')
  }

  # eliminate saturated lambda values, if any
  ind <- !is.na(iter)
  iter <- iter[ind]
  converged <- converged[ind]
  lambda <- lambda[ind]
  loss <- loss[ind]
  if (warn & sum(iter) == max.iter) warning("Maximum number of iterations reached")

  ret <- structure(list(
    n = prep$n,
    p = prep$p,
    std_X_n = prep$std_X_n,
    std_X_p = prep$std_X_p,
    y = prep$y,
    s = prep$s,
    U = prep$U,
    # rot_X = rot_X,
    rot_y = rot_y,
    stdrot_X = stdrot_X,
    stdrot_X_details = stdrot_X_details,
    lambda = lambda,
    stdrot_scale_beta = stdrot_scale_beta,
    std_scale_beta = std_scale_beta,
    linear.predictors = linear.predictors,
    eta = prep$eta,
    penalty.factor = penalty.factor,
    iter = iter,
    converged = converged,
    loss = loss,
    penalty = penalty,
    gamma = gamma,
    alpha = alpha,
    ns = prep$ns,
    penalty = penalty,
    nlambda = nlambda,
    eps = eps,
    max.iter = max.iter,
    warn = warn,
    init = init,
    trace = prep$trace))

  if (exists("rot_X")){
    ret$rot_X <- rot_X
  }

  return(ret)
}

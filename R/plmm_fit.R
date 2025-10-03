#' PLMM fit: a function that fits a PLMM using the values returned by plmm_prep()
#'
#' @param prep A list as returned from \code{plmm_prep}
#' @param y    The original (not centered) outcome vector. Need this for intercept estimate
#' @param std_X_details A list with components 'center' (values used to center X), 'scale' (values used to scale X), and 'ns' (indices for nonsignular columns of X)
#' @param eta_star The ratio of variances (passed from plmm())
#' @param penalty_factor A multiplicative factor for the penalty applied to each coefficient. If supplied, penalty_factor must be a numeric vector of length equal to the number of columns of X. The purpose of penalty_factor is to apply differential penalization if some coefficients are thought to be more likely than others to be in the model. In particular, penalty_factor can be 0, in which case the coefficient is always in the model without shrinkage.
#' @param fbm_flag Logical: is std_X an FBM object? Passed from `plmm()`.
#' @param penalty The penalty to be applied to the model. Either "MCP" (the default), "SCAD", or "lasso".
#' @param gamma The tuning parameter of the MCP/SCAD penalty (see details). Default is 3 for MCP and 3.7 for SCAD.
#' @param alpha Tuning parameter for the Mnet estimator which controls the relative contributions from the MCP/SCAD penalty and the ridge, or L2 penalty. alpha=1 is equivalent to MCP/SCAD penalty, while alpha=0 would be equivalent to ridge regression. However, alpha=0 is not supported; alpha may be arbitrarily small, but not exactly 0.
#' @param lambda_min The smallest value for lambda, as a fraction of lambda.max. Default is .001 if the number of observations is larger than the number of covariates and .05 otherwise.
#' @param nlambda Length of the sequence of lambda. Default is 100.
#' @param lambda A user-specified sequence of lambda values. By default, a sequence of values of length nlambda is computed, equally spaced on the log scale.
#' @param eps Convergence threshold. The algorithm iterates until the RMSD for the change in linear predictors for each coefficient is less than eps. Default is \code{1e-4}.
#' @param max_iter Maximum number of iterations (total across entire path). Default is 10000.
#' @param init Initial values for coefficients. Default is 0 for all columns of X.
#' @param warn Return warning messages for failures to converge and model saturation? Default is TRUE.
#' @param ... Additional arguments that can be passed to `biglasso::biglasso_simple_path()`
#'
#' @keywords internal
#'

plmm_fit <- function(prep,
                     y,
                     std_X_details,
                     penalty_factor,
                     fbm_flag,
                     penalty,
                     gamma = 3,
                     alpha = 1,
                     lambda_min,
                     nlambda = 100,
                     lambda,
                     eps = 1e-04,
                     max_iter = 10000,
                     init = NULL,
                     warn = TRUE,
                     ...) {

  # error checking ------------------------------------------------------------
  if (penalty == "MCP" && gamma <= 1) stop("gamma must be greater than 1 for the MC penalty", call. = FALSE)
  if (penalty == "SCAD" && gamma <= 2) stop("gamma must be greater than 2 for the SCAD penalty", call. = FALSE)
  if (nlambda < 2) stop("nlambda must be at least 2", call. = FALSE)
  if (alpha <= 0) stop("alpha must be greater than 0; choose a small positive number instead", call. = FALSE)

  if (prep$trace) {
    cat("Beginning rotation ('preconditioning').\n")
  }

  # rotate data ----------------------------------------------------------------
  if (!fbm_flag) {
    w <- (prep$eta * prep$s + (1 - prep$eta))^(-1/2)
    wUt <- sweep(x = t(prep$U), MARGIN = 1, STATS = w, FUN = "*")
    rot_X <- wUt %*% prep$std_X
    rot_y <- wUt %*% prep$centered_y # remember: we're using the centered outcome vector

    # re-standardize rot_X
    stdrot_info <- standardize_in_memory(rot_X)
    stdrot_X <- stdrot_info$std_X
    stdrot_X_details <- stdrot_info$std_X_details

  } else {
    rot_res <- rotate_filebacked(prep)
    stdrot_X <- rot_res$stdrot_X
    rot_y <- rot_res$rot_y
    stdrot_X_details <- list(center = rot_res$stdrot_X_center,
                             scale = rot_res$stdrot_X_scale)
  }


  if (prep$trace) {
    (cat("Rotation (preconditiong) finished at ",
         format(Sys.time(), "%Y-%m-%d %H:%M:%S\n")))
  }

  # set up lambda -------------------------------------------------------
  if (missing(lambda)) {
    if (prep$trace) cat("Setting up lambda/preparing for model fitting.\n")
    lambda <- setup_lambda(X = stdrot_X,
                           y = rot_y,
                           alpha = alpha,
                           nlambda = nlambda,
                           lambda_min = lambda_min,
                           penalty_factor = penalty_factor)
    user.lambda <- FALSE
  } else {
    # make sure (if user-supplied sequence) is in DESCENDING order
    if (length(lambda) > 1 && max(diff(lambda)) > 0) {
      stop("\nUser-supplied lambda sequence must be in descending (largest -> smallest) order")
    }
    nlambda <- length(lambda)
    user.lambda <- TRUE
  }

  # placeholders for results ---------------------------------
  # setting up 'init' as below is not needed when this is called from plmm or cv_plmm
  # as those user-facing functions set up 'init' as a vector of zeroes
  if (is.null(init)) {
    init <- rep(0, ncol(stdrot_X))
  }

  if (!fbm_flag) {
    r <- drop(rot_y - stdrot_X %*% init)
    stdrot_scale_beta <- matrix(NA, nrow = ncol(stdrot_X), ncol = nlambda)
  } else {
    r <- rot_y - stdrot_X %*% as.matrix(init) # again, using bigalgebra method here
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
    if (prep$trace) {
      pb <- utils::txtProgressBar(min = 0, max = nlambda, style = 3)
    }
    for (ll in 1:nlambda) {
      lam <- lambda[ll]
      res <- ncvreg::ncvfit(X = stdrot_X,
                            y = rot_y,
                            init = init,
                            r = r,
                            xtx = xtx,
                            penalty = penalty,
                            gamma = gamma,
                            alpha = alpha,
                            lambda = lam,
                            eps = eps,
                            max.iter = max_iter,
                            penalty.factor = penalty_factor,
                            warn = warn)

      stdrot_scale_beta[, ll] <- init <- res$beta
      iter[ll] <- res$iter
      converged[ll] <- res$iter < max_iter
      loss[ll] <- res$loss
      r <- res$resid
      if (prep$trace) {
        utils::setTxtProgressBar(pb, ll)
      }
    }
    if (prep$trace) close(pb)

    # reverse the POST-ROTATION standardization on estimated betas
    std_scale_beta <- matrix(0,
                             nrow = nrow(stdrot_scale_beta) + 1,
                             ncol = ncol(stdrot_scale_beta))

    stdrot_unscale <- ifelse(stdrot_X_details$scale < 1e-3, 1, stdrot_X_details$scale)
    bb <-  stdrot_scale_beta / (stdrot_unscale)
    std_scale_beta[-1, ] <- bb
    std_scale_beta[1, ] <- mean(y) - crossprod(stdrot_X_details$center, bb)

    # calculate linear predictors on the scale of std_X
    std_Xbeta <- prep$std_X %*% bb
    std_Xbeta <- sweep(std_Xbeta, 2, std_scale_beta[1, ], "+")

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
      max.iter = max_iter,
      penalty.factor = penalty_factor,
      ...)

    stdrot_scale_beta <- res$beta

    iter <- res$iter
    converged <- iter < max_iter
    loss <- res$loss
    r <- res$resid

    # reverse the POST-ROTATION standardization on estimated betas
    # NB: the intercept of a PLMM is always the mean of y. We prove this in our methods work.
    std_scale_beta <- Matrix::sparseMatrix(i = rep(1, ncol(stdrot_scale_beta)),
                                           j = seq_len(ncol(stdrot_scale_beta)),
                                           x = mean(y),
                                           dims = c(nrow(stdrot_scale_beta) + 1,
                                                    ncol = ncol(stdrot_scale_beta)))
    bb <-  stdrot_scale_beta / stdrot_X_details$scale
    std_scale_beta[-1, ] <- bb
    std_scale_beta[1, ] <- mean(y) - crossprod(stdrot_X_details$center, bb)

    # calculate linear predictors on the scale of std_X
    std_Xbeta <- prep$std_X %*% bb
    std_Xbeta <- sweep(std_Xbeta, 2, std_scale_beta[1, ], "+")
  }

  if (prep$trace) {
    cat("Model fitting finished at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  }

  # eliminate saturated lambda values, if any
  ind <- !is.na(iter)
  iter <- iter[ind]
  converged <- converged[ind]
  lambda <- lambda[ind]
  loss <- loss[ind]
  if (warn && sum(iter) == max_iter) warning("Maximum number of iterations reached")

  ret <- structure(list(
    y = y,
    std_scale_beta = std_scale_beta,
    std_Xbeta = std_Xbeta,
    centered_y = prep$centered_y, # the centered outcome vector
    s = prep$s,
    U = prep$U,
    lambda = lambda,
    penalty = penalty,
    penalty_factor = penalty_factor,
    iter = iter,
    converged = converged,
    loss = loss,
    eta = prep$eta,
    gamma = gamma,
    alpha = alpha,
    nlambda = nlambda,
    eps = eps,
    max_iter = max_iter,
    warn = warn,
    trace = prep$trace))

  if (fbm_flag) {
    ret$std_X <- bigmemory::describe(prep$std_X)
  }

  return(ret)
}

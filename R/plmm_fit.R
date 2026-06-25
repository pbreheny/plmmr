#' PLMM fit: A function that fits a PLMM using the values returned by `plmm_prep()`
#'
#' @param prep A list as returned from `plmm_prep`
#' @param y    The original (not centered) outcome vector. Need this for intercept estimate
#' @param std_X_details A list with components `center` (values used to center X), `scale` (values used to scale X), and `ns` (indices for nonsingular columns of X)
#' @param fbm_flag Logical: is std_X a filebacked `big.matrix` object? Passed from `plmm()`.
#' @param penalty The penalty to be applied to the model. Either "MCP" (the default), "SCAD", or "lasso".
#' @param gamma The tuning parameter of the MCP/SCAD penalty (see details). Default is 3 for MCP and 3.7 for SCAD.
#' @param alpha Tuning parameter for the Mnet estimator which controls the relative contributions from the MCP/SCAD penalty and the ridge, or L2 penalty. `alpha = 1` is equivalent to MCP/SCAD penalty, while `alpha = 0` would be equivalent to ridge regression.
#'              However, `alpha = 0` is not supported; alpha may be arbitrarily small, but not exactly 0.
#' @param lambda_min The smallest value for lambda, as a fraction of the maximum lambda. Default is .001 if the number of observations is larger than the number of covariates and .05 otherwise.
#' @param nlambda Length of the sequence of lambda. Default is 100.
#' @param lambda A user-specified sequence of lambda values. By default, a sequence of values of length `nlambda` is computed, equally spaced on the log scale.
#' @param eps Convergence threshold. The algorithm iterates until the RMSE for the change in linear predictors for each coefficient is less than `eps`. Default is `1e-4`.
#' @param max_iter Maximum number of iterations (total across entire path). Default is 10000.
#' @param dfmax Maximum number of non-zero coefficients that may enter the model. Default is NULL (no maximum).
#' @param init Initial values for coefficients. Default is 0 for all columns of X.
#' @param warn Return warning messages for failures to converge and model saturation? Default is TRUE.
#' @param restandardize Should the X matrix be restandardized after rotation? Default is TRUE.
#' @param ... Additional arguments that can be passed to `biglasso::biglasso_simple_path()`
#'
#' @return A list which includes 21 items:
#'  * `y`: The outcome vector used in model fitting.
#'  * `std_scale_beta`: The matrix of estimated coefficients on the standardized scale. Rows are predictors (with the first row being the intercept), and columns are values of `lambda`.
#'  * `std_Xbeta`: A matrix of the linear predictors on the scale of the standardized design matrix. Rows are predictors, columns are values of `lambda`.
#'              **Note**: `std_Xbeta` will not include rows for the intercept or for constant features.
#'  * `centered_y`: The centered outcome vector.
#'  * `s`: a vector of the non-zero eigenvalues of the relatedness matrix K (note: K is the kinship matrix for genetic/genomic data; see the article on notation for details)
#'  * `U`: a matrix of the eigenvectors of K associated with `s`
#'  * `lambda`: A numeric vector of the tuning parameter values used in model fitting.
#'  * `penalty`: A character string indicating the penalty with which the model was fit (e.g., 'MCP')
#'  * `penalty_factor`: A vector of indicators corresponding to each predictor, where 1 = predictor was penalized.
#'  * `iter`: An integer vector with the number of iterations needed in model fitting for each value of `lambda`
#'  * `converged`: A vector of logical values indicating whether the model fitting converged at each value of `lambda`
#'  * `loss`: A vector with the numeric values of the loss at each value of `lambda` (calculated on the ~rotated~ scale)
#'  * `eta`: A double between 0 and 1 representing the estimated proportion of the variance in the outcome attributable to population/correlation structure.
#'  * `gamma`: A numeric value indicating the tuning parameter used for the SCAD or MCP penalties. Not relevant for lasso models.
#'  * `alpha`: A numeric value indicating the elastic net tuning parameter.
#'  * `nlambda` Length of the sequence of lambda.
#'  * `eps`: Convergence threshold. The algorithm iterates until the RMSE for the change in linear predictors for each coefficient is less than `eps`
#'  * `max_iter`: Maximum number of iterations (total across entire path)
#'  * `warn`: Return warning messages for failures to converge and model saturation?
#'  * `trace`: If set to TRUE, inform the user of progress by announcing the beginning of each step of the modeling process
#'  * `std_X`: If design matrix is filebacked, the descriptor for the filebacked data is returned using `bigmemory::describe()`.
#'
#' @keywords internal
#'
plmm_fit <- function(
  prep,
  y,
  std_X_details,
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
  dfmax = NULL,
  warn = TRUE,
  restandardize = TRUE,
  ...
) {
  # error checking ------------------------------------------------------------
  if (penalty == "MCP" && gamma <= 1) {
    stop("gamma must be greater than 1 for the MC penalty", call. = FALSE)
  }
  if (penalty == "SCAD" && gamma <= 2) {
    stop("gamma must be greater than 2 for the SCAD penalty", call. = FALSE)
  }
  if (nlambda < 2) {
    stop("nlambda must be at least 2", call. = FALSE)
  }
  if (alpha <= 0) {
    stop(
      "alpha must be greater than 0; choose a small positive number instead",
      call. = FALSE
    )
  }

  if (prep$trace) {
    cat("Beginning rotation ('preconditioning').\n")
  }

  # rotate data ----------------------------------------------------------------
  if (!fbm_flag) {
    w <- (prep$eta * prep$s + (1 - prep$eta))^(-1 / 2)
    wUt <- sweep(x = t(prep$U), MARGIN = 1, STATS = w, FUN = "*")
    rot_X <- prep$U %*% wUt %*% prep$std_X
    if (prep$incpt_flag) {
      rot_y <- prep$U %*% wUt %*% y
    } else {
      rot_y <- prep$U %*% wUt %*% prep$centered_y # remember: we're using the centered outcome vector
    }

    if (restandardize) {
      # scale rot_X
      stdrot_info <- standardize_in_memory(rot_X, tocenter = FALSE)
      stdrot_X <- stdrot_info$std_X
      stdrot_X_details <- stdrot_info$std_X_details
    } else {
      stdrot_X <- rot_X
    }
  } else {
    rot_res <- rotate_filebacked(prep, tocenter = FALSE, restandardize = restandardize)
    stdrot_X <- rot_res$stdrot_X
    rot_y <- rot_res$rot_y
    stdrot_X_details <- list(
      center = rot_res$stdrot_X_center,
      scale = rot_res$stdrot_X_scale
    )
  }

  if (prep$trace) {
    (cat(
      "Rotation (preconditioning) finished at ",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S\n")
    ))
  }

  # set up lambda -------------------------------------------------------
  if (missing(lambda)) {
    if (prep$trace) {
      cat("Setting up lambda/preparing for model fitting.\n")
    }
    lambda <- setup_lambda(
      X = stdrot_X,
      y = rot_y,
      alpha = alpha,
      nlambda = nlambda,
      lambda_min = lambda_min,
      penalty_factor = prep$penalty_factor
    )
  } else {
    # make sure (if user-supplied sequence) is in DESCENDING order
    if (length(lambda) > 1 && max(diff(lambda)) > 0) {
      stop(
        "\nUser-supplied lambda sequence must be in descending (largest -> smallest) order"
      )
    }
    nlambda <- length(lambda)
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

  if (restandardize) {
    # population var is 1 since we re-standardized (pass this arg to fit function)
    xtx <- rep(1, ncol(stdrot_X))
  }

  # main attraction -----------------------------------------------------------
  if (prep$trace) {
    cat("Beginning model fitting.\n")
  }
  if (!fbm_flag) {
    # set up progress bar -- this can take a while
    if (prep$trace) {
      pb <- utils::txtProgressBar(min = 0, max = nlambda, style = 3)
    }
    for (ll in 1:nlambda) {
      lam <- lambda[ll]
      if (restandardize) {
        res <- ncvreg::ncvfit(
          X = stdrot_X,
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
          penalty.factor = prep$penalty_factor,
          warn = warn
        )
      } else {
        res <- ncvreg::ncvfit(
          X = stdrot_X,
          y = rot_y,
          init = init,
          r = r,
          penalty = penalty,
          gamma = gamma,
          alpha = alpha,
          lambda = lam,
          eps = eps,
          max.iter = max_iter,
          penalty.factor = prep$penalty_factor,
          warn = warn
        )
      }

      stdrot_scale_beta[, ll] <- init <- res$beta
      iter[ll] <- res$iter
      converged[ll] <- res$iter < max_iter
      loss[ll] <- res$loss
      r <- res$resid
      if (prep$trace) {
        utils::setTxtProgressBar(pb, ll)
      }
      if (sum(res$beta != 0) > dfmax) {
        iter[ll:nlambda] <- NA
        break
      }
    }
    if (prep$trace) {
      close(pb)
    }

    ind <- !is.na(iter)
    iter <- iter[ind]
    converged <- converged[ind]
    loss <- loss[ind]
    lambda <- lambda[ind]
    stdrot_scale_beta <- stdrot_scale_beta[, ind, drop = FALSE]

    # adjust dimensions of beta matrix based on whether the intercept was included
    std_scale_beta <- matrix(0,
                             nrow = nrow(stdrot_scale_beta) + 1 * !(prep$incpt_flag),
                             ncol = ncol(stdrot_scale_beta))
  } else {
    if (restandardize) {
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
        penalty.factor = prep$penalty_factor,
        dfmax = dfmax,
        warn = warn,
        ...)
    } else {
      res <- biglasso::biglasso_path(
        X = stdrot_X,
        y = rot_y,
        r = r,
        init = init,
        penalty = penalty,
        lambda = lambda,
        alpha = alpha,
        gamma = gamma,
        eps = eps,
        max.iter = max_iter,
        penalty.factor = prep$penalty_factor,
        dfmax = dfmax,
        warn = warn,
        ...)
    }

    stdrot_scale_beta <- res$beta

    iter <- res$iter
    converged <- iter < max_iter
    loss <- res$loss
    r <- res$resid
    lambda <- res$lambda

    std_scale_beta <- Matrix::sparseMatrix(
      i = rep(1, ncol(stdrot_scale_beta)),
      j = seq_len(ncol(stdrot_scale_beta)),
      x = mean(y),
      dims = c(nrow(stdrot_scale_beta) + 1, ncol = ncol(stdrot_scale_beta))
    )
  }

  if (restandardize) {
    # reverse the POST-ROTATION standardization on estimated betas
    stdrot_unscale <- ifelse(
      stdrot_X_details$scale < 1e-3,
      1,
      stdrot_X_details$scale
    )
    bb <- stdrot_scale_beta / stdrot_unscale

    if (prep$incpt_flag) {
      std_scale_beta <- bb
      # calculate linear predictors on the scale of std_X
      std_Xbeta <- prep$std_X %*% std_scale_beta
    } else {
      std_scale_beta[-1, ] <- bb
      std_scale_beta[1, ] <- mean(y)

      # calculate linear predictors on the scale of std_X
      std_Xbeta <- prep$std_X %*% bb
      std_Xbeta <- sweep(std_Xbeta, 2, std_scale_beta[1, ], "+")
    }
  } else {
    if (prep$incpt_flag) {
      std_scale_beta <- stdrot_scale_beta
      # calculate linear predictors on the scale of std_X
      std_Xbeta <- prep$std_X %*% std_scale_beta
    } else {
      std_scale_beta[-1, ] <- stdrot_scale_beta
      std_scale_beta[1, ] <- mean(y)

      # calculate linear predictors on the scale of std_X
      std_Xbeta <- prep$std_X %*% stdrot_scale_beta
      std_Xbeta <- sweep(std_Xbeta, 2, std_scale_beta[1, ], "+")
    }
  }

  if (prep$trace) {
    cat(
      "Model fitting finished at ",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      "\n"
    )
  }

  # eliminate saturated lambda values, if any
  ind <- !is.na(iter)
  iter <- iter[ind]
  converged <- converged[ind]
  lambda <- lambda[ind]
  loss <- loss[ind]
  if (warn && sum(iter) == max_iter) {
    warning("Maximum number of iterations reached")
  }

  ret <- structure(list(
    y = y,
    std_scale_beta = std_scale_beta,
    std_Xbeta = std_Xbeta,
    centered_y = prep$centered_y, # the centered outcome vector
    s = prep$s,
    U = prep$U,
    lambda = lambda,
    penalty = penalty,
    penalty_factor = prep$penalty_factor,
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
    trace = prep$trace
  ))

  if (fbm_flag) {
    ret$std_X <- bigmemory::describe(prep$std_X)
  }

  ret
}

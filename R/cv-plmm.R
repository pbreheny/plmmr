#' Cross-validation for plmm
#'
#' Performs k-fold cross validation for lasso-, MCP-, or SCAD-penalized
#'  linear mixed models over a grid of values for the regularization parameter `lambda`.
#'
#' @param design          The first argument must be one of three things:
#'                          (1) `plmm_design` object (as created by `create_design()`)
#'                          (2) a string with the file path to a design object (the file path must end in '.rds')
#'                          (3) a `matrix` or `data.frame` object representing the design matrix of interest
#' @param y               Optional: In the case where `design` is a `matrix` or `data.frame`, the user must also supply
#'                        a numeric outcome vector as the `y` argument. In this case, `design` and `y` will be passed
#'                        internally to `create_design(X = design, y = y)`.
#' @param K               Similarity matrix used to rotate the data. This should either be (1) a known matrix that reflects the covariance of y, (2) an estimate (Default is \eqn{\frac{1}{p}(XX^T)}), or (3) a list with components 's' and 'u', as returned by choose_k().
#' @param diag_K          Logical: should K be a diagonal matrix? This would reflect observations that are unrelated, or that can be treated as unrelated. Defaults to FALSE.
#'                        Note: plmm() does not check to see if a matrix is diagonal. If you want to use a diagonal K matrix, you must set diag_K = TRUE.
#' @param eta_star        Optional argument to input a specific eta term rather than estimate it from the data. If K is a known covariance matrix that is full rank, this should be 1.
#' @param penalty         The penalty to be applied to the model. Either "lasso" (the default), "SCAD", or "MCP".
#' @param gamma           The tuning parameter of the MCP/SCAD penalty (see details). Default is 3 for MCP and 3.7 for SCAD.
#' @param alpha           Tuning parameter for the Mnet estimator which controls the relative contributions from the MCP/SCAD penalty and the ridge, or L2 penalty. alpha=1 is equivalent to MCP/SCAD penalty, while alpha=0 would be equivalent to ridge regression. However, alpha=0 is not supported; alpha may be arbitrarily small, but not exactly 0.
#' @param lambda_min      The smallest value for lambda, as a fraction of lambda.max. Default is .001 if the number of observations is larger than the number of covariates and .05 otherwise.
#' @param nlambda         Length of the sequence of lambda. Default is 100.
#' @param lambda          A user-specified sequence of lambda values. By default, a sequence of values of length nlambda is computed, equally spaced on the log scale.
#' @param eps             Convergence threshold. The algorithm iterates until the RMSD for the change in linear predictors for each coefficient is less than eps. Default is \code{1e-4}.
#' @param max_iter        Maximum number of iterations (total across entire path). Default is 10000.
#' @param init            Initial values for coefficients. Default is 0 for all columns of X.
#' @param warn            Return warning messages for failures to converge and model saturation? Default is TRUE.
#' @param type            A character argument indicating what should be returned from predict.plmm(). If type == 'lp', predictions are
#'                        based on the linear predictor, X beta. If type == 'blup', predictions are based on the sum of the linear predictor
#'                        and the estimated random effect (BLUP). Defaults to 'blup', as this has shown to be a superior prediction method
#'                        in many applications.
#' @param cluster         Option for **in-memory data only**: cv_plmm() can be run in parallel across a cluster using the parallel package.
#'                        The cluster must be set up in advance using parallel::makeCluster(). The cluster must then be passed to cv_plmm().
#'                        **Note**: this option is not yet implemented for filebacked data.
#' @param nfolds          The number of cross-validation folds. Default is 5.
#' @param fold            Which fold each observation belongs to. By default, the observations are randomly assigned.
#' @param seed            You may set the seed of the random number generator in order to obtain reproducible results.
#' @param trace           If set to TRUE, inform the user of progress by announcing the beginning of each CV fold. Default is FALSE.
#' @param save_rds        Optional: if a filepath and name *without* the '.rds' suffix is specified (e.g., `save_rds = "~/dir/my_results"`), then the model results are saved to the provided location (e.g., "~/dir/my_results.rds").
#'                        Defaults to NULL, which does not save the result.
#'                        **Note**: Along with the model results, two '.rds' files ('loss' and 'yhat') will be created in the same directory as 'save_rds'.
#'                        These files contain the loss and predicted outcome values in each fold; both files will be updated during after prediction within each fold.
#' @param return_fit      Optional: a logical value indicating whether the fitted model should be returned as a `plmm` object in the current (assumed interactive) session. Defaults to TRUE.
#' @param ...             Additional arguments to `plmm_fit`
#'
#' @returns A list that includes 15 items:
#'
#' * type: The type of prediction used ('lp' or 'blup')
#' * cve: A numeric vector with the cross validation error (CVE) at each value of `lambda`
#' * cvse: A numeric vector with the estimated standard error associated with each value of `cve`
#' * fold: A numeric `n` length vector of integers indicating the fold to which each observation was assigned
#' * lambda: A numeric vector of `lambda` values
#' * fit: The overall fit of the object, including all predictors; this is a list as returned by `plmm()`
#' * min: The index corresponding to the value of `lambda` that minimizes `cve`
#' * lambda_min: The `lambda` value at which `cve` is minimized
#' * min1se: The index corresponding to the value of `lambda` within 1 standard error of
#'   that which minimizes `cve`
#' * lambda1se: The largest value of lambda such that `cve` is within 1 standard error of the minimum
#' * null.dev: A numeric value representing the deviance for the intercept-only model. If you have supplied
#'   your own `lambda` sequence, this quantity may not be meaningful.
#' * Y: A matrix with the predicted outcome (\eqn{\hat{y}}) values at each value of `lambda`.
#'   Rows are observations, columns are values of `lambda`.
#' * bias: A numeric value with the estimated bias of the minimized CVE.
#' * loss: A matrix with the loss values at each value of lambda. Rows are observations,
#'   columns are values of `lambda`.
#' * estimated_Sigma: An n x n matrix representing the estimated covariance matrix.
#' @export
#'
#' @examples
#' admix_design <- create_design(X = admix$X, y = admix$y)
#' cv_fit <- cv_plmm(design = admix_design)
#' print(summary(cv_fit))
#' plot(cv_fit)
#'
#' # Note: for examples with filebacked data, see the filebacking vignette
#' # https://pbreheny.github.io/plmmr/articles/filebacking.html
#'
#'
cv_plmm <- function(design,
                    y = NULL,
                    K = NULL,
                    diag_K = NULL,
                    eta_star = NULL,
                    penalty = "lasso",
                    type = "blup",
                    gamma,
                    alpha = 1,
                    lambda_min, # passed to internal function setup_lambda()
                    nlambda = 100,
                    lambda,
                    eps = 1e-04,
                    max_iter = 10000,
                    warn = TRUE,
                    init = NULL,
                    cluster,
                    nfolds = 5,
                    seed,
                    fold = NULL,
                    trace = FALSE,
                    save_rds = NULL,
                    return_fit = TRUE,
                    ...) {

  # check filepaths for saving results ------------------------------
  if (!is.null(save_rds)) {
    save_rds <- tools::file_path_sans_ext(save_rds)
    # ^^ internally, we need to take off the extension from the file name

    # start the log
    logfile <- create_log(outfile = save_rds)
  }

  # create a design if needed -------------
  if (inherits(design, "data.frame") || inherits(design, "matrix")) {
    # error check: if 'design' is matrix/data.frame, user must specify 'y'
    if (is.null(y)) {
      stop("If you supply a matrix or data frame as 'design', you
                         must also specify a numeric vector as the outcome of your
                         model to the 'y' argument.")
    }

    design <- create_design_in_memory(X = design, y = y)


  } else if (!is.null(y)) {
    stop("If you are supplying a plmm_design object or filepath to
                          the 'design' argument, that design already has a 'y' --
                          please do not specify a 'y' argument here in plmm()")
  }

  # run data checks ------------------------------
  checked_data <- plmm_checks(design = design,
                              K = K,
                              diag_K = diag_K,
                              eta_star = eta_star,
                              penalty = penalty,
                              init = init,
                              gamma = gamma,
                              alpha = alpha,
                              trace = trace,
                              ...)

  if (!is.null(save_rds)) {
    cat("\nInput data passed all checks at ",
        pretty_time(),
        file = logfile, append = TRUE)
  }

  # prep  ------------------------
  prep_args <- c(list(std_X = checked_data$std_X,
                      std_X_n = checked_data$std_X_n,
                      std_X_p = checked_data$std_X_p,
                      n = checked_data$n,
                      p = checked_data$p,
                      centered_y = checked_data$centered_y,
                      K = checked_data$K,
                      diag_K = checked_data$diag_K,
                      eta_star = checked_data$eta_star,
                      fbm_flag = checked_data$fbm_flag,
                      trace = trace))

  prep <- do.call("plmm_prep", prep_args)

  if (!is.null(save_rds)) {
    cat("\nEigendecomposition finished at ",
        pretty_time(),
        file = logfile, append = TRUE)
  }

  # full model fit ----------------------------------
  fit_args <- c(list(
    prep = prep,
    y = checked_data$y,
    std_X_details = checked_data$std_X_details,
    penalty_factor = checked_data$penalty_factor,
    fbm_flag = checked_data$fbm_flag,
    penalty = checked_data$penalty,
    gamma = checked_data$gamma,
    alpha = alpha,
    nlambda = nlambda,
    max_iter = max_iter,
    eps = eps,
    warn = warn))

  if (!missing(lambda_min)) {
    fit_args$lambda_min <- lambda_min
  }

  fit <- do.call("plmm_fit", fit_args)

  if (!is.null(save_rds)) {
    cat("\nFull model fit finished at",
        pretty_time(),
        file = logfile, append = TRUE)
  }

  fit_to_return <- plmm_format(fit = fit,
                               p = checked_data$p,
                               std_X_details = checked_data$std_X_details,
                               fbm_flag = checked_data$fbm_flag,
                               plink_flag = checked_data$plink_flag)

  if (!is.null(save_rds)) {
    cat("\nFormatting for full model finished at",
        pretty_time(),
        file = logfile, append = TRUE)
  }

  if (!is.null(save_rds)) {
    saveRDS(fit_to_return, paste0(save_rds, ".rds"))
    cat("Results saved to:", paste0(save_rds, ".rds"), "at",
        pretty_time(),
        file = logfile, append = TRUE)
  }
  gc()

  # set up arguments for cv ---------------------------
  cv_args <- fit_args
  cv_args$warn <- FALSE
  cv_args$lambda <- fit$lambda
  cv_args$eta_star <- fit$eta # note: here we use the same eta estimate in each fold of CV
  cv_args$plink_flag <- checked_data$plink_flag

  estimated_Sigma <- NULL
  if (type == "blup") {
    estimated_Sigma <- construct_variance(eta = fit$eta, K = prep$K)
  }

  # initialize objects to hold CV results
  n <- checked_data$std_X_n
  E <- Y <- matrix(NA, nrow = n, ncol = length(fit$lambda))

  # set up folds for cross validation
  if (!missing(seed)) {
    original_seed <- .GlobalEnv$.Random.seed
    on.exit(.GlobalEnv$.Random.seed <- original_seed)
    set.seed(seed)
  }

  #sde <- sqrt(.Machine$double.eps)

  if (is.null(fold)) {
    if (trace) {
      cat("'Fold' argument is either NULL or missing; assigning folds randomly (by default).
          \nTo specify folds for each observation, supply a vector with fold assignments.\n")
    }
    fold <- sample(1:n %% nfolds)
    fold[fold == 0] <- nfolds
  } else {
    nfolds <- max(fold)
  }

  # set up cluster if user-specified ---------------------------------------
  if (!missing(cluster)) {

    # error check: parallelization is not yet implemented for filebacked data
    if (checked_data$fbm_flag) stop("Parallelization is not yet implemented for filebacked data.
                                    You cannot specify a 'cluster' argument if data are stored filebacked.\n.")

    # check type of 'cluster' argument
    if (!inherits(cluster, "cluster")) stop("cluster is not of class 'cluster'; see ?makeCluster", call. = FALSE)

    # check if variables are defined
    if (!exists("fold") || !exists("type") || !exists("cv_args")) {
      stop("One or more required variables (fold, type, cv_args) are not defined", call. = FALSE)
    }

    parallel::clusterExport(cl = cluster,
                            varlist = c("fold", "type", "cv_args"),
                            envir = environment())
    parallel::clusterCall(cluster, function() library(plmmr))
    fold_results <- parallel::parLapply(cl = cluster,
                                        X = 1:max(fold),
                                        fun = cvf,
                                        design = design,
                                        fold = fold,
                                        type = type,
                                        cv_args = cv_args)

    # stop the cluster when done
    parallel::stopCluster(cluster)

  }

  # carry out CV -------------------------------------
  if (trace) cat("\nStarting cross validation\n")

  if (!is.null(save_rds)) {
    cat("\nCross validation started at: ",
        pretty_time(),
        file = logfile,
        append = TRUE)
  }

  # set up progress bar -- this can take a while
  if(trace) {
    pb <- utils::txtProgressBar(min = 0, max = nfolds, style = 3)
  }
  for (i in 1:nfolds) {
    # case 1: user-specified cluster
    if (!missing(cluster)) {
      res <- fold_results[[i]] # refers to lines from above
      if (trace) {
        utils::setTxtProgressBar(pb, i)
        }
    } else {
      # case 2: cluster NOT user specified
      if (!is.null(save_rds)) {
        cat("Started fold", i, "at", pretty_time(),
            file = logfile, append = TRUE)
      }

      res <- cvf(i = i,
                 fold = fold,
                 type = type,
                 cv_args = cv_args)
      if (trace) {
        utils::setTxtProgressBar(pb, i)
        }
    }

    if(trace) close(pb)

    # update E and Y
    E[fold == i, 1:res$nl] <- res$loss

    if (!is.matrix(res$yhat)) {
      res$yhat <- as.matrix(res$yhat)
    }
    Y[fold == i, 1:res$nl] <- res$yhat

    if (!is.null(save_rds)) {
      saveRDS(E, paste0(save_rds, "_loss.rds"))
      cat("Loss saved to:", paste0(save_rds, "_loss.rds"), "at", pretty_time(),
          file = logfile, append = TRUE)

      saveRDS(Y, paste0(save_rds, "_yhat.rds"))
      cat("Predicted outcomes saved to:", paste0(save_rds, "_yhat.rds"), "at", pretty_time(),
          file = logfile, append = TRUE)
    }
    gc()
  }

  # post-process results -----------------------------------------
  # eliminate saturated lambda values, if any
  ind <- which(apply(is.finite(E), 2, all)) # index for lambda values to keep
  E <- E[, ind, drop = FALSE]
  Y <- Y[, ind]
  lambda <- fit$lambda[ind]

  # return min lambda idx
  cve <- apply(E, 2, mean)
  cvse <- apply(E, 2, stats::sd) / sqrt(nrow(Y))
  min <- which.min(cve)

  # return lambda 1se idx
  l.se <- cve[min] - cvse[min]
  u.se <- cve[min] + cvse[min]
  within1se <- which(cve >= l.se & cve <= u.se)
  min1se <- which.max(lambda %in% lambda[within1se])

  # bias correction
  e <- sapply(1:nfolds, function(i) apply(E[fold == i, , drop = FALSE], 2, mean))
  Bias <- mean(e[min, ] - apply(e, 2, min))

  val <- list(type = type,
              cve = cve,
              cvse = cvse,
              fold = fold,
              lambda = lambda,
              fit = fit_to_return,
              min = min,
              lambda_min = lambda[min],
              min1se = min1se,
              lambda.1se = lambda[min1se],
              null.dev = mean(plmm_loss(checked_data$y, rep(mean(checked_data$y), n))),
              Y = Y,
              bias = Bias,
              loss = E)

  if (type == "blup") {
    val$estimated_Sigma <- estimated_Sigma
  }

  # save output
  if (!is.null(save_rds)) {
    saveRDS(val, paste0(save_rds, ".rds"))
    cat("Results saved to:", paste0(save_rds, ".rds"), "at",
        pretty_time(),
        file = logfile, append = TRUE)
  }

  # create a failsafe -- if save_rds is NULL, make sure return_fit = TRUE
  if (is.null(save_rds) && !return_fit) {
    cat("You accidentally left save_rds = NULL and return_fit = FALSE;
        to prevent you from losing your work, plmm() is returning the output as if return_fit = TRUE")

    return_fit <- TRUE
  }

  # release pointer
  gc()

  # return ----------------------------------------------------------------
  if (return_fit) {
    return(structure(val, class = "cv_plmm"))
  }



}

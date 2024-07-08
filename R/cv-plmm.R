#' Cross-validation for plmm
#'
#' Performs k-fold cross validation for lasso-, MCP-, or SCAD-penalized
#'  linear mixed models over a grid of values for the regularization parameter `lambda`.
#'
#' @param X               Design matrix for model fitting. May include clinical covariates and other non-SNP data.
#' @param y               Continuous outcome vector. Defaults to NULL, assuming that the outcome is the 6th column in the .fam PLINK file data. Can also be a user-supplied numeric vector.
#' @param col_names       Optional vector of column names for design matrix. Defaults to NULL.
#' @param non_genomic     Optional vector specifying which columns of the design matrix represent features that are *not* genomic, as these features are excluded from the empirical estimation of genomic relatedness.
#'                        For cases where X is a filepath to an object created by `process_plink()`, this is handled automatically via the arguments to `process_plink()`.
#'                        For all other cases, 'non_genomic' defaults to NULL (meaning `plmm()` will assume that all columns of `X` represent genomic features).
#' @param K               Similarity matrix used to rotate the data. This should either be (1) a known matrix that reflects the covariance of y, (2) an estimate (Default is \eqn{\frac{1}{p}(XX^T)}), or (3) a list with components 'd' and 'u', as returned by choose_k().
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
#' @param convex          (future idea; not yet incorporated) Calculate index for which objective function ceases to be locally convex? Default is TRUE.
#' @param dfmax           (future idea; not yet incorporated) Upper bound for the number of nonzero coefficients. Default is no upper bound. However, for large data sets, computational burden may be heavy for models with a large number of nonzero coefficients.
#' @param penalty_factor  A multiplicative factor for the penalty applied to each coefficient.
#'                        If supplied, penalty_factor must be a numeric vector of length equal to the number of columns of X.
#'                        The purpose of penalty_factor is to apply differential penalization if some coefficients are thought to be more likely than others to be in the model.
#'                        In particular, penalty_factor can be 0, in which case the coefficient is always in the model without shrinkage.
#' @param init            Initial values for coefficients. Default is 0 for all columns of X.
#' @param warn            Return warning messages for failures to converge and model saturation? Default is TRUE.
#' @param type            A character argument indicating what should be returned from predict.plmm(). If type == 'lp', predictions are
#'                        based on the linear predictor, X beta. If type == 'blup', predictions are based on the sum of the linear predictor
#'                        and the estimated random effect (BLUP). Defaults to 'blup', as this has shown to be a superior prediction method
#'                        in many applications.
#' @param cluster         cv_plmm() can be run in parallel across a cluster using the parallel package. The cluster must be set up in
#'                        advance using parallel::makeCluster(). The cluster must then be passed to cv_plmm().
#' @param nfolds          The number of cross-validation folds. Default is 5.
#' @param fold            Which fold each observation belongs to. By default, the observations are randomly assigned.
#' @param seed            You may set the seed of the random number generator in order to obtain reproducible results.
#' @param returnY         Should cv_plmm() return the linear predictors from the cross-validation folds? Default is FALSE; if TRUE,
#'                        this will return a matrix in which the element for row i, column j is the fitted value for observation i from
#'                        the fold in which observation i was excluded from the fit, at the jth value of lambda.
#' @param returnBiasDetails Logical: should the cross-validation bias (numeric value) and loss (n x p matrix) be returned? Defaults to FALSE.
#' @param trace           If set to TRUE, inform the user of progress by announcing the beginning of each CV fold. Default is FALSE.
#' @param save_rds        Optional: if a filepath and name *without* the '.rds' suffix is specified (e.g., `save_rds = "~/dir/my_results"`), then the model results are saved to the provided location (e.g., "~/dir/my_results.rds").
#'                        Defaults to NULL, which does not save the result.
#' @param save_fold_res   Optional: a logical value indicating whether the results (loss and predicted values) from each CV fold should be saved?
#'                        If TRUE, then two '.rds' files will be saved ('loss' and 'yhat') will be created in the same directory as 'save_rds'.
#'                        Both files will be updated after each fold is done.
#'                        Defaults to FALSE.
#' @param return_fit      Optional: a logical value indicating whether the fitted model should be returned as a `plmm` object in the current (assumed interactive) session. Defaults to TRUE.
#' @param compact_save    Optional: if TRUE, three separate .rds files will saved: one with the 'beta_vals', one with 'K', and one with everything else (see below).
#'                        Defaults to FALSE. **Note**: you must specify `save_rds` for this argument to be called.
#' @param ...             Additional arguments to `plmm_fit`
#'
#' @returns a list with 11 items:
#'
#' * type: the type of prediction used ('lp' or 'blup')
#' * cve: numeric vector with the cross validation error (CVE) at each value of `lambda`
#' * cvse: numeric vector with the estimated standard error associated with each value of for `cve`
#' * fold: numeric `n` length vector of integers indicating the fold to which each observation was assigned
#' * lambda: numeric vector of `lambda` values
#' * fit: the overall fit of the object, including all predictors; this is a
#'  list as returned by `plmm()`
#' * min: The index corresponding to the value of `lambda` that minimizes `cve`
#' * lambda_min: The `lambda` value at which `cve` is minmized
#' * min1se: The index corresponding to the value of `lambda` within
#' standard error of that which minimizes `cve`
#' * lambda1se: largest value of lambda such that error is within 1 standard error of the minimum.
#' * null.dev: numeric value representing the deviance for the
#'  intercept-only model. If you have supplied your own `lambda` sequence,
#'  this quantity may not be meaningful.
#' @export
#'
#' @examples
#' cv_fit <- cv_plmm(X = cbind(admix$race,admix$X), y = admix$y,
#'  non_genomic = 1, penalty_factor = c(0, rep(1, ncol(admix$X))))
#' print(summary(cv_fit))
#' plot(cv_fit)
#'
#' # Note: for examples with filebacked data, see the filebacking vignette
#' # https://pbreheny.github.io/plmmr/articles/filebacking.html
#'
#'
cv_plmm <- function(X,
                    y = NULL,
                    col_names = NULL,
                    non_genomic = NULL,
                    K = NULL,
                    diag_K = NULL,
                    eta_star = NULL,
                    penalty = "lasso",
                    penalty_factor = NULL,
                    type = 'blup',
                    gamma,
                    alpha = 1,
                    lambda_min, # passed to internal function setup_lambda()
                    nlambda = 100,
                    lambda,
                    eps = 1e-04,
                    max_iter = 10000,
                    convex = TRUE,
                    dfmax = NULL,
                    warn = TRUE,
                    init = NULL,
                    cluster,
                    nfolds=5,
                    seed,
                    fold = NULL,
                    returnY=FALSE,
                    returnBiasDetails = FALSE,
                    trace=FALSE,
                    save_rds = NULL,
                    save_fold_res = FALSE,
                    return_fit = TRUE,
                    compact_save = FALSE,
                    ...) {

  # check filepaths for saving results ------------------------------
  if (save_fold_res & is.null(save_rds)) {
    stop("You have set 'save_fold_res = TRUE', but no argument was supplied to 'save_rds'.
         \nPlease specify a filepath (as a string) to 'save_rds'")
  }

  if (compact_save & is.null(save_rds)) {
    stop("You have set 'compact_save = TRUE', but no argument was supplied to 'save_rds'.
          \nPlease specify a filepath (as a string) to 'save_rds'")
  }

  save_rds <- check_for_file_extension(save_rds)
  # ^^ internally, we need to take off the extension from the file name

  # start the log -----------------------
  logfile <- create_log(outfile = ifelse(!is.null(save_rds),
                                         save_rds,
                                         "./cv-plmm"))

  # run data checks ------------------------------
  checked_data <- plmm_checks(X,
                              col_names = col_names,
                              non_genomic = non_genomic,
                              y = y,
                              K = K,
                              diag_K = diag_K,
                              eta_star = eta_star,
                              penalty = penalty,
                              penalty_factor = penalty_factor,
                              init = init,
                              dfmax = dfmax,
                              gamma = gamma,
                              alpha = alpha,
                              trace = trace,
                              ...)

  cat("\nInput data passed all checks at ",
      pretty_time(),
      file = logfile, append = TRUE)
  # prep  ------------------------
  prep_args <- c(list(std_X = checked_data$std_X,
                      std_X_n = checked_data$std_X_n,
                      std_X_p = checked_data$std_X_p,
                      genomic = checked_data$genomic,
                      n = checked_data$n,
                      p = checked_data$p,
                      centered_y = checked_data$centered_y,
                      K = checked_data$K,
                      diag_K = checked_data$diag_K,
                      eta_star = eta_star,
                      fbm_flag = checked_data$fbm_flag,
                      trace = trace))

  prep <- do.call('plmm_prep', prep_args)
  cat("\nEigendecomposition finished at ",
      pretty_time(),
      file = logfile, append = TRUE)

  # full model fit ----------------------------------
  fit_args <- c(list(
    prep = prep,
    y = checked_data$y,
    std_X_details = checked_data$std_X_details,
    eta_star = eta_star,
    penalty_factor = checked_data$penalty_factor,
    fbm_flag = checked_data$fbm_flag,
    penalty = checked_data$penalty,
    gamma = checked_data$gamma,
    alpha = alpha,
    nlambda = nlambda,
    max_iter = max_iter,
    eps = eps,
    warn = warn,
    convex = convex,
    dfmax = dfmax))

  if (!missing(lambda_min)){
    fit_args$lambda_min <- lambda_min
  }
  fit <- do.call('plmm_fit', fit_args)

  cat("\nFull model fit finished at",
      pretty_time(),
      file = logfile, append = TRUE)

  if (is.null(col_names)){
    if (!is.null(checked_data$dat)) {
      col_names <- checked_data$dat$X_colnames
    }
  }

  fit_to_return <- plmm_format(fit = fit,
                               p = checked_data$p,
                               std_X_details = checked_data$std_X_details,
                               feature_names = col_names,
                               fbm_flag = checked_data$fbm_flag,
                               non_genomic = checked_data$non_genomic)

  cat("\nFormatting for full model finished at",
      pretty_time(),
      file = logfile, append = TRUE)

  if (!is.null(save_rds)) {
    if (compact_save) {
      # go ahead and save the most important results from the full model fit

      saveRDS(fit_to_return$beta_vals, paste0(save_rds, "_coefficients.rds"))
      cat("Coefficients (estimated beta values) saved to:", paste0(save_rds, "_coefficients.rds"), "at",
          pretty_time(),
          file = logfile, append = TRUE)

      saveRDS(fit_to_return$K, paste0(save_rds, "_K.rds"))
      cat("K (eigendecomposition) saved to:", paste0(save_rds, "_K.rds"), "at",
          pretty_time(),
          file = logfile, append = TRUE)
    }
  }


  # set up arguments for cv ---------------------------
  cv_args <- fit_args
  cv_args$warn <- FALSE
  cv_args$lambda <- fit$lambda
  cv_args$non_genomic <- checked_data$non_genomic

  estimated_Sigma <- NULL
  if (type == 'blup') {
    estimated_Sigma <- construct_variance(eta = fit$eta, K = prep$K)
  }

  # initialize objects to hold CV results
  n <- checked_data$std_X_n
  E <- Y <- matrix(NA, nrow=n, ncol=length(fit$lambda))


  # set up folds for cross validation
  if (!missing(seed)) {
    original_seed <- .GlobalEnv$.Random.seed
    on.exit(.GlobalEnv$.Random.seed <- original_seed)
    set.seed(seed)
  } else {
    # TODO: determine how, if at all, the random seed should be documented
    # cat("Random seed:",
    #     .GlobalEnv$.Random.seed[1],
    #     "\n",
    #     file = logfile,
    #     append = TRUE)
  }

  sde <- sqrt(.Machine$double.eps)

  if (is.null(fold)) {
    if (trace) {
      cat("'Fold' argument is either NULL or missing; assigning folds randomly (by default).
          \nTo specify folds for each observation, supply a vector with fold assignments.\n")
    }
    fold <- sample(1:n %% nfolds)
    fold[fold==0] <- nfolds
  } else {
    nfolds <- max(fold)
  }

  # set up cluster if user-specified ---------------------------------------
  if (!missing(cluster)) {
    if (!inherits(cluster, "cluster")) stop("cluster is not of class 'cluster'; see ?makeCluster", call.=FALSE)
    parallel::clusterExport(cluster, c("X", "y", "K", "fold", "type", "cv_args", "estimated_Sigma"), envir=environment())
    parallel::clusterCall(cluster, function() library(plmmr))
    fold.results <- parallel::parLapply(cl=cluster, X=1:max(fold), fun=cvf, X=X, y=y,
                                        fold=fold, type=type, cv_args=cv_args,
                                        estimated_Sigma = estimated_Sigma)
  }

  # carry out CV -------------------------------------
  if (trace) cat("\nStarting cross validation\n")
  cat("\nCross validation started at: ",
      pretty_time(),
      file = logfile,
      append = TRUE)

  # set up progress bar -- this can take a while
  if(trace){pb <- utils::txtProgressBar(min = 0, max = nfolds, style = 3)}
  for (i in 1:nfolds) {
    # case 1: user-specified cluster
    if (!missing(cluster)) {
      res <- fold.results[[i]] # refers to lines from above
      if (trace) {utils::setTxtProgressBar(pb, i)}

    } else {
      # case 2: cluster NOT user specified (this is the typical use case)
      cat("Started fold", i, "at", pretty_time(),
          file = logfile, append = TRUE)
      res <- cvf(i = i,
                 fold = fold,
                 type = type,
                 cv_args = cv_args,
                 estimated_Sigma = estimated_Sigma)
      if (trace) {utils::setTxtProgressBar(pb, i)}

    }

    if(trace) close(pb)

    # update E and Y
    E[fold==i, 1:res$nl] <- res$loss

    if (!is.matrix(res$yhat)) {
      res$yhat <- as.matrix(res$yhat)
    }
    Y[fold==i, 1:res$nl] <- res$yhat

    if (save_fold_res) {
      saveRDS(E, paste0(save_rds, "_loss.rds"))
      cat("Loss saved to:", paste0(save_rds, "_loss.rds"), "at", pretty_time(),
          file = logfile, append = TRUE)

      saveRDS(Y, paste0(save_rds, "_yhat.rds"))
      cat("Predicted outcomes saved to:", paste0(save_rds, "_yhat.rds"), "at", pretty_time(),
          file = logfile, append = TRUE)
    }

  }

  # post-process results -----------------------------------------
  # eliminate saturated lambda values, if any
  ind <- which(apply(is.finite(E), 2, all)) # index for lambda values to keep
  E <- E[, ind, drop=FALSE]
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
  e <- sapply(1:nfolds, function(i) apply(E[fold==i, , drop=FALSE], 2, mean))
  Bias <- mean(e[min,] - apply(e, 2, min))

  val <- list(type=type,
              cve=cve,
              cvse=cvse,
              fold=fold,
              lambda=lambda,
              fit=fit_to_return,
              min=min,
              lambda_min=lambda[min],
              min1se = min1se,
              lambda.1se = lambda[min1se],
              null.dev=mean(plmm_loss(checked_data$y, rep(mean(checked_data$y), n))))

  if (returnY) val$Y <- Y
  if (returnBiasDetails){
    val$Bias <- Bias
    val$Loss <- E
  }

  if (type == "blup"){
    val$estimated_Sigma = estimated_Sigma
  }

  # handle output
  if (!is.null(save_rds)){
    if (compact_save) {
      # save the rest of the output across multiple files
      saveRDS(val[c(1:5, 7:11)], paste0(save_rds, "_cv_details.rds"))
      cat("CV details (CVE, fold assignments, etc) saved to:",
          paste0(save_rds, "_cv_details.rds"),
          "at",
          pretty_time(),
          file = logfile, append = TRUE)

      saveRDS(estimated_Sigma, paste0(save_rds, "_estimated_Sigma.rds"))
      cat("Estimated variance matrix saved to:", paste0(save_rds, "_estimated_Sigma.rds"),
          "at", pretty_time(), file = logfile, append = TRUE)

      saveRDS(fit_to_return$linear_predictors, paste0(save_rds, "_linear_predictors.rds"))
      cat("Linear predictors (on rotated scale) saved to:", paste0(save_rds, "_linear_predictors.rds"), "at",
          pretty_time(),
          file = logfile, append = TRUE)

      saveRDS(fit_to_return[c(2:3, 5:12)], paste0(save_rds, "_full_fit_details.rds"))
      cat("All other results (loss, # of iterations, ...) saved to:", paste0(save_rds, "_full_fit_details.rds"), "at",
          pretty_time(),
          file = logfile, append = TRUE)

    } else {
      # save all output in one file (default)
      saveRDS(val, paste0(save_rds, ".rds"))
      cat("Results saved to:", paste0(save_rds, ".rds"), "at",
          pretty_time(),
          file = logfile, append = TRUE)
    }


  }

  if (is.null(save_rds) & !return_fit){
    cat("\nYou accidentally left save_rds NULL while setting return_fit = FALSE;
        to prevent you from losing your work, I am saving the output as cv_plmm_results.rds
        in your current working directory (current folder).
        \nNext time, make sure to specify your own filepath to the save_rds argument.")

    rdsfile <- paste0(getwd(),"/cv_plmm_results.rds")
    saveRDS(val, rdsfile)
    cat("Results saved in", rdsfile, "at", pretty_time(),
        file = logfile, append = TRUE)
  }

  if (return_fit){
    return(structure(val, class="cv_plmm"))
  }

}

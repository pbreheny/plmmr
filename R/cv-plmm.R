#' Cross-validation for plmm
#'
#' Performs k-fold cross validation for lasso-, MCP-, or SCAD-penalized 
#'  linear mixed models over a grid of values for the regularization parameter `lambda`.
#'
#' @param X              Design matrix for model fitting. May include clinical covariates and other non-SNP data.
#' @param y              Continuous outcome vector. Defaults to NULL, assuming that the outcome is the 6th column in the .fam PLINK file data. Can also be a user-supplied numeric vector. 
#' @param std_needed     Logical: does the supplied X need to be standardized? Defaults to TRUE. For data processed from PLINK files, standardization happens in `process_plink()`. For data supplied as a matrix, standardization happens here in `plmm()`. If you know your data are already standardized, set `std_needed = FALSE` -- this would be an atypical case. **Note**: failing to standardize data will lead to incorrect analyses. 
#'                       By default, X will be standardized internally. For data processed from PLINK files, standardization happens in `process_plink()`. For data supplied as a matrix, standardization happens here in `plmm()`. If you know your data are already standardized, set `std_needed = FALSE` -- this would be an atypical case. **Note**: failing to standardize data will lead to incorrect analyses. 
#' @param col_names      Optional vector of column names for design matrix. Defaults to NULL.
#' @param k An integer specifying the number of singular values to be used in the approximation of the rotated design matrix. This argument is passed to `RSpectra::svds()`. Defaults to `min(n, p) - 1`, where n and p are the dimensions of the _standardized_ design matrix.
#' @param K Similarity matrix used to rotate the data. This should either be (1) a known matrix that reflects the covariance of y, (2) an estimate (Default is \eqn{\frac{1}{p}(XX^T)}), or (3) a list with components 'd' and 'u', as returned by choose_k().
#' @param diag_K Logical: should K be a diagonal matrix? This would reflect observations that are unrelated, or that can be treated as unrelated. Defaults to FALSE. 
#'  Note: plmm() does not check to see if a matrix is diagonal. If you want to use a diagonal K matrix, you must set diag_K = TRUE.
#' @param eta_star Optional argument to input a specific eta term rather than estimate it from the data. If K is a known covariance matrix that is full rank, this should be 1.
#' @param penalty The penalty to be applied to the model. Either "MCP" (the default), "SCAD", or "lasso".
#' @param gamma The tuning parameter of the MCP/SCAD penalty (see details). Default is 3 for MCP and 3.7 for SCAD.
#' @param alpha Tuning parameter for the Mnet estimator which controls the relative contributions from the MCP/SCAD penalty and the ridge, or L2 penalty. alpha=1 is equivalent to MCP/SCAD penalty, while alpha=0 would be equivalent to ridge regression. However, alpha=0 is not supported; alpha may be arbitrarily small, but not exactly 0.
#' @param lambda.min The smallest value for lambda, as a fraction of lambda.max. Default is .001 if the number of observations is larger than the number of covariates and .05 otherwise.
#' @param nlambda Length of the sequence of lambda. Default is 100. 
#' @param lambda A user-specified sequence of lambda values. By default, a sequence of values of length nlambda is computed, equally spaced on the log scale.
#' @param eps Convergence threshold. The algorithm iterates until the RMSD for the change in linear predictors for each coefficient is less than eps. Default is \code{1e-4}.
#' @param max.iter Maximum number of iterations (total across entire path). Default is 10000.
#' @param convex Calculate index for which objective function ceases to be locally convex? Default is TRUE.
#' @param dfmax Upper bound for the number of nonzero coefficients. Default is no upper bound. However, for large data sets, computational burden may be heavy for models with a large number of nonzero coefficients.
#' @param penalty.factor A multiplicative factor for the penalty applied to each coefficient. If supplied, penalty.factor must be a numeric vector of length equal to the number of columns of X. The purpose of penalty.factor is to apply differential penalization if some coefficients are thought to be more likely than others to be in the model. In particular, penalty.factor can be 0, in which case the coefficient is always in the model without shrinkage.
#' @param init Initial values for coefficients. Default is 0 for all columns of X. 
#' @param warn Return warning messages for failures to converge and model saturation? Default is TRUE.
#' @param type           A character argument indicating what should be returned from predict.plmm(). If type == 'lp', predictions are 
#'                      based on the linear predictor, X beta. If type == 'blup', predictions are based on the sum of the linear predictor 
#'                      and the estimated random effect (BLUP). Defaults to 'blup', as this has shown to be a superior prediction method 
#'                      in many applications.
#' @param cluster        cv.plmm() can be run in parallel across a cluster using the parallel package. The cluster must be set up in 
#'                      advance using parallel::makeCluster(). The cluster must then be passed to cv.plmm().
#' @param nfolds         The number of cross-validation folds. Default is 10.
#' @param fold           Which fold each observation belongs to. By default, the observations are randomly assigned.
#' @param seed           You may set the seed of the random number generator in order to obtain reproducible results.
#' @param returnY        Should cv.plmm() return the linear predictors from the cross-validation folds? Default is FALSE; if TRUE, 
#'                      this will return a matrix in which the element for row i, column j is the fitted value for observation i from 
#'                      the fold in which observation i was excluded from the fit, at the jth value of lambda.
#' @param returnBiasDetails Logical: should the cross-validation bias (numeric value) and loss (n x p matrix) be returned? Defaults to FALSE.
#' @param trace          If set to TRUE, inform the user of progress by announcing the beginning of each CV fold. Default is FALSE.
#' @param ...            Additional arguments to `plmm_fit`

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
#' * lambda.min: The `lambda` value at which `cve` is minmized 
#' * min1se: The index corresponding to the value of `lambda` within 
#' standard error of that which minimizes `cve`  
#' * lambda1se: largest value of lambda such that error is within 1 standard error of the minimum.
#' * null.dev: numeric value representing the deviance for the
#'  intercept-only model. If you have supplied your own `lambda` sequence,
#'  this quantity may not be meaningful.
#' @export
#' 
#' @examples 
#' cv_fit <- cv.plmm(X = admix$X, y = admix$y, seed = 321)
#' \donttest{
#' cv_s <- summary.cv.plmm(cv_fit, lambda = "1se")
#' print(cv_s)
#' plot(cv_fit)
#' 
#' # filebacked example (file path is specific to current machine)
#' # since this dataset is < 100Mb, have to specify returnX = FALSE for 
#' # `get_data()` to return an FBM
#' my_fb_data <- paste0(get_example_data(parent = TRUE), "/penncath_lite")
#' 
#' cv_fb_fit <- cv.plmm(X = my_fb_data, type = 'blup', returnX = FALSE,
#'  trace = TRUE, nfolds = 3)
#' 
#' plot(cv_fb_fit)
#' summary(cv_fb_fit)
#' }
#' 
#' 
#' 
cv.plmm <- function(X, 
                    y = NULL,
                    std_needed = TRUE,
                    col_names = NULL,
                    k = NULL,
                    K = NULL,
                    diag_K = NULL,
                    eta_star = NULL,
                    penalty = "MCP",
                    penalty.factor = NULL,
                    type = 'blup',
                    gamma,
                    alpha = 1,
                    lambda.min, # passed to internal function setup_lambda()
                    nlambda = 100,
                    lambda,
                    eps = 1e-04,
                    max.iter = 10000,
                    convex = TRUE,
                    dfmax = NULL,
                    warn = TRUE,
                    init = NULL,
                    cluster,
                    nfolds=10,
                    seed,
                    fold = NULL,
                    returnY=FALSE,
                    returnBiasDetails = FALSE,
                    trace=FALSE,
                    ...) {
  # run checks ------------------------------
  checked_data <- plmm_checks(X = X,
                              y = y,
                              std_needed = std_needed,
                              trace = trace,
                              ...)  
  
  # prep  ------------------------
  prep.args <- c(list(std_X = checked_data$std_X,
                      std_X_n = checked_data$std_X_n,
                      std_X_p = checked_data$std_X_p,
                      n = checked_data$n,
                      p = checked_data$p,
                      y = checked_data$y,
                      K = checked_data$K,
                      k = checked_data$k,
                      diag_K = checked_data$diag_K,
                      fbm_flag = checked_data$fbm_flag,
                      trace = trace,
                      ...)) # ... additional arguments to plmm_prep()

  prep <- do.call('plmm_prep', prep.args)

  # full model fit ----------------------------------
  fit.args <- c(list(prep = prep,
                     penalty = penalty,
                     std_X_details = checked_data$std_X_details,
                     eta_star = eta_star,
                     penalty.factor = checked_data$penalty.factor,
                     fbm_flag = checked_data$fbm_flag,
                     gamma = checked_data$gamma,
                     alpha = alpha,
                     nlambda = nlambda,
                     eps = eps,
                     max.iter = max.iter,
                     warn = warn,
                     convex = convex,
                     # TODO: figure out if/when to include dfmax... (for now, it is not used)
                     dfmax = checked_data$dfmax,
                     init = checked_data$init),
                list(...))
  
  fit <- do.call('plmm_fit', fit.args)
  
  if (is.null(col_names)){
    if (!is.null(checked_data$dat)) {
      col_names <- checked_data$dat$map$marker.ID
    }
  }
  fit_to_return <- plmm_format(fit = fit,
                               std_X_details = checked_data$std_X_details,
                               snp_names = col_names,
                               fbm_flag = checked_data$fbm_flag)
  
  # set up arguments for cv ---------------------------
  cv.args <- fit.args
  cv.args$warn <- FALSE
  cv.args$lambda <- fit$lambda 
  
  estimated_V <- NULL 
  if (type == 'blup') {
   estimated_V <- construct_variance(eta = fit$eta, K = prep$K)
  }
 
  # initialize objects to hold CV results 
  n <- length(fit$y)
  E <- Y <- matrix(NA, nrow=fit$std_X_n, ncol=length(fit$lambda))

  
  # set up folds for cross validation 
  if (!missing(seed)) {
    original_seed <- .GlobalEnv$.Random.seed
    on.exit(.GlobalEnv$.Random.seed <- original_seed)
    set.seed(seed)
  }
  
  sde <- sqrt(.Machine$double.eps)
  
  if(is.null(fold)) {
    if(trace){
      cat("'Fold' argument is either NULL or missing; assigning folds randomly (by default). 
          \nTo specify folds for each observation, supply a vector with fold assignments.")
    }
    fold <- sample(1:n %% nfolds)
    fold[fold==0] <- nfolds
  } else {
    nfolds <- max(fold)
  }

  
  # set up cluster if user-specified ---------------------------------------
  if (!missing(cluster)) {
    if (!inherits(cluster, "cluster")) stop("cluster is not of class 'cluster'; see ?makeCluster", call.=FALSE)
    parallel::clusterExport(cluster, c("X", "y", "K", "fold", "type", "cv.args", "estimated_V"), envir=environment())
    parallel::clusterCall(cluster, function() library(plmmr))
    fold.results <- parallel::parLapply(cl=cluster, X=1:max(fold), fun=cvf, X=X, y=y,
                                        fold=fold, type=type, cv.args=cv.args, 
                                        estimated_V = estimated_V)
  }
  
  # carry out CV -------------------------------------
  if (trace) cat("\nStarting cross validation\n")  
  # set up progress bar -- this can take a while
  if(trace){pb <- utils::txtProgressBar(min = 0, max = nfolds, style = 3)}
  for (i in 1:nfolds) {
    # case 1: user-specified cluster
    if (!missing(cluster)) {
      res <- fold.results[[i]] # refers to lines from above
      if (trace) {utils::setTxtProgressBar(pb, i)}
    } else {
      # case 2: cluster NOT user specified (this is the typical use case)
      res <- cvf(i = i,
                 fold = fold,
                 type = type,
                 cv.args = cv.args,
                 estimated_V = estimated_V)
      if (trace) {utils::setTxtProgressBar(pb, i)}

    }
    
    # update E and Y
    E[fold==i, 1:res$nl] <- res$loss

    if (!is.matrix(res$yhat)) {
      res$yhat <- as.matrix(res$yhat)
    }
    Y[fold==i, 1:res$nl] <- res$yhat
  }

  # post-process results -----------------------------------------
  # eliminate saturated lambda values, if any
  ind <- which(apply(is.finite(E), 2, all)) # index for lambda values to keep
  E <- E[, ind, drop=FALSE]
  Y <- Y[, ind]
  lambda <- fit$lambda[ind]

  # return min lambda idx
  cve <- apply(E, 2, mean)
  cvse <- apply(E, 2, stats::sd) / sqrt(n)
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
              lambda.min=lambda[min],
              min1se = min1se,
              lambda.1se = lambda[min1se],
              null.dev=mean(loss.plmm(checked_data$y, rep(mean(checked_data$y), n))))
  if (returnY) val$Y <- Y
  if (returnBiasDetails){
    val$Bias <- Bias
    val$Loss <- E
  }
  structure(val, class="cv.plmm")
}

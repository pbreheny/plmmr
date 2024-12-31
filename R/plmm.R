#' Fit a linear mixed model via non-convex penalized maximum likelihood.
#' @param design                  The first argument must be one of three things:
#'                                  (1) `plmm_design` object (as created by `create_design()`)
#'                                  (2) a string with the file path to a design object (the file path must end in '.rds')
#'                                  (3) a `matrix` or `data.frame` object representing the design matrix of interest
#' @param y                       Optional: In the case where `design` is a `matrix` or `data.frame`, the user must also supply
#'                                a numeric outcome vector as the `y` argument. In this case, `design` and `y` will be passed
#'                                internally to `create_design(X = design, y = y)`.
#' @param K                       Similarity matrix used to rotate the data. This should either be:
#'                                  (1) a known matrix that reflects the covariance of y,
#'                                  (2) an estimate (Default is \eqn{\frac{1}{p}(XX^T)}), or
#'                                  (3) a list with components 'd' and 'U', as returned by a previous `plmm()` model fit on the same data.
#' @param diag_K                  Logical: should K be a diagonal matrix? This would reflect observations that are unrelated, or that can be treated as unrelated. Defaults to FALSE.
#'                                Note: plmm() does not check to see if a matrix is diagonal. If you want to use a diagonal K matrix, you must set diag_K = TRUE.
#' @param eta_star                Optional argument to input a specific eta term rather than estimate it from the data. If K is a known covariance matrix that is full rank, this should be 1.
#' @param penalty                 The penalty to be applied to the model. Either "lasso" (the default), "SCAD", or "MCP".
#' @param init                    Initial values for coefficients. Default is 0 for all columns of X.
#' @param gamma                   The tuning parameter of the MCP/SCAD penalty (see details). Default is 3 for MCP and 3.7 for SCAD.
#' @param alpha                   Tuning parameter for the Mnet estimator which controls the relative contributions from the MCP/SCAD penalty and the ridge, or L2 penalty. alpha=1 is equivalent to MCP/SCAD penalty, while alpha=0 would be equivalent to ridge regression. However, alpha=0 is not supported; alpha may be arbitrarily small, but not exactly 0.
#' @param lambda_min              The smallest value for lambda, as a fraction of lambda.max. Default is .001 if the number of observations is larger than the number of covariates and .05 otherwise.
#' @param nlambda                 Length of the sequence of lambda. Default is 100.
#' @param lambda                  A user-specified sequence of lambda values. By default, a sequence of values of length nlambda is computed, equally spaced on the log scale.
#' @param eps                     Convergence threshold. The algorithm iterates until the RMSD for the change in linear predictors for each coefficient is less than eps. Default is \code{1e-4}.
#' @param max_iter                Maximum number of iterations (total across entire path). Default is 10000.
#' @param warn                    Return warning messages for failures to converge and model saturation? Default is TRUE.
#' @param trace                   If set to TRUE, inform the user of progress by announcing the beginning of each step of the modeling process. Default is FALSE.
#' @param save_rds                Optional: if a filepath and name *without* the '.rds' suffix is specified (e.g., `save_rds = "~/dir/my_results"`), then the model results are saved to the provided location (e.g., "~/dir/my_results.rds").
#'                                Accompanying the RDS file is a log file for documentation, e.g., "~/dir/my_results.log".
#'                                Defaults to NULL, which does not save any RDS or log files.
#' @param return_fit              Optional: a logical value indicating whether the fitted model should be returned as a `plmm` object in the current (assumed interactive) session.
#'                                Defaults to TRUE for in-memory data, and defaults to FALSE for filebacked data.
#' @param ...                     Additional optional arguments to `plmm_checks()`
#'
#' @returns A list which includes 19 items:
#'  * beta_vals: the matrix of estimated coefficients on the original scale. Rows are predictors, columns are values of `lambda`
#'  * std_Xbeta: a matrix of the linear predictors on the scale of the standardized design matrix. Rows are predictors, columns are values of `lambda`.
#'              **Note**: std_Xbeta will not include rows for the intercept or for constant features.
#'  * std_X_details: a list with 3 items: the center & scale values used to center/scale the data, and a vector ('ns') of the nonsingular columns
#'                  of the original data. Nonsingular columns cannot be standardized (by definition), and so were removed from analysis.
#'  * std_X: if design matrix is filebacked, the descriptor for the filebacked data is returned using \code{bigmemory::describe()}
#'  * y: the outcome vector used in model fitting.
#'  * p: the total number of columns in the design matrix (including singular columns).
#'  * plink_flag: a logical flag: did the data come from PLINK files?
#'  * lambda: a numeric vector of the lasso tuning parameter values used in model fitting.
#'  * eta: a number (double) between 0 and 1 representing the estimated proportion of the variance in the outcome attributable to population/correlation structure
#'  * std_Xbeta: the matrix resulting from the product of `stdrot_X` and the estimated coefficients on the ~rotated~ scale.
#'  * penalty: character string indicating the penalty with which the model was fit (e.g., 'MCP')
#'  * gamma: numeric value indicating the tuning parameter used for the SCAD or lasso penalties was used. Not relevant for lasso models.
#'  * alpha: numeric value indicating the elastic net tuning parameter.
#'  * loss: vector with the numeric values of the loss at each value of `lambda` (calculated on the ~rotated~ scale)
#'  * penalty_factor: vector of indicators corresponding to each predictor, where 1 = predictor was penalized.
#'  * ns_idx: vector with the indices of predictors which were non-singular features (i.e., features which had variation).
#'  * iter: numeric vector with the number of iterations needed in model fitting for each value of `lambda`
#'  * converged: vector of logical values indicating whether the model fitting converged at each value of `lambda`
#'  * K: a list with 2 elements, `s` and `U` ---
#'    * s: a vector of the eigenvalues of the relatedness matrix; see `relatedness_mat()` for details.
#'    * U: a matrix of the eigenvectors of the relatedness matrix
#' @export
#'
#' @examples
#' # using admix data
#' admix_design <- create_design(X = admix$X, y = admix$y)
#' fit_admix1 <- plmm(design = admix_design)
#' s1 <- summary(fit_admix1, idx = 50)
#' print(s1)
#' plot(fit_admix1)
#'
#' # Note: for examples with large data that are too big to fit in memory,
#' # see the article "PLINK files/file-backed matrices" on our website
#' # https://pbreheny.github.io/plmmr/articles/filebacking.html
#'
plmm <- function(design,
                 y = NULL,
                 K = NULL,
                 diag_K = NULL,
                 eta_star = NULL,
                 penalty = "lasso",
                 init = NULL,
                 gamma,
                 alpha = 1,
                 lambda_min, # passed to internal function setup_lambda()
                 nlambda = 100,
                 lambda,
                 eps = 1e-04,
                 max_iter = 10000,
                 warn = TRUE,
                 trace = FALSE,
                 save_rds = NULL,
                 return_fit = NULL,
                 ...) {

  # check filepaths for saving results, if requested  --------------------------

  if (!is.null(save_rds)) {
    save_rds <- check_for_file_extension(save_rds)
    # ^^ internally, we need to take off the extension from the file name

    # start the log file
    logfile <- create_log(outfile = ifelse(!is.null(save_rds),
                                           save_rds,
                                           "./plmm"))

  }


  # create a design if needed -------------
  if (inherits(design, "data.frame") | inherits(design, 'matrix')) {
    # error check: if 'design' is matrix/data.frame, user must specify 'y'
    if (is.null(y)) {
      stop("If you supply a matrix or data frame as 'design', you
                         must also specify a numeric vector as the outcome of your
                         model to the 'y' argument.")
    }

    design <- create_design_in_memory(X = design, y = y)


  } else {
    if (!is.null(y)) stop("If you are supplying a plmm_design object or filepath to
                          the 'design' argument, that design already has a 'y' --
                          please do not specify a 'y' argument here in plmm()")
  }

  # run checks ------------------------------
  checked_data <- plmm_checks(design,
                              K = K,
                              diag_K = diag_K,
                              eta_star = eta_star,
                              penalty = penalty,
                              init = init,
                              gamma = gamma,
                              alpha = alpha,
                              trace = trace,
                              ...)

  # set defaults for returning fit
  if (is.null(return_fit)) {
    if (checked_data$fbm_flag) {
      return_fit <- FALSE
    } else {
      return_fit <- TRUE
    }
  }

  # prep (SVD)-------------------------------------------------
  if(trace){cat("Input data passed all checks at ",
                pretty_time())}
  if (!is.null(save_rds)) {
  cat("Input data passed all checks at ",
      pretty_time(),
      "\n", file = logfile, append = TRUE)
  }

  the_prep <- plmm_prep(std_X = checked_data$std_X,
                        std_X_n = checked_data$std_X_n,
                        std_X_p = checked_data$std_X_p,
                        n = checked_data$n,
                        p = checked_data$p,
                        centered_y = checked_data$centered_y,
                        K = checked_data$K,
                        diag_K = checked_data$diag_K,
                        fbm_flag = checked_data$fbm_flag,
                        trace = trace)

  if (trace)(cat("Eigendecomposition finished at ",
                 pretty_time()))

  if (!is.null(save_rds)) {
    cat("Eigendecomposition finished at ",
        pretty_time(),
        "\n", file = logfile, append = TRUE)
  }


  # rotate & fit -------------------------------------------------------------
  the_fit <- plmm_fit(prep = the_prep,
                      y = checked_data$y,
                      std_X_details = checked_data$std_X_details,
                      eta_star = eta_star,
                      penalty_factor = checked_data$penalty_factor,
                      fbm_flag = checked_data$fbm_flag,
                      penalty = checked_data$penalty,
                      gamma = checked_data$gamma,
                      alpha = alpha,
                      lambda_min = lambda_min,
                      nlambda = nlambda,
                      lambda = lambda,
                      eps = eps,
                      max_iter = max_iter,
                      warn = warn)

  if (trace) cat("Beta values are estimated -- almost done!\n")

  # format results ---------------------------------------------------
  if(trace){cat("Formatting results (backtransforming coefs. to original scale).\n")}
  the_final_product <- plmm_format(fit = the_fit,
                                   p = checked_data$p,
                                   std_X_details = checked_data$std_X_details,
                                   fbm_flag = checked_data$fbm_flag,
                                   plink_flag = checked_data$plink_flag)

  if (trace)(cat("Model ready at ",
                 pretty_time()))

  if (!is.null(save_rds)){
    cat("Model ready at ",
        pretty_time(),
        file = logfile, append = TRUE)
  }

  # handle output
  if (!is.null(save_rds)){
    # save all output in one file (default); *not* including std_X
    saveRDS(the_final_product[c(1:3, 5:19)],
            paste0(save_rds, ".rds"))
    cat("Results saved to:", paste0(save_rds, ".rds"), "at",
        pretty_time(),
        file = logfile, append = TRUE)
  }


  # create a failsafe -- if save_rds is NULL, make sure return_fit = TRUE
  if (is.null(save_rds) & !return_fit){
    cat("You accidentally left save_rds = NULL and return_fit = FALSE;
        to prevent you from losing your work, plmm() is returning the output as if return_fit = TRUE")

    return_fit <- TRUE
  }

  if (return_fit){
    return(the_final_product)
  }

}

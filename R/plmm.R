#' Fit a linear mixed model via non-convex penalized maximum likelihood.
#' @param design                  A `plmm_design` object (as created by `create_design()`) or a string with the file path to a design object (the file path must end in '.rds').
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
#' @param dfmax                   (**Future idea; not yet incorporated**): Upper bound for the number of nonzero coefficients. Default is no upper bound. However, for large data sets, computational burden may be heavy for models with a large number of nonzero coefficients.
#' @param lambda_min              The smallest value for lambda, as a fraction of lambda.max. Default is .001 if the number of observations is larger than the number of covariates and .05 otherwise.
#' @param nlambda                 Length of the sequence of lambda. Default is 100.
#' @param lambda                  A user-specified sequence of lambda values. By default, a sequence of values of length nlambda is computed, equally spaced on the log scale.
#' @param eps                     Convergence threshold. The algorithm iterates until the RMSD for the change in linear predictors for each coefficient is less than eps. Default is \code{1e-4}.
#' @param max_iter                Maximum number of iterations (total across entire path). Default is 10000.
#' @param convex                  (**Future idea; not yet incorporated**): Calculate index for which objective function ceases to be locally convex? Default is TRUE.
#' @param warn                    Return warning messages for failures to converge and model saturation? Default is TRUE.
#' @param trace                   If set to TRUE, inform the user of progress by announcing the beginning of each step of the modeling process. Default is FALSE.
#' @param save_rds                Optional: if a filepath and name *without* the '.rds' suffix is specified (e.g., `save_rds = "~/dir/my_results"`), then the model results are saved to the provided location (e.g., "~/dir/my_results.rds").
#'                                Defaults to NULL, which does not save the result.
#' @param return_fit              Optional: a logical value indicating whether the fitted model should be returned as a `plmm` object in the current (assumed interactive) session.
#'                                Defaults to TRUE for in-memory data, and defaults to FALSE for filebacked data.
#' @param compact_save            Optional: if TRUE, three separate .rds files will saved: one with the 'beta_vals', one with 'K', and one with everything else (see below).
#'                                Defaults to FALSE. **Note**: you must specify `save_rds` for this argument to be called.
#' @param ...                     Additional optional arguments to `plmm_checks()`
#'
#' @returns A list which includes:
#'  * beta_vals: the matrix of estimated coefficients on the original scale. Rows are predictors, columns are values of `lambda`
#'  * rotated_scale_beta_vals: the matrix of estimated coefficients on the ~rotated~ scale. This is the scale on which the model was fit.
#'  * lambda: a numeric vector of the lasso tuning parameter values used in model fitting.
#'  * eta: a number (double) between 0 and 1 representing the estimated proportion of the variance in the outcome attributable to population/correlation structure
#'  * linear_predictors: the matrix resulting from the product of `stdrot_X` and the estimated coefficients on the ~rotated~ scale.
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
#' admix_design <- create_design(X = admix$X, outcome_col = admix$y)
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
                 K = NULL,
                 diag_K = NULL,
                 eta_star = NULL,
                 penalty = "lasso",
                 init = NULL,
                 gamma,
                 alpha = 1,
                 dfmax = NULL,
                 lambda_min, # passed to internal function setup_lambda()
                 nlambda = 100,
                 lambda,
                 eps = 1e-04,
                 max_iter = 10000,
                 convex = TRUE,
                 warn = TRUE,
                 trace = FALSE,
                 save_rds = NULL,
                 compact_save = FALSE,
                 return_fit = NULL,
                 ...) {

  # check filepaths for saving results ------------------------------

  if (compact_save & is.null(save_rds)) {
    stop("You have set 'compact_save = TRUE', but no argument was supplied to 'save_rds'.
          \nPlease specify a filepath (as a string) to 'save_rds'")
  }

  if (!is.null(save_rds)) {
    save_rds <- check_for_file_extension(save_rds)
    # ^^ internally, we need to take off the extension from the file name
  }

  # start the log -----------------------
  logfile <- create_log(outfile = ifelse(!is.null(save_rds),
                                         save_rds,
                                         "./plmm"))

  # run checks ------------------------------
  checked_data <- plmm_checks(design,
                              K = K,
                              diag_K = diag_K,
                              eta_star = eta_star,
                              penalty = penalty,
                              init = init,
                              dfmax = dfmax,
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

  cat("Input data passed all checks at ",
      pretty_time(),
      "\n", file = logfile, append = TRUE)

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

  cat("Eigendecomposition finished at ",
      pretty_time(),
      "\n", file = logfile, append = TRUE)

    # rotate & fit -------------------------------------------------------------
  if (is.null(dfmax)) dfmax <- checked_data$p + 1
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
                                   fbm_flag = checked_data$fbm_flag)

  if (trace)(cat("Model ready at ",
                 pretty_time()))
  cat("Model ready at ",
      pretty_time(),
      file = logfile, append = TRUE)

  # handle output
  if (!is.null(save_rds)){
    if (compact_save) {
      # save output across multiple files
      saveRDS(the_final_product$beta_vals, paste0(save_rds, "_coefficients.rds"))
      cat("Coefficients (estimated beta values) saved to:", paste0(save_rds, "_coefficients.rds"), "at",
          pretty_time(),
          file = logfile, append = TRUE)

      saveRDS(the_final_product$K, paste0(save_rds, "_K.rds"))
      cat("K (eigendecomposition) saved to:", paste0(save_rds, "_K.rds"), "at",
          pretty_time(),
          file = logfile, append = TRUE)

      saveRDS(the_final_product$linear_predictors, paste0(save_rds, "_linear_predictors.rds"))
      cat("Linear predictors (on rotated scale) saved to:", paste0(save_rds, "_linear_predictors.rds"), "at",
          pretty_time(),
          file = logfile, append = TRUE)

      saveRDS(the_final_product[c(2:3, 5:12)], paste0(save_rds, "_details.rds"))
      cat("All other results (loss, # of iterations, ...) saved to:", paste0(save_rds, "_details.rds"), "at",
          pretty_time(),
          file = logfile, append = TRUE)

    } else {
      # save all output in one file (default)
      saveRDS(the_final_product, paste0(save_rds, ".rds"))
      cat("Results saved to:", paste0(save_rds, ".rds"), "at",
          pretty_time(),
          file = logfile, append = TRUE)
    }


  }

  if (is.null(save_rds) & !return_fit){
    cat("You accidentally left save_rds NULL while setting return_fit = FALSE;
        to prevent you from losing your work, I am saving the output as plmm_results.rds
        in your current working directory (current folder).\n
        Next time, make sure to specify your own filepath to the save_rds argument.\n",
        file = logfile, append = TRUE)

    rdsfile <- paste0(getwd(),"/plmm_results.rds")
    saveRDS(the_final_product, rdsfile)
    cat("Results saved to:", rdsfile, file = logfile, append = TRUE)
  }

  if (return_fit){
    return(the_final_product)
  }

}

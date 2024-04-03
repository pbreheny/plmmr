#' Fit a linear mixed model with non-convex regularization
#'
#' This function allows you to fit a linear mixed model via non-convex penalized maximum likelihood.
#' NB: this function is simply a wrapper for plmm_prep -> plmm_fit -> plmm_format
#' @param X Design matrix object or a string with the file path to a design matrix. If a string, string will be passed to `get_data()`. 
#' * Note: X may include clinical covariates and other non-SNP data, but no missing values are allowed.
#' @param fbm Logical: should X be treated as filebacked? Relevant only when X is a string to be passed to `get_data()`. Defaults to NULL, using the default setttings of `get_data()` to determine whether X should be stored in memory.
#' @param std_needed Logical: does the supplied X need to be standardized? Defaults to FALSE, since `process_plink()` standardizes the design matrix by default. 
#' By default, X will be standardized internally. For data processed from PLINK files, standardization happens in `process_plink()`. For data supplied as a matrix, standardization happens here in `plmm()`. If you know your data are already standardized, set `std_needed = FALSE` -- this would be an atypical case. **Note**: failing to standardize data will lead to incorrect analyses. 
#' @param col_names Optional vector of column names for design matrix. Defaults to NULL.
#' @param y Continuous outcome vector. Defaults to NULL, assuming that the outcome is the 6th column in the .fam PLINK file data. Can also be a user-supplied numeric vector. 
#' @param k An integer specifying the number of singular values to be used in the approximation of the rotated design matrix. This argument is passed to `RSpectra::svds()`. Defaults to `min(n, p) - 1`, where n and p are the dimensions of the _standardized_ design matrix.
#' @param K Similarity matrix used to rotate the data. This should either be (1) a known matrix that reflects the covariance of y, (2) an estimate (Default is \eqn{\frac{1}{p}(XX^T)}), or (3) a list with components 'd' and 'u', as returned by choose_k().
#' @param diag_K Logical: should K be a diagonal matrix? This would reflect observations that are unrelated, or that can be treated as unrelated. Defaults to FALSE. 
#'  Note: plmm() does not check to see if a matrix is diagonal. If you want to use a diagonal K matrix, you must set diag_K = TRUE.
#' @param eta_star Optional argument to input a specific eta term rather than estimate it from the data. If K is a known covariance matrix that is full rank, this should be 1.
#' @param penalty The penalty to be applied to the model. Either "MCP" (the default), "SCAD", or "lasso".
#' @param gamma The tuning parameter of the MCP/SCAD penalty (see details). Default is 3 for MCP and 3.7 for Spenncath.
#' @param alpha Tuning parameter for the Mnet estimator which controls the relative contributions from the MCP/Spenncath penalty and the ridge, or L2 penalty. alpha=1 is equivalent to MCP/Spenncath penalty, while alpha=0 would be equivalent to ridge regression. However, alpha=0 is not supported; alpha may be arbitrarily small, but not exactly 0.
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
#' @param trace If set to TRUE, inform the user of progress by announcing the beginning of each step of the modeling process. Default is FALSE.
#' @returns A list which includes: 
#'  * `beta_vals`: the matrix of estimated coefficients on the original scale. Rows are predictors, columns are values of `lambda`
#'  * `rotated_scale_beta_vals`: the matrix of estimated coefficients on the ~rotated~ scale. This is the scale on which the model was fit. 
#'  * `lambda`: a numeric vector of the lasso tuning parameter values used in model fitting. 
#'  * `eta`: a number (double) between 0 and 1 representing the estimated proportion of the variance in the outcome attributable to population/correlation structure.
#'  * `s`: a vectof of the eigenvalues of relatedness matrix `K`; see `relatedness_mat()` for details.
#'  * `U`: a matrix of the eigenvalues of relatedness matrix `K`
#'  * `rot_y`: the vector of outcome values on the rotated scale. This is the scale on which the model was fit. 
#'  * `linear.predictors`: the matrix resulting from the product of `stdrot_X` and the estimated coefficients on the ~rotated~ scale.
#'  * `penalty`: character string indicating the penalty with which the model was fit (e.g., 'MCP')
#'  * `gamma`: numeric value indicating the tuning parameter used for the SCAD or lasso penalties was used. Not relevant for lasso models.
#'  * `alpha`: numeric value indicating the elastic net tuning parameter. 
#'  * `convex.min`: NULL (This is an option we will add in the future!)
#'  * `loss`: vector with the numeric values of the loss at each value of `lambda` (calculated on the ~rotated~ scale)
#'  * `penalty.factor`: vector of indicators corresponding to each predictor, where 1 = predictor was penalized. 
#'  * `ns_idx`: vector with the indicies of predictors which were constant features (i.e., had no variation).
#'  * `p`: the number of features 
#'  * `n`: the number of observations (instances)
#'  * `iter`: numeric vector with the number of iterations needed in model fitting for each value of `lambda`
#'  * `converged`: vector of logical values indicating whether the model fitting converged at each value of `lambda`

#' @export
#' 
#' @examples 
#' # using admix data 
#' fit_admix1 <- plmm(X = admix$X, y = admix$y, std_needed = TRUE)
#' s1 <- summary(fit_admix1, idx = 99)
#' print(s1)
#' plot(fit_admix1)
#' 
#' # an example with p > n:
#' fit_admix2 <- plmm(X = admix$X[1:50, ], y = admix$y[1:50])
#' s2 <- summary(fit_admix2, idx = 99)
#' print(s2)
#' plot(fit_admix2) # notice: the default penalty is MCP
#' 
#' # now use PLINK data files
#' \dontrun{
#' 
#' # file-backed example
#' plmm(X = "~/tmp_files/penncath_lite", # adjust this line to fit current machine
#'  fbm = TRUE, trace = TRUE) -> fb_fit
#' 
#' penncath_mid <- process_plink(prefix = "penncath_mid", dataDir = plink_example(path="penncath_mid.fam", parent=T))
#' penncath_clinical <- read.csv(plink_example(path="penncath_clinical.csv"))
#' # for the sake of illustration, I use a simple mean imputation for the outcome 
#' penncath_clinical$hdl_impute <- ifelse(is.na(penncath_clinical$hdl), mean(penncath_clinical$hdl, na.rm = T), penncath_clinical$hdl)
#' 
#' # fit with no 'k' specified
#' fit_plink1 <- plmm(X = penncath_mid$X, y = penncath_clinical$hdl_impute, trace = TRUE)
#' summary(fit_plink1, idx = 5)
#' # Runs in ~219 seconds (3.65 mins) on my 2015 MacBook Pro
#' 
#' # fit with 'k = 5' specified (so using RSpectra::svds())
#' fit_plink2 <- plmm(X = penncath_mid$X, y = penncath_clinical$hdl_impute, k = 5, trace = TRUE)
#' # Runs in ~44 seconds on my 2015 MacBook Pro
#' summary(fit_plink2, idx = 5);summary(fit_plink2, idx = 95)
#' 
#' 
#' # case where X is an FBM
#' lite <- get_data("../temp_files/penncath_lite", fbm = TRUE)
#' clinical <- read.csv("../temp_files/penncath_clinical.csv")
#' hdl <- ifelse(is.na(clinical$hdl), mean(clinical$hdl, na.rm = TRUE), clinical$hdl)
#' fit_fbm <- plmm(X = lite, y = hdl, k = 1200)
#' }
#' 
plmm <- function(X,
                 fbm = NULL,
                 std_needed = NULL,
                 col_names = NULL,
                 y = NULL,
                 k = NULL, 
                 K = NULL,
                 diag_K = NULL,
                 eta_star = NULL,
                 penalty = c("MCP", "SCAD", "lasso"), # TODO: think about making lasso default
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
                 penalty.factor = NULL,
                 init = NULL,
                 trace = FALSE) {
  
# run checks ------------------------------

  checked_data <- plmm_checks(X,
                              fbm = fbm,
                              std_needed = std_needed,
                              col_names = col_names,
                              y = y,
                              k = k, 
                              K = K,
                              diag_K = diag_K,
                              eta_star = eta_star,
                              penalty = penalty, # TODO: think about making lasso default
                              penalty.factor = penalty.factor,
                              init = init,
                              dfmax = dfmax,
                              gamma = gamma,
                              alpha = alpha,
                              trace = trace)  
  
  # prep (SVD)-------------------------------------------------
  if(trace){cat("\nInput data passed all checks.")}
    the_prep <- plmm_prep(std_X = checked_data$std_X,
                          std_X_n = checked_data$std_X_n,
                          std_X_p = checked_data$std_X_p,
                          n = checked_data$n,
                          p = checked_data$p,
                          y = checked_data$y,
                          K = checked_data$K,
                          k = checked_data$k,
                          diag_K = checked_data$diag_K,
                          fbm_flag = checked_data$fbm_flag,
                          trace = trace)
  
  # rotate & fit -------------------------------------------------------------
  if(trace){cat("\nBeginning model fitting.\n")}
  
  the_fit <- plmm_fit(prep = the_prep,
                      std_X_details = checked_data$std_X_details,
                      eta_star = eta_star,
                      penalty.factor = checked_data$penalty.factor,
                      fbm_flag = checked_data$fbm_flag,
                      penalty = checked_data$penalty,
                      gamma = checked_data$gamma,
                      alpha = alpha,
                      lambda.min = lambda.min,
                      nlambda = nlambda,
                      lambda = lambda,
                      eps = eps,
                      max.iter = max.iter,
                      warn = warn,
                      convex = convex,
                      # TODO: figure out if/when to include dfmax... (for now, it is not used)
                      dfmax = checked_data$dfmax,
                      init = checked_data$init)
  
  if(trace){cat("\nBeta values are estimated -- almost done!")}
  # format results ---------------------------------------------------
  if(trace){cat("\nFormatting results (backtransforming coefs. to original scale).\n")}

  if (is.null(col_names)){
    if (!is.null(checked_data$dat)) {
      col_names <- checked_data$dat$map$marker.ID
    }
  }
    the_final_product <- plmm_format(fit = the_fit,
                                     std_X_details = checked_data$std_X_details,
                                     snp_names = col_names,
                                     fbm_flag = checked_data$fbm_flag)
    
  
  
  return(the_final_product)
  
  
}

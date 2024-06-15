#' PLMM format: a function to format the output of a model constructed with \code{plmm_fit}
#'
#' @param fit A list of parameters describing the output of a model constructed with \code{plmm_fit}
#' @param p The number of features in the original data (including constant features)
#' @param std_X_details A list with 3 items:
#'  * 'center': the centering values for the columns of `X`
#'  * 'scale': the scaling values for the non-singular columns of `X`
#'  * 'ns': indicies of nonsingular columns in `std_X`
#' @param fbm_flag Logical: is the corresponding design matrix filebacked? Passed from `plmm()`.
#' @param feature_names A vector of names for features, passed internally if such names are included with the data `X` passed to `plmm()`
#' @param non_genomic Optional vector specifying which columns of the design matrix represent features that are *not* genomic, as these features are excluded from the empirical estimation of genomic relatedness.
#' For cases where X is a filepath to an object created by `process_plink()`, this is handled automatically via the arguments to `process_plink()`.
#' For all other cases, 'non_genomic' defaults to NULL (meaning `plmm()` will assume that all columns of `X` represent genomic features).
#'
#' @returns A list with the components:
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
#'  * `loss`: vector with the numeric values of the loss at each value of `lambda` (calculated on the ~rotated~ scale)
#'  * `penalty.factor`: vector of indicators corresponding to each predictor, where 1 = predictor was penalized.
#'  * `ns_idx`: vector with the indicies of predictors which were constant features (i.e., had no variation).
#'  * `p`: the number of features
#'  * `n`: the number of observations (instances)
#'  * `iter`: numeric vector with the number of iterations needed in model fitting for each value of `lambda`
#'  * `converged`: vector of logical values indicating whether the model fitting converged at each value of `lambda`
#'
#' @keywords internal

plmm_format <- function(fit, p, std_X_details, fbm_flag, feature_names = NULL, non_genomic = NULL){

  # get beta values back in original scale; reverse the PRE-ROTATION standardization
  og_scale_beta <- untransform(
    std_scale_beta = fit$std_scale_beta,
    p = p,
    std_X_details = std_X_details,
    fbm_flag = fbm_flag,
    non_genomic = non_genomic)

  # give the matrix of beta_values readable names
  # SNPs (or covariates) on the rows, lambda values on the columns

  if (is.null(feature_names)) {
    feature_names <- paste("Var", 1:(p + length(non_genomic)), sep="")
  }

  varnames <- c("(Intercept)", feature_names) # add intercept label
  dimnames(og_scale_beta) <- list(varnames, lam_names(fit$lambda))
  colnames(fit$linear.predictors) <- lam_names(fit$lambda)


  # output
  structure(list(
    beta_vals = og_scale_beta,
    # rotated_scale_beta_vals = fit$std_scale_beta,
    lambda = fit$lambda,
    eta = fit$eta,
    # rot_y = fit$rot_y,
    linear.predictors = fit$linear.predictors,
    penalty = fit$penalty,
    gamma = fit$gamma,
    alpha = fit$alpha,
    loss = fit$loss,
    penalty.factor = fit$penalty.factor,
    ns_idx = c(1, 1 + fit$ns), # PAY ATTENTION HERE!
    # p = fit$p,
    # n = fit$n,
    # std_X_n = fit$std_X_n,
    # std_X_p = fit$std_X_p,
    iter = fit$iter,
    converged = fit$converged,
    K = list(s = fit$s, U = fit$U)),
    class = "plmm")
}

#' PLMM format: a function to format the output of a model constructed with `plmm_fit()`
#'
#' @param fit A list of parameters describing the output of a model constructed with `plmm_fit()`
#' @param p The number of features in the original data (including constant features)
#' @param std_X_details A list with 3 items:
#'  * `center`: the centering values for the columns of `X`
#'  * `scale`: the scaling values for the non-singular columns of `X`
#'  * `ns`: indices of nonsingular columns in `std_X`
#' @param fbm_flag Logical: is the corresponding design matrix filebacked? Passed from `plmm()`.
#' @param plink_flag Logical: did these data come from PLINK files?
#'                    **Note**: This flag matters because of how non-genomic features
#'                    are handled for PLINK files -- in data from PLINK files,
#'                    unpenalized columns are *not* counted in the `p` argument. For delimited
#'                    files, `p` does include unpenalized columns. This difference has
#'                    implications for how the `untransform()` function determines the
#'                    appropriate dimensions for the estimated coefficient matrix it returns.
#'
#' @return A list with 18 components:
#'  * `beta_vals`: the matrix of estimated coefficients on the original scale. Rows are predictors, columns are values of `lambda`
#'  * `std_Xbeta`: A matrix of the linear predictors on the scale of the standardized design matrix. Rows are predictors, columns are values of `lambda`.
#'                **Note**: `std_Xbeta` will not include rows for the intercept or for constant features.
#'  * `std_X_details`: A list with 9 items:
#'    - `center`: The center values used to center the columns of the design matrix
#'    - `scale`: The scaling values used to scale the columns of the design matrix
#'    - `ns`: An integer vector of the nonsingular columns of the original data
#'    - `unpen`: An integer vector of indices of the unpenalized features, if any were specified in the design
#'    - `unpen_colnames`: A character vector of the column names of any unpenalized features.
#'    - `X_colnames`: A character vector with the column names of all features in the original design matrix
#'    - `X_rownames`: A character vector with the row names of all features in the original design matrix; if none were provided, these are named 'row1', 'row2', etc.
#'    - `std_X_colnames`: A subset of `X_colnames` representing only nonsingular columns (i.e., the columns indexed by `ns`)
#'    - `std_X_rownames`: A subset of `X_rownames` representing rows that passed QC filtering & and are represented in both the genotype and phenotype data sets (this only applies to PLINK data)
#'  * `y`: The original outcome vector.
#'  * `p`: The total number of columns in the design matrix (including any singular columns, excluding the intercept).
#'  * `plink_flag`: Logical - did the data come from PLINK files?
#'  * `lambda`: a numeric vector of the lasso tuning parameter values used in model fitting.
#'  * `eta`: a number (double) between 0 and 1 representing the estimated proportion of the variance in the outcome attributable to population/correlation structure.
#'  * `penalty`: character string indicating the penalty with which the model was fit (e.g., 'MCP')
#'  * `gamma`: numeric value indicating the tuning parameter used for the SCAD or lasso penalties was used. Not relevant for lasso models.
#'  * `alpha`: numeric value indicating the elastic net tuning parameter.
#'  * `loss`: vector with the numeric values of the loss at each value of `lambda` (calculated on the ~rotated~ scale)
#'  * `penalty_factor`: vector of indicators corresponding to each predictor, where 1 = predictor was penalized.
#'  * `ns_idx`: vector with the indices of predictors which were nonsingular features (i.e., had variation).
#'  * `iter`: numeric vector with the number of iterations needed in model fitting for each value of `lambda`
#'  * `converged`: vector of logical values indicating whether the model fitting converged at each value of `lambda`
#'  * `K`: a list with 2 elements, `s` and `U` ---
#'    - `s`: a vector of the non-zero eigenvalues of the relatedness matrix K (note: K is the kinship matrix for genetic/genomic data; see the article on notation for details)
#'    - `U`: a matrix of the eigenvectors of K associated with `s`
#'  * `std_X`: If design matrix is filebacked, the descriptor for the filebacked data is returned using `bigmemory::describe()`.
#'
#' @keywords internal
#'
plmm_format <- function(fit, p, std_X_details, fbm_flag, plink_flag) {
  # get beta values back in original scale; reverse the PRE-ROTATION standardization
  og_scale_beta <- untransform(
    std_scale_beta = fit$std_scale_beta,
    p = p,
    std_X_details = std_X_details,
    fbm_flag = fbm_flag,
    plink_flag = plink_flag)

  # give the matrix of beta_values readable names
  # features on the rows, lambda values on the columns
  colnames(og_scale_beta) <- colnames(fit$std_Xbeta) <- lam_names(fit$lambda)

  # output (18 items)
  out <- list(
    beta_vals = og_scale_beta,
    std_Xbeta = fit$std_Xbeta,
    std_X_details = std_X_details,
    y = fit$y,
    p = p, # (need to hold onto the total number of features in the original data)
    plink_flag = plink_flag,
    lambda = fit$lambda,
    eta = fit$eta,
    penalty = fit$penalty,
    gamma = fit$gamma,
    alpha = fit$alpha,
    loss = fit$loss,
    penalty_factor = fit$penalty_factor,
    ns_idx = c(1, 1 + std_X_details$ns), # NOTE: this indexing is *very* important
    iter = fit$iter,
    converged = fit$converged,
    K = list(s = fit$s, U = fit$U))

  if (fbm_flag) {
    out$std_X <- fit$std_X
  }

  structure(out, class = "plmm")
}

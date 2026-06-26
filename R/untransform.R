#' Untransform coefficient values back to the original scale
#'
#' This function unwinds the initial standardization of the data to obtain
#' coefficient values on their original scale. It is called by `plmm_format()`.
#'
#' @param std_scale_beta The estimated coefficients on the standardized scale
#' @param p  The number of columns in the original design matrix
#' @param std_X_details A list with 3 elements describing the standardized design matrix BEFORE rotation; this should have elements `scale`, `center`, and `ns`
#' @param fbm_flag Logical: is the corresponding design matrix filebacked?
#' @param plink_flag Logical: did these data come from PLINK files?
#'                    **Note**: This flag matters because of how non-genomic features
#'                    are handled for PLINK files -- in data from PLINK files,
#'                    unpenalized columns are *not* counted in the `p` argument. For delimited
#'                    files, `p` does include unpenalized columns. This difference has
#'                    implications for how the `untransform()` function determines the
#'                    appropriate dimensions for the estimated coefficient matrix it returns.
#' @param use_names Logical: should names be added? Defaults to TRUE. Set to FALSE inside of `cvf()` helper, as `ns` will vary within CV folds.
#'
#' @return a matrix of estimated coefficients, `untransformed_beta`, that is on the scale of the original data.
#'
#' @keywords internal
#'
untransform <- function(std_scale_beta, p, std_X_details, fbm_flag, plink_flag, use_names = TRUE) {
  if (is.null(std_X_details$X_colnames)) {
    use_names <- FALSE
  }

  ns <- std_X_details$ns
  scale <- std_X_details$scale
  center <- std_X_details$center

  a <- std_scale_beta[1, , drop = FALSE]
  b <- std_scale_beta[-1, , drop = FALSE]

  nrow_out <- p + 1L + if (plink_flag) length(std_X_details$unpen) else 0L

  if (fbm_flag) {
    untransformed_beta <- Matrix::Matrix(
      0,
      nrow = nrow_out,
      ncol = ncol(std_scale_beta),
      sparse = TRUE
    )
  } else {
    untransformed_beta <- matrix(0, nrow = nrow_out, ncol = ncol(std_scale_beta))
  }

  # When scale has more entries than ns (e.g. CV folds where some columns became
  # constant), index scale/center to ns. When b has more rows than ns, subset b.
  s <- if (length(scale) == length(ns)) scale else scale[ns]
  ctr <- if (length(scale) == length(ns)) center else center[ns]
  b2 <- sweep(if (nrow(b) == length(ns)) b else b[ns, , drop = FALSE], 1, s, "/")

  untransformed_beta[ns + 1L, ] <- b2
  untransformed_beta[1L, ] <- a - crossprod(ctr, b2)

  if (use_names) {
    unpen_names <- if (plink_flag) std_X_details$unpen_colnames else NULL
    rownames(untransformed_beta) <- c("(Intercept)", unpen_names, std_X_details$X_colnames)
  }

  untransformed_beta
}

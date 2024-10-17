#' Untransform coefficient values back to the original scale
#'
#' This function unwinds the initial standardization of the data to obtain
#' coefficient values on their original scale. It is called by plmm_format().
#'
#' @param std_scale_beta The estimated coefficients on the standardized scale
#' @param p  The number of columns in the original design matrix
#' @param std_X_details A list with 3 elements describing the standardized design matrix BEFORE rotation; this should have elements 'scale', 'center', and 'ns'
#' @param fbm_flag Logical: is the corresponding design matrix filebacked?
#' @param use_names Logical: should names be added? Defaults to TRUE. Set to FALSE inside of `cvf()` helper, as 'ns' will vary within CV folds.
#' @keywords internal
#'
#' @returns a matrix of estimated coeffcients, 'beta_vals', that is on the scale of the original data.


untransform <- function(std_scale_beta, p, std_X_details, fbm_flag, use_names = TRUE) {
  if (is.null(std_X_details$X_colnames)) use_names <- FALSE

  if (fbm_flag) {
    untransform_filebacked(std_scale_beta = std_scale_beta,
                           p = p,
                           std_X_details = std_X_details,
                           use_names = use_names)
  } else {
    untransform_in_memory(
      std_scale_beta = std_scale_beta,
      p = p,
      std_X_details = std_X_details,
      use_names = use_names)
  }
}





#' An internal function to adjust the dimensions of a matrix of estimated coefficients
#'  returned by plmm_fit().
#'
#'  This function is designed for use in BLUP prediction.
#'  The objective here is to get a matrix of estimated beta coefficients that
#'  are on the standardized scale, but have the dimension of the original/training data.
#'  We do this by adding rows of 0s to the std_scale_beta matrix corresponding to
#'  the singular features of X.
#'
#' @param std_scale_beta  A matrix of estimated beta coefficients on the scale of the standardized original/training data
#'                        Note: the rows of this matrix represent the *nonsingular* columns of the design matrix
#' @param p               The number of columns in the original/training design matrix
#' @param std_X_details   A list with 3 elements describing the standardized design matrix BEFORE rotation; this should have elements 'scale', 'center', and 'ns'
#' @param fbm_flag        Logical: was this model fit filebacked?
#' @param plink_flag Logical: did these data come from PLINK files?
#'                    **Note**: This flag matters because of how non-genomic features
#'                    are handled for PLINK files -- in data from PLINK files,
#'                    unpenalized columns are *not* counted in the `p` argument. For delimited
#'                    files, `p` does include unpenalized columns. This difference has
#'                    implications for how the `untransform()` function determines the
#'                    appropriate dimensions for the estimated coefficient matrix it returns.
#' @returns std_scale_b_og_dim: a matrix of estimated beta coefs. that is still on the scale of std_X, but has the dimension of X
#' @keywords internal
#'
adjust_beta_dimension <- function(std_scale_beta, p, std_X_details,
                                  fbm_flag, plink_flag){

  ## Note: for in-memory matrix data & delimited data, unpenalized columns are included in 'p',
  ##     the number of columns in X;
  ##     This differs from the PLINK data, where unpenalized columns are
  ##     *appended* to the columns of 'X'.
  if (fbm_flag){
    a <- std_scale_beta[1, , drop = FALSE] # this is the intercept
    b <- std_scale_beta[-1, , drop=FALSE]
    if (plink_flag) {
      std_scale_b_og_dim <- Matrix::Matrix(0,
                                           nrow = (p + length(std_X_details$unpen) + 1), # + 1 is for the intercept; see note below
                                           ncol = ncol(std_scale_beta), # this is the number of lambda values
                                           sparse = TRUE)
    } else {
      std_scale_b_og_dim <- Matrix::Matrix(0,
                                           nrow = (p + 1), # + 1 is for the intercept; see note below
                                           ncol = ncol(std_scale_beta), # this is the number of lambda values
                                           sparse = TRUE)
    }

    if (nrow(b) != (nrow(std_scale_b_og_dim)-1)) {
      std_scale_b_og_dim[std_X_details$ns+1,] <- b
    } else {
      std_scale_b_og_dim[-1,] <- b
    }
    std_scale_b_og_dim[1,] <- a
  } else {

    # initialize beta with zeros; nrow = # of predictors, ncol = # of lambda values
    # this will create columns of zeros for betas corresponding to singular columns
    std_scale_b_og_dim <- matrix(0,
                                 nrow = (p + 1), # + 1 is for the intercept
                                 ncol = ncol(std_scale_beta))
    a <- std_scale_beta[1, , drop = FALSE] # this is the intercept
    b <- std_scale_beta[-1, , drop=FALSE]
    std_scale_b_og_dim[std_X_details$ns+1,] <- b
    std_scale_b_og_dim[1,] <- a
  }
  return(std_scale_b_og_dim)
}
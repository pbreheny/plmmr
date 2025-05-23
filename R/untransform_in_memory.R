##' Untransform coefficient values back to the original scale *In memory*
#'
#' This function unwinds the initial standardization of the data to obtain
#' coefficient values on their original scale. It is called by plmm_format().
#'
#' @param std_scale_beta The estimated coefficients on the standardized scale
#' @param p  The number of columns in the original design matrix
#' @param std_X_details A list with 3 elements describing the standardized design matrix BEFORE rotation; this should have elements 'scale', 'center', and 'ns'
#' @param use_names Logical: should names be added? Defaults to TRUE. Set to FALSE inside of `cvf()` helper, as 'ns' will vary within CV folds.
#' @keywords internal
#'
#' @returns a matrix of estimated coeffcients, 'beta_vals', that is on the scale of the original data.

untransform_in_memory <- function(std_scale_beta, p, std_X_details, use_names = TRUE) {

  # goal: reverse the PRE-ROTATION standardization #
  # partition the values from Step 1 into intercept and non-intercept parts
  a <- std_scale_beta[1, , drop = FALSE] # this is the intercept
  b <- std_scale_beta[-1, , drop = FALSE]

  # initialize beta with zeros; nrow = # of predictors, ncol = # of lambda values
  # this will create columns of zeros for betas corresponding to singular columns
  untransformed_beta <- matrix(0,
                               nrow = (p + 1), # + 1 is for the intercept; see note below
                               ncol = ncol(std_scale_beta)) # again, # lambda values

  ## Note: for in-memory matrix data, unpenalized columns are included in 'p',
  ##     the number of columns in X;
  ##     This differs from the filebacked data, where unpenalized columns are
  ##     *appended* to the columns of 'X'.


  # next, unscale the beta values for non-singular, non-intercept columns
  # NB: this requires the details of standardization (centering/scaling values
  # and indices of nonsingular columns). The details of how these details are
  # passed around varies depending on whether data are stored filebacked,
  # hence, the division into cases below:
  # in-memory cases
  if (length(std_X_details$ns) == length(std_X_details$scale)) {
    # case 1: ns and center/scale values have same length
    untransformed_b2 <- sweep(x = b,
                              MARGIN = 1,
                              STATS = std_X_details$scale,
                              FUN = "/")

    # fill in the un-transformed values
    untransformed_beta[std_X_details$ns + 1, ] <- untransformed_b2 # again, the + 1 is for the intercept
    untransformed_beta[1, ] <- a - crossprod(std_X_details$center,
                                             untransformed_b2)
  } else if ((length(std_X_details$ns) != length(std_X_details$scale)) &&
             (length(std_X_details$ns) == nrow(b))) {
    # case 2: ns and center/scale values **do not** have same length, but ns is
    #   equal to the number of rows of 'b'
    untransformed_b2 <- sweep(x = b,
                              MARGIN = 1,
                              STATS = std_X_details$scale[std_X_details$ns],
                              FUN = "/")

    # fill in the un-transformed values
    untransformed_beta[std_X_details$ns + 1, ] <- untransformed_b2 # again, the + 1 is for the intercept
    untransformed_beta[1, ] <- a - crossprod(std_X_details$center[std_X_details$ns],
                                             untransformed_b2)
  } else {
    # case 3: ns and center/scale values **do not** have same length, and ns is
    #  *NOT* equal to the number of rows of 'b'. This is a scenario that may arise
    #   in cross-validation.

    untransformed_b2 <- sweep(x = b[std_X_details$ns, ],
                              MARGIN = 1,
                              STATS = std_X_details$scale[std_X_details$ns],
                              FUN = "/")
    # fill in the un-transformed values
    untransformed_beta[std_X_details$ns + 1, ] <- untransformed_b2 # again, the + 1 is for the intercept
    untransformed_beta[1, ] <- a - crossprod(std_X_details$center[std_X_details$ns],
                                             untransformed_b2)
  }

  if (use_names) {

    rownames(untransformed_beta) <- c("(Intercept)",
                                      std_X_details$X_colnames)

  }


  # Final step: return un-transformed beta values
  return(untransformed_beta)
}

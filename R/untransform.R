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


untransform <- function(std_scale_beta, p, std_X_details, fbm_flag, use_names = TRUE) {

  # goal: reverse the PRE-ROTATION standardization #
  # partition the values from Step 1 into intercept and non-intercept parts
  a <- std_scale_beta[1, , drop = FALSE] # this is the intercept
  b <- std_scale_beta[-1, , drop=FALSE]

  # initialize beta with zeros; nrow = # of predictors, ncol = # of lambda values
  # this will create columns of zeros for betas corresponding to singular columns
  if (fbm_flag) {
    untransformed_beta <- Matrix::Matrix(0,
                                 nrow = (p + length(std_X_details$unpen) + 1), # + 1 is for the intercept
                                 ncol = ncol(std_scale_beta), # this is the number of lambda values
                                 sparse = TRUE)
  } else {
    untransformed_beta <- matrix(0,
                                 nrow = (p + length(std_X_details$unpen) + 1), # + 1 is for the intercept
                                 ncol = ncol(std_scale_beta)) # again, # lambda values
  }


  # next, unscale the beta values for non-singular, non-intercept columns
  # NB: this requires the details of standardization (centering/scaling values
  # and indices of nonsingular columns). The details of how these details are
  # passed around varies depending on whether data are stored filebacked,
  # hence, the division into cases below:
  if (fbm_flag){
    # filebacked cases
    if (length(std_X_details$ns) == length(std_X_details$scale)){
      # case 1: ns and center/scale values have same length
      untransformed_b2 <- sweep(x = b,
                                MARGIN = 1,
                                STATS = std_X_details$scale,
                                FUN = "/")

      # fill in the un-transformed values
      untransformed_beta[std_X_details$ns+1,] <- untransformed_b2 # again, the + 1 is for the intercept
      cp <- apply(X = untransformed_b2, 2, function(c){crossprod(std_X_details$center, c)})
      untransformed_beta[1,] <- a - cp
    } else {
      # case 2: ns and center/scale values **do not** have same length
      # (this will often be the case in cross-validation, where features can
      # become constant in a given fold)
      untransformed_b2 <- sweep(x = b[std_X_details$ns,],
                                MARGIN = 1,
                                STATS = std_X_details$scale[std_X_details$ns],
                                FUN = "/")

      # fill in the un-transformed values
      untransformed_beta[std_X_details$ns+1,] <- untransformed_b2 # again, the + 1 is for the intercept
      cp <- apply(X = untransformed_b2, 2, function(c){crossprod(std_X_details$center[std_X_details$ns], c)})
      untransformed_beta[1,] <- a - cp
    }

  } else {
    # in-memory cases
    if (length(std_X_details$ns) == length(std_X_details$scale)){
      # case 1: ns and center/scale values have same length
      untransformed_b2 <- sweep(x = b,
                                MARGIN = 1,
                                STATS = std_X_details$scale,
                                FUN = "/")

      # fill in the un-transformed values
      untransformed_beta[std_X_details$ns+1,] <- untransformed_b2 # again, the + 1 is for the intercept
      untransformed_beta[1,] <- a - crossprod(std_X_details$center,
                                              untransformed_b2)
    } else {
      # case 2: ns and center/scale values **do not** have same length
      untransformed_b2 <- sweep(x = b,
                                MARGIN = 1,
                                STATS = std_X_details$scale[std_X_details$ns],
                                FUN = "/")

      # fill in the un-transformed values
      untransformed_beta[std_X_details$ns+1,] <- untransformed_b2 # again, the + 1 is for the intercept
      untransformed_beta[1,] <- a - crossprod(std_X_details$center[std_X_details$ns],
                                              untransformed_b2)
    }

  }

  if (use_names) {
    if (!is.null(std_X_details$unpen_colnames)) {
      rownames(untransformed_beta) <- c("(Intercept)",
                                        std_X_details$unpen_colnames,
                                        std_X_details$X_colnames)
    } else {
      rownames(untransformed_beta) <- c("(Intercept)",
                                        std_X_details$X_colnames)
    }
  }


  # Final step: return un-transformed beta values
  return(untransformed_beta)
}





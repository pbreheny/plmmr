#' Untransform coefficient values back to the original scale 
#' 
#' This function unwinds the initial standardization of the data to obtain 
#' coefficient values on their original scale. It is called by plmm_format().
#' 
#' @param untransformed_b1 The estimated coefficients on the standardized scale
#' @param p  The number of columns in the original design matrix
#' @param std_X_details A list with 3 elements describing the standardized design matrix BEFORE rotation; this should have elements 'scale', 'center', and 'ns'
#' @param fbm_flag Logical: is the corresponding design matrix filebacked?
#' @param non_genomic Optional vector specifying which columns of the design matrix represent features that are *not* genomic, as these features are excluded from the empirical estimation of genomic relatedness. 
#' For cases where X is a filepath to an object created by `process_plink()`, this is handled automatically via the arguments to `process_plink()`.
#' For all other cases, 'non_genomic' defaults to NULL (meaning `plmm()` will assume that all columns of `X` represent genomic features).
#' @keywords internal


untransform <- function(untransformed_b1, p, std_X_details, fbm_flag, non_genomic){

  # goal: reverse the PRE-ROTATION standardization #
  # partition the values from Step 1 into intercept and non-intercept parts
  a <- untransformed_b1[1, , drop = FALSE] # this is the intercept 
  b <- untransformed_b1[-1, , drop=FALSE]
  # initialize beta with zeros; nrow = # of predictors, ncol = # of lambda values
  # this will create columns of zeros for betas corresponding to singular columns
  if (fbm_flag) {
    untransformed_beta <- Matrix::Matrix(0,
                                 nrow = (p + length(non_genomic) + 1), # + 1 is for the intercept
                                 ncol = ncol(untransformed_b1),
                                 sparse = TRUE)
  } else {
    untransformed_beta <- matrix(0,
                                 nrow = (p + length(non_genomic) + 1), # + 1 is for the intercept
                                 ncol = ncol(untransformed_b1))
  }
  

  # next, unscale the beta values for non-singular, non-intercept columns
  # NB: this requires the details of standardization (centering/scaling values
  # and indices of nonsingular columns). The details of how these details are 
  # passed around varies depending on whether data are stored filebacked, 
  # hence, the division into cases below:
  if (fbm_flag){
    untransformed_b2 <- sweep(x = b,
                              MARGIN = 1,
                              STATS = std_X_details$scale,
                              FUN = "/")

    # fill in the un-transformed values
    untransformed_beta[std_X_details$ns+1,] <- untransformed_b2 # again, the + 1 is for the intercept
    cp <- apply(X = untransformed_b2, 2, function(c){crossprod(std_X_details$center, c)})
    untransformed_beta[1,] <- a - cp
    # ... the bigger question is -- do we even need an intercept? 
  } else {
    # in-memory case 
    untransformed_b2 <- sweep(x = b,
                              MARGIN = 1,
                              STATS = std_X_details$scale[std_X_details$ns],
                              FUN = "/")
    
    # fill in the un-transformed values
    untransformed_beta[std_X_details$ns+1,] <- untransformed_b2 # again, the + 1 is for the intercept
    untransformed_beta[1,] <- a - crossprod(std_X_details$center[std_X_details$ns],
                                            untransformed_b2)
    
  }
  
  if (!is.null(std_X_details$X_colnames)) {
    rownames(untransformed_beta) <- c("(Intercept)",std_X_details$X_colnames)
  }
  
  # Final step: return un-transformed beta values 
  return(untransformed_beta)
}





#' A helper function to standardize matrices
#'
#' @param X a matrix
#'
#' @returns a list with the standardized matrix and its details
#' @keywords internal
#'
#' @details
#' This function is adapted from https://github.com/pbreheny/ncvreg/blob/master/R/std.R
#' NOTE: this function returns a matrix **in memory**. For standardizing filebacked
#' data, use `big_std()`  -- see src/big_standardize.cpp
standardize_matrix <- function(X){

  if (typeof(X) == 'integer') storage.mode(X) <- 'double'
  if (!inherits(X, "matrix")) {
    if (is.numeric(X)) {
      X <- matrix(as.double(X), ncol=1)
    } else {
      tmp <- try(X <- stats::model.matrix(~0+., data=X), silent=TRUE)
      if (inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix", call.=FALSE)
    }
}

    standardization <- .Call("in_mem_std", X, PACKAGE = 'plmmr')
    dimnames(standardization[[1]]) <- dimnames(X)
    ns <- which(standardization[[3]] > 1e-6)
    std_X <- standardization[[1]]
    std_X_details <- list()

    # difference from ncvreg::std(): instead of removing singular columns,
    #   this version fills constant columns with 0s to preserve the dimension
    #   of the matrix.
    std_X[,-ns] <- 0

    attr(std_X, "center") <- std_X_details$center <- standardization[[2]]
    attr(std_X, "scale") <- std_X_details$scale <- standardization[[3]]
    attr(std_X, "nonsingular") <- std_X_details$ns <- ns

    res <- list(std_X = std_X[,],
                std_X_details = std_X_details)

    return(res)
}

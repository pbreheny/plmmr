#'  a version of cbind() for file-backed matrices
#'
#' @param A in-memory data
#' @param B file-backed data
#' @param C file-backed placeholder for combined data
#' @param quiet Logical
#'
#' @return C, filled in with all column values of A and B combined
#' @keywords internal
big_cbind <- function(A, B, C, quiet){
  options(bigmemory.typecast.warning=FALSE)
  if (!quiet) {
    pb <- txtProgressBar(min = 0, max = ncol(A) + ncol(B), style = 3)
  }

  # fill in first few columns of C with A
  for (j in 1:ncol(A)) {
    C[, j] <- A[, j]
    if (!quiet) {
      setTxtProgressBar(pb, j)
    }
  }
  # fill in the rest of the columns of C with B
  for (j in 1:ncol(B)) {
    C[, j + ncol(A)] <- B[, j]
    if (!quiet) {
      setTxtProgressBar(pb, j)
    }
  }

  if (!quiet) {
    close(pb)
  }

  return(C)
}

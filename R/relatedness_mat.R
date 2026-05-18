#' Calculate a relatedness matrix
#'
#' Given a matrix of genotypes, this function estimates the genetic relatedness matrix (GRM,
#' also known as the RRM, see Hayes et al. 2009, \doi{10.1017/S0016672308009981}) among
#' the subjects: \eqn{\frac{1}{p}(XX^T)}, where X is standardized.
#'
#' @param X An n x p numeric matrix of genotypes (from *fully-imputed* data). Can be a filebacked `big.matrix` object.
#' Note: This matrix should *not* include non-genetic features.
#' @param std Logical: should `X` be standardized? If you set this to FALSE, you should have a good reason for doing so,
#' as standardization is a best practice.
#' @export
#'
#' @return An n x n numeric matrix capturing the genomic relatedness of the
#' samples represented in `X`. In our notation, we call this matrix `K` for 'kinship';
#' this is also known as the GRM or RRM.
#'
#' @examples
#' RRM <- relatedness_mat(X = admix$X)
#' RRM[1:5, 1:5]
relatedness_mat <- function(X, std = TRUE) {
  p <- ncol(X)

  if (inherits(X, "big.matrix")) {
    if (std) {
      Xdesc <- standardize_filebacked(X, outfile = nullfile(), quiet = TRUE)
      X <- bigmemory::attach.big.matrix(Xdesc$std_X)
      p <- p - sum(Xdesc$std_X_scale < 1e-3) # Don't count singular columns
    }

    XX <- bigalgebra::dgemm(TRANSA = "N",
                            TRANSB = "T",
                            A = X,
                            B = X)[,]
  } else if (std) {
    std_X <- ncvreg::std(X)
    XX <- tcrossprod(std_X)
    p <- ncol(std_X)
  } else {
    XX <- tcrossprod(X)
  }

  # scale by the number of non-singular columns in the design matrix
  XX / p
}

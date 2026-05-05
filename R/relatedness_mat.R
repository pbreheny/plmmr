#' Calculate a relatedness matrix
#'
#' Given a matrix of genotypes, this function estimates the genetic relatedness matrix (GRM,
#' also known as the RRM, see Hayes et al. 2009, \doi{10.1017/S0016672308009981}) among
#' the subjects: \eqn{\frac{1}{p}(XX^T)}, where X is standardized.
#'
#' @param X An n x p numeric matrix of genotypes (from *fully-imputed* data).
#' Note: This matrix should *not* include non-genetic features. Can be a filebacked `big.matrix` object.
#' @param std Logical: should `X` be standardized? If you set this to FALSE (which can only be done
#' if data are stored in memory), you should have a good reason for doing so, as standardization
#' is a best practice.
#' @param ns Optional vector of values indicating the indices of nonsingular features.
#' **Note**: If a filebacked `big.matrix` is passed to `X` along with a non-null `ns`,
#' a temporary deep copy will be created.
#'
#' @export
#'
#' @return An n x n numeric matrix capturing the genomic relatedness of the
#' samples represented in `X`. In our notation, we call this matrix `K` for 'kinship';
#' this is also known as the GRM or RRM.
#'
#' @examples
#' RRM <- relatedness_mat(X = admix$X)
#' RRM[1:5, 1:5]
relatedness_mat <- function(X, std = TRUE, ns = NULL) {
  p <- ncol(X)

  if (length(ns) %in% 1:(ncol(X) - 1)) {
    if(inherits(X, "big.matrix")) {
      X <- subset_filebacked(X,
                             new_file = "temp_X",
                             complete_samples = NULL,
                             ns = ns,
                             rds_dir = tempdir(),
                             outfile = nullfile(),
                             quiet = TRUE)$subset_X

      on.exit({
        gc()
        unlink(list.files(tempdir(), pattern = paste0("temp_X", "\\.(bk|desc)"),
                          full.names = TRUE),
               force = TRUE)
      })
    } else {
      X <- X[, ns]
    }
  }

  if (inherits(X, "big.matrix")) {
    if(std) {
      Xdesc <- standardize_filebacked(X, outfile = nullfile(), quiet = TRUE)
      X <- bigmemory::attach.big.matrix(Xdesc$std_X)
    }

    XX <- bigalgebra::dgemm(TRANSA = "N",
                            TRANSB = "T",
                            A = X,
                            B = X)[,]
  } else if (std) {
    XX <- tcrossprod(ncvreg::std(X))
  } else {
    XX <- tcrossprod(X)
  }

  # scale by p, the number of columns in the design matrix (including constant features)
  XX / p
}

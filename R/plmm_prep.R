#' PLMM prep: a function to run checks, eigendecomposition, and rotation prior to fitting a PLMM model
#'
#' This is an internal function for `plmm()`
#'
#' @param std_X Column standardized design matrix. May include clinical covariates and other non-SNP data.
#' @param std_X_n The number of observations in `std_X` (integer)
#' @param std_X_p The number of features in `std_X` (integer)
#' @param n The number of instances in the *original* design matrix X. This should not be altered by standardization.
#' @param p The number of features in the *original* design matrix X, including constant features
#' @param centered_y Continuous outcome vector, centered.
#' @param penalty_factor A multiplicative factor for the penalty applied to each coefficient.
#' @param K Similarity matrix used to rotate the data. This should either be:
#'            (1) a known matrix that reflects the covariance of y,
#'            (2) an estimate (Default is \eqn{\frac{1}{p}(XX^T)}), or
#'            (3) a list with components `s` and `U`, as returned by a previous `plmm()` model fit on the same data.
#' @param eta Optional argument to input a specific eta term rather than estimate it from the data. If K is a known covariance matrix that is full rank, this should be 1.
#' @param fbm_flag Logical: is `std_X` a filebacked `big.matrix` object? This is set internally by `plmm()`.
#' @param trace If set to TRUE, inform the user of progress by announcing the beginning of each step of the modeling process. Default is FALSE.
#' @param ... Not used
#'
#' @return List with these components:
#' * `std_X`: Standardized design matrix. If design matrix is filebacked, the descriptor for the filebacked data is returned using `bigmemory::describe()`.
#' * `centered_y`: Vector of centered outcomes
#' * `K`: Similarity matrix
#' * `s`: Vector of the non-zero eigenvalues of `K`
#' * `U`: Matrix of eigenvectors of `K` associated with `s` (same as left singular values of X).
#' * `eta`: The numeric value of the estimated eta parameter
#' * `penalty_factor` A multiplicative factor for the penalty applied to each coefficient.
#' * `incpt_flag` Logical: Does the model require fitting an intercept?
#' * `trace`: If set to TRUE, inform the user of progress by announcing the beginning of each step of the modeling process
#'
#' @keywords internal
#'
plmm_prep <- function(std_X,
                      std_X_n,
                      std_X_p,
                      n,
                      p,
                      centered_y,
                      penalty_factor,
                      K = NULL,
                      eta = NULL,
                      fbm_flag,
                      trace = NULL,
                      ...) {

  ## coercion
  U <- s <- eta <- NULL

  # First: handle the cases where no decomposition is needed ------------------

  # case 1: K is user-supplied diagonal matrix (like a weighted lm())
  flag1 <- !is.null(K) & ifelse(!(inherits(K, "matrix") || inherits(K, "dcGMatrix")),
                                FALSE,
                                Matrix::isDiagonal(K))
  if (flag1) {
    if (trace) {
      cat("Using supplied diagonal matrix for K, similar to a lm() with weights.\n")
    }
    s <- sort(diag(K), decreasing = TRUE)
    U <- diag(nrow = n)[, order(diag(K), decreasing = TRUE)]
  }
  # case 2: K is a user-supplied list
  flag2 <- !is.null(K) & inherits(K, "list")
  if (flag2) {
    if (trace) {
      cat("K is a list; will pass U,s components from list to model fitting.\n")
    }
    s <- K$s # no need to adjust singular values by p
    U <- K$U
  }

  # otherwise, need to do eigendecomposition -----------------------------
  if (sum(c(flag1, flag2)) == 0) {
    if (trace) {
      cat("Starting decomposition.\n")
    }
    # set default K: if not specified and not diagonal, use realized relatedness matrix
    if (is.null(K) && is.null(s)) {
      # NB: the is.null(s) keeps you from overwriting the 3 preceding special cases

      if (trace) cat("Calculating the eigendecomposition of K\n")
      eigen_res <- eigen_K(std_X)
      K <- eigen_res$K
      s <- eigen_res$s
      U <- eigen_res$U

    } else {
      # last case: K is a user-supplied matrix
      eigen_res <- eigen(K, symmetric = TRUE)
      nz <- eigen_res$values > 1e-4
      s <- eigen_res$values[nz]
      U <- eigen_res$vectors[, nz, drop = FALSE]
    }

  }

  # error check: what if the combination of args. supplied was none of the cases above?
  if (is.null(s) || is.null(U)) {
    stop("\nSomething is wrong in the eigendecomposition.
    \nThe combination of supplied arguments does not match any cases handled in
         \n plmm_prep(), the internal function called by plmm() to do this step of the modeling process.
         \n Re-examine the supplied arguments -- here are some common mistakes:
         \n Is the K argument you supplied something other than a list, a matrix, or a filepath to an RDS file with one of those objects?
         \n \tDid you supply a list to K? Check its element names -- they must be 's' and 'U'.")
  }

  # check if an intercept is required for model fitting based on U - not required if column means of U are 0
  incpt_flag <- !isTRUE(all.equal(rep(0, ncol(U)), colSums(U)))

  if (incpt_flag) {
    if (!fbm_flag) {
      std_X <- cbind(1, std_X)
      penalty_factor <- c(0, penalty_factor)
    } else {
      incpt <- as.matrix(rep(1, std_X_n), ncol = std_X_n)

      filename <- bigmemory::describe(std_X)@description$filename |> tools::file_path_sans_ext()
      dirname <- bigmemory::describe(std_X)@description$dirname

      if (trace) {
        cat("\nNote: Due to the format of your provided K, the model will require an intercept in the design matrix. ",
	    "This will be written to a filebacked big.matrix in the same directory where you created your design as ",
            filename, "_incpt.bk/desc. \n", sep = "")
      }

      # Overwrite old backing files
      list.files(dirname, pattern = paste0(filename, "_incpt"), full.names = TRUE) |>
        unlink(force = TRUE)

      tmp_matrix <- big.matrix(nrow = std_X_n,
                               ncol = std_X_p + 1,
                               type = "double",
                               backingfile = paste0(filename, "_incpt.bk"),
                               backingpath = dirname,
                               descriptorfile = paste0(filename, "_incpt.desc"))

      std_X <- big_cbind(A = incpt,
                         B = std_X,
                         C = tmp_matrix,
                         quiet = TRUE)
      penalty_factor <- c(0, penalty_factor)
    }
  }

  # estimate eta if needed; otherwise, use the user-supplied value (this option is mainly used for simulation studies)
  if (is.null(eta)) {
    eta <- estimate_eta(n = std_X_n, s = s, U = U, y = centered_y, incpt_flag = incpt_flag)
  }

  # return values to be passed into plmm_fit():
  list(
    std_X = std_X,
    centered_y = centered_y,
    K = K, # Note: need this for CV (see call to construct_variance() within cv_plmm())
    s = s,
    U = U,
    eta = eta, # carry eta over to fit
    penalty_factor = penalty_factor,
    incpt_flag = incpt_flag,
    trace = trace)
}

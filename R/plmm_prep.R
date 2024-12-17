#' PLMM prep: a function to run checks, SVD, and rotation prior to fitting a PLMM model
#' This is an internal function for \code{cv_plmm}
#'
#' @param std_X Column standardized design matrix. May include clinical covariates and other non-SNP data.
#' @param std_X_n The number of observations in std_X (integer)
#' @param std_X_p The number of features in std_X (integer)
#' @param n The number of instances in the *original* design matrix X. This should not be altered by standardization.
#' @param p The number of features in the *original* design matrix X, including constant features
#' @param centered_y Continuous outcome vector, centered.
#' @param k An integer specifying the number of singular values to be used in the approximation of the rotated design matrix. This argument is passed to `RSpectra::svds()`. Defaults to `min(n, p) - 1`, where n and p are the dimensions of the **standardized** design matrix.
#' @param K Similarity matrix used to rotate the data. This should either be a known matrix that reflects the covariance of y, or an estimate (Default is \eqn{\frac{1}{p}(XX^T)}, where X is standardized). This can also be a list, with components d and u (as returned by choose_k)
#' @param diag_K Logical: should K be a diagonal matrix? This would reflect observations that are unrelated, or that can be treated as unrelated. Passed from `plmm()`.
#' @param eta_star Optional argument to input a specific eta term rather than estimate it from the data. If K is a known covariance matrix that is full rank, this should be 1.
#' @param fbm_flag Logical: is std_X an FBM type object? This is set internally by `plmm()`.
#' @param trace If set to TRUE, inform the user of progress by announcing the beginning of each step of the modeling process. Default is FALSE.
#' @param ... Not used yet
#'
#' @returns List with these components:
#' * centered_y: The vector of centered outcomes
#' * std_X: standardized design matrix
#' * K: a list with 2 elements. (1) s: vector with the eigenvalues of K,
#'  and (2) U: the eigenvectors of K (same as left singular values of X).
#' * eta: the numeric value of the estimated eta parameter
#' * trace: logical.
#'
#' @keywords internal

plmm_prep <- function(std_X,
                      std_X_n,
                      std_X_p,
                      n,
                      p,
                      centered_y,
                      k = NULL,
                      K = NULL,
                      diag_K = NULL,
                      eta_star = NULL,
                      fbm_flag,
                      penalty_factor = rep(1, ncol(std_X)),
                      trace = NULL,
                      ...){


  ## coersion
  U <- s <- eta <- NULL

  if(is.null(k)){
    trunc_flag <- FALSE
  }

  if(!is.null(k)){

    if(k < min(std_X_n, std_X_p)){
      trunc_flag <- TRUE
    } else {
      stop("\nk value is out of bounds.
         \nIf specified, k must be an integer in the range from 1 to min(nrow(std_X), ncol(std_X)).
         \nwhere std_X does not include singular columns. For help detecting singularity,
         \nsee ncvreg::std()")
    }

  }

  # set default: if diag_K not specified, set to false
  if(is.null(diag_K)){diag_K <- FALSE}

  # handle the cases where no decomposition is needed:
  # case 1: K is the identity matrix
  flag1 <- diag_K & is.null(K)
  if(flag1){
    if (trace) {(cat("Using identity matrix for K.\n"))}
    s <- rep(1, n)
    U <- diag(nrow = n)
  }
  # case 2: K is user-supplied diagonal matrix (like a weighted lm())
  flag2 <- diag_K & !is.null(K) & ('matrix' %in% class(K))
  if(flag2){
    if (trace) {(cat("Using supplied diagonal matrix for K, similar to a lm() with weights.\n"))}
    s <- sort(diag(K), decreasing = TRUE)
    U <- diag(nrow = n)[,order(diag(K), decreasing = TRUE)]
  }
  # case 3: K is a user-supplied list
  flag3 <- !is.null(K) & ('list' %in% class(K))
  if( flag3) {
    if (trace) {cat("K is a list; will pass U,s components from list to model fitting.\n")}
    s <- K$s # no need to adjust singular values by p
    if ('FBM' %in% class(K$U)){
      U <- K$U[,]
    } else {
      U <- K$U
    }
  }

  # otherwise, need to do eigendecomposition:
  if (sum(c(flag1, flag2, flag3)) == 0) {
    if (trace) {cat("Starting decomposition.\n")}
    # set default K: if not specified and not diagonal, use realized relatedness matrix
    if (is.null(K) & is.null(s)) {
      # NB: the is.null(s) keeps you from overwriting the 3 preceding special cases

      if (trace) cat("Calculating the eigendecomposition of K\n")
      eigen_res <- eigen_K(std_X, fbm_flag = fbm_flag)
      K <- eigen_res$K
      s <- eigen_res$s
      U <- eigen_res$U

    } else {
      # last case: K is a user-supplied matrix
      eigen_res <- eigen(K)
      s <- eigen_res$values*(1/std_X_p)
      # note: our definition of the RRM averages over the number of features used to calculate K
      U <- eigen_res$vectors
    }

  }

  # error check: what if the combination of args. supplied was none of the SVD cases above?
  if (is.null(s) | is.null(U)){
    stop("\nSomething is wrong in the SVD/eigendecomposition.
    \nThe combination of supplied arguments does not match any cases handled in
         \n svd_X(), the internal function called by plmm() via plmm_prep().
         \n Re-examine the supplied arguments -- here are some common mistakes:
         \n \tDid you supply a list to K? Check its element names -- they must be 's' and 'U'.
         \n \t \t *If you used choose_k(), make sure you are supplying the 'svd_K' element returned from that function's result as the 'K' here.*
         \n \tDid you intend to set diag_K = TRUE?
         \n \tDid you set diag_K = TRUE and specifiy a k value at the same time? This combination of arguments is incompatible.")
  }

  # estimate eta if needed
  if (is.null(eta_star)) {
    eta <- estimate_eta(n = std_X_n, s = s, U = U, y = centered_y)
  } else {
    # otherwise, use the user-supplied value (this option is mainly for simulation)
    eta <- eta_star
  }

  # return values to be passed into plmm_fit():
  ret <- structure(list(
    std_X = std_X,
    centered_y = centered_y,
    K = K, # Note: need this for CV (see call to construct_variance() within cv_plmm())
    s = s,
    U = U,
    eta = eta, # carry eta over to fit
    trace = trace))

  return(ret)
}

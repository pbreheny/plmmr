#' plmm_checks
#'
#' @param X Design matrix object or a string with the file path to a design matrix. If a string, string will be passed to `get_data()`.
#' * Note: X may include clinical covariates and other non-SNP data, but no missing values are allowed.
#' @param y Continuous outcome vector. Defaults to NULL, assuming that the outcome is the 6th column in the .fam PLINK file data. Can also be a user-supplied numeric vector.
#' @param col_names Optional vector of column names for design matrix. Defaults to NULL.
#' @param non_genomic Optional vector specifying which columns of the design matrix represent features that are *not* genomic, as these features are excluded from the empirical estimation of genomic relatedness.
#' For cases where X is a filepath to an object created by `process_plink()`, this is handled automatically via the arguments to `process_plink()`.
#' For all other cases, 'non_genomic' defaults to NULL (meaning `plmm()` will assume that all columns of `X` represent genomic features).
#' @param K Similarity matrix used to rotate the data. This should either be (1) a known matrix that reflects the covariance of y, (2) an estimate (Default is \eqn{\frac{1}{p}(XX^T)}), or (3) a list with components 'd' and 'u', as returned by choose_k().
#' @param diag_K Logical: should K be a diagonal matrix? This would reflect observations that are unrelated, or that can be treated as unrelated. Defaults to FALSE.
#'  Note: plmm() does not check to see if a matrix is diagonal. If you want to use a diagonal K matrix, you must set diag_K = TRUE.
#' @param eta_star Optional argument to input a specific eta term rather than estimate it from the data. If K is a known covariance matrix that is full rank, this should be 1.
#' @param penalty The penalty to be applied to the model. Either "MCP" (the default), "SCAD", or "lasso".
#' @param penalty_factor A multiplicative factor for the penalty applied to each coefficient. If supplied, penalty_factor must be a numeric vector of length equal to the number of columns of X.
#' The purpose of penalty_factor is to apply differential penalization if some coefficients are thought to be more likely than others to be in the model.
#' In particular, penalty_factor can be 0, in which case the coefficient is always in the model without shrinkage.
#' @param init Initial values for coefficients. Default is 0 for all columns of X.
#' @param gamma The tuning parameter of the MCP/SCAD penalty (see details). Default is 3 for MCP and 3.7 for SCAD.
#' @param alpha Tuning parameter for the Mnet estimator which controls the relative contributions from the MCP/SCAD penalty and the ridge, or L2 penalty. alpha=1 is equivalent to MCP/SCAD penalty, while alpha=0 would be equivalent to ridge regression. However, alpha=0 is not supported; alpha may be arbitrarily small, but not exactly 0.
#' @param dfmax Option to be added soon: Upper bound for the number of nonzero coefficients. Default is no upper bound. However, for large data sets, computational burden may be heavy for models with a large number of nonzero coefficients.
#' @param trace If set to TRUE, inform the user of progress by announcing the beginning of each step of the modeling process. Default is FALSE.
#' @param save_rds  Optional: if a filepath and name is specified (e.g., `save_rds = "~/dir/my_results.rds"`), then the model results are saved to the provided location. Defaults to NULL, which does not save the result.
#' @param return_fit Optional: a logical value indicating whether the fitted model should be returned as a `plmm` object in the current (assumed interactive) session. Defaults to TRUE.
#' @param ... Additional arguments to `get_data()`
#'
#' @keywords internal
#'
plmm_checks <- function(X,
                        y = NULL,
                        col_names = NULL,
                        non_genomic = NULL,
                        K = NULL,
                        diag_K = NULL,
                        eta_star = NULL,
                        penalty = "lasso",
                        penalty_factor = NULL,
                        init = NULL,
                        gamma,
                        alpha = 1,
                        dfmax = NULL,
                        trace = FALSE,
                        save_rds = NULL,
                        return_fit = TRUE,
                        ...){

  # check X types -------------------------------------------------
  if (!any(class(X) %in% c("character", "matrix"))) {
    stop("\nThe X argument must be either (1) a numeric matrix or (2) a character
         string specifying a filepath to an RDS object that you created using
         process_delim() or process_plink().")
  }
  # read in X -----------------------------------------------------
  if("character" %in% class(X)){
    # case 1: X is a filebacked matrix from process_delim() or process_plink()
    dat <- get_data(path = X, trace = trace, ...)
    X <- std_X <- dat$std_X
    std_X_n <- dat$std_X_n
    std_X_p <- dat$std_X_p
    col_names <- dat$X_colnames
    genomic <- index_std_X(std_X_p = std_X_p, non_genomic = dat$non_gen)

    # create a list that captures the centering/scaling for std_X;
    # will need this later, see `untransform()`
    std_X_details <- list(
      center = dat$std_X_center,
      scale = dat$std_X_scale,
      ns = dat$ns)

    if ('colnames' %in% names(dat) | 'std_X_colnames' %in% names(dat)){
      std_X_details$X_colnames <- dat$colnames
      std_X_details$X_rownames <-  dat$rownames
      std_X_details$std_X_rownames <- dat$std_X_rownames
      std_X_details$std_X_colnames <-  dat$std_X_colnames
    } else if (!missing(col_names)){
      std_X_details$X_colnames <- col_names
      std_X_details$std_X_colnames <- col_names[std_X_details$ns]
    }

    if(inherits(std_X, "big.matrix")){
      fbm_flag <- TRUE
    } else {
      fbm_flag <- FALSE
    }
  } else {
    # case 2: X is a matrix in-memory
    if (!inherits(X, "matrix")) {
      tmp <- try(X <- stats::model.matrix(~0+., data=X), silent=TRUE)
      if (inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix", call.=FALSE)
    }
    if (typeof(X)=="integer") storage.mode(X) <- "double"
    if (typeof(X)=="character") stop("if X is a matrix, it must be a numeric matrix", call.=FALSE)

    fbm_flag <- FALSE

    # designate the dimensions of the original design matrix
    n <- nrow(X)
    p <- ncol(X) - length(non_genomic)
    # Note: p only includes genomic markers (in GWAS context);
    # This parallels the behavior of `process_plink()` for filebacked data

    # set default column names of X
    if (is.null(col_names) & !is.null(attr(X, "dimnames")[[2]])) {
      col_names <- attr(X, "dimnames")[[2]]
    }

    # handle standardization (for in-memory matrix case)
      std_res <- standardize_matrix(X, penalty_factor)
      std_X <- std_res$std_X
      std_X_details <- std_res$std_X_details
      penalty_factor <- std_res$penalty_factor

    # designate dimensions of the standardized data
      std_X_n <- nrow(std_X)
      std_X_p <- ncol(std_X)
      genomic <- index_std_X(std_X_p = std_X_p, non_genomic = non_genomic)


  } # this bracket closes case 2 (the matrix case )

  #  check y types & read y -------------------------------
  # if y is null, use .fam file
  if(is.null(y)){
    # default: use data from 6th column of .fam file
    if ("fam" %in% names(dat)){
      # case 1: if handle_missing_phen was set to 'prune', use 'complete_phen' for subsetting
      if (length(dat$complete_phen) > 1) {
        y <- dat$fam$affection[dat$complete_phen]
      } else {
        # case 2: use all dat$fam$affection values
        y <- dat$fam$affection
      }

    } else {
      stop("If the data did not come from process_plink(), you must specify a
           'y' argument\n")
    }

  }

  if (!is.double(y)) {
    op <- options(warn=2)
    on.exit(options(op))
    y <- tryCatch(
      error = function(cond) stop("y must be numeric or able to be coerced to numeric\n", call.=FALSE),
      as.double(y))
    options(op)
  }

  # set up defaults --------------------------------------------------
  if(is.null(dfmax)){dfmax <- std_X_p + 1}

  # default: penalize everything except the intercept, which we will add later
  if(is.null(penalty_factor)){penalty_factor <- rep(1, std_X_p)}

  # set default init
  if(is.null(init)){init <- rep(0, std_X_p)}


  # set default gamma (gamma not used in 'lasso' option)
  if (missing(gamma)) gamma <- switch(penalty, SCAD = 3.7, MCP = 3, lasso = 1)

  # error checking ------------------------------------------------------------
  if(!fbm_flag){
    # error check for matrix X
    if (length(y) != std_X_n) stop("X and y do not have the same number of observations", call.=FALSE)

    if (any(is.na(y)) | any(is.na(std_X))) stop("Missing data (NA's) detected in outcome, 'y'.
                                            \nTake actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to plmm", call.=FALSE)
    if (length(penalty_factor)!=std_X_p) stop("Dimensions of penalty_factor and X do not match", call.=FALSE)
  } else {
    #  error checking for filebacked X
    if (length(y) != std_X_n) stop("X and y do not have the same number of observations", call.=FALSE)
    if (any(is.na(y))) stop("Missing data (NA's) detected in the outcome.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to plmm", call.=FALSE)
    if (length(penalty_factor)!=std_X_p) stop("Dimensions of penalty_factor and X do not match", call.=FALSE)

    # if penalty factor has 0s, these must be contiguous
    # (necessary for subsetting later -- see setup_lambda(), for example)
    if (any(penalty_factor < 1e-8)){
      pf <- which(penalty_factor < 1e-8)
      if (!identical(pf, 1:length(pf))) {
        stop("It looks like you are trying to make some covariates *not* be penalized, as
           you have some penalty_factor values set to 0.\n
           If this is your intention, you must make all 'n' unpenalized covariates appear
           as the first 'n' columns in your design matrix.\n This is needed for subsetting later
             (for now, subsetting filebacked matrices requires a contiguous submatrix).\n")
      }
    }


  }

  # check K types -------------------------------------------------------
  if (!is.null(K)){
    # first, check type/class:
    if (!inherits(K, "matrix") & !is.list(K)) {
      tmp <- try(K <- stats::model.matrix(~0+., data=K), silent=TRUE)
      if (inherits(tmp, "try-error")) stop("K must be either (1) able to be coerced to a matrix or (2) be a list.", call.=FALSE)
    }
    if (typeof(K)=="integer") storage.mode(std_X) <- "double" # change K to X
    if (typeof(K)=="character") stop("K must be a numeric matrix", call.=FALSE)
    if (is.list(K)) {
      if(!('s' %in% names(K) & 'U' %in% names(K))){stop('Components s and U not both found in list supplied for K.')}
    }
    # last check: look at dimensions
    if (is.matrix(K)){
      if (dim(K)[1] != std_X_n || dim(K)[2] != std_X_n) stop("Dimensions of K and X do not match", call.=FALSE)
    } else if (is.list(K)) {
      if (std_X_n != nrow(K$U)) stop("Dimensions of K and X do not match.")
    }

  }

  # return list for model preparation ---------------------------------

  ret <- list(
    std_X = std_X,
    std_X_details = std_X_details,
    std_X_n = std_X_n,
    std_X_p = std_X_p,
    genomic = genomic,
    std_X_n = std_X_n,
    std_X_p = std_X_p,
    col_names = col_names,
    y = y,
    centered_y = y - mean(y),
    K = K,
    diag_K = diag_K,
    fbm_flag = fbm_flag,
    penalty = penalty,
    penalty_factor = penalty_factor,
    gamma = gamma,
    init = init
  )

  if (exists('dat')){
    ret$n <- dat$n
    ret$p <- dat$p
    ret$non_genomic <- dat$non_gen
    ret$dat <- dat
  } else {
    ret$n <- n
    ret$p <- p
    ret$non_genomic <- non_genomic
  }

  if (trace & !is.null(save_rds)){cat("Your results will be saved to ", save_rds, "\n")}

  return(ret)

}

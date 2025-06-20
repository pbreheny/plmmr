#' plmm_checks
#' @param design The design object, as created by `create_design()`
#' @param K Similarity matrix used to rotate the data. This should either be (1) a known matrix that reflects the covariance of y, (2) an estimate (Default is \eqn{\frac{1}{p}(XX^T)}), or (3) a list with components 'd' and 'u', as returned by choose_k().
#' @param diag_K Logical: should K be a diagonal matrix? This would reflect observations that are unrelated, or that can be treated as unrelated. Defaults to FALSE.
#'  Note: plmm() does not check to see if a matrix is diagonal. If you want to use a diagonal K matrix, you must set diag_K = TRUE.
#' @param eta_star Optional argument to input a specific eta term rather than estimate it from the data. If K is a known covariance matrix that is full rank, this should be 1.
#' @param penalty The penalty to be applied to the model. Either "MCP" (the default), "SCAD", or "lasso".
#' @param init Initial values for coefficients. Default is 0 for all columns of X.
#' @param gamma The tuning parameter of the MCP/SCAD penalty (see details). Default is 3 for MCP and 3.7 for SCAD.
#' @param alpha Tuning parameter for the Mnet estimator which controls the relative contributions from the MCP/SCAD penalty and the ridge, or L2 penalty. alpha=1 is equivalent to MCP/SCAD penalty, while alpha=0 would be equivalent to ridge regression. However, alpha=0 is not supported; alpha may be arbitrarily small, but not exactly 0.
#' @param trace If set to TRUE, inform the user of progress by announcing the beginning of each step of the modeling process. Default is FALSE.
#' @param save_rds  Optional: if a filepath and name is specified (e.g., `save_rds = "~/dir/my_results.rds"`), then the model results are saved to the provided location. Defaults to NULL, which does not save the result.
#' @param return_fit Optional: a logical value indicating whether the fitted model should be returned as a `plmm` object in the current (assumed interactive) session. Defaults to TRUE.
#' @param ... Additional arguments to `get_data()`
#'
#' @keywords internal
#'
#' @returns A list of parameters to pass on to model fitting. The list includes the standardized design matrix, the outcome, and meta-data
#'
plmm_checks <- function(design,
                        K = NULL,
                        diag_K = NULL,
                        eta_star = NULL,
                        penalty = "lasso",
                        init = NULL,
                        gamma,
                        alpha = 1,
                        trace = FALSE,
                        save_rds = NULL,
                        return_fit = TRUE,
                        ...) {

  # read in X -----------------------------------------------------
  # TODO: maybe need to add a 'catch' here for incorrect list input...

  if (!(inherits(design, "plmm_design") || inherits(design, "character"))) {
    stop("Input to 'design' in plmm() must be either a plmm_design object output from create_design()
         or an .rds filepath to a saved list output from create_design()")
  }
  if ("character" %in% class(design)) {
    design <- get_data(path = design, trace = trace, ...)
    std_X <- design$std_X
  } else if (inherits(design$std_X, "big.matrix")) {
    std_X <- bigmemory::attach.big.matrix(design$std_X)
  } else if (inherits(design$std_X, "matrix")) {
    std_X <- design$std_X
  }

  std_X_n <- design$std_X_n
  std_X_p <- design$std_X_p

  # create a list that captures the centering/scaling for std_X;
  # will need this later, see `untransform()`
  std_X_details <- list(
    center = design$std_X_center,
    scale = design$std_X_scale,
    ns = design$ns,
    unpen = design$unpen,
    unpen_colnames = design$unpen_colnames,
    X_colnames = design$X_colnames,
    X_rownames = design$X_rownames,
    std_X_rownames = design$std_X_rownames,
    std_X_colnames =  design$std_X_colnames)

  if (inherits(std_X, "big.matrix")) {
    fbm_flag <- TRUE
  } else {
    fbm_flag <- FALSE
  }

  if (inherits(design$y, "matrix") || inherits(design$y, "data.frame")) {
    y <- unlist(design$y[, 1]) # design$y will have one column
  } else {
    y <- design$y
  }

  penalty_factor <- design$penalty_factor

  # set up defaults --------------------------------------------------
  # set default init
  if (is.null(init)) {
    init <- rep(0, std_X_p)
    }

  # set default gamma (gamma not used in 'lasso' option)
  if (missing(gamma)) gamma <- switch(penalty, SCAD = 3.7, MCP = 3, lasso = 1)

  # error checking design matrix  ---------------------------------------------
  if (length(y) != std_X_n) stop("X and y do not have the same number of observations", call.=FALSE)
  if (length(penalty_factor) != std_X_p) stop("Dimensions of penalty_factor and X do not match; something is off in the supplied design", call.=FALSE)

  # check K types -------------------------------------------------------
  if (!is.null(K)) {
    # first, check type/class:
    if (!inherits(K, "matrix") && !is.list(K)) {
      tmp <- try(K <- stats::model.matrix(~0+., data = K), silent = TRUE)
      if (inherits(tmp, "try-error")) stop("K must be either (1) able to be coerced to a matrix or (2) be a list.", call.=FALSE)
    }
    if (typeof(K) == "integer") storage.mode(std_X) <- "double" # change K to X
    if (typeof(K) == "character") stop("K must be a numeric matrix", call. = FALSE)
    if (is.list(K)) {
      if(!("s" %in% names(K) & "U" %in% names(K))) {
        stop("Components s and U not both found in list supplied for K.")
        }
    }
    # last check: look at dimensions
    if (is.matrix(K)) {
      if (dim(K)[1] != std_X_n || dim(K)[2] != std_X_n) stop("Dimensions of K and X do not match", call.=FALSE)
    } else if (is.list(K)) {
      if (std_X_n != nrow(K$U)) stop("Dimensions of K and X do not match.")
    }

  }

  # for in in-memory designs, set plink flag to FALSE
  if (is.null(design$is_plink)) {
    plink_flag <- FALSE
  } else if (design$is_plink) {
    plink_flag <- TRUE
  } else {
    plink_flag <- FALSE
  }

  # warn about bigalgebra issues if this is user's first call to filebacked analysis
  # Check if the warning has been shown
  if (!.plmmr_env$warning_shown && fbm_flag) {

    warning("
Note: version 1.1.2 of bigalgebra produces occasional 'stack imbalance'
warnings; see
https://github.com/fbertran/bigalgebra/issues/2
These warnings don't affect the results of plmm(), although you may wish to
download the previous version of the package to avoid these warnings:

remotes::install_version(\"bigalgebra\", version = \"1.1.1\")")

    # Update the state
    .plmmr_env$warning_shown <- TRUE
  }

  # return list for model preparation ---------------------------------
  ret <- list(
    std_X = std_X,
    std_X_details = std_X_details,
    std_X_n = std_X_n,
    std_X_p = std_X_p,
    y = y,
    y_name = colnames(design$y),
    centered_y = y - mean(y),
    K = K,
    diag_K = diag_K,
    fbm_flag = fbm_flag,
    plink_flag = plink_flag,
    penalty = penalty,
    penalty_factor = penalty_factor,
    gamma = gamma,
    init = init,
    n = design$n,
    p = design$p
  )

  if (trace && !is.null(save_rds)) {
    cat("Your results will be saved to ", save_rds, "\n")
  }

  return(ret)

}

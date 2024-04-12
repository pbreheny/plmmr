#' plmm_checks
#'
#' @param X Design matrix object or a string with the file path to a design matrix. If a string, string will be passed to `get_data()`. 
#' * Note: X may include clinical covariates and other non-SNP data, but no missing values are allowed.
#' @param y Continuous outcome vector. Defaults to NULL, assuming that the outcome is the 6th column in the .fam PLINK file data. Can also be a user-supplied numeric vector. 
#' @param std_needed Logical: does the supplied X need to be standardized? Defaults to NULL, since by default, X will be standardized internally. For data processed from PLINK files, standardization happens in `process_plink()`. For data supplied as a matrix, standardization happens here in `plmm()`. If you know your data are already standardized, leave `std_needed = FALSE` -- this would be an atypical case. **Note**: failing to standardize data will lead to incorrect analyses. 
#' @param col_names Optional vector of column names for design matrix. Defaults to NULL.
#' @param k An integer specifying the number of singular values to be used in the approximation of the rotated design matrix. This argument is passed to `RSpectra::svds()`. Defaults to `min(n, p) - 1`, where n and p are the dimensions of the _standardized_ design matrix.
#' @param K Similarity matrix used to rotate the data. This should either be (1) a known matrix that reflects the covariance of y, (2) an estimate (Default is \eqn{\frac{1}{p}(XX^T)}), or (3) a list with components 'd' and 'u', as returned by choose_k().
#' @param diag_K Logical: should K be a diagonal matrix? This would reflect observations that are unrelated, or that can be treated as unrelated. Defaults to FALSE. 
#'  Note: plmm() does not check to see if a matrix is diagonal. If you want to use a diagonal K matrix, you must set diag_K = TRUE.
#' @param eta_star Optional argument to input a specific eta term rather than estimate it from the data. If K is a known covariance matrix that is full rank, this should be 1.
#' @param penalty The penalty to be applied to the model. Either "MCP" (the default), "SCAD", or "lasso".
#' @param penalty.factor A multiplicative factor for the penalty applied to each coefficient. If supplied, penalty.factor must be a numeric vector of length equal to the number of columns of X. The purpose of penalty.factor is to apply differential penalization if some coefficients are thought to be more likely than others to be in the model. In particular, penalty.factor can be 0, in which case the coefficient is always in the model without shrinkage.
#' @param init Initial values for coefficients. Default is 0 for all columns of X. 
#' @param gamma The tuning parameter of the MCP/SCAD penalty (see details). Default is 3 for MCP and 3.7 for SCAD.
#' @param alpha Tuning parameter for the Mnet estimator which controls the relative contributions from the MCP/SCAD penalty and the ridge, or L2 penalty. alpha=1 is equivalent to MCP/SCAD penalty, while alpha=0 would be equivalent to ridge regression. However, alpha=0 is not supported; alpha may be arbitrarily small, but not exactly 0.
#' @param dfmax Upper bound for the number of nonzero coefficients. Default is no upper bound. However, for large data sets, computational burden may be heavy for models with a large number of nonzero coefficients.
#' @param trace If set to TRUE, inform the user of progress by announcing the beginning of each step of the modeling process. Default is FALSE.
#' @param ... Additional arguments to `get_data()`
#'
#' @keywords internal
#'
plmm_checks <- function(X,
                        y = NULL,
                        std_needed = NULL,
                        col_names = NULL,
                        k = NULL, 
                        K = NULL,
                        diag_K = NULL,
                        eta_star = NULL,
                        penalty = c("MCP", "SCAD", "lasso"), # TODO: think about making lasso default
                        penalty.factor = NULL,
                        init = NULL,
                        gamma,
                        alpha = 1,
                        dfmax = NULL,
                        trace = FALSE,
                        ...){
  # check X types -------------------------------------------------
  ## string with a filepath -----------------------------
  if("character" %in% class(X)){
    dat <- get_data(path = X, trace = trace, ...)
    std_X <- dat$std_X
    if("FBM.code256" %in% class(std_X) | "FBM" %in% class(std_X)){
      fbm_flag <- TRUE
    } else {
      fbm_flag <- FALSE
    }
  }
  
  ## FBM object ----------------
  if("FBM.code256" %in% class(X) | "FBM" %in% class(X)){stop("To analyze data from a file-backed X matrix, meta-data must also be supplied. 
                                       For the 'X' arg, you need to supply either (1) a character string representing the filepath to the .rds object or (2) a list as returned by get_data(). 
                                       If you don't have an .rds object yet, see process_plink() for preparing your data.")}
  ## list -----------------------
  if("list" %in% class(X)){
    # check for X element 
    if(!("std_X" %in% names(X))){stop("The list supplied for the X argument does not have a design matrix element named 'std_X'. Rename as needed and try again.")}
    if("FBM.code256" %in% class(X$std_X) | "FBM" %in% class(X$std_X)){
      dat <- X
      std_X <- dat$std_X
      # designate the dimensions of the standardized design matrix, with only ns columns
      std_X_n <- std_X$nrow
      std_X_p <- std_X$ncol
      fbm_flag <- TRUE
    } else {
      std_X <- X$std_X
      std_X_n <- nrow(std_X)
      std_X_p <- ncol(std_X)
      # TODO: add case to handle X passed as list where X is not an FBM
      fbm_flag <- FALSE
    }
    
  }
  
  ### set fbm flag ---------------------------
  # if FBM flag is not 'on' by now, set it 'off' & turn on standardization
  if (!exists('fbm_flag')){
    fbm_flag <- FALSE
  }
  
  ## matrix -------------------------------
  if (!fbm_flag){
    
    if (!inherits(X, "matrix")) {
      tmp <- try(X <- stats::model.matrix(~0+., data=X), silent=TRUE)
      if (inherits(tmp, "try-error")) stop("X must be a matrix or able to be coerced to a matrix", call.=FALSE)
    }
    if (typeof(X)=="integer") storage.mode(X) <- "double"
    if (typeof(X)=="character") stop("if X is a matrix, it must be a numeric matrix", call.=FALSE)
    
    # designate the dimensions of the original design matrix 
    n <- nrow(X)
    p <- ncol(X) 
    
    # set default column names of X
    if (is.null(col_names) & !is.null(attr(X, "dimnames")[[2]])) {
      col_names <- attr(X, "dimnames")[[2]]
    }
 
    # handle standardization ---------------------------------------------
    if (is.null(std_needed) & !fbm_flag){
      std_X <- ncvreg::std(X)
      std_X_details <- list(center = attr(std_X, 'center'),
                            scale = attr(std_X, 'scale'),
                            ns = attr(std_X, 'nonsingular'))
    } else if(std_needed){
      std_X <- ncvreg::std(X)
      std_X_details <- list(center = attr(std_X, 'center'),
                            scale = attr(std_X, 'scale'),
                            ns = attr(std_X, 'nonsingular'))
      
      # designate the dimensions of the standardized design matrix, with only ns columns
    } else if (!std_needed) {
      std_X <- X
      # TODO: may need to change this default setting 
      
      if (trace) { 
        cat("\nYou have left std_needed = FALSE; this means you are assuming that 
            all columns in X have mean 0 and variance 1. This also means 
            you are assuming that none of the columns in X are constant features. 
            \n If this is not what you want to assume, leave std_needed = NULL (default).")  
        
      }
      
      std_X_details <- list(center = rep(0, ncol(X)),
                            scale = rep(1, ncol(X)),
                            ns = 1:ncol(X))
    }
    
  }
  # designate dimensions of the 
  std_X_n <- nrow(std_X)
  std_X_p <- ncol(std_X)
  
  #  check y types -------------------------------
  # if y is null, use .fam file 
  if(is.null(y)){
    # default: uses bigSNP naming convention as would be returned in get_data()
    y <- dat$fam$affection 
  }
  
  if (!is.double(y)) {
    op <- options(warn=2)
    on.exit(options(op))
    y <- tryCatch(
      error = function(cond) stop("y must be numeric or able to be coerced to numeric", call.=FALSE),
      as.double(y))
    options(op)
  }  
  
  # set up defaults --------------------------------------------------
  if(is.null(dfmax)){dfmax <- std_X_p + 1}
  
  # default: penalize everything except the intercept, which we will add later
  if(is.null(penalty.factor)){penalty.factor <- rep(1, std_X_p)}
  
  # set default init
  if(is.null(init)){init <- rep(0, std_X_p)}
  
  # coercion for penalty
  penalty <- match.arg(penalty)
  # for now, filebacked data are limited to lasso only (until biglasso gets another upgrade)
  # if (fbm_flag & penalty != "lasso") {
  #   if (trace) cat("\nFor now, filebacked data must be analyzed with penalty=lasso 
  #                  \n (we are working to expand this). Changing the penalty to 
  #                  lasso.")
  #   
  #   penalty <- 'lasso'
  # }
  
  # set default gamma
  if (missing(gamma)) gamma <- switch(penalty, SCAD = 3.7, 3)
  
  # error checking ------------------------------------------------------------
  if(!fbm_flag){
    # error check for matrix X
    if (length(y) != std_X_n) stop("X and y do not have the same number of observations", call.=FALSE)
    if (any(is.na(y)) | any(is.na(std_X))) stop("Missing data (NA's) detected.  
                                            \nTake actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to plmm", call.=FALSE)
    if (length(penalty.factor)!=std_X_p) stop("Dimensions of penalty.factor and X do not match", call.=FALSE)
  } else {
    #  error checking for FBM X 
    if (length(y) != std_X_n) stop("X and y do not have the same number of observations", call.=FALSE)
    if (any(is.na(y))) stop("Missing data (NA's) detected in the outcome.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to plmm", call.=FALSE)
    if (length(penalty.factor)!=std_X_p) stop("Dimensions of penalty.factor and X do not match", call.=FALSE)
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

  # create a list that captures the centering/scaling for std_X; will need this 
  # later, see `untransform()`
  if(fbm_flag){
    std_X_details <- list(
      center = dat$std_X_center,
      scale = dat$std_X_scale,
      ns = dat$ns
    )
  } 
  
  
  # return list for model preparation ---------------------------------
  
  ret <- list(
    std_X = std_X,
    std_X_details = std_X_details,
    std_X_n = std_X_n,
    std_X_p = std_X_p,
    col_names = col_names,
    y = y,
    K = K,
    diag_K = diag_K,
    fbm_flag = fbm_flag,
    penalty = penalty,
    penalty.factor = penalty.factor,
    gamma = gamma,
    init = init
  )
  
  if (exists('dat')){
    ret$n <- dat$n
    ret$p <- dat$p
    ret$dat <- dat
  } else {
    ret$n <- n
    ret$p <- p
  }
  
  return(ret)
  
}
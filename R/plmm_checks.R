#' plmm_checks
#'
#' @param X 
#' @param fbm 
#' @param std_needed 
#' @param col_names 
#' @param y 
#' @param k 
#' @param K 
#' @param diag_K 
#' @param eta_star 
#' @param penalty 
#' @param penalty.factor 
#' @param init 
#' @param gamma 
#' @param alpha 
#' @param dfmax 
#' @param trace 
#'
#' @return
#' @keywords internal
#'
plmm_checks <- function(X,
                        fbm = NULL,
                        std_needed = NULL,
                        col_names = NULL,
                        y = NULL,
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
                        trace = FALSE){
  
  # check X types -------------------------------------------------
  ## string with a filepath -----------------------------
  if("character" %in% class(X)){
    if(is.null(fbm)){
      dat <- get_data(path = X, trace = trace)
    } else if (fbm){
      dat <- get_data(path = X, returnX = FALSE, trace = trace)
    } else if (!fbm){
      dat <- get_data(path = X, returnX = TRUE, trace = trace)
    }
    
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
      std_X <- obj$std_X
      std_X_n <- std_X$nrow
      std_X_p <- std_X$ncol
      # TODO: add case to handle X passed as list where X is not an FBM
      fbm_flag <- FALSE
    }
    
  }
  
  ### set fbm flag ---------------------------
  # if FBM flag is not 'on' by now, set it 'off' & turn on standardization
  if (!exists('fbm_flag')){
    fbm_flag <- FALSE
    std_needed <- TRUE
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
    
    # handle standardization ---------------------------------------------
    if(std_needed){
      std_X <- ncvreg::std(X)
      std_X_details <- list(center = attr(std_X, 'center'),
                            scale = attr(std_X, 'scale'),
                            ns = attr(std_X, 'nonsingular'))
      
      # designate the dimensions of the standardized design matrix, with only ns columns
    } else if (is.null(std_needed) & !fbm_flag){
      std_X <- ncvreg::std(X)
      std_X_details <- list(center = attr(std_X, 'center'),
                            scale = attr(std_X, 'scale'),
                            ns = attr(std_X, 'nonsingular'))
    } else if (!std_needed) {
      std_X <- X
      # TODO: may need to change this default setting 
      
      if (trace) { 
        cat("\nYou have set std_needed = FALSE; this means you are assuming that 
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
  if (fbm_flag & penalty != "lasso") {
    if (trace) cat("\nFor now, filebacked data must be analyzed with penalty=lasso 
                   \n (we are working to expand this). Changing the penalty to 
                   lasso.")
    
    penalty <- 'lasso'
  }
  
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
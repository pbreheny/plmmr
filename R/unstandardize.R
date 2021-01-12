#' Unstandardize coefficient values
#'
#' This function allows you to unstandardize coefficient values based on attributes of X.
#' @param b p x nlambda matrix of standardized coefficient path values.
#' @param X Design matrix which *includes* the intercept column if present. May include clinical covariates and other non-SNP data. Must include the attributes 'center', 'scale', and 'nonsingular' corresponding to non-intercept variables.
#' @param intercept Logical flag for whether an intercept row is included in b, that is, a coefficient which accounts for the mean of the original y value on its rotated scale.
#' @param centerY Logical flag for whether the outcome variable, y, should be centered and an intercept coefficient estimated based on its mean. Defaults to FALSE.
#' @param a Mean out outcome variable to be used in unstandardization. Required if centerY = TRUE. Defaults to NULL.
#' @export

unstandardize <- function(b, X, intercept = TRUE, centerY = FALSE, a = NULL) {
  if (intercept & centerY) stop('It looks like you are fitting the intercept twice', call.=FALSE)
  if (centerY & is.null(a)) stop('Must supply a value if centerY = TRUE', call. = FALSE)
  ns <- attr(X, 'nonsingular')
  if (intercept){
    bbb <- matrix(0, nrow = nrow(b), ncol = ncol(b))
    a <- b[1, , drop = FALSE]
    bb <- b[-1, , drop=FALSE]
    bbb[1 + ns,] <- bb / attr(X, 'scale')[ns]
    bbb[1,] <- a - crossprod(attr(X, 'center')[ns], bb)
  } else if (centerY){
    bbb <- matrix(0, nrow = nrow(b) + 1, ncol = ncol(b))
    a <- a
    bb <- b[, , drop=FALSE]
    bbb[1 + ns,] <- bb / attr(X, 'scale')[ns]
    bbb[1,] <- a - crossprod(attr(X, 'center')[ns], bb)
  } else {
    bbb <- matrix(0, nrow = nrow(b), ncol = ncol(b))
    # a better logical check would be better...
    if (nrow(b) == ncol(X)){
      bbb[ns, ] <- b/attr(X, 'scale')[ns]
    } else {
    # this means that the intercept was already calculated using the rotated outcome mean
    a <- b[1, , drop = FALSE]
    bb <- b[-1, , drop=FALSE]
    bbb[1 + ns,] <- bb / attr(X, 'scale')[ns]
    bbb[1,] <- a - crossprod(attr(X, 'center')[ns], bb)
    }
  }
  return(bbb)
}

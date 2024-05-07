#' Extract residuals from a PLMM fit
#' 
#' Currently, only deviance residuals are supported.
#' 
#' @param object   Object of class `plmm`.
#' @param lambda   Optional numeric value(s) of the regularization parameter at which residuals are requested (numeric vector). For values of lambda not in the sequence of fitted models, linear interpolation is used.
#' @param which    Index of the penalty parameter at which residuals are requested (default = all indices). If `lambda` is specified, this take precedence over `which`.
#' @param drop     By default, if a single value of lambda is supplied, a vector of residuals is returned (logical; default=`TRUE`). Set `drop=FALSE` if you wish to have the function always return a matrix (see [drop()]).
#' @param unrotate Logical: should residuals be 'unrotated', i.e. transformed back to the original scale? 
#' @param ...      Not used.
#' 
#' @rdname residuals.plmm
#' 
#' @returns a numeric matrix with the residuals on the **transformed scale**, where
#' rows correspond to observations (e.g., samples) and columns correspond to 
#' values of `lambda` at which the model was fit. 
#' 
#' @export
#' 
#' @examples
#' \donttest{
#' fit <- plmm(admix$X, admix$y)
#' residuals.plmm(fit)[1:5, 1:5]
#' head(residuals.plmm(fit, which = 50))
#' }
#' 
#' @returns 
#'
#' 
#' 

residuals.plmm <- function(object,
                           lambda,
                           which=1:length(object$lambda),
                           drop=TRUE, 
                           unrotate = FALSE,
                           ...) {
  cat("\nNote: Residuals are returned on the transformed scale 
      \n(i.e., the scale on which the model was fit)")
  # Calculate matrix of residuals
  R <- matrix(nrow = nrow(object$linear.predictors), ncol = ncol(object$linear.predictors))
  for(j in 1:ncol(R)){
    R[,j] <- object$rot_y - object$linear.predictors[j]
  }
  if(unrotate){
    stop("\nThis option is still underdevelopment. Don't try to unrotate residuals yet.")
    # w <- (object$eta * object$s + (1 - object$eta))^(-1/2)
    # W_inv <- sweep(object$U, 2, 1/(w), "*")
    # WR <- W_inv%*%R
    # TODO: if this is worth pursuing, I need to finish working this through

  } 
  
  # Interpolate and return
  if (!missing(lambda)) {
    ind <- approx(object$lambda, seq(object$lambda), lambda)$y
    l <- floor(ind)
    r <- ceiling(ind)
    q <- ind %% 1
    out <- (1-q)*R[, l, drop=FALSE] + q*R[, r, drop=FALSE]
    colnames(out) <- round(lambda, 4)
  } else {
    out <- R[, which, drop=FALSE]
  }
  if (drop) return(drop(out)) else return(out)
}

# from ncvreg
coef.plmm <- function (object, lambda, which = 1:length(object$lambda), drop = TRUE,  ...){
  if (!missing(lambda)) {
    if (max(lambda) > max(object$lambda) | min(lambda) <
        min(object$lambda)) {
      stop("Supplied lambda value(s) are outside the range of the model fit.",
           call. = FALSE)
    }
    ind <- approx(object$lambda, seq(object$lambda), lambda)$y
    l <- floor(ind)
    r <- ceiling(ind)
    w <- ind%%1
    beta <- (1 - w) * object$beta[, l, drop = FALSE] + w *
      object$beta[, r, drop = FALSE]
    colnames(beta) <- lamNames(lambda)
  }
  else beta <- object$beta[, which, drop = FALSE]
  if (drop)
    return(drop(beta))
  else return(beta)
}

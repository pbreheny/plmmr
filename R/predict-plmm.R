## from ncvreg
predict.plmm <- function(object, X, type=c("response", "coefficients", "vars", "nvars"),
                           lambda, which=1:length(object$lambda), ...) {
  type <- match.arg(type)
  beta <- coef.plmm(object, lambda=lambda, which=which, drop=FALSE)
  if (type=="coefficients") return(beta)
  alpha <- beta[1,]
  beta <- beta[-1, , drop=FALSE]
  if (type=="nvars") return(apply(beta!=0, 2, sum))
  if (type=="vars") return(drop(apply(beta!=0, 2, FUN=which)))
  eta <- sweep(X %*% beta, 2, alpha, "+")
  if (type=="link" || object$family=="gaussian") return(drop(eta))
}

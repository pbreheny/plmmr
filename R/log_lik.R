#' Evaluate the negative log-likelihood of an intercept-only Gaussian plmm model
#'
#' This function allows you to evaluate the negative log-likelihood of a linear mixed model under the assumption of a null model in order to estimate the variance parameter, eta.
#'
#' @param eta Estimated proportion of the variance in the outcome attributable to population/correlation structure
#' @param n The number of observations
#' @param s The non-zero eigenvalues of K, the realized relationship matrix
#' @param U The eigenvectors of K associated with s
#' @param y Continuous outcome vector
#' @param incpt_flag Logical: Does the model require fitting an intercept? Passed from `estimate_eta`.
#'
#' @return the value of the log-likelihood of the PLMM, evaluated with the supplied parameters
#'
#' @keywords internal
#'
log_lik <- function(eta, n, s, U, y, incpt_flag) {
  # first, the constant (comes from 1st term in derivation)
  constant <- n * log(2 * pi)

  # we will need the sum of the nonzero values from the diagonal matrix of weights
  w2 <- ((eta * s) + (1 - eta))

  # get w2 on the log scale
  sum_det_log <- sum(log(w2))

  # rotate y
  w <- w2^(-1 / 2)
  wUt <- sweep(x = t(U), MARGIN = 1, STATS = w, FUN = "*")

  rot_y <- U %*% wUt %*% y

  if (incpt_flag) {
    incpt <- rep(1, n)
    rot_incpt <- U %*% wUt %*% incpt
    # calculate the term representing the \hat\beta(\eta) MLE
    incpt_crossprod <- crossprod(rot_incpt)
    incpt_y_crossprod <- crossprod(rot_incpt, rot_y)
    hat_beta_mle <- drop(incpt_y_crossprod / incpt_crossprod)
    # using hat_beta_mle, calculate the quadratic term from the log likelihood
    rot_e <- rot_y - (rot_incpt * hat_beta_mle)
  } else {
    rot_e <- rot_y
  }

  # rotated intercept is the only term in the null model
  rot_sqe <- crossprod(rot_e) # mse = sq. error
  quad_term <- (1 / n) * rot_sqe / sum(w2)

  # Return the negative log-likelihood
  nLL <- 0.5 * (constant + n * log(quad_term) + sum_det_log + n)
  drop(nLL)
}

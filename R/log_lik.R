#' Evaluate the negative log-likelihood of a null Gaussian penalizedLMM model
#'
#' This function allows you to evaluate the negtive log-likelihood of a linear mixed model under the assumption of a null model in order to estimate the variance parameter, eta.
#' @param eta The proportion of variance in the outcome that is attributable to causal SNP effects. In other words, SNR. Sometimes referred to as the narrow-sense heritability.
#' @param Uy The the continuous outcome, y, rotated by the eigenvectors of the similarity matrix, K.
#' @param S The eigenvalues of the similarity matrix, K.
#' @keywords internal
#' 
#' @examples 
#' \dontrun{
#' admix$K <- relatedness_mat(admix$X) # create an estimated covariance matrix 
#' ev <- eigen(admix$K)
#' U <- ev$vectors
#' fit <- plmm(X = admix$X, y = admix$y, K = admix$K)
#' (log_lik(eta = fit$eta, Uy = U%*%admix$y, S = ev$values ))
#' }

log_lik <- function(eta, Uy, S){

  n <- dim(Uy)[1]

  # evaluate log determinant
  Sd <- eta * S + (1 - eta)
  ldet <- sum(log(Sd))  # log of product = sum of the logs

  # evaluate the variance
  Sdi <- 1/Sd 
  Uy <- as.vector(Uy)
  ss <- (1/n) * sum(Uy*Uy*Sdi)

  # evaluate the negative log likelihood
  # NB: keep constant here to be consistent with log_lik.lm() method
  nLL <- 0.5*(n*log(2*pi) + ldet + n + n*log(ss))
  # TODO: double check the derivation for this. Do we need the factor of n 
  #   in the last term? 

  return(nLL)

}
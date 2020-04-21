
#' Generate data with population structure
#'
#' This function allows you to simulate structured genetic data (SNP) with an unobserved environmental confounding effect.
#' @param n Number of observations/samples to simulate.
#' @param p Number of SNPs to simulate.
#' @param p1 Number of SNPs that are causal. Defaults to floor(p/2).
#' @param nJ Number of observations in each subpopulation. The length of nJ corresponds to the number of subpopulations.
#' @param structureX Type of structure to simulate.
#' @param Fst The desired final FST of the admixed individuals. Ranges from 0 to 1. A high Fst implies greater differentiation among populations. Defaults to 0.1 if structureX = 1d_linear and 0.2 if structureX = indep_subpops. Otherwise defaults to NULL.
#' @param inbr Indicates whether the desired inbreeding is homogeneous or heterogeneous. Defaults to heterogeneous.
#' @param structureGamma The desired structure of the environmental confounding effect. Defaults to halfandhalf_decreasing_heterogeneous.
#' @param eta The desired proportion of variance in the outcome that is attributable to causal SNP effects. In other words, SNR.
#' @param xi The desired proportion of the non-signal variance in the outcome that is attributable to unobserved environmental confounding effects.
#' @param standardizeX Should the generated X matrix be standardized? Defaults to TRUE.
#' @keywords
#' @export
#' @examples

genDataPS <- function(n = 197, p = 98, p1 = floor(p/2), nJ = c(47, 50, 50, 50),
                      structureX = c("admixture", "indep_subpops", "1d_linear", "independent"),
                      Fst = NULL,
                      inbr = c("heterogeneous", "homogeneous"),
                      structureGamma = c("halfandhalf_decreasing_heterogeneous"),
                      eta, xi, standardizeX = TRUE){

  structureX <- match.arg(structureX)
  J <- length(nJ)

  # dimension checks
  if (structureX == "admixture"){
    if ((n != 197 || p != 98 || !missing(nJ) & length(nJ) != 4)){
      warning("Admixture structureX is specified which will override supplied values of n, p, J, and nJ.")
    }
    n <- 197
    p <- 98
    J <- 4
    # Admixture row order: AA (47), EU (50), JP (50), AF (50)
    nJ <- c(47, 50, 50, 50)
  }

  if (missing(nJ)) nJ <- c(rep(n%/%J, J-1), n%/%J + n%%J)
  if (!missing(nJ) & sum(nJ) != n) stop("sum(nJ) must equal n")

  # Generate appropriate scaling factors depending on if partitioning var(y) or not
  scaleXbeta <- sqrt(eta)
  scaleZgamma <- sqrt((1 - eta) * xi)
  scaleEps <- sqrt((1 - eta) * (1 - xi))

  # Gen X
  X <- genXps(n, nJ, p, structureX, Fst, inbr, standardizeX)

  # Gen (scaled) beta and standardized Xbeta (scale by eta)
  j <- 1:p
  s <- rep(1, times = p)
  b <- (j <= p1)
  b <- b * s
  covX <- var(X) * (n - 1)/n
  beta <- b * scaleXbeta / sqrt(drop(crossprod(b, covX) %*% b)) # this scales Xbeta
  Xbeta <- X%*%beta
  # Xbeta <- X%*%beta - mean(X%*%beta) # this centers Xbeta to mean 0
  # older way...
  # Xbeta <- X%*%b
  # Xbeta_sd <- sqrt(varp(Xbeta))
  # Xbeta <- Xbeta / Xbeta_sd * scaleXbeta
  # beta <- b / Xbeta_sd * scaleXbeta

  # Gen Z
  mlist <- vector("list", J)
  for (jj in 1:J){
    mlist[[jj]] <- matrix(rep(1, nJ[jj]), ncol = 1)
  }
  Z <- as.matrix(Matrix::bdiag(mlist))

  # Gen (scaled) gamma and standardized Zgamma (scale by eta and xi)
  if (!is.numeric(structureGamma) && structureGamma == "banded"){
    gamma_unscaled <- (1:n)/n
    Z <- diag(n)
  } else {
    gamma_unscaled <- genGammaUnscaled(structureGamma, J)
  }
  env <- Z %*% gamma_unscaled
  # env <- drop(as.matrix(std(env))) * scaleZgamma
  env_sd <- sqrt(varp(env))
  env <- env/env_sd * scaleZgamma
  # gamma <- (gamma_unscaled - attributes(env)$center)/attributes(env)$scale * scaleZgamma
  gamma <- gamma_unscaled/env_sd * scaleZgamma
  # Gen Y
  y <- Xbeta + env + scaleEps * drop(std(as.matrix(rnorm(n))))

  # name things
  w <- 1 + floor(log10(p))
  vlab <- paste0("V", formatC(1:p, format = "d", width = w, flag = "0"))
  colnames(X) <- names(beta) <- vlab

  # return
  list(X = X, y = y, env = env, beta = beta, Z = Z, gamma = gamma, structureX = structureX)

}



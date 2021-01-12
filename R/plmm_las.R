#' Fit a penalizedLMM with lasso regularization using penalizedLMM using various intercept and standardization options
#'
#' This function allows you to fit a linear mixed model via penalized maximum likelihood with null model variance component estimation.
#' @param X Design matrix.
#' @param y Continuous outcome vector.
#' @param X_for_K X matrix used to compute the similarity matrix, K. For multi-chromosome analysis this may be supplied in order to perform a leave-one-chromosome-out correction. The objective here is to adjust for population stratification and unobserved confounding without rotating out the causal SNP effects.
#' @param p1 Number of causal SNPs. Lambda will be selected such that <= p1 variables enter the model.
#' @param intercept Logical flag for whether an intercept should be included.
#' @param centerY Logical flag for whether the outcome variable, y, should be centered and an intercept coefficient estimated based on its mean. Defaults to FALSE.
#' @param centerRtY Logical flag for whether the rotated outcome variable, SUy, should be centered and an intercept coefficient estimated based on its mean. Defaults to FALSE.
#' @param standardizeX Logical flag for X variable standardization, prior to data transformation. The coefficients are always returned on the original scale. Default is TRUE. If variables are in the same units already, or manually standardized, you might not wish to standardize.
#' @param standardizeRtX Logical flag for transformed X variable standardization. The coefficients are always returned on the original scale. Default is TRUE. If variables are in the same units already, or manually standardized, you might not wish to standardize, but this is generally recommended.
#' @param standardizeX Should standardization be performed on the original X matrix? Defaults to FALSE.
#' @param standardizeRtX Should standardization be performed on the rotated X matrix? Defaults to FALSE.
#' @importFrom zeallot %<-%
#' @export



plmm_las <- function(X, y, X_for_K = X, p1, intercept = TRUE,
                       centerY = FALSE, centerRtY = FALSE,
                       standardizeX = FALSE, standardizeRtX = FALSE) {
  fit <- plmm(X, y, X_for_K, penalty = 'lasso',
              intercept = intercept,
              centerY = centerY, centerRtY = centerRtY,
              standardizeX = standardizeX, standardizeRtX = standardizeRtX)
  sel <- predict.plmm(fit, type = "nvar")
  coef <- coef(fit, min(fit$lambda[sel <= p1])) # this returns the intercept
  names(coef) <- colnames(X)
  return(list(fit = fit,
              nonzero = length(which(coef != 0)) - 1,
              coef = coef,
              delta = (1/fit$eta) - 1,
              eta = fit$eta))
}



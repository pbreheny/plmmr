
#' Generate and X matrix with population structure
#'
#' This function allows you to simulate structured genetic data (SNP).
#' @param n Number of observations/samples to simulate.
#' @param nJ Number of observations in each subpopulation. The length of nJ corresponds to the number of subpopulations.
#' @param p Number of SNPs to simulate.
#' @param structureX Type of structure to simulate.
#' @param Fst The desired final FST of the admixed individuals. Ranges from 0 to 1. A high Fst implies greater differentiation among populations. Defaults to 0.1 if structureX = 1d_linear and 0.2 if structureX = indep_subpops. Otherwise defaults to NULL.
#' @param inbr Indicates whether the desired inbreeding is homogeneous or heterogeneous. Defaults to heterogeneous.
#' @param standardizeX Should the generated X matrix be standardized? Defaults to TRUE.
#' @param plot Should a plot of the kinship matrix be generated? Defaults to FALSE.
#' @export

genXps <- function(n, nJ, p,
                   structureX = c("admixture", "indep_subpops", "1d_linear", "independent"),
                   Fst = NULL,
                   inbr = c("homogeneous", "heterogeneous"),
                   standardizeX = TRUE, plot = FALSE){
  structureX <- match.arg(structureX)
  if (structureX == "admixture"){
    dat <- utils::read.delim("https://s3.amazonaws.com/pbreheny-data-sets/admixture.txt")
    XX <- as.matrix(dat[,-1])
    polymorphic <- apply(XX, 2, stats::sd) != 0
    X <- XX[, polymorphic]
  } else if (structureX == "independent"){
    pi <- 0.2
    X <- matrix(stats::rbinom(n*p, 2, pi), nrow = n, ncol = p)
  } else if (structureX == "indep_subpops"){
    subpops <- rep(c(paste0("S", 1:length(nJ))), times = nJ)
    admix_proportions <- bnpsd::admix_prop_indep_subpops(subpops)
    if (is.null(Fst)){
      Fst <- 0.2
    }
    if (inbr == "heterogeneous"){
      inbr_subpops <- 1:length(nJ)
    } else {
      inbr_subpops <- rep(1, length(nJ))
    }
    inbr_subpops <- inbr_subpops / popkin::fst(inbr_subpops) * Fst
    coancestry <- bnpsd::coanc_admix(admix_proportions, inbr_subpops)
    weights <- popkin::weights_subpops(subpops) # function from `popkin` package
    X <- t(bnpsd::draw_all_admix(admix_proportions = admix_proportions,
                          inbr_subpops = inbr_subpops,
                          m_loci = p)$X)
    ### plot things
    if (plot == TRUE) {
      # fname <- paste0("n=", n, "-p=", p, "-nJ=", length(nJ), "-structureX=", structureX, "-Fst=", Fst, "-inbr=", inbr, ".pdf")
      # pdf(file = fname, width = 15, height =15)
      # col_subpops <- brewer.pal(length(nJ), "Paired")
      # par(mar = c(2.5, 2.5, 0, 0) + 1, lab = c(2, 1, 7), mgp = c(1.5, 0.5, 0))
      # barplot(inbr_subpops, col = col_subpops, names.arg = colnames(admix_proportions), xlab = 'Subpopulation', ylab = 'Inbr')
      # par(mar = c(2.5, 2.5, 0, 0) + 1, lab = c(2, 2, 7))
      # barplot(t(admix_proportions), col = col_subpops, border = NA, space = 0, ylab = 'Admix prop')
      # mtext('Individuals', 1)
      popkin::plot_popkin(coancestry, mar = c(1, 1, 1, 3))
      # dev.off()
    }

  } else if (structureX == "1d_linear"){
    subpops <- rep(c(paste0("S", 1:length(nJ))), times = nJ)
    if (inbr == "heterogeneous"){
      inbr_subpops <- 1:length(nJ)
    } else {
      inbr_subpops <- rep(1, length(nJ))
    }
    if (is.null(Fst)){
      Fst <- 0.1
    }
    # making this smaller makes blocks more distinct
    # bigger means more blurred together and indistinguishable
    # 0.35 - 0.5 seems reasonable
    bias_coeff <- 0.5
    tmp <- bnpsd::admix_prop_1d_linear(n_ind = n,
                                k_subpops = length(nJ),
                                bias_coeff = bias_coeff,
                                coanc_subpops = inbr_subpops,
                                fst = Fst)
    admix_proportions <- tmp$admix_proportions
    inbr_subpops <- tmp$coanc_subpops
    coancestry <- bnpsd::coanc_admix(admix_proportions, inbr_subpops)
    X <- t(bnpsd::draw_all_admix(admix_proportions = admix_proportions,
                          inbr_subpops = inbr_subpops,
                          m_loci = p)$X)
    ### plot things
    if (plot == TRUE) {
      # fname <- paste0("n=", n, "-p=", p, "-nJ=", length(nJ), "-structureX=", structureX, "-Fst=", Fst, "-inbr=", inbr, ".pdf")
      # pdf(file = fname, width = 15, height =15)
      # col_subpops <- brewer.pal(length(nJ), "Paired")
      # par(mar = c(2.5, 2.5, 0, 0) + 1, lab = c(2, 1, 7), mgp = c(1.5, 0.5, 0))
      # barplot(inbr_subpops, col = col_subpops, names.arg = colnames(admix_proportions), xlab = 'Subpopulation', ylab = 'Inbr')
      # par(mar = c(2.5, 2.5, 0, 0) + 1, lab = c(2, 2, 7))
      # barplot(t(admix_proportions), col = col_subpops, border = NA, space = 0, ylab = 'Admix prop')
      # mtext('Individuals', 1)
      popkin::plot_popkin(coancestry, mar = c(1, 1, 1, 3))
      # dev.off()
    }
  }

  # randomize column order so there's no LD
  X <- X[, sample(1:ncol(X), ncol(X), replace = FALSE)]
  # make numeric, not integer
  X <- apply(X, 2, as.numeric)
  # standardize
  if (standardizeX == TRUE) X <- ncvreg::std(X)
  return(X)
}

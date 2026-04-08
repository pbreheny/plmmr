#' Admix: Semi-simulated SNP data
#'
#' A dataset containing the 100 SNPs, a demographic variable representing ancestry,
#' and a simulated outcome.
#' @format A list with 3 components:
#' \describe{
#'   \item{X}{SNP matrix (197 observations of 100 SNPs)}
#'   \item{y}{197 x 1 matrix of simulated (continuous) outcomes}
#'   \item{ancestry}{vector with ancestry categorization: 0 = African, 1 = African American, 2 = European, 3 = Japanese}
#' }
#' @source \url{https://hastie.su.domains/CASI/}
"admix"

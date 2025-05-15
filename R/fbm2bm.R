#' Functions to convert between FBM and big.matrix type objects
#'
#' @param fbm An FBM object; see `bigstatsr::FBM()` for details
#' @param desc Logical: is the descriptor file desired (as opposed to the filebacked big matrix)? Defaults to FALSE.
#'
#' @return a `big.matrix` - see `bigmemory::filebacked.big.matrix()` for details
#'
#' @keywords internal
#'
fbm2bm <- function(fbm, desc = FALSE) {
  if (desc) {
    desc <- fbm$bm.desc()
    return(bigmemory::attach.big.matrix(desc))

  } else {
    return(fbm$bm())
  }

}

#' Functions to convert between FBM and big.matrix type objects
#'
#' @param fbm An FBM object; see `bigstatsr::FBM()` for details
#' @param desc Logical: is the descriptor file desired (as opposed to the filebacked big matrix)? Defaults to FALSE.
#'
#' @return a `big.matrix` - see `bigmemory::filebacked.big.matrix()` for details
#'
#' @export
#'
#' @examples
#' colon2 <- read.delim(get_example_data("colon2.txt"), header = TRUE)
#' class(colon2)
#' bigstatsr::as_FBM(colon2) -> colon2_fbm
#' fbm2bm(colon2_fbm) -> colon2_bm
#' class(colon2_bm)
#'
fbm2bm <- function(fbm, desc = FALSE){
  if(desc){
    desc <- fbm$bm.desc()
    return(bigmemory::attach.big.matrix(desc))

  } else {
    return(fbm$bm())
  }

}

#' Create an FBM that points to the same backingfile as a big.matrix
#' Credit to Florian Prive for this function
#' @param bm A `big.matrix` object
#'
#' @return An `FBM` object pointing to the same backingfile
#'
#' @keywords internal
#' @details
#' For an extended example, see the tutorial on the `bigstatsr` package site:
#' https://privefl.github.io/bigstatsr/articles/bigstatsr-and-bigmemory.html
#'
bm2fbm <- function(bm){
  warning("\nThis function is under construction - don't rely on this yet!")
  bigstatsr::FBM(nrow = nrow(bm), ncol = ncol(bm), type = typeof(bm),
      backingfile = file.path(dir.name(bm), bigstatsr::sub_bk(file.name(bm))),
      create_bk = FALSE)
}


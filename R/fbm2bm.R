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
#' colon2 <- read.delim(find_example_data("colon2.txt"), header = TRUE)
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

bm2fbm <- function(bm){
  f <- paste0(bigmemory::dir.name(bm), bigmemory::file.name(bm))
  bk_file_path_sans_extension <- bigstatsr::sub_bk(f)
  fbm <- bigstatsr::FBM(nrow = nrow(bm), ncol = ncol(bm),
                          create_bk = FALSE,
                          backingfile = bk_file_path_sans_extension)
  return(fbm)
}


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
#' \dontrun{
#' pen <- get_data("../temp_files/penncath_lite", fbm = T)
#' pen_bm <- fbm2bm(pen$std_X)
#' fit <- biglasso::biglasso(X = pen_bm, y = pen$fam$affection, family = "gaussian")
#' }
fbm2bm <- function(fbm, desc = FALSE){
  if(desc){
    desc <- fbm$bm.desc()
    return(bigmemory::attach.big.matrix(desc))
    
  } else {
    return(fbm$bm())
  }

}


# bm2fbm <- function(bm, type = "integer"){
#   warning("\nThis function is under construction")
#   # Credit to Florian Prive for this function - see https://privefl.github.io/bigstatsr/articles/bigstatsr-and-bigmemory.html
#   fbm <- bigstatsr::FBM(nrow = nrow(bm),
#                         ncol = ncol(bm), 
#                         type = type,
#                         backingfile = file.path(dir.name(bm),
#                                                 bigstatsr::sub_bk(file.name(bm))),
#                         create_bk = F)
#   
#   return(fbm)
# }

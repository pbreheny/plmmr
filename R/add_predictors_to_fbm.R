#' Title
#'
#' @param X 
#' @param add_predictor_ext 
#' @param id_var 
#' @param rownames_X 
#' @param non_gen 
#' @param quiet 
#'
#' @return
#' @export
#'
add_predictors_to_fbm <- function(X, add_predictor_ext, id_var, rownames_X, non_gen, quiet){
  
  if (!is.null(add_predictor_ext)) {
    if (!quiet) {
      cat("\nAdding predictors from external data.")
    }
    if (is.vector(add_predictor_ext)) {
      ### vector case -------------------------------
      # make sure types match
      if (!is.numeric(add_predictor_ext)) {
        stop("\nThe vector supplied to the 'add_predictor_ext' argument must be numeric.")
      }
      names(add_predictor_ext) <- as.numeric(names(add_predictor_ext))
      
      if (var(add_predictor_ext) == 0) {
        stop("\nThe supplied argument to add_predictor_ext is constant (no variation).
             This would not be a meaningful predictor.")
      }
      
      # check for alignment 
      if (is.null(names(add_predictor_ext)) | 
          length(intersect(rownames_X, names(add_predictor_ext))) == 0) {
        stop("\nYou supplied an argument to 'add_predictor_ext', but the names of this
         vector either (a) do not exist or (b) do not align with either of the ID columns in the PLINK fam file.
         \nPlease create or align the names of this vector - alignment is essential for accurate analysis.")
      }
      browser()
      # TODO: pickup here
      add_predictor_ext <- align_ids(rownames_X = rownames_X,
                                     id_var = id_var,
                                     quiet = quiet,
                                     add_predictor = add_predictor_ext)
      
      # save non_gen: an index marking the first column as non-genomic predictor
      non_gen <- 1
      
      obj$geno_plus_predictors <- bigstatsr::FBM(init = 0,
                                                 nrow = nrow(obj$fam),
                                                 ncol = obj$genotypes$ncol + length(non_gen)) 
      # fill in new matrix 
      bigstatsr::big_apply(obj$genotypes,
                           a.FUN = function(X, ind, res){
                             res[,1:length(non_gen)] <- add_predictor_ext
                             res[,ind+length(non_gen)] <- X[,ind]
                           },
                           a.combine = cbind,
                           res = obj$geno_plus_predictors)
      
      # adjust colnames
      obj$colnames <- c(deparse(substitute(add_predictor_ext)), obj$colnames)
      
    } else if (is.matrix(add_predictor_ext) | is.data.frame(add_predictor_ext)) {
      ### matrix case --------------------------------
      if (is.data.frame(add_predictor_ext)) {
        add_predictor_ext <- as.matrix(add_predictor_ext)
      }
      # make sure types match
      if (!is.numeric(add_predictor_ext[,1])) {
        stop("\nThe matrix supplied to the 'add_predictor_ext' argument must have numeric values only.")
      }
      
      if (any(apply(add_predictor_ext, 2, var) == 0)) {
        stop("\nThe matrix supplied to the 'add_predictor_ext' argument has at least one
             constant column (a column that does not vary over the given samples).")
      }
      
      # check for alignment 
      if (is.null(rownames(add_predictor_ext)) | 
          length(intersect(rownames_X, rownames(add_predictor_ext))) == 0) {
        stop("\nYou supplied an argument to 'add_predictor_ext', but the row names of this
         matrix either (a) do not exist or (b) do not align with either of the ID columns in the PLINK fam file.
         \nPlease create or align the names of this matrix - alignment is essential for accurate analysis.")
      }
      
      add_predictor_ext <- align_famfile_ids(id_var = id_var,
                                             quiet = quiet,
                                             add_predictor = add_predictor_ext,
                                             rownames_X = rownames_X)
      
      # save non_gen: an index marking added columns as non-genomic predictors
      non_gen <- 1:ncol(add_predictor_ext)
      
      obj$geno_plus_predictors <- bigstatsr::FBM(init = 0,
                                                 nrow = nrow(obj$fam),
                                                 ncol = obj$genotypes$ncol + length(non_gen)) 
      # fill in new matrix 
      bigstatsr::big_apply(obj$genotypes,
                           a.FUN = function(X, ind, res){
                             res[,1:length(non_gen)] <- add_predictor_ext
                             res[,ind+length(non_gen)] <- X[,ind]
                           },
                           a.combine = cbind,
                           res = obj$geno_plus_predictors)
      
      # adjust colnames if applicable 
      if (!is.null(colnames(add_predictor_ext))){
        obj$colnames <- c(colnames(add_predictor_ext), obj$colnames)
      }
      
    }
    
  }
  
  return(list(obj = obj, non_gen = non_gen))
  
}
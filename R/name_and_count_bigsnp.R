#' A helper function to label and summarize the contents of a `bigSNP`
#'
#' @param obj a `bigSNP` object, possibly subset by `add_external_phenotype()`
#' @param id_var String specifying which column of the PLINK `.fam` file has the unique sample identifiers. Options are "IID" (default) and "FID".
#' @param quiet Logical: should messages be printed to the console? Defaults to TRUE
#'
#' @return a list with components:
#' * counts: column-wise summary of the minor allele counts in 'genotypes'
#' * obj: a modified `bigSNP` list with additional components
#' * X: the `obj$genotypes` as its own FBM
#' * pos: the `obj$map$physical.pos` vector
#' @keywords internal

name_and_count_bigsnp <- function(obj, id_var, quiet){

  X <- obj$genotypes

  # set object names
  obj$colnames <- obj$map$marker.ID

  if (id_var == "FID"){
    obj$rownames <- og_plink_ids <- as.character(obj$fam$family.ID)
  } else if (id_var == "IID") {
    obj$rownames <- og_plink_ids <- as.character(obj$fam$sample.ID)
  }

  chr_range <- range(obj$map$chromosome)
  # save the dimensions of the *original* (pre-standardized) design matrix
  obj$n <- X$nrow
  obj$p <- X$ncol

  if(!quiet){
    cat("\nThere are", obj$n, "observations and",
        obj$p, "genomic features in the specified data files, representing chromosomes",
        chr_range[1], "-", chr_range[2])
  }

  # save these counts
  counts <- bigstatsr::big_counts(X) # NB: this is a matrix

  return(list(counts = counts,
              obj = obj,
              og_plink_ids = og_plink_ids,
              chr = obj$map$chromosome,
              X = X,
              pos = obj$map$physical.pos,
              chr_range = chr_range))

}
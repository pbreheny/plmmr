#' A helper function to label and summarize the contents of a `bigSNP`
#'
#' @param obj a `bigSNP` object, possibly subset by `add_external_phenotype()`
#' @param id_var String specifying which column of the PLINK `.fam` file has the unique sample identifiers. Options are "IID" (default) and "FID".
#' @param quiet Logical: should console messages be silenced? Defaults to FALSE
#' @param outfile The string with the name of the `.log` file
#'
#' @return a list with 7 components:
#' * `na_counts`: vector of missing SNP counts in `genotypes`
#' * `obj`: a modified `bigSNP` list with additional components
#' * `og_plink_ids`: either the IID or FID column from `.fam`, determined by `id_var`
#' * `chr`: p-length containing the chromosomes for each SNP
#' * `X`: the `obj$genotypes` as its own FBM
#' * `pos`: vector of physical positions of the SNPs
#' * `chr_range`: vector containing the minimum and maximum values of `chr`.
#'    Character strings are treated as the **maximum**.
#'
#' @keywords internal
#'
name_and_count_bigsnp <- function(obj, id_var, quiet, outfile) {

  X <- obj$genotypes

  # set object names
  obj$colnames <- obj$map$marker.ID

  if (id_var == "FID") {
    obj$rownames <- as.character(obj$fam$family.ID)
  } else if (id_var == "IID") {
    obj$rownames <- as.character(obj$fam$sample.ID)
  }

  chr_range <- range(obj$map$chromosome)
  # save the dimensions of the *original* (pre-standardized) design matrix
  obj$n <- X$nrow
  obj$p <- X$ncol

  if (!quiet) {
    cat("\nThere are", obj$n, "observations and",
        obj$p, "genomic features in the specified data files, representing chromosomes",
        chr_range[1], "-", chr_range[2], "\n")
  }
  cat("\nThere are", obj$n, "observations and",
      obj$p, "genomic features in the specified data files, representing chromosomes",
      chr_range[1], "-", chr_range[2], "\n",
      file = outfile, append = TRUE)

  # save these counts
  counts <- bigstatsr::big_counts(X) # NB: this is a matrix

  list(na_counts = counts[4, ],
       obj = obj,
       og_plink_ids = obj$rownames,
       chr = obj$map$chromosome,
       X = X,
       pos = obj$map$physical.pos,
       chr_range = chr_range)
}

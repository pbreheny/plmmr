#' A function to align genotype and phenotype data
#'
#' @param geno      An object created by `process_plink()`
#' @param rds_dir   The path to the directory in which you want to create the new '.rds' and '.bk' files.
#' @param geno_id   A character string indicating the ID column name in the 'fam'
#'                  element of the genotype data list. Defaults to 'sample.ID', equivalent to 'IID' in PLINK. The other option is 'family.ID', equivalent to 'FID' in PLINK.
#' @param pheno     A data frame with at least two columns: and ID column and a phenotype column
#' @param pheno_id  A string specifying the name of the ID column in `pheno`
#' @param pheno_col A string specifying the name of the phenotype column in `pheno`. This column will be used as the default `y` argument to 'plmm()'.
#' @param outfile   A string with the name of the filepath for the log file
#' @keywords        internal
#'
#' @returns       The indices of
#' @param quiet   Logical: should messages be printed to the console? Defaults to FALSE (which leaves the print messages on...)
add_external_phenotype <- function(geno, rds_dir, geno_id = "sample.ID",
                                      pheno, pheno_id, pheno_col, outfile, quiet){

  # check to make sure IDs overlap
  if (inherits(pheno, 'matrix')) {
    overlap <- intersect(geno$fam[[geno_id]], pheno[,pheno_id])
  } else if (inherits(pheno, 'numeric')) {
    overlap <- intersect(geno$fam[[geno_id]], pheno[[pheno_id]])
  }

  if (length(overlap) < 10) {
    stop("\nThe amount of overlap between the supplied IDs is less than 10 observations.
         This seems really small -- are you sure you chose the right variable names?
         \nAre the ID columns in the genotype and phenotype data of the same type (e.g., both characters)?")
  }

  id_to_keep <- which(geno$fam[[geno_id]] %in% overlap) # indices of IDs in the order of the PLINK data

  # case 1: subset geno data to keep only IDs with corresponding rows in pheno data
  if (length(id_to_keep) < nrow(geno$fam)){

    if (!quiet) {
      cat("Based on the 'id_var' argument you supplied, a total of", length(overlap),
          "samples are in both your genotype and phenotype data. We will subset our analysis to include only these samples.\n")
    }

    cat("Based on the 'id_var' argument you supplied, a total of", length(overlap),
        "samples are in both your genotype and phenotype data. We will subset our analysis to include only these samples.\n",
        file = outfile, append = TRUE)

  } else {
    id_to_keep <- 1:nrow(geno$X)
  }

  # subset and order pheno data to match genotype data
  subset_pheno <- data.table::merge.data.table(x = data.table::as.data.table(geno$fam),
                                               by.x = geno_id,
                                               y = data.table::as.data.table(pheno),
                                               by.y = pheno_id)

  # add new phenotype value to fam file
  geno_keep$fam[6] <- subset_pheno[[pheno_col]]

  return(geno_keep)

}


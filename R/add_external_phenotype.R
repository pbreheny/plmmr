#' A function to align genotype and phenotype data
#'
#' @param geno      An object created by `process_plink()`
#' @param geno_id   A character string indicating the ID column name in the 'fam'
#'                  element of the genotype data list. Defaults to 'sample.ID', equivalent to 'IID' in PLINK. The other option is 'family.ID', equivalent to 'FID' in PLINK.
#' @param pheno     A data frame with at least two columns: and ID column and a phenotype column
#' @param pheno_id  A string specifying the name of the ID column in `pheno`
#' @param pheno_col A string specifying the name of the phenotype column in `pheno`. This column will be used as the default `y` argument to 'plmm()'.
#' @param outfile   A string with the name of the filepath for the log file
#' @keywords        internal
#'
#' @returns         A `bigSNP` object that has had the phenotype values from the external file assigned as the values in the 6th column of the 'fam' object.
#'                  Results are written to an RDS file, then reattached with `bigsnpr;:snp_attach()` and returned
#'                  By default, the RDS file in `geno` is overwritten. You may change this behavior by specifying a new
#'                  filename to `rdsfile`
#' @param quiet   Logical: should messages be printed to the console? Defaults to FALSE (which leaves the print messages on...)
add_external_phenotype <- function(geno, geno_id = "sample.ID",
                                      pheno, pheno_id, pheno_col, outfile, quiet){
browser()
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

  id_to_keep <- which(geno$fam[[geno_id]] %in% overlap) # indices of IDs

  # case 1: subset geno data to keep only IDs with corresponding rows in pheno data
  if (length(id_to_keep) < nrow(geno$fam)){

    if (!quiet) {
      cat("Based on the 'id_var' argument you supplied, a total of", length(overlap),
          "samples are in both your genotype and phenotype data. We will subset our analysis to include only these samples.\n")
    }

    cat("Based on the 'id_var' argument you supplied, a total of", length(overlap),
        "samples are in both your genotype and phenotype data. We will subset our analysis to include only these samples.\n",
        file = outfile, append = TRUE)

    geno_keep_file <- bigsnpr::snp_subset(geno,
                                          ind.row = id_to_keep)
    # TODO: should the subset created in the lines above be given a .bk file in a
    #   temporary directory instead?
    geno_keep_file <- bigmemory::big.matrix(nrow = length(id_to_keep),
                                            ncol = ncol(obj$X),
                                            backingfile = )

    geno_keep <- bigsnpr::snp_attach(geno_keep_file)
    cat("\nMerging the genotype data and phenotype information; new RDS is ", geno_keep$genotypes$rds)
    geno_keep$fam <- geno$fam[id_to_keep,]
  } else {
    # case 2: geno and pheno data represent the same samples
    geno_keep <- geno
  }

  # subset and order pheno data to match genotype data
  subset_pheno <- data.table::merge.data.table(x = data.table::as.data.table(geno_keep$fam),
                                               by.x = geno_id,
                                               y = data.table::as.data.table(pheno),
                                               by.y = pheno_id)

  # add new phenotype value to fam file
  geno_keep$fam[6] <- subset_pheno[[pheno_col]]

  return(geno_keep)

}


#' A function to align genotype and phenotype data 
#'
#' @param geno      An object created by `name_and_count_bigsnp()`
#' @param geno_id   A character string indicating the ID column name in the 'fam' 
#'                  element of the genotype data list. Defautls to 'sample.ID', equivalent to 'IID' in PLINK.
#' @param pheno     A data frame with at least two columns: and ID column and a phenotype column 
#' @param pheno_id  A string specifying the name of the ID column in `pheno`
#' @param pheno_col A string specifying the name of the phenotype column in `pheno`. This column will be used as the default `y` argument to 'plmm()'.
#' 
#' @keywords        internal
#' 
#' @returns         Nothing is returned; instead, the results are written to an RDS file. 
#'                  By default, the RDS file in `geno` is overwritten. You may change this behavior by specifying a new 
#'                  filename to `rdsfile` 
#'                  
add_external_phenotype <- function(geno, geno_id = "sample.ID",
                                      pheno, pheno_id, pheno_col){
  
  # check to make sure IDs overlap
  overlap <- intersect(geno$fam[[geno_id]], pheno[[pheno_id]])
  if (length(overlap) < 10) {
    stop("\nThe amount of overlap between the supplied IDs is less than 10 observations.
         This seems really small -- are you sure you chose the right variable names?
         \nAre the ID columns in the genotype and phenotype data of the same type (e.g., both characters)?")
  }
  
  id_to_keep <- which(geno$fam[[geno_id]] %in% overlap) # indices of IDs
  # case 1: subset geno data to keep only IDs with corresponding rows in pheno data
  if (length(id_to_keep) < nrow(geno$fam)){
    geno_keep_file <- bigsnpr::snp_subset(geno, ind.row = id_to_keep)
    geno_keep <- bigsnpr::snp_attach(geno_keep_file)
    cat("\nMerging the genotype data and phenotype information; new RDS is ", geno_keep$genotypes$rds)
    X <- geno_keep$genotypes
    geno_keep$fam <- geno$fam[id_to_keep,]
  } else {
    # case 2: geno and pheno data represent the same samples
    geno_keep <- geno
  }

  # subset and order pheno data to match genotype data
  subset_pheno <- dplyr::inner_join(geno_keep$fam, pheno, by = c("sample.ID" = "IID"))
  
  # add new phenotype value to fam file 
  geno_keep$fam$affection <- subset_pheno[[pheno_col]]

  return(geno_keep)
  
}

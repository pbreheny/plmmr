
#' Preprocess PLINK files
#'
#' This function allows you to preprocess PLINK bed/bim/fam files for use with penalizedLMM functions. Unreliable SNPs are removed and missing values are imptued using either SnpStats, or if not tagged, the HWE mean value.
#' @param group Character argument that is the prefix of your bed/bim/fam files.
#' @keywords
#' @export
#' @examples

preprocess <- function(group){

  require(snpStats)
  require(magrittr)

  cat("\nPreprocessing", group, "data:\n")

  obj <- read.plink(file.path(dataDir, group), na.strings = "-9") # `na.strings = "-9"` makes chr and position show up as 0 instead of NA in the map obj

  # Only keep polygenic SNPs and those on chr1:22
  genotypes <- obj$genotypes[, col.summary(obj$genotypes)$P.AA != 1 &
                               col.summary(obj$genotypes)$P.AB != 1 &
                               col.summary(obj$genotypes)$P.BB != 1 & obj$map$chromosome %in% 1:22]

  cat("\nRemoving ", dim(obj$geno)[2] - dim(genotypes)[2], "SNPs that are not polygenic or outside of chr 1:22\n")

  # Now how many SNPs have missing data?
  missing <- table(col.summary(genotypes)$Call.rate == 1)
  cat("\n", ifelse(length(missing) == 2, missing[1], 0), "SNPs with missing data that we will attempt to impute\n")

  # Impute missing values
  rules <- snp.imputation(genotypes, minA=0)
  out <- impute.snps(rules, genotypes, as.numeric=FALSE)

  # How many SNPs have missing data after imputation?
  missing <- table(col.summary(out)$Call.rate == 1)
  cat("\n", ifelse(length(missing) == 2, missing[1], 0), "SNPs with missing data after imputation\n")

  # How many SNPs cannot be imputed and still have >= 50% missing values?
  missing <- table(col.summary(out)$Call.rate <= 0.5)
  cat("\n", ifelse(length(missing) == 2, missing[2], 0), "of these SNPs have call rates <= 50% and will be removed\n")

  # Throw out SNPs that have >= 50% missingness, even after imputation
  out2 <- out[, col.summary(out)$Call.rate > 0.5]

  # How many SNPs will be imputed with mean values?
  missing <- table(col.summary(out2)$Call.rate == 1)
  cat("\n", ifelse(length(missing) == 2, missing[1], 0), "SNPs still have missing data that will be imputed with the HWE expected value\n")

  # keep snp order for later
  order_snps <- colnames(out2)

  # which SNPs have missingness
  to_impute <- which(col.summary(out2)$Call.rate < 1)
  miss <- out2[, to_impute]

  imputed_mean <- apply(as(miss, "numeric"), 2, function(s){
    these <- which(is.na(s))
    s[these] <- mean(s, na.rm = TRUE)
    s <- mean2g(s)
    return(s)
  })

  out3 <- out2
  out3@.Data[, to_impute] <- imputed_mean

  # make sure there are no missing values remaining
  missing <- 0
  chunks <- ceiling(nrow(out3) / 100)
  start <- 1
  for (i in 1:chunks){
    stop <- min(i*100, nrow(out3))
    missing <- missing + sum(is.na(out3[start:stop]))
    start <- stop + 1
  }

  # keep corresponding map data
  map <- obj$map[colnames(out3),]

  # keep corresponding fam data
  fam <- obj$fam[rownames(out3),]

  # done...
  cat("\nPreprocessing", group, "data DONE!\n\n",
      nrow(out3), "subjects,", ncol(out3), "SNPs,", missing, "missing values\n")

  return(list(genotypes = out3,
              map = map,
              fam = fam))
}


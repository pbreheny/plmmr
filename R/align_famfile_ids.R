#' A helper function to support `process_plink()`
#' This function is called by another helper, `add_predictors_to_bigsnp()`
#' @param id_var String specifying which variable in the PLINK file was the unique id: 'FID' or 'IID'
#' @param quiet Logical: should a message be printed?
#' @param add_predictor External data to include in design matrix. This is the add_predictors... arg in `process_plink()`
#' @param og_plink_ids Character vector with the PLINK ids (FID or IID) from the *original* data (i.e., the data before any subsetting from handling missing phenotypes)
#'
#' @return An object of the same type as add_predictor
#' @keywords internal

align_famfile_ids <- function(id_var, quiet, add_predictor, og_plink_ids){
  if (id_var == "FID") {
    if (!quiet){cat("\nAligning external data with the PLINK .fam data by FID\n")}
  } else {
    if (!quiet){cat("\nAligning external data with the PLINK .fam data by IID\n")}
  }

  # check alignment between the geno/pheno data we've already merged in step 1
  #   and the supplied predictor file. Need to subset predictors to match pheno file
  if (length(og_plink_ids) < nrow(add_predictor)) {
    add_predictor <- add_predictor[rownames(add_predictor) %in% og_plink_ids,]
  } else if (length(og_plink_ids) > nrow(add_predictor)) {
    stop("There are more rows (samples) in the supplied PLINK data than there are rows in the 'add_predictor_ext' data.
         For now, this is not supported by plmmr. You need to subset your PLINK data to represent
         only the rows represented in your 'add_predictor_ext' file.
         There are at least two ways to do this: use the PLINK software directly, or use methods from the R package 'bigsnpr'.
         If you don't have a lot of background in computing, I'd recommend the 'bigsnpr' approach.")
  }

  ordered_ids <- order(rownames(add_predictor), og_plink_ids)
  if (is.vector(add_predictor)) {
    add_predictor <- add_predictor[ordered_ids]
  } else {
    add_predictor <- add_predictor[ordered_ids,]
  }
  # rownames(add_predictor) <- [ordered_ids]
  return(add_predictor)
}


#' A helper function to support `process_delim()`
#' @param rownames_X A character vector with the rownames or IDs of the observations in the filebacked matrix of data
#' @param id_var character vector with the same values as `rownames_X`
#' @param quiet Logical: should a message be printed?
#' @param add_predictor External data to include in design matrix. This is the add_predictors... arg in `process_delim()`
#'
#' @return An object of the same type as add_predictor
#' @keywords internal
#'
align_ids <- function(rownames_X, id_var, quiet, add_predictor){
  # check alignment between the data we've already merged in step 1
  #   and the supplied predictor file. Need to subset predictors to match main dataframe
  if (length(rownames_X) != nrow(add_predictor)) {
    add_predictor <- add_predictor[rownames(add_predictor) %in% rownames_X,]
  }
  ordered_ids <- order(id_var, rownames_X)
  if (is.vector(add_predictor)) {
    add_predictor <- add_predictor[ordered_ids]
  } else {
    add_predictor <- add_predictor[ordered_ids,]
  }
  return(add_predictor)
}

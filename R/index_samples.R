#' A function to align genotype and phenotype data
#'
#' @param obj     An object created by `process_plink()`
#' @param rds_dir   The path to the directory in which you want to create the new '.rds' and '.bk' files.
#' @param indiv_id   A character string indicating the ID column name in the 'fam'
#'                  element of the genotype data list. Defaults to 'sample.ID', equivalent to 'IID' in PLINK. The other option is 'family.ID', equivalent to 'FID' in PLINK.
#' @param add_outcome    A data frame with at least two columns: and ID column and a phenotype column
#' @param outcome_id  A string specifying the name of the ID column in `pheno`
#' @param outcome_col A string specifying the name of the phenotype column in `pheno`. This column will be used as the default `y` argument to 'plmm()'.
#' @param na_outcome_vals A vector of numeric values used to code NA values in the outcome. Defaults to `c(-9, NA_integer)` (the -9 matches PLINK conventions).
#' @param outfile   A string with the name of the filepath for the log file
#' @param quiet     Logical: should messages be printed to the console? Defaults to FALSE (which leaves the print messages on...

#' @keywords        internal
#'
#' @returns a list with two items:
#' * a data.table with rows corresponding to the samples for which both genotype and phenotype are available.
#' * a numeric vector with indices indicating which samples were 'complete' (i.e., which samples from add_outcome had corresponding data in the PLINK files)
index_samples <- function(obj,
                          rds_dir,
                          indiv_id,
                          add_outcome,
                          outcome_id,
                          outcome_col,
                          na_outcome_vals,
                          outfile,
                          quiet) {

  # first, determine which IDs have both feature data and outcome data ---------
  if (inherits(add_outcome, "matrix") || inherits(add_outcome, "data.frame")) {
    overlap <- intersect(indiv_id, add_outcome[, outcome_id])
  } else if (inherits(add_outcome, "numeric")) {
    overlap <- intersect(indiv_id, add_outcome[[outcome_id]])
  }

  # check to make sure IDs overlap
  if (length(overlap) < 10) {
    stop("\nThe amount of overlap between the supplied IDs is less than 10 observations.
         This seems really small -- are you sure you chose the right variable names?")
  }

  id_in_both <- which(indiv_id %in% overlap)

  if (length(id_in_both) < nrow(obj$X)) {

    if (!quiet) {
      cat("Based on the 'id_var' argument you supplied, a total of", length(overlap),
          "samples are in both your processed data and your outcome data. We will subset our analysis to include only these samples.\n")
    }

    cat("Based on the 'id_var' argument you supplied, a total of", length(overlap),
        "samples are in both your processed data and your outcome data. We will subset our analysis to include only these samples.\n",
        file = outfile, append = TRUE)

  }

  # finally, subset & sort add_outcome -----------------------------------------
  indiv_id_df <- data.table::as.data.table(as.character(indiv_id))
  colnames(indiv_id_df) <- "ID"
  complete_samples <- data.table::merge.data.table(x = indiv_id_df,
                                                   by.x = "ID",
                                                   y = data.table::as.data.table(add_outcome),
                                                   by.y = outcome_id,
                                                   sort = FALSE)
  # NB: here, 'complete samples' means samples that appear in both the outcome data
  # *and* the feature data

  # keep the indices
  return(list(complete_samples = complete_samples,
              outcome_idx = which(indiv_id %in% add_outcome[, outcome_id])))

}

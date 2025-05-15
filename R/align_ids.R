#' A helper function to support `process_plink()`
#' @param id_var String specifying the variable name of the ID column
#' @param quiet Logical: should a message be printed?
#' @param add_predictor External data to include in design matrix. This is the add_predictors... arg in `process_plink()`
#' @param og_ids Character vector with the PLINK ids (FID or IID) from the *original* data (i.e., the data before any subsetting from handling missing phenotypes)
#'
#' @return A matrix with the same dimensions as add_predictor
#' @keywords internal
align_ids <- function(id_var, quiet, add_predictor, og_ids) {

  if (is.numeric(add_predictor[, id_var])) {
    add_predictor[, id_var] <- as.character(add_predictor[, id_var])
  }

  # check for alignment
  if (length(intersect(og_ids, add_predictor[, id_var])) == 0) {
    stop("You supplied an argument to 'add_predictor', but the IDs do not align with either of the ID columns in the feature data file.
         \nPlease create or align the names of this vector - alignment is essential for accurate analysis.\n",
         call. = FALSE)
  }

  # check alignment between the geno/pheno data we've already merged in step 1
  #   and the supplied predictor file.
  if (length(og_ids) > nrow(add_predictor)) {
    stop("There are more rows (samples) in the supplied PLINK data than there are rows in the 'add_predictor' data.
         For now, this is not supported by plmmr. You need to subset your PLINK data to represent
         only the rows represented in your 'add_predictor' file.
         There are at least two ways to do this: use the PLINK software directly, or use methods from the R package 'bigsnpr'.
         If you don't have a lot of background in computing, I'd recommend the 'bigsnpr' approach.",
         call. = FALSE)
  }

  if (!quiet) {
    cat("\nAligning external data with the feature data by", id_var, "\n")
  }
  # make original IDs into a table with which the added predictors can be merged
  og_ids <- data.table::as.data.table(data.frame(ID = og_ids))
  add_predictor <- data.table::data.table(add_predictor)
  new_add_predictor <- data.table::merge.data.table(x = og_ids,
                                                    y = add_predictor,
                                                    by.x = "ID",
                                                    by.y = id_var,
                                                    all = FALSE,
                                                    sort = FALSE)

  # for downstream calls, output *must* be a numeric matrix
  new_add_predictor_mat <- as.matrix(new_add_predictor[, -1])
  rownames(new_add_predictor_mat) <- as.matrix(new_add_predictor[, "ID", with = FALSE]) |>
    as.character()

  return(new_add_predictor_mat)
}

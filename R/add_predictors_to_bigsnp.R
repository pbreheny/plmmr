#' A helper function to add predictors to a filebacked matrix of data
#'
#' @param obj               A `bigSNP` object
#' @param add_predictor_fam Optional: if you want to include "sex" (the 5th column of `.fam` file) in the analysis, specify 'sex' here.
#' @param add_predictor_ext Optional: add additional covariates/predictors/features from an external file (i.e., not a PLINK file).
#' @param id_var            String specifying which column of the PLINK `.fam` file has the unique sample identifiers. Options are "IID" (default) and "FID".
#' @param og_plink_ids      Character vector passed from `name_and_count_bigsnp()`
#' @param rds_dir           The path to the directory in which you want to create the new '.rds' and '.bk' files. Defaults to `data_dir`(from `process_plink()` call)
#' @param quiet             Logical: should messages be printed to the console? Defaults to FALSE (which leaves the print messages on...)
#'
#' @return A list of 2 components:
#' * 'obj' - a `bigSNP` object with an added element representing the matrix that includes the additional predictors as the first few columns
#' * 'non_gen' - an integer vector that ranges from 1 to the number of added predictors. Example: if 2 predictors are added, non_gen = 1:2
#' @keywords internal
#'
add_predictors <- function(obj,
                           add_predictor_ext,
                           id_var,
                           og_ids,
                           rds_dir,
                           rds_file,
                           quiet){


  # add additional covariates -----------------------
  # first, set up some indices; even if no additional args are used, these NULL
  #   values are important for checks downstream
  non_gen <- NULL

  if (!is.null(add_predictor_ext)) {
    if (!quiet) {
      cat("Adding predictors from external data.\n")
    }
    # vector case ----------------------------------------------------------------
    if (is.vector(add_predictor_ext)) {

      # make sure types match
      if (!is.numeric(add_predictor_ext)) {
        stop("The vector supplied to the 'add_predictor_ext' argument must be numeric.\n")
      }
      names(add_predictor_ext) <- as.numeric(names(add_predictor_ext))

      if (var(add_predictor_ext) == 0) {
        stop("The supplied argument to add_predictor_ext is constant (no variation).
             This would not be a meaningful predictor.\n")
      }

      # check for alignment
      if (is.null(names(add_predictor_ext)) |
          length(intersect(og_plink_ids, names(add_predictor_ext))) == 0) {
        stop("You supplied an argument to 'add_predictor_ext', but the names of this
         vector either (a) do not exist or (b) do not align with either of the ID columns in the PLINK fam file.
         \nPlease create or align the names of this vector - alignment is essential for accurate analysis.\n")
      }

      if (!quiet) cat('Aligning IDs between fam and predictor files\n')
      add_predictor_ext <- align_famfile_ids(id_var = id_var,
                                             quiet = quiet,
                                             add_predictor = add_predictor_ext,
                                             og_plink_ids = og_plink_ids)

      # save non_gen: an index marking the first column as non-genomic predictor
      non_gen <- 1

      obj$geno_plus_predictors <- bigstatsr::FBM(type = 'double',
                                                 nrow = nrow(obj$fam),
                                                 ncol = obj$genotypes$ncol + length(non_gen))
      # fill in new matrix
      obj$geno_plus_predictors <- big_cbind(A = as.matrix(add_predictor_ext),
                                            B = geno_bm,
                                            C = obj$geno_plus_predictors,
                                            quiet = quiet)

      # adjust colnames
      obj$colnames <- c(deparse(substitute(add_predictor_ext)), obj$colnames)

    } else if (is.matrix(add_predictor_ext) | is.data.frame(add_predictor_ext)) {
      # matrix case --------------------------------------------------------
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
          length(intersect(og_ids, rownames(add_predictor_ext))) == 0) {
        stop("\nYou supplied an argument to 'add_predictor_ext', but the row names of this
         matrix either (a) do not exist or (b) do not align with either of the ID columns in the PLINK fam file.
         \nPlease create or align the names of this matrix - alignment is essential for accurate analysis.")
      }

      if (!quiet) cat('Aligning IDs between fam and predictor files\n')
      add_predictor_ext <- align_famfile_ids(id_var = id_var,
                                             quiet = quiet,
                                             add_predictor = add_predictor_ext,
                                             og_ids = og_ids)

      # save non_gen: an index marking added columns as non-genomic predictors
      non_gen <- 1:ncol(add_predictor_ext)


      design_matrix <- big.matrix(nrow = nrow(obj$fam), # TODO: think about whether this should be nrow(obj$X)...
                                  ncol = ncol(obj$X) + length(non_gen),
                                  type = 'double',
                                  backingfile = "unstd_design_matrix.bk",
                                  backingpath = rds_dir,
                                  descriptorfile = "unstd_design_matrix.desc")

      design_matrix <- big_cbind(A = add_predictor_ext,
                                 B = obj$X,
                                 C = design_matrix,
                                 quiet = quiet)

      ret <- list(design_matrix = design_matrix, non_gen = non_gen)
      # adjust colnames if applicable
      if (!is.null(colnames(add_predictor_ext))){
        ret$colnames <- c(colnames(add_predictor_ext), obj$colnames)
      }

    }


    return(ret)

  }
}


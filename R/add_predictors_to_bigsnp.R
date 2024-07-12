#' A helper function to add predictors to a filebacked matrix of data
#'
#' @param obj               A `bigSNP` object
#' @param add_predictor_fam Optional: if you want to include "sex" (the 5th column of `.fam` file) in the analysis, specify 'sex' here.
#' @param add_predictor_ext Optional: add additional covariates/predictors/features from an external file (i.e., not a PLINK file).
#' @param id_var            String specifying which column of the PLINK `.fam` file has the unique sample identifiers. Options are "IID" (default) and "FID".
#' @param og_plink_ids      Character vector passed from `name_and_count_bigsnp()`
#' @param rds_dir           The path to the directory in which you want to create the new '.rds' and '.bk' files. Defaults to `data_dir`(from `process_plink()` call)
#' @param quiet             Logical: should messages be printed to the console? Defaults to TRUE
#'
#' @return A list of 2 components:
#' * 'obj' - a `bigSNP` object with an added element representing the matrix that includes the additional predictors as the first few columns
#' * 'non_gen' - an integer vector that ranges from 1 to the number of added predictors. Example: if 2 predictors are added, non_gen = 1:2
#' @keywords internal
#'
add_predictors_to_bigsnp <- function(obj, add_predictor_fam, add_predictor_ext,
                                     id_var, og_plink_ids,
                                     rds_dir, quiet){

  # genotypes need to have type 'double' from now on, in order to merge
  geno_bm <- fbm2bm(bigstatsr::big_copy(X = obj$genotypes, type = "double"))

  # add additional covariates -----------------------
  # first, set up some indices; even if no additional args are used, these NULL
  #   values are important for checks downstream
  non_gen <- NULL
  ## covariates from .fam file ---------------------------
  if (!is.null(add_predictor_fam)) {
    if (!quiet) {
      cat("\nAdding predictors from .fam file.")
    }
    if (add_predictor_fam == "sex"){

      # add space for extra column
      obj$geno_plus_predictors <- bigstatsr::FBM(init = 0,
                                                 nrow = nrow(obj$fam),
                                                 ncol = obj$genotypes$ncol + 1)
      # fill in new matrix
      obj$geno_plus_predictors <- big_cbind(A = as.matrix(obj$fam[add_predictor_fam]),
                                            B = geno_bm,
                                            C = obj$geno_plus_predictors,
                                            quiet = quiet)

      # adjust colnames
      obj$colnames <- c(add_predictor_fam, obj$colnames)

      # save non_gen: an index marking the first column as non-genomic predictor
      non_gen <- 1

    }
  }

  ## covariates from external file   ----------------------------------
  if (!is.null(add_predictor_ext)) {
    if (!quiet) {
      cat("\nAdding predictors from external data.")
    }
    if (is.vector(add_predictor_ext)) {
      ### vector case -------------------------------
      # make sure types match
      if (!is.numeric(add_predictor_ext)) {
        stop("\nThe vector supplied to the 'add_predictor_ext' argument must be numeric.")
      }
      names(add_predictor_ext) <- as.numeric(names(add_predictor_ext))

      if (var(add_predictor_ext) == 0) {
        stop("\nThe supplied argument to add_predictor_ext is constant (no variation).
             This would not be a meaningful predictor.")
      }

      # check for alignment
      if (is.null(names(add_predictor_ext)) |
          length(intersect(og_plink_ids, names(add_predictor_ext))) == 0) {
        stop("\nYou supplied an argument to 'add_predictor_ext', but the names of this
         vector either (a) do not exist or (b) do not align with either of the ID columns in the PLINK fam file.
         \nPlease create or align the names of this vector - alignment is essential for accurate analysis.")
      }

      add_predictor_ext <- align_famfile_ids(id_var = id_var,
                                             quiet = quiet,
                                             add_predictor = add_predictor_ext,
                                             og_plink_ids = og_plink_ids)

      # save non_gen: an index marking the first column as non-genomic predictor
      non_gen <- 1

      obj$geno_plus_predictors <- bigstatsr::FBM(init = 0,
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
      ### matrix case --------------------------------
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
          length(intersect(og_plink_ids, rownames(add_predictor_ext))) == 0) {
        stop("\nYou supplied an argument to 'add_predictor_ext', but the row names of this
         matrix either (a) do not exist or (b) do not align with either of the ID columns in the PLINK fam file.
         \nPlease create or align the names of this matrix - alignment is essential for accurate analysis.")
      }

      add_predictor_ext <- align_famfile_ids(id_var = id_var,
                                             quiet = quiet,
                                             add_predictor = add_predictor_ext,
                                             og_plink_ids = og_plink_ids)

      # save non_gen: an index marking added columns as non-genomic predictors
      non_gen <- 1:ncol(add_predictor_ext)

      obj$geno_plus_predictors <- bigstatsr::FBM(init = 0,
                                                 type = "double",
                                                 #backingfile = file.path(rds_dir, "combined_data"
                                                 nrow = nrow(obj$fam),
                                                 ncol = obj$genotypes$ncol + length(non_gen)) |> fbm2bm()

      obj$geno_plus_predictors <- big_cbind(A = add_predictor_ext,
                                          B = geno_bm,
                                          C = obj$geno_plus_predictors,
                                          quiet = quiet)


      # adjust colnames if applicable
      if (!is.null(colnames(add_predictor_ext))){
        obj$colnames <- c(colnames(add_predictor_ext), obj$colnames)
      }

    }


  return(list(obj = obj, non_gen = non_gen))

  }
}


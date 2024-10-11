#' A function to create a design with an in-memory X matrix
#'
#' @param X             A numeric matrix in which rows correspond to observations (e.g., samples) and columns correspond to features.
#' @param outcome_col   A numeric vector representing the outcome for the model.
#'                      **Note**: it is the responsibility of the user to ensure that the outcome_col and X have the same row order!
#' @param unpen         An optional character vector with the names of columns to mark as unpenalized (i.e., these features would always be included in a model).
#'                      **Note**: if you choose to use this option, X must have column names.
#'
#' @return A list with elements including a standardized X and model design information
#' @keywords internal
#'
create_design_in_memory <- function(X, outcome_col, unpen = NULL){

  # standardize X
  std_X <- ncvreg::std(X)
  design <- list(std_X = std_X,
                 std_X_n = nrow(std_X),
                 std_X_p = ncol(std_X),
                 ns = attr(std_X, 'nonsingular'),
                 std_X_center = attr(std_X, 'center'),
                 std_X_scale = attr(std_X, 'scale'))

  # format meta-data
    design$X_colnames <- colnames(X)
    design$X_rownames <- rownames(X)
    design$n <- nrow(X)
    design$p <- ncol(X)
    design$y <- outcome_col
    design$std_X_rownames <- rownames(X)
    design$std_X_colnames <- colnames(X)[design$ns]

# handle unpenalized columns
  if (is.null(unpen)) {
    design$penalty_factor <- rep(1, design$std_X_p)
  } else {
    # save indices for unpenalized covariates
    design$unpen_colnames <- unpen
    design$penalty_factor <- rep(1, ncol(design$std_X))
    design$unpen <- which_unpen <- which(unpen %in% design$std_X_colnames)
    design$penalty_factor[which_unpen] <- 0
  }
  return(structure(design, class = "plmm_design"))
}
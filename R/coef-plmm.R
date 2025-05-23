#' Coef method for "plmm" class
#'
#' @param object An object of class "plmm."
#' @param lambda A numeric vector of lambda values.
#' @param which Vector of lambda indices for which coefficients to return.
#' @param drop Logical.
#' @param ... Additional arguments.
#'
#' @rdname coef.plmm
#'
#' @returns Either a numeric matrix (if model was fit on data stored in memory)
#' or a sparse matrix (if model was fit on data stored filebacked). Rownames are
#' feature names, columns are values of `lambda`.
#'
#' @export
#'
#' @examples
#' admix_design <- create_design(X = admix$X, y = admix$y)
#' fit <- plmm(design = admix_design)
#' coef(fit)[1:10, 41:45]


coef.plmm <- function(object, lambda, which = seq_along(object$lambda), drop = TRUE, ...) {
  # error check for supplied lambda value
  if (!missing(lambda)) {
    if (max(lambda) > max(object$lambda) || min(lambda) < min(object$lambda)) {
      stop("Supplied lambda value(s) are outside the range of the model fit.",
           call. = FALSE)
    }

    ind <- stats::approx(object$lambda, seq(object$lambda), lambda)$y
    l <- floor(ind)
    r <- ceiling(ind)
    w <- ind %% 1
    beta_vals <- (1 - w) * object$beta_vals[, l, drop = FALSE] +
      w * object$beta_vals[, r, drop = FALSE]

    # format dim. names
    if (is.null(dim(beta_vals))) {
      # case 1: beta_vals is a vector
      names(beta_vals) <- lam_names(lambda)
    } else {
      # case 2: beta_vals is a matrix
      colnames(beta_vals) <- lam_names(lambda)
    }

  } else {
    beta_vals <- object$beta_vals[, which, drop = FALSE]
  }

  if (drop) {
    return(drop(beta_vals))
  } else {
    return(beta_vals)
  }
}

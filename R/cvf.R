#' Cross-validation internal function for cv_plmm
#'
#' Internal function for cv_plmm which calls plmm on a fold subset of the original data.
#'
#' @param i Fold number to be excluded from fit.
#' @param fold n-length vector of fold-assignments.
#' @param type A character argument indicating what should be returned from predict.plmm. If \code{type == 'lp'} predictions are based on the linear predictor, \code{$X beta$}. If \code{type == 'individual'} predictions are based on the linear predictor plus the estimated random effect (BLUP).
#' @param cv_args List of additional arguments to be passed to plmm.
#' @param estimated_V Estimated variance-covariance matrix using all observations when computing BLUP; NULL if type = "lp" in cv_plmm.
#' @param ... Optional arguments to `predict_within_cv`
#'
#' @keywords internal

cvf <- function(i, fold, type, cv_args, estimated_V, ...) {

  # save the 'prep' object from the plmm_prep() in cv_plmm
  full_cv_prep <- cv_args$prep
  y <- cv_args$y
  # subset std_X, U, and y to match fold indices
  #   (and in so doing, leave out the ith fold)
  if (cv_args$fbm_flag) {
    cv_args$prep$std_X <- bigstatsr::big_copy(full_cv_prep$std_X,
                                              ind.row = which(fold!=i)) |> fbm2bm()
    # re-scale data & check for singularity
    train_data <- .Call("big_std",
                           cv_args$prep$std_X@address,
                           as.integer(bigstatsr::nb_cores()),
                           PACKAGE = "plmmr")
    cv_args$prep$std_X@address <- train_data$std_X
    cv_args$std_X_details$center <- train_data$std_X_center
    cv_args$std_X_details$scale <- train_data$std_X_scale
    cv_args$std_X_details$ns <- which(train_data$std_X_scale > 1e-3)
    singular <- train_data$std_X_scale < 1e-3

    # do not fit a model on singular features!
    if (sum(singular) >= 1) cv_args$penalty.factor[singular] <- Inf

  } else {
    cv_args$prep$std_X <- full_cv_prep$std_X[fold!=i, ,drop=FALSE]
    # Note: subsetting the data into test/train sets may cause low variance features
    #   to become constant features in the training data. The following lines address this issue

    # re-scale data & check for singularity
    cv_args$prep$std_X <- ncvreg::std(cv_args$prep$std_X)
    cv_args$std_X_details$center <- attr(cv_args$prep$std_X, "center")
    cv_args$std_X_details$scale <- attr(cv_args$prep$std_X, "scale")
    cv_args$std_X_details$ns <- attr(cv_args$prep$std_X, "nonsingular")

    # do not fit a model on these singular features!
    cv_args$penalty.factor <- cv_args$penalty.factor[cv_args$std_X_details$ns]

  }

  cv_args$prep$U <- full_cv_prep$U[fold!=i, , drop=FALSE]
  cv_args$prep$centered_y <- full_cv_prep$centered_y[fold!=i] |> scale(scale=FALSE)
  cv_args$y <- y[fold!=i]

  # extract test set (comes from cv prep on full data)
  if (cv_args$fbm_flag){
    test_X <- bigstatsr::big_copy(full_cv_prep$std_X,
                                  ind.row = which(fold==i)) |> fbm2bm()

  } else {
    test_X <- full_cv_prep$std_X[fold==i, cv_args$std_X_details$ns, drop=FALSE]
  }
  test_y <- y[fold==i]

  # NB: we are assuming that the eta is the same across the training and testing data.

  # fit a plmm within each fold at each value of lambda
  # lambda stays the same for each fold; comes from the overall fit in cv_plmm()
  if (cv_args$prep$trace) {
    cat("Fitting model in fold ", i, ":\n")
  }

  fit.i <- do.call("plmm_fit", cv_args)

  if(type == "lp"){
    yhat <- predict_within_cv(fit = fit.i,
                              oldX = cv_args$prep$std_X,
                              newX = test_X,
                              type = 'lp',
                              fbm = cv_args$fbm_flag)
  }

  if (type == 'blup'){
    # estimated_V here comes from the overall fit in cv_plmm.R, an n*n matrix
    V21 <- estimated_V[fold==i, fold!=i, drop = FALSE]
    V11 <- estimated_V[fold!=i, fold!=i, drop = FALSE]

    yhat <- predict_within_cv(fit = fit.i,
                              oldX = cv_args$prep$std_X,
                              newX = test_X,
                              type = 'blup',
                              fbm = cv_args$fbm_flag,
                              V11 = V11,
                              V21 = V21, ...)

  }

  loss <- sapply(1:ncol(yhat), function(ll) plmm_loss(test_y, yhat[,ll]))
  list(loss=loss, nl=length(fit.i$lambda), yhat=yhat)
}

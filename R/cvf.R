#' Cross-validation internal function for cv_plmm
#'
#' Internal function for cv_plmm which calls plmm on a fold subset of the original data.
#'
#' @param i Fold number to be excluded from fit.
#' @param fold n-length vector of fold-assignments.
#' @param type A character argument indicating what should be returned from predict.plmm. If \code{type == 'lp'} predictions are based on the linear predictor, \code{$X beta$}. If \code{type == 'individual'} predictions are based on the linear predictor plus the estimated random effect (BLUP).
#' @param cv.args List of additional arguments to be passed to plmm.
#' @param estimated_V Estimated variance-covariance matrix using all observations when computing BLUP; NULL if type = "lp" in cv_plmm.
#' @param ... Optional arguments to `predict_within_cv`
#'
#' @keywords internal

cvf <- function(i, fold, type, cv.args, estimated_V, ...) {

  # save the 'prep' object from the plmm_prep() in cv_plmm
  full_cv_prep <- cv.args$prep
  # subset std_X, U, and y to match fold indices
  #   (and in so doing, leave out the ith fold)
  if (cv.args$fbm_flag) {
    cv.args$prep$std_X <- bigstatsr::big_copy(full_cv_prep$std_X,
                                              ind.row = which(fold!=i)) |> fbm2bm()

    # check for singularity
    train_data_sd <- .Call("big_sd",
                           cv.args$prep$std_X@address,
                           as.integer(bigstatsr::nb_cores()),
                           PACKAGE = "plmmr")
    singular <- train_data_sd$sd_vals < 1e-3

    # do not fit a model on these singular features!
    if (sum(singular) >= 1) cv.args$penalty.factor[singular] <- Inf

  } else {
    cv.args$prep$std_X <- full_cv_prep$std_X[fold!=i, ,drop=FALSE]
    # Note: subsetting the data into test/train sets may cause low variance features
    #   to become constant features in the training data. The following lines address this issue

    # check 1
    singular <- apply(cv.args$prep$std_X, 2, sd) < 1e-3

    # check 2
    # balance_ratios <- apply(cv.args$prep$std_X, 2, balance_ratio)
    # singular <- balance_ratios > 19 # per Krstajic et al. 2014

    # check 3
    # bm_std_X <- bigstatsr::as_FBM(cv.args$prep$std_X) |> fbm2bm()
    # train_data_sd <- .Call("big_sd",
    #                        bm_std_X@address,
    #                        as.integer(bigstatsr::nb_cores()),
    #                        PACKAGE = "plmmr")
    # singular <- train_data_sd$sd_vals < 1e-3
    #
    # # do not fit a model on these singular features!
    if (sum(singular) >= 1) cv.args$penalty.factor[singular] <- Inf
  }

  cv.args$prep$U <- full_cv_prep$U[fold!=i, , drop=FALSE]
  cv.args$prep$y <- full_cv_prep$y[fold!=i]# |> scale(scale=FALSE)

  # extract test set (comes from cv prep on full data)
  if (cv.args$fbm_flag){
    test_X <- bigstatsr::big_copy(full_cv_prep$std_X,
                                  ind.row = which(fold==i)) |> fbm2bm()

  } else {
    test_X <- full_cv_prep$std_X[fold==i, , drop=FALSE]
  }
  test_y <- full_cv_prep$y[fold==i]

  # NB: we are assuming that the eta is the same across the training and testing data.

  # fit a plmm within each fold at each value of lambda
  # lambda stays the same for each fold; comes from the overall fit in cv_plmm()
  if (cv.args$prep$trace) {
    cat("Fitting model in fold ", i, ":\n")
  }

  fit.i <- do.call("fit_within_cv", cv.args)
  # checks
  # stdrot_beta_max_check <- apply(fit.i$stdrot_scale_beta, 1, max)
  # if (any(abs(stdrot_beta_max_check) > 10)) browser()
  #
  # stdrot_beta_min_check <- apply(fit.i$stdrot_scale_beta, 1, min)
  # if (any(abs(stdrot_beta_min_check) > 10)) browser()
  #
  # std_beta_max_check <- apply(fit.i$std_scale_beta, 1, max)
  # if (any(abs(std_beta_max_check) > 10)) browser()
  #
  # std_beta_min_check <- apply(fit.i$std_scale_beta, 1, min)
  # if (any(abs(std_beta_min_check) > 10)) browser()


  if(type == "lp"){
    yhat <- predict_within_cv(fit = fit.i,
                              oldX = cv.args$prep$std_X,
                              newX = test_X,
                              type = 'lp',
                              fbm = cv.args$fbm_flag)
  }

  if (type == 'blup'){
    # estimated_V here comes from the overall fit in cv_plmm.R, an n*n matrix
    V21 <- estimated_V[fold==i, fold!=i, drop = FALSE]
    V11 <- estimated_V[fold!=i, fold!=i, drop = FALSE]

    yhat <- predict_within_cv(fit = fit.i,
                              oldX = cv.args$prep$std_X,
                              newX = test_X,
                              type = 'blup',
                              fbm = cv.args$fbm_flag,
                              V11 = V11,
                              V21 = V21, ...)

  }

  loss <- sapply(1:ncol(yhat), function(ll) plmm_loss(test_y, yhat[,ll]))
  list(loss=loss, nl=length(fit.i$lambda), yhat=yhat)
}

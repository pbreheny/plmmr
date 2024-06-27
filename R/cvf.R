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
  # make list to hold the data for this particular fold:
  fold_args <- list(std_X_details = list(),
                    fbm_flag = cv_args$fbm_flag,
                    penalty = cv_args$penalty,
                    penalty.factor = cv_args$penalty.factor,
                    gamma = cv_args$gamma,
                    alpha = cv_args$alpha,
                    nlambda = cv_args$nlambda,
                    lambda.min = cv_args$lambda.min,
                    max.iter = cv_args$max.iter,
                    eps = cv_args$eps,
                    warn = cv_args$warn,
                    convex = cv_args$convex,
                    dfmax = cv_args$dfmax,
                    lambda = cv_args$lambda)
  # subset std_X, U, and y to match fold indices ------------------------
  #   (and in so doing, leave out the ith fold)
  if (cv_args$fbm_flag) {

    fold_args$std_X <- train_X <- bigstatsr::big_copy(full_cv_prep$std_X,
                                              ind.row = which(fold!=i)) |> fbm2bm()
    # re-scale data & check for singularity
    train_data <- .Call("big_std",
                           fold_args$std_X@address,
                           as.integer(bigstatsr::nb_cores()),
                           PACKAGE = "plmmr")
    fold_args$std_X@address <- train_data$std_X
    fold_args$std_X_details$center <- train_data$std_X_center
    fold_args$std_X_details$scale <- train_data$std_X_scale
    fold_args$std_X_details$ns <- which(train_data$std_X_scale > 1e-3)
    singular <- train_data$std_X_scale < 1e-3

    # do not fit a model on singular features!
    if (sum(singular) >= 1) fold_args$penalty.factor[singular] <- Inf

  } else {
    fold_args$std_X <- train_X <- full_cv_prep$std_X[fold!=i, ,drop=FALSE]
    # Note: subsetting the data into test/train sets may cause low variance features
    #   to become constant features in the training data. The following lines address this issue

    # re-scale data & check for singularity
    fold_args$std_X <- ncvreg::std(fold_args$std_X) # notice: singular columns are *removed* here
    fold_args$std_X_details$center <- attr(fold_args$std_X, "center")
    fold_args$std_X_details$scale <- attr(fold_args$std_X, "scale")
    fold_args$std_X_details$ns <- attr(fold_args$std_X, "nonsingular")

    # do not fit a model on these singular features!
    fold_args$penalty.factor <- fold_args$penalty.factor[fold_args$std_X_details$ns]

  }
  fold_args$centered_y <- full_cv_prep$centered_y[fold!=i] |> scale(scale=FALSE) |> drop()
  fold_args$y <- y[fold!=i]

  # extract test set --------------------------------------
  # this comes from cv prep on full data
  if (cv_args$fbm_flag){
    test_X <- bigstatsr::big_copy(full_cv_prep$std_X,
                                  ind.row = which(fold==i)) |> fbm2bm()

  } else {
    test_X <- full_cv_prep$std_X[fold==i, fold_args$std_X_details$ns, drop=FALSE]
  }
  test_y <- y[fold==i]

  # decomposition for current fold ------------------------------
  if (cv_args$prep$trace) {
    cat("Beginning eigendecomposition in fold ", i, ":\n")
  }
  fold_prep <- plmm_prep(std_X = fold_args$std_X,
                         std_X_n = nrow(fold_args$std_X),
                         std_X_p = ncol(fold_args$std_X),
                         n = nrow(full_cv_prep$std_X),
                         p = ncol(full_cv_prep$std_X),
                         centered_y = fold_args$centered_y,
                         fbm_flag = fold_args$fbm_flag,
                         penalty.factor = fold_args$penalty.factor,
                         trace = cv_args$prep$trace)

  fold_args$prep <- fold_prep

  # fit a plmm within each fold at each value of lambda
  # lambda stays the same for each fold; comes from the overall fit in cv_plmm()
  if (cv_args$prep$trace) {
    cat("Fitting model in fold ", i, ":\n")
  }

  fit.i <- plmm_fit(prep = fold_prep,
                    y = fold_args$y,
                    std_X_details = fold_args$std_X_details,
                    penalty.factor = fold_args$penalty.factor,
                    fbm_flag = fold_args$fbm_flag,
                    penalty = fold_args$penalty,
                    gamma = fold_args$gamma,
                    alpha = fold_args$alpha,
                    lambda.min = fold_args$lambda.min,
                    nlambda = fold_args$nlambda,
                    lambda = fold_args$lambda,
                    eps = fold_args$eps,
                    max.iter = fold_args$max.iter,
                    warn = fold_args$warn,
                    convex = fold_args$convex,
                    dfmax = ncol(train_X) + 1)

  format.i <- plmm_format(fit = fit.i,
              p =  ncol(train_X),
              std_X_details = fold_args$std_X_details,
              use_feature_names = FALSE, # no need for names in internal CV fits
              # TODO: figure out how to track non_genomic features in CV
              non_genomic = NULL,
              # Note: must keep non_genomic length 0 for untransform() to work correctly within each fold
              fbm_flag = fold_args$fbm_flag)

  if(type == "lp"){
    yhat <- predict_within_cv(fit = fit.i,
                              trainX = train_X,
                              testX = test_X,
                              og_scale_beta = format.i$beta_vals,
                              type = 'lp',
                              fbm = cv_args$fbm_flag)
  }

  if (type == 'blup'){
    # estimated_V here comes from the overall fit in cv_plmm.R, an n*n matrix
    V21 <- estimated_V[fold==i, fold!=i, drop = FALSE]
    V11 <- estimated_V[fold!=i, fold!=i, drop = FALSE]

    yhat <- predict_within_cv(fit = fit.i,
                              trainX = train_X,
                              trainY = fold_args$y,
                              testX = test_X,
                              og_scale_beta = format.i$beta_vals,
                              std_X_details = fold_args$std_X_details,
                              type = 'blup',
                              fbm = cv_args$fbm_flag,
                              V11 = V11,
                              V21 = V21, ...)

  }

  loss <- sapply(1:ncol(yhat), function(ll) plmm_loss(test_y, yhat[,ll]))
  list(loss=loss, nl=length(fit.i$lambda), yhat=yhat)
}

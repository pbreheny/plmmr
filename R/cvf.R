#' Cross-validation internal function for cv_plmm
#'
#' Internal function for cv_plmm which calls plmm on a fold subset of the original data.
#'
#' @param i Fold number to be excluded from fit.
#' @param fold n-length vector of fold-assignments.
#' @param type A character argument indicating what should be returned from predict.plmm. If \code{type == 'lp'} predictions are based on the linear predictor, \code{$X beta$}. If \code{type == 'individual'} predictions are based on the linear predictor plus the estimated random effect (BLUP).
#' @param cv_args List of additional arguments to be passed to plmm.
#' @param estimated_Sigma Estimated variance-covariance matrix using all observations when computing BLUP; NULL if type = "lp" in cv_plmm.
#' @param ... Optional arguments to `predict_within_cv`
#'
#' @keywords internal
#'
#' @returns a list with three elements:
#' * a numeric vector with the loss at each value of lambda
#' * a numeric value indicating the number of lambda values used
#' * a numeric value with the predicted outcome (y hat) values at each lambda
#'
cvf <- function(i, fold, type, cv_args, estimated_Sigma, ...) {

  # save the 'prep' object from the plmm_prep() in cv_plmm
  full_cv_prep <- cv_args$prep

  # save outcome information -- will need to subset this into test and train sets
  y <- cv_args$y

  # make list to hold the data for this particular fold:
  fold_args <- list(std_X_details = list(),
                    fbm_flag = cv_args$fbm_flag,
                    penalty = cv_args$penalty,
                    penalty_factor = cv_args$penalty_factor,
                    gamma = cv_args$gamma,
                    alpha = cv_args$alpha,
                    nlambda = cv_args$nlambda,
                    lambda_min = cv_args$lambda_min,
                    max_iter = cv_args$max_iter,
                    eps = cv_args$eps,
                    warn = cv_args$warn,
                    convex = cv_args$convex,
                    dfmax = cv_args$dfmax,
                    lambda = cv_args$lambda)

  # subset std_X, U, and y to match fold indices ------------------------
  #   (and in so doing, leave out the ith fold)
  if (cv_args$fbm_flag) {

    # designate the training set
    train_X <- bigmemory::deepcopy(full_cv_prep$std_X,
                                   rows = which(fold!=i),
                                   type = "double",
                                   backingfile = paste0("train_fold",i,".bk"),
                                   descriptorfile = paste0("train_fold",i,".desc"),
                                   backingpath = bigmemory::dir.name(full_cv_prep$std_X))

    # TODO should these fold-specific data files be created in a tempdir? Then,
    #   the entire tempdir could be deleted at the end of this function....

    # create the copy of the training data to be standardized
    fold_args$std_X <- bigmemory::deepcopy(full_cv_prep$std_X,
                                           rows = which(fold!=i),
                                           type = "double",
                                           backingfile = paste0("std_train_fold",i,".bk"),
                                           descriptorfile = paste0("std_train_fold",i,".desc"),
                                           backingpath = bigmemory::dir.name(full_cv_prep$std_X))

    # re-scale data & check for singularity
    train_data <- .Call("big_std",
                        fold_args$std_X@address,
                        as.integer(count_cores()),
                        PACKAGE = "plmmr")

    fold_args$std_X@address <- train_data$std_X
    fold_args$std_X_details$center <- train_data$std_X_center
    fold_args$std_X_details$scale <- train_data$std_X_scale
    fold_args$std_X_details$ns <- which(train_data$std_X_scale > 1e-3)
    singular <- train_data$std_X_scale < 1e-3

    # do not fit a model on singular features!
    if (sum(singular) >= 1) fold_args$penalty_factor[singular] <- Inf

  } else {

    # subset training data
    train_X <- full_cv_prep$std_X[fold!=i, ,drop=FALSE]

    # Note: subsetting the data into test/train sets may cause low variance features
    #   to become constant features in the training data. The following lines address this issue

    # re-standardize training data & check for singularity

    std_info <- standardize_in_memory(train_X)
    fold_args$std_X <- std_info$std_X
    fold_args$std_X_details <- std_info$std_X_details
    # std_train_X <- ncvreg::std(train_X)
    # fold_args$std_X <- std_train_X
    # fold_args$std_X_details$center <- attr(std_train_X,'center')
    # fold_args$std_X_details$scale <- attr(std_train_X,'scale')
    # fold_args$std_X_details$ns <- attr(std_train_X,'nonsingular')

    # do not fit a model on these (near) singular features!

    singular <- fold_args$std_X_details$scale < 1e-3
    if (sum(singular) >= 1) fold_args$penalty_factor[singular] <- Inf
    # fold_args$penalty_factor <- fold_args$penalty_factor[fold_args$std_X_details$ns]

  }

  # subset outcome vector to include outcomes for training data only
  fold_args$y <- y[fold!=i]

  # center the training outcome
  fold_args$centered_y <- fold_args$y |>
    scale(scale=FALSE) |>
    drop()

  # extract test set --------------------------------------
  # this comes from cv prep on full data
  if (cv_args$fbm_flag){
    test_X <- bigmemory::deepcopy(full_cv_prep$std_X,
                                  rows = which(fold==i),
                                  type = "double",
                                  backingfile = paste0("test_fold",i,".bk"),
                                  descriptorfile = paste0("test_fold",i,".desc"),
                                  backingpath = bigmemory::dir.name(full_cv_prep$std_X))

  } else {
    test_X <- full_cv_prep$std_X[fold==i, , drop=FALSE]

    # use center/scale values from train_X to standardize test_X
    std_test_X <- scale(test_X,
                        center = fold_args$std_X_details$center,
                        scale = fold_args$std_X_details$scale)
  }

  # subset outcome for test set
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
                         penalty_factor = fold_args$penalty_factor,
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
                    penalty_factor = fold_args$penalty_factor,
                    fbm_flag = fold_args$fbm_flag,
                    penalty = fold_args$penalty,
                    gamma = fold_args$gamma,
                    alpha = fold_args$alpha,
                    lambda_min = fold_args$lambda_min,
                    nlambda = fold_args$nlambda,
                    lambda = fold_args$lambda,
                    eps = fold_args$eps,
                    max_iter = fold_args$max_iter,
                    warn = fold_args$warn,
                    convex = fold_args$convex,
                    dfmax = ncol(train_X) + 1)


  # first, get beta hat back on the scale of the training data
  # og_betas.i <- untransform(
  #   std_scale_beta = fit.i$std_scale_beta,
  #   p = ncol(train_X),
  #   std_X_details = fold_args$std_X_details,
  #   fbm_flag = fold_args$fbm_flag,
  #   use_names = FALSE)

  if(type == "lp"){
    # yhat <- predict_within_cv(fit = fit.i,
    #                           trainX = train_X,
    #                           testX = test_X,
    #                           og_scale_beta = og_betas.i,
    #                           type = 'lp',
    #                           fbm = cv_args$fbm_flag)

    # # a working idea... what if predictions were on the scale of the standardized data?
    yhat <- predict_within_cv(fit = fit.i,
                              trainX = NULL,
                              testX = std_test_X,
                              og_scale_beta = fit.i$std_scale_beta,
                              type = 'lp',
                              fbm = cv_args$fbm_flag)
  }

  if (type == 'blup'){
    # estimated_Sigma here comes from the overall fit in cv_plmm.R, an n*n matrix
    # Sigma_11 <- estimated_Sigma[fold!=i, fold!=i, drop = FALSE]
    # Sigma_21 <- estimated_Sigma[fold==i, fold!=i, drop = FALSE]

    # explicit calculation of Sigma_11 and Sigma_21
    Sigma_11 <- construct_variance(K = fold_prep$K, eta = fit.i$eta)
    Sigma_21 <- fit.i$eta*(1/ncol(train_X))*tcrossprod(std_test_X, fold_args$std_X)

    # yhat <- predict_within_cv(fit = fit.i,
    #                           trainX = train_X,
    #                           trainY = fold_args$y,
    #                           testX = test_X,
    #                           og_scale_beta = og_betas.i,
    #                           std_X_details = fold_args$std_X_details,
    #                           type = 'blup',
    #                           fbm = cv_args$fbm_flag,
    #                           Sigma_11 = Sigma_11,
    #                           Sigma_21 = Sigma_21, ...)

    yhat <- predict_within_cv(fit = fit.i,
                              trainX = train_X,
                              trainY = fold_args$y,
                              testX = std_test_X,
                              og_scale_beta = fit.i$std_scale_beta,
                              std_X_details = fold_args$std_X_details,
                              type = 'blup',
                              fbm = cv_args$fbm_flag,
                              Sigma_11 = Sigma_11,
                              Sigma_21 = Sigma_21, ...)

  }

  # cleanup -----------------------------------------------------------------
  # delete files created in cross-validation, if data is filebacked
  if (cv_args$fbm_flag) {
    list.files(path = bigmemory::dir.name(full_cv_prep$std_X),
               pattern = paste0("fold",i),
               full.names = TRUE) |> file.remove()
    gc() # release the pointer
  }

  # return -----------------------------------------------------------------
  loss <- sapply(1:ncol(yhat), function(ll) plmm_loss(test_y, yhat[,ll]))
  list(loss=loss, nl=length(fit.i$lambda), yhat=yhat)
}

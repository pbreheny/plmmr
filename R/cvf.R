#' Cross-validation internal function for `cv_plmm()`
#'
#' Internal function for `cv_plmm()` which calls `plmm()` on a fold subset of the original data.
#'
#' @param i       Fold number to be excluded from fit.
#' @param fold    n-length vector of fold-assignments.
#' @param type    A character argument indicating what should be returned from `predict.plmm()`. If `type = 'lp'` predictions are based on the linear predictor, \eqn{X \beta}.
#'                If `type = 'blup'`, predictions are based on the linear predictor plus the estimated random effect (BLUP).
#' @param cv_args List of additional arguments to be passed to plmm.
#' @param ...     Optional arguments to `predict_within_cv()`
#'
#' @return A list with three elements:
#' * `loss`: a numeric vector with the loss at each value of lambda
#' * `nl`: a numeric value indicating the number of lambda values used
#' * `yhat`: a numeric value with the predicted outcome values at each lambda
#'
#' @keywords internal
#'
cvf <- function(i, fold, type, cv_args, ...) {
  # save the 'prep' object from the plmm_prep() in cv_plmm
  full_cv_prep <- cv_args$prep

  # save outcome information -- will need to subset this into test and train sets
  y <- cv_args$y

  # make list to hold the data for this particular fold:
  fold_args <- list(
    std_X_details = list(),
    fbm_flag = cv_args$fbm_flag,
    plink_flag = cv_args$plink_flag,
    penalty = cv_args$penalty,
    gamma = cv_args$gamma,
    alpha = cv_args$alpha,
    nlambda = cv_args$nlambda,
    lambda_min = cv_args$lambda_min,
    max_iter = cv_args$max_iter,
    eps = cv_args$eps,
    dfmax = cv_args$dfmax,
    warn = cv_args$warn,
    lambda = cv_args$lambda,
    penalty_factor = full_cv_prep$penalty_factor
  )

  # If the K provided by the user was non-null, subset variance blocks
  if (!is.null(full_cv_prep$K)) {
    fold_args$K <- full_cv_prep$K[fold != i, fold != i]
    Sigma_21 <- full_cv_prep$K[fold == i, fold != i]
  }

  # subset std_X, U, and y to match fold indices ------------------------
  #   (and in so doing, leave out the ith fold)
  if (cv_args$fbm_flag) {
    # create the copy of the training data to be standardized
    fold_args$std_X <- bigmemory::deepcopy(
      full_cv_prep$std_X,
      rows = which(fold != i),
      type = "double",
      backingfile = paste0("std_train_fold", i, ".bk"),
      descriptorfile = paste0("std_train_fold", i, ".desc"),
      backingpath = bigmemory::dir.name(full_cv_prep$std_X)
    )

    # re-scale data & check for singularity
    std_trainX_info <- .Call(
      "big_std",
      fold_args$std_X@address,
      as.integer(count_cores()),
      to_center = TRUE,
      NULL,
      NULL,
      PACKAGE = "plmmr"
    )

    fold_args$std_X@address <- std_trainX_info$std_X
    fold_args$std_X_details$center <- std_trainX_info$std_X_center
    fold_args$std_X_details$scale <- std_trainX_info$std_X_scale
    fold_args$std_X_details$ns <- which(std_trainX_info$std_X_scale > 1e-3)
    singular <- std_trainX_info$std_X_scale < 1e-3

    if (any(singular)) {
      fold_args$penalty_factor[singular] <- Inf
    }
  } else {
    # subset training data
    train_X <- full_cv_prep$std_X[fold != i, , drop = FALSE]

    # Note: subsetting the data into test/train sets may cause low variance features
    #   to become constant features in the training data. The following lines address this issue

    # re-standardize training data & check for singularity
    std_info <- standardize_in_memory(train_X)
    fold_args$std_X <- std_info$std_X
    fold_args$std_X_details <- std_info$std_X_details

    # do not fit a model on these (near) singular features!
    singular <- fold_args$std_X_details$scale < 1e-3
    if (any(singular)) {
      fold_args$penalty_factor[singular] <- Inf
    }
  }

  # subset outcome vector to include outcomes for training data only
  fold_args$y <- cv_args$y[fold != i]

  # center the training outcome
  fold_args$centered_y <- fold_args$y |> scale(scale = FALSE) |> drop()
  # extract test set --------------------------------------
  # this comes from cv prep on full data
  if (cv_args$fbm_flag) {
    test_X <- bigmemory::deepcopy(
      full_cv_prep$std_X,
      rows = which(fold == i),
      type = "double",
      backingfile = paste0("test_fold", i, ".bk"),
      descriptorfile = paste0("test_fold", i, ".desc"),
      backingpath = bigmemory::dir.name(full_cv_prep$std_X)
    )
  } else {
    test_X <- full_cv_prep$std_X[fold == i, , drop = FALSE]
  }

  # subset outcome for test set
  test_y <- y[fold == i]

  # decomposition for current fold ------------------------------
  if (cv_args$prep$trace) {
    cat("Beginning eigendecomposition in fold ", i, ":\n")
  }

  fold_prep <- plmm_prep(
    std_X = fold_args$std_X,
    std_X_n = nrow(fold_args$std_X),
    std_X_p = ncol(fold_args$std_X),
    centered_y = fold_args$centered_y,
    fbm_flag = fold_args$fbm_flag,
    eta = cv_args$eta,
    penalty_factor = fold_args$penalty_factor,
    K = fold_args$K,
    trace = cv_args$prep$trace
  )
  fold_args$prep <- fold_prep
  # fit a plmm within each fold at each value of lambda
  # lambda stays the same for each fold; comes from the overall fit in cv_plmm()
  if (cv_args$prep$trace) {
    cat("** Fitting model in fold ", i, "\n", sep = "")
  }

  fit.i <- plmm_fit(
    prep = fold_prep,
    y = fold_args$y,
    std_X_details = fold_args$std_X_details,
    fbm_flag = fold_args$fbm_flag,
    penalty = fold_args$penalty,
    gamma = fold_args$gamma,
    alpha = fold_args$alpha,
    lambda_min = fold_args$lambda_min,
    nlambda = fold_args$nlambda,
    lambda = fold_args$lambda,
    eps = fold_args$eps,
    max_iter = fold_args$max_iter,
    dfmax = fold_args$dfmax,
    warn = fold_args$warn
  )

  format.i <- plmm_format(
    fit = fit.i,
    p = ncol(full_cv_prep$std_X),
    std_X_details = fold_args$std_X_details,
    fbm_flag = fold_args$fbm_flag,
    plink_flag = fold_args$plink_flag
  )

  # prediction ---------------------------------------------------
  # note: predictions are on the scale of the standardized training data
  if (type == "lp") {
    yhat <- predict_within_cv(fit = format.i, testX = test_X, type = "lp", fbm = cv_args$fbm_flag)
  }

  if (type == "blup") {
    if (cv_args$fbm_flag) {
      # we will need a copy of the testing data that is standardized
      std_test_X <- bigmemory::deepcopy(
        full_cv_prep$std_X,
        rows = which(fold == i),
        type = "double",
        backingfile = paste0("std_test_fold", i, ".bk"),
        descriptorfile = paste0("std_test_fold", i, ".desc"),
        backingpath = bigmemory::dir.name(full_cv_prep$std_X)
      )

      # use center/scale values from train_X to standardize test_X
      if (any(singular)) {
        fold_args$std_X_details$scale[singular] <- 1
      }
      std_test_info <- .Call(
        "big_std",
        std_test_X@address,
        as.integer(count_cores()),
        tocenter = TRUE,
        fold_args$std_X_details$center,
        fold_args$std_X_details$scale,
        PACKAGE = "plmmr"
      )
      std_test_X@address <- std_test_info$std_X

      const <- (fit.i$eta / ncol(fold_args$std_X))
      XXt <- bigalgebra::dgemm(TRANSA = "N", TRANSB = "T", A = std_test_X, B = fold_args$std_X)
      Sigma_21 <- const * XXt
      Sigma_21 <- Sigma_21[,] # convert to in-memory matrix
    } else if (!exists("Sigma_21")) {
      # don't rescale columns that were singular features in std_train_X;
      #   these features will have an estimated beta of 0 anyway
      if (any(singular)) {
        fold_args$std_X_details$scale[singular] <- 1
      }
      # use center/scale values from train_X to standardize test_X
      std_test_X <- scale(
        test_X,
        center = fold_args$std_X_details$center,
        scale = fold_args$std_X_details$scale
      )
      Sigma_21 <- (fit.i$eta / ncol(fold_args$std_X)) * tcrossprod(std_test_X, fold_args$std_X)
    }

    yhat <- predict_within_cv(
      fit = format.i,
      testX = test_X,
      type = "blup",
      fbm = cv_args$fbm_flag,
      Sigma_21 = Sigma_21
    )
  }

  # cleanup -----------------------------------------------------------------
  # delete files created in cross-validation, if data is filebacked
  if (cv_args$fbm_flag) {
    gc() # release the pointer
    list.files(
      path = bigmemory::dir.name(full_cv_prep$std_X),
      pattern = paste0("fold", i),
      full.names = TRUE
    ) |>
      unlink(force = TRUE)
  }

  # return -----------------------------------------------------------------
  loss <- sapply(seq_len(ncol(yhat)), function(ll) {
    plmm_loss(test_y, yhat[, ll])
  })
  list(loss = loss, nl = length(fit.i$lambda), yhat = yhat)
}

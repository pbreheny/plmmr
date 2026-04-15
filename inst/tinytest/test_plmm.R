# Tabitha Peter - tests for plmmr

# Test 0: Case where K = identity ----------------------------------------------

# set up
lambda0 <- c(1, 0.1, 0.01, 0.001)
admix_design <- create_design(X = admix$X, y = admix$y)
plmm0 <- plmm(
  design = admix_design,
  K = diag(nrow(admix$X)),
  lambda = lambda0,
  penalty = "lasso",
  trace = TRUE,
  eps = 1e-15)

lasso0 <- glmnet::glmnet(
  x = admix$X,
  y = admix$y,
  family = "gaussian",
  lambda = lambda0,
  thres = 1e-15)

A0 <- coef(plmm0)
dimnames(A0) <- NULL
B0 <- as.matrix(coef(lasso0))
dimnames(B0) <- NULL

# test 0 - implementation
expect_equivalent(A0, B0, tolerance = 0.01)


# Test 1 Case where K is diagonal and lambda is 0 ------------------------------

# use 'lm' here instead of glmnet

K_diagonal <- diag(x = (rnorm(n = nrow(admix$X))^2),
                   nrow = nrow(admix$X))

plmm1 <- plmm(
  design = admix_design,
  K = K_diagonal,
  # TODO: Need to fix plmm so that lambda can be a single value
  lambda = c(0.001, 0),
  penalty = "lasso",
  eps = 1e-15)

v1 <- diag(K_diagonal) * plmm1$eta + 1
print(summary(plmm1, lambda = 0))

A1 <- plmm1$beta_vals[,"0.0000"]
names(A1) <- NULL

lm1 <- lm(admix$y ~ admix$X,
          weights = 1/v1)

B1 <- lm1$coefficients

names(B1) <- NULL
B1 <- ifelse(is.na(B1), 0, B1)

# test 1: implementation
expect_equivalent(A1, B1, tolerance = 0.01)


# Test 2: Case where K is diagonal and lambda != 0 -----------------------------

lambda2 <- c(1, 0.1, 0.01)

plmm2 <- plmm(
  design = admix_design,
  K = K_diagonal,
  lambda = lambda2,
  penalty = "lasso",
  eps = 1e-15)

v2 <- diag(K_diagonal) * plmm2$eta + 1

lasso2 <- glmnet::glmnet(
  x = admix$X,
  y = admix$y,
  family = "gaussian",
  lambda = lambda2,
  # weights are by INVERSE variance
  weights = 1/v2,
  thres = 1e-15)


A2 <- as.matrix(plmm2$beta_vals[2:10, ])
dimnames(A2) <- NULL
B2 <- as.matrix(lasso2$beta[1:9, ]) # NB: glmnet() does not return intercept values
dimnames(B2) <- NULL

# test 2 - implementation
expect_equivalent(A2, B2, tolerance = 0.01)


# Test 3: show that monomorphic SNPs are given beta values of 0s ---------------

monomorphic <- apply(admix$X[,1:15], 2, var) == 0
monomorphic_snps <- paste0("Snp", which(monomorphic))
# NB: SNPs 8 and 14 are monomorphic
fit3 <- plmm(design = admix_design)
expect_equivalent(matrix(0,
                         nrow = length(monomorphic_snps),
                         ncol = length(fit3$lambda)
),
fit3$beta_vals[monomorphic_snps,])


# Test 4: make sure in-memory and filebacked computations match ----------------

local({
  temp_dir <- withr::local_tempdir() # using a temp dir -- change to fit your preference

  # process delimited files
  colon_dat <- process_delim(
    data_file = "colon2.txt",
    data_dir = find_example_data(parent = TRUE),
    rds_dir = temp_dir,
    rds_prefix = "processed_colon2",
    sep = "\t",
    overwrite = TRUE,
    header = TRUE)

  # prepare outcome data
  colon_outcome <- read.delim(find_example_data(path = "colon2_outcome.txt"))

  # filebacked
  fb_design <- create_design(
    data_file = colon_dat,
    rds_dir = temp_dir,
    new_file = "std_colon2",
    add_outcome = colon_outcome,
    outcome_id = "ID",
    outcome_col = "y",
    logfile = "fb_design",
    overwrite = TRUE)

  fb_fit <- plmm(
    design = fb_design,
    trace = TRUE,
    return_fit = TRUE)

  # in-memory
  colon_path <- find_example_data("colon2.txt")
  colon_X <- read.delim(colon_path)

  in_mem_design <- create_design(X = colon_X, y = colon_outcome$y)

  fit <- plmm(
    design = in_mem_design,
    lambda = fb_fit$lambda,
    K = fb_fit$fit$K,
    trace = TRUE,
    return_fit = TRUE)

  # check: these results match
  b1 <- fb_fit$beta_vals |> as.matrix()
  b2 <- fit$beta_vals
  expect_equivalent(b1, b2, tolerance = 0.025) # allowing slightly higher tolerance here; small coefficients
})


# Test 5: make sure predict method is working ----------------------------------

# plmm_fit <- plmm(design = admix_design,
#                  penalty = 'lasso',
#                  lambda = c(0.1, 0.01))
# plmm_pred <- predict(object = plmm_fit, newX = admix$X, type = "lp")
#
# # use glmnet as gold standard
# glmnet_fit <- glmnet::glmnet(admix$X, admix$y, lambda = c(0.1, 0.01))
# glmnet_pred <- predict(glmnet_fit, newx = admix$X, type = "response")
#
# test <- cbind(admix$y, plmm_pred, glmnet_pred)
# colnames(test) <- c('y',
#                     'y_hat_plmm0.1',
#                     'y_hat_plmm0.01',
#                     'y_hat_glmnet0.1',
#                     'y_hat_glmnet0.01')
#
# glmnet_spe0.1 <- crossprod(glmnet_pred[,1] - admix$y)/nrow(admix$X)
# plmmr_spe0.1 <- crossprod(plmm_pred[,1] - admix$y)/nrow(admix$X)
#
# glmnet_spe0.01 <- crossprod(glmnet_pred[,2] - admix$y)/nrow(admix$X)
# plmmr_spe0.01 <- crossprod(plmm_pred[,2] - admix$y)/nrow(admix$X)

# examine the values to see how much the two sets of predictions differ...
# mean(abs(test[,2] - test[,4]))
# test[1:10,]

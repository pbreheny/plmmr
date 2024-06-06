# Tabitha Peter - tests for plmmr

# Test 0: Case where K = identity -------------------------------------------

# set up
lambda0 <- c(1, 0.1, 0.01, 0.001)

plmm0 <- plmm(X = admix$X,
              y = admix$y,
              diag_K = TRUE,
              lambda = lambda0,
              penalty = "lasso",
              trace = TRUE)

lasso0 <- glmnet::glmnet(x = admix$X,
                 y = admix$y,
                 family = "gaussian",
                 lambda = lambda0)

A0 <- as.matrix(plmm0$beta_vals[2:10, ])
dimnames(A0) <- NULL
B0 <- as.matrix(lasso0$beta[1:9, ]) # NB: glmnet() does not return intercept values
dimnames(B0) <- NULL

# test 0 - implementation
expect_equivalent(A0, B0, tolerance = 0.01)


# Test 1 Case where K is diagonal and lambda is 0 ---------------------------

# use 'lm' here instead of glmnet

K_diagonal <- diag(x = (rnorm(n = nrow(admix$X))^2),
                   nrow = nrow(admix$X))

plmm1 <- plmm(X = admix$X,
              y = admix$y,
              K = K_diagonal,
              diag_K = TRUE,
              # FIXME: Need to fix plmm so that lambda can be a single value
              lambda = c(0.001, 0),
              penalty = "lasso")

v1 <- diag(K_diagonal)*plmm1$eta + 1
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

# check
# head(data.frame(A1, B1))

# Test 2: Case where K is diagonal and lambda != 0 -----------------------------

lambda2 <- c(1, 0.1, 0.01)

plmm2 <- plmm(X = admix$X,
              y = admix$y,
              diag_K = TRUE,
              K = K_diagonal,
              lambda = lambda2,
              penalty = "lasso")

v2 <- diag(K_diagonal)*plmm2$eta + 1

lasso2 <- glmnet::glmnet(x = admix$X,
                 y = admix$y,
                 family = "gaussian",
                 lambda = lambda2,
                 # weights are by INVERSE variance
                 weights = 1/v2)


A2 <- as.matrix(plmm2$beta_vals[2:10, ])
dimnames(A2) <- NULL
B2 <- as.matrix(lasso2$beta[1:9, ]) # NB: glmnet() does not return intercept values
dimnames(B2) <- NULL

# test 2 - implementation
expect_equivalent(A2, B2, tolerance = 0.1)

# Test 3: show that monomorphic SNPs are given beta values of 0s -------------
monomorphic <- apply(admix$X[,1:15], 2, var) == 0
monomorphic_snps <- paste0("Snp", which(monomorphic))
# NB: SNPs 8 and 14 are monomorphic
fit3 <- plmm(X = admix$X[,1:15], y = admix$y)
tinytest::expect_equivalent(matrix(0,
                 nrow = length(monomorphic_snps),
                 ncol = length(fit3$lambda)
                 ),
          fit3$beta_vals[monomorphic_snps,])

# Test 4: make sure in-memory and filebacked computations match ---------------
if (interactive()) {
  # filebacked fit
  process_plink(data_dir = find_example_data(parent = TRUE),
                prefix = "penncath_lite",
                gz = TRUE,
                outfile = "process_penncath",
                overwrite = TRUE,
                impute_method = "mode",
                keep_bigSNP = TRUE)

  my_fb_data <- paste0(find_example_data(parent = TRUE), "/penncath_lite")
  fb_fit <- plmm(X = my_fb_data,
                 returnX = FALSE,
                 trace = TRUE)

  # in memory fit
  pen <- bigsnpr::snp_attach(paste0(find_example_data(parent = TRUE), "/penncath_lite.rds"))
  inmem_fit <- plmm(X = pen$genotypes[pen$complete_phen,],
                    y = pen$fam$affection[pen$complete_phen],
                    trace = TRUE)

  tinytest::expect_equivalent(inmem_fit$beta_vals, as.matrix(fb_fit$beta_vals),
                              tolerance = 0.01)
}


# Test 5: make sure predict method is working -------------------
plmm_fit <- plmm(admix$X,
                 admix$y,
                 # K = relatedness_mat(admix$X),
                 penalty = 'lasso',
                 lambda = c(0.1, 0.01))
plmm_pred <- predict(object = plmm_fit, newX = admix$X, type = "lp")

# use glmnet as gold standard
glmnet_fit <- glmnet::glmnet(admix$X, admix$y, lambda = c(0.1, 0.01))
glmnet_pred <- predict(glmnet_fit, newx = admix$X, type = "response")

cbind(admix$y, plmm_pred, glmnet_pred) -> test
colnames(test) <- c('y',
                 'y_hat_plmm0.1',
                 'y_hat_plmm0.01',
                 'y_hat_glmnet0.1',
                 'y_hat_glmnet0.01')
# test[1:10,] # in plmm method, all rows of X have same predicted value!
if(abs(mean(test[,2] - test[,4])) > 5) stop("PLMM and GLMNET predictions are far off for the test model.")
# NB: the 5 above is chosen arbitrarily, based on my experience with the admix data


# Test 6: is resid. method working ----------------------------------------------
R <- residuals(object = plmm(admix$X, admix$y, penalty = "lasso",
                                  diag_K = TRUE, lambda = lambda0))

ncv_fit <- ncvreg::ncvreg(X = admix$X, y = admix$y, penalty = "lasso", lambda = lambda0)
ncv_R <- matrix(nrow = nrow(ncv_fit$linear.predictors), ncol = ncol(ncv_fit$linear.predictors))
for(j in 1:ncol(ncv_R)){
  ncv_R[,j] <- ncv_fit$y - ncv_fit$linear.predictors[j]
}

tinytest::expect_equivalent(R, ncv_R)

# Test 7: make sure plmm() runs in-memory and filebacked ---------------------

if (interactive()) {
  # process data
  pen_clinic <- read.csv(paste0(find_example_data(parent = TRUE), "/penncath_clinical.csv"))
  extdata <- pen_clinic[,3:4]
  rownames(extdata) <- pen_clinic$FamID # This is important!

  # create a new temporary directory
  temp_dir <- paste0(tempdir(), sample(LETTERS, 1))

  process_plink(data_dir = find_example_data(parent = TRUE),
                rds_dir = temp_dir, # using a temporary directory
                prefix = "penncath_lite",
                id_var = "FID", # this is KEY!
                outfile = "process_penncath",
                impute_method = "mode",
                add_predictor_ext = extdata)

  # filebacked
  dat_plus_newvars <- paste0(temp_dir, "/std_penncath_lite")
  pen2 <- readRDS(paste0(temp_dir, "/std_penncath_lite.rds"))
  foo <- plmm(X = dat_plus_newvars,
                           penalty.factor = c(0, 0, rep(1, ncol(pen2$std_X) - 2)),
                           returnX = FALSE,
                           trace = TRUE)

  # NB: returnX = FALSE is needed to pass to get_data(); otherwise, this
  #   will run in-memory because of the small size of this test data set
  foo_nz <- which(foo$beta_vals[,4] != 0)

  # in memory
  foo2 <- plmm(X = dat_plus_newvars,
                           penalty.factor = c(0, 0, rep(1, ncol(pen2$std_X) - 2)),
                           trace = TRUE)
  foo2_nz <- which(foo2$beta_vals[,4] != 0)

  # look at head of values (checks SNP names)
  foo$beta_vals[,4] |> head()
  foo2$beta_vals[,4] |> head()

  # look at just nonzero values
  foo$beta_vals[foo_nz,4] |> head()
  foo2$beta_vals[foo2_nz,4] |> head()

  tinytest::expect_equivalent(foo$beta_vals[,3], foo2$beta_vals[,3],
                              tolerance = 0.001)
  tinytest::expect_equivalent(foo$beta_vals[,4], foo2$beta_vals[,4],
                              tolerance = 0.001)

}

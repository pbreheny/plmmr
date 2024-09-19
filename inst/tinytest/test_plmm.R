# Tabitha Peter - tests for plmmr

# Test 0: Case where K = identity -------------------------------------------

# set up
lambda0 <- c(1, 0.1, 0.01, 0.001)
admix_design <- create_design(X = admix$X, outcome_col = admix$y)
plmm0 <- plmm(design = admix_design,
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
tinytest::expect_equivalent(A0, B0, tolerance = 0.01)


# Test 1 Case where K is diagonal and lambda is 0 ---------------------------

# use 'lm' here instead of glmnet

K_diagonal <- diag(x = (rnorm(n = nrow(admix$X))^2),
                   nrow = nrow(admix$X))

plmm1 <- plmm(design = admix_design,
              K = K_diagonal,
              diag_K = TRUE,
              # TODO: Need to fix plmm so that lambda can be a single value
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
tinytest::expect_equivalent(A1, B1, tolerance = 0.01)

# check
# head(data.frame(A1, B1))

# Test 2: Case where K is diagonal and lambda != 0 -----------------------------

lambda2 <- c(1, 0.1, 0.01)

plmm2 <- plmm(design = admix_design,
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
tinytest::expect_equivalent(A2, B2, tolerance = 0.1)

# Test 3: show that monomorphic SNPs are given beta values of 0s -------------
monomorphic <- apply(admix$X[,1:15], 2, var) == 0
monomorphic_snps <- paste0("Snp", which(monomorphic))
# NB: SNPs 8 and 14 are monomorphic
fit3 <- plmm(design = admix_design)
tinytest::expect_equivalent(matrix(0,
                 nrow = length(monomorphic_snps),
                 ncol = length(fit3$lambda)
                 ),
          fit3$beta_vals[monomorphic_snps,])

# Test 4: make sure in-memory and filebacked computations match ---------------
if (interactive()) {
  # process delimited files
  temp_dir <- tempdir() # using a temp dir -- change to fit your preference
  colon_dat <- process_delim(data_file = "colon2.txt",
                             data_dir = find_example_data(parent = TRUE),
                             rds_dir = temp_dir,
                             rds_prefix = "processed_colon2",
                             sep = "\t",
                             overwrite = TRUE,
                             header = TRUE)
  # prepare outcome data
  colon_outcome <- read.delim(find_example_data(path = "colon2_outcome.txt"))

  # create a design
  colon_design <- create_design(data_file = colon_dat,
                                rds_dir = temp_dir,
                                new_file = "std_colon2",
                                add_outcome = colon_outcome,
                                outcome_id = "ID",
                                outcome_col = "y",
                                logfile = "colon_design",
                                overwrite = TRUE)
  # filebacked
  fb_fit <- plmm(design = colon_design, trace = TRUE, return_fit = TRUE)

  # in-memory
  colon_X <- read.delim(file = "inst/extdata/colon2.txt")
  in_mem_design <- create_design(X = colon_X, outcome_col = colon_outcome$y)
  fit <- plmm(design = in_mem_design,
              # make sure to use the same K!
              K = fb_fit$K,
              trace = TRUE)

  # check: these results match
  b1 <- fb_fit$beta_vals |> as.matrix()
  b2 <- fit$beta_vals
  tinytest::expect_equivalent(b1, b2, tolerance = 0.001) # passes

}


# Test 5: make sure predict method is working -------------------
plmm_fit <- plmm(design = admix_design,
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

if(abs(mean(test[,2] - test[,4])) > 0.01) stop("PLMM and GLMNET predictions do not align well.")
# examine the values to see how much the two sets of predictions differ...
# test[1:10,]

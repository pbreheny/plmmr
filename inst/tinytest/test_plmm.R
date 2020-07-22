
set.seed(7)
Data <- genData(5, 2, 1) # can change this to check more high-dimensional data
X <- Data$X
y <- Data$y

dont_run <- FALSE

### no rotation checks ------------------------------------------------------###

# std(X)
plmm1 <- plmm(ncvreg::std(X),
              y,
              penalty = "lasso",
              alpha = 1,
              nlambda = 100,
              eps = 1e-04,
              max.iter = 1000,
              convex = FALSE,
              warn = TRUE,
              standardize = FALSE,
              rotation = FALSE,
              returnX = FALSE)

ncv1 <- ncvreg::ncvreg(ncvreg::std(X), y, "gaussian", penalty = "lasso", lambda = plmm1$lambda)
glm1 <- glmnet::glmnet(ncvreg::std(X), y, "gaussian", standardize = TRUE, lambda = plmm1$lambda)

expect_equivalent(plmm1$beta, coef(ncv1), tol = 1e-3)
expect_equivalent(plmm1$beta, as.matrix(coef(glm1)), tol = 1e-3)

# X
plmm2 <- plmm(X,
              y,
              penalty = "lasso",
              alpha = 1,
              nlambda = 100,
              eps = 1e-04,
              max.iter = 1000,
              convex = FALSE,
              warn = TRUE,
              standardize = TRUE,
              rotation = FALSE,
              returnX = FALSE)

ncv2 <- ncvreg::ncvreg(X, y, "gaussian", penalty = "lasso", lambda = plmm2$lambda)
glm2 <- glmnet::glmnet(X, y, "gaussian", standardize = TRUE, lambda = plmm2$lambda)

expect_equivalent(plmm2$beta, coef(ncv2), tol = 1e-3)
expect_equivalent(plmm2$beta, as.matrix(coef(glm2)), tol = 1e-3)

### add check for unpenalized covars here with additional equiv.check - modify penalty.factor for glmnet


### rotation checks ---------------------------------------------------------###

# std(X)
plmm3 <- plmm(ncvreg::std(X),
              y,
              penalty = "lasso",
              alpha = 1,
              lambda = 0, # compare to ols solutions
              eps = 1e-4,
              max.iter = 1e3,
              convex = FALSE,
              warn = TRUE,
              standardize = FALSE,
              rotation = TRUE,
              returnX = TRUE)

if (nrow(X) > ncol(X)){
  ols_soln <- solve(t(plmm3$X) %*% plmm3$X) %*% t(plmm3$X) %*% plmm3$y
  expect_equivalent(plmm3$beta, ols_soln, tol = 1e-3)
}

# X
plmm4 <- plmm(X,
              y,
              penalty = "lasso",
              alpha = 1,
              lambda = 0, # compare to ols solutions
              eps = 1e-4,
              max.iter = 1e3,
              convex = FALSE,
              warn = TRUE,
              standardize = FALSE,
              rotation = TRUE,
              returnX = TRUE)

if (nrow(X) > ncol(X)){
  ols_soln <- lm(y ~ X)
  expect_equivalent(as.numeric(plmm4$beta), coef(ols_soln), tol = 1e-3)
}


### I don't think this is going to work in general for glmnet - it's not the same for a manual int with int = F
# https://stackoverflow.com/questions/49495494/glmnet-is-different-with-intercept-true-compared-to-intercept-false-and-with-pen
# glm3 <- glmnet::glmnet(plmm3$X, plmm3$y, "gaussian", standardize = FALSE,
#                        intercept = FALSE,
#                        penalty.factor = rep(c(0, 1), c(1, ncol(X))), lambda = 0)
# coef(glm3)[-1,, drop = FALSE]
#
# expect_equivalent(as.matrix(plmm3$beta), as.matrix(coef(glm3)[-1,]), tol = 1e-3)

if (dont_run){
  fit <- glmnet::glmnet(X, y, "gaussian", standardize = FALSE, intercept = TRUE, nlambda = 5)
  fit1 <- glmnet::glmnet(cbind(1, X), y, "gaussian", standardize = FALSE, intercept = FALSE, penalty.factor = c(0, rep(1, ncol(X))), nlambda = 5)
  coef(fit)
  coef(fit1)
}


### can't do a rotated and standardized version check - standardization occurs
### for the unrotated version of X, glmnet would give the standardization for the
### rotated version, which will not be the same.
### Using the unrotated test to check that standardization is working as it should be

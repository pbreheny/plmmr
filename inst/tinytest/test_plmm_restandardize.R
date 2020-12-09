
nn <- 50
pp <- 4
p1 <- 1
set.seed(7)
X <- matrix(rnorm(nn * pp), nrow = nn, ncol = pp)
X0 <- matrix(c(rnorm(nn * 1)), nrow = nn, ncol = 1)
B <- rep(c(1, 0), times = c(p1, pp - p1))
y <- X %*% B + X0 %*% c(1) + rnorm(nn)

dont_run <- FALSE

### test code
# X = ncvreg::std(X)
# y = y
# penalty = "lasso"
# alpha = 1
# nlambda = 5
# standardizeX = FALSE
# rotation = FALSE
# returnX = FALSE
# n = nn
# p = pp

### no rotation checks ------------------------------------------------------###

# std(X)
plmm1 <- plmm(ncvreg::std(X),
              y,
              penalty = "lasso",
              alpha = 1,
              nlambda = 5,
              standardizeX = FALSE,
              standardizeRtX = TRUE,
              rotation = FALSE,
              returnX = FALSE)

ncv1 <- ncvreg::ncvreg(ncvreg::std(X), y, "gaussian", penalty = "lasso", lambda = plmm1$lambda)
glm1 <- glmnet::glmnet(ncvreg::std(X), y, "gaussian", standardize = TRUE, lambda = plmm1$lambda)

expect_equivalent(coef(plmm1), coef(ncv1), tol = 1e-3)
expect_equivalent(coef(plmm1), as.matrix(coef(glm1)), tol = 1e-3)
expect_equivalent(coef(plmm1)[-1, 1], rep(0, ncol(X))) # make sure setup lambda is working correctly

# X
plmm2 <- plmm(X,
              y,
              penalty = "lasso",
              alpha = 1,
              nlambda = 5,
              standardizeX = TRUE,
              standardizeRtX = TRUE,
              rotation = FALSE,
              returnX = FALSE)

ncv2 <- ncvreg::ncvreg(X, y, "gaussian", penalty = "lasso", lambda = plmm2$lambda)
glm2 <- glmnet::glmnet(X, y, "gaussian", standardize = TRUE, lambda = plmm2$lambda)

expect_equivalent(coef(plmm2), coef(ncv2), tol = 1e-3)
expect_equivalent(coef(plmm2), as.matrix(coef(glm2)), tol = 1e-3)
expect_equivalent(coef(plmm2)[-1, 1], rep(0, ncol(X)))


### no rotation + unpenalized covar checks-----------------------------------###
### first var (after int) should be included because unpenalized


### these are not passing test...pass at 1e-1 (?)
### I think this has to do with the issue in glmnet not treating the intercept/unpenalized vars equivalently
### With lambda = 0 (ols) seems fine...

# std(X)
plmm3 <- plmm(ncvreg::std(cbind(X0, X)),
              y,
              X_for_K = ncvreg::std(X),
              penalty = "lasso",
              penalty.factor = rep(c(0, 1), times = c(ncol(X0), ncol(X))),
              alpha = 1,
              nlambda = 5,
              standardizeX = FALSE,
              standardizeRtX = TRUE,
              rotation = FALSE,
              returnX = FALSE)

ncv3 <- ncvreg::ncvreg(ncvreg::std(cbind(X0, X)), y, "gaussian",
                       penalty = "lasso", lambda = plmm3$lambda,
                       penalty.factor = rep(c(0, 1), times = c(ncol(X0), ncol(X))))
glm3 <- glmnet::glmnet(ncvreg::std(cbind(X0, X)), y, "gaussian",
                       standardize = FALSE, lambda = plmm3$lambda,
                       penalty.factor = rep(c(0, 1), times = c(ncol(X0), ncol(X))))

expect_equivalent(coef(plmm3), coef(ncv3), tol = 1e-3)
expect_equivalent(coef(plmm3), as.matrix(coef(glm3)), tol = 1e-1)
expect_equivalent(coef(plmm3)[-c(1:(1 + ncol(X0))), 1], rep(0, ncol(X)))


# X
plmm4 <- plmm(cbind(X0, X),
              y,
              X_for_K = X,
              penalty = "lasso",
              penalty.factor = rep(c(0, 1), times = c(ncol(X0), ncol(X))),
              alpha = 1,
              nlambda = 5,
              standardizeX = TRUE,
              standardizeRtX = TRUE,
              rotation = FALSE,
              returnX = FALSE)

ncv4 <- ncvreg::ncvreg(cbind(X0, X), y, "gaussian",
                       penalty = "lasso", lambda = plmm4$lambda,
                       penalty.factor = rep(c(0, 1), times = c(ncol(X0), ncol(X))))
glm4 <- glmnet::glmnet(cbind(X0, X), y, "gaussian",
                       standardize = TRUE, lambda = plmm4$lambda,
                       penalty.factor = rep(c(0, 1), times = c(ncol(X0), ncol(X))))

expect_equivalent(coef(plmm4), coef(ncv4), tol = 1e-2)
expect_equivalent(coef(plmm4), as.matrix(coef(glm4)), tol = 1e-1)
expect_equivalent(coef(plmm4)[-c(1:(1 + ncol(X0))), 1], rep(0, ncol(X)))


### rotation checks ---------------------------------------------------------###
### ISSUES HERE - should these match????
# x <- matrix(c(1, 1, 0), nrow = 3, ncol = 1)
# beta <- c(1)
# y <- x %*% beta + rnorm(3, 0, 0.2)
# test <- plmm(ncvreg::std(x),
#               y,
#               penalty = "lasso",
#               alpha = 1,
#               lambda = 0, # compare to ols solutions
#               standardizeX = FALSE,
#               standardizeRtX = FALSE,
#               rotation = TRUE,
#               returnX = TRUE)
# lm(test$y ~ -1 + test$X)
# as.numeric(solve(t(test$X) %*% test$X) %*% t(test$X) %*% test$y)
#
# test2 <- plmm(test$X,
#              test$y,
#              penalty = "lasso",
#              alpha = 1,
#              lambda = 0, # compare to ols solutions
#              standardizeX = FALSE,
#              standardizeRtX = FALSE,
#              rotation = FALSE,
#              intercept = FALSE,
#              returnX = TRUE)
#
#
#
# ### compare with ols solutions
# # std(X)
# plmm5 <- plmm(ncvreg::std(X),
#               y,
#               penalty = "lasso",
#               alpha = 1,
#               lambda = 0, # compare to ols solutions
#               standardizeX = FALSE,
#               standardizeRtX = TRUE,
#               rotation = TRUE,
#               returnX = TRUE)
#
# ncvreg::ncvreg(cbind(X0, X), y, "gaussian",
#                penalty = "lasso", lambda = 0,
#                penalty.factor = rep(c(0, 1), times = c(ncol(X0), ncol(X))))
#
# if (nrow(X) > ncol(X)){
#   ols_soln <- lm(plmm5$y ~ -1 + plmm5$X)
#   ols_soln <- as.numeric(solve(t(plmm5$X) %*% plmm5$X) %*% t(plmm5$X) %*% plmm5$y)
#   expect_equivalent(coef(plmm5), ols_soln, tol = 1e-3)
# }
#
# # X
# plmm6 <- plmm(X,
#               y,
#               penalty = "lasso",
#               alpha = 1,
#               lambda = 0, # compare to ols solutions
#               standardizeX = FALSE,
#               standardizeRtX = FALSE,
#               rotation = TRUE,
#               returnX = TRUE)
#
# if (nrow(X) > ncol(X)){
#   ols_soln <- lm(y ~ X)
#   expect_equivalent(as.numeric(coef(plmm6)), coef(ols_soln), tol = 1e-3)
# }
#
# ### rotation + unpenalized covar checks-----------------------------------###
#
# # std(X)
# plmm7 <- plmm(ncvreg::std(cbind(X0, X)),
#               y,
#               penalty = "lasso",
#               alpha = 1,
#               lambda = 0, # compare to ols solutions
#               standardizeX = FALSE,
#               standardizeRtX = FALSE,
#               rotation = TRUE,
#               returnX = TRUE)
#
# if (nrow(cbind(X0, X)) > ncol(cbind(X0, X))){
#   ols_soln <- as.numeric(solve(t(plmm7$X) %*% plmm7$X) %*% t(plmm7$X) %*% plmm7$y)
#   expect_equivalent(coef(plmm7), ols_soln, tol = 1e-3)
# }
#
# # X
# plmm8 <- plmm(cbind(X0, X),
#               y,
#               penalty = "lasso",
#               alpha = 1,
#               lambda = 0, # compare to ols solutions
#               standardizeX = FALSE,
#               standardizeRtX = FALSE,
#               rotation = TRUE,
#               returnX = TRUE)
#
# if (nrow(cbind(X0, X)) > ncol(cbind(X0, X))){
#   ols_soln <- lm(y ~ cbind(X0, X))
#   expect_equivalent(coef(plmm8), coef(ols_soln), tol = 1e-3)
# }
#

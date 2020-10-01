

nn <- 50
pp <- 4
p1 <- 1
set.seed(7)
X <- matrix(rnorm(nn * pp), nrow = nn, ncol = pp)
X0 <- matrix(c(rnorm(nn * 1)), nrow = nn, ncol = 1)
B <- rep(c(1, 0), times = c(p1, pp - p1))
y <- X %*% B + X0 %*% c(1) + rnorm(nn)

### 1-fold cv.plmm should match plmm ----------------------------------------###

plmm0 <- plmm(X,
              y,
              penalty = "lasso",
              alpha = 1,
              nlambda = 5,
              standardize = TRUE,
              rotation = TRUE,
              returnX = TRUE)

cv_plmm0 <- cv.plmm(X,
                   y,
                   penalty = "lasso",
                   alpha = 1,
                   lambda = plmm0$lambda,
                   standardize = TRUE,
                   rotation = TRUE,
                   returnX = TRUE,
                   nfolds = 2)

cv_plmm1 <- cv.plmm(X,
                    y,
                    penalty = "lasso",
                    alpha = 1,
                    nlambda = 5,
                    standardize = TRUE,
                    rotation = TRUE,
                    returnX = TRUE,
                    nfolds = 2)

expect_equivalent(coef(plmm0), coef(cv_plmm0$fit), tol = 1e-5) # same overall fit?
expect_equivalent(plmm0$lambda, cv_plmm1$lambda, tol = 1e-5) # are they both computing lambda sequences the same way



cv_plmm2 <- cv.plmm(X,
                    y,
                    penalty = "lasso",
                    alpha = 1,
                    nlambda = 100,
                    standardize = TRUE,
                    rotation = TRUE,
                    returnX = TRUE,
                    nfolds = 10)
plot(cv_plmm2)
plot(cv_plmm2, type = 'all')

### need to get cv_plmm working with (1) lambda 1se and (2) parallel package for bigger data sets
# (3) havea a testing option for it to return the fold-specific coefficient values




### is cv getting the correct coefficients?
# plmm1 <- plmm(plmm0$X[1:(nn/2),],
#               plmm0$y[1:(nn/2)],
#               penalty = "lasso",
#               alpha = 1,
#               lambda = plmm0$lambda,
#               standardize = FALSE,
#               intercept = FALSE,
#               rotation = FALSE,
#               returnX = TRUE)
#
# plmm2 <- plmm(plmm0$X[((nn/2) + 1):nn,],
#               plmm0$y[((nn/2) + 1):nn],
#               penalty = "lasso",
#               alpha = 1,
#               lambda = plmm0$lambda,
#               standardize = FALSE,
#               intercept = FALSE,
#               rotation = FALSE,
#               returnX = TRUE)

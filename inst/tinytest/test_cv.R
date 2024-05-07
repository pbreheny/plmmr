# Test 1 - sanity check -----------------------------------------
# 'blup' and 'response' should not give the same result 
cv1 <- cv.plmm(X = admix$X, y = admix$y, seed = 321, type = 'lp')
cv2 <- cv.plmm(X = admix$X, y = admix$y, seed = 321, type='blup')
if(identical(cv1$cve, cv2$cve))stop("BLUP and response types give same result")


# Test 2 (think of what to add here) -------------------------
# cv3 <- cv.plmm(X = admix$X, y = admix$y, k = 70, trace = TRUE)
# cv4 <- cv.plmm(X = admix$X, y = admix$y, k = 70, type = "blup", trace = TRUE)
# identical(cv3$cve, cv4$cve)


# Test 3: BLUP method for prediction -------------------------

# look at predict.plmm(..., type = 'blup')

# #################################
# # used to check BLUP calculation #
# #################################
# library(zeallot)
# aggregate X and newX to compute V
# X_all <- rbind(std_X, std_newX)
# c(S_all, U_all) %<-% svd(X_all, nv = 0) # D, U
# S_all <- S_all^2 / p

# assuming newX has the same eta as X
# eta_all <- object$eta
# Vhat_all <- eta_all * tcrossprod(U_all %*% diag(S_all), U_all) + (1-eta_all)*diag(nrow(U_all))
# V21_check <- Vhat_all[-c(1:n), 1:n, drop = FALSE]
# V11_check <- Vhat_all[1:n, 1:n, drop = FALSE]

# ranef_check <- V21_check %*% chol2inv(chol(V11_check)) %*% (drop(y) - cbind(1, std_X) %*% beta_vals)
# print(eta)

# blup_check <- Xb + ranef_check

# frob_distance <- (Matrix::norm(blup_check) - Matrix::norm(blup))/Matrix::norm(blup_check)
# if (abs(frob_distance) > 0.01) stop("\nThe two calculations of the BLUP do not match")


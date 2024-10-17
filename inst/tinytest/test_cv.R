admix_design <- create_design(X = admix$X, y = admix$y)
# Test 1 - sanity check -----------------------------------------
# 'blup' and 'response' should not give the same result
cv1 <- cv_plmm(design = admix_design, seed = 321, type = 'lp')
cv2 <- cv_plmm(design = admix_design, seed = 321, type='blup')
if(identical(cv1$cve, cv2$cve))stop("BLUP and response types give same result")


# Test 2 (think of what to add here) -------------------------

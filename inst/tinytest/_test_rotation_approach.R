# TKP 
# May 2024
# Goal: figure out how to make rotate_filebacked() run more efficiently

admix_fbm <- as_FBM(admix$X)
K <- relatedness_mat(admix_fbm, fbm = TRUE)
eig <- eigen(K[,])
U <- eig$vectors
w <- eig$values

# question: can I alter the way I am constructing WUt?
wUt <- sweep(x = t(prep$U), MARGIN = 1, STATS = w, FUN = "*")
wU <- sweep(x = U, MARGIN = 2, STATS = w, FUN = "*")
identical(t(wU), wUt) # yes! this works 

# question: can I multiply a matrix wUt with a filebacked X? 
# rot_X <- FBM(nrow = nrow(admix_fbm), ncol = ncol(admix_fbm), init = 0)
rot_X <- crossprod(wU, admix_fbm) # works, but creates a matrix in-memory... :( 

# question: ok, since bigstatsr multiplication methods don't return an FBM (which is what I need),
# can I do the rotation in C and make the C code return a big.matrix?








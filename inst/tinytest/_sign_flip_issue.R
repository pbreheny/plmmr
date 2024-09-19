# TKP
# Sept. 2024
# Objective: using the 'colon2' data as an example, I want to illustrate
#   that the in-memory method for eigendecomposition may produce calculations
#   which differ from the filebacked method up to the +/- signs of the vectors of U

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
# filebacked fit
fb_fit <- plmm(design = colon_design, trace = TRUE, return_fit = TRUE)

# in-memory fit
colon_X <- read.delim(file = "inst/extdata/colon2.txt")
in_mem_design <- create_design(X = colon_X, outcome_col = colon_outcome$y)
fit <- plmm(design = in_mem_design,
            # make sure to use the same K
            # K = fb_fit$K,
            trace = TRUE)

# check: these results match
b1 <- fb_fit$beta_vals |> as.matrix()
b2 <- fit$beta_vals
tinytest::expect_equivalent(b1, b2) # does not pass

# compare S -- they align, because the eigenvalues must be positive
tinytest::expect_equivalent(fb_fit$K$s, fit$K$s)
# ... but the columns of U do not match
tinytest::expect_equivalent(fb_fit$K$U, fit$K$U)

# plot filebacked K
fb_K <- tcrossprod(fb_fit$K$U%*%diag(fb_fit$K$s), fb_fit$K$U)
corrplot::corrplot(fb_K, is.corr = F)

# plot in-memory K and compare to filebacked version
K <- tcrossprod(fit$K$U%*%diag(fit$K$s), fit$K$U)
corrplot::corrplot(K, is.corr = F, col = corrplot::COL2("PuOr", 200))

# examine signs
fb_U_sign <- apply(fb_fit$K$U, 2, FUN = sign)
U_sign <- apply(fit$K$U, 2, FUN = sign)
tinytest::expect_equivalent(U_sign, fb_U_sign)

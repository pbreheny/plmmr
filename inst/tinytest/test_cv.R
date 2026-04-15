# Test 1 - make sure in-memory and filebacked LOOCV match ----------------------

local({
  lambda <- exp(seq(log(1), log(0.05), length.out = 25))
  # process delimited files
  temp_dir <- withr::local_tempdir() # using a temp dir -- change to fit your preference
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
  n <- nrow(colon_outcome)

  # create a design
  fb_design <- create_design(
    data_file = colon_dat,
    rds_dir = temp_dir,
    new_file = "std_colon2",
    add_outcome = colon_outcome,
    outcome_id = "ID",
    outcome_col = "y",
    logfile = "fb_design",
    overwrite = TRUE)
  # filebacked
  fb_fit <- cv_plmm(
    design = fb_design,
    lambda = lambda,
    trace = TRUE,
    return_fit = TRUE,
    nfolds = n,
    eps = 1e-15,
    warn = FALSE) # need for small epsilon

  # in-memory
  colon_path <- find_example_data("colon2.txt")
  colon_X <- read.delim(colon_path)

  in_mem_design <- create_design(X = colon_X, y = colon_outcome$y)

  fit <- cv_plmm(
    design = in_mem_design,
    lambda = lambda,
    K = fb_fit$K,
    trace = TRUE,
    return_fit = TRUE,
    nfolds = n,
    eps = 1e-15,
    warn = FALSE)

  # check: these results match
  b1 <- coef(fb_fit) |> as.numeric()
  b2 <- coef(fit)
  expect_equivalent(b1, b2, tolerance = 0.01)
})

# Test 1 - make sure in-memory and filebacked LOOCV match ----------------------

local({
  temp_dir <- withr::local_tempdir() # using a temp dir -- change to fit your preference
  # log-transform the colon data
  colon_path <- find_example_data("colon2.txt")
  colon_X <- read.delim(colon_path)
  colon_X[,-1] <- log(colon_X[,-1])
  write.table(colon_X,
              file.path(temp_dir, "colon2_log.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  # process delimited files
  colon_dat <- process_delim(
    data_file = "colon2_log.txt",
    data_dir = temp_dir,
    rds_dir = temp_dir,
    rds_prefix = "processed_colon2",
    sep = "\t",
    overwrite = TRUE,
    header = TRUE)

  # prepare outcome data
  colon_outcome <- read.delim(find_example_data(path = "colon2_outcome.txt"))
  n <- nrow(colon_outcome)

  # filebacked
  fb_design <- create_design(
    data_file = colon_dat,
    rds_dir = temp_dir,
    new_file = "log_colon2",
    add_outcome = colon_outcome,
    outcome_id = "ID",
    outcome_col = "y",
    logfile = "fb_design",
    overwrite = TRUE)

  fb_fit <- cv_plmm(
    design = fb_design,
    trace = TRUE,
    return_fit = TRUE,
    nfolds = n,
    warn = FALSE) # need for small epsilon

  # in-memory
  in_mem_design <- create_design(X = colon_X, y = colon_outcome$y)

  fit <- cv_plmm(
    design = in_mem_design,
    lambda = fb_fit$lambda,
    K = fb_fit$fit$K,
    trace = TRUE,
    return_fit = TRUE,
    fold = fb_fit$fold,
    warn = FALSE)

  # check: these results match
  b1 <- coef(fb_fit) |> as.numeric()
  b2 <- coef(fit)
  expect_equivalent(b1, b2, tolerance = 0.01)
})

local({
  temp_dir <- withr::local_tempdir()
  unzip_example_data(outdir = temp_dir)

  # test 1: can all 3 pieces be combined correctly?  ----------------------
  # use external pheno and predictor files, no NAs

  ## create design matrix-------------------------------------------------------
  penncath_pheno <- read.csv(find_example_data(path = "penncath_clinical.csv"))

  predictors <- penncath_pheno[, c("FamID", "age", "tg")]
  predictors$tg <- ifelse(is.na(predictors$tg), mean(predictors$tg, na.rm = TRUE), predictors$tg)
  predictors$FamID <- as.character(predictors$FamID)

  # note: CAD has no missing values in phen file
  phen <- cbind(penncath_pheno$FamID, penncath_pheno$CAD) |> as.data.frame()

  colnames(phen) <- c("FamID", "CAD")
  phen$FamID <- as.character(phen$FamID)

  penncath_lite <- process_plink(
    data_dir = temp_dir,
    data_prefix = "penncath_lite",
    rds_dir = temp_dir,
    rds_prefix = "process_penncath",
    id_var = "FID",
    quiet = TRUE,
    overwrite = TRUE,
    parallel = FALSE
  )

  X <- create_design(
    data_file = penncath_lite,
    feature_id = "FID",
    rds_dir = temp_dir,
    new_file = "std_penncath_lite",
    add_outcome = phen,
    outcome_id = "FamID",
    outcome_col = "CAD",
    add_predictor = predictors,
    predictor_id = 'FamID',
    overwrite = TRUE,
    quiet = TRUE
  )

  res <- readRDS(X)

  ### checks  -------------------------------------------------------------------
  # are external files and PLINK fam file aligned on ID var?
  plink <- readRDS(file.path(temp_dir, "process_penncath.rds"))
  expect_identical(as.character(plink$fam$family.ID), res$X_rownames)

  # are row & column names aligned?
  expect_identical(c(res$unpen_colnames, res$X_colnames)[res$ns], res$std_X_colnames)
  expect_identical(as.character(plink$fam$family.ID[res$outcome_idx]), res$std_X_rownames)

  # test 2: is alignment working? -----------------------------
  # shuffle the IDs here, to test alignment
  samp <- sample(1:nrow(predictors))
  predictors_shuff <- predictors[samp, ]
  phen_shuff <- phen[samp, ]

  ### create design ---------------------------------------------------------------
  X <- create_design(
    data_file = penncath_lite,
    rds_dir = temp_dir,
    new_file = "std_penncath_lite",
    add_outcome = phen_shuff,
    outcome_id = "FamID",
    outcome_col = "CAD",
    add_predictor = predictors_shuff,
    predictor_id = 'FamID',
    feature_id = "FID",
    overwrite = TRUE,
    quiet = TRUE
  )

  res <- readRDS(X)

  ### checks  -------------------------------------------------------------------
  # are external files and PLINK fam file aligned on ID var?
  expect_identical(as.character(plink$fam$family.ID), res$X_rownames)

  # are row & column names (IDs) aligned?
  expect_identical(c(res$unpen_colnames, res$X_colnames)[res$ns], res$std_X_colnames)
  expect_identical(as.character(plink$fam$family.ID[res$outcome_idx]), res$std_X_rownames)

  # test 3: more obs. in pheno than geno ------------------------------------------

  penncath_pheno <- read.csv(find_example_data(path = "penncath_clinical.csv"))

  # subset geno data - take just 1000 samples
  imputed_dat <- readRDS(file.path(temp_dir, "process_penncath.rds"))
  imputed_dat$X <- bigmemory::sub.big.matrix(imputed_dat$X, firstRow = 1, lastRow = 1000) |>
    bigmemory::describe()
  imputed_dat$fam <- imputed_dat$fam[1:1000, ]
  saveRDS(imputed_dat, file.path(temp_dir, 'imputed_data_n1000.rds'))

  ### create design ---------------------------------------------------------------
  X <- create_design(
    data_file = file.path(temp_dir, 'imputed_data_n1000.rds'),
    rds_dir = temp_dir,
    new_file = "std_penncath_n1000",
    feature_id = "FID",
    add_outcome = phen,
    outcome_id = "FamID",
    outcome_col = "CAD",
    add_predictor = predictors,
    predictor_id = 'FamID',
    overwrite = TRUE,
    quiet = TRUE
  )

  res <- readRDS(X)

  ### checks  -------------------------------------------------------------------
  # are external files and PLINK fam file aligned on ID var?
  expect_identical(as.character(plink$fam$family.ID[1:1000]), res$X_rownames)

  # are row & column names (IDs) aligned?
  expect_identical(res$std_X_rownames, res$X_rownames[res$outcome_idx])
  expect_identical(res$std_X_colnames, c(res$unpen_colnames, res$X_colnames)[res$ns])

  # test 4: handle the case with no unpenalized predictors ----------------------------

  X <- create_design(
    data_file = penncath_lite,
    feature_id = "FID",
    rds_dir = temp_dir,
    new_file = "std_penncath_lite",
    add_outcome = phen,
    outcome_id = "FamID",
    outcome_col = "CAD",
    overwrite = TRUE,
    quiet = TRUE
  )

  res <- readRDS(X)

  ## checks ----------------------------------------------------------------------
  # are external files and PLINK fam file aligned on ID var?
  expect_identical(as.character(plink$fam$family.ID), res$X_rownames)

  # are row & column names (IDs) aligned?
  expect_identical(res$std_X_rownames, res$X_rownames[res$outcome_idx])
  expect_identical(res$std_X_colnames, c(res$unpen_colnames, res$X_colnames)[res$ns])

  # test 5: confirm model runs without error -------------------------------------

  expect_silent(plmm(X))
})

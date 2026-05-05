local({
  temp_dir <- withr::local_tempdir()
  extdata <- file.path("..", "extdata")
  unzip_example_data(outdir = extdata)

  # test 1: can all 3 pieces be combined correctly?  ----------------------
  # use external pheno and predictor files, no NAs

  ## process plink -------------------------------------------------------------
  penncath_lite <- process_plink(data_dir = extdata,
                                 data_prefix = "penncath_lite",
                                 rds_dir = temp_dir,
                                 rds_prefix = "imputed_penncath_lite",
                                 id_var = "FID",
                                 quiet = FALSE,
                                 overwrite = TRUE)

  ## create design matrix-------------------------------------------------------
  penncath_pheno <- read.csv(file.path(extdata, "penncath_clinical.csv"))

  predictors <- penncath_pheno |>
    dplyr::select(FamID, age, tg) |>
    dplyr::mutate(tg = dplyr::if_else(is.na(tg), mean(tg, na.rm = TRUE), tg),
                  FamID = as.character(FamID))

  phen <- cbind(penncath_pheno$FamID, penncath_pheno$CAD) |>
    as.data.frame() |>
    tidyr::drop_na() |>
    as.matrix()
  colnames(phen) <- c("FamID", "CAD")
  penncath_lite <- process_plink(data_dir = extdata,
                                 data_prefix = "penncath_lite",
                                 rds_dir = temp_dir,
                                 rds_prefix = "process_pencath",
                                 id_var = "FID",
                                 add_phen = phen,
                                 pheno_id = "FamID",
                                 pheno_name = "CAD",
                                 add_predictor_ext = predictors,
                                 quiet = FALSE,
                                 overwrite = TRUE,
                                 logfile = "./test_process_plink")

  phen <- data.frame(FamID = penncath_pheno$FamID, CAD = penncath_pheno$CAD) |>
    dplyr::mutate(FamID = as.character(FamID))
  # note: CAD has no missing values in phen file

  X <- create_design(data_file = penncath_lite,
                     feature_id = "FID",
                     rds_dir = extdata,
                     new_file = "std_penncath_lite",
                     add_outcome = phen,
                     outcome_id = "FamID",
                     outcome_col = "CAD",
                     add_predictor = predictors,
                     predictor_id = 'FamID',
                     overwrite = TRUE,
                     logfile = "design")

  res <- readRDS(X)

  ### checks  -------------------------------------------------------------------
  # are external files and PLINK fam file aligned on ID var?
  plink <- readRDS(file.path(temp_dir, "imputed_penncath_lite.rds"))
  expect_identical(as.character(plink$fam$family.ID), res$X_rownames)

  # are row & column names aligned?
  expect_identical(c(res$unpen_colnames, res$X_colnames)[res$ns], res$std_X_colnames)
  expect_identical(as.character(plink$fam$family.ID[res$outcome_idx]), res$std_X_rownames)

  # test 2: is alignment working? -----------------------------
  # shuffle the IDs here, to test alignment
  penncath_pheno <- penncath_pheno[sample(1:nrow(penncath_pheno)),]

  predictors <- penncath_pheno |>
    dplyr::select(FamID, age, tg) |>
    dplyr::mutate(tg = dplyr::if_else(is.na(tg), mean(tg, na.rm = TRUE), tg),
                  FamID = as.character(FamID))

  phen <- data.frame(FamID = penncath_pheno$FamID, CAD = penncath_pheno$CAD) |>
    dplyr::mutate(FamID = as.character(FamID))

  ### create design ---------------------------------------------------------------
  X <- create_design(data_file = penncath_lite,
                     rds_dir = temp_dir,
                     new_file = "std_penncath_lite",
                     add_outcome = phen,
                     outcome_id = "FamID",
                     outcome_col = "CAD",
                     add_predictor = predictors,
                     predictor_id = 'FamID',
                     feature_id = "FID",
                     overwrite = TRUE,
                     logfile = "design")

  res <- readRDS(X)

  ### checks  -------------------------------------------------------------------
  # are external files and PLINK fam file aligned on ID var?
  plink <- readRDS(file.path(temp_dir, 'imputed_penncath_lite.rds'))
  expect_identical(as.character(plink$fam$family.ID), res$X_rownames)

  # are row & column names (IDs) aligned?
  expect_identical(c(res$unpen_colnames, res$X_colnames)[res$ns], res$std_X_colnames)
  expect_identical(as.character(plink$fam$family.ID[res$outcome_idx]), res$std_X_rownames)

  # test 3: more obs. in pheno than geno ------------------------------------------

  penncath_pheno <- read.csv(file.path(extdata, "penncath_clinical.csv"))

  # subset geno data - take just 1000 samples
  imputed_dat <- readRDS(file.path(temp_dir, 'imputed_penncath_lite.rds'))
  imputed_dat$X <- bigmemory::sub.big.matrix(imputed_dat$X, firstRow = 1, lastRow = 1000) |> bigmemory::describe()
  imputed_dat$fam <- imputed_dat$fam[1:1000,]
  saveRDS(imputed_dat, file.path(temp_dir, 'imputed_data_n1000.rds'))

  predictors <- penncath_pheno |>
    dplyr::select(FamID, age, tg) |>
    dplyr::mutate(tg = dplyr::if_else(is.na(tg), mean(tg, na.rm = T), tg),
                  FamID = as.character(FamID))

  phen <- data.frame(FamID = penncath_pheno$FamID, CAD = penncath_pheno$CAD) |>
    dplyr::mutate(FamID = as.character(FamID))

  ### create design ---------------------------------------------------------------
  X <- create_design(data_file = file.path(temp_dir, 'imputed_data_n1000.rds'),
                     rds_dir = temp_dir,
                     new_file = "std_penncath_n1000",
                     feature_id = "FID",
                     add_outcome = phen,
                     outcome_id = "FamID",
                     outcome_col = "CAD",
                     add_predictor = predictors,
                     predictor_id = 'FamID',
                     overwrite = TRUE,
                     logfile = "design_penncath_n1000")

  res <- readRDS(X)

  ### checks  -------------------------------------------------------------------
  # are external files and PLINK fam file aligned on ID var?
  expect_identical(as.character(plink$fam$family.ID[1:1000]), res$X_rownames)

  # are row & column names (IDs) aligned?
  expect_identical(res$std_X_rownames, res$X_rownames[res$outcome_idx])
  expect_identical(res$std_X_colnames, c(res$unpen_colnames, res$X_colnames)[res$ns])

  # test 4: handle the case with no unpenalized predictors ----------------------------

  penncath_pheno <- read.csv(file.path(extdata, "penncath_clinical.csv"))

  phen <- data.frame(FamID = penncath_pheno$FamID, CAD = penncath_pheno$CAD) |>
    dplyr::mutate(FamID = as.character(FamID))

  X <- create_design(data_file = penncath_lite,
                     feature_id = "FID",
                     rds_dir = temp_dir,
                     new_file = "std_penncath_lite",
                     add_outcome = phen,
                     outcome_id = "FamID",
                     outcome_col = "CAD",
                     overwrite = TRUE,
                     logfile = "design")

  res <- readRDS(X)

  ## checks ----------------------------------------------------------------------
  # are external files and PLINK fam file aligned on ID var?
  expect_identical(as.character(plink$fam$family.ID), res$X_rownames)

  # are row & column names (IDs) aligned?
  expect_identical(res$std_X_rownames, res$X_rownames[res$outcome_idx])
  expect_identical(res$std_X_colnames, c(res$unpen_colnames,res$X_colnames)[res$ns])

  # test 5: confirm model runs without error -------------------------------------

  expect_silent(plmm(X))
})
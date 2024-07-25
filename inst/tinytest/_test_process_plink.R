if (interactive()){
  # tests with penncath_lite data (ships with package)
  # setup --------------------------------------
  devtools::load_all(".")
  library(dplyr)
  library(tidyr)
  library(bigsnpr)
  unzip_example_data(outdir = 'inst/extdata')

  # test 1: can all 3 pieces be combined correctly?  ----------------------
  # use external pheno and predictor files, no NAs

  ## process plink -------------------------------------------------------------
  penncath_lite <- process_plink(data_dir = "inst/extdata",
                                 data_prefix = "penncath_lite",
                                 rds_prefix = "imputed_penncath_lite",
                                 id_var = "FID",
                                 quiet = FALSE,
                                 overwrite = TRUE)

  imputed_dat <- readRDS(penncath_lite)

  # check out imputed data
  # str(imputed_dat)
  # any(is.na(imputed_dat$genotypes[,]))

  ## create design matrix-------------------------------------------------------
  penncath_pheno <- read.csv("inst/extdata/penncath_clinical.csv")
  str(penncath_pheno)

  predictors <- penncath_pheno |>
    dplyr::select(FamID, age, tg) |>
    dplyr::mutate(tg = dplyr::if_else(is.na(tg), mean(tg, na.rm = T), tg),
                  FamID = as.character(FamID))

  phen <- data.frame(FamID = penncath_pheno$FamID, CAD = penncath_pheno$CAD) |>
    mutate(FamID = as.character(FamID))
  # note: CAD has no missing values in phen file

  X <- create_design(dat_file = penncath_lite,
                     feature_id = "FID",
                     rds_dir = "inst/extdata",
                     new_file = "std_penncath_lite",
                     add_outcome = phen,
                     outcome_id = "FamID",
                     outcome_col = "CAD",
                     add_predictor = predictors,
                     predictor_id = 'FamID',
                     overwrite = TRUE,
                     logfile = "design")

  res <- readRDS(X)
  str(res)

}

### checks  -------------------------------------------------------------------
# are external files and PLINK fam file aligned on ID var?
plink <- readRDS('inst/extdata/imputed_penncath_lite.rds')
tinytest::expect_identical(as.character(plink$fam$family.ID), res$X_rownames)

# does final fam[,6] have the expected outcome?
table(plink$fam$affection, penncath_pheno$CAD)

# are row & column names aligned?
str(res$X_colnames[res$ns]); str(res$X_rownames)
str(res$std_X_colnames); str(res$std_X_rownames)
tinytest::expect_identical(c(res$non_gen_colnames,res$X_colnames)[res$ns], res$std_X_colnames)
tinytest::expect_identical(plink$rownames[res$outcome_idx], res$std_X_rownames)

# are .bk files 'cleaned up' and labeled correctly?
list.files("inst/extdata", pattern = "*.bk")

# clear example
rm(X); rm(res); rm(predictors); rm(phen)
list.files("inst/extdata", pattern = "std_penncath_lite.*", full.names = T) |> file.remove()

# test 2: is alignment working? -----------------------------
# shuffle the IDs here, to test alignment
penncath_pheno <- penncath_pheno[sample(1:nrow(penncath_pheno)),]

predictors <- penncath_pheno |>
  dplyr::select(FamID, age, tg) |>
  dplyr::mutate(tg = dplyr::if_else(is.na(tg), mean(tg, na.rm = T), tg),
                FamID = as.character(FamID))

phen <- data.frame(FamID = penncath_pheno$FamID, CAD = penncath_pheno$CAD) |>
  mutate(FamID = as.character(FamID))

### create design ---------------------------------------------------------------
X <- create_design(dat_file = penncath_lite,
                   rds_dir = "inst/extdata",
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
str(res)

### checks  -------------------------------------------------------------------
# are external files and PLINK fam file aligned on ID var?
plink <- readRDS('inst/extdata/imputed_penncath_lite.rds')
tinytest::expect_identical(as.character(plink$fam$family.ID), res$X_rownames)

# does final fam[,6] have the expected outcome? (same number of affected?)
table(plink$fam[,6])
table(penncath_pheno$CAD)

# are row & column names (IDs) aligned?
str(res$X_colnames[res$ns]); str(res$X_rownames)
str(res$std_X_colnames); str(res$std_X_rownames)
tinytest::expect_identical(c(res$non_gen_colnames,res$X_colnames)[res$ns], res$std_X_colnames)
tinytest::expect_identical(plink$rownames[res$outcome_idx], res$std_X_rownames)

# are .bk files 'cleaned up' and labeled correctly?
list.files("inst/extdata", pattern = "*.bk")

# clear example
rm(X); rm(res); rm(predictors); rm(phen)
list.files("inst/extdata", pattern = "std_penncath_lite.*", full.names = T) |> file.remove()

# test 3: more obs. in pheno than geno ------------------------------------------
penncath_pheno <- read.csv("inst/extdata/penncath_clinical.csv")

# subset geno data - take just 1000 samples
imputed_dat <- readRDS('inst/extdata/imputed_penncath_lite.rds')
imputed_dat$X <- bigmemory::sub.big.matrix(imputed_dat$X, firstRow = 1, lastRow = 1000) |> bigmemory::describe()
imputed_dat$fam <- imputed_dat$fam[1:1000,]
saveRDS(imputed_dat,'inst/extdata/imputed_data_n1000.rds')

predictors <- penncath_pheno |>
  dplyr::select(FamID, age, tg) |>
  dplyr::mutate(tg = dplyr::if_else(is.na(tg), mean(tg, na.rm = T), tg),
                FamID = as.character(FamID))


phen <- data.frame(FamID = penncath_pheno$FamID, CAD = penncath_pheno$CAD) |>
  mutate(FamID = as.character(FamID))

### create design ---------------------------------------------------------------
X <- create_design(dat = 'inst/extdata/imputed_data_n1000.rds',
                   rds_dir = 'inst/extdata',
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
str(res)

### checks  -------------------------------------------------------------------
# does final fam[,6] have the expected outcome?
tinytest::expect_identical(as.character(plink$fam$family.ID), res$X_rownames)

# are row & column names (IDs) aligned?
str(res$X_colnames); str(res$X_rownames)
str(res$std_X_colnames); str(res$std_X_rownames)
tinytest::expect_identical(res$std_X_rownames, res$X_rownames[res$outcome_idx])
tinytest::expect_identical(res$std_X_colnames, c(res$non_gen_colnames,res$X_colnames)[res$ns])

# are .bk files 'cleaned up' and labeled correctly?
list.files("inst/extdata", pattern = "*.bk")

# clear example
rm(X); rm(res); rm(predictors); rm(phen)

# test 4: handle the case with no predictors ----------------------------

penncath_pheno <- read.csv("inst/extdata/penncath_clinical.csv")

phen <- data.frame(FamID = penncath_pheno$FamID, CAD = penncath_pheno$CAD) |>
  mutate(FamID = as.character(FamID))

X <- create_design(dat_file = penncath_lite,
                   feature_id = "FID",
                   rds_dir = "inst/extdata",
                   new_file = "std_penncath_lite",
                   add_outcome = phen,
                   outcome_id = "FamID",
                   outcome_col = "CAD",
                   overwrite = TRUE,
                   logfile = "design")

res <- readRDS(X)
str(res)

## chekcs ----------------------------------------------------------------------
# does final fam[,6] have the expected outcome?
tinytest::expect_identical(as.character(plink$fam$family.ID), res$X_rownames)

# are row & column names (IDs) aligned?
str(res$X_colnames); str(res$X_rownames)
str(res$std_X_colnames); str(res$std_X_rownames)
tinytest::expect_identical(res$std_X_rownames, res$X_rownames[res$outcome_idx])
tinytest::expect_identical(res$std_X_colnames, c(res$non_gen_colnames,res$X_colnames)[res$ns])

# are .bk files 'cleaned up' and labeled correctly?
list.files("inst/extdata", pattern = "*.bk")

# clear example
rm(X); rm(res); rm(predictors); rm(phen)


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
                                 prefix = "penncath_lite",
                                 id_var = "FID",
                                 quiet = FALSE,
                                 overwrite = TRUE)

  imputed_dat <- readRDS(penncath_lite)
  str(imputed_dat)
  any(is.na(imputed_dat$genotypes[,]))

  ## create design matrix-------------------------------------------------------
  penncath_pheno <- read.csv("inst/extdata/penncath_clinical.csv")
  str(penncath_pheno)

  predictors <- penncath_pheno |>
    dplyr::select(FamID, age, tg) |>
    dplyr::mutate(tg = dplyr::if_else(is.na(tg), mean(tg, na.rm = T), tg)) |>
    tibble::column_to_rownames('FamID') |>
    as.matrix()
  colnames(predictors) <- c("age", "tg")

  phen <- cbind(penncath_pheno$FamID, penncath_pheno$CAD) |>
    as.data.frame() |>
    as.matrix()
  colnames(phen) <- c("FamID", "CAD") # CAD has no missing values in phen file


  X <- create_design(dat = penncath_lite,
                     rds_dir = "inst/extdata",
                     new_file = "std_penncath_lite",
                     is_bigsnp = TRUE,
                     add_phen = phen,
                     pheno_id = "FamID",
                     pheno_name = "CAD",
                     add_predictor_ext = predictors,
                     id_var = "FID",
                     overwrite = TRUE,
                     outfile = "design")

  res <- readRDS(X)
  str(res)

}

### checks  -------------------------------------------------------------------
# are external files and PLINK fam file aligned on ID var?
tinytest::expect_identical(as.character(res$fam$family.ID), res$X_rownames)

# does final fam[,6] have the expected outcome?
table(res$fam$affection)
table(penncath_pheno$CAD)

# are row & column names aligned?
str(res$X_colnames[res$ns]); str(res$X_rownames)
str(res$std_X_colnames); str(res$std_X_rownames)
tinytest::expect_identical(res$X_colnames[res$ns], res$std_X_colnames)
tinytest::expect_identical(res$X_rownames[res$complete_phen], res$std_X_rownames)

# are .bk files 'cleaned up' and labeled correctly?
list.files("inst/extdata", pattern = "*.bk")

# clear example
rm(X); rm(res); rm(predictors); rm(phen)
list.files("inst/extdata", pattern = "std_penncath_lite.*", full.names = T) |> file.remove()

# test 2: is alignment working? -----------------------------
# shuffle the IDs here, to test alignment
penncath_pheno <- penncath_pheno[sample(1:nrow(penncath_pheno)),]

predictors <- penncath_pheno |>
  dplyr::select(age, tg) |>
  dplyr::mutate(tg = dplyr::if_else(is.na(tg), mean(tg, na.rm = T), tg)) |>
  as.matrix()
colnames(predictors) <- c("age", "tg")
rownames(predictors) <- penncath_pheno$FamID

phen <- cbind(penncath_pheno$FamID, penncath_pheno$CAD) |>
  as.data.frame() |>
  as.matrix()
colnames(phen) <- c("FamID", "CAD") # CAD has no missing values in phen file

### create design ---------------------------------------------------------------
X <- create_design(dat = penncath_lite,
                   rds_dir = "inst/extdata",
                   new_file = "std_penncath_lite",
                   is_bigsnp = TRUE,
                   add_phen = phen,
                   pheno_id = "FamID",
                   pheno_name = "CAD",
                   add_predictor_ext = predictors,
                   id_var = "FID",
                   overwrite = TRUE,
                   outfile = NULL)

res <- readRDS(X)
str(res)

### checks  -------------------------------------------------------------------
# are external files and PLINK fam file aligned on ID var?
tinytest::expect_identical(as.character(res$fam$family.ID), res$X_rownames)

# does final fam[,6] have the expected outcome?
table(res$fam[,6])
table(penncath_pheno$CAD)

# are row & column names (IDs) aligned?
str(res$X_colnames); str(res$X_rownames)
str(res$std_X_colnames); str(res$std_X_rownames)
tinytest::expect_identical(res$std_X_rownames, res$X_rownames[res$complete_phen])
tinytest::expect_identical(res$std_X_colnames, res$X_colnames[res$ns])

# are .bk files 'cleaned up' and labeled correctly?
list.files("inst/extdata", pattern = "*.bk")

# clear example
rm(X); rm(res); rm(predictors); rm(phen)
list.files("inst/extdata", pattern = "std_penncath_lite.*", full.names = T) |> file.remove()

# test 3: what if there are missing outcomes? ----------------------------------

# shuffle the IDs here, to test alignment
penncath_pheno <- penncath_pheno[sample(1:nrow(penncath_pheno)),]

predictors <- penncath_pheno |>
  dplyr::select(age, tg) |>
  dplyr::mutate(tg = dplyr::if_else(is.na(tg), mean(tg, na.rm = T), tg)) |>
  as.matrix()
colnames(predictors) <- c("age", "tg")
rownames(predictors) <- penncath_pheno$FamID

# use 'hdl' for outcome - missing in some observations
summary(penncath_pheno$hdl)
phen <- cbind(penncath_pheno$FamID, penncath_pheno$hdl) |>
  as.data.frame() |>
  as.matrix()
colnames(phen) <- c("FamID", "hdl")

### create design --------------------------------------------------------------
X <- create_design(dat = penncath_lite,
                   rds_dir = "inst/extdata",
                   new_file = "std_penncath_lite",
                   is_bigsnp = TRUE,
                   add_phen = phen,
                   pheno_id = "FamID",
                   pheno_name = "hdl",
                   na_phenotype_vals = c(NA_integer_),
                   add_predictor_ext = predictors,
                   id_var = "FID",
                   overwrite = TRUE,
                   outfile = NULL)

res <- readRDS(X)
str(res)

### checks  -------------------------------------------------------------------
# are external files and PLINK fam file aligned on ID var?
tinytest::expect_identical(as.character(res$fam$family.ID), res$X_rownames)

# does final fam[,6] have the expected outcome?
summary(res$fam[,6])
summary(penncath_pheno$hdl)

# are row & column names (IDs) aligned?
str(res$X_colnames); str(res$X_rownames)
str(res$std_X_colnames); str(res$std_X_rownames)
tinytest::expect_identical(res$std_X_rownames, res$X_rownames[res$complete_phen])
tinytest::expect_identical(res$std_X_colnames, res$X_colnames[res$ns])

# are .bk files 'cleaned up' and labeled correctly?
list.files("inst/extdata", pattern = "*.bk")

# clear example
rm(X); rm(res); rm(predictors); rm(phen)
list.files("inst/extdata", pattern = "std_penncath_lite.*", full.names = T) |> file.remove()

# test 4: more obs. in pheno than geno ------------------------------------------
penncath_pheno <- read.csv("inst/extdata/penncath_clinical.csv")

# subset geno data
n1000 <- sample(1:nrow(imputed_dat$genotypes), 1000) |> sort()
geno_n1000 <- snp_subset(x = imputed_dat,
                         ind.row = n1000)

bigsnp_n1000 <- snp_attach(geno_n1000)
bigsnp_n1000$colnames <- imputed_dat$colnames
bigsnp_n1000$rownames <- imputed_dat$rownames[n1000]
bigsnp_n1000$n <- 1000
bigsnp_n1000$p <- ncol(bigsnp_n1000$genotypes)
snp_save(bigsnp_n1000)

predictors <- penncath_pheno |>
  dplyr::select(FamID, age, tg) |>
  dplyr::mutate(tg = dplyr::if_else(is.na(tg), mean(tg, na.rm = T), tg)) |>
  tibble::column_to_rownames('FamID') |>
  as.matrix()
colnames(predictors) <- c("age", "tg")

phen <- cbind(penncath_pheno$FamID, penncath_pheno$CAD) |>
  as.data.frame() |>
  as.matrix()
colnames(phen) <- c("FamID", "CAD") # CAD has no missing values in phen file

### create design ---------------------------------------------------------------
X <- create_design(dat = geno_n1000,
                   rds_dir = 'inst/extdata',
                   new_file = "std_penncath_n1000",
                   is_bigsnp = T,
                   add_phen = phen,
                   pheno_id = "FamID",
                   pheno_name = "CAD",
                   add_predictor_ext = predictors,
                   id_var = "FID",
                   overwrite = TRUE,
                   outfile = "design_penncath_n1000")

res <- readRDS(X)
str(res)

### checks  -------------------------------------------------------------------
# are external files and PLINK fam file aligned on ID var?
tinytest::expect_identical(as.character(res$fam$family.ID), res$X_rownames)

# does final fam[,6] have the expected outcome?
table(res$fam[,6], penncath_pheno$CAD[n1000], useNA = 'ifany')

# are row & column names (IDs) aligned?
str(res$X_colnames); str(res$X_rownames)
str(res$std_X_colnames); str(res$std_X_rownames)
tinytest::expect_identical(res$std_X_rownames, res$X_rownames[res$complete_phen])
tinytest::expect_identical(res$std_X_colnames, res$X_colnames[res$ns])

# are .bk files 'cleaned up' and labeled correctly?
list.files("inst/extdata", pattern = "*.bk")

# clear example
rm(X); rm(res); rm(predictors); rm(phen)


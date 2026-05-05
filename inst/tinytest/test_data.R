# Test 1 - Example datasets can be unzipped ------------------------------------
tmp <- tempdir()
expect_message(unzip_example_data(tmp), pattern = "Unzipped")

# Test 2 - Overwriting filebacked data succeeds ----------------------------

expect_silent(local({
  # using PennCath data that ships with the package
  temp_dir <- withr::local_tempdir()
  unzip_example_data(outdir = temp_dir)

  plink_data <- process_plink(data_dir = temp_dir,
                              data_prefix = "penncath_lite",
                              rds_dir = temp_dir,
                              rds_prefix = "imputed_penncath_lite",
                              impute_method = "mode",
                              parallel = FALSE)

  # get outcome data
  penncath_pheno <- read.csv(find_example_data(path = 'penncath_clinical.csv'))

  phen <- data.frame(FamID = as.character(penncath_pheno$FamID),
                     CAD = penncath_pheno$CAD)

  # prepare a data.frame of the predictors for which we want to adjust:
  other_predictors <- penncath_pheno[,c('FamID', 'sex', 'age')]
  other_predictors$FamID <- as.character(other_predictors$FamID)

  pen_design <- create_design(data_file = plink_data,
                              feature_id = "FID",
                              rds_dir = temp_dir,
                              new_file = "std_penncath_lite",
                              add_outcome = phen,
                              outcome_id = "FamID",
                              outcome_col = "CAD",
                              add_predictor = other_predictors,
                              predictor_id = 'FamID',
                              logfile = "design")

  # Run again
  plink_data <- process_plink(data_dir = temp_dir,
                              data_prefix = "penncath_lite",
                              rds_dir = temp_dir,
                              rds_prefix = "imputed_penncath_lite",
                              impute_method = "mode",
                              parallel = FALSE,
                              overwrite = TRUE)

  pen_design <- create_design(data_file = plink_data,
                              feature_id = "FID",
                              rds_dir = temp_dir,
                              new_file = "std_penncath_lite",
                              add_outcome = phen,
                              outcome_id = "FamID",
                              outcome_col = "CAD",
                              add_predictor = other_predictors,
                              predictor_id = 'FamID',
                              logfile = "design",
                              overwrite = TRUE)
}))

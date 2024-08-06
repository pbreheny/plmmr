# TKP
# August 2024
# Goal: analyze the penncath_mid data -- test how plmmr methods scale up

# process the data -----------------------------------------------------
pen_mid <- process_plink(data_dir = "inst/extdata/",
                         data_prefix = "penncath_mid",
                         rds_dir = "inst/extdata/",
                         rds_prefix = "imputed_penncath_mid",
                         id_var = "FID",
                         logfile = "imputed_penncath_mid")

# external data ---------------------------------------------------------
pen_phen <- read.csv("inst/extdata/penncath_clinical.csv") |>
  dplyr::mutate(FamID = as.character(FamID)) # ID variable must be a character
extdata <- pen_phen |>
  dplyr::select(FamID, sex, age)

# create a design ------------------------------------------------------
pen_design <- create_design(dat_file = pen_mid,
                            feature_id = 'FID',
                            rds_dir = "inst/extdata/",
                            new_file = "std_penncath_mid",
                            add_outcome = pen_phen,
                            outcome_id = "FamID",
                            outcome_col = "CAD",
                            add_predictor = extdata,
                            predictor_id = "FamID",
                            overwrite = TRUE,
                            logfile = 'pen_design')

# fit a model -----------------------------------------------------------
fit <- plmm(design = pen_design, trace = TRUE)


# CV -------------------------------------------------------------------




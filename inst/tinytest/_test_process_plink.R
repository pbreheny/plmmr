if (interactive()){
# test1: penncath_lite data (ships with package) ----------------------------
## setup --------------------------------------
devtools::load_all(".")
library(dplyr)
library(tidyr)
library(bigsnpr)

## read in phenotype data --------------------------
pen_clinic <- read.csv(paste0(find_example_data(parent = TRUE), "/penncath_clinical.csv"))
extdata <- pen_clinic[,3:4]
rownames(extdata) <- pen_clinic$FamID # This is important!

## process PLINK data ------------------------------
temp_dir <- paste0(tempdir(), sample(LETTERS, 1))
process_plink(data_dir = find_example_data(parent = TRUE),
              rds_dir = temp_dir, # using a temporary directory
              prefix = "penncath_lite",
              id_var = "FID", # this is KEY!
              outfile = "process_penncath",
              impute_method = "mode",
              add_predictor_fam = "sex",
              add_phen = pen_clinic,
              pheno_id = "FamID",
              pheno_name = "CAD")

# check this out:
pen2 <- readRDS(paste0(temp_dir, "/std_penncath_lite.rds"))
pen2$std_X_colnames |> head() # std_X includes our non-genomic covariates

# note: dim(pen2$std_X) is different from the vignette; this is because
# here we are testing the functionality of passing in an external phenotype
# file. The vignette uses the 6th column of the fam file as the outcome,
# and it recognizes the -9 values as missing

# test 2: use external pheno and predictor files ----------------------------
unzip_example_data(outdir = 'inst/extdata')

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
  tidyr::drop_na() |>
  as.matrix()
colnames(phen) <- c("FamID", "CAD")

system.time(

  penncath_lite <- process_plink(data_dir = "inst/extdata/",
                                 prefix = "penncath_lite",
                                 id_var = "FID",
                                 add_phen = phen,
                                 pheno_id = "FamID",
                                 pheno_name = "CAD",
                                 add_predictor_ext = predictors,
                                 quiet = FALSE,
                                 overwrite = TRUE,
                                 outfile = "./test_process_plink")


)

}


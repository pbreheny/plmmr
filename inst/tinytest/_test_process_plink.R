# test1: penncath_lite data (ships with package) ----------------------------
## setup --------------------------------------
devtools::load_all(".")
library(dplyr)
library(tidyr)
library(bigsnpr)

## read in phenotype data --------------------------
pen_clinic <- read.csv(paste0(get_example_data(parent = TRUE), "/penncath_clinical.csv"))
extdata <- pen_clinic[,3:4]
rownames(extdata) <- pen_clinic$FamID # This is important! 

## process PLINK data ------------------------------
temp_dir <- paste0(tempdir(), sample(LETTERS, 1))
process_plink(data_dir = get_example_data(parent = TRUE),
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

# test 2: chr10_lite data (from TKP's OFC data) ---------------------
# Note: filepaths here will be machine-dependent
thesis_path <- "~/thesis/" # change this for current machine; use TKP's thesis dir
pheno <- readRDS(paste0(thesis_path, "plmm_tests/pheno_complete.rds"))
str(pheno)

# create a binary outcome
pheno <- pheno |>
  dplyr::mutate(ofc_bin = dplyr::case_when(
    CleftType == 0 ~ 0, # 0 = unaffected
    CleftType > 1 ~ 1, # 1 = indeterminate phenotype , > 1 = cleft
    .default = NA_integer_
  ),
  # make ID a character (for compatibility with genotype data)
  IID = as.character(IID)) |>
  # remove missing phenotype
  dplyr::filter(!is.na(ofc_bin))

# process plink files 
process_plink(data_dir = paste0(thesis_path,"plmm_tests/"),
              prefix = "chr10_lite",
              add_predictor_fam = "sex",
              add_phen = pheno,
              pheno_id = "IID",
              pheno_name = "ofc_bin",
              gz = FALSE,
              outfile = "./filebacked_ch10_lite")

# checks 
chr10_lite_rds <- readRDS(paste0(thesis_path,"plmm_tests/std_chr10_lite.rds"))
str(chr10_lite_rds)
chr10_lite_rds$std_X[1:10, 1:8]
table(chr10_lite_rds$fam$affection, useNA = "ifany")
std_X <- chr10_lite_rds$std_X[,]
colMeans(std_X) |> summary() # columns have mean zero...
apply(std_X, 2, var) |> summary() # ... & variance 1

# run checks 
checked_data <- plmmr:::plmm_checks(X = paste0(thesis_path,
                                               "plmm_tests/std_chr10_lite"),
                                    returnX = FALSE)
str(checked_data)

# try fitting a model
# after 1st test run, I save the 'prep' for future testing - this is about 1.8 GB
chr10_lite_prep <- readRDS(paste0(thesis_path, "plmm_tests/std_chr10_lite_prep.rds"))


chr10_lite_fit <- plmm(X = paste0(thesis_path, 
                             "plmm_tests/std_chr10_lite"),
                  penalty.factor = c(0, rep(1, ncol(checked_data$std_X) - 1)),
                  # give 'sex' a penalty factor of 0 -- we always want this in the model.
                  # The 'intercept' will get a 0 penalty factor by default.
                  returnX = FALSE, # run filebacked,
                  K = list(U = chr10_lite_prep$U, s = chr10_lite_prep$s), 
                  save_rds = paste0(thesis_path, 
                                    "plmm_tests/std_chr10_lite_fit.rds"),
                  trace = T)







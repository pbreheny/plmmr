# Objective of this script: run plmmr::process_plink() from command line

# directory setup -----------------------------------------
# set directory to 'cleft' LSS drive
if (Sys.info()["user"] == "tabpeter") {
  # case 1: TKP is working on her desktop machine
  if (grepl("DENT", Sys.info()["nodename"])) cleft_drive  <- "/mnt/cleft"
  # case 2: TKP is working on Argon
  if (grepl("argon", Sys.info()["nodename"])) cleft_drive  <- "/Shared/cleft"

} else if (Sys.info()["user"] == "pbreheny") {
  cleft_drive <- glue::glue('~/lss/cleft') # TODO: check to make sure this is right
}

rds_file <- "plmm_analysis_ofc/phs000774"
bfile <- "stage1_filtered"

# read in external data --------------------------
all_pheno <- readRDS(file.path(cleft_drive, "data", "phs000774", "pheno", "whole",
                               "ofc_whole_pheno.rds"))

pheno <- all_pheno[,c('IID', 'OFC')] |> as.matrix()
colnames(pheno) <- c('IID', 'OFC')

# str(pheno) # see what's here
sex_and_site_matrix <- read.delim(file.path(cleft_drive, "data", "phs000774",
                                            "pheno", "whole",
                                            "sex_and_site_matrix.txt"))

# str(sex_and_site_matrix) # IIDs are rownames here

# process PLINK data ------------------------------
# Note: the call below to process_plink() created 'stage1_filtered.rds' & 'stage1_filtered.bk'
# imputed_dat <- plmmr::process_plink(data_dir = file.path(cleft_drive, "data", "phs000774",
#                                                          "qc", "whole"),
#                                     rds_dir = file.path(cleft_drive, rds_file),
#                                     prefix = bfile,
#                                     outfile = paste0(file.path(cleft_drive, rds_file),"/" , bfile))
#
# the_imputed_data <- readRDS(imputed_dat)

imputed_dat <- '/mnt/cleft/plmm_analysis_ofc/phs000774/stage1_filtered.rds'

# checks
# str(the_imputed_data)
# the_imputed_data$genotypes[,1:1500] -> foo
# any(is.na(the_imputed_data))

# create design ------------------------------------
design <- plmmr::create_design(dat = imputed_dat,
                               rds_dir = file.path(cleft_drive, rds_file),
                               new_file = "yet_another_test_design",
                               is_bigsnp = TRUE,
                               na_phenotype_vals = NA_integer_,
                               add_predictor_ext = sex_and_site_matrix,
                               add_phen = pheno, # add in phenotype info from external file
                               pheno_id = "IID",
                               pheno_name = "OFC", # name of phenotype column
                               outfile = "yet_another_test_design",
                               quiet = FALSE)

# 'browser' will be called before standardize_bigsnp()

# for argon ------------------------------------------------------------------
# to run on argon:
# qlogin -q BIOSTAT -t 1-5 -b y Rscript plmm_analysis_ofc/scripts/cleft_process_plink.R

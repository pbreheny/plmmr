## code to prepare `admix_raw` dataset goes here
# load libraries
library(dplyr)
library(plmmr)
# read in data for examples
admix_raw <- read.delim("https://s3.amazonaws.com/pbreheny-data-sets/admixture.txt")
# str(admix_raw) # includes 197 obs., 100 SNPs, and racial category


# Simulate Y, an outcome representing a continuous phenotype
# Note: Y will depend on race as well as some SNP information, in order to illustrate
#   population structure (as is common in GWAS studies)

# make the ancestry variable into a numeric value
table(admix_raw$Race) # groups are approx. equal size
admix_raw$Race <- admix_raw$Race |> as.factor() |> as.numeric() - 1
# check:
table(admix_raw$Race) # 0 = African, 1 = African American, 2 = European, 3 = Japanese

# determine which 2 SNPs will be significant
set.seed(522)
signif_snps <- sample(x = 1:100, size = 2) # SNPs 17 & 85 are chosen

# determine the impact of each significant SNP on the outcome
true_beta <- rep(0, ncol(admix_raw))
true_beta[1] <- 2 # make ancestry have a relatively strong relationship
true_beta[signif_snps] <- c(1, -1) # make two SNPs have an association with the outcome

# simulate some random noise
set.seed(522)
noise <- rnorm(n = nrow(admix_raw))

# create an outcome variable
true_X <- as.matrix(admix_raw)
y <- true_X%*%true_beta + noise

# create objects to export to user level
X <- admix_raw |>
  dplyr::select(-c(Race)) |>
  as.matrix()

# create a list with the data needed for analyses
# note: 'ancestry' is a more appropriate word choice than 'Race'
admix <- list(X = X,
              y = y,
              ancestry = admix_raw$Race
)

usethis::use_data(admix, overwrite = TRUE)

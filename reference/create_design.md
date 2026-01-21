# a function to create a design for PLMM modeling

a function to create a design for PLMM modeling

## Usage

``` r
create_design(data_file = NULL, rds_dir = NULL, X = NULL, y = NULL, ...)
```

## Arguments

- data_file:

  For **filebacked data** (data from
  [`process_plink()`](https://pbreheny.github.io/plmmr/reference/process_plink.md)
  or
  [`process_delim()`](https://pbreheny.github.io/plmmr/reference/process_delim.md)),
  this is the filepath to the processed data. Defaults to NULL (this
  argument does not apply for in-memory data).

- rds_dir:

  For **filebacked data**, this is the filepath to the directory/folder
  where you want the design to be saved. **Note**: do not include/append
  the name you want for the to-be-created file – the name is the
  argument `new_file`, passed to
  [`create_design_filebacked()`](https://pbreheny.github.io/plmmr/reference/create_design_filebacked.md).
  Defaults to NULL (this argument does not apply for in-memory data).

- X:

  For **in-memory data (data in a matrix or data frame)**, this is the
  design matrix. Defaults to NULL (this argument does not apply for
  filebacked data).

- y:

  For **in-memory data**, this is the numeric vector representing the
  outcome. Defaults to NULL (this argument does not apply for filebacked
  data). **Note**: it is the responsibility of the user to ensure that
  the rows in X and the corresponding elements of y have the same row
  order, i.e., observations must be in the same order in both the design
  matrix and in the outcome vector.

- ...:

  Additional arguments to pass to
  [`create_design_filebacked()`](https://pbreheny.github.io/plmmr/reference/create_design_filebacked.md)
  or
  [`create_design_in_memory()`](https://pbreheny.github.io/plmmr/reference/create_design_in_memory.md).
  See the documentation for those helper functions for details.

## Value

A filepath to an object of class `plmm_design`, which is a named list
with the design matrix, outcome, penalty factor vector, and other
details needed for fitting a model. This list is stored as an .rds file
for filebacked data, so in the filebacked case a string with the path to
that file is returned. For in-memory data, the list itself is returned.

## Details

This function is a wrapper for the other `create_design...()` inner
functions; all arguments included here are passed along to the
`create_design...()` inner function that matches the type of the data
being supplied. Note which arguments are optional and which ones are
not.

Additional arguments for **all filebacked** data:

- **new_file** User-specified filename (*without .bk/.rds extension*)
  for the to-be-created .rds/.bk files. Must be different from any
  existing .rds/.bk files in the same folder.

- **feature_id** Optional: A string specifying the column in the data X
  (the feature data) with the row IDs (e.g., identifiers for each
  row/sample/participant/, etc.). No duplicates allowed. - for PLINK
  data: a string specifying an ID column of the PLINK `.fam` file.
  Options are "IID" (default) and "FID" - for all other filebacked data:
  a character vector of unique identifiers (IDs) for each row of the
  feature data (i.e., the data processed with
  [`process_delim()`](https://pbreheny.github.io/plmmr/reference/process_delim.md)) -
  if left NULL (default), X is assumed to have the same row-order as
  add_outcome. **Note**: if this assumption is made in error,
  calculations downstream will be incorrect. Pay close attention here.

- **add_outcome** A data frame or matrix with two columns: and ID column
  and a column with the outcome value (to be used as 'y' in the final
  design). IDs must be characters, outcome must be numeric.

- **outcome_id** A string specifying the name of the ID column in
  'add_outcome'

- **outcome_col** A string specifying the name of the phenotype column
  in 'add_outcome'

- **na_outcome_vals** Optional: a vector of numeric values used to code
  NA values in the outcome. Defaults to `c(-9, NA_integer)` (the -9
  matches PLINK conventions).

- **overwrite** Optional: logical - should existing .rds files be
  overwritten? Defaults to FALSE.

- **logfile** Optional: name of the '.log' file to be written –
  **Note:** do not append a `.log` to the filename; this is done
  automatically.

- **quiet** Optional: logical - should messages to be printed to the
  console be silenced? Defaults to FALSE

Additional arguments specific to **PLINK** data:

- **add_predictor** Optional (for PLINK data only): a matrix or data
  frame to be used for adding additional **unpenalized**
  covariates/predictors/features from an external file (i.e., not a
  PLINK file). This matrix must have one column that is an ID column;
  all other columns aside the ID will be used as covariates in the
  design matrix. Columns must be named.

- **predictor_id** Optional (for PLINK data only): A string specifying
  the name of the column in 'add_predictor' with sample IDs. Required if
  'add_predictor' is supplied. The names will be used to subset and
  align this external covariate with the supplied PLINK data.

Additional arguments specific to **delimited file** data:

- **unpen** Optional: an character vector with the names of columns to
  mark as unpenalized (i.e., these features would always be included in
  a model). **Note**: if you choose to use this option, your delimited
  file **must** have column names.

Additional arguments for **in-memory** data:

- **unpen** Optional: an character vector with the names of columns to
  mark as unpenalized (i.e., these features would always be included in
  a model). **Note**: if you choose to use this option, X must have
  column names.

## Examples

``` r
## Example 1: matrix data in-memory ##
admix_design <- create_design(X = admix$X, y = admix$y, unpen = "Snp1")

## Example 2: delimited data ##
# process delimited data
temp_dir <- tempdir()
colon_dat <- process_delim(data_file = "colon2.txt",
 data_dir = find_example_data(parent = TRUE), overwrite = TRUE,
 rds_dir = temp_dir, rds_prefix = "processed_colon2", sep = "\t", header = TRUE)
#> There are 62 observations and 2001 features in the specified data files.
#> At this time, plmmr::process_delim() does not not handle missing values in delimited data.
#>       Please make sure you have addressed missingness before you proceed.
#> 
#> process_plink() completed 
#> Processed files now saved as /tmp/RtmpZ6kWm2/processed_colon2.rds

# prepare outcome data
colon_outcome <- read.delim(find_example_data(path = "colon2_outcome.txt"))

# create a design
colon_design <- create_design(data_file = colon_dat, rds_dir = temp_dir, new_file = "std_colon2",
add_outcome = colon_outcome, outcome_id = "ID", outcome_col = "y", unpen = "sex",
overwrite = TRUE, logfile = "test.log")
#> No feature_id supplied; will assume data X are in same row-order as add_outcome.
#> There are 0 constant features in the data
#> Subsetting data to exclude constant features (e.g., monomorphic SNPs)
#> Column-standardizing the design matrix...
#> Standardization completed at 2026-01-21 16:54:33
#> Done with standardization. File formatting in progress

# look at the results
colon_rds <- readRDS(colon_design)
str(colon_rds)
#> List of 18
#>  $ X_colnames    : chr [1:2001] "sex" "Hsa.3004" "Hsa.13491" "Hsa.13491.1" ...
#>  $ X_rownames    : chr [1:62] "row1" "row2" "row3" "row4" ...
#>  $ n             : num 62
#>  $ p             : num 2001
#>  $ is_plink      : logi FALSE
#>  $ outcome_idx   : int [1:62] 1 2 3 4 5 6 7 8 9 10 ...
#>  $ y             : int [1:62] 1 0 1 0 1 0 1 0 1 0 ...
#>  $ std_X_rownames: chr [1:62] "row1" "row2" "row3" "row4" ...
#>  $ unpen         : int 1
#>  $ unpen_colnames: chr "sex"
#>  $ ns            : int [1:2001] 1 2 3 4 5 6 7 8 9 10 ...
#>  $ std_X_colnames: chr [1:2001] "sex" "Hsa.3004" "Hsa.13491" "Hsa.13491.1" ...
#>  $ std_X         :Formal class 'big.matrix.descriptor' [package "bigmemory"] with 1 slot
#>   .. ..@ description:List of 13
#>   .. .. ..$ sharedType: chr "FileBacked"
#>   .. .. ..$ filename  : chr "std_colon2.bk"
#>   .. .. ..$ dirname   : chr "/tmp/RtmpZ6kWm2/"
#>   .. .. ..$ totalRows : int 62
#>   .. .. ..$ totalCols : int 2001
#>   .. .. ..$ rowOffset : num [1:2] 0 62
#>   .. .. ..$ colOffset : num [1:2] 0 2001
#>   .. .. ..$ nrow      : num 62
#>   .. .. ..$ ncol      : num 2001
#>   .. .. ..$ rowNames  : NULL
#>   .. .. ..$ colNames  : chr [1:2001] "sex" "Hsa.3004" "Hsa.13491" "Hsa.13491.1" ...
#>   .. .. ..$ type      : chr "double"
#>   .. .. ..$ separated : logi FALSE
#>  $ std_X_n       : num 62
#>  $ std_X_p       : num 2001
#>  $ std_X_center  : num [1:2001] 1.47 7015.79 4966.96 4094.73 3987.79 ...
#>  $ std_X_scale   : num [1:2001] 0.499 3067.926 2171.166 1803.359 2002.738 ...
#>  $ penalty_factor: num [1:2001] 0 1 1 1 1 1 1 1 1 1 ...
#>  - attr(*, "class")= chr "plmm_design"

## Example 3: PLINK data ##
# \donttest{
# process PLINK data
temp_dir <- tempdir()
unzip_example_data(outdir = temp_dir)
#> Unzipped files are saved in /tmp/RtmpZ6kWm2 

plink_data <- process_plink(data_dir = temp_dir,
  data_prefix = "penncath_lite",
  rds_dir = temp_dir,
  rds_prefix = "imputed_penncath_lite",
  # imputing the mode to address missing values
  impute_method = "mode",
  # overwrite existing files in temp_dir
  # (you can turn this feature off if you need to)
  overwrite = TRUE,
  # turning off parallelization - leaving this on causes problems knitting this vignette
  parallel = FALSE)
#> 
#> Preprocessing penncath_lite data:
#> Creating penncath_lite.rds
#> 
#> There are 1401 observations and 4367 genomic features in the specified data files, representing chromosomes 1 - 22 
#> There are a total of 3514 SNPs with missing values
#> Of these, 13 are missing in at least 50% of the samples
#> 
#> Imputing the missing (genotype) values using mode method
#> 
#> process_plink() completed
#> Processed files now saved as /tmp/RtmpZ6kWm2/imputed_penncath_lite.rds

# get outcome data
penncath_pheno <- read.csv(find_example_data(path = 'penncath_clinical.csv'))

outcome <- data.frame(FamID = as.character(penncath_pheno$FamID),
                  CAD = penncath_pheno$CAD)

unpen_predictors <- data.frame(FamID = as.character(penncath_pheno$FamID),
                               sex = penncath_pheno$sex,
                               age = penncath_pheno$age)


# create design where sex and age are always included in the model
pen_design <- create_design(data_file = plink_data,
  feature_id = "FID",
  rds_dir = temp_dir,
  new_file = "std_penncath_lite",
  add_outcome = outcome,
  outcome_id = "FamID",
  outcome_col = "CAD",
  add_predictor = unpen_predictors,
  predictor_id = "FamID",
  logfile = "design",
  # again, overwrite if needed; use with caution
  overwrite = TRUE)
#> 
#> Aligning external data with the feature data by FamID 
#> Adding predictors from external data.
#> Aligning IDs between fam and predictor files
#> Column-wise combining data sets
#>   |                                                                              |                                                                      |   0%  |                                                                              |                                                                      |   1%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |==                                                                    |   4%  |                                                                              |===                                                                   |   4%  |                                                                              |===                                                                   |   5%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |=====                                                                 |   8%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |=======                                                               |  11%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |=========                                                             |  14%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  16%  |                                                                              |============                                                          |  17%  |                                                                              |============                                                          |  18%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |==============                                                        |  21%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  22%  |                                                                              |================                                                      |  23%  |                                                                              |================                                                      |  24%  |                                                                              |=================                                                     |  24%  |                                                                              |=================                                                     |  25%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |=====================                                                 |  31%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |=======================                                               |  34%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |==========================                                            |  38%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  39%  |                                                                              |============================                                          |  40%  |                                                                              |============================                                          |  41%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |==============================                                        |  44%  |                                                                              |===============================                                       |  44%  |                                                                              |===============================                                       |  45%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |=================================                                     |  48%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |===================================                                   |  51%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |=====================================                                 |  54%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  56%  |                                                                              |========================================                              |  57%  |                                                                              |========================================                              |  58%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |==========================================                            |  61%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  62%  |                                                                              |============================================                          |  63%  |                                                                              |============================================                          |  64%  |                                                                              |=============================================                         |  64%  |                                                                              |=============================================                         |  65%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |===============================================                       |  68%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |=================================================                     |  71%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |===================================================                   |  74%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  76%  |                                                                              |======================================================                |  77%  |                                                                              |======================================================                |  78%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  79%  |                                                                              |========================================================              |  80%  |                                                                              |========================================================              |  81%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |==========================================================            |  84%  |                                                                              |===========================================================           |  84%  |                                                                              |===========================================================           |  85%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |=============================================================         |  88%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |===============================================================       |  91%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |=================================================================     |  94%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |====================================================================  |  98%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================|  99%  |                                                                              |======================================================================| 100%
#> There are 62 constant features in the data
#> Subsetting data to exclude constant features (e.g., monomorphic SNPs)
#> Column-standardizing the design matrix...
#> Standardization completed at 2026-01-21 16:54:36
#> Done with standardization. File formatting in progress

# examine the design - notice the components of this object
pen_design_rds <- readRDS(pen_design)

# }

```

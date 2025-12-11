# A function to create a design matrix, outcome, and penalty factor to be passed to a model fitting function

A function to create a design matrix, outcome, and penalty factor to be
passed to a model fitting function

## Usage

``` r
create_design_filebacked(
  data_file,
  rds_dir,
  obj,
  new_file,
  add_outcome,
  outcome_id,
  outcome_col,
  na_outcome_vals = c(-9, NA_integer_),
  feature_id = NULL,
  add_predictor = NULL,
  predictor_id = NULL,
  unpen = NULL,
  logfile = NULL,
  overwrite = FALSE,
  quiet = FALSE
)
```

## Arguments

- data_file:

  A filepath to rds file of processed data (data from
  [`process_plink()`](https://pbreheny.github.io/plmmr/reference/process_plink.md)
  or
  [`process_delim()`](https://pbreheny.github.io/plmmr/reference/process_delim.md))

- rds_dir:

  The path to the directory in which you want to create the new '.rds'
  and '.bk' files.

- obj:

  The RDS object read in by
  [`create_design()`](https://pbreheny.github.io/plmmr/reference/create_design.md)

- new_file:

  User-specified filename (*without .bk/.rds extension*) for the
  to-be-created .rds/.bk files. Must be different from any existing
  .rds/.bk files in the same folder.

- add_outcome:

  A data frame or matrix with two columns: and ID column and a column
  with the outcome value (to be used as 'y' in the final design). IDs
  must be characters, outcome must be numeric.

- outcome_id:

  A string specifying the name of the ID column in 'add_outcome'

- outcome_col:

  A string specifying the name of the phenotype column in 'add_outcome'

- na_outcome_vals:

  A vector of numeric values used to code NA values in the outcome.
  Defaults to `c(-9, NA_integer)` (the -9 matches PLINK conventions).

- feature_id:

  A string specifying the column in the data X (the feature data) with
  the row IDs (e.g., identifiers for each row/sample/participant/,
  etc.). No duplicates allowed. - for PLINK data: a string specifying an
  ID column of the PLINK `.fam` file. Options are "IID" (default) and
  "FID" - for all other filebacked data: a character vector of unique
  identifiers (IDs) for each row of the feature data (i.e., the data
  processed with
  [`process_delim()`](https://pbreheny.github.io/plmmr/reference/process_delim.md)) -
  if left NULL (default), X is assumed to have the same row-order as
  add_outcome. **Note**: if this assumption is made in error,
  calculations downstream will be incorrect. Pay close attention here.

- add_predictor:

  Optional (for PLINK data only): a matrix or data frame to be used for
  adding additional **unpenalized** covariates/predictors/features from
  an external file (i.e., not a PLINK file). This matrix must have one
  column that is an ID column; all other columns aside the ID will be
  used as covariates in the design matrix. Columns must be named.

- predictor_id:

  Optional (for PLINK data only): A string specifying the name of the
  column in 'add_predictor' with sample IDs. **Required** if
  'add_predictor' is supplied. The names will be used to subset and
  align this external covariate with the supplied PLINK data.

- unpen:

  Optional (for delimited file data only): an optional character vector
  with the names of columns to mark as unpenalized (i.e., these features
  would always be included in a model). **Note**: if you choose to use
  this option, X must have column names.

- logfile:

  Optional: name of the '.log' file to be written – **Note:** do not
  append a `.log` to the filename; this is done automatically.

- overwrite:

  Logical: should existing .rds files be overwritten? Defaults to FALSE.

- quiet:

  Logical: should messages to be printed to the console be silenced?
  Defaults to FALSE

## Value

A filepath to the created .rds file containing all the information for
model fitting, including a standardized X and model design information

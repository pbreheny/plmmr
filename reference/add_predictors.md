# A helper function to add predictors to a filebacked matrix of data

A helper function to add predictors to a filebacked matrix of data

## Usage

``` r
add_predictors(obj, add_predictor, id_var, rds_dir, quiet)
```

## Arguments

- obj:

  A `bigSNP` object

- add_predictor:

  Optional: add additional covariates/predictors/features from an
  external file (i.e., not a PLINK file).

- id_var:

  String specifying which column of the PLINK `.fam` file has the unique
  sample identifiers.

- rds_dir:

  The path to the directory in which you want to create the new '.rds'
  and '.bk' files. Defaults to `data_dir`(from
  [`process_plink()`](https://pbreheny.github.io/plmmr/reference/process_plink.md)
  call)

- quiet:

  Logical: should messages be printed to the console? Defaults to FALSE
  (which leaves the print messages on...)

## Value

A list of 2 components:

- 'obj' - a `bigSNP` object with an added element representing the
  matrix that includes the additional predictors as the first few
  columns

- 'non_gen' - an integer vector that ranges from 1 to the number of
  added predictors. Example: if 2 predictors are added, unpen= 1:2

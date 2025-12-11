# A helper function to support `process_plink()`

A helper function to support
[`process_plink()`](https://pbreheny.github.io/plmmr/reference/process_plink.md)

## Usage

``` r
align_ids(id_var, quiet, add_predictor, og_ids)
```

## Arguments

- id_var:

  String specifying the variable name of the ID column

- quiet:

  Logical: should a message be printed?

- add_predictor:

  External data to include in design matrix. This is the
  add_predictors... arg in
  [`process_plink()`](https://pbreheny.github.io/plmmr/reference/process_plink.md)

- og_ids:

  Character vector with the PLINK ids (FID or IID) from the *original*
  data (i.e., the data before any subsetting from handling missing
  phenotypes)

## Value

A matrix with the same dimensions as add_predictor

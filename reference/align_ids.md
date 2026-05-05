# A helper function to support `create_design_filebacked()`

A helper function to support
[`create_design_filebacked()`](https://pbreheny.github.io/plmmr/reference/create_design_filebacked.md)

## Usage

``` r
align_ids(id_var, add_predictor, og_ids, outfile, quiet)
```

## Arguments

- id_var:

  String specifying the variable name of the ID column

- add_predictor:

  External data to include in design matrix. This is the `add_predictor`
  arg in
  [`create_design_filebacked()`](https://pbreheny.github.io/plmmr/reference/create_design_filebacked.md)

- og_ids:

  Character vector with the PLINK ids (FID or IID) from the *original*
  data (i.e., the data before any subsetting from handling missing
  phenotypes)

- outfile:

  A string with the name of the filepath for the log file

- quiet:

  Logical: should console messages be silenced? Defaults to FALSE

## Value

A matrix with the same dimensions as add_predictor

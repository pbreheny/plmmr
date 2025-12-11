# Read in processed data This function is intended to be called after either `process_plink()` or `process_delim()` has been called once.

Read in processed data This function is intended to be called after
either
[`process_plink()`](https://pbreheny.github.io/plmmr/reference/process_plink.md)
or
[`process_delim()`](https://pbreheny.github.io/plmmr/reference/process_delim.md)
has been called once.

## Usage

``` r
get_data(path, returnX = FALSE, trace = TRUE)
```

## Arguments

- path:

  The file path to the RDS object containing the processed data. Do not
  add the '.rds' extension to the path.

- returnX:

  Logical: Should the design matrix be returned as a numeric matrix that
  will be stored in memory. By default, this will be FALSE.

- trace:

  Logical: Should trace messages be shown? Default is TRUE.

## Value

A list with these components:

- std_X, the column-standardized design matrix as either (1) a numeric
  matrix or (2) a filebacked matrix (FBM). See
  [`bigstatsr::FBM()`](https://privefl.github.io/bigstatsr/reference/FBM-class.html)
  and `bigsnpr::bigSnp-class` documentation for details.

- (if PLINK data) fam, a data frame containing the pedigree information
  (like a .fam file in PLINK)

- (if PLINK data) map, a data frame containing the feature information
  (like a .bim file in PLINK)

- ns: A vector indicating the which columns of X contain nonsingular
  features (i.e., features with variance != 0.

- center: A vector of values for centering each column in X

- scale: A vector of values for scaling each column in X

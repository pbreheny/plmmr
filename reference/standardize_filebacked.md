# A helper function to standardize a filebacked matrix

A helper function to standardize a filebacked matrix

## Usage

``` r
standardize_filebacked(X, outfile, quiet, tocenter = TRUE)
```

## Arguments

- X:

  A `big.matrix` object that has been subset &/or had any additional
  predictors appended as columns

- outfile:

  Optional: the name (character string) of the logfile to be written.

- quiet:

  Logical: should console messages be silenced? Defaults to FALSE

- tocenter:

  Should the matrix be centered in addition to scaled? Defaults to TRUE.

## Value

A list with a component called `std_X` - this is a filebacked
`big.matrix` with column-standardized data. List also includes several
other indices/meta-data on the standardized matrix

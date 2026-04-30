# A helper function to count constant features

A helper function to count constant features

## Usage

``` r
count_constant_features(fbm, outfile, quiet)
```

## Arguments

- fbm:

  A filebacked `big.matrix`

- outfile:

  String specifying name of log file

- quiet:

  Logical: should console messages be silenced? Defaults to FALSE

## Value

A numeric vector with the indices of the non-singular columns of the
matrix associated with `fbm`

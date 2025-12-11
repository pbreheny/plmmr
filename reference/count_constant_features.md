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

  Logical: should a message be printed to the console

## Value

ns A numeric vector with the indices of the non-singular columns of the
matrix associated with `counts`

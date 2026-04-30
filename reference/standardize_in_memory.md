# A helper function to standardize matrices

A helper function to standardize matrices

## Usage

``` r
standardize_in_memory(X, tocenter = TRUE)
```

## Arguments

- X:

  A matrix

- tocenter:

  Should the matrix be centered in addition to scaled? Defaults to TRUE.

## Value

A list containing the standardized `X` matrix and associated metadata

## Details

This function is adapted from
https://github.com/pbreheny/ncvreg/blob/master/R/std.R NOTE: this
function returns a matrix **in memory**. For standardizing filebacked
data, use
[`standardize_filebacked()`](https://pbreheny.github.io/plmmr/reference/standardize_filebacked.md)
– see `src/big_standardize.cpp`

# A helper function to standardize matrices

A helper function to standardize matrices

## Usage

``` r
standardize_in_memory(X, tocenter = TRUE)
```

## Arguments

- X:

  a matrix

## Value

a list with the standardized matrix, vectors with the centering/scaling
values, and a vector with the indices of nonsingular columns

## Details

This function is adapted from
https://github.com/pbreheny/ncvreg/blob/master/R/std.R NOTE: this
function returns a matrix **in memory**. For standardizing filebacked
data, use `big_std()` – see src/big_standardize.cpp

# A version of `cbind()` for file-backed matrices

A version of [`cbind()`](https://rdrr.io/r/base/cbind.html) for
file-backed matrices

## Usage

``` r
big_cbind(A, B, C, quiet)
```

## Arguments

- A:

  in-memory data

- B:

  file-backed data

- C:

  file-backed placeholder for combined data

- quiet:

  Logical: should console messages be silenced? Defaults to FALSE

## Value

C, filled in with all column values of A and B combined

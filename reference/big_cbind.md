# A version of `cbind()` for file-backed matrices

A version of [`cbind()`](https://rdrr.io/r/base/cbind.html) for
file-backed matrices

## Usage

``` r
big_cbind(A, B, C)
```

## Arguments

- A:

  in-memory data

- B:

  file-backed data

- C:

  file-backed placeholder for combined data

## Value

C, filled in with all column values of A and B combined

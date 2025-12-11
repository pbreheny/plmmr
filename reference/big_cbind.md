# a version of cbind() for file-backed matrices

a version of cbind() for file-backed matrices

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

  Logical

## Value

C, filled in with all column values of A and B combined

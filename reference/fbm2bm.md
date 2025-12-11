# Functions to convert between FBM and big.matrix type objects

Functions to convert between FBM and big.matrix type objects

## Usage

``` r
fbm2bm(fbm, desc = FALSE)
```

## Arguments

- fbm:

  An FBM object; see
  [`bigstatsr::FBM()`](https://privefl.github.io/bigstatsr/reference/FBM-class.html)
  for details

- desc:

  Logical: is the descriptor file desired (as opposed to the filebacked
  big matrix)? Defaults to FALSE.

## Value

a `big.matrix` - see
[`bigmemory::filebacked.big.matrix()`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html)
for details

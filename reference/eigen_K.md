# A function to take the eigendecomposition of K

Note: This is faster than taking SVD of X when p \>\> n

## Usage

``` r
eigen_K(std_X)
```

## Arguments

- std_X:

  The *standardized* design matrix.

## Value

A list with three elements:

- `s`: The non-zero eigenvalues of K

- `U`: The eigenvectors of K associated with s

- `K`: The fully computed K matrix

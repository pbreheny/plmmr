# A function to take the eigendecomposition of K Note: This is faster than taking SVD of X when p \>\> n

A function to take the eigendecomposition of K Note: This is faster than
taking SVD of X when p \>\> n

## Usage

``` r
eigen_K(std_X, fbm_flag)
```

## Arguments

- std_X:

  The *standardized* design matrix, stored as big.matrix object.

- fbm_flag:

  Logical: is std_X an FBM obejct? Passed from
  [`plmm()`](https://pbreheny.github.io/plmmr/reference/plmm.md).

## Value

A list with the eigenvectors and eigenvalues of K

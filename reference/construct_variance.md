# a function to create the estimated variance matrix from a PLMM fit

a function to create the estimated variance matrix from a PLMM fit

## Usage

``` r
construct_variance(fit, K = NULL, eta = NULL)
```

## Arguments

- fit:

  An object returned by
  [`plmm()`](https://pbreheny.github.io/plmmr/reference/plmm.md)

- K:

  An optional matrix

- eta:

  An optional numeric value between 0 and 1; if `fit` is not supplied,
  then this option must be specified.

## Value

Sigma_hat, a matrix representing the estimated variance

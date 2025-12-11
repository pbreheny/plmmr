# Estimate eta (to be used in rotating the data) This function is called internally by `plmm()`

Estimate eta (to be used in rotating the data) This function is called
internally by
[`plmm()`](https://pbreheny.github.io/plmmr/reference/plmm.md)

## Usage

``` r
estimate_eta(n, s, U, y, eta_star)
```

## Arguments

- n:

  The number of observations

- s:

  The singular values of K, the realized relationship matrix

- U:

  The left-singular vectors of the *standardized* design matrix

- y:

  Continuous outcome vector.

## Value

a numeric value with the estimated value of eta, the variance parameter

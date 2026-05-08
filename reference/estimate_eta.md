# Estimate eta (to be used in rotating the data)

This function is called internally by
[`plmm()`](https://pbreheny.github.io/plmmr/reference/plmm.md)

## Usage

``` r
estimate_eta(n, s, U, y, incpt_flag)
```

## Arguments

- n:

  The number of observations

- s:

  The non-zero eigenvalues of K, the realized relationship matrix

- U:

  The eigenvectors of K associated with s

- y:

  Continuous outcome vector

- incpt_flag:

  Logical: Does the model require fitting an intercept?

## Value

a numeric value with the estimated value of eta, the variance parameter

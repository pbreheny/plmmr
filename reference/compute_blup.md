# a function to compute the BLUP

a function to compute the BLUP

## Usage

``` r
compute_blup(fit, Xb, Sigma_21, idx)
```

## Arguments

- fit:

  An object returned by
  [`plmm()`](https://pbreheny.github.io/plmmr/reference/plmm.md)

- Xb:

  Linear predictor

- Sigma_21:

  Covariance matrix between the training and the testing data. Extracted
  from `estimated_Sigma` that is generated using all observations

- idx:

  Vector of indices of the penalty parameter `lambda` at which
  predictions are required. By default, all indices are returned.

## Value

Sigma_hat, a matrix representing the estimated variance

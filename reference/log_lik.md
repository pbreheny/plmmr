# Evaluate the negative log-likelihood of an intercept-only Gaussian plmm model

This function allows you to evaluate the negative log-likelihood of a
linear mixed model under the assumption of a null model in order to
estimate the variance parameter, eta.

## Usage

``` r
log_lik(eta, n, s, U, y, rot_y = NULL)
```

## Arguments

- eta:

  The proportion of variance in the outcome that is attributable to
  causal SNP effects. In other words, signal-to-noise ratio.

- n:

  The number of observations

- s:

  The singular values of K, the realized relationship matrix

- U:

  The left-singular vectors of the *standardized* design matrix

- y:

  Continuous outcome vector.

- rot_y:

  Optional: if y has already been rotated, then this can be supplied.

## Value

the value of the log-likelihood of the PLMM, evaluated with the supplied
parameters

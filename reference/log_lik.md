# Evaluate the negative log-likelihood of an intercept-only Gaussian plmm model

This function allows you to evaluate the negative log-likelihood of a
linear mixed model under the assumption of a null model in order to
estimate the variance parameter, eta.

## Usage

``` r
log_lik(eta, n, s, U, y, incpt_flag)
```

## Arguments

- eta:

  Estimated proportion of the variance in the outcome attributable to
  population/correlation structure

- n:

  The number of observations

- s:

  The non-zero eigenvalues of K, the realized relationship matrix

- U:

  The eigenvectors of K associated with s

- y:

  Continuous outcome vector

- incpt_flag:

  Logical: Does the model require fitting an intercept? Passed from
  `estimate_eta`.

## Value

the value of the log-likelihood of the PLMM, evaluated with the supplied
parameters

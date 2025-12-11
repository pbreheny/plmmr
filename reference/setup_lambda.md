# Compute sequence of lambda values

This function allows you compute a sequence of lambda values for plmm
models.

## Usage

``` r
setup_lambda(
  X,
  y,
  alpha,
  lambda_min,
  nlambda,
  penalty_factor,
  intercept = TRUE
)
```

## Arguments

- X:

  Rotated and standardized design matrix which *includes* the intercept
  column if present. May include clinical covariates and other non-SNP
  data. This can be either a 'matrix' or 'FBM' object.

- y:

  Continuous outcome vector.

- alpha:

  Tuning parameter for the Mnet estimator which controls the relative
  contributions from the MCP/SCAD penalty and the ridge, or L2 penalty.
  alpha=1 is equivalent to MCP/SCAD penalty, while alpha=0 would be
  equivalent to ridge regression. However, alpha=0 is not supported;
  alpha may be arbitrarily small, but not exactly 0.

- lambda_min:

  The smallest value for lambda, as a fraction of lambda.max. Default is
  .001 if the number of observations is larger than the number of
  covariates and .05 otherwise. A value of lambda_min = 0 is not
  supported.

- nlambda:

  The desired number of lambda values in the sequence to be generated.

- penalty_factor:

  A multiplicative factor for the penalty applied to each coefficient.
  If supplied, penalty_factor must be a numeric vector of length equal
  to the number of columns of X. The purpose of penalty_factor is to
  apply differential penalization if some coefficients are thought to be
  more likely than others to be in the model. In particular,
  penalty_factor can be 0, in which case the coefficient is always in
  the model without shrinkage.

- intercept:

  Logical: does X contain an intercept column? Defaults to TRUE.

## Value

a numeric vector of lambda values, equally spaced on the log scale

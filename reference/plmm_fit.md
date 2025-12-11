# PLMM fit: a function that fits a PLMM using the values returned by plmm_prep()

PLMM fit: a function that fits a PLMM using the values returned by
plmm_prep()

## Usage

``` r
plmm_fit(
  prep,
  y,
  std_X_details,
  penalty_factor,
  fbm_flag,
  penalty,
  gamma = 3,
  alpha = 1,
  lambda_min,
  nlambda = 100,
  lambda,
  eps = 1e-04,
  max_iter = 10000,
  init = NULL,
  warn = TRUE,
  ...
)
```

## Arguments

- prep:

  A list as returned from `plmm_prep`

- y:

  The original (not centered) outcome vector. Need this for intercept
  estimate

- std_X_details:

  A list with components 'center' (values used to center X), 'scale'
  (values used to scale X), and 'ns' (indices for nonsingular columns of
  X)

- penalty_factor:

  A multiplicative factor for the penalty applied to each coefficient.
  If supplied, penalty_factor must be a numeric vector of length equal
  to the number of columns of X. The purpose of penalty_factor is to
  apply differential penalization if some coefficients are thought to be
  more likely than others to be in the model. In particular,
  penalty_factor can be 0, in which case the coefficient is always in
  the model without shrinkage.

- fbm_flag:

  Logical: is std_X an FBM object? Passed from
  [`plmm()`](https://pbreheny.github.io/plmmr/reference/plmm.md).

- penalty:

  The penalty to be applied to the model. Either "MCP" (the default),
  "SCAD", or "lasso".

- gamma:

  The tuning parameter of the MCP/SCAD penalty (see details). Default is
  3 for MCP and 3.7 for SCAD.

- alpha:

  Tuning parameter for the Mnet estimator which controls the relative
  contributions from the MCP/SCAD penalty and the ridge, or L2 penalty.
  alpha=1 is equivalent to MCP/SCAD penalty, while alpha=0 would be
  equivalent to ridge regression. However, alpha=0 is not supported;
  alpha may be arbitrarily small, but not exactly 0.

- lambda_min:

  The smallest value for lambda, as a fraction of lambda.max. Default is
  .001 if the number of observations is larger than the number of
  covariates and .05 otherwise.

- nlambda:

  Length of the sequence of lambda. Default is 100.

- lambda:

  A user-specified sequence of lambda values. By default, a sequence of
  values of length nlambda is computed, equally spaced on the log scale.

- eps:

  Convergence threshold. The algorithm iterates until the RMSD for the
  change in linear predictors for each coefficient is less than eps.
  Default is `1e-4`.

- max_iter:

  Maximum number of iterations (total across entire path). Default is
  10000.

- init:

  Initial values for coefficients. Default is 0 for all columns of X.

- warn:

  Return warning messages for failures to converge and model saturation?
  Default is TRUE.

- ...:

  Additional arguments that can be passed to
  `biglasso::biglasso_simple_path()`

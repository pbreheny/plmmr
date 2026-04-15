# plmm_checks

plmm_checks

## Usage

``` r
plmm_checks(
  design,
  K = NULL,
  eta = NULL,
  penalty = "lasso",
  init = NULL,
  gamma,
  alpha = 1,
  trace = FALSE,
  save_rds = NULL,
  return_fit = TRUE,
  ...
)
```

## Arguments

- design:

  The design object, as created by
  [`create_design()`](https://pbreheny.github.io/plmmr/reference/create_design.md)

- K:

  Similarity matrix used to rotate the data. This should either be (1) a
  known matrix that reflects the covariance of y, (2) an estimate
  (Default is \\\frac{1}{p}(XX^T)\\), or (3) a list with components 'd'
  and 'U', as returned by a previous
  [`plmm()`](https://pbreheny.github.io/plmmr/reference/plmm.md) model
  fit on the same data.

- eta:

  Optional argument to input a specific eta term rather than estimate it
  from the data. If K is a known covariance matrix that is full rank,
  this should be 1.

- penalty:

  The penalty to be applied to the model. Either "MCP" (the default),
  "SCAD", or "lasso".

- init:

  Initial values for coefficients. Default is 0 for all columns of X.

- gamma:

  The tuning parameter of the MCP/SCAD penalty (see details). Default is
  3 for MCP and 3.7 for SCAD.

- alpha:

  Tuning parameter for the Mnet estimator which controls the relative
  contributions from the MCP/SCAD penalty and the ridge, or L2 penalty.
  alpha=1 is equivalent to MCP/SCAD penalty, while alpha=0 would be
  equivalent to ridge regression. However, alpha=0 is not supported;
  alpha may be arbitrarily small, but not exactly 0.

- trace:

  If set to TRUE, inform the user of progress by announcing the
  beginning of each step of the modeling process. Default is FALSE.

- save_rds:

  Optional: if a filepath and name is specified (e.g.,
  `save_rds = "~/dir/my_results.rds"`), then the model results are saved
  to the provided location. Defaults to NULL, which does not save the
  result.

- return_fit:

  Optional: a logical value indicating whether the fitted model should
  be returned as a `plmm` object in the current (assumed interactive)
  session. Defaults to TRUE.

- ...:

  Additional arguments to
  [`get_data()`](https://pbreheny.github.io/plmmr/reference/get_data.md)

## Value

A list of parameters to pass on to model fitting. The list includes the
standardized design matrix, the outcome, and meta-data

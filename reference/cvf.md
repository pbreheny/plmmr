# Cross-validation internal function for cv_plmm

Internal function for cv_plmm which calls plmm on a fold subset of the
original data.

## Usage

``` r
cvf(i, fold, type, cv_args, ...)
```

## Arguments

- i:

  Fold number to be excluded from fit.

- fold:

  n-length vector of fold-assignments.

- type:

  A character argument indicating what should be returned from
  predict.plmm. If `type == 'lp'` predictions are based on the linear
  predictor, `$X beta$`. If `type == 'individual'` predictions are based
  on the linear predictor plus the estimated random effect (BLUP).

- cv_args:

  List of additional arguments to be passed to plmm.

- ...:

  Optional arguments to `predict_within_cv`

## Value

a list with three elements:

- a numeric vector with the loss at each value of lambda

- a numeric value indicating the number of lambda values used

- a numeric value with the predicted outcome (y hat) values at each
  lambda

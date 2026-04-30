# Cross-validation internal function for `cv_plmm()`

Internal function for
[`cv_plmm()`](https://pbreheny.github.io/plmmr/reference/cv_plmm.md)
which calls
[`plmm()`](https://pbreheny.github.io/plmmr/reference/plmm.md) on a fold
subset of the original data.

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
  [`predict.plmm()`](https://pbreheny.github.io/plmmr/reference/predict.plmm.md).
  If `type = 'lp'` predictions are based on the linear predictor, \\X
  \beta\\. If `type = 'individual'` predictions are based on the linear
  predictor plus the estimated random effect (BLUP).

- cv_args:

  List of additional arguments to be passed to plmm.

- ...:

  Optional arguments to
  [`predict_within_cv()`](https://pbreheny.github.io/plmmr/reference/predict_within_cv.md)

## Value

A list with three elements:

- `loss`: a numeric vector with the loss at each value of lambda

- `nl`: a numeric value indicating the number of lambda values used

- `yhat`: a numeric value with the predicted outcome values at each
  lambda

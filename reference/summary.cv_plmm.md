# A summary function for cv_plmm objects

A summary function for cv_plmm objects

## Usage

``` r
# S3 method for class 'cv_plmm'
summary(object, lambda = "min", ...)
```

## Arguments

- object:

  A `cv_plmm` object

- lambda:

  The regularization parameter value at which inference should be
  reported. Can choose a numeric value, 'min', or '1se'. Defaults to
  'min.'

- ...:

  Not used

## Value

The return value is an object with S3 class `summary.cv_plmm`. The class
has its own print method and contains the following list elements:

- `lambda_min`: The lambda value at the minimum cross validation error

- `lambda.1se`: The maximum lambda value within 1 standard error of the
  minimum cross validation error

- `penalty`: The penalty applied to the fitted model

- `nvars`: The number of non-zero coefficients at the selected lambda
  value

- `cve`: The cross validation error at all folds

- `min`: The minimum cross validation error

- `fit`: The `plmm` fit used in the cross validation

if `return_bias_details = TRUE`, two more items are returned:

- `bias`: The mean bias of the cross validation

- `loss`: The loss at each value of `lambda`

## Examples

``` r
admix_design <- create_design(X = admix$X, y = admix$y)
cv_fit <- cv_plmm(design = admix_design)
summary(cv_fit)
#> lasso-penalized model with n=197 and p=101
#> At minimum cross-validation error (lambda=0.2632):
#> -------------------------------------------------
#>   Nonzero coefficients: 5
#>   Cross-validation error (deviance): 2.57
#>   Scale estimate (sigma): 1.602
```

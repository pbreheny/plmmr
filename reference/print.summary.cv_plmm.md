# Print method for summary.cv_plmm objects

Print method for summary.cv_plmm objects

## Usage

``` r
# S3 method for class 'summary.cv_plmm'
print(x, digits, ...)
```

## Arguments

- x:

  An object of class `summary.cv_plmm`

- digits:

  The number of digits to use in formatting output

- ...:

  Not used

## Value

Nothing is returned; instead, a message is printed to the console
summarizing the results of the cross-validated model fit.

## Examples

``` r
admix_design <- create_design(X = admix$X, y = admix$y)
cv_fit <- cv_plmm(design = admix_design)
print(summary(cv_fit))
#> lasso-penalized model with n=197 and p=101
#> At minimum cross-validation error (lambda=0.6204):
#> -------------------------------------------------
#>   Nonzero coefficients: 0
#>   Cross-validation error (deviance): 2.95
#>   Scale estimate (sigma): 1.717
```

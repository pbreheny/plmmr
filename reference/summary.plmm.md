# A summary method for the plmm objects

A summary method for the plmm objects

## Usage

``` r
# S3 method for class 'plmm'
summary(object, lambda, idx, eps = 1e-05, ...)
```

## Arguments

- object:

  An object of class `plmm`

- lambda:

  The regularization parameter value at which inference should be
  reported.

- idx:

  Alternatively, `lambda` may be specified by an index; `idx=10` means:
  report inference for the 10th value of `lambda` along the
  regularization path. If both `lambda` and `idx` are specified,
  `lambda` takes precedence.

- eps:

  If lambda is given, eps is the tolerance for difference between the
  given lambda value and a lambda value from the object. Defaults to
  0.0001 (1e-5)

- ...:

  Not used

## Value

The return value is an object with S3 class `summary.plmm`. The class
has its own print method and contains the following list elements:

- `penalty`: The penalty used by `plmm` (e.g. SCAD, MCP, lasso)

- `n`: Number of instances/observations

- `std_X_n`: the number of observations in the standardized data; the
  only time this would differ from 'n' is if data are from PLINK and the
  external data does not include all the same samples

- `p`: Number of regression coefficients (not including the intercept)

- `converged`: Logical indicator for whether the model converged

- `lambda`: The `lambda` value at which inference is being reported

- `lambda_char`: A formatted character string indicating the lambda
  value

- `nvars`: The number of nonzero coefficients (again, not including the
  intercept) at that value of `lambda`

- `nonzero`: The column names indicating the nonzero coefficients in the
  model at the specified value of `lambda`

## Examples

``` r
admix_design <- create_design(X = admix$X, y = admix$y)
fit <- plmm(design = admix_design)
summary(fit, idx = 97)
#> lasso-penalized regression model with n=197, p=101 at lambda=0.00053
#> -------------------------------------------------
#> The model converged 
#> -------------------------------------------------
#> # of non-zero coefficients:  98 
#> -------------------------------------------------
```

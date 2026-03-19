# Loss method for "plmm" class

Loss method for "plmm" class

## Usage

``` r
plmm_loss(y, yhat)
```

## Arguments

- y:

  Observed outcomes (response) vector

- yhat:

  Predicted outcomes (response) vector

## Value

A numeric vector of the squared-error loss values for the given observed
and predicted outcomes

## Examples

``` r
admix_design <- create_design(X = admix$X, y = admix$y)
fit <- plmm(design = admix_design, K = relatedness_mat(admix$X))
yhat <- predict(object = fit, newX = admix$X, type = 'lp', lambda = 0.05)
head(plmm_loss(yhat = yhat, y = admix$y))
#>            [,1]
#> [1,] 0.61339955
#> [2,] 0.03271895
#> [3,] 0.71148771
#> [4,] 0.24866116
#> [5,] 1.72930587
#> [6,] 2.34539978
```

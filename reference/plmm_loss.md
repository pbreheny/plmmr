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
#> [1,] 0.54465392
#> [2,] 0.02489129
#> [3,] 1.10632045
#> [4,] 0.82823115
#> [5,] 0.79800605
#> [6,] 0.84456943
```

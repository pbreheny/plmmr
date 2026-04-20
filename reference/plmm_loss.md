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
#> [1,] 0.78965760
#> [2,] 0.07225795
#> [3,] 1.21246359
#> [4,] 1.48951064
#> [5,] 0.67003851
#> [6,] 0.60754400
```

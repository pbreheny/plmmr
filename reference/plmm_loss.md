# Loss method for `plmm` class

Loss method for `plmm` class

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
fit <- plmm(design = admix_design)
yhat <- predict(object = fit, newX = admix$X, type = 'lp', lambda = 0.05)
head(plmm_loss(yhat = yhat, y = admix$y))
#>            [,1]
#> [1,] 0.78944147
#> [2,] 0.07330212
#> [3,] 1.22563765
#> [4,] 1.51621494
#> [5,] 0.66752505
#> [6,] 0.60393682
```

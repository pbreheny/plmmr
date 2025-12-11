# A function to print the summary of a `plmm` model

A function to print the summary of a `plmm` model

## Usage

``` r
# S3 method for class 'summary.plmm'
print(x, ...)
```

## Arguments

- x:

  A `summary.plmm` object

- ...:

  Not used

## Value

Nothing is returned; instead, a message is printed to the console
summarizing the results of the model fit.

## Examples

``` r
lam <- rev(seq(0.01, 1, length.out=20)) |> round(2) # for sake of example
admix_design <- create_design(X = admix$X, y = admix$y)
fit <- plmm(design = admix_design, lambda = lam)
fit2 <- plmm(design = admix_design, penalty = "SCAD", lambda = lam)
print(summary(fit, idx = 18))
#> lasso-penalized regression model with n=197, p=101 at lambda=0.1100
#> -------------------------------------------------
#> The model converged 
#> -------------------------------------------------
#> # of non-zero coefficients:  50 
#> -------------------------------------------------
print(summary(fit2, idx = 18))
#> SCAD-penalized regression model with n=197, p=101 at lambda=0.1100
#> -------------------------------------------------
#> The model converged 
#> -------------------------------------------------
#> # of non-zero coefficients:  50 
#> -------------------------------------------------
```

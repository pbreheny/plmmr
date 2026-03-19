# Coef method for "plmm" class

Coef method for "plmm" class

## Usage

``` r
# S3 method for class 'plmm'
coef(object, lambda, which = seq_along(object$lambda), drop = TRUE, ...)
```

## Arguments

- object:

  An object of class "plmm."

- lambda:

  A numeric vector of lambda values.

- which:

  Vector of lambda indices for which coefficients to return.

- drop:

  Logical.

- ...:

  Additional arguments.

## Value

Either a numeric matrix (if model was fit on data stored in memory) or a
sparse matrix (if model was fit on data stored filebacked). Rownames are
feature names, columns are values of `lambda`.

## Examples

``` r
admix_design <- create_design(X = admix$X, y = admix$y)
fit <- plmm(design = admix_design)
coef(fit)[1:10, 41:45]
#>                 0.1849     0.1794     0.1741     0.1689    0.1638
#> (Intercept)  4.5317512  4.5781075  4.6268174  4.6741033  4.725200
#> Snp1        -0.3589324 -0.3729881 -0.3867645 -0.4001278 -0.413018
#> Snp2         0.0000000  0.0000000  0.0000000  0.0000000  0.000000
#> Snp3         1.0848873  1.1587859  1.2298792  1.2988576  1.364996
#> Snp4         0.0000000  0.0000000  0.0000000  0.0000000  0.000000
#> Snp5         0.0000000  0.0000000  0.0000000  0.0000000  0.000000
#> Snp6         0.0000000  0.0000000  0.0000000  0.0000000  0.000000
#> Snp7         0.0000000  0.0000000  0.0000000  0.0000000  0.000000
#> Snp8         0.0000000  0.0000000  0.0000000  0.0000000  0.000000
#> Snp9         0.0000000  0.0000000  0.0000000  0.0000000  0.000000
```

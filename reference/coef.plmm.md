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
#>                 0.1856     0.1801     0.1747     0.1695     0.1644
#> (Intercept)  4.5226834  4.5705611  4.6194859  4.6672104  4.7185077
#> Snp1        -0.3594119 -0.3735471 -0.3873291 -0.4007035 -0.4136123
#> Snp2         0.0000000  0.0000000  0.0000000  0.0000000  0.0000000
#> Snp3         1.0622440  1.1362918  1.2078129  1.2771748  1.3437531
#> Snp4         0.0000000  0.0000000  0.0000000  0.0000000  0.0000000
#> Snp5         0.0000000  0.0000000  0.0000000  0.0000000  0.0000000
#> Snp6         0.0000000  0.0000000  0.0000000  0.0000000  0.0000000
#> Snp7         0.0000000  0.0000000  0.0000000  0.0000000  0.0000000
#> Snp8         0.0000000  0.0000000  0.0000000  0.0000000  0.0000000
#> Snp9         0.0000000  0.0000000  0.0000000  0.0000000  0.0000000
```

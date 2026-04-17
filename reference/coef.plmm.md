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
#>                 0.02632    0.02454     0.02289     0.02135     0.01991
#> (Intercept)  6.52428850  6.5546722  6.58301479  6.60771523  6.63250278
#> Snp1        -0.75533945 -0.7622179 -0.76862980 -0.77458003 -0.78006802
#> Snp2         0.17396681  0.1809709  0.18749609  0.19357674  0.19920799
#> Snp3         2.93136526  2.9599321  2.98654974  3.01160079  3.03501122
#> Snp4         0.10959159  0.1149682  0.11998422  0.12470256  0.12918443
#> Snp5         0.26914495  0.2944826  0.31809271  0.34007337  0.36054725
#> Snp6        -0.06785054 -0.0717795 -0.07544617 -0.07886263 -0.08203357
#> Snp7         0.10025728  0.1051686  0.10974791  0.11400506  0.11793499
#> Snp8         0.00000000  0.0000000  0.00000000  0.00000000  0.00000000
#> Snp9         0.20099143  0.2056513  0.20999860  0.21407593  0.21787670
```

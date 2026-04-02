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
#>                 0.02629     0.02451     0.02286     0.02132     0.01988
#> (Intercept)  6.52434653  6.55474052  6.58309378  6.60786351  6.63231243
#> Snp1        -0.75495057 -0.76186423 -0.76830950 -0.77429085 -0.77980403
#> Snp2         0.17412955  0.18112740  0.18764712  0.19372261  0.19935811
#> Snp3         2.93730059  2.96546704  2.99170567  3.01640229  3.03952138
#> Snp4         0.10863048  0.11405949  0.11912434  0.12388886  0.12841085
#> Snp5         0.26942185  0.29474498  0.31834369  0.34031413  0.36076774
#> Snp6        -0.06796699 -0.07188781 -0.07554629 -0.07895494 -0.08212034
#> Snp7         0.10021628  0.10513369  0.10971862  0.11398080  0.11791722
#> Snp8         0.00000000  0.00000000  0.00000000  0.00000000  0.00000000
#> Snp9         0.20155652  0.20617793  0.21048898  0.21453268  0.21830578
```

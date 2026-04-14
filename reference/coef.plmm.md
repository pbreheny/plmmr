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
#>                 0.02632     0.02454     0.02289     0.02135     0.01991
#> (Intercept)  6.52429230  6.55467658  6.58301968  6.60772387  6.63249260
#> Snp1        -0.75531802 -0.76219841 -0.76861216 -0.77456411 -0.78005347
#> Snp2         0.17397609  0.18097985  0.18750469  0.19358504  0.19921649
#> Snp3         2.93169702  2.96024149  2.98683792  3.01186920  3.03526329
#> Snp4         0.10953812  0.11491760  0.11993637  0.12465729  0.12914142
#> Snp5         0.26916094  0.29449773  0.31810720  0.34008725  0.36055998
#> Snp6        -0.06785712 -0.07178562 -0.07545183 -0.07886785 -0.08203847
#> Snp7         0.10025514  0.10516677  0.10974640  0.11400383  0.11793410
#> Snp8         0.00000000  0.00000000  0.00000000  0.00000000  0.00000000
#> Snp9         0.20102298  0.20568071  0.21002598  0.21410143  0.21790065
```

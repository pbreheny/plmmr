# Coef method for "cv_plmm" class

Coef method for "cv_plmm" class

## Usage

``` r
# S3 method for class 'cv_plmm'
coef(object, lambda, which = object$min, ...)
```

## Arguments

- object:

  An object of class "cv_plmm."

- lambda:

  A numeric vector of lambda values.

- which:

  Vector of lambda indices for which coefficients to return. Defaults to
  lambda index with minimum CVE.

- ...:

  Additional arguments (not used).

## Value

Returns a named numeric vector. Values are the coefficients of the model
at the specified value of either `lambda` or `which`. Names are the
values of `lambda`.

## Examples

``` r
cv_fit <- cv_plmm(admix$X, admix$y, return_fit = TRUE)
head(coef(cv_fit))
#> (Intercept)        Snp1        Snp2        Snp3        Snp4        Snp5 
#>    4.145401    0.000000    0.000000    0.000000    0.000000    0.000000 
```

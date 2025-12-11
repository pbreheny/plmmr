# Calculate a relatedness matrix

Given a matrix of genotypes, this function estimates the genetic
relatedness matrix (GRM, also known as the RRM, see Hayes et al. 2009,
[doi:10.1017/S0016672308009981](https://doi.org/10.1017/S0016672308009981)
) among the subjects: XX'/p, where X is standardized.

## Usage

``` r
relatedness_mat(X, std = TRUE, fbm = FALSE, ns = NULL, ...)
```

## Arguments

- X:

  An n x p numeric matrix of genotypes (from *fully-imputed* data).
  Note: This matrix should *not* include non-genetic features.

- std:

  Logical: should X be standardized? If you set this to FALSE (which can
  only be done if data are stored in memory), you should have a good
  reason for doing so, as standardization is a best practice.

- fbm:

  Logical: is X stored as an FBM? Defaults to FALSE

- ns:

  Optional vector of values indicating the indices of nonsingular
  features

- ...:

  Other optional arguments to `bigstatsr::bigapply()` (like
  `ncores = ...`)

## Value

An n x n numeric matrix capturing the genomic relatedness of the samples
represented in `X`. In our notation, we call this matrix K for
'kinship'; this is also known as the GRM or RRM.

## Examples

``` r
RRM <- relatedness_mat(X = admix$X)
RRM[1:5, 1:5]
#>             [,1]        [,2]        [,3]        [,4]        [,5]
#> [1,]  0.81268908 -0.09098097 -0.07888910  0.06770613  0.08311777
#> [2,] -0.09098097  0.81764801  0.20480021  0.02112812 -0.02640295
#> [3,] -0.07888910  0.20480021  0.82177986 -0.02864226  0.18693970
#> [4,]  0.06770613  0.02112812 -0.02864226  0.89327266 -0.03541470
#> [5,]  0.08311777 -0.02640295  0.18693970 -0.03541470  0.79589686
```

# Calculate a relatedness matrix

Given a matrix of genotypes, this function estimates the genetic
relatedness matrix (GRM, also known as the RRM, see Hayes et al. 2009,
[doi:10.1017/S0016672308009981](https://doi.org/10.1017/S0016672308009981)
) among the subjects: \\\frac{1}{p}(XX^T)\\, where X is standardized.

## Usage

``` r
relatedness_mat(X, std = TRUE)
```

## Arguments

- X:

  An n x p numeric matrix of genotypes (from *fully-imputed* data). Can
  be a filebacked `big.matrix` object. Note: This matrix should *not*
  include non-genetic features.

- std:

  Logical: should `X` be standardized? If you set this to FALSE, you
  should have a good reason for doing so, as standardization is a best
  practice.

## Value

An n x n numeric matrix capturing the genomic relatedness of the samples
represented in `X`. In our notation, we call this matrix `K` for
'kinship'; this is also known as the GRM or RRM.

## Examples

``` r
RRM <- relatedness_mat(X = admix$X)
RRM[1:5, 1:5]
#>             [,1]        [,2]        [,3]        [,4]        [,5]
#> [1,]  0.82927457 -0.09283772 -0.08049908  0.06908789  0.08481405
#> [2,] -0.09283772  0.83433471  0.20897980  0.02155931 -0.02694179
#> [3,] -0.08049908  0.20897980  0.83855088 -0.02922680  0.19075479
#> [4,]  0.06908789  0.02155931 -0.02922680  0.91150272 -0.03613745
#> [5,]  0.08481405 -0.02694179  0.19075479 -0.03613745  0.81213966
```

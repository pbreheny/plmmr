# PLMM prep: a function to run checks, SVD, and rotation prior to fitting a PLMM model This is an internal function for `cv_plmm`

PLMM prep: a function to run checks, SVD, and rotation prior to fitting
a PLMM model This is an internal function for `cv_plmm`

## Usage

``` r
plmm_prep(
  std_X,
  std_X_n,
  std_X_p,
  n,
  p,
  centered_y,
  K = NULL,
  eta = NULL,
  fbm_flag,
  trace = NULL,
  ...
)
```

## Arguments

- std_X:

  Column standardized design matrix. May include clinical covariates and
  other non-SNP data.

- std_X_n:

  The number of observations in std_X (integer)

- std_X_p:

  The number of features in std_X (integer)

- n:

  The number of instances in the *original* design matrix X. This should
  not be altered by standardization.

- p:

  The number of features in the *original* design matrix X, including
  constant features

- centered_y:

  Continuous outcome vector, centered.

- K:

  Similarity matrix used to rotate the data. This should either be: (1)
  a known matrix that reflects the covariance of y, (2) an estimate
  (Default is \\\frac{1}{p}(XX^T)\\), or (3) a list with components 's'
  and 'U', as returned by a previous
  [`plmm()`](https://pbreheny.github.io/plmmr/reference/plmm.md) model
  fit on the same data.

- eta:

  Optional argument to input a specific eta term rather than estimate it
  from the data. If K is a known covariance matrix that is full rank,
  this should be 1.

- fbm_flag:

  Logical: is std_X an FBM type object? This is set internally by
  [`plmm()`](https://pbreheny.github.io/plmmr/reference/plmm.md).

- trace:

  If set to TRUE, inform the user of progress by announcing the
  beginning of each step of the modeling process. Default is FALSE.

- ...:

  Not used yet

## Value

List with these components:

- centered_y: The vector of centered outcomes

- std_X: standardized design matrix

- K: a list with 2 elements. (1) s: vector with the eigenvalues of K,
  and (2) U: the eigenvectors of K (same as left singular values of X).

- eta: the numeric value of the estimated eta parameter

- trace: logical.

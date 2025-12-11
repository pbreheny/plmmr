# Untransform coefficient values back to the original scale

This function unwinds the initial standardization of the data to obtain
coefficient values on their original scale. It is called by
plmm_format().

## Usage

``` r
untransform(
  std_scale_beta,
  p,
  std_X_details,
  fbm_flag,
  plink_flag,
  use_names = TRUE
)
```

## Arguments

- std_scale_beta:

  The estimated coefficients on the standardized scale

- p:

  The number of columns in the original design matrix

- std_X_details:

  A list with 3 elements describing the standardized design matrix
  BEFORE rotation; this should have elements 'scale', 'center', and 'ns'

- fbm_flag:

  Logical: is the corresponding design matrix filebacked?

- plink_flag:

  Logical: did these data come from PLINK files? **Note**: This flag
  matters because of how non-genomic features are handled for PLINK
  files – in data from PLINK files, unpenalized columns are *not*
  counted in the `p` argument. For delimited files, `p` does include
  unpenalized columns. This difference has implications for how the
  `untransform()` function determines the appropriate dimensions for the
  estimated coefficient matrix it returns.

- use_names:

  Logical: should names be added? Defaults to TRUE. Set to FALSE inside
  of [`cvf()`](https://pbreheny.github.io/plmmr/reference/cvf.md)
  helper, as 'ns' will vary within CV folds.

## Value

a matrix of estimated coeffcients, 'beta_vals', that is on the scale of
the original data.

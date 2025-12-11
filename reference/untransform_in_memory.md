# Untransform coefficient values back to the original scale *In memory*

This function unwinds the initial standardization of the data to obtain
coefficient values on their original scale. It is called by
plmm_format().

## Usage

``` r
untransform_in_memory(std_scale_beta, p, std_X_details, use_names = TRUE)
```

## Arguments

- std_scale_beta:

  The estimated coefficients on the standardized scale

- p:

  The number of columns in the original design matrix

- std_X_details:

  A list with 3 elements describing the standardized design matrix
  BEFORE rotation; this should have elements 'scale', 'center', and 'ns'

- use_names:

  Logical: should names be added? Defaults to TRUE. Set to FALSE inside
  of [`cvf()`](https://pbreheny.github.io/plmmr/reference/cvf.md)
  helper, as 'ns' will vary within CV folds.

## Value

a matrix of estimated coeffcients, 'beta_vals', that is on the scale of
the original data.

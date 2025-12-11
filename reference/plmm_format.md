# PLMM format: a function to format the output of a model constructed with `plmm_fit`

PLMM format: a function to format the output of a model constructed with
`plmm_fit`

## Usage

``` r
plmm_format(fit, p, std_X_details, fbm_flag, plink_flag)
```

## Arguments

- fit:

  A list of parameters describing the output of a model constructed with
  `plmm_fit`

- p:

  The number of features in the original data (including constant
  features)

- std_X_details:

  A list with 3 items: \* 'center': the centering values for the columns
  of `X` \* 'scale': the scaling values for the non-singular columns of
  `X` \* 'ns': indicesof nonsingular columns in `std_X`

- fbm_flag:

  Logical: is the corresponding design matrix filebacked? Passed from
  [`plmm()`](https://pbreheny.github.io/plmmr/reference/plmm.md).

- plink_flag:

  Logical: did these data come from PLINK files? **Note**: This flag
  matters because of how non-genomic features are handled for PLINK
  files – in data from PLINK files, unpenalized columns are *not*
  counted in the `p` argument. For delimited files, `p` does include
  unpenalized columns. This difference has implications for how the
  [`untransform()`](https://pbreheny.github.io/plmmr/reference/untransform.md)
  function determines the appropriate dimensions for the estimated
  coefficient matrix it returns.

## Value

A list with the components:

- `beta_vals`: the matrix of estimated coefficients on the original
  scale. Rows are predictors, columns are values of `lambda`

- `lambda`: a numeric vector of the lasso tuning parameter values used
  in model fitting.

- `eta`: a number (double) between 0 and 1 representing the estimated
  proportion of the variance in the outcome attributable to
  population/correlation structure.

- `s`: a vectof of the eigenvalues of relatedness matrix `K`; see
  [`relatedness_mat()`](https://pbreheny.github.io/plmmr/reference/relatedness_mat.md)
  for details.

- `U`: a matrix of the eigenvalues of relatedness matrix `K`

- `rot_y`: the vector of outcome values on the rotated scale. This is
  the scale on which the model was fit.

- `linear_predictors`: the matrix resulting from the product of
  `stdrot_X` and the estimated coefficients on the ~rotated~ scale.

- `penalty`: character string indicating the penalty with which the
  model was fit (e.g., 'MCP')

- `gamma`: numeric value indicating the tuning parameter used for the
  SCAD or lasso penalties was used. Not relevant for lasso models.

- `alpha`: numeric value indicating the elastic net tuning parameter.

- `loss`: vector with the numeric values of the loss at each value of
  `lambda` (calculated on the ~rotated~ scale)

- `penalty_factor`: vector of indicators corresponding to each
  predictor, where 1 = predictor was penalized.

- `ns_idx`: vector with the indices of predictors which were nonsingular
  features (i.e., had variation).

- `iter`: numeric vector with the number of iterations needed in model
  fitting for each value of `lambda`

- `converged`: vector of logical values indicating whether the model
  fitting converged at each value of `lambda`

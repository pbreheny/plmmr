# Predict method to use in cross-validation (within `cvf()`)

Predict method to use in cross-validation (within
[`cvf()`](https://pbreheny.github.io/plmmr/reference/cvf.md))

## Usage

``` r
predict_within_cv(fit, testX, type, fbm = FALSE, Sigma_21 = NULL)
```

## Arguments

- fit:

  A list with the components returned by `plmm_fit`.

- testX:

  A design matrix used for computing predicted values (i.e, the test
  data).

- type:

  A character argument indicating what type of prediction should be
  returned. Passed from
  [`cvf()`](https://pbreheny.github.io/plmmr/reference/cvf.md). Options
  are "lp," "coefficients," "vars," "nvars," and "blup." See details.

- fbm:

  Logical: is `trainX` a filebacked `big.matrix` object? If so, this
  function expects that `testX` is also an FBM. The two X matrices must
  be stored the same way.

- Sigma_21:

  Covariance matrix between the training and the testing data. Required
  if `type == 'blup'`.

## Value

A numeric vector of predicted values

## Details

- `lp` (linear predictor): uses the product of `testX` and the beta
  coefficients of `fit` to predict new values of the outcome. This does
  not incorporate the correlation structure of the data.

- `blup` (acronym for Best Linear Unbiased Predictor): adds to the `lp`
  a value that represents the estimated random effect. This addition is
  a way of incorporating the estimated correlation structure of data
  into our prediction of the outcome.

- `coefficients`: returns the estimated beta-hat

- `vars`: returns the *indices* of variables (e.g., SNPs) with nonzero
  coefficients at each value of lambda. EXCLUDES intercept.

- `nvars`: returns the *number* of variables (e.g., SNPs) with nonzero
  coefficients at each value of lambda. EXCLUDES intercept.

Note: the main difference between this function and the
[`predict.plmm()`](https://pbreheny.github.io/plmmr/reference/predict.plmm.md)
method is that here in CV, the standardized testing data (`std_test_X`),
`Sigma_11`, and `Sigma_21` are calculated in
[`cvf()`](https://pbreheny.github.io/plmmr/reference/cvf.md) instead of
the function defined here.

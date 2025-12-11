# Predict method to use in cross-validation (within `cvf`)

Predict method to use in cross-validation (within `cvf`)

## Usage

``` r
predict_within_cv(
  fit,
  testX,
  type,
  fbm = FALSE,
  Sigma_11 = NULL,
  Sigma_21 = NULL
)
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
  [`cvf()`](https://pbreheny.github.io/plmmr/reference/cvf.md), Options
  are "lp," "coefficients," "vars," "nvars," and "blup." See details.

- fbm:

  Logical: is trainX an FBM object? If so, this function expects that
  testX is also an FBM. The two X matrices must be stored the same way.

- Sigma_11:

  Variance-covariance matrix of the training data. Extracted from
  `estimated_Sigma` that is generated using all observations. Required
  if `type == 'blup'`.

- Sigma_21:

  Covariance matrix between the training and the testing data. Extracted
  from `estimated_Sigma` that is generated using all observations.
  Required if `type == 'blup'`.

## Value

A numeric vector of predicted values

## Details

Define beta-hat as the coefficients estimated at the value of lambda
that minimizes cross-validation error (CVE). Then options for `type` are
as follows:

- 'lp' (default): uses the linear predictor (i.e., product of test data
  and estimated coefficients) to predict test values of the outcome.
  Note that this approach does not incorporate the correlation structure
  of the data.

- 'blup' (acronym for Best Linear Unbiased Predictor): adds to the 'lp'
  a value that represents the estimated random effect. This addition is
  a way of incorporating the estimated correlation structure of data
  into our prediction of the outcome.

Note: the main difference between this function and the
[`predict.plmm()`](https://pbreheny.github.io/plmmr/reference/predict.plmm.md)
method is that here in CV, the standardized testing data (std_test_X),
Sigma_11, and Sigma_21 are calculated in
[`cvf()`](https://pbreheny.github.io/plmmr/reference/cvf.md) instead of
the function defined here.

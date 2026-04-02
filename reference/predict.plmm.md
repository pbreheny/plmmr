# Predict method for plmm class

Predict method for plmm class

## Usage

``` r
# S3 method for class 'plmm'
predict(
  object,
  newX,
  type = c("blup", "coefficients", "vars", "nvars", "lp"),
  X = NULL,
  lambda,
  idx = seq_along(object$lambda),
  ...
)
```

## Arguments

- object:

  An object of class `plmm`.

- newX:

  Matrix of values at which predictions are to be made (not used for
  `type="coefficients"` or for some of the `type` settings in
  `predict`). This can be either a FBM object or a 'matrix' object.
  **Note**: Columns of this argument must be named!

- type:

  A character argument indicating what type of prediction should be
  returned. Options are "lp," "coefficients," "vars," "nvars," and
  "blup." See details.

- X:

  Optional: if `type = 'blup'` and the model was fit in-memory, the
  design matrix used to fit the model represented in `object` must be
  supplied. When supplied, this design matrix will be standardized using
  the center/scale values in `object$std_X_details`, so please **do
  not** standardize this matrix before supplying here. **Note**: If the
  model was fit file-backed, then the filepath to the .bk file with this
  standardized design matrix is returned as 'std_X' in the fit supplied
  to 'object'.

- lambda:

  A numeric vector of regularization parameter `lambda` values at which
  predictions are requested.

- idx:

  Vector of indices of the penalty parameter `lambda` at which
  predictions are required. By default, all indices are returned.

- ...:

  Additional optional arguments

## Value

Depends on the `type` - see Details

## Details

Define beta-hat as the coefficients estimated at the value of lambda
that minimizes cross-validation error (CVE). Then options for `type` are
as follows:

- 'response' (default): uses the product of newX and beta-hat to predict
  new values of the outcome. This does not incorporate the correlation
  structure of the data. For the stats folks out there, this is simply
  the linear predictor.

- 'blup' (acronym for Best Linear Unbiased Predictor): adds to the
  'response' a value that represents the esetimated random effect. This
  addition is a way of incorporating the estimated correlation structure
  of data into our prediction of the outcome.

- 'coefficients': returns the estimated beta-hat

- 'vars': returns the *indices* of variables (e.g., SNPs) with nonzero
  coefficients at each value of lambda. EXCLUDES intercept.

- 'nvars': returns the *number* of variables (e.g., SNPs) with nonzero
  coefficients at each value of lambda. EXCLUDES intercept.

## Examples

``` r
set.seed(123)
train_idx <- sample(1:nrow(admix$X), 100)
# Note: ^ shuffling is important here! Keeps test and train groups comparable.
train <- list(X = admix$X[train_idx,], y = admix$y[train_idx])
train_design <- create_design(X = train$X, y = train$y)

test <- list(X = admix$X[-train_idx,], y = admix$y[-train_idx])
fit <- plmm(design = train_design)

# make predictions for all lambda values
 pred1 <- predict(object = fit, newX = test$X, type = "lp")
 pred2 <- predict(object = fit, newX = test$X, type = "blup", X = train$X)

# look at mean squared prediction error
mspe <- apply(pred1, 2, function(c){crossprod(test$y - c)/length(c)})
min(mspe)
#> [1] 2.819134

mspe_blup <- apply(pred2, 2, function(c){crossprod(test$y - c)/length(c)})
min(mspe_blup) # BLUP is better
#> [1] 2.127884

# compare the MSPE of our model to a null model, for reference
# null model = intercept only -> y_hat is always mean(y)
crossprod(mean(test$y) - test$y)/length(test$y)
#>          [,1]
#> [1,] 6.381748

```

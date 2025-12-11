# Package index

## Fitting penalized linear mixed models (PLMMs - the featured presentation)

- [`plmm()`](https://pbreheny.github.io/plmmr/reference/plmm.md) : Fit a
  linear mixed model via penalized maximum likelihood.

## Loss and cross-validation for PLMMs

- [`cv_plmm()`](https://pbreheny.github.io/plmmr/reference/cv_plmm.md) :
  Cross-validation for plmm
- [`plmm_loss()`](https://pbreheny.github.io/plmmr/reference/plmm_loss.md)
  : Loss method for "plmm" class

## Coefficient methods for PLMMs

- [`coef(`*`<cv_plmm>`*`)`](https://pbreheny.github.io/plmmr/reference/coef.cv_plmm.md)
  : Coef method for "cv_plmm" class
- [`coef(`*`<plmm>`*`)`](https://pbreheny.github.io/plmmr/reference/coef.plmm.md)
  : Coef method for "plmm" class

## Data (pre)processing and wrangling

- [`create_design()`](https://pbreheny.github.io/plmmr/reference/create_design.md)
  : a function to create a design for PLMM modeling

- [`find_example_data()`](https://pbreheny.github.io/plmmr/reference/find_example_data.md)
  : A function to help with accessing example PLINK files

- [`process_plink()`](https://pbreheny.github.io/plmmr/reference/process_plink.md)
  :

  Preprocess PLINK files using the `bigsnpr` package

- [`process_delim()`](https://pbreheny.github.io/plmmr/reference/process_delim.md)
  : A function to read in large data files as an FBM

- [`relatedness_mat()`](https://pbreheny.github.io/plmmr/reference/relatedness_mat.md)
  : Calculate a relatedness matrix

- [`unzip_example_data()`](https://pbreheny.github.io/plmmr/reference/unzip_example_data.md)
  :

  For Linux/Unix and MacOS only, here is a companion function to unzip
  the .gz files that ship with the `plmmr` package

## Plotting, summarizing, and formatting

- [`summary(`*`<plmm>`*`)`](https://pbreheny.github.io/plmmr/reference/summary.plmm.md)
  : A summary method for the plmm objects

- [`summary(`*`<cv_plmm>`*`)`](https://pbreheny.github.io/plmmr/reference/summary.cv_plmm.md)
  : A summary function for cv_plmm objects

- [`plot(`*`<plmm>`*`)`](https://pbreheny.github.io/plmmr/reference/plot.plmm.md)
  : Plot method for plmm class

- [`plot(`*`<cv_plmm>`*`)`](https://pbreheny.github.io/plmmr/reference/plot.cv_plmm.md)
  : Plot method for cv_plmm class

- [`print(`*`<summary.cv_plmm>`*`)`](https://pbreheny.github.io/plmmr/reference/print.summary.cv_plmm.md)
  : Print method for summary.cv_plmm objects

- [`print(`*`<summary.plmm>`*`)`](https://pbreheny.github.io/plmmr/reference/print.summary.plmm.md)
  :

  A function to print the summary of a `plmm` model

## Prediction

- [`predict(`*`<plmm>`*`)`](https://pbreheny.github.io/plmmr/reference/predict.plmm.md)
  : Predict method for plmm class

## Data sets

- [`admix`](https://pbreheny.github.io/plmmr/reference/admix.md) :
  Admix: Semi-simulated SNP data

<!-- badges: start -->

[![GitHub version](https://img.shields.io/static/v1?label=GitHub&message=3.2.0.0&color=blue&logo=github)](https://github.com/pbreheny/plmmr) [![R-CMD-check](https://github.com/pbreheny/plmmr/workflows/R-CMD-check/badge.svg)](https://github.com/pbreheny/plmmr/actions) [![Codecov test coverage](https://codecov.io/gh/pbreheny/plmmr/branch/master/graph/badge.svg)](https://app.codecov.io/gh/pbreheny/plmmr?branch=master)

<!-- badges: end -->

## plmmr <img src="man/figures/plmmr_hex_sticker.png" align="right" width="150"/>

The `plmmr` (**p**enalized **l**inear **m**ixed **m**odels in **R**) package contains functions that fit penalized linear mixed models to correct for unobserved confounding effects.

## Installation

To install the latest version of the package:

``` r
devtools::install_github("pbreheny/plmmr")
```

For a description of the motivation of the functions in this package (along with examples) refer to the second module of [this GWAS data tutorial](https://pbreheny.github.io/adv-gwas-tutorial/index.html)

## Minimal example

``` r
library(plmmr)
X <- rnorm(100*20) |> matrix(100, 20)
y <- rnorm(100)
fit <- plmm(X, y) 
plot(fit)

cvfit <- cv_plmm(X, y)
plot(cvfit)
summary(cvfit)
```

## So how fast is `plmmr`? And how well does it scale?

To illustrate these important questions, I created a separate [GitHub repository](https://github.com/tabpeter/demo_plmmr/tree/master) that has all the scripts for a `plmmr` workflow using publicly-available genome-wide association (GWAS) data. The main takeaway: using GWAS data from a study with 1,400 samples and 800,000 SNPs, a full `plmmr` analysis will run in about half an hour using a single core on a laptop.

Three smaller datasets ship with `plmmr`, and tutorials walking through how to analyze these data sets are documented in the [documentation site](https://pbreheny.github.io/plmmr/). While these datasets are useful for didactic purposes, they are not large enough to really highlight the computational scalability of `plmmr` -- this is what motivated the creation of the separate repository for a GWAS workflow.

## Note on branches

The branches of this repo are organized in the following way:

-   `master` is the main (or 'head') branch.

-   `gh_pages` is where we are keeping all the documentation for `plmmr`

-   `gwas_scale` is an **archived** branch that contains the development version of the package I used to run my dissertation analysis. Will delete this eventually.

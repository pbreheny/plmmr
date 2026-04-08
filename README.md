<!-- badges: start -->
[![GitHub version](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/pbreheny/plmmr/master/.version.json&style=flat&logo=github)](https://github.com/pbreheny/plmmr)
[![CRAN version](https://img.shields.io/cran/v/plmmr?logo=R)](https://cran.r-project.org/package=plmmr)
[![R-CMD-check](https://github.com/pbreheny/plmmr/workflows/R-CMD-check/badge.svg)](https://github.com/pbreheny/plmmr/actions) 
[![Codecov test coverage](https://codecov.io/gh/pbreheny/plmmr/branch/main/graph/badge.svg)](https://app.codecov.io/gh/pbreheny/plmmr?branch=main)
<!-- badges: end -->

## plmmr <img src="man/figures/plmmr_hex_sticker.png" align="right" width="150"/>

The `plmmr` (**p**enalized **l**inear **m**ixed **m**odels in **R**) package contains functions that fit penalized linear mixed models to correct for unobserved confounding effects.

Three small datasets ship with `plmmr`, and tutorials walking through how to analyze these data sets are documented in the [plmmr website](https://pbreheny.github.io/plmmr/).

## Installation

To install the latest version of the package from GitHub, use this:

``` r
devtools::install_github("pbreheny/plmmr")
```

You can also install `plmmr` from CRAN: 

```r
install.packages('plmmr')
```

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

These questions are addressed in our [manuscript describing plmmr](https://doi.org/10.1093/bib/bbaf672), along with its accompanying [GitHub repository](https://github.com/tabpeter/reproduce_plmmr_manuscript/). However, using GWAS data from a study with 1,400 samples and 800,000 SNPs, a full `plmmr` analysis will run in about half an hour using a single core on a laptop.

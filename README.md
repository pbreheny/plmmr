<!-- badges: start -->
[![GitHub version](https://img.shields.io/static/v1?label=GitHub&message=3.0.0&color=blue&logo=github)](https://github.com/pbreheny/plmmr)
[![R-CMD-check](https://github.com/pbreheny/plmmr/workflows/R-CMD-check/badge.svg)](https://github.com/pbreheny/plmmr/actions)
[![Codecov test coverage](https://codecov.io/gh/pbreheny/plmmr/branch/master/graph/badge.svg)](https://app.codecov.io/gh/pbreheny/plmmr?branch=master)
<!-- badges: end -->

## plmmr <img src="man/figures/plmmr_hex_sticker.png" align="right" alt="" width="150" />

The `plmmr` (**p**enalized **l**inear **m**ixed **m**odels in **R**) package contains functions that fit penalized linear mixed models to correct for unobserved confounding effects.

🚧🛠️ **Note**: this package is still under construction 🛠️🚧, as both the package and its underlying methodology are the core of my (Tabitha's 👷‍♀️) in-progress dissertation work. Will keep this page updated as I make progress. The dream is for this package to be able to fit penalized regression models in GWAS-scale data. 

## Installation 

To install the latest version of the package: 

```r
devtools::install_github("pbreheny/plmmr")
```

For a description of the motivation of the functions in this package (along with examples) refer to the second module of [this GWAS data tutorial](https://pbreheny.github.io/adv-gwas-tutorial/index.html)

## Note on branches 

The branches of this repo are organized in the following way: 

  - `master` is the main (or 'head') branch.

  - `gwas_scale` is where I've been working to improve `plmmr` functionality for analyzing genome-wide scale data. This will bring some notable user-facing changes; stay tuned for more updates. 
    
  - `gh_pages` is where we are keeping all the documentation for `plmmr`
  
  - `sign_flip` is an **archived** branch where we have examined the issues caused by +/- signs being flipped as part of truncated SVD.

  

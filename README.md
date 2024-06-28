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

  - `master` is the main branch with all the latest updates.

  - `scale_up` is where we are working to improve `plmm()`'s ability to scale up to larger datasets. In particular, I am in the process of converting more functions to `C++`/adjusting how `plmmr` interacts with [file-backed data](https://pbreheny.github.io/plmmr/articles/filebacking.html). Stay tuned for more on this.
    
  - `gh_pages` is where we are keeping all the documentation for `plmmr`
  
  - `estimate_eta` is an **archived** branch where we worked through an alternative approach for estimating $\eta$. We're keeping this around for reference as we are writing. 
  
  - `sign_flip` is an **archived** branch where we have examined the issues caused by +/- signs being flipped as part of truncated SVD.

  

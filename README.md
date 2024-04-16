<!-- badges: start -->
[![GitHub version](https://img.shields.io/static/v1?label=GitHub&message=2.2.1&color=blue&logo=github)](https://github.com/pbreheny/plmm)
[![R-CMD-check](https://github.com/pbreheny/plmm/workflows/R-CMD-check/badge.svg)](https://github.com/pbreheny/plmm/actions)
<!-- badges: end -->

## Welcome 

The `plmmr` package contains functions that fit penalized linear mixed models to correct for unobserved confounding effects. Documentation for this package is in progress. 


## Installation 

To install the latest version of the package: 

```r
devtools::install_github("pbreheny/plmmr")
```

For a description of the motivation of the functions in this package (along with examples) refer to the second module of [this GWAS data tutorial](https://pbreheny.github.io/adv-gwas-tutorial/index.html)

## Latest changes 

The newest features of `plmmr` are: 
  - version 2.2.0: parameter $\eta$ is no longer estimated in each fold of cross-validation

  - version 2.1.0: A new function `mfdr()` for inference on model coefficients. 
  
  - version 2.0.3: An `xgboost` method is now available in `process_plink()`. Check out the documentation for details. This option should be regarded as being in 'beta-testing' mode.
  
## Note on branches 

The branches of this repo are organized in the following way: 

  - `master` is the main branch with all the latest updates

  - `integrate_fbm` is the development branch where we are working to extend our current methods to analyze data from a design matrix stored as a file-backed object (a Filebacked Big Matrix, or FBM). See package [bigstatsr](https://privefl.github.io/bigstatsr/) for more info on these objects. When we merge this into `master`, we will have `plmmr` version 3.0

  - `gh_pages` is where we are keeping all the documentation for `plmmr`
  
  - `estimate_eta` is an **archived** branch where we worked through an alternative approach for estimating $\eta$. We're keeping this around for reference as we are writing. 
  
  - `sign_flip` is an **archived** branch where we have examined the issues caused by +/- signs being flipped as part of truncated SVD.

  

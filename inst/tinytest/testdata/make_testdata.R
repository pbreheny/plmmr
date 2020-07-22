# install.packages("remotes")
# library(remotes)
# install_github("pbreheny/hdrm")
library(hdrm)
library(magrittr)
# lib_loc <- "C:/Users/areisett/Documents/R/win-library/3.6" ### pointing to previous version libraries
# to_install <- unname(installed.packages(lib.loc = lib_loc)[, "Package"])
# to_install
# install.packages(pkgs = to_install)

# should path be relative to project dir or tinytest folder?
# make sure to save in the right place...
cd <- getwd() %>%
  strsplit(., '/') %>%
  .[[1]] %>%
  .[length(.)]

if (cd == 'tinytest'){
  rt <- 'testdata'
} else if (cd == 'penalizedLMM'){
  rt <- paste0('tests/tinytest')
}

set.seed(7)
Data <- genData(5, 2, 1) # can change this to check more high-dimensional data
saveRDS(Data, file = file.path(rt, 'testdata.rds'))

### readRDS is weirdly not working on this...

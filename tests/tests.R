
library(tinytest)
# remotes::install_github('pbreheny/ncvreg')
# remotes::install_github('pbreheny/hdrm')
library(ncvreg)
library(hdrm)
library(glmnet)
library(zeallot)
for (f in list.files(path = 'R', pattern = '*.R')) source(paste0('R/', f))
setwd('tests/tinytest')
run_test_dir('.', pattern = '*.R') # is there a way to make the display output from this better?

# this test is taken from plmm_analysis_ofc/interactive/filebacked_chr10_lite.R
# Objective: troubleshoot the CV implementation -- something throws an error
#  somewhere in cvf() that made the session crash when I tried this the first time.

# read in data from previous full model fit -----------------------------------
data_dir <- "/mnt/cleft"  # adjust to fit specific machine
pheno <- readRDS(file.path(data_dir, "data/phs000774/pheno/whole/pheno_complete.rds"))
str(pheno) # see what's here

# create a binary outcome (1 = cleft diagnosis, 0 = no cleft)
pheno <- pheno |>
  dplyr::mutate(ofc_bin = dplyr::case_when(
    CleftType == 0 ~ 0, # 0 = unaffected
    CleftType > 1 ~ 1, # 1 = indeterminate phenotype , > 1 = cleft
    .default = NA_integer_
  ),
  # make ID a character (for compatibility with genotype data)
  IID = as.character(IID)) |>
  # remove missing phenotype
  dplyr::filter(!is.na(ofc_bin))

res <- paste0(file.path(data_dir, "plmm_analysis_ofc/chr10_lite/"), "std_chr10_lite")
# read file created by process_plink()
chr10_rds <- readRDS(paste0(res, ".rds"))
str(chr10_rds)
std_X <- bigmemory::attach.big.matrix(chr10_rds$std_X)
n <- nrow(std_X)
p <- ncol(std_X) # gets number of features in standardized data
chr10_fit <- readRDS(file.path(data_dir, "plmm_analysis_ofc/chr10_lite/chr10_lite_fit.rds"))

checked_data <- plmm_checks(X = paste0(file.path(data_dir,
                                                 "plmm_analysis_ofc/chr10_lite/"),
                                       "std_chr10_lite"),
                            returnX = FALSE, # KEY line: forces thsi to run filebacked
                            col_names = col_names,
                            non_genomic = chr10_rds$non_gen,
                            K = chr10_fit$K,
                            penalty = 'lasso',
                            penalty_factor =  c(0, rep(1, p - 1)),
                            gamma = 1,
                            alpha = 1,
                            trace = TRUE)

prep <- list(
  std_X = std_X,
  centered_y = checked_data$centered_y,
  K = chr10_fit$K, # Note: need this for CV (see call to construct_variance() within cv_plmm())
  s = chr10_fit$K$s,
  U = chr10_fit$K$U,
  trace = TRUE,
  eta = chr10_fit$eta # carry eta over to fit
)

# set up args for CV ------------------------------------------------------
# assign folds
nfolds <- 3
set.seed(52242)
fold <- sample(1:n) %% nfolds
fold[fold==0] <- nfolds

cv_args <- list(
  prep = prep,
  save_intermed = file.path(data_dir, "plmm_analysis_ofc/chr10_lite/filebacked_cv_intermed_res.rds"),
  seed = 52242,
  returnBiasDetails = TRUE,
  fold = fold, # for testing
  trace = TRUE,
  warn = FALSE,
  lambda = chr10_fit$lambda,
  non_genomic = chr10_rds$non_gen,
  type = "blup",
  y = checked_data$y,
  std_X_details = checked_data$std_X_details,
  penalty_factor = checked_data$penalty_factor,
  fbm_flag = checked_data$fbm_flag,
  penalty = checked_data$penalty,
  gamma = checked_data$gamma,
  alpha = 1,
  nlambda = 100,
  eps = 1e-04,
  max_iter = 10000,
  convex = TRUE,
  dfmax = NULL,
  warn = TRUE
)

# use same V as in in-memory test
estimated_V <- readRDS(paste0(file.path(data_dir, "plmm_analysis_ofc", "chr10_lite"), '/estimated_V.rds'))

# initialize objects to hold CV results
E <- Y <- matrix(NA, nrow=n, ncol=length(chr10_fit$lambda))
sde <- sqrt(.Machine$double.eps)
pb <- utils::txtProgressBar(min = 0, max = nfolds, style = 3)

# carry out CV -----------------------------------------------------------
for (i in 1:nfolds) {
  utils::setTxtProgressBar(pb, i)

  res <- cvf(i = i,
             save_intermed = file.path(data_dir, "plmm_analysis_ofc/chr10_lite/filebacked_cvf_intermed.rds"), # new; added for debugging purposes
             fold = fold,
             type = cv_args$type,
             cv_args = cv_args,
             estimated_V = estimated_V)
  utils::setTxtProgressBar(pb, i)
  close(pb)

  # update E and Y
  E[fold==i, 1:res$nl] <- res$loss

  if (!is.matrix(res$yhat)) {
    res$yhat <- as.matrix(res$yhat)
  }
  Y[fold==i, 1:res$nl] <- res$yhat

  # save results along the way
  val_intermed <- list(type=cv_args$type,
                       res = res,
                       loss = E,
                       yhat = Y,
                       current_fold = i)
  saveRDS(val_intermed, cv_args$save_intermed)

}

# June 30, 2024 - that runs, in all 3 folds!
# post-process results -----------------------------------------
## function to post-processs ---------------------------------------------------

postprocess_cvf_res <- function(rds){
  # eliminate saturated lambda values, if any
  ind <- which(apply(is.finite(rds$loss), 2, all)) # index for lambda values to keep
  rds$loss <- rds$loss[, ind, drop=FALSE]
  rds$yhat <- rds$yhat[, ind]
  lambda <- cv_args$lambda[ind]

  # return min lambda idx
  cve <- apply(rds$loss, 2, mean)
  cvse <- apply(rds$loss, 2, stats::sd) / sqrt(nrow(rds$yhat))
  min <- which.min(cve)

  # return lambda 1se idx
  l.se <- cve[min] - cvse[min]
  u.se <- cve[min] + cvse[min]
  within1se <- which(cve >= l.se & cve <= u.se)
  min1se <- which.max(lambda %in% lambda[within1se])
  # bias correction
  e <- sapply(1:nfolds, function(i) apply(rds$loss[fold==i, , drop=FALSE], 2, mean))
  Bias <- mean(e[min,] - apply(e, 2, min))

  cv_results <- structure(
    list(type=cv_args$type,
         cve=cve,
         cvse=cvse,
         fold=fold,
         lambda=lambda,
         fit=chr10_fit,
         min=min,
         lambda_min=lambda[min],
         min1se = min1se,
         lambda.1se = lambda[min1se],
         null.dev=mean(plmm_loss(checked_data$y, rep(mean(checked_data$y), n)))),
    class="cv_plmm")

  return(cv_results)
}

rds <- readRDS("/mnt/cleft/plmm_analysis_ofc/chr10_lite/filebacked_cv_intermed_res.rds")
fb_res <- postprocess_cvf_res(rds)
# look at results ---------------------------------------------------------------
plot(fb_res)
summary(fb_res)
best_betas <- chr10_fit$beta_vals[,fb_res$min]
best_betas[abs(best_betas) > 1e-4]
coef(fb_res, lambda = fb_res$lambda_min)[abs(best_betas > 1e-4)]

# compare to in-memory fit -----------------------------------------------------
in_mem_fit <- readRDS('/mnt/cleft/plmm_analysis_ofc/chr10_lite/cv_intermed_res.rds')
inmem_res <- postprocess_cvf_res(in_mem_fit)
summary(inmem_res)
# different (but very similar) lambdas are chosen
fb_res$lambda_min; inmem_res$lambda_min
fb_res$min; inmem_res$min

inmem_betas <- chr10_fit$beta_vals[,inmem_res$min]
inmem_betas[abs(inmem_betas) > 1e-4]

# look at two chosen sets of coefficients, side by side:
bb <- chr10_fit$beta_vals[,fb_res$min:inmem_res$min]
bb[abs(bb[,1]) > 1e-4 | abs(bb[,2]) > 1e-4,]

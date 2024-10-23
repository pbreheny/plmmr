# plmmr 4.1.0 (2024-10-23)

- **Restore plmm(X,y) syntax**: Where version 4.0.0 required that `create_design()` always be called prior to `plmm()` or `cv_plmm()`; this update restores the X,y syntax consistent with other packages (e.g., `glmnet`, `ncvreg`). Note that this syntax is only available for the case where the design matrix is stored in-memory as a `matrix` or `data.frame` object. The `create_design()` function is still required for cases where the design matrix/dataset is stored in an external file. 

- **Bug fix**: The 4.0.0 version of `create_design()` required `X` to have column names, and errored out with an uninformative message if no names were supplied (see issue 61). This is now fixed -- column names are not required unless the user wants to specify an argument to `unpen`. 

- **Argument name change**: In `create_design()`, the argument to specify an outcome in the in-memory case has been renamed to `y`; this makes the syntax consistent, e.g., `create_design(X, y)`. Note again that this change is relevant to in-memory data only. 

* **Internal:** Fixed LTO type mismatch bug.

# plmmr 4.0.0 (2024-10-07)

- **Major re-structuring of preprocessing pipeline:** Data from external files must now be processed with `process_plink()` or `process_delim()`. All data (including in-memory data) must be prepared for analysis via `create_design()`. This change ensures that data are funneled into a uniform format for analysis.

- **Documentation updated:** The vignettes for the package are now all revised to include examples of the complete pipeline with the new `create_design()` syntax. There is an article for each type of data input (matrix/data.frame, delimited file, and PLINK).

- **CRAN:** The package is on CRAN now.

# plmmr 3.2.0 (2024-09-02)

- **bigsnpr now in Suggests, not Imports:** The essential filebacking support is now all done with `bigmemory` and `bigalgebra`. The `bigsnpr` package is used only for processing PLINK files. 

- **dev branch gwas_scale** has a version of the pipeline that runs completely file-backed. 

# plmmr 3.1.0 (2024-07-13)

- **Enhancement:** To make `plmmr` have better functionality for writing scripts, the functions `process_plink()`, `plmmm()`, and `cv_plmm()` now (optionally) write '.log' files, as in PLINK.

- **Enhancement:** In cases where users are working with large datasets, it may not be practical or desirable for all the results returned by `plmmm()` or `cv_plmm()` to be saved in a single '.rds' file. There is now an option in both of these model fitting functions called 'compact_save', which gives users the option to save the output in multiple, smaller '.rds' files.

- **Argument removed:** Argument `std_needed` is no longer available in `plmm()` and `cv_plmm()` functions.

# plmmr 3.0.0 (2024-06-27)

- **Bug fix:** Cross-validation implementation issues fixed. Previously, the full set of eigenvalues were used inside CV folds, which is not ideal as it involves information from outside the fold. Now, the entire modeling process is cross-validated: the standardization, the eigendecomposition of the relatedness matrix, the model fitting, and the backtransformation onto the original scale for prediction.

- **Computational speedup:** The standardization and rotation of filebacked data are now much faster; `bigalgebra` and `bigmemory` are now used for these computations.

- **Internal:** On the standardized scale, the intercept of the PLMM is the mean of the outcome. This derivation considerably simplifies the handling of the intercept internally during model fitting.

# plmmr 2.2.1 (2024-03-16)

- **Name change:** Changed package name to `plmmr`; note that `plmm()`, `cv_plmm()`, and other functions starting with `plmm_` have not changed names.

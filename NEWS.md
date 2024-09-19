# plmmr 3.2.0.0

## Most recent changes

## 3.2.0.0

These changes were borne out of the troubleshooting of the real data analysis elements of my dissertation. This is the version that has a integrated system for analyzing data from in-memory matrices/data frames, delimited files, and PLINK files. Checkout the highlights from these changes below: 

-   **Major re-structuring of preprocessing pipeline**: Data from external files must now be processed with `process_plink()` or `process_delim()`. All data (including in-memory data) must be prepared for analysis via `create_design()`. This change ensures that data are funneled into a uniform format for analysis. 

-  **bigsnpr now in Suggests, not Imports**: The essential filebacking support is now all done with `bigmemory` and `bigalgebra`. The `bigsnpr` package is used only for processing PLINK files. 

## 3.1.0

-   **Enhancement**: To make `plmmr` have better functionality for writing scripts, the functions `process_plink()`, `plmmm()`, and `cv_plmm()` now write '.log' files each time they are called. The design of these log files is based on the behavior of [PLINK software](https://www.cog-genomics.org/plink/).

-   **Enhancement**: In cases where users are working with large datasets, it may not be practical or desirable for all the results returned by `plmmm()` or `cv_plmm()` to be saved in a single '.rds' file. There is now an option in both of these model fitting functions called 'compact_save', which gives users the option to save the output in multiple, smaller '.rds' files.

-   **Argument removed**: Previously, the functions `plmmm()` and `cv_plmm()` had an argument 'std_needed' that allowed users to 'turn off' standardization. The more I have worked with this and thought about it, the more I realized that the package isn't mature enough yet to handle this complex, nuanced option. This option has been removed for now -- we may add this back in later.

## 3.0.0

-   **Bug fix**: As of June 27, 2024, we have addressed some issues with our cross-validation implementation. Previously, we were using all eigenvalues and the estimated $\hat\eta$ in each CV fold -- but this is not consistent with the best practices for CV implementation, as information from outside a given fold should not inform predictions. These best practices are summarized in the *Elements of Statistical Learning* by Hastie et al., section 7.10 (in 2nd edition, available [online](https://hastie.su.domains/Papers/ESLII.pdf) from author).

As of this update, all modeling steps are carried out in each fold: the standardization, the eigendecomposition of the relatedness matrix, the model fitting, and the backtransformation onto the original scale for prediction. There may be a way to make the eigendecomposition step faster -- this is a question we are actively studying.

-   **Computational speedup**: The standardization and rotation of filebacked data are now much faster; we have moved toward using methods from `bigalgebra` and `bigmemory` for these computations.

-   **Methods development** (for the nerds): We have derived that on the standardized scale, the intercept of our PLMM is the mean of the outcome. With this in mind, we have updated the way we handle the intercept in model fitting -- which has made our code in internal functions more 'readable'.


# plmmr 2.2.1

-   **Notable changes**

    -   Changed package name to `plmmr` - note that `plmm()`, `cv_plmm()`, and other functions starting with `plmm_` have not changed names.

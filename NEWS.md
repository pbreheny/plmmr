# plmmr 3.1.0

## Most recent changes

## 3.1.0

-   **Enhancement**: To make `plmmr` have better functionality for writing scripts, the functions `process_plink()`, `plmmm()`, and `cv_plmm()` now write '.log' files each time they are called. The design of these log files is based on the behavior of [PLINK software](https://www.cog-genomics.org/plink/).

-   **Enhancement**: In cases where users are working with large datasets, it may not be practical or desirable for all the results returned by `plmmm()` or `cv_plmm()` to be saved in a single '.rds' file. There is now an option in both of these model fitting functions called 'compact_save', which gives users the option to save the output in multiple, smaller '.rds' files.

-   **Argument removed**: Previously, the functions `plmmm()` and `cv_plmm()` had an argument 'std_needed' that allowed users to 'turn off' standardization. The more I have worked with this and thought about it, the more I realized that the package isn't mature enough yet to handle this complex, nuanced option. This option has been removed for now -- we may add this back in later.

## 3.0.0

-   **Bug fix**: As of June 27, 2024, we have addressed some issues with our cross-validation implementation. Previously, we were using all eigenvalues and the estimated $\hat\eta$ in each CV fold -- but this is not consistent with the best practices for CV implementation, as information from outside a given fold should not inform predictions.[^news-1] As of this update, all modeling steps are carried out in each fold: the standardization, the eigendecomposition of the relatedness matrix, the model fitting, and the backtransformation onto the original scale for prediction. There may be a way to make the eigendecomposition step faster -- this is a question we are actively studying.

-   **Computational speedup**: The standardization and rotation of filebacked data are now much faster; we have moved toward using methods from `bigalgebra` and `bigmemory` for these computations.

-   **Methods development** (for the nerds): We have derived that on the standardized scale, the intercept of our PLMM is the mean of the outcome. With this in mind, we have updated the way we handle the intercept in model fitting -- which has made our code in internal functions more 'readable'.

[^news-1]: See *Elements of Statistical Learning* by Hastie et al., section 7.10 (in 2nd edition). Available [online](https://hastie.su.domains/Papers/ESLII.pdf) from author.

# plmmr 2.2.1

-   **Notable changes**

    -   Changed package name to `plmmr` - note that `plmm()`, `cv_plmm()`, and other functions starting with `plmm_` have not changed names.

# plmmr 3.0.0

## Most recent changes
  - **Bug fix**: As of June 27, 2024, we have addressed some issues with our cross-validation implementation. Previously, we were using all eigenvalues and the estimated $\hat\eta$ in each CV fold -- but this is not consistent with the best practices for CV implementation, as information from outside a given fold should not inform predictions.[^news-1] As of this update, all modeling steps are carried out in each fold: the standardization, the eigendecomposition of the relatedness matrix, the model fitting, and the backtransformation onto the original scale for prediction. There may be a way to make the eigendecompsition step faster -- this is a question we are actively studying. 
    
  - **Computational speedup**: The standardization and rotation of filebacked data are now much faster; we have moved toward using methods from `bigalgebra` and `bigmemory` for these computations. 
    
  - **Methods developement** (for the nerds): We have derived that on the standardized scale, the intercept of our PLMM is the mean of the outcome. With this in mind, we have updated the way we handle the intercept in model fitting -- which has made our code in internal functions more 'readable'. 

  - **CRAN**. This is soon-to-be our initial CRAN submission!

[^news-1]: See *Elements of Statistical Learning* by Hastie et al., section 7.10 (in 2nd edition). Available [online](https://hastie.su.domains/Papers/ESLII.pdf) from author.

# plmmr 2.2.1

-   **Notable changes**

    -   Changed package name to `plmmr` - note that `plmm()`, `cv_plmm()`, and other functions starting with `plmm_` have not changed names.

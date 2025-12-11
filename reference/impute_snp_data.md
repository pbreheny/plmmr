# A function to impute SNP data

A function to impute SNP data

## Usage

``` r
impute_snp_data(
  obj,
  X,
  impute,
  impute_method,
  parallel,
  outfile,
  quiet,
  seed = as.numeric(Sys.Date()),
  ...
)
```

## Arguments

- obj:

  a `bigSNP` object (as created by
  [`read_plink_files()`](https://pbreheny.github.io/plmmr/reference/read_plink_files.md))

- X:

  A matrix of genotype data as returned by `name_and_count_bigsnp`

- impute:

  Logical: should data be imputed? Default to TRUE.

- impute_method:

  If 'impute' = TRUE, this argument will specify the kind of imputation
  desired. Options are:

  - mode (default): Imputes the most frequent call. See
    [`bigsnpr::snp_fastImputeSimple()`](https://privefl.github.io/bigsnpr/reference/snp_fastImputeSimple.html)
    for details.

  - random: Imputes sampling according to allele frequencies.

  - mean0: Imputes the rounded mean.

  - mean2: Imputes the mean rounded to 2 decimal places.

  - xgboost: Imputes using an algorithm based on local XGBoost models.
    See
    [`bigsnpr::snp_fastImpute()`](https://privefl.github.io/bigsnpr/reference/snp_fastImpute.html)
    for details. Note: this can take several minutes, even for a
    relatively small data set.

- parallel:

  Logical: should the computations within this function be run in
  parallel? Defaults to TRUE. See
  [`count_cores()`](https://pbreheny.github.io/plmmr/reference/count_cores.md)
  and
  [`?bigparallelr::assert_cores`](https://rdrr.io/pkg/bigparallelr/man/assert_cores.html)
  for more details. In particular, the user should be aware that too
  much parallelization can make computations *slower*.

- outfile:

  Optional: the name (character string) of the prefix of the logfile to
  be written. Defaults to 'process_plink', i.e. you will get
  'process_plink.log' as the outfile.

- quiet:

  Logical: should messages be printed to the console? Defaults to TRUE

- seed:

  Numeric value to be passed as the seed for
  `impute_method = 'xgboost'`. Defaults to `as.numeric(Sys.Date())`

- ...:

  Optional: additional arguments to
  [`bigsnpr::snp_fastImpute()`](https://privefl.github.io/bigsnpr/reference/snp_fastImpute.html)
  (relevant only if impute_method = "xgboost")

## Value

Nothing is returned, but the `obj$genotypes` is overwritten with the
imputed version of the data

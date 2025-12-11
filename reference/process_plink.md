# Preprocess PLINK files using the `bigsnpr` package

Preprocess PLINK files using the `bigsnpr` package

## Usage

``` r
process_plink(
  data_dir,
  data_prefix,
  rds_dir = data_dir,
  rds_prefix = NULL,
  logfile = NULL,
  impute = TRUE,
  impute_method = "mode",
  id_var = "IID",
  parallel = TRUE,
  quiet = FALSE,
  overwrite = FALSE,
  ...
)
```

## Arguments

- data_dir:

  The path to the bed/bim/fam data files, *without* a trailing "/"
  (e.g., use `data_dir = '~/my_dir'`, **not** `data_dir = '~/my_dir/'`)

- data_prefix:

  The prefix (as a character string) of the bed/fam data files (e.g.,
  `data_prefix = 'mydata'`)

- rds_dir:

  The path to the directory in which you want to create the new '.rds'
  and '.bk' files. Defaults to `data_dir`

- rds_prefix:

  String specifying the user's preferred filename for the to-be-created
  .rds file (will be create inside `rds_dir` folder). If no rds_prefix
  is provided, the processed data files will be returned in memory.
  Note: 'rds_prefix' cannot be the same as 'data_prefix'

- logfile:

  Optional: the name (character string) of the prefix of the logfile to
  be written in 'rds_dir'. Default to NULL (no log file written). Note:
  if you supply a file path in this argument, it will error out with a
  "file not found" error. Only supply the string; e.g., if you want
  my_log.log, supply 'my_log', the my_log.log file will appear in
  rds_dir.

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

- id_var:

  String specifying which column of the PLINK `.fam` file has the unique
  sample identifiers. Options are "IID" (default) and "FID"

- parallel:

  Logical: should the computations within this function be run in
  parallel? Defaults to TRUE. See
  [`count_cores()`](https://pbreheny.github.io/plmmr/reference/count_cores.md)
  and
  [`?bigparallelr::assert_cores`](https://rdrr.io/pkg/bigparallelr/man/assert_cores.html)
  for more details. In particular, the user should be aware that too
  much parallelization can make computations *slower*.

- quiet:

  Logical: should messages to be printed to the console be silenced?
  Defaults to FALSE

- overwrite:

  Logical: if existing `.bk`/`.rds` files exist for the specified
  directory/prefix, should these be overwritten? Defaults to FALSE. Set
  to TRUE if you want to change the imputation method you're using, etc.
  **Note**: If there are multiple `.rds` files with names that start
  with "std_prefix\_...", **this will error out**. To protect users from
  accidentally deleting files with saved results, only one `.rds` file
  can be removed with this option.

- ...:

  Optional: additional arguments to
  [`bigsnpr::snp_fastImpute()`](https://privefl.github.io/bigsnpr/reference/snp_fastImpute.html)
  (relevant only if impute_method = "xgboost")

## Value

The filepath to the '.rds' object created; see details for explanation.

## Details

Three files are created in the location specified by `rds_dir`:

- 'rds_prefix.rds': This is a list with three items: (1) `X`: the
  filebacked
  [`bigmemory::big.matrix`](https://rdrr.io/pkg/bigmemory/man/big.matrix.html)
  object pointing to the imputed genotype data. This matrix has type
  'double', which is important for downstream operations in
  [`create_design()`](https://pbreheny.github.io/plmmr/reference/create_design.md) (2)
  `map`: a data.frame with the PLINK 'bim' data (i.e., the variant
  information) (3) `fam`: a data.frame with the PLINK 'fam' data (i.e.,
  the pedigree information)

- 'prefix.bk': This is the backingfile that stores the numeric data of
  the genotype matrix

- 'rds_prefix.desc'" This is the description file, as needed by the

Note that `process_plink()` need only be run once for a given set of
PLINK files; in subsequent data analysis/scripts,
[`get_data()`](https://pbreheny.github.io/plmmr/reference/get_data.md)
will access the '.rds' file.

For an example, see vignette on processing PLINK files

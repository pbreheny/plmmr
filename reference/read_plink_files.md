# A function to read in PLINK files using `bigsnpr` methods

A function to read in PLINK files using `bigsnpr` methods

## Usage

``` r
read_plink_files(
  data_dir,
  data_prefix,
  rds_dir,
  rds_prefix,
  outfile,
  parallel,
  overwrite,
  quiet
)
```

## Arguments

- data_dir:

  The path to the bed/bim/fam data files, *without* a trailing "/"
  (e.g., use `data_dir = '~/my_dir'`, **not** `data_dir = '~/my_dir/'`)

- data_prefix:

  The prefix (as a character string) of the bed/fam data files (e.g.,
  `prefix = 'mydata'`)

- rds_dir:

  The path to the directory in which you want to create the new `.rds`
  and `.bk` files. Defaults to `data_dir`

- rds_prefix:

  String specifying the user's preferred filename for the to-be-created
  `.rds` file (will be create inside `rds_dir` folder). If no rds_prefix
  is provided, the processed data files will be returned in memory.
  Note: `rds_prefix` cannot be the same as `data_prefix`

- outfile:

  Optional: the name (character string) of the prefix of the logfile to
  be written. Defaults to NULL (no log written).

- parallel:

  Logical: should the computations within this function be run in
  parallel? Defaults to TRUE. See
  [`count_cores()`](https://pbreheny.github.io/plmmr/reference/count_cores.md)
  and
  [`?bigparallelr::assert_cores`](https://rdrr.io/pkg/bigparallelr/man/assert_cores.html)
  for more details. In particular, the user should be aware that too
  much parallelization can make computations *slower*.

- overwrite:

  Logical: if existing `.bk`/`.rds` files exist for the specified
  directory/prefix, should these be overwritten? Defaults to FALSE. Set
  to TRUE if you want to change the imputation method you're using, etc.

- quiet:

  Logical: should messages be printed to the console? Defaults to TRUE

## Value

`.rds` and `.bk` files are created in `data_dir`, and `obj` (a `bigSNP`
object) is returned. See `bigsnpr` documentation for more info on the
`bigSNP` class.

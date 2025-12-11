# A helper function to standardize a filebacked matrix

A helper function to standardize a filebacked matrix

## Usage

``` r
standardize_filebacked(
  X,
  new_file,
  rds_dir,
  non_gen,
  complete_outcome,
  id_var,
  outfile,
  quiet,
  overwrite,
  tocenter = TRUE
)
```

## Arguments

- X:

  A list that includes: (1) subset_X: a `big.matrix` object that has
  been subset &/or had any additional predictors appended as columns (2)
  ns: a numeric vector indicating the indices of nonsingular columns in
  subset_X

- new_file:

  The new_file (as a character string) of the bed/fam data files (e.g.,
  `new_file = 'mydata'`)

- rds_dir:

  The path to the directory in which you want to create the new '.rds'
  and '.bk' files. Defaults to `data_dir`

- outfile:

  Optional: the name (character string) of the new_file of the logfile
  to be written. Defaults to 'process_plink', i.e. you will get
  'process_plink.log' as the outfile.

- quiet:

  Logical: should messages be printed to the console? Defaults to FALSE
  (which leaves the print messages on...)

- overwrite:

  Logical: if existing `.bk`/`.rds` files exist for the specified
  directory/new_file, should these be overwritten?

## Value

A list with a new component of `obj` called 'std_X' - this is an FBM
with column-standardized data. List also includes several other
indices/meta-data on the standardized matrix

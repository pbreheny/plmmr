# A helper function to subset `big.matrix` objects

A helper function to subset `big.matrix` objects

## Usage

``` r
subset_filebacked(X, new_file, complete_samples, ns, rds_dir, outfile, quiet)
```

## Arguments

- X:

  A filebacked `big.matrix` with the to-be-standardized design matrix

- new_file:

  Optional user-specified new_file for the to-be-created .rds/.bk files.

- complete_samples:

  Numeric vector with indicesmarking the rows of the original data which
  have a non-missing entry in the 6th column of the `.fam` file

- ns:

  Numeric vector with the indices of the non-singular columns This
  vector is created in `handle_missingness()`

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

## Value

A list with two components. First, a `big.matrix` object, 'subset_X',
representing a design matrix wherein:

- rows are subset according to user's specification in
  `handle_missing_phen`

- columns are subset so that no constant features remain – this is
  important for standardization downstream The list also includes the
  integer vector 'ns' which marks which columns of the original matrix
  were 'non-singular' (i.e. *not* constant features). The 'ns' index
  plays an important role in
  [`plmm_format()`](https://pbreheny.github.io/plmmr/reference/plmm_format.md)
  and
  [`untransform()`](https://pbreheny.github.io/plmmr/reference/untransform.md)
  (both helper functions in model fitting)

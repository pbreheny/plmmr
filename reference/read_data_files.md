# A function to read in a large file as a numeric file-backed matrix (`FBM`) Note: this function is a wrapper for `bigstatsr::big_read()`

A function to read in a large file as a numeric file-backed matrix
(`FBM`) Note: this function is a wrapper for
[`bigstatsr::big_read()`](https://privefl.github.io/bigstatsr/reference/big_read.html)

## Usage

``` r
read_data_files(
  data_file,
  data_dir,
  rds_dir,
  rds_prefix,
  outfile,
  overwrite,
  quiet,
  ...
)
```

## Arguments

- data_file:

  The name of the file to read, not including its directory. Directory
  should be specified in `data_dir`

- data_dir:

  The path to the directory where 'file' is

- rds_dir:

  The path to the directory in which you want to create the new '.rds'
  and '.bk' files. Defaults to `data_dir`

- rds_prefix:

  String specifying the user's preferred filename for the to-be-created
  .rds/.bk files (will be create insie `rds_dir` folder) Note:
  'rds_prefix' cannot be the same as 'data_file'

- outfile:

  Optional: the name (character string) of the prefix of the logfile to
  be written. Defaults to 'process_plink', i.e. you will get
  'process_plink.log' as the outfile.

- overwrite:

  Logical: if existing `.bk`/`.rds` files exist for the specified
  directory/prefix, should these be overwritten? Defaults to FALSE. Set
  to TRUE if you want to change the imputation method you're using, etc.

- quiet:

  Logical: should messages be printed to the console? Defaults to TRUE

- ...:

  Optional: other arguments to be passed to
  [`bigmemory::read.big.matrix()`](https://rdrr.io/pkg/bigmemory/man/write.big.matrix.html).
  Note: 'sep' is an option to pass here.

## Value

'.rds', '.bk', and '.desc' files are created in `data_dir`, and `obj` (a
filebacked `bigmemory big.matrix` object) is returned. See `bigmemory`
documentation for more info on the `big.matrix` class.

# A function to read in large data files as an FBM

A function to read in large data files as an FBM

## Usage

``` r
process_delim(
  data_dir,
  data_file,
  feature_id,
  rds_dir = data_dir,
  rds_prefix,
  logfile = NULL,
  overwrite = FALSE,
  quiet = FALSE,
  ...
)
```

## Arguments

- data_dir:

  The directory to the file.

- data_file:

  The file to be read in, without the filepath. This should be a file of
  numeric values. Example: use `data_file = "myfile.txt"`, not
  `data_file = "~/mydirectory/myfile.txt"` Note: if your file has
  headers/column names, set 'header = TRUE' – this will be passed into
  [`bigmemory::read.big.matrix()`](https://rdrr.io/pkg/bigmemory/man/write.big.matrix.html).

- feature_id:

  A string specifying the column in the data X (the feature data) with
  the row IDs (e.g., identifiers for each row/sample/participant/,
  etc.). No duplicates allowed.

- rds_dir:

  The directory where the user wants to create the '.rds' and '.bk'
  files Defaults to `data_dir`

- rds_prefix:

  String specifying the user's preferred filename for the to-be-created
  .rds file (will be create inside `rds_dir` folder) Note: 'rds_prefix'
  cannot be the same as 'data_prefix'

- logfile:

  Optional: the name (character string) of the prefix of the logfile to
  be written. Defaults to 'process_delim', i.e. you will get
  'process_delim.log' as the outfile.

- overwrite:

  Optional: the name (character string) of the prefix of the logfile to
  be written. Defaults to 'process_plink', i.e. you will get
  'process_plink.log' as the outfile. **Note**: If there are multiple
  `.rds` files with names that start with "std_prefix\_...", **this will
  error out**. To protect users from accidentally deleting files with
  saved results, only one `.rds` file can be removed with this option.

- quiet:

  Logical: should the messages printed to the console be silenced?
  Defaults to FALSE.

- ...:

  Optional: other arguments to be passed to
  [`bigmemory::read.big.matrix()`](https://rdrr.io/pkg/bigmemory/man/write.big.matrix.html).
  Note: 'sep' is an option to pass here, as is 'header'.

## Value

The file path to the newly created '.rds' file

## Examples

``` r
temp_dir <- tempdir()
colon_dat <- process_delim(data_file = "colon2.txt",
 data_dir = find_example_data(parent = TRUE), overwrite = TRUE,
 rds_dir = temp_dir, rds_prefix = "processed_colon2", sep = "\t", header = TRUE)
#> 
#> Overwriting existing files:processed_colon2.bk/.rds/.desc
#> There are 62 observations and 2001 features in the specified data files.
#> At this time, plmmr::process_delim() does not not handle missing values in delimited data.
#>       Please make sure you have addressed missingness before you proceed.
#> 
#> process_plink() completed 
#> Processed files now saved as /tmp/RtmpLOJR73/processed_colon2.rds

colon2 <- readRDS(colon_dat)
str(colon2)
#> List of 3
#>  $ X:Formal class 'big.matrix.descriptor' [package "bigmemory"] with 1 slot
#>   .. ..@ description:List of 13
#>   .. .. ..$ sharedType: chr "FileBacked"
#>   .. .. ..$ filename  : chr "processed_colon2.bk"
#>   .. .. ..$ dirname   : chr "/tmp/RtmpLOJR73/"
#>   .. .. ..$ totalRows : int 62
#>   .. .. ..$ totalCols : int 2001
#>   .. .. ..$ rowOffset : num [1:2] 0 62
#>   .. .. ..$ colOffset : num [1:2] 0 2001
#>   .. .. ..$ nrow      : num 62
#>   .. .. ..$ ncol      : num 2001
#>   .. .. ..$ rowNames  : NULL
#>   .. .. ..$ colNames  : chr [1:2001] "sex" "Hsa.3004" "Hsa.13491" "Hsa.13491.1" ...
#>   .. .. ..$ type      : chr "double"
#>   .. .. ..$ separated : logi FALSE
#>  $ n: num 62
#>  $ p: num 2001
#>  - attr(*, "class")= chr "processed_delim"
```

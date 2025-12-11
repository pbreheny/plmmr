# A function to help with accessing example PLINK files

A function to help with accessing example PLINK files

## Usage

``` r
find_example_data(path, parent = FALSE)
```

## Arguments

- path:

  Argument (string) specifying a path (filename) for an external data
  file in `extdata/`

- parent:

  If `path=TRUE` and the user wants the name of the parent directory
  where that file is located, set `parent=TRUE`. Defaults to FALSE.

## Value

If `path=NULL`, a character vector of file names is returned. If path is
given, then a character string with the full file path

## Examples

``` r
find_example_data(parent = TRUE)
#> [1] "/home/runner/work/_temp/Library/plmmr/extdata"
```

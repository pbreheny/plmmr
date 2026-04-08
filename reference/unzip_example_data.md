# Companion function to unzip the .gz files that ship with the `plmmr` package.

Companion function to unzip the .gz files that ship with the `plmmr`
package.

## Usage

``` r
unzip_example_data(outdir)
```

## Arguments

- outdir:

  The file path to the directory to which the .gz files should be
  written.

## Value

Nothing is returned; the PLINK files that ship with the `plmmr` package
are stored in the directory specified by 'outdir'.

## Details

For an example of this function, look at
`vignette('plink_files', package = "plmmr")`.

# A helper function to label and summarize the contents of a `bigSNP`

A helper function to label and summarize the contents of a `bigSNP`

## Usage

``` r
name_and_count_bigsnp(obj, id_var, quiet, outfile)
```

## Arguments

- obj:

  a `bigSNP` object, possibly subset by `add_external_phenotype()`

- id_var:

  String specifying which column of the PLINK `.fam` file has the unique
  sample identifiers. Options are "IID" (default) and "FID".

- quiet:

  Logical: should messages be printed to the console? Defaults to TRUE

- outfile:

  The string with the name of the .log file

## Value

a list with components:

- counts: column-wise summary of the minor allele counts in 'genotypes'

- obj: a modified `bigSNP` list with additional components

- X: the `obj$genotypes` as its own FBM

- pos: the `obj$map$physical.pos` vector

# A helper function to label and summarize the contents of a `bigSNP`

A helper function to label and summarize the contents of a `bigSNP`

## Usage

``` r
name_and_count_bigsnp(obj, id_var, outfile, quiet)
```

## Arguments

- obj:

  a `bigSNP` object, possibly subset by `add_external_phenotype()`

- id_var:

  String specifying which column of the PLINK `.fam` file has the unique
  sample identifiers. Options are "IID" (default) and "FID".

- outfile:

  The string with the name of the `.log` file

- quiet:

  Logical: should console messages be silenced? Defaults to FALSE

## Value

a list with 7 components:

- `na_counts`: vector of missing SNP counts in `genotypes`

- `obj`: a modified `bigSNP` list with additional components

- `og_plink_ids`: either the IID or FID column from `.fam`, determined
  by `id_var`

- `chr`: p-length containing the chromosomes for each SNP

- `X`: the `obj$genotypes` as its own FBM

- `pos`: vector of physical positions of the SNPs

- `chr_range`: vector containing the minimum and maximum values of
  `chr`. Character strings are treated as the **maximum**.

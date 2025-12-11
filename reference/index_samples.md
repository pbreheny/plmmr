# A function to align genotype and phenotype data

A function to align genotype and phenotype data

## Usage

``` r
index_samples(
  obj,
  rds_dir,
  indiv_id,
  add_outcome,
  outcome_id,
  outcome_col,
  na_outcome_vals,
  outfile,
  quiet
)
```

## Arguments

- obj:

  An object created by
  [`process_plink()`](https://pbreheny.github.io/plmmr/reference/process_plink.md)

- rds_dir:

  The path to the directory in which you want to create the new '.rds'
  and '.bk' files.

- indiv_id:

  A character string indicating the ID column name in the 'fam' element
  of the genotype data list. Defaults to 'sample.ID', equivalent to
  'IID' in PLINK. The other option is 'family.ID', equivalent to 'FID'
  in PLINK.

- add_outcome:

  A data frame with at least two columns: and ID column and a phenotype
  column

- outcome_id:

  A string specifying the name of the ID column in `pheno`

- outcome_col:

  A string specifying the name of the phenotype column in `pheno`. This
  column will be used as the default `y` argument to 'plmm()'.

- na_outcome_vals:

  A vector of numeric values used to code NA values in the outcome.
  Defaults to `c(-9, NA_integer)` (the -9 matches PLINK conventions).

- outfile:

  A string with the name of the filepath for the log file

- quiet:

  Logical: should messages be printed to the console? Defaults to FALSE
  (which leaves the print messages on...

## Value

a list with two items:

- a data.table with rows corresponding to the samples for which both
  genotype and phenotype are available.

- a numeric vector with indices indicating which samples were 'complete'
  (i.e., which samples from add_outcome had corresponding data in the
  PLINK files)

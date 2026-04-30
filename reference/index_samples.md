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

  The path to the directory in which you want to create the new `.rds`
  and `.bk` files.

- indiv_id:

  A character string indicating the ID column name in the 'fam' element
  of the genotype data list. Defaults to 'sample.ID', equivalent to
  'IID' in PLINK. The other option is 'family.ID', equivalent to 'FID'
  in PLINK.

- add_outcome:

  A data frame with at least two columns: an ID column and a phenotype
  column

- outcome_id:

  A string specifying the name of the ID column in `add_outcome`

- outcome_col:

  A string specifying the name of the phenotype column in `add_outcome`.
  This column will be used as the default `y` argument to
  [`plmm()`](https://pbreheny.github.io/plmmr/reference/plmm.md).

- na_outcome_vals:

  A vector of numeric values used to code NA values in the outcome.
  Defaults to `c(-9, NA_integer)` (the -9 matches PLINK conventions).

- outfile:

  A string with the name of the filepath for the log file

- quiet:

  Logical: should console messages be silenced? Defaults to FALSE

## Value

a list with two items:

- `complete_samples`: a data.table with rows corresponding to the
  samples for which both genotype and phenotype are available.

- `outcome_idx`: a numeric vector with indices indicating which samples
  were 'complete' (i.e., which samples from add_outcome had
  corresponding data in the PLINK files)

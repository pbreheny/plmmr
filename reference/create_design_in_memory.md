# A function to create a design with an in-memory X matrix

A function to create a design with an in-memory X matrix

## Usage

``` r
create_design_in_memory(X, y, unpen = NULL)
```

## Arguments

- X:

  A numeric matrix in which rows correspond to observations (e.g.,
  samples) and columns correspond to features.

- y:

  A numeric vector representing the outcome for the model. **Note**: it
  is the responsibility of the user to ensure that the outcome_col and X
  have the same row order!

- unpen:

  An optional character vector with the names of columns to mark as
  unpenalized (i.e., these features would always be included in a
  model). **Note**: if you choose to use this option, X must have column
  names.

## Value

A list with elements including a standardized X and model design
information

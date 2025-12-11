# helper function to implement SCAD penalty

helper function to implement SCAD penalty

## Usage

``` r
SCAD(z, l1, l2, gamma, v)
```

## Arguments

- z:

  solution over active set at each feature

- l1:

  upper bound

- l2:

  lower bound

- gamma:

  The tuning parameter of the SCAD penalty

- v:

  the 'xtx' term

## Value

numeric vector of the SCAD-penalized coefficient estimates within the
given bounds

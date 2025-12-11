# helper function to implement MCP penalty The helper functions to implement each penalty.

helper function to implement MCP penalty The helper functions to
implement each penalty.

## Usage

``` r
MCP(z, l1, l2, gamma, v)
```

## Arguments

- z:

  a vector representing the solution over active set at each feature

- l1:

  upper bound (on beta)

- l2:

  lower bound (on beta)

- gamma:

  The tuning parameter of the MCP penalty

- v:

  the 'xtx' term

## Value

numeric vector of the MCP-penalized coefficient estimates within the
given bounds

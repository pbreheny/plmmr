# A helper function to count the number of cores available on the current machine

A helper function to count the number of cores available on the current
machine

## Usage

``` r
count_cores()
```

## Value

A number of cores to use; if `parallel` is installed, this will be
[`parallel::detectCores()`](https://rdrr.io/r/parallel/detectCores.html).
Otherwise, this returns a 1.

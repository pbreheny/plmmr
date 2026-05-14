# A function to rotate filebacked data

A function to rotate filebacked data

## Usage

``` r
rotate_filebacked(prep, y, tocenter = TRUE, ...)
```

## Arguments

- prep:

  The object returned by
  [`plmm_prep()`](https://pbreheny.github.io/plmmr/reference/plmm_prep.md)

- y:

  The continuous outcome vector.

- tocenter:

  Should the matrix be centered in addition to scaled? Defaults to TRUE

- ...:

  Not used

## Value

a list with 4 items:

- `stdrot_X`: `X` on the rotated and re-standardized scale

- `rot_y`: `y` on the rotated scale (a numeric vector)

- `stdrot_X_center`: numeric vector of values used to center `rot_X`

- `stdrot_X_scale`: numeric vector of values used to scale `rot_X`

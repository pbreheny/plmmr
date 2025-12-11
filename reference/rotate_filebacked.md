# A function to rotate filebacked data

A function to rotate filebacked data

## Usage

``` r
rotate_filebacked(prep, tocenter = TRUE, ...)
```

## Value

a list with 4 items:

- stdrot_X: `X` on the rotated and re-standardized scale

- rot_y: `y` on the rotated scale (a numeric vector)

- stdrot_X_center: numeric vector of values used to center `rot_X`

- stdrot_X_scale: numeric vector of values used to scale `rot_X`

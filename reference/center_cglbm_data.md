# Center reported matrix columns under censoring

Estimate one center per column from exact and one-sided censored
reports, using either a residual variance pooled across columns or one
variance per column. The function subtracts centers from exact values,
report bounds, latent simulation values, and fixed reporting-design
limits; it never divides by a scale.

## Usage

``` r
center_cglbm_data(
  data,
  variance = c("pooled", "column"),
  max_iter = 5000L,
  tolerance = 1e-11,
  score_tolerance = 1e-06,
  parameter_tolerance = 1e-07
)
```

## Arguments

- data:

  A `cglbm_data` object or numeric matrix interpreted as exact data.

- variance:

  `"pooled"` for one residual variance or `"column"` for one variance
  per column.

- max_iter:

  Maximum EM iterations.

- tolerance:

  Relative log-likelihood convergence tolerance.

- score_tolerance:

  Maximum absolute score used in convergence checks.

- parameter_tolerance:

  Maximum relative parameter change used in convergence checks.

## Value

A centered `cglbm_data` object. Estimation details are stored in
`metadata$centering`.

## Examples

``` r
x <- matrix(c(1, 2, 3, 4, 10, 11, 12, 13), nrow = 4)
centered <- center_cglbm_data(exact_cglbm_data(x), variance = "pooled")
centered$metadata$centering$centers
#> [1]  2.5 11.5
```

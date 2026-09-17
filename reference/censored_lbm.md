# Fit a censored Gaussian latent block model

Jointly estimate row memberships, column memberships, class proportions,
and Gaussian block parameters from exact, censored, and missing reports.
The five profiled numerical kernels use the compiled Rcpp backend;
starts, validity rules, convergence, and model control remain in R.

## Usage

``` r
censored_lbm(
  data,
  g,
  m,
  algorithm = "SS",
  S = c("Skl", "Sk", "Sl", "S0"),
  P = c("free", "equal"),
  n_starts = 20L,
  starts = NULL,
  seed = 1L,
  max_iter = 400L,
  tol = 1e-07,
  sigma2_min = 1e-08,
  effective_weight_min = 1e-08,
  monotone_tolerance = 1e-07,
  max_start_attempts = NULL,
  require_map_nonempty = TRUE,
  verbose = FALSE
)

# S3 method for class 'cglbm_fit'
print(x, ...)

# S3 method for class 'cglbm_fit'
summary(object, ...)

# S3 method for class 'cglbm_fit'
predict(object, type = c("labels", "memberships", "parameters"), ...)
```

## Arguments

- data:

  A
  [`new_cglbm_data()`](https://stagill.github.io/CGLBM/reference/new_cglbm_data.md)
  object or a numeric matrix interpreted as exact data.

- g, m:

  Positive numbers of row and column classes.

- algorithm:

  Membership configuration: `"HH"`, `"HS"`, `"SH"`, or `"SS"`; the first
  letter refers to rows and the second to columns.

- S:

  Variance structure: block-specific (`"Skl"`), shared by row class
  (`"Sk"`), shared by column class (`"Sl"`), or globally shared
  (`"S0"`).

- P:

  Class-proportion structure: estimated freely (`"free"`) or fixed
  equally (`"equal"`).

- n_starts:

  Number of valid starts requested.

- starts:

  Optional ordered list of starts, each containing integer vectors `z`
  and `w`.

- seed:

  Integer random seed.

- max_iter:

  Maximum GEM iterations per start.

- tol:

  Relative objective convergence tolerance.

- sigma2_min:

  Positive variance floor.

- effective_weight_min:

  Minimum admissible row, column, and observed block mass.

- monotone_tolerance:

  Allowed numerical decrease in a retained GEM substep.

- max_start_attempts:

  Maximum attempted starts.

- require_map_nonempty:

  Require every declared class in the final MAP partitions.

- verbose:

  Emit per-attempt progress messages.

- x:

  A `cglbm_fit` object.

- ...:

  Additional arguments, currently unused.

- object:

  A `cglbm_fit` object.

- type:

  One of `"labels"`, `"memberships"`, or `"parameters"`.

## Value

A `cglbm_fit` object containing memberships, MAP labels, fitted
parameters, objective diagnostics, and Fisher-information summaries.

`summary.cglbm_fit()` returns a data frame with one row per block;
`predict.cglbm_fit()` returns labels, memberships, or parameters.

## Examples

``` r
y <- matrix(c(-2.1, -1.9, 1.9, 2.1,
              -2.0, -1.8, 1.8, 2.0), nrow = 4)
dat <- exact_cglbm_data(y)
fit <- censored_lbm(
  dat, g = 2, m = 1, algorithm = "HH",
  n_starts = 1,
  starts = list(list(z = c(1L, 1L, 2L, 2L), w = rep(1L, 2))),
  max_iter = 30, max_start_attempts = 1, seed = 1
)
predict(fit, type = "labels")
#> $z
#> [1] 1 1 2 2
#> 
#> $w
#> [1] 1 1
#> 
summary(fit)
#>   K L fitted_mean standard_error fitted_variance effective_block_weight
#> 1 1 1       -1.95      0.0559017          0.0125                      4
#> 2 2 1        1.95      0.0559017          0.0125                      4
#>   exact_count left_censored_count right_censored_count nominal_cell_count
#> 1           4                   0                    0                  4
#> 2           4                   0                    0                  4
#>   effective_sample_size
#> 1                     4
#> 2                     4
```

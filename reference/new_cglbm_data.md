# Construct reported CGLBM data

Create the matrix-valued data representation used by
[`censored_lbm()`](https://stagill.github.io/CGLBM/reference/censored_lbm.md).
Each cell is recorded as `"exact"`, `"left"`, `"right"`, `"interval"`,
or `"missing"`, with report bounds on the same analysis scale as the
values. Missing cells are represented by the whole real line.

## Usage

``` r
new_cglbm_data(
  value,
  status = NULL,
  lower = NULL,
  upper = NULL,
  latent_truth = NULL,
  metadata = list()
)

exact_cglbm_data(y, metadata = list())

apply_detection_limits(
  y,
  lower_limit = -Inf,
  upper_limit = Inf,
  metadata = list()
)

# S3 method for class 'cglbm_data'
print(x, ...)
```

## Arguments

- value:

  Numeric matrix of exact values. Non-exact and missing cells may be
  `NA`.

- status:

  Optional character matrix of report types.

- lower, upper:

  Optional numeric matrices of report bounds.

- latent_truth:

  Optional latent matrix retained only for simulation diagnostics; it is
  never used in fitting.

- metadata:

  Optional list of application metadata. A fixed reporting rule may be
  supplied as `metadata$reporting_design$lower_limit` and
  `metadata$reporting_design$upper_limit`.

- y:

  Numeric matrix of latent or observed values; `NA` cells are recorded
  as missing.

- lower_limit, upper_limit:

  Scalar, column-specific vector, or matrix of fixed detection limits.

- x:

  A `cglbm_data` object.

- ...:

  Additional arguments, currently unused.

## Value

An object of class `cglbm_data`.

## Examples

``` r
x <- matrix(c(-2, -0.5, 0.5, 2), nrow = 2)
dat <- apply_detection_limits(x, lower_limit = -1, upper_limit = 1)
dat
#> CGLBM reported-data object
#>   dimensions: 2 x 2
#>   exact=2, left=1, right=1, interval=0
```

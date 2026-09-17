# Complete censored reports before co-clustering

Completion comparators used in the accompanying manuscript. Boundary and
half/double substitution are deterministic. Truncated-normal and
multiple imputation require explicit seeds. Completed matrices are
returned as exact `cglbm_data` objects.

## Usage

``` r
complete_boundary(data)

complete_halfdouble(data, analysis_scale = c("log", "log2"))

complete_truncnorm(data, seed)

complete_mi(data, seed, m = 10L, burnin = 50L, thin = 20L)

pool_mi_fits(fits)
```

## Arguments

- data:

  A `cglbm_data` object.

- analysis_scale:

  Scale used for half/double substitution.

- seed:

  Integer random seed.

- m:

  Number of multiple-imputation matrices.

- burnin:

  Number of initial sampler iterations discarded.

- thin:

  Iterations between retained imputations.

- fits:

  List of fitted CGLBM objects from multiple imputations.

## Value

`complete_boundary()`, `complete_halfdouble()`, and
`complete_truncnorm()` return an exact `cglbm_data` object.
`complete_mi()` returns a list with `completed` numeric matrices,
parameter `draws`, and sampler `diagnostics`. `pool_mi_fits()` returns a
pooled `cglbm_fit` after fitting the completed matrices separately.

## Examples

``` r
x <- matrix(c(-2, -0.5, 0.5, 2, -1.5, 0, 1, 2.5), nrow = 4)
dat <- apply_detection_limits(x, lower_limit = -1, upper_limit = 2)
complete_boundary(dat)
#> CGLBM reported-data object
#>   dimensions: 4 x 2
#>   exact=8, left=0, right=0, interval=0
complete_halfdouble(dat, analysis_scale = "log2")
#> CGLBM reported-data object
#>   dimensions: 4 x 2
#>   exact=8, left=0, right=0, interval=0
```

# Getting started with CGLBM

`CGLBM` co-clusters matrix rows and columns while using exact values and
censoring intervals in one reported-data likelihood.

``` r

library(CGLBM)

z0 <- rep(1:2, each = 4)
w0 <- rep(1:2, each = 4)
mu0 <- matrix(c(-2, 0.5, 0.5, 2.5), 2, 2)
y <- matrix(0, 8, 8)
for (i in seq_len(8)) for (j in seq_len(8)) {
  y[i, j] <- mu0[z0[i], w0[j]] + 0.02 * (i - j)
}

reported <- apply_detection_limits(
  y,
  lower_limit = -2,
  upper_limit = 2.5
)
reported
#> CGLBM reported-data object
#>   dimensions: 8 x 8
#>   exact=44, left=10, right=10, interval=0
```

Fit one model with a fixed start. The `S` argument controls the block
variance structure and `P` controls whether row and column class
proportions are free or equal.

``` r

fit <- censored_lbm(
  reported,
  g = 2,
  m = 2,
  algorithm = "HH",
  S = "Skl",
  P = "free",
  n_starts = 1,
  starts = list(list(z = z0, w = w0)),
  max_start_attempts = 1,
  max_iter = 50,
  seed = 1
)
fit
#> Censored Gaussian LBM fit (HH: hard rows / hard columns)
#>   model: g=2, m=2, S=Skl, P=free
#>   objective: 63.498409
#>   iterations: 26; converged: TRUE
#>   valid starts: 1/1 attempts
#>   MAP labels contain all classes: TRUE
predict(fit, type = "labels")
#> $z
#> [1] 1 1 1 1 2 2 2 2
#> 
#> $w
#> [1] 1 1 1 1 2 2 2 2
summary(fit)
#>   K L fitted_mean standard_error fitted_variance effective_block_weight
#> 1 1 1   -2.010528    0.014498566     0.001682111                     16
#> 2 1 2    0.420000    0.007905694     0.001000000                     16
#> 3 2 1    0.580000    0.007905694     0.001000000                     16
#> 4 2 2    2.510528    0.014498566     0.001682111                     16
#>   exact_count left_censored_count right_censored_count nominal_cell_count
#> 1           6                  10                    0                 16
#> 2          16                   0                    0                 16
#> 3          16                   0                    0                 16
#> 4           6                   0                   10                 16
#>   effective_sample_size
#> 1               11.9297
#> 2               16.0000
#> 3               16.0000
#> 4               11.9297
```

The summary reports fitted block means and variances, plug-in standard
errors, effective observed block weights, report counts, and
information-equivalent sample sizes.

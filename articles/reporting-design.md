# Reporting design and block-mean information

Fisher-information summaries use the fixed reporting rule, not only the
realized report status.
[`apply_detection_limits()`](https://stagill.github.io/CGLBM/reference/new_cglbm_data.md)
retains the complete lower and upper limits, including at cells that
happened to be reported exactly.

``` r

library(CGLBM)

y <- matrix(c(-2, -0.5, 0.5, 2, -1.5, 0, 1, 2.5), nrow = 4)
reported <- apply_detection_limits(y, lower_limit = -1, upper_limit = 2)
reported$metadata$reporting_design
#> $lower_limit
#> [1] -1
#> 
#> $upper_limit
#> [1] 2

fit <- censored_lbm(
  reported,
  g = 1,
  m = 1,
  algorithm = "SS",
  n_starts = 1,
  max_start_attempts = 1,
  max_iter = 50,
  seed = 3
)
summary(fit)
#>   K L fitted_mean standard_error fitted_variance effective_block_weight
#> 1 1 1   0.3582451      0.7925361        4.447167                      8
#>   exact_count left_censored_count right_censored_count nominal_cell_count
#> 1           4                   2                    2                  8
#>   effective_sample_size
#> 1                7.0924
```

For externally assembled reports, provide the complete rule through
`metadata$reporting_design`. Missing or inconsistent design metadata
does not invalidate the fit, but information-derived fields are returned
as `NA` with a reason. A finite interval-valued report is supported by
the likelihood, while its Fisher-information contribution is not
implemented in the first release.

# Selecting row and column class numbers

[`select_cglbm_model()`](https://stagill.github.io/CGLBM/reference/select_cglbm_model.md)
compares candidate values of `K` and `L` under one fixed membership
configuration, variance structure, and class-proportion structure.

``` r

library(CGLBM)

z0 <- rep(1:2, each = 4)
w0 <- rep(1:2, each = 4)
mu0 <- matrix(c(-2, 0.5, 0.5, 2.5), 2, 2)
y <- matrix(0, 8, 8)
for (i in seq_len(8)) for (j in seq_len(8)) {
  y[i, j] <- mu0[z0[i], w0[j]] + 0.02 * (i - j)
}
reported <- apply_detection_limits(y, -2, 2.5)

selection <- select_cglbm_model(
  reported,
  K = 1:2,
  L = 1:2,
  algorithm = "SS",
  S = "Skl",
  P = "free",
  n_starts = 1,
  max_start_attempts = 10,
  max_iter = 50,
  seed = 11
)
selection
#> CGLBM model selection by ICL
#>   selected model: K=2, L=2
#>   ICL: 44.783437
#>   eligible candidates: 4/4
selection$table[, c("K", "L", "eligible", "icl")]
#>   K L eligible        icl
#> 1 1 1     TRUE -123.01466
#> 2 1 2     TRUE -115.24910
#> 3 2 1     TRUE -112.07063
#> 4 2 2     TRUE   44.78344
```

Candidates without a valid fit are excluded. `S` and `P` remain fixed
during one search; the function does not automatically search over those
structures.

# CGLBM

<!-- badges: start -->
[![R-CMD-check](https://github.com/StaGill/CGLBM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/StaGill/CGLBM/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/StaGill/CGLBM/actions/workflows/pkgdown.yaml/badge.svg)](https://stagill.github.io/CGLBM/)
<!-- badges: end -->

Documentation: <https://stagill.github.io/CGLBM/>

`CGLBM` fits censored Gaussian latent block models for co-clustering matrix
rows and columns when cells may be exact, left-censored, right-censored,
interval-censored, or missing. The reported-data likelihood is optimized
jointly over row memberships, column memberships, class proportions, and
Gaussian block parameters. The package also provides ICL-based class-number
selection, block-mean standard errors from Fisher information, completion
comparators, and censoring-aware column centering.

## Installation

```r
# install.packages("remotes")
remotes::install_github("StaGill/CGLBM")
```

## Quick start

```r
library(CGLBM)

y <- matrix(
  c(-2.1, -1.9,  1.9,  2.1,
    -2.0, -1.8,  1.8,  2.0),
  nrow = 4
)

dat <- exact_cglbm_data(y)

fit <- censored_lbm(
  dat,
  g = 2,
  m = 1,
  algorithm = "HH",
  S = "Skl",
  P = "free",
  n_starts = 1,
  starts = list(list(
    z = c(1L, 1L, 2L, 2L),
    w = rep(1L, 2)
  )),
  max_start_attempts = 1,
  seed = 1
)

predict(fit, type = "labels")
summary(fit)
```

For detection-limit data, construct reports and preserve the complete fixed
reporting rule with `apply_detection_limits()`:

```r
reported <- apply_detection_limits(
  y,
  lower_limit = -1.5,
  upper_limit = 1.5
)
```

## Model structures

The fitting and model-selection functions support all fixed combinations of:

- membership configuration: `HH`, `HS`, `SH`, or `SS`
- variance structure: `Skl`, `Sk`, `Sl`, or `S0`
- class proportions: `free` or `equal`

`select_cglbm_model()` searches only the supplied candidate values of `K` and
`L`; `S` and `P` remain fixed during one search.

## Documentation

See `vignette("getting-started", package = "CGLBM")` for the core workflow,
`vignette("model-selection", package = "CGLBM")` for ICL-based class-number
selection, and `vignette("reporting-design", package = "CGLBM")` for fixed
reporting rules and block-mean information summaries.

## Reproducibility project

The installable package and the manuscript reproduction project are maintained
as separate directories. The reproduction project depends on this package and
retains the formal simulation manifests, deterministic seed schedules,
manuscript-ready CSVs, rendering scripts, and private-data boundary checks.
The CARE-SSc participant-level data are not distributed with either project.

## Manuscript setting

The accompanying manuscript *Model-Based Co-Clustering for Censored Data* uses
block-specific variances and free row and column class proportions
(`S = "Skl"`, `P = "free"`). The package additionally supports shared-variance
structures and equal class proportions.

## License

MIT (c) 2026 Kejun Fang.

# Changelog

## CGLBM 0.1.0

- Initial package release.
- Fits censored Gaussian latent block models under HH, HS, SH, and SS
  membership configurations.
- Supports exact, left-censored, right-censored, finite-interval, and
  missing reports.
- Supports `Skl`, `Sk`, `Sl`, and `S0` variance structures with free or
  equal class proportions.
- Provides ICL-based class-number selection, block-mean
  Fisher-information summaries, completion comparators, and
  censoring-aware column centering.
- Uses a package-native Rcpp backend for the five profiled numerical
  kernels.

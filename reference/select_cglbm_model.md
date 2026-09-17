# Select row and column class numbers by ICL

Fit a grid of candidate class-number pairs under one fixed membership,
variance, and class-proportion structure, then maximize the implemented
BIC approximation to the integrated classification likelihood.

## Usage

``` r
cglbm_icl_bic(fit, data)

select_cglbm_model(
  data,
  K = 1:5,
  L = 1:5,
  algorithm = "SS",
  n_starts = 20L,
  seed = 1L,
  max_iter = 400L,
  tol = 1e-07,
  sigma2_min = 1e-08,
  effective_weight_min = 1e-08,
  monotone_tolerance = 1e-07,
  max_start_attempts = 200L,
  verbose = FALSE,
  S = "Skl",
  P = "free"
)

# S3 method for class 'cglbm_model_selection'
print(x, ...)
```

## Arguments

- fit:

  A fitted `cglbm_fit` object.

- data:

  A `cglbm_data` object or numeric matrix.

- K, L:

  Candidate row and column class numbers.

- algorithm, S, P:

  Fixed model structures passed to
  [`censored_lbm()`](https://stagill.github.io/CGLBM/reference/censored_lbm.md).

- n_starts:

  Number of valid starts requested for each candidate.

- seed:

  Base integer seed used to construct deterministic candidate seeds.

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

  Maximum attempted starts for each candidate.

- verbose:

  Emit per-candidate progress messages.

- x:

  A `cglbm_model_selection` object.

- ...:

  Additional arguments, currently unused.

## Value

`cglbm_icl_bic()` returns one numeric ICL score.

A `cglbm_model_selection` object containing the candidate table,
retained fits, and selected fit.

## Examples

``` r
y <- matrix(seq_len(24), nrow = 6)
sel <- select_cglbm_model(
  exact_cglbm_data(y), K = 1, L = 1,
  n_starts = 1, max_start_attempts = 1, max_iter = 10, seed = 1
)
sel
#> CGLBM model selection by ICL
#>   selected model: K=1, L=1
#>   ICL: -83.666139
#>   eligible candidates: 1/1
```

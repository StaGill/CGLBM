#' CGLBM: Censored Gaussian Latent Block Models
#'
#' CGLBM fits Gaussian latent block models to exact, censored, and missing
#' matrix-valued reports. It provides hard and soft co-clustering algorithms,
#' fixed-structure ICL model selection, Fisher-information summaries,
#' completion comparators, and censoring-aware column centering.
#'
#' @useDynLib CGLBM, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats dnorm integrate pnorm qnorm rgamma rnorm runif sd var
#' @importFrom utils tail
#' @keywords internal
"_PACKAGE"

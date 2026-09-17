make_completion_fixture <- function() {
  status <- matrix(c(
    "exact", "left", "right",
    "interval", "exact", "left",
    "right", "exact", "exact"
  ), nrow = 3, byrow = TRUE)
  value <- matrix(c(
    0.2, NA, NA,
    NA, 0.8, NA,
    NA, 1.4, 2.0
  ), nrow = 3, byrow = TRUE)
  lower <- matrix(NA_real_, 3, 3)
  upper <- matrix(NA_real_, 3, 3)
  upper[1, 2] <- 0
  lower[1, 3] <- 1
  lower[2, 1] <- -0.5; upper[2, 1] <- 0.5
  upper[2, 3] <- 1.2
  lower[3, 1] <- 1.5
  new_cglbm_data(value, status, lower, upper)
}

test_that("extreme-tail truncated draws remain random and inside support", {
  draw <- getFromNamespace("comparison_rtruncnorm_inverse", "CGLBM")
  set.seed(20260909)
  samples <- list(
    draw(64L, 10, 11, 0, 1),
    draw(64L, -40, -39, 0, 1),
    draw(64L, 10, Inf, 0, 1),
    draw(64L, -Inf, -10, 0, 1)
  )
  expect_true(all(samples[[1]] > 10 & samples[[1]] < 11))
  expect_true(all(samples[[2]] > -40 & samples[[2]] < -39))
  expect_true(all(samples[[3]] > 10))
  expect_true(all(samples[[4]] < -10))
  expect_true(all(vapply(samples, function(x) {
    all(is.finite(x)) && length(unique(signif(x, 12))) > 1L
  }, logical(1))))
})

test_that("multiple imputation retains the manuscript schedule and prior", {
  dat <- make_completion_fixture()
  set.seed(2468)
  before <- .Random.seed
  a <- complete_mi(dat, seed = 20260908L)
  after <- .Random.seed
  b <- complete_mi(dat, seed = 20260908L)

  expect_identical(after, before)
  expect_length(a$completed, 10)
  expect_equal(a$diagnostics$save_iterations, 50L + seq_len(10L) * 20L)
  expect_equal(a$diagnostics$prior$kappa0, 0.01, tolerance = 0)
  expect_equal(a$diagnostics$prior$alpha0, 2.01, tolerance = 0)
  expect_equal(a$diagnostics$prior$beta0,
               1.01 * a$diagnostics$prior$s20, tolerance = 1e-14)
  expect_true(isTRUE(a$diagnostics$exact_preserved))
  expect_true(isTRUE(a$diagnostics$interval_valid))
  expect_equal(a$completed, b$completed, tolerance = 0)
})

test_that("assignment-based label alignment preserves physical blocks", {
  one_hot <- getFromNamespace("one_hot", "CGLBM")
  align <- getFromNamespace("comparison_align_fit", "CGLBM")
  make_fit <- function(z, w, mu, variance) {
    K <- nrow(mu); L <- ncol(mu)
    structure(list(
      z = as.integer(z), w = as.integer(w),
      tau = one_hot(z, K), omega = one_hot(w, L),
      params = list(
        pi = tabulate(z, nbins = K) / length(z),
        rho = tabulate(w, nbins = L) / length(w),
        mu = mu, var = variance
      ),
      model = list(g = K, m = L, S = "Skl", P = "free"),
      objective = 0
    ), class = "cglbm_fit")
  }

  mu <- matrix(c(1, 2, 3, 4), 2, 2)
  variance <- matrix(c(1.1, 1.2, 1.3, 1.4), 2, 2)
  a <- make_fit(c(1, 1, 2, 2), c(1, 2, 1), mu, variance)
  b <- make_fit(
    c(2, 2, 1, 1), c(2, 1, 2),
    mu[c(2, 1), c(2, 1), drop = FALSE],
    variance[c(2, 1), c(2, 1), drop = FALSE]
  )
  aligned <- align(b, a, 2L, 2L)$fit
  expect_identical(aligned$z, a$z)
  expect_identical(aligned$w, a$w)
  expect_equal(aligned$params$mu, a$params$mu, tolerance = 0)
  expect_equal(aligned$params$var, a$params$var, tolerance = 0)

  pooled <- pool_mi_fits(list(a, b))
  expect_equal(pooled$mi_diagnostics$pairwise_cari[1, 2], 1, tolerance = 0)
  expect_equal(pooled$params$mu, mu, tolerance = 0)
  expect_equal(pooled$params$var, variance, tolerance = 0)
})

test_that("alignment and CARI retain the polynomial-time implementations", {
  align <- getFromNamespace("comparison_align_fit", "CGLBM")
  align_text <- paste(deparse(body(align)), collapse = "\n")
  expect_true(grepl("clue::solve_LSAP", align_text, fixed = TRUE))
  expect_false(grepl("comparison_permutations", align_text, fixed = TRUE))

  cari_text <- paste(deparse(body(comparison_cari)), collapse = "\n")
  expect_true(grepl("kronecker", cari_text, fixed = TRUE))
  expect_false(grepl("paste(", cari_text, fixed = TRUE))

  expanded_reference <- function(z_a, w_a, z_b, w_b) {
    labels_a <- paste(
      rep(z_a, times = length(w_a)),
      rep(w_a, each = length(z_a)),
      sep = ":"
    )
    labels_b <- paste(
      rep(z_b, times = length(w_b)),
      rep(w_b, each = length(z_b)),
      sep = ":"
    )
    choose2 <- function(x) x * (x - 1) / 2
    n <- length(labels_a)
    if (n < 2L) return(1)
    tab <- table(labels_a, labels_b)
    A <- sum(choose2(rowSums(tab)))
    B <- sum(choose2(colSums(tab)))
    expected <- A * B / choose2(n)
    denominator <- (A + B) / 2 - expected
    if (abs(denominator) < 1e-12) return(1)
    as.numeric((sum(choose2(tab)) - expected) / denominator)
  }

  cases <- list(
    list(
      z_a = c(1L, 1L, 2L, 2L),
      w_a = c(1L, 2L, 1L),
      z_b = c(1L, 2L, 1L, 2L),
      w_b = c(2L, 2L, 1L)
    ),
    list(
      z_a = c(1L, 1L, 1L, 2L, 2L, 3L),
      w_a = c(1L, 1L, 2L, 3L),
      z_b = c(2L, 2L, 2L, 1L, 1L, 3L),
      w_b = c(3L, 3L, 1L, 2L)
    )
  )
  for (case in cases) {
    expect_equal(
      do.call(comparison_cari, case),
      do.call(expanded_reference, case),
      tolerance = 1e-12
    )
  }
})

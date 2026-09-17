test_that("seed-locked likelihood and truncated-moment references are preserved", {
  logprob <- getFromNamespace("log_normal_interval_prob", "CGLBM")
  moments <- getFromNamespace("truncated_normal_moments", "CGLBM")

  expect_equal(logprob(300, 306), -45006.622732118663,
               tolerance = 1e-9)

  central <- moments(-1, 1, 0, 1)
  central_variance <- 1 - 2 * dnorm(1) / (pnorm(1) - pnorm(-1))
  expect_equal(central$first, 0, tolerance = 1e-15)
  expect_equal(central$variance, central_variance, tolerance = 1e-14)

  translated <- moments(1e8, Inf, 1e8, 1)
  expect_equal(translated$first, 1e8 + sqrt(2 / pi), tolerance = 3e-8)
  expect_equal(translated$variance, 1 - 2 / pi, tolerance = 1e-12)

  deep <- moments(300, 306, 0, 1)
  reflected <- moments(-306, -300, 0, 1)
  expect_equal(deep$first, 300.00333325926337415, tolerance = 2e-10)
  expect_equal(deep$variance, 1.1110370438949582e-05, tolerance = 1e-12)
  expect_equal(reflected$first, -deep$first, tolerance = 2e-10)
  expect_equal(reflected$variance, deep$variance, tolerance = 1e-12)

  very_deep <- moments(1000, 1010, 0, 1)
  expect_equal(very_deep$first, 1000.00099999800001, tolerance = 2e-10)
  expect_equal(very_deep$variance, 9.999940000499995e-07,
               tolerance = 1e-12)

  narrow <- moments(300, 300.000001, 0, 1)
  expect_gt(narrow$first, 300)
  expect_lt(narrow$first, 300.000001)
  expect_gte(narrow$variance, 0)
  expect_lte(narrow$variance, (1e-6)^2 / 4)
  expect_equal(narrow$variance, (1e-6)^2 / 12, tolerance = 3e-5)
})

test_that("report-likelihood dispatch retains fixed numerical references", {
  cell_loglik <- getFromNamespace("cell_loglik_matrix", "CGLBM")
  dat <- new_cglbm_data(
    value = matrix(c(0, NA_real_, NA_real_, NA_real_), nrow = 1),
    status = matrix(c("exact", "left", "right", "interval"), nrow = 1),
    lower = matrix(c(NA_real_, -Inf, 0, -1), nrow = 1),
    upper = matrix(c(NA_real_, 0, Inf, 1), nrow = 1)
  )
  expected <- c(
    dnorm(0, log = TRUE),
    log(0.5),
    log(0.5),
    -0.38171514630212616
  )
  expect_equal(as.numeric(cell_loglik(dat, mean = 0, variance = 1)),
               expected, tolerance = 1e-15)
})

test_that("completion methods preserve the caller RNG state", {
  x <- matrix(c(-2, -0.5, 0.5, 2, -1.5, 0, 1, 2.5), nrow = 4)
  dat <- apply_detection_limits(x, lower_limit = -1, upper_limit = 2)

  set.seed(2718)
  state_before <- .Random.seed
  invisible(complete_truncnorm(dat, seed = 31415))
  expect_identical(.Random.seed, state_before)

  invisible(complete_mi(dat, seed = 1618, m = 2, burnin = 2, thin = 1))
  expect_identical(.Random.seed, state_before)
})

test_that("membership ties and numerical-zero masses retain locked behavior", {
  update <- getFromNamespace("update_membership_from_scores", "CGLBM")
  validate <- getFromNamespace("validate_membership_masses", "CGLBM")
  spec <- getFromNamespace("algorithm_spec", "CGLBM")

  tie_scores <- matrix(c(1, 1, 2, 2), nrow = 2, byrow = TRUE)
  hard <- update(tie_scores, "H")
  soft <- update(tie_scores, "S")
  expect_identical(max.col(hard, ties.method = "first"), c(1L, 1L))
  expect_equal(rowSums(soft), c(1, 1), tolerance = 1e-12)
  expect_equal(soft[, 1], c(0.5, 0.5), tolerance = 1e-12)

  tau <- matrix(c(1 - 2e-9, 2e-9, 1 - 2e-9, 2e-9),
                nrow = 2, byrow = TRUE)
  check <- validate(
    tau = tau,
    omega = matrix(1, nrow = 1, ncol = 1),
    spec = spec("SS"),
    effective_weight_min = 1e-8
  )
  expect_false(check$ok)
  expect_identical(check$reason, "effective row-class mass too small")
})

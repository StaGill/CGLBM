test_that("column centering subtracts centers without scaling", {
  x <- matrix(c(1, 2, 3, 4, 10, 11, 12, 13), nrow = 4)
  pooled <- center_cglbm_data(exact_cglbm_data(x), variance = "pooled")
  column <- center_cglbm_data(exact_cglbm_data(x), variance = "column")
  expect_equal(colMeans(pooled$value), c(0, 0), tolerance = 1e-10)
  expect_equal(colMeans(column$value), c(0, 0), tolerance = 1e-10)
  expect_equal(pooled$value,
               sweep(x, 2, pooled$metadata$centering$centers, "-"),
               tolerance = 1e-12)
  expect_length(pooled$metadata$centering$variance, 1)
  expect_length(column$metadata$centering$variance, 2)
})

test_that("centering shifts reports and fixed detection limits together", {
  x <- matrix(c(-2, -0.5, 0.5, 2, -1.5, 0, 1, 2.5), nrow = 4)
  dat <- apply_detection_limits(x, lower_limit = -1, upper_limit = 2)
  centers <- center_cglbm_data(dat, variance = "pooled")
  cvec <- centers$metadata$centering$centers
  design <- centers$metadata$reporting_design
  expected_lower <- sweep(matrix(-1, nrow(x), ncol(x)), 2, cvec, "-")
  expected_upper <- sweep(matrix(2, nrow(x), ncol(x)), 2, cvec, "-")
  expect_equal(design$lower_limit, expected_lower, tolerance = 1e-10)
  expect_equal(design$upper_limit, expected_upper, tolerance = 1e-10)
})

test_that("censored pooled centering preserves the frozen CARE-style fixture", {
  status <- matrix(
    c("exact", "exact", "left", "exact", "right", "exact",
      "left", "exact", "exact", "exact", "exact", "exact"),
    nrow = 4, ncol = 3, byrow = TRUE
  )
  value <- matrix(
    c(1, 5, NA, 2, NA, 9, NA, 6, 10, 3, 7, 11),
    nrow = 4, ncol = 3, byrow = TRUE
  )
  lower <- matrix(
    c(1, 5, -Inf, 2, 7, 9, -Inf, 6, 10, 3, 7, 11),
    nrow = 4, ncol = 3, byrow = TRUE
  )
  upper <- matrix(
    c(1, 5, 8, 2, Inf, 9, 1.5, 6, 10, 3, 7, 11),
    nrow = 4, ncol = 3, byrow = TRUE
  )
  dat <- new_cglbm_data(value, status, lower, upper)
  before <- dat
  centered <- center_cglbm_data(dat, variance = "pooled")
  info <- centered$metadata$centering

  expect_identical(dat, before)
  expect_true(info$converged)
  expect_identical(info$iterations, 16L)
  expect_lte(info$max_score, 1e-6)
  expect_true(all(diff(info$trace) >= -1e-8))
  expect_equal(
    info$centers,
    c(1.6623236342689385, 6.432287265298132, 9.361267136696938),
    tolerance = 2e-10
  )
  expect_equal(info$sigma, 1.1370416321540893, tolerance = 2e-10)
  expect_equal(info$auxiliary_loglik, -16.715437454957655, tolerance = 2e-10)
  expect_equal(centered$upper[1, 3], -1.361267136696938, tolerance = 2e-10)
  expect_equal(centered$lower[2, 2], 0.5677127347018676, tolerance = 2e-10)
  expect_equal(centered$upper[3, 1], -0.16232363426893848, tolerance = 2e-10)
})

test_that("deterministic completion methods use the documented substitutions", {
  x <- matrix(c(-2, -0.5, 0.5, 2, -1.5, 0, 1, 2.5), nrow = 4)
  dat <- apply_detection_limits(x, lower_limit = -1, upper_limit = 2)
  boundary <- complete_boundary(dat)
  halfdouble <- complete_halfdouble(dat, analysis_scale = "log2")
  expect_true(all(boundary$status == "exact"))
  expect_true(all(halfdouble$status == "exact"))
  expect_true(all(boundary$value[dat$status == "left"] == -1))
  expect_true(all(boundary$value[dat$status == "right"] == 2))
  expect_true(all(halfdouble$value[dat$status == "left"] == -2))
  expect_true(all(halfdouble$value[dat$status == "right"] == 3))
})

test_that("stochastic completion is seed reproducible and respects bounds", {
  x <- matrix(c(-2, -0.5, 0.5, 2, -1.5, 0, 1, 2.5), nrow = 4)
  dat <- apply_detection_limits(x, lower_limit = -1, upper_limit = 2)
  a <- complete_truncnorm(dat, seed = 77)
  b <- complete_truncnorm(dat, seed = 77)
  expect_equal(a$value, b$value, tolerance = 0)
  expect_true(all(a$value[dat$status == "left"] < -1))
  expect_true(all(a$value[dat$status == "right"] > 2))

  mi <- complete_mi(dat, seed = 88, m = 2, burnin = 2, thin = 1)
  expect_length(mi$completed, 2)
  expect_true(all(vapply(mi$completed, is.matrix, logical(1))))
  expect_true(isTRUE(mi$diagnostics$exact_preserved))
  expect_true(isTRUE(mi$diagnostics$interval_valid))
})

test_that("CARI is invariant to separate label permutations", {
  z <- c(1, 1, 2, 2)
  w <- c(1, 2, 1)
  expect_equal(comparison_cari(z, w, 3 - z, 3 - w), 1)
})

test_that("MI fits pool without reusing single-fit information", {
  x <- matrix(c(-2, -0.5, 0.5, 2, -1.5, 0, 1, 2.5), nrow = 4)
  dat <- apply_detection_limits(x, lower_limit = -1, upper_limit = 2)
  imputed <- complete_mi(dat, seed = 99, m = 2, burnin = 2, thin = 1)
  fits <- lapply(seq_along(imputed$completed), function(index) {
    censored_lbm(
      imputed$completed[[index]], g = 1, m = 1,
      n_starts = 1, max_start_attempts = 1, max_iter = 20,
      seed = 100 + index
    )
  })
  pooled <- pool_mi_fits(fits)
  expect_s3_class(pooled, "cglbm_fit")
  expect_null(pooled[["information"]])
  expect_warning(out <- summary(pooled), "MI-pooled")
  expect_true(is.na(out$standard_error))
})

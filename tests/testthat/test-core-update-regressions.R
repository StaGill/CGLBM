test_that("Gaussian updates retain translation invariance and the variance floor", {
  update <- getFromNamespace("weighted_gaussian_em_step", "CGLBM")
  model <- list(g = 1L, m = 1L, S = "Skl", P = "free")
  tau <- matrix(1, 2, 1)
  omega <- matrix(1, 1, 1)

  run_exact <- function(values, old_mean = mean(values), sigma2_min = 1e-12) {
    update(
      data = exact_cglbm_data(matrix(values, ncol = 1)),
      tau = tau,
      omega = omega,
      old_params = list(
        pi = 1, rho = 1,
        mu = matrix(old_mean, 1, 1),
        var = matrix(1, 1, 1)
      ),
      model = model,
      sigma2_min = sigma2_min,
      effective_weight_min = 0,
      iteration = 1L
    )
  }

  base <- run_exact(c(0, 1))
  shifted <- run_exact(c(1e8, 1e8 + 1), old_mean = 1e8 + 0.5)
  expect_true(base$ok)
  expect_true(shifted$ok)
  expect_equal(as.numeric(base$mu), 0.5, tolerance = 1e-12)
  expect_lt(abs(as.numeric(shifted$mu) - (1e8 + 0.5)), 3e-8)
  expect_equal(as.numeric(base$var), 0.25, tolerance = 1e-12)
  expect_equal(as.numeric(shifted$var), 0.25, tolerance = 1e-12)
  expect_equal(as.numeric(base$block_ss), 0.5, tolerance = 1e-12)
  expect_equal(as.numeric(shifted$block_ss), 0.5, tolerance = 1e-12)

  floored <- run_exact(c(4, 4), old_mean = 4, sigma2_min = 1e-8)
  expect_true(floored$ok)
  expect_equal(floored$var, matrix(1e-8, 1, 1), tolerance = 0)
  expect_true(isTRUE(floored$variance_floor_history$activated[[1L]]))
})

test_that("mixed exact and censored updates are translation invariant", {
  update <- getFromNamespace("weighted_gaussian_em_step", "CGLBM")
  model <- list(g = 1L, m = 1L, S = "Skl", P = "free")
  tau <- matrix(1, 2, 1)
  omega <- matrix(1, 1, 1)

  run_mixed <- function(shift) {
    dat <- new_cglbm_data(
      value = matrix(c(shift, NA_real_), ncol = 1),
      status = matrix(c("exact", "right"), ncol = 1),
      lower = matrix(c(NA_real_, shift + 1), ncol = 1),
      upper = matrix(c(NA_real_, Inf), ncol = 1)
    )
    update(
      data = dat, tau = tau, omega = omega,
      old_params = list(
        pi = 1, rho = 1,
        mu = matrix(shift, 1, 1),
        var = matrix(1, 1, 1)
      ),
      model = model, sigma2_min = 1e-12,
      effective_weight_min = 0, iteration = 1L
    )
  }

  base <- run_mixed(0)
  shifted <- run_mixed(1e8)
  expect_true(base$ok)
  expect_true(shifted$ok)
  expect_lt(abs((as.numeric(shifted$mu) - 1e8) - as.numeric(base$mu)), 3e-8)
  expect_lt(abs(as.numeric(shifted$var) - as.numeric(base$var)), 1e-8)
})

test_that("complete one-sided censoring is rejected before a Gaussian update", {
  update <- getFromNamespace("weighted_gaussian_em_step", "CGLBM")
  model <- list(g = 1L, m = 1L, S = "Skl", P = "free")
  tau <- matrix(1, 2, 1)
  omega <- matrix(1, 1, 1)
  params <- list(pi = 1, rho = 1, mu = matrix(0, 1, 1), var = matrix(1, 1, 1))

  left <- new_cglbm_data(
    value = matrix(NA_real_, 2, 1),
    status = matrix("left", 2, 1),
    lower = matrix(-Inf, 2, 1),
    upper = matrix(c(0, 1), 2, 1)
  )
  right <- new_cglbm_data(
    value = matrix(NA_real_, 2, 1),
    status = matrix("right", 2, 1),
    lower = matrix(c(0, 1), 2, 1),
    upper = matrix(Inf, 2, 1)
  )

  left_result <- update(left, tau, omega, params, model, 1e-8, 0, 1L)
  right_result <- update(right, tau, omega, params, model, 1e-8, 0, 1L)
  expect_false(left_result$ok)
  expect_match(left_result$reason, "complete left censoring", fixed = TRUE)
  expect_false(right_result$ok)
  expect_match(right_result$reason, "complete right censoring", fixed = TRUE)
})

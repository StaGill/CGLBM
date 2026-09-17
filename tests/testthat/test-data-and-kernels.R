test_that("reported-data constructors preserve all report types", {
  value <- matrix(c(0, NA, NA, NA, NA), nrow = 1)
  status <- matrix(c("exact", "left", "right", "interval", "missing"),
                   nrow = 1)
  lower <- matrix(c(NA, -Inf, 0, -1, -Inf), nrow = 1)
  upper <- matrix(c(NA, 0, Inf, 1, Inf), nrow = 1)
  dat <- new_cglbm_data(value, status, lower, upper)
  expect_s3_class(dat, "cglbm_data")
  expect_identical(as.character(dat$status), as.character(status))
  expect_true(is.na(dat$value[1, 5]))
  expect_equal(dat$lower[1, 5], -Inf)
  expect_equal(dat$upper[1, 5], Inf)
})

test_that("stable interval probabilities retain seed-locked references", {
  logprob <- getFromNamespace("log_normal_interval_prob", "CGLBM")
  expect_equal(logprob(-Inf, 0), log(0.5), tolerance = 1e-15)
  expect_equal(logprob(0, Inf), log(0.5), tolerance = 1e-15)
  expect_equal(logprob(-1, 1), -0.38171514630212616, tolerance = 1e-15)
  expect_equal(logprob(10, 11), -53.23131022558312, tolerance = 1e-12)
})

test_that("compiled kernels agree with the retained R references", {
  fixture <- make_cglbm_fixture(censor = TRUE, missing = TRUE)
  dat <- fixture$data
  model <- list(g = 2L, m = 2L, S = "Skl", P = "free")
  params <- list(
    pi = c(0.5, 0.5), rho = c(0.5, 0.5),
    mu = matrix(c(-1.8, 0.4, 0.4, 2.2), 2, 2),
    var = matrix(c(0.8, 1.1, 1.2, 0.9), 2, 2)
  )
  tau <- matrix(c(0.9, 0.1, 0.8, 0.2, 0.7, 0.3, 0.6, 0.4,
                  0.4, 0.6, 0.3, 0.7, 0.2, 0.8, 0.1, 0.9), 8, 2, byrow = TRUE)
  omega <- tau

  cpp_log <- getFromNamespace("all_cell_loglik", "CGLBM")(dat, params, model)
  r_log <- getFromNamespace(".cglbm_r_all_cell_loglik", "CGLBM")(dat, params, model)
  expect_equal(cpp_log, r_log, tolerance = 1e-12)

  cpp_row <- getFromNamespace("row_membership_scores", "CGLBM")(
    omega, params, cpp_log, model)
  r_row <- getFromNamespace(".cglbm_r_row_membership_scores", "CGLBM")(
    omega, params, r_log, model)
  expect_equal(cpp_row, r_row, tolerance = 1e-12)

  cpp_col <- getFromNamespace("col_membership_scores", "CGLBM")(
    tau, params, cpp_log, model)
  r_col <- getFromNamespace(".cglbm_r_col_membership_scores", "CGLBM")(
    tau, params, r_log, model)
  expect_equal(cpp_col, r_col, tolerance = 1e-12)

  cpp_obj <- getFromNamespace("membership_objective", "CGLBM")(
    tau, omega, params, cpp_log, model)
  r_obj <- getFromNamespace(".cglbm_r_membership_objective", "CGLBM")(
    tau, omega, params, r_log, model)
  expect_equal(cpp_obj, r_obj, tolerance = 1e-11)

  cpp_step <- getFromNamespace("weighted_gaussian_em_step", "CGLBM")(
    dat, tau, omega, params, model, sigma2_min = 1e-8,
    effective_weight_min = 1e-8, iteration = 1L)
  r_step <- getFromNamespace(".cglbm_r_weighted_gaussian_em_step", "CGLBM")(
    dat, tau, omega, params, model, sigma2_min = 1e-8,
    effective_weight_min = 1e-8, iteration = 1L)
  expect_identical(cpp_step$ok, r_step$ok)
  expect_equal(cpp_step$mu, r_step$mu, tolerance = 1e-10)
  expect_equal(cpp_step$var, r_step$var, tolerance = 1e-10)
  expect_equal(cpp_step$block_weight, r_step$block_weight, tolerance = 1e-12)
  expect_equal(cpp_step$block_ss, r_step$block_ss, tolerance = 1e-9)
})

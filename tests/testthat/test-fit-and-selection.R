test_that("all membership configurations fit a fixed start", {
  fixture <- make_cglbm_fixture(censor = FALSE)
  for (algorithm in c("HH", "HS", "SH", "SS")) {
    fit <- fit_fixture(fixture$data, fixture$start, algorithm = algorithm)
    expect_s3_class(fit, "cglbm_fit")
    expect_true(is.finite(fit$objective))
    expect_true(all(fit$substep_increments >= -1e-8))
    expect_setequal(unique(fit$z), 1:2)
    expect_setequal(unique(fit$w), 1:2)
  }
})

test_that("all eight fixed S/P combinations are supported", {
  fixture <- make_cglbm_fixture(censor = FALSE)
  for (S in c("Skl", "Sk", "Sl", "S0")) {
    for (P in c("free", "equal")) {
      fit <- fit_fixture(fixture$data, fixture$start, S = S, P = P)
      expect_identical(fit$model$S, S)
      expect_identical(fit$model$P, P)
      expect_true(all(is.finite(fit$params$mu)))
      expect_true(all(is.finite(fit$params$var)))
      out <- summary(fit)
      expect_equal(nrow(out), 4)
      expect_named(out, c(
        "K", "L", "fitted_mean", "standard_error", "fitted_variance",
        "effective_block_weight", "exact_count", "left_censored_count",
        "right_censored_count", "nominal_cell_count",
        "effective_sample_size"
      ))
    }
  }
})

test_that("fixed seeds reproduce the same fitted result", {
  fixture <- make_cglbm_fixture(censor = TRUE)
  a <- fit_fixture(fixture$data, fixture$start, algorithm = "SS", seed = 813L)
  b <- fit_fixture(fixture$data, fixture$start, algorithm = "SS", seed = 813L)
  expect_identical(a$z, b$z)
  expect_identical(a$w, b$w)
  expect_equal(a$objective, b$objective, tolerance = 0)
  expect_equal(a$params$mu, b$params$mu, tolerance = 0)
  expect_equal(a$params$var, b$params$var, tolerance = 0)
})

test_that("generalized ICL parameter counts match each S/P combination", {
  penalty <- getFromNamespace("cglbm_icl_penalty", "CGLBM")
  K <- 3L; L <- 2L; n <- 120L; m <- 60L
  q <- c(Skl = K * L, Sk = K, Sl = L, S0 = 1L)
  for (S in names(q)) for (P in c("free", "equal")) {
    expected <- (K * L + q[[S]]) / 2 * log(as.numeric(n) * m)
    if (P == "free") {
      expected <- expected + (K - 1) / 2 * log(n) + (L - 1) / 2 * log(m)
    }
    expect_equal(penalty(K, L, n, m, S = S, P = P), expected,
                 tolerance = 1e-12)
  }
  expect_true(is.finite(penalty(3, 3, 60000L, 40000L, "Skl", "free")))
})

test_that("model selection returns a consistent selected fit", {
  fixture <- make_cglbm_fixture(censor = TRUE)
  selection <- select_cglbm_model(
    fixture$data, K = 1:2, L = 1:2, algorithm = "SS",
    S = "Sk", P = "equal", n_starts = 1,
    max_start_attempts = 10, max_iter = 50, seed = 919L
  )
  expect_s3_class(selection, "cglbm_model_selection")
  expect_true(selection$success)
  expect_true(selection$table$eligible[selection$table$K == selection$selected_K &
                                       selection$table$L == selection$selected_L])
  expect_identical(selection$selected_fit$model$S, "Sk")
  expect_identical(selection$selected_fit$model$P, "equal")
  expect_equal(cglbm_icl_bic(selection$selected_fit, fixture$data),
               selection$selected_icl, tolerance = 1e-10)
})

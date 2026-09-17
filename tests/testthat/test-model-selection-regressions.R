test_that("classification scoring and ICL agree with direct Gaussian calculations", {
  components <- getFromNamespace("cglbm_classification_loglik_components", "CGLBM")
  class_loglik <- getFromNamespace("cglbm_classification_loglik", "CGLBM")
  penalty <- getFromNamespace("cglbm_icl_penalty", "CGLBM")

  y <- matrix(c(-1, 0, 1, 2), nrow = 2, byrow = TRUE)
  dat <- exact_cglbm_data(y)
  model <- list(g = 1L, m = 1L, S = "Skl", P = "free")
  params <- list(
    pi = 1, rho = 1,
    mu = matrix(0.25, 1, 1),
    var = matrix(1.5, 1, 1)
  )
  z <- rep(1L, nrow(y)); w <- rep(1L, ncol(y))
  direct <- sum(dnorm(y, mean = 0.25, sd = sqrt(1.5), log = TRUE))
  score <- components(dat, z, w, params, model)
  expect_equal(score, direct, tolerance = 1e-12)

  fit <- structure(list(
    z = z, w = w, params = params, model = model,
    data_dimensions = c(n = nrow(y), d = ncol(y)), map_nonempty = TRUE
  ), class = "cglbm_fit")
  expect_equal(class_loglik(fit, dat), score, tolerance = 1e-12)
  expect_equal(cglbm_icl_bic(fit, dat),
               score - penalty(1, 1, nrow(y), ncol(y), "Skl", "free"),
               tolerance = 1e-12)
})

test_that("large candidate dimensions do not overflow parameter counts", {
  penalty <- getFromNamespace("cglbm_icl_penalty", "CGLBM")
  grid <- getFromNamespace("cglbm_candidate_grid", "CGLBM")

  value <- penalty(3L, 3L, 60000L, 40000L, "Skl", "free")
  expected <-
    (3 - 1) / 2 * log(60000) +
    (3 - 1) / 2 * log(40000) +
    3 * 3 * log(as.double(60000) * as.double(40000))
  expect_true(is.finite(value))
  expect_equal(value, expected, tolerance = 1e-12)

  candidate <- grid(50000L, 50000L, 60000L, 60000L,
                    S = "Skl", P = "free")
  expected_count <-
    (as.double(50000) - 1) +
    (as.double(50000) - 1) +
    2 * as.double(50000) * as.double(50000)
  expect_equal(nrow(candidate), 1)
  expect_true(is.finite(candidate$parameter_count[[1L]]))
  expect_equal(candidate$parameter_count[[1L]], expected_count, tolerance = 0)
})

test_that("classification scores reject incomplete labels unless explicitly permitted", {
  components <- getFromNamespace("cglbm_classification_loglik_components", "CGLBM")
  dat <- exact_cglbm_data(matrix(0, nrow = 3, ncol = 2))
  model <- list(g = 2L, m = 2L, S = "Skl", P = "free")
  params <- list(
    pi = c(0.5, 0.5), rho = c(0.5, 0.5),
    mu = matrix(0, 2, 2), var = matrix(1, 2, 2)
  )
  z <- c(1L, 1L, 1L); w <- c(1L, 2L)
  expect_error(components(dat, z, w, params, model), "declared classes")
  expect_true(is.finite(components(
    dat, z, w, params, model, require_all_classes = FALSE
  )))
})

test_that("missing cells do not alter a clear fixed-start partition", {
  full <- make_cglbm_fixture(censor = TRUE, missing = FALSE)
  miss <- make_cglbm_fixture(censor = TRUE, missing = TRUE)
  fit_full <- fit_fixture(full$data, full$start)
  fit_miss <- fit_fixture(miss$data, miss$start)
  expect_equal(comparison_cari(fit_full$z, fit_full$w,
                               fit_miss$z, fit_miss$w), 1)
  expect_lt(max(abs(fit_full$params$mu - fit_miss$params$mu)), 0.1)
})

test_that("exact-data mean SE reduces to sqrt(variance / block weight)", {
  y <- matrix(c(-1, 0, 1, 2, -0.5, 0.5, 1.5, 2.5), nrow = 4)
  fit <- censored_lbm(
    exact_cglbm_data(y), g = 1, m = 1,
    n_starts = 1, max_start_attempts = 1, max_iter = 20, seed = 1
  )
  out <- summary(fit)
  expect_equal(out$standard_error,
               sqrt(out$fitted_variance / out$effective_block_weight),
               tolerance = 1e-10)
  expect_equal(out$effective_sample_size, out$effective_block_weight,
               tolerance = 1e-12)
})

test_that("fixed detection rules yield available information summaries", {
  fixture <- make_cglbm_fixture(censor = TRUE)
  fit <- fit_fixture(fixture$data, fixture$start, algorithm = "SS")
  out <- summary(fit)
  expect_true(all(is.finite(out$standard_error)))
  expect_true(all(out$effective_sample_size <= out$effective_block_weight + 1e-10))
  expect_true(all(out$effective_sample_size > 0))
})

test_that("fixed-rule information does not depend on the realized status", {
  cell_information <- getFromNamespace(".cglbm_cell_information", "CGLBM")

  left_rule <- apply_detection_limits(
    matrix(c(-1, 1), nrow = 1),
    lower_limit = 0,
    upper_limit = Inf
  )
  expect_identical(as.vector(left_rule$status), c("left", "exact"))
  left <- cell_information(left_rule, mean = 0, variance = 1)
  A0 <- 0.5 + 1 / pi
  expect_equal(as.numeric(left$A), rep(A0, 2), tolerance = 1e-14)
  expect_equal(as.numeric(left$I11), rep(A0, 2), tolerance = 1e-14)

  right_rule <- apply_detection_limits(
    matrix(c(-1, 1), nrow = 1),
    lower_limit = -Inf,
    upper_limit = 0
  )
  right <- cell_information(right_rule, mean = 0, variance = 1)
  expect_equal(as.numeric(right$I11), as.numeric(left$I11), tolerance = 1e-14)
  expect_equal(as.numeric(right$I12), -as.numeric(left$I12), tolerance = 1e-14)
  expect_equal(as.numeric(right$I22), as.numeric(left$I22), tolerance = 1e-14)
})

test_that("shared-variance mean SE uses the joint information matrix", {
  attach_information <- getFromNamespace(".cglbm_attach_information", "CGLBM")
  one_hot <- getFromNamespace("one_hot", "CGLBM")

  y <- matrix(
    c(-2, 0, 1, 2, -0.5, 0.3, 1.2, -1, 0.4, 2, -2, 0.2,
      1, -0.2, 2, -1, 0.7, 0.9, -0.1, 1, 2, 0, -1, 0.6),
    nrow = 4, ncol = 6
  )
  y[1, 2] <- NA_real_
  dat <- apply_detection_limits(y, lower_limit = -0.7, upper_limit = 1.3)
  z <- c(1L, 1L, 2L, 2L)
  w <- rep(1:3, each = 2L)
  fit <- structure(
    list(
      model = list(g = 2L, m = 3L, S = "S0", P = "free"),
      z = z,
      w = w,
      tau = one_hot(z, 2L),
      omega = one_hot(w, 3L),
      params = list(
        mu = matrix(c(-0.8, 0.6, 0.2, -0.9, 1.1, 0.7), 2L, 3L),
        var = matrix(1.05, 2L, 3L),
        pi = c(0.5, 0.5),
        rho = rep(1 / 3, 3L)
      ),
      data_dimensions = c(n = 4L, d = 6L),
      map_nonempty = TRUE
    ),
    class = "cglbm_fit"
  )

  fit <- attach_information(fit, dat)
  expect_true(all(is.finite(fit$information$mean_standard_error)))
  expect_length(fit$information$variance_groups, 1L)
  expect_identical(fit$information$variance_groups[[1L]]$status, "ok")

  separate <- matrix(NA_real_, 2L, 3L)
  for (k in 1:2) for (l in 1:3) {
    separate[k, l] <- sqrt(
      solve(fit$information$block_information[k, l, , ])[1, 1]
    )
  }
  expect_gt(
    max(abs(separate - fit$information$mean_standard_error)),
    1e-4
  )
})

test_that("finite interval reports produce explicit unavailable information", {
  dat <- new_cglbm_data(
    value = matrix(c(0, NA, 1, 2), nrow = 2),
    status = matrix(c("exact", "interval", "exact", "exact"), nrow = 2),
    lower = matrix(c(NA, -0.5, NA, NA), nrow = 2),
    upper = matrix(c(NA, 0.5, NA, NA), nrow = 2),
    metadata = list(reporting_design = list(lower_limit = -1, upper_limit = 1))
  )
  fit <- censored_lbm(
    dat, g = 1, m = 1, n_starts = 1,
    max_start_attempts = 1, max_iter = 20, seed = 1
  )
  expect_warning(out <- summary(fit), "finite interval-censored")
  expect_true(is.na(out$standard_error))
  expect_true(is.na(out$effective_sample_size))
})

test_that("one percent missing at random preserves clear exact and censored fits", {
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(
    20260916L,
    kind = "Mersenne-Twister",
    normal.kind = "Inversion",
    sample.kind = "Rejection"
  )
  z <- rep(1:3, each = 40L)
  w <- rep(1:3, each = 20L)
  block_means <- outer(c(-3, 0, 3), c(-1, 0, 1), "+")
  y <- block_means[z, w] + matrix(rnorm(120L * 60L, sd = 0.8), 120L, 60L)

  set.seed(
    20260917L,
    kind = "Mersenne-Twister",
    normal.kind = "Inversion",
    sample.kind = "Rejection"
  )
  mask <- matrix(FALSE, 120L, 60L)
  mask[sample.int(length(mask), size = 72L)] <- TRUE
  expect_equal(sum(mask), 0.01 * length(mask), tolerance = 0)

  adjusted_rand <- getFromNamespace("comparison_adjusted_rand", "CGLBM")
  align_fit <- getFromNamespace("comparison_align_fit", "CGLBM")
  starts <- list(list(z = z, w = w))
  controls <- list(
    g = 3L,
    m = 3L,
    algorithm = "SS",
    n_starts = 1L,
    starts = starts,
    seed = 20260918L,
    max_iter = 400L,
    tol = 1e-8,
    max_start_attempts = 1L
  )

  for (kind in c("exact", "censored")) {
    intact <- if (kind == "exact") {
      exact_cglbm_data(y)
    } else {
      apply_detection_limits(y, lower_limit = -3, upper_limit = 3)
    }
    intact$latent_truth <- NULL

    status <- intact$status
    status[mask] <- "missing"
    masked <- new_cglbm_data(
      value = intact$value,
      status = status,
      lower = intact$lower,
      upper = intact$upper,
      metadata = list(reporting_design = intact$metadata$reporting_design)
    )

    fit_full <- do.call(censored_lbm, c(list(data = intact), controls))
    fit_missing <- do.call(censored_lbm, c(list(data = masked), controls))

    expect_true(isTRUE(fit_full$converged), info = paste(kind, "intact fit"))
    expect_true(isTRUE(fit_missing$converged), info = paste(kind, "missing fit"))

    row_ari <- adjusted_rand(fit_full$z, fit_missing$z)
    column_ari <- adjusted_rand(fit_full$w, fit_missing$w)
    cari <- comparison_cari(
      fit_full$z, fit_full$w, fit_missing$z, fit_missing$w
    )
    aligned <- align_fit(fit_missing, fit_full, K = 3L, L = 3L)$fit
    mean_difference <- max(abs(aligned$params$mu - fit_full$params$mu))

    expect_true(row_ari >= 0.95, info = paste(kind, "row ARI"))
    expect_true(column_ari >= 0.95, info = paste(kind, "column ARI"))
    expect_true(cari >= 0.95, info = paste(kind, "CARI"))
    expect_true(mean_difference <= 0.10, info = paste(kind, "block means"))
    expect_true(
      all(fit_missing$substep_increments >= -1e-7),
      info = paste(kind, "substep monotonicity")
    )
  }
})

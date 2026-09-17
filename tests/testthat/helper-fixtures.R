make_cglbm_fixture <- function(censor = TRUE, missing = FALSE) {
  z <- rep(1:2, each = 4)
  w <- rep(1:2, each = 4)
  mu <- matrix(c(-2, 0.5, 0.5, 2.5), 2, 2)
  y <- matrix(0, 8, 8)
  for (i in seq_len(8)) for (j in seq_len(8)) {
    y[i, j] <- mu[z[i], w[j]] + 0.02 * (i - j)
  }
  if (missing) y[1, 8] <- NA_real_
  data <- if (censor) {
    apply_detection_limits(y, lower_limit = -2, upper_limit = 2.5)
  } else {
    exact_cglbm_data(y)
  }
  list(data = data, y = y, z = z, w = w,
       start = list(z = z, w = w))
}

fit_fixture <- function(data, start, algorithm = "HH", S = "Skl", P = "free",
                        seed = 1L) {
  censored_lbm(
    data, g = 2, m = 2, algorithm = algorithm, S = S, P = P,
    n_starts = 1, starts = list(start), seed = seed,
    max_iter = 50, max_start_attempts = 1,
    monotone_tolerance = 1e-8
  )
}

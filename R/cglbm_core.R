# Core functions for the Censored Gaussian Latent Block Model (CGLBM).
# The model combines exact observations and interval-censored reports in one
# Gaussian latent block likelihood. HH, HS, SH, and SS specify hard or soft
# membership domains for the row and column partitions.
#
# Boundary and midpoint values are used only to initialize block parameters.
# Likelihood evaluations and Gaussian updates use exact observations and the
# reported censoring intervals.

.cglbm_or <- function(x, y) if (is.null(x)) y else x

.cglbm_stop <- function(...) stop(sprintf(...), call. = FALSE)
.cglbm_assert <- function(condition, ...) {
  if (!isTRUE(condition)) .cglbm_stop(...)
  invisible(TRUE)
}

as_numeric_matrix <- function(x, name = deparse(substitute(x))) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (length(dim(x)) != 2L) .cglbm_stop("%s must be a two-dimensional matrix.", name)
  x
}

expand_to_matrix <- function(x, n, d, name) {
  if (length(x) == 1L) return(matrix(as.numeric(x), nrow = n, ncol = d))
  if (is.atomic(x) && is.null(dim(x)) && length(x) == d) {
    return(matrix(rep(as.numeric(x), each = n), nrow = n, ncol = d))
  }
  x <- as_numeric_matrix(x, name)
  if (!identical(dim(x), c(n, d))) {
    .cglbm_stop("%s must be scalar, length d, or an n x d matrix.", name)
  }
  x
}

normalize_labels <- function(x) {
  x <- as.vector(x)
  match(x, sort(unique(x)))
}

one_hot <- function(labels, n_classes) {
  labels <- as.integer(labels)
  if (anyNA(labels) || any(labels < 1L | labels > n_classes)) {
    .cglbm_stop("Labels must be integers between 1 and n_classes.")
  }
  out <- matrix(0, nrow = length(labels), ncol = n_classes)
  out[cbind(seq_along(labels), labels)] <- 1
  out
}

sample_nonempty_labels <- function(n, g) {
  .cglbm_assert(n >= g, "Cannot create %d nonempty classes from %d items.", g, n)
  labels <- c(seq_len(g), if (n > g) sample.int(g, n - g, replace = TRUE) else integer())
  sample(labels, size = n, replace = FALSE)
}

generate_start_partitions <- function(n, d, g, m, n_starts, seed = 1L) {
  n_starts <- as.integer(n_starts)
  .cglbm_assert(n_starts >= 1L, "n_starts must be at least 1.")
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else NULL
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(
  as.integer(seed),
  kind = "Mersenne-Twister",
  normal.kind = "Inversion")
  lapply(seq_len(n_starts), function(s) list(
    z = sample_nonempty_labels(n, g),
    w = sample_nonempty_labels(d, m),
    start_id = s
  ))
}

logsumexp <- function(x) {
  if (any(is.infinite(x) & x > 0, na.rm = TRUE)) return(Inf)
  x <- x[is.finite(x)]
  if (!length(x)) return(-Inf)
  mx <- max(x)
  mx + log(sum(exp(x - mx)))
}

row_log_softmax <- function(scores) {
  scores <- as.matrix(scores)
  out <- matrix(NA_real_, nrow(scores), ncol(scores))
  for (i in seq_len(nrow(scores))) {
    mx <- max(scores[i, ])
    if (!is.finite(mx)) .cglbm_stop("All variational scores are non-finite for item %d.", i)
    out[i, ] <- scores[i, ] - (mx + log(sum(exp(scores[i, ] - mx))))
  }
  out
}
row_softmax <- function(scores) exp(row_log_softmax(scores))

xlogx_sum <- function(x) {
  idx <- x > 0
  sum(x[idx] * log(x[idx]))
}

logdiffexp <- function(log_x, log_y) {
  log_x <- as.numeric(log_x)
  log_y <- rep_len(as.numeric(log_y), length(log_x))
  if (any(log_y > log_x, na.rm = TRUE)) .cglbm_stop("logdiffexp requires log_x >= log_y.")
  out <- rep(NA_real_, length(log_x))
  same <- is.finite(log_x) & is.finite(log_y) & log_x == log_y
  out[same] <- -Inf
  idx <- !same & is.finite(log_x)
  out[idx] <- log_x[idx] + log(-expm1(log_y[idx] - log_x[idx]))
  out[is.infinite(log_x) & log_x < 0] <- -Inf
  out
}


.log_standard_normal_interval_prob_numerical <- function(a, b, rel_tol = 1e-12) {
  if (!is.finite(a) || !is.finite(b) || a >= b) {
    .cglbm_stop("Numerical interval integration requires finite endpoints with a < b.")
  }
  center <- a / 2 + b / 2
  half_width <- (b - a) / 2
  anchor <- min(max(0, a), b)
  log_anchor <- dnorm(anchor, log = TRUE)
  scaled <- tryCatch(
    stats::integrate(
      function(t) {
        z <- center + half_width * t
        exp(dnorm(z, log = TRUE) - log_anchor)
      },
      lower = -1, upper = 1, rel.tol = rel_tol,
      subdivisions = 200L, stop.on.error = TRUE
    )$value,
    error = function(e) NA_real_
  )
  if (!is.finite(scaled) || scaled <= 0) {
    .cglbm_stop("Failed to compute a positive Gaussian interval probability.")
  }
  log_anchor + log(half_width) + log(scaled)
}

log_normal_interval_prob <- function(lower, upper, mean = 0, sd = 1) {
  target_length <- max(length(lower), length(upper), length(mean), length(sd))
  lower <- rep_len(as.numeric(lower), target_length)
  upper <- rep_len(as.numeric(upper), target_length)
  mean <- rep_len(as.numeric(mean), target_length)
  sd <- rep_len(as.numeric(sd), target_length)
  if (any(!is.finite(sd) | sd <= 0)) .cglbm_stop("All standard deviations must be positive and finite.")
  if (any(lower >= upper, na.rm = TRUE)) .cglbm_stop("Every censored interval must satisfy lower < upper.")

  a <- (lower - mean) / sd
  b <- (upper - mean) / sd
  out <- rep(NA_real_, target_length)

  left_inf <- is.infinite(a) & a < 0
  right_inf <- is.infinite(b) & b > 0
  out[left_inf & right_inf] <- 0

  idx <- left_inf & !right_inf
  if (any(idx)) out[idx] <- pnorm(b[idx], log.p = TRUE)
  idx <- right_inf & !left_inf
  if (any(idx)) out[idx] <- pnorm(a[idx], lower.tail = FALSE, log.p = TRUE)

  finite <- is.finite(a) & is.finite(b)
  neg <- finite & b <= 0
  pos <- finite & a >= 0
  cross <- finite & !neg & !pos
  if (any(neg)) {
    out[neg] <- logdiffexp(pnorm(b[neg], log.p = TRUE), pnorm(a[neg], log.p = TRUE))
  }
  if (any(pos)) {
    out[pos] <- logdiffexp(
      pnorm(a[pos], lower.tail = FALSE, log.p = TRUE),
      pnorm(b[pos], lower.tail = FALSE, log.p = TRUE)
    )
  }
  if (any(cross)) out[cross] <- log(pnorm(b[cross]) - pnorm(a[cross]))

  standardized_width <- b - a
  finite_scale <- pmax(1, abs(a), abs(b))
  narrow_finite <- finite & standardized_width <= 1e-5 * finite_scale
  if (any(narrow_finite)) {
    for (ii in which(narrow_finite)) {
      out[ii] <- .log_standard_normal_interval_prob_numerical(a[ii], b[ii])
    }
  }

  bad_finite <- finite & !is.finite(out)
  if (any(bad_finite)) {
    for (ii in which(bad_finite)) {
      out[ii] <- .log_standard_normal_interval_prob_numerical(a[ii], b[ii])
    }
  }
  if (anyNA(out) || any(!is.finite(out))) {
    .cglbm_stop("Failed to compute one or more positive Gaussian interval probabilities.")
  }
  out
}

.truncated_normal_positive_finite_tail_moments <- function(a, b,
                                                           rel_tol = 1e-10) {
  if (!is.finite(a) || !is.finite(b) || a <= 0 || a >= b) {
    .cglbm_stop(
      "Positive finite-tail moments require 0 < a < b with finite endpoints."
    )
  }

  upper_t <- a * (b - a)
  if (!is.finite(upper_t) || upper_t <= 0) {
    .cglbm_stop("Failed to construct a finite-tail integration interval.")
  }

  cutoff <- max(40, -log(rel_tol) + 10)
  if (upper_t <= cutoff) {
    # Scaling the finite interval to [0, 1] keeps very short tail intervals
    # numerically visible to the quadrature rule.
    scaled_kernel <- function(u) {
      t <- upper_t * u
      exp(-t - t^2 / (2 * a^2))
    }
    integrate_unit_interval <- function(f) {
      tryCatch(
        stats::integrate(
          f, lower = 0, upper = 1, rel.tol = rel_tol, abs.tol = rel_tol,
          subdivisions = 200L, stop.on.error = TRUE
        )$value,
        error = function(e) NA_real_
      )
    }

    mass <- integrate_unit_interval(scaled_kernel)
    if (!is.finite(mass) || mass <= 0) {
      .cglbm_stop("Failed to compute a finite-tail probability mass.")
    }
    mean_u <- integrate_unit_interval(function(u) u * scaled_kernel(u)) / mass
    if (!is.finite(mean_u)) {
      .cglbm_stop("Failed to compute a finite-tail first moment.")
    }
    variance_u <- integrate_unit_interval(
      function(u) (u - mean_u)^2 * scaled_kernel(u)
    ) / mass
    mean_t <- upper_t * mean_u
    variance_t <- upper_t^2 * variance_u
  } else {
    # The omitted transformed tail is negligible at the requested tolerance
    # for the probability mass and the first two centered moments.
    cutoff_kernel <- function(t) exp(-t - t^2 / (2 * a^2))
    integrate_cutoff_interval <- function(f) {
      tryCatch(
        stats::integrate(
          f, lower = 0, upper = cutoff, rel.tol = rel_tol,
          abs.tol = rel_tol, subdivisions = 200L, stop.on.error = TRUE
        )$value,
        error = function(e) NA_real_
      )
    }

    mass <- integrate_cutoff_interval(cutoff_kernel)
    if (!is.finite(mass) || mass <= 0) {
      .cglbm_stop("Failed to compute a finite-tail probability mass.")
    }
    mean_t <- integrate_cutoff_interval(function(t) t * cutoff_kernel(t)) / mass
    if (!is.finite(mean_t)) {
      .cglbm_stop("Failed to compute a finite-tail first moment.")
    }
    variance_t <- integrate_cutoff_interval(
      function(t) (t - mean_t)^2 * cutoff_kernel(t)
    ) / mass
  }

  standardized_mean <- a + mean_t / a
  standardized_variance <- variance_t / a^2
  if (!is.finite(standardized_variance) || standardized_variance < 0) {
    .cglbm_stop("Failed to compute a finite-tail variance.")
  }
  list(mean = standardized_mean, variance = standardized_variance)
}

.truncated_normal_moments_numerical <- function(a, b, rel_tol = 1e-10) {
  left_tail <- is.infinite(a) && a < 0 && is.finite(b) && b < 0
  right_tail <- is.finite(a) && a > 0 && is.infinite(b) && b > 0
  if (left_tail || right_tail) {
    endpoint <- if (left_tail) b else a
    scale <- abs(endpoint)
    one_sided_kernel <- function(t) exp(-t - t^2 / (2 * scale^2))
    integrate_one_sided_tail <- function(f) {
      tryCatch(
        stats::integrate(
          f, lower = 0, upper = Inf, rel.tol = rel_tol,
          subdivisions = 200L, stop.on.error = TRUE
        )$value,
        error = function(e) NA_real_
      )
    }

    mass <- integrate_one_sided_tail(one_sided_kernel)
    if (!is.finite(mass) || mass <= 0) {
      .cglbm_stop("Failed to compute a truncated-normal probability mass.")
    }
    mean_t <- integrate_one_sided_tail(function(t) t * one_sided_kernel(t)) / mass
    if (!is.finite(mean_t)) {
      .cglbm_stop("Failed to compute a truncated-normal first moment.")
    }
    variance_t <- integrate_one_sided_tail(
      function(t) (t - mean_t)^2 * one_sided_kernel(t)
    ) / mass
    direction <- if (left_tail) -1 else 1
    standardized_mean <- endpoint + direction * mean_t / scale
    standardized_variance <- variance_t / scale^2
    if (!is.finite(standardized_variance) || standardized_variance < 0) {
      .cglbm_stop("Failed to compute a truncated-normal variance.")
    }
    return(list(mean = standardized_mean, variance = standardized_variance))
  }

  if (is.finite(a) && is.finite(b)) {
    if (a >= b) .cglbm_stop("Finite truncated-normal endpoints must satisfy a < b.")

    if (a >= 8) {
      return(.truncated_normal_positive_finite_tail_moments(a, b, rel_tol))
    }
    if (b <= -8) {
      reflected <- .truncated_normal_positive_finite_tail_moments(
        -b, -a, rel_tol
      )
      return(list(mean = -reflected$mean, variance = reflected$variance))
    }

    center <- a / 2 + b / 2
    half_width <- (b - a) / 2
    anchor <- min(max(0, a), b)
    log_anchor <- dnorm(anchor, log = TRUE)
    finite_interval_kernel <- function(t) {
      z <- center + half_width * t
      exp(dnorm(z, log = TRUE) - log_anchor)
    }
    integrate_finite_interval <- function(f, abs_tol = rel_tol) {
      tryCatch(
        stats::integrate(
          f, lower = -1, upper = 1, rel.tol = rel_tol, abs.tol = abs_tol,
          subdivisions = 200L, stop.on.error = TRUE
        )$value,
        error = function(e) NA_real_
      )
    }

    mass <- integrate_finite_interval(finite_interval_kernel)
    if (!is.finite(mass) || mass <= 0) {
      .cglbm_stop("Failed to compute a truncated-normal probability mass.")
    }
    first_t <- integrate_finite_interval(function(t) t * finite_interval_kernel(t)) / mass
    if (!is.finite(first_t)) {
      .cglbm_stop("Failed to compute a truncated-normal first moment.")
    }
    standardized_mean <- center + half_width * first_t
    variance_t <- integrate_finite_interval(
      function(t) (t - first_t)^2 * finite_interval_kernel(t)
    ) / mass
    standardized_variance <- half_width^2 * variance_t
    if (!is.finite(standardized_variance) || standardized_variance < 0) {
      .cglbm_stop("Failed to compute a truncated-normal variance.")
    }
    return(list(mean = standardized_mean, variance = standardized_variance))
  }

  anchor <- min(max(0, a), b)
  log_anchor <- dnorm(anchor, log = TRUE)
  general_kernel <- function(z) exp(dnorm(z, log = TRUE) - log_anchor)
  integrate_general_interval <- function(f) {
    tryCatch(
      stats::integrate(
        f, lower = a, upper = b, rel.tol = rel_tol,
        subdivisions = 200L, stop.on.error = TRUE
      )$value,
      error = function(e) NA_real_
    )
  }

  mass <- integrate_general_interval(general_kernel)
  if (!is.finite(mass) || mass <= 0) {
    .cglbm_stop("Failed to compute a truncated-normal probability mass.")
  }
  centered_first <- integrate_general_interval(
    function(z) (z - anchor) * general_kernel(z)
  )
  if (!is.finite(centered_first)) {
    .cglbm_stop("Failed to compute a truncated-normal first moment.")
  }
  standardized_mean <- anchor + centered_first / mass
  centered_second <- integrate_general_interval(
    function(z) (z - standardized_mean)^2 * general_kernel(z)
  )
  if (!is.finite(centered_second) || centered_second < 0) {
    .cglbm_stop("Failed to compute a truncated-normal variance.")
  }
  list(mean = standardized_mean, variance = centered_second / mass)
}

truncated_normal_moments <- function(lower, upper, mean = 0, variance = 1,
                                     numerical_tolerance = 1e-10) {
  numerical_tolerance <- as.numeric(numerical_tolerance)
  if (length(numerical_tolerance) != 1L || !is.finite(numerical_tolerance) ||
      numerical_tolerance <= 0) {
    .cglbm_stop("numerical_tolerance must be positive and finite.")
  }
  target_length <- max(length(lower), length(upper), length(mean), length(variance))
  lower <- rep_len(as.numeric(lower), target_length)
  upper <- rep_len(as.numeric(upper), target_length)
  mean <- rep_len(as.numeric(mean), target_length)
  variance <- rep_len(as.numeric(variance), target_length)
  if (any(!is.finite(variance) | variance <= 0)) {
    .cglbm_stop("All variances must be positive and finite.")
  }
  if (any(lower >= upper, na.rm = TRUE)) {
    .cglbm_stop("Every censored interval must satisfy lower < upper.")
  }

  sd <- sqrt(variance)
  a <- (lower - mean) / sd
  b <- (upper - mean) / sd
  log_delta <- log_normal_interval_prob(lower, upper, mean, sd)

  ratio_a <- numeric(target_length)
  ratio_b <- numeric(target_length)
  finite_a <- is.finite(a)
  finite_b <- is.finite(b)
  ratio_a[finite_a] <- exp(dnorm(a[finite_a], log = TRUE) - log_delta[finite_a])
  ratio_b[finite_b] <- exp(dnorm(b[finite_b], log = TRUE) - log_delta[finite_b])

  lambda <- ratio_a - ratio_b
  a_term <- numeric(target_length)
  b_term <- numeric(target_length)
  a_term[finite_a] <- a[finite_a] * ratio_a[finite_a]
  b_term[finite_b] <- b[finite_b] * ratio_b[finite_b]
  chi <- a_term - b_term

  standardized_mean <- lambda
  standardized_second <- 1 + chi
  standardized_variance <- standardized_second - standardized_mean^2
  first <- mean + sd * standardized_mean
  cond_var <- variance * standardized_variance
  second <- first^2 + cond_var

  finite_interval <- is.finite(lower) & is.finite(upper)
  standardized_width <- b - a
  finite_scale <- pmax(1, abs(a), abs(b))
  narrow_finite <- finite_interval & standardized_width <= 1e-3 * finite_scale
  deep_finite_tail <- finite_interval & (a >= 8 | b <= -8)
  extreme_one_sided <-
    (is.infinite(a) & a < 0 & is.finite(b) & b <= -30) |
    (is.finite(a) & a >= 30 & is.infinite(b) & b > 0)
  support_violation <- rep(FALSE, target_length)
  mean_outside <- rep(FALSE, target_length)
  if (any(finite_interval)) {
    width <- upper[finite_interval] - lower[finite_interval]
    variance_bound <- width^2 / 4
    variance_slack <- 10 * .Machine$double.eps * pmax(1, width^2)
    support_violation[finite_interval] <-
      cond_var[finite_interval] > variance_bound + variance_slack

    endpoint_scale <- pmax(
      1, abs(lower[finite_interval]), abs(upper[finite_interval])
    )
    mean_slack <- 10 * .Machine$double.eps * endpoint_scale
    mean_outside[finite_interval] <-
      first[finite_interval] < lower[finite_interval] - mean_slack |
      first[finite_interval] > upper[finite_interval] + mean_slack
  }

  # Direct standardized formulas avoid cancellation from a large location shift.
  # Numerical centered calculations are still needed for narrow intervals and
  # deep same-tail truncation, where the standardized variance itself is a
  # difference between large, nearly equal terms.
  bad <- !is.finite(first) | !is.finite(second) | !is.finite(cond_var) |
    cond_var < 0 | narrow_finite | deep_finite_tail | extreme_one_sided |
    support_violation | mean_outside
  if (any(bad)) {
    rel_tol <- max(as.numeric(numerical_tolerance), 1e-12)
    for (ii in which(bad)) {
      stable <- .truncated_normal_moments_numerical(a[ii], b[ii], rel_tol)
      first[ii] <- mean[ii] + sd[ii] * stable$mean
      cond_var[ii] <- variance[ii] * stable$variance
      second[ii] <- first[ii]^2 + cond_var[ii]
    }
  }

  if (any(!is.finite(first)) || any(!is.finite(second)) || any(cond_var < 0)) {
    .cglbm_stop("Failed to compute valid truncated-normal moments.")
  }
  list(first = first, second = second, variance = cond_var, log_probability = log_delta)
}

#' Construct reported CGLBM data
#'
#' Create the matrix-valued data representation used by [censored_lbm()].
#' Each cell is recorded as `"exact"`, `"left"`, `"right"`, `"interval"`,
#' or `"missing"`, with report bounds on the same analysis scale as the
#' values. Missing cells are represented by the whole real line.
#'
#' @param value Numeric matrix of exact values. Non-exact and missing cells
#'   may be `NA`.
#' @param status Optional character matrix of report types.
#' @param lower,upper Optional numeric matrices of report bounds.
#' @param latent_truth Optional latent matrix retained only for simulation
#'   diagnostics; it is never used in fitting.
#' @param metadata Optional list of application metadata. A fixed reporting
#'   rule may be supplied as `metadata$reporting_design$lower_limit` and
#'   `metadata$reporting_design$upper_limit`.
#'
#' @return An object of class `cglbm_data`.
#' @export
#' @examples
#' x <- matrix(c(-2, -0.5, 0.5, 2), nrow = 2)
#' dat <- apply_detection_limits(x, lower_limit = -1, upper_limit = 1)
#' dat
new_cglbm_data <- function(value, status = NULL, lower = NULL, upper = NULL,
                           latent_truth = NULL, metadata = list()) {
  value <- as_numeric_matrix(value, "value")
  n <- nrow(value); d <- ncol(value)
  if (is.null(status)) {
    status <- matrix("exact", n, d)
    status[is.na(value)] <- "missing"
  } else {
    status <- as.matrix(status)
    if (!identical(dim(status), c(n, d))) .cglbm_stop("status must have the same dimensions as value.")
    status <- tolower(as.character(status)); dim(status) <- c(n, d)
  }
  allowed <- c("exact", "left", "right", "interval", "missing")
  if (any(!status %in% allowed)) .cglbm_stop("status entries must be exact, left, right, interval, or missing.")
  if (is.null(lower)) lower <- matrix(NA_real_, n, d)
  if (is.null(upper)) upper <- matrix(NA_real_, n, d)
  lower <- expand_to_matrix(lower, n, d, "lower")
  upper <- expand_to_matrix(upper, n, d, "upper")

  exact <- status == "exact"; left <- status == "left"
  right <- status == "right"; interval <- status == "interval"
  missing <- status == "missing"
  if (any(!is.finite(value[exact]))) .cglbm_stop("Every exact cell must have a finite value.")
  lower[exact] <- value[exact]; upper[exact] <- value[exact]
  value[!exact] <- NA_real_
  lower[missing] <- -Inf; upper[missing] <- Inf
  lower[left] <- -Inf
  if (any(!is.finite(upper[left]))) .cglbm_stop("Every left-censored cell needs a finite upper endpoint.")
  upper[right] <- Inf
  if (any(!is.finite(lower[right]))) .cglbm_stop("Every right-censored cell needs a finite lower endpoint.")
  if (any(interval)) {
    if (any(!is.finite(lower[interval]) | !is.finite(upper[interval]))) {
      .cglbm_stop("Finite interval-censored cells need two finite endpoints.")
    }
    if (any(lower[interval] >= upper[interval])) .cglbm_stop("Finite intervals require lower < upper.")
  }
  if (!is.null(latent_truth)) {
    latent_truth <- as_numeric_matrix(latent_truth, "latent_truth")
    if (!identical(dim(latent_truth), c(n, d))) .cglbm_stop("latent_truth must match value dimensions.")
  }
  structure(list(
    value = value, status = status, lower = lower, upper = upper,
    n = n, d = d, latent_truth = latent_truth, metadata = metadata
  ), class = "cglbm_data")
}

#' @rdname new_cglbm_data
#' @param y Numeric matrix interpreted as an uncensored Gaussian matrix;
#'   `NA` cells are missing.
#' @export
exact_cglbm_data <- function(y, metadata = list()) {
  y <- as_numeric_matrix(y, "y")
  # This constructor explicitly declares an uncensored numerical matrix.
  # Historical censoring metadata from completion must not change that model.
  metadata$reporting_design <- list(lower_limit = -Inf, upper_limit = Inf)
  new_cglbm_data(y, latent_truth = y, metadata = metadata)
}

#' @rdname new_cglbm_data
#' @param y Numeric matrix of latent or observed values; `NA` cells are
#'   recorded as missing.
#' @param lower_limit,upper_limit Scalar, column-specific vector, or matrix of
#'   fixed detection limits.
#' @export
apply_detection_limits <- function(y, lower_limit = -Inf, upper_limit = Inf, metadata = list()) {
  y <- as_numeric_matrix(y, "y")
  n <- nrow(y); d <- ncol(y)
  lo <- expand_to_matrix(lower_limit, n, d, "lower_limit")
  hi <- expand_to_matrix(upper_limit, n, d, "upper_limit")
  missing <- is.na(y)
  observed <- !missing
  if (any(is.na(lo[observed]) | is.na(hi[observed]) | lo[observed] >= hi[observed])) {
    .cglbm_stop("Every lower detection limit must be below its upper limit.")
  }
  left <- matrix(FALSE, n, d); right <- matrix(FALSE, n, d)
  left[observed] <- y[observed] <= lo[observed]
  right[observed] <- y[observed] >= hi[observed]
  if (any(left & right)) .cglbm_stop("A cell cannot be both left- and right-censored.")
  status <- matrix("exact", n, d); status[left] <- "left"; status[right] <- "right"
  status[missing] <- "missing"
  value <- y; value[left | right] <- NA_real_
  lower <- matrix(NA_real_, n, d); upper <- matrix(NA_real_, n, d)
  upper[left] <- lo[left]; lower[right] <- hi[right]
  metadata$censoring <- list(
    lower_limit = lower_limit, upper_limit = upper_limit,
    left_count = sum(left), right_count = sum(right), exact_count = sum(!(left | right | missing)),
    total_count = length(y), left_rate = mean(left), right_rate = mean(right),
    total_rate = mean(left | right)
  )
  if (any(missing)) {
    metadata$censoring$missing_count <- sum(missing)
    metadata$censoring$missing_rate <- mean(missing)
  }
  # Preserve the entire rule, including limits at realized exact cells.
  metadata$reporting_design <- list(lower_limit = lower_limit, upper_limit = upper_limit)
  new_cglbm_data(value, status, lower, upper, latent_truth = y, metadata = metadata)
}

#' @rdname new_cglbm_data
#' @param x A `cglbm_data` object.
#' @param ... Additional arguments, currently unused.
#' @export
print.cglbm_data <- function(x, ...) {
  counts <- table(factor(x$status, levels = c("exact", "left", "right", "interval", "missing")))
  cat("CGLBM reported-data object\n")
  cat(sprintf("  dimensions: %d x %d\n", x$n, x$d))
  cat(sprintf("  exact=%d, left=%d, right=%d, interval=%d",
              counts[[1]], counts[[2]], counts[[3]], counts[[4]]))
  if (counts[[5]] > 0L) cat(sprintf(", missing=%d", counts[[5]]))
  cat("\n")
  invisible(x)
}

initialization_pseudo_matrix <- function(data, z = NULL, w = NULL) {
  .cglbm_assert(inherits(data, "cglbm_data"), "data must be a cglbm_data object.")
  out <- data$value
  left <- data$status == "left"; right <- data$status == "right"
  interval <- data$status == "interval"
  out[left] <- data$upper[left]
  out[right] <- data$lower[right]
  out[interval] <- data$lower[interval] / 2 + data$upper[interval] / 2
  missing <- data$status == "missing"
  if (any(missing)) {
    # Before a partition is supplied, treat the observed matrix as one block.
    # With initial labels, use only that block's nonmissing pseudo-values:
    # exact values, one-sided endpoints, and finite-interval midpoints.
    if (xor(is.null(z), is.null(w))) .cglbm_stop("Supply both z and w for block initialization.")
    if (is.null(z)) {
      z <- rep.int(1L, data$n); w <- rep.int(1L, data$d)
    }
    if (length(z) != data$n || length(w) != data$d || anyNA(z) || anyNA(w)) {
      .cglbm_stop("Initialization labels must match the reported-data dimensions and have no missing values.")
    }
    for (k in sort(unique(z))) for (l in sort(unique(w))) {
      rows <- z == k; cols <- w == l
      block <- out[rows, cols, drop = FALSE]
      absent <- missing[rows, cols, drop = FALSE]
      if (!any(!absent)) {
        .cglbm_stop("No observed cells in initialization block (%d,%d).", k, l)
      }
      if (any(absent)) {
        block[absent] <- mean(block[!absent])
        out[rows, cols] <- block
      }
    }
  }
  if (any(!is.finite(out))) .cglbm_stop("Initialization pseudo-matrix contains non-finite values.")
  out
}

validate_model <- function(model, n, d) {
  if (!all(c("g", "m", "S", "P") %in% names(model))) .cglbm_stop("model must contain g, m, S, and P.")
  model$g <- as.integer(model$g); model$m <- as.integer(model$m)
  model$S <- as.character(model$S); model$P <- tolower(as.character(model$P))
  .cglbm_assert(model$g >= 1L && model$g <= n, "g must be between 1 and n.")
  .cglbm_assert(model$m >= 1L && model$m <= d, "m must be between 1 and d.")
  .cglbm_assert(model$S %in% c("Skl", "Sk", "Sl", "S0"), "S must be Skl, Sk, Sl, or S0.")
  .cglbm_assert(model$P %in% c("free", "equal"), "P must be free or equal.")
  model
}

variance_group_id <- function(g, m, S) {
  if (S == "Skl") return(matrix(seq_len(g * m), nrow = g, ncol = m))
  if (S == "Sk") return(matrix(rep(seq_len(g), times = m), nrow = g, ncol = m))
  if (S == "Sl") return(matrix(rep(seq_len(m), each = g), nrow = g, ncol = m))
  if (S == "S0") return(matrix(1L, nrow = g, ncol = m))
  .cglbm_stop("Unknown variance structure: %s", S)
}

variance_floor_history_row <- function(phase, iteration, group_id,
                                       raw_variance, sigma2_min) {
  data.frame(
    phase = as.character(phase),
    iteration = as.integer(iteration),
    group_id = as.integer(group_id),
    raw_variance = as.numeric(raw_variance),
    sigma2_min = as.numeric(sigma2_min),
    activated = is.finite(raw_variance) && raw_variance <= sigma2_min,
    stringsAsFactors = FALSE
  )
}

summarize_variance_floor_history <- function(history) {
  if (is.null(history) || !is.data.frame(history) || !nrow(history)) {
    history <- data.frame(
      phase = character(), iteration = integer(), group_id = integer(),
      raw_variance = numeric(), sigma2_min = numeric(), activated = logical(),
      stringsAsFactors = FALSE
    )
  }
  active <- history$activated %in% TRUE
  finite_raw <- history$raw_variance[is.finite(history$raw_variance)]
  list(
    bound_at_initialization = any(active & history$phase == "initialization"),
    bound_during_update = any(active & history$phase == "update"),
    activation_count = as.integer(sum(active)),
    min_unconstrained_variance = if (length(finite_raw)) min(finite_raw) else NA_real_,
    history = history,
    activation_table = history[active, , drop = FALSE]
  )
}

initialize_parameters <- function(data, model, z, w, sigma2_min,
                                  record_variance_floor = FALSE) {
  model <- validate_model(model, data$n, data$d)
  z <- as.integer(z); w <- as.integer(w)
  if (length(z) != data$n || length(w) != data$d) .cglbm_stop("Initialization labels have incorrect lengths.")
  if (anyNA(z) || anyNA(w) || any(z < 1L | z > model$g) || any(w < 1L | w > model$m)) {
    .cglbm_stop("Initialization labels must be integers in 1:g and 1:m.")
  }
  if (length(unique(z)) != model$g || length(unique(w)) != model$m) {
    .cglbm_stop("Initialization labels must contain every declared class.")
  }
  y0 <- initialization_pseudo_matrix(data, z, w)
  missing <- data$status == "missing"
  has_missing <- any(missing)
  mu <- matrix(NA_real_, model$g, model$m)
  weight <- matrix(0, model$g, model$m)
  ss <- matrix(0, model$g, model$m)
  for (k in seq_len(model$g)) for (l in seq_len(model$m)) {
    vals <- y0[z == k, w == l, drop = FALSE]
    # Filled pseudo-values do not become observations, even when computing
    # the initial variance. A block with no observed cells was rejected above.
    if (has_missing) vals <- vals[!missing[z == k, w == l, drop = FALSE]]
    mu[k, l] <- mean(vals); weight[k, l] <- length(vals)
    ss[k, l] <- sum((vals - mu[k, l])^2)
  }
  gid <- variance_group_id(model$g, model$m, model$S)
  var <- matrix(NA_real_, model$g, model$m)
  floor_history <- list()
  for (gr in sort(unique(as.vector(gid)))) {
    idx <- gid == gr
    raw_variance <- sum(ss[idx]) / sum(weight[idx])
    var[idx] <- max(raw_variance, sigma2_min)
    floor_history[[length(floor_history) + 1L]] <- variance_floor_history_row(
      "initialization", 0L, gr, raw_variance, sigma2_min
    )
  }
  if (model$P == "free") {
    pi <- tabulate(z, nbins = model$g) / data$n
    rho <- tabulate(w, nbins = model$m) / data$d
  } else {
    pi <- rep(1 / model$g, model$g); rho <- rep(1 / model$m, model$m)
  }
  out <- list(pi = pi, rho = rho, mu = mu, var = var)
  if (isTRUE(record_variance_floor)) {
    attr(out, "variance_floor_history") <- do.call(rbind, floor_history)
  }
  out
}

conditional_moment_matrices <- function(data, mean, variance) {
  m1 <- matrix(NA_real_, data$n, data$d); m2 <- m1
  exact <- data$status == "exact"; cens <- !exact & data$status != "missing"
  m1[exact] <- data$value[exact]; m2[exact] <- data$value[exact]^2
  if (any(cens)) {
    mom <- truncated_normal_moments(data$lower[cens], data$upper[cens], mean, variance)
    m1[cens] <- mom$first; m2[cens] <- mom$second
  }
  list(first = m1, second = m2)
}

cell_loglik_matrix <- function(data, mean, variance) {
  sd <- sqrt(variance)
  if (!is.finite(sd) || sd <= 0) .cglbm_stop("Variance must be positive and finite.")
  out <- matrix(NA_real_, data$n, data$d)
  exact <- data$status == "exact"; cens <- !exact
  out[exact] <- dnorm(data$value[exact], mean = mean, sd = sd, log = TRUE)
  if (any(cens)) out[cens] <- log_normal_interval_prob(
    data$lower[cens], data$upper[cens], mean = mean, sd = sd
  )
  if (any(!is.finite(out))) .cglbm_stop("A non-finite cell log-likelihood was produced.")
  out
}

.cglbm_r_all_cell_loglik <- function(data, params, model) {
  out <- array(NA_real_, dim = c(data$n, data$d, model$g, model$m))
  for (k in seq_len(model$g)) for (l in seq_len(model$m)) {
    out[, , k, l] <- cell_loglik_matrix(data, params$mu[k, l], params$var[k, l])
  }
  out
}

complete_censoring_in_one_direction <- function(data, weights) {
  weights <- as_numeric_matrix(weights, "weights")
  if (!identical(dim(weights), c(data$n, data$d))) {
    .cglbm_stop("weights must match the reported-data dimensions.")
  }
  active <- is.finite(weights) & weights > 0 & data$status != "missing"
  if (!any(active)) return("")

  active_status <- data$status[active]
  if (all(active_status == "left") && all(is.finite(data$upper[active]))) {
    return("complete left censoring in one direction")
  }
  if (all(active_status == "right") && all(is.finite(data$lower[active]))) {
    return("complete right censoring in one direction")
  }
  ""
}

complete_censoring_at_common_threshold <- function(data, weights) {
  complete_censoring_in_one_direction(data, weights)
}

.cglbm_r_weighted_gaussian_em_step <- function(data, tau, omega, old_params, model,
                                      sigma2_min, effective_weight_min = 0,
                                      iteration = NA_integer_,
                                      cached_complete_censoring_reason = NULL) {
  g <- model$g; m <- model$m
  R <- matrix(0, g, m)
  mu_new <- matrix(NA_real_, g, m)
  block_ss <- matrix(NA_real_, g, m)
  missing <- data$status == "missing"
  has_missing <- any(missing)
  for (k in seq_len(g)) for (l in seq_len(m)) {
    weights <- tcrossprod(tau[, k], omega[, l])
    # Zero the missing weights before the block mass, anchor and active-cell
    # moments are formed. Never multiply missing values or moments by zero.
    if (has_missing) weights[missing] <- 0
    R[k, l] <- sum(weights)
    if (!is.finite(R[k, l]) || R[k, l] <= effective_weight_min) {
      return(list(ok = FALSE, reason = sprintf("effective block weight too small at (%d,%d)", k, l)))
    }
    boundary_reason <- if (is.null(cached_complete_censoring_reason)) {
      complete_censoring_in_one_direction(data, weights)
    } else {
      cached_complete_censoring_reason
    }
    if (nzchar(boundary_reason)) {
      return(list(
        ok = FALSE,
        reason = sprintf("%s in block (%d,%d)", boundary_reason, k, l)
      ))
    }
    active <- weights > 0
    active_weights <- weights[active]
    active_status <- data$status[active]
    active_value <- data$value[active]
    active_lower <- data$lower[active]
    active_upper <- data$upper[active]
    conditional_mean <- numeric(length(active_weights))
    conditional_variance <- numeric(length(active_weights))

    exact_active <- active_status == "exact"
    cens_active <- !exact_active
    if (any(exact_active)) {
      conditional_mean[exact_active] <- active_value[exact_active]
    }
    if (any(cens_active)) {
      mom <- truncated_normal_moments(
        active_lower[cens_active], active_upper[cens_active],
        old_params$mu[k, l], old_params$var[k, l]
      )
      conditional_mean[cens_active] <- mom$first
      conditional_variance[cens_active] <- mom$variance
    }
    if (any(!is.finite(conditional_mean)) ||
        any(!is.finite(conditional_variance)) ||
        any(conditional_variance < 0)) {
      return(list(ok = FALSE, reason = "invalid conditional moments"))
    }

    anchor <- conditional_mean[which.max(active_weights)]
    mu_new[k, l] <- anchor +
      sum(active_weights * (conditional_mean - anchor)) / R[k, l]
    centered_mean <- conditional_mean - mu_new[k, l]
    # This is algebraically identical to the raw second-moment update, but it
    # avoids subtracting quantities of order mu^2 when the data are translated.
    block_ss[k, l] <- sum(
      active_weights * (conditional_variance + centered_mean^2)
    )
  }
  if (any(!is.finite(block_ss)) || any(block_ss < 0)) {
    return(list(ok = FALSE, reason = "invalid expected residual sum of squares"))
  }
  gid <- variance_group_id(g, m, model$S)
  var_new <- matrix(NA_real_, g, m)
  floor_history <- list()
  for (gr in sort(unique(as.vector(gid)))) {
    idx <- gid == gr; denom <- sum(R[idx])
    if (!is.finite(denom) || denom <= effective_weight_min) {
      return(list(ok = FALSE, reason = sprintf("variance-group weight too small for group %d", gr)))
    }
    raw_variance <- sum(block_ss[idx]) / denom
    var_new[idx] <- max(raw_variance, sigma2_min)
    floor_history[[length(floor_history) + 1L]] <- variance_floor_history_row(
      "update", iteration, gr, raw_variance, sigma2_min
    )
  }
  list(
    ok = TRUE, mu = mu_new, var = var_new, block_weight = R, block_ss = block_ss,
    variance_floor_history = do.call(rbind, floor_history)
  )
}

# -----------------------------------------------------------------------------
# Fisher information for block means
# -----------------------------------------------------------------------------
# Expected information is conditional on the FIXED reporting rule and fitted
# memberships. A realized exact/left/right status does not select a different
# expected-information formula under that rule. See R/README.md for the input
# contract and the two-detection-limit extension of the manuscript's lemma.

.cglbm_left_information_components <- function(alpha) {
  alpha <- as.numeric(alpha)
  if (anyNA(alpha)) .cglbm_stop("Standardized censoring thresholds cannot be NA.")
  A <- B <- C <- numeric(length(alpha))
  log_phi <- dnorm(alpha, log = TRUE)
  phi <- exp(log_phi)
  exact_limit <- phi == 0 & alpha < 0
  A[exact_limit] <- 1
  C[exact_limit] <- 2
  middle <- phi > 0
  if (any(middle)) {
    a <- alpha[middle]
    r <- exp(log_phi[middle] - pnorm(a, log.p = TRUE))
    survival <- pnorm(a, lower.tail = FALSE)
    A[middle] <- survival + phi[middle] * (a + r)
    B[middle] <- phi[middle] * (a^2 + 1 + a * r)
    C[middle] <- 2 * survival + phi[middle] * (a^3 + a + a^2 * r)
  }
  tolerance <- 100 * .Machine$double.eps
  A[A < 0 & A > -tolerance] <- 0
  B[B < 0 & B > -tolerance] <- 0
  C[C < 0 & C > -tolerance] <- 0
  A[A > 1 & A < 1 + tolerance] <- 1
  C[C > 2 & C < 2 + tolerance] <- 2
  bad <- !is.finite(A) | !is.finite(B) | !is.finite(C) |
    A < 0 | A > 1 | B < 0 | C < 0 | C > 2
  # Numerical unavailability of information is not a failure of the fit.
  A[bad] <- B[bad] <- C[bad] <- NA_real_
  list(A = A, B = B, C = C)
}

.cglbm_reporting_design <- function(data) {
  n <- data$n; d <- data$d
  reason <- matrix("reporting_design_unavailable", n, d)
  lo <- hi <- matrix(NA_real_, n, d)
  metadata <- data$metadata
  spec <- if (is.list(metadata)) metadata$reporting_design else NULL
  source <- "reporting_design"
  if (is.null(spec)) {
    spec <- if (is.list(metadata)) metadata$censoring else NULL
    source <- "censoring_metadata"
  }
  if (is.list(spec) && !is.null(spec$lower_limit) && !is.null(spec$upper_limit)) {
    # Catch only optional metadata conversion/dimension errors here. Invalid
    # fitted parameters and programming errors elsewhere are not swallowed.
    limits <- tryCatch({
      if (!is.numeric(spec$lower_limit) || !is.numeric(spec$upper_limit)) {
        .cglbm_stop("Reporting limits must be numeric.")
      }
      list(lower = expand_to_matrix(spec$lower_limit, n, d, "lower_limit"),
           upper = expand_to_matrix(spec$upper_limit, n, d, "upper_limit"))
    }, error = function(e) NULL)
    if (!is.null(limits)) {
      lo <- limits$lower; hi <- limits$upper
      known <- !is.na(lo) & !is.na(hi) & lo < hi
      reason[known] <- "ok"
      same_endpoint <- function(x, y) {
        is.finite(x) & is.finite(y) &
          abs(x - y) <= 64 * .Machine$double.eps * pmax(1, abs(x), abs(y))
      }
      # The rule must be on the same analysis scale as the reports. Do not
      # use limits known to contradict centered/subset reports. A check cannot
      # establish that otherwise compatible caller-supplied limits are correct.
      exact <- known & data$status == "exact"
      left <- known & data$status == "left"
      right <- known & data$status == "right"
      reason[exact] <- ifelse(
        (data$value[exact] >= lo[exact] | same_endpoint(data$value[exact], lo[exact])) &
          (data$value[exact] <= hi[exact] | same_endpoint(data$value[exact], hi[exact])),
        "ok", "reporting_design_inconsistent")
      reason[left] <- ifelse(same_endpoint(data$upper[left], lo[left]),
                            "ok", "reporting_design_inconsistent")
      reason[right] <- ifelse(same_endpoint(data$lower[right], hi[right]),
                             "ok", "reporting_design_inconsistent")
    }
  }
  # A finite interval-valued report is not the exact-between-two-limits rule.
  reason[data$status == "interval"] <- "finite_interval_unsupported"
  reason[data$status == "missing"] <- "missing"
  list(lower = lo, upper = hi, reason = reason, source = source)
}

.cglbm_exact_region_information <- function(a, b) {
  # Stable fallback for cancellation in exact-region score integrals. Map a
  # narrow interval to [-1,1]; this does not infer latent values or alter EM.
  center <- a / 2 + b / 2
  half_width <- b / 2 - a / 2
  out <- numeric(3L)
  for (component in seq_len(3L)) {
    integral <- stats::integrate(function(t) {
      z <- center + half_width * t
      density <- dnorm(z)
      value <- numeric(length(z))
      active <- density > 0
      x <- z[active]
      polynomial <- switch(component, x^2, x * (x^2 - 1), (x^2 - 1)^2)
      value[active] <- half_width * density[active] * polynomial
      value
    }, lower = -1, upper = 1, subdivisions = 200L,
    rel.tol = 1e-11, abs.tol = 1e-14, stop.on.error = FALSE)
    if (!identical(integral$message, "OK") || !is.finite(integral$value)) {
      return(rep(NA_real_, 3L))
    }
    out[component] <- integral$value
  }
  out
}

.cglbm_detection_information_components <- function(a, b) {
  .cglbm_assert(length(a) == length(b) && !anyNA(a) && !anyNA(b),
                "Standardized detection limits are incompatible.")
  a <- as.numeric(a); b <- as.numeric(b)
  A <- B <- C <- rep(NA_real_, length(a))
  exact <- a == -Inf & b == Inf
  A[exact] <- 1; B[exact] <- 0; C[exact] <- 2
  left <- is.finite(a) & b == Inf
  if (any(left)) {
    one <- .cglbm_left_information_components(a[left])
    A[left] <- one$A; B[left] <- one$B; C[left] <- one$C
  }
  right <- a == -Inf & is.finite(b)
  if (any(right)) {
    one <- .cglbm_left_information_components(-b[right])
    A[right] <- one$A; B[right] <- -one$B; C[right] <- one$C
  }
  # Both endpoints can overflow to the same signed infinity when a finite
  # block parameter is far beyond both limits. The limiting information is 0.
  saturated <- (a == Inf & b == Inf) | (a == -Inf & b == -Inf)
  A[saturated] <- B[saturated] <- C[saturated] <- 0
  both <- is.finite(a) & is.finite(b) & a < b
  if (any(both)) {
    x <- a[both]; y <- b[both]
    px <- dnorm(x); py <- dnorm(y)
    # Exact-region mass evaluated in the stable tail.
    use_right_tail <- x >= 0
    large <- pnorm(y, log.p = TRUE)
    small <- pnorm(x, log.p = TRUE)
    large[use_right_tail] <- pnorm(x[use_right_tail], lower.tail = FALSE, log.p = TRUE)
    small[use_right_tail] <- pnorm(y[use_right_tail], lower.tail = FALSE, log.p = TRUE)
    q <- numeric(length(x))
    finite_mass <- is.finite(large)
    q[finite_mass] <- exp(large[finite_mass]) *
      (-expm1(small[finite_mass] - large[finite_mass]))
    # Avoid Inf*0 in remote tails. Where phi is positive, |z| is small enough
    # for these low-degree products to be representable.
    x1 <- x2 <- x3 <- y1 <- y2 <- y3 <- numeric(length(x))
    ix <- px > 0; iy <- py > 0
    x1[ix] <- x[ix] * px[ix]; x2[ix] <- (x[ix]^2 + 1) * px[ix]
    x3[ix] <- (x[ix]^3 + x[ix]) * px[ix]
    y1[iy] <- y[iy] * py[iy]; y2[iy] <- (y[iy]^2 + 1) * py[iy]
    y3[iy] <- (y[iy]^3 + y[iy]) * py[iy]
    exact_A <- q + x1 - y1
    exact_B <- x2 - y2
    exact_C <- 2 * q + x3 - y3
    narrow <- (y - x) < 1e-4 * (1 + pmin(abs(x), abs(y)))
    cancellation <- exact_A < 1e-8 * (q + abs(x1) + abs(y1)) |
      exact_C < 1e-8 * (2 * q + abs(x3) + abs(y3))
    needs_quadrature <- narrow | cancellation |
      !is.finite(exact_A) | !is.finite(exact_B) | !is.finite(exact_C)
    for (j in which(needs_quadrature)) {
      integrated <- .cglbm_exact_region_information(x[j], y[j])
      exact_A[j] <- integrated[1L]
      exact_B[j] <- integrated[2L]
      exact_C[j] <- integrated[3L]
    }
    # Add the two censoring atoms to the exact-region score integrals.
    atom_left <- atom_right <- numeric(length(x))
    atom_left[ix] <- exp(2 * dnorm(x[ix], log = TRUE) - pnorm(x[ix], log.p = TRUE))
    atom_right[iy] <- exp(2 * dnorm(y[iy], log = TRUE) -
                            pnorm(y[iy], lower.tail = FALSE, log.p = TRUE))
    left_B <- left_C <- right_B <- right_C <- numeric(length(x))
    left_B[ix] <- x[ix] * atom_left[ix]
    left_C[ix] <- x[ix]^2 * atom_left[ix]
    right_B[iy] <- y[iy] * atom_right[iy]
    right_C[iy] <- y[iy]^2 * atom_right[iy]
    A[both] <- exact_A + atom_left + atom_right
    B[both] <- exact_B + left_B + right_B
    C[both] <- exact_C + left_C + right_C
  }
  tolerance <- 100 * .Machine$double.eps
  A[!is.na(A) & A < 0 & A > -tolerance] <- 0
  C[!is.na(C) & C < 0 & C > -tolerance] <- 0
  A[!is.na(A) & A > 1 & A < 1 + tolerance] <- 1
  C[!is.na(C) & C > 2 & C < 2 + tolerance] <- 2
  bad <- !is.finite(A) | !is.finite(B) | !is.finite(C) |
    A < 0 | A > 1 | C < 0 | C > 2
  bad[is.na(bad)] <- TRUE
  A[bad] <- B[bad] <- C[bad] <- NA_real_
  list(A = A, B = B, C = C)
}

.cglbm_cell_information <- function(data, mean, variance, design = NULL) {
  .cglbm_assert(inherits(data, "cglbm_data"), "data must be a cglbm_data object.")
  .cglbm_assert(length(mean) == 1L && is.finite(mean), "The fitted block mean must be finite.")
  .cglbm_assert(length(variance) == 1L && is.finite(variance) && variance > 0,
                "The fitted block variance must be positive and finite.")
  if (is.null(design)) design <- .cglbm_reporting_design(data)
  A <- B <- C <- matrix(NA_real_, data$n, data$d)
  reason <- design$reason
  known <- reason == "ok"
  if (any(known)) {
    components <- .cglbm_detection_information_components(
      (design$lower[known] - mean) / sqrt(variance),
      (design$upper[known] - mean) / sqrt(variance))
    A[known] <- components$A; B[known] <- components$B; C[known] <- components$C
    bad <- known & (!is.finite(A) | !is.finite(B) | !is.finite(C))
    reason[bad] <- "numerical_information_unavailable"
  }
  missing <- data$status == "missing"
  A[missing] <- B[missing] <- C[missing] <- 0
  # Keep the dimensionless components for a unit-scaled 2x2 inversion.
  list(A = A, B = B, C = C,
       I11 = A / variance, I12 = (B / variance) / (2 * sqrt(variance)),
       I22 = ((C / variance) / variance) / 4, reason = reason)
}

.cglbm_mean_information_se <- function(A, B, C, variance) {
  # I = D J D, D = diag(1/sigma, 1/(2*sigma^2)), J = [A B; B C].
  # This is sqrt((I^{-1})[1,1]), including the variance nuisance parameter.
  # Equilibrating J by its diagonal avoids rejecting a valid matrix merely
  # because the mean and variance parameters have different physical units.
  if (!all(is.finite(c(A, B, C))) || A <= 0 || C <= 0) {
    return(list(se = NA_real_, reason = "singular_information"))
  }
  correlation <- (B / sqrt(A)) / sqrt(C)
  determinant <- (1 - correlation) * (1 + correlation)
  if (!is.finite(determinant) || determinant <= 100 * .Machine$double.eps) {
    return(list(se = NA_real_, reason = "singular_information"))
  }
  se <- (sqrt(variance) / sqrt(A)) / sqrt(determinant)
  if (!is.finite(se) || se <= 0) {
    return(list(se = NA_real_, reason = "numerical_information_unavailable"))
  }
  list(se = se, reason = "ok")
}

# Independent mean coordinates coupled to ONE shared variance. The matrix is
# J = [diag(A), B; t(B), sum(C)], with I = D J D and
# D = diag(rep(1/sqrt(variance), length(A)), 1/(2*variance)).
# All components have already been summed with the fitted gamma weights.
.cglbm_shared_information_inverse <- function(A, B, C, variance) {
  q <- length(A)
  unavailable <- function(reason, information = matrix(NA_real_, q + 1L, q + 1L)) {
    list(information = information, covariance = matrix(NA_real_, q + 1L, q + 1L),
         se = rep(NA_real_, q), status = reason)
  }
  if (!q || length(B) != q || length(C) != q ||
      !all(is.finite(c(A, B, C))) || !is.finite(sum(C))) {
    return(unavailable("numerical_information_unavailable"))
  }
  J <- matrix(0, q + 1L, q + 1L)
  diag(J) <- c(A, sum(C))
  J[seq_len(q), q + 1L] <- B
  J[q + 1L, seq_len(q)] <- B
  # Divide successively instead of forming sigma^4 or other extreme products.
  unit <- c(rep(sqrt(variance), q), 2 * variance)
  I <- sweep(sweep(J, 1L, unit, "/"), 2L, unit, "/")
  if (any(!is.finite(I)) || any(!is.finite(unit))) {
    return(unavailable("numerical_information_unavailable"))
  }
  if (any(A <= 0) || any(C < 0) || sum(C) <= 0 || any(diag(I) <= 0)) {
    return(unavailable("singular_information", I))
  }
  root_diag <- sqrt(diag(J))
  correlation <- sweep(sweep(J, 1L, root_diag, "/"), 2L, root_diag, "/")
  # Equilibration removes parameter units. The Schur complement must be
  # positive; no pseudoinverse or independent 2x2 inversion is substituted.
  links <- correlation[seq_len(q), q + 1L]
  schur <- 1 - sum(links^2)
  if (!is.finite(schur) || schur <= 100 * .Machine$double.eps) {
    return(unavailable("singular_information", I))
  }
  # chol/Cholesky inversion may fail on a numerically singular matrix. This
  # narrow handler covers only that numerical operation, not fitted inputs.
  inverse <- tryCatch(chol2inv(chol(correlation)), error = function(e) NULL)
  if (is.null(inverse)) return(unavailable("singular_information", I))
  covariance <- sweep(sweep(inverse, 1L, root_diag, "/"), 2L, root_diag, "/")
  covariance <- sweep(sweep(covariance, 1L, unit, "*"), 2L, unit, "*")
  if (any(!is.finite(covariance)) || any(diag(covariance) <= 0)) {
    return(unavailable("numerical_information_unavailable", I))
  }
  list(information = I, covariance = covariance,
       se = sqrt(diag(covariance)[seq_len(q)]), status = "ok")
}

# Rebuild joint groups from block totals, also after label permutations. Keep
# block_information as LOCAL 2x2 contributions, not a shared-model covariance.
.cglbm_refresh_shared_information <- function(information) {
  S <- information$variance_structure
  v <- information$block_variance
  g <- nrow(v); m <- ncol(v)
  gid <- variance_group_id(g, m, S)
  information$variance_group_id <- gid
  information$information_status <- information$block_information_status
  information$mean_standard_error[,] <- NA_real_
  groups <- vector("list", max(gid))
  for (gr in seq_len(max(gid))) {
    blocks <- arrayInd(which(gid == gr), dim(gid))
    q <- nrow(blocks)
    variance <- v[blocks[1L, 1L], blocks[1L, 2L]]
    # Sharing is an exact constraint in the EM update. A different variance
    # here is a malformed fitted object, not an unavailable information case.
    if (any(v[blocks] != variance)) {
      .cglbm_stop("Fitted block variances violate the declared shared variance structure.")
    }
    components <- matrix(NA_real_, q, 3L)
    for (h in seq_len(q)) {
      components[h, ] <- information$block_dimensionless[
        blocks[h, 1L], blocks[h, 2L], ]
    }
    if (any(!is.finite(components))) {
      answer <- list(information = matrix(NA_real_, q + 1L, q + 1L),
                     covariance = matrix(NA_real_, q + 1L, q + 1L),
                     se = rep(NA_real_, q), status = "shared_variance_group_unavailable")
      # Do not discard a known block's per-cell information or ESS just
      # because another member makes their JOINT covariance unavailable.
      usable <- information$information_status[blocks] %in% c("ok", "singular_information")
      status <- information$information_status[blocks]
      status[usable] <- answer$status
      information$information_status[blocks] <- status
    } else {
      answer <- .cglbm_shared_information_inverse(
        components[, 1L], components[, 2L], components[, 3L], variance)
      information$information_status[blocks] <- answer$status
    }
    information$mean_standard_error[blocks] <- answer$se
    parameter_names <- c(sprintf("mu[%d,%d]", blocks[, 1L], blocks[, 2L]),
                          sprintf("variance[%d]", gr))
    dimnames(answer$information) <- dimnames(answer$covariance) <-
      list(parameter_names, parameter_names)
    groups[[gr]] <- list(blocks = blocks, variance = variance,
                         information = answer$information, covariance = answer$covariance,
                         status = answer$status)
  }
  information$variance_groups <- groups
  information
}

.cglbm_attach_information <- function(fit, data) {
  .cglbm_assert(inherits(fit, "cglbm_fit"), "fit must be a cglbm_fit object.")
  .cglbm_assert(inherits(data, "cglbm_data"), "data must be a cglbm_data object.")
  .cglbm_assert(length(fit$model$S) == 1L &&
                  fit$model$S %in% c("Skl", "Sk", "Sl", "S0"),
                "S must be Skl, Sk, Sl, or S0.")
  shared <- !identical(fit$model$S, "Skl")
  g <- fit$model$g; m <- fit$model$m
  .cglbm_assert(identical(dim(fit$tau), c(data$n, g)) &&
                  identical(dim(fit$omega), c(data$d, m)) &&
                  length(fit$z) == data$n && length(fit$w) == data$d,
                "The fitted memberships do not match the reported-data dimensions.")
  .cglbm_assert(identical(dim(fit$params$mu), c(g, m)) &&
                  identical(dim(fit$params$var), c(g, m)),
                "The fitted block parameters do not match the declared model.")
  .cglbm_assert(all(is.finite(fit$params$mu)), "The fitted block means must be finite.")
  .cglbm_assert(all(is.finite(fit$params$var)) && all(fit$params$var > 0),
                "The fitted block variance must be positive and finite.")
  .cglbm_assert(all(is.finite(fit$tau)) && all(fit$tau >= 0) &&
                  all(is.finite(fit$omega)) && all(fit$omega >= 0),
                "The fitted membership weights must be finite and nonnegative.")
  design <- .cglbm_reporting_design(data)
  per_cell_A <- array(NA_real_, dim = c(data$n, data$d, g, m))
  block_information <- array(NA_real_, dim = c(g, m, 2L, 2L))
  block_dimensionless <- if (shared) array(NA_real_, c(g, m, 3L)) else NULL
  mean_standard_error <- effective_sample_size <- matrix(NA_real_, g, m)
  effective_block_weight <- matrix(0, g, m)
  exact_count <- left_censored_count <- right_censored_count <-
    nominal_cell_count <- matrix(0L, g, m)
  interval_block <- matrix(FALSE, g, m)
  information_status <- matrix("ok", g, m)
  observed <- data$status != "missing"
  interval <- data$status == "interval"
  for (k in seq_len(g)) for (l in seq_len(m)) {
    cell <- .cglbm_cell_information(data, fit$params$mu[k, l], fit$params$var[k, l], design)
    per_cell_A[, , k, l] <- cell$A
    weights <- tcrossprod(fit$tau[, k], fit$omega[, l])
    active <- weights > 0 & observed
    effective_block_weight[k, l] <- sum(weights[observed])
    map_block <- tcrossprod(as.numeric(fit$z == k), as.numeric(fit$w == l)) > 0
    exact_count[k, l] <- sum(map_block & data$status == "exact")
    left_censored_count[k, l] <- sum(map_block & data$status == "left")
    right_censored_count[k, l] <- sum(map_block & data$status == "right")
    nominal_cell_count[k, l] <- sum(map_block)
    interval_block[k, l] <- any(weights[interval] > 0)
    if (interval_block[k, l]) {
      information_status[k, l] <- "finite_interval_unsupported"
      next
    }
    if (any(cell$reason[active] != "ok")) {
      information_status[k, l] <- cell$reason[active][which(cell$reason[active] != "ok")[1L]]
      next
    }
    A <- sum(weights[active] * cell$A[active])
    B <- sum(weights[active] * cell$B[active])
    C <- sum(weights[active] * cell$C[active])
    if (shared) block_dimensionless[k, l, ] <- c(A, B, C)
    effective_sample_size[k, l] <- A
    variance <- fit$params$var[k, l]
    I11 <- A / variance
    I12 <- (B / variance) / (2 * sqrt(variance))
    I22 <- ((C / variance) / variance) / 4
    information_matrix <- matrix(c(I11, I12, I12, I22), 2L, 2L)
    if (any(!is.finite(information_matrix))) {
      information_status[k, l] <- "numerical_information_unavailable"
      next
    }
    block_information[k, l, , ] <- information_matrix
    if (I11 <= 0 || I22 <= 0) {
      information_status[k, l] <- "singular_information"
      next
    }
    uncertainty <- .cglbm_mean_information_se(A, B, C, variance)
    mean_standard_error[k, l] <- uncertainty$se
    information_status[k, l] <- uncertainty$reason
  }
  fit$information <- list(
    method = "expected_fisher_fixed_reporting_rule",
    per_cell_A = per_cell_A, block_information = block_information,
    mean_standard_error = mean_standard_error, effective_block_weight = effective_block_weight,
    exact_count = exact_count, left_censored_count = left_censored_count,
    right_censored_count = right_censored_count, nominal_cell_count = nominal_cell_count,
    effective_sample_size = effective_sample_size, interval_block = interval_block,
    information_status = information_status
  )
  if (shared) {
    fit$information$variance_structure <- fit$model$S
    fit$information$block_dimensionless <- block_dimensionless
    fit$information$block_variance <- fit$params$var
    fit$information$block_information_status <- information_status
    fit$information <- .cglbm_refresh_shared_information(fit$information)
  }
  fit$information_unavailable_reason <- NULL
  fit
}

.cglbm_reorder_information <- function(information, row_order, col_order) {
  if (is.null(information)) return(NULL)
  out <- information
  if (!is.null(information$per_cell_A)) {
    out$per_cell_A <- information$per_cell_A[, , row_order, col_order, drop = FALSE]
  }
  if (!is.null(information$block_information)) {
    out$block_information <- information$block_information[row_order, col_order, , , drop = FALSE]
  }
  for (field in c("mean_standard_error", "effective_block_weight", "exact_count",
                   "left_censored_count", "right_censored_count", "nominal_cell_count",
                   "effective_sample_size", "interval_block", "information_status")) {
    if (!is.null(information[[field]])) {
      out[[field]] <- information[[field]][row_order, col_order, drop = FALSE]
    }
  }
  if (!is.null(information$variance_structure)) {
    out$block_dimensionless <- information$block_dimensionless[row_order, col_order, , drop = FALSE]
    out$block_variance <- information$block_variance[row_order, col_order, drop = FALSE]
    out$block_information_status <- information$block_information_status[row_order, col_order, drop = FALSE]
    out <- .cglbm_refresh_shared_information(out)
  }
  out
}

# -----------------------------------------------------------------------------
# Membership configurations
# -----------------------------------------------------------------------------
# The first letter refers to rows and the second to columns.
#   HH = hard rows, hard columns
#   HS = hard rows, soft columns
#   SH = soft rows, hard columns
#   SS = soft rows, soft columns
#
# All four use the same censored Gaussian LBM objective and differ only in the
# admissible row and column membership domains.

normalize_algorithm <- function(algorithm) {
  if (length(algorithm) != 1L || is.na(algorithm)) {
    .cglbm_stop("algorithm must be one of HH, HS, SH, or SS.")
  }
  key <- gsub("[^a-z]", "", tolower(as.character(algorithm)))
  aliases <- c(
    hh = "HH", hardhard = "HH", classification = "HH", cgem = "HH",
    hs = "HS", hardsoft = "HS",
    sh = "SH", softhard = "SH", nadifgem = "SH",
    ss = "SS", softsoft = "SS", variational = "SS", vem = "SS"
  )
  if (!key %in% names(aliases)) {
    .cglbm_stop("Unknown algorithm '%s'; use HH, HS, SH, or SS.", algorithm)
  }
  unname(aliases[[key]])
}

algorithm_spec <- function(algorithm) {
  code <- normalize_algorithm(algorithm)
  list(
    code = code,
    row_mode = substr(code, 1L, 1L),
    col_mode = substr(code, 2L, 2L),
    label = switch(code,
      HH = "hard rows / hard columns",
      HS = "hard rows / soft columns",
      SH = "soft rows / hard columns",
      SS = "soft rows / soft columns"
    )
  )
}

membership_entropy <- function(x) {
  # Returns positive Shannon entropy.  It is exactly zero for one-hot rows.
  -xlogx_sum(x)
}

.cglbm_r_membership_objective <- function(tau, omega, params, logL, model) {
  tau <- as.matrix(tau); omega <- as.matrix(omega)
  if (!identical(dim(tau), c(dim(logL)[1], model$g))) {
    .cglbm_stop("tau has incompatible dimensions.")
  }
  if (!identical(dim(omega), c(dim(logL)[2], model$m))) {
    .cglbm_stop("omega has incompatible dimensions.")
  }
  if (any(params$pi < 0) || any(params$rho < 0)) return(-Inf)

  row_term <- 0
  for (k in seq_len(model$g)) {
    mass <- sum(tau[, k])
    if (params$pi[k] == 0 && mass > 0) return(-Inf)
    if (params$pi[k] > 0) row_term <- row_term + mass * log(params$pi[k])
  }
  col_term <- 0
  for (l in seq_len(model$m)) {
    mass <- sum(omega[, l])
    if (params$rho[l] == 0 && mass > 0) return(-Inf)
    if (params$rho[l] > 0) col_term <- col_term + mass * log(params$rho[l])
  }

  n <- dim(logL)[1]; d <- dim(logL)[2]
  data_term <- 0
  for (k in seq_len(model$g)) for (l in seq_len(model$m)) {
    ll_kl <- matrix(logL[, , k, l, drop = FALSE], nrow = n, ncol = d)
    data_term <- data_term + sum(tau[, k] * as.vector(ll_kl %*% omega[, l]))
  }
  as.numeric(row_term + col_term + data_term +
               membership_entropy(tau) + membership_entropy(omega))
}

classification_objective <- function(z, w, params, logL, model) {
  membership_objective(one_hot(z, model$g), one_hot(w, model$m), params, logL, model)
}

variational_objective <- function(tau, omega, params, logL, model) {
  membership_objective(tau, omega, params, logL, model)
}

.cglbm_r_row_membership_scores <- function(omega, params, logL, model) {
  n <- dim(logL)[1]; d <- dim(logL)[2]
  scores <- matrix(NA_real_, nrow = n, ncol = model$g)
  for (k in seq_len(model$g)) {
    scores[, k] <- log(params$pi[k])
    for (l in seq_len(model$m)) {
      ll_kl <- matrix(logL[, , k, l, drop = FALSE], nrow = n, ncol = d)
      scores[, k] <- scores[, k] + as.vector(ll_kl %*% omega[, l])
    }
  }
  scores
}

.cglbm_r_col_membership_scores <- function(tau, params, logL, model) {
  n <- dim(logL)[1]; d <- dim(logL)[2]
  scores <- matrix(NA_real_, nrow = d, ncol = model$m)
  for (l in seq_len(model$m)) {
    scores[, l] <- log(params$rho[l])
    for (k in seq_len(model$g)) {
      ll_kl <- matrix(logL[, , k, l, drop = FALSE], nrow = n, ncol = d)
      scores[, l] <- scores[, l] + as.vector(t(ll_kl) %*% tau[, k])
    }
  }
  scores
}

update_membership_from_scores <- function(scores, mode = c("H", "S")) {
  mode <- match.arg(mode)
  if (mode == "S") return(row_softmax(scores))
  one_hot(max.col(scores, ties.method = "first"), ncol(scores))
}

classification_row_update <- function(w, params, logL, model) {
  max.col(row_membership_scores(one_hot(w, model$m), params, logL, model),
          ties.method = "first")
}

classification_col_update <- function(z, params, logL, model) {
  max.col(col_membership_scores(one_hot(z, model$g), params, logL, model),
          ties.method = "first")
}

variational_tau_update <- function(omega, params, logL, model) {
  row_softmax(row_membership_scores(omega, params, logL, model))
}

variational_omega_update <- function(tau, params, logL, model) {
  row_softmax(col_membership_scores(tau, params, logL, model))
}

check_nonempty_hard_labels <- function(z, w, model) {
  length(unique(z)) == model$g && length(unique(w)) == model$m
}

check_nonempty_dimension <- function(membership, n_classes) {
  labels <- max.col(membership, ties.method = "first")
  length(unique(labels)) == n_classes
}

validate_membership_masses <- function(tau, omega, spec, effective_weight_min) {
  row_mass <- colSums(tau)
  col_mass <- colSums(omega)
  if (spec$row_mode == "H" && any(row_mass < 1 - 1e-12)) {
    return(list(ok = FALSE, reason = "empty hard row class"))
  }
  if (spec$col_mode == "H" && any(col_mass < 1 - 1e-12)) {
    return(list(ok = FALSE, reason = "empty hard column class"))
  }
  if (any(!is.finite(row_mass)) || any(row_mass <= effective_weight_min)) {
    return(list(ok = FALSE, reason = "effective row-class mass too small"))
  }
  if (any(!is.finite(col_mass)) || any(col_mass <= effective_weight_min)) {
    return(list(ok = FALSE, reason = "effective column-class mass too small"))
  }
  list(ok = TRUE, row_mass = row_mass, col_mass = col_mass)
}

fit_membership_start <- function(data, model, start, algorithm, max_iter, tol,
                                 sigma2_min, effective_weight_min,
                                 monotone_tolerance,
                                 cached_complete_censoring_reason = NULL) {
  spec <- algorithm_spec(algorithm)
  z0 <- as.integer(start$z); w0 <- as.integer(start$w)
  if (!check_nonempty_hard_labels(z0, w0, model)) {
    return(list(ok = FALSE, reason = "initial hard labels are empty"))
  }

  tau <- one_hot(z0, model$g)
  omega <- one_hot(w0, model$m)
  params <- initialize_parameters(
    data, model, z0, w0, sigma2_min, record_variance_floor = TRUE
  )
  variance_floor_history <- attr(params, "variance_floor_history")
  attr(params, "variance_floor_history") <- NULL
  floor_diagnostics <- function() {
    summarize_variance_floor_history(variance_floor_history)
  }
  logL <- all_cell_loglik(data, params, model)
  objective <- membership_objective(tau, omega, params, logL, model)
  trace <- objective
  substeps <- list()
  membership_deltas <- numeric()
  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    old_tau <- tau; old_omega <- omega
    old_params <- params; old_logL <- logL; old_obj <- objective

    row_scores <- row_membership_scores(old_omega, old_params, old_logL, model)
    tau_new <- update_membership_from_scores(row_scores, spec$row_mode)
    mass_check <- validate_membership_masses(
      tau_new, old_omega, spec, effective_weight_min
    )
    if (!isTRUE(mass_check$ok)) {
      return(list(ok = FALSE, reason = mass_check$reason, iteration = iter,
                  variance_floor_diagnostics = floor_diagnostics()))
    }
    obj_row <- membership_objective(tau_new, old_omega, old_params, old_logL, model)

    col_scores <- col_membership_scores(tau_new, old_params, old_logL, model)
    omega_new <- update_membership_from_scores(col_scores, spec$col_mode)
    mass_check <- validate_membership_masses(
      tau_new, omega_new, spec, effective_weight_min
    )
    if (!isTRUE(mass_check$ok)) {
      return(list(ok = FALSE, reason = mass_check$reason, iteration = iter,
                  variance_floor_diagnostics = floor_diagnostics()))
    }
    obj_col <- membership_objective(tau_new, omega_new, old_params, old_logL, model)

    params_mid <- old_params
    if (model$P == "free") {
      params_mid$pi <- colMeans(tau_new)
      params_mid$rho <- colMeans(omega_new)
    } else {
      params_mid$pi <- rep(1 / model$g, model$g)
      params_mid$rho <- rep(1 / model$m, model$m)
    }
    if (any(!is.finite(params_mid$pi)) || any(params_mid$pi <= 0) ||
        any(!is.finite(params_mid$rho)) || any(params_mid$rho <= 0)) {
      return(list(ok = FALSE, reason = "zero or non-finite class proportion", iteration = iter,
                  variance_floor_diagnostics = floor_diagnostics()))
    }
    obj_prop <- membership_objective(tau_new, omega_new, params_mid, old_logL, model)

    upd <- weighted_gaussian_em_step(
      data = data,
      tau = tau_new,
      omega = omega_new,
      old_params = old_params,
      model = model,
      sigma2_min = sigma2_min,
      effective_weight_min = effective_weight_min,
      iteration = iter,
      cached_complete_censoring_reason = cached_complete_censoring_reason
    )
    if (!isTRUE(upd$ok)) {
      return(list(ok = FALSE, reason = upd$reason, iteration = iter,
                  variance_floor_diagnostics = floor_diagnostics()))
    }
    variance_floor_history <- rbind(
      variance_floor_history, upd$variance_floor_history
    )
    params_new <- params_mid
    params_new$mu <- upd$mu
    params_new$var <- upd$var
    logL_new <- all_cell_loglik(data, params_new, model)
    obj_new <- membership_objective(tau_new, omega_new, params_new, logL_new, model)

    inc <- c(
      row = obj_row - old_obj,
      column = obj_col - obj_row,
      proportions = obj_prop - obj_col,
      gaussian = obj_new - obj_prop,
      total = obj_new - old_obj
    )
    substeps[[iter]] <- inc
    if (any(inc < -monotone_tolerance)) {
      return(list(
        ok = FALSE,
        reason = sprintf("%s objective decreased beyond numerical tolerance", spec$code),
        iteration = iter,
        increments = inc,
        variance_floor_diagnostics = floor_diagnostics()
      ))
    }

    membership_delta <- max(
      max(abs(tau_new - old_tau)),
      max(abs(omega_new - old_omega))
    )
    membership_deltas <- c(membership_deltas, membership_delta)

    tau <- tau_new; omega <- omega_new
    params <- params_new; logL <- logL_new
    objective <- obj_new; trace <- c(trace, objective)

    objective_small <- abs(objective - old_obj) <= tol * (1 + abs(old_obj))
    hard_stable <- TRUE
    if (spec$row_mode == "H") {
      hard_stable <- hard_stable && identical(
        max.col(tau, ties.method = "first"),
        max.col(old_tau, ties.method = "first")
      )
    }
    if (spec$col_mode == "H") {
      hard_stable <- hard_stable && identical(
        max.col(omega, ties.method = "first"),
        max.col(old_omega, ties.method = "first")
      )
    }
    if (objective_small && hard_stable) {
      converged <- TRUE
      break
    }
  }

  z_map <- max.col(tau, ties.method = "first")
  w_map <- max.col(omega, ties.method = "first")
  row_map_nonempty <- length(unique(z_map)) == model$g
  col_map_nonempty <- length(unique(w_map)) == model$m
  list(
    ok = TRUE,
    algorithm = spec$code,
    row_mode = spec$row_mode,
    col_mode = spec$col_mode,
    algorithm_label = spec$label,
    tau = tau,
    omega = omega,
    z = z_map,
    w = w_map,
    params = params,
    objective = objective,
    objective_trace = trace,
    substep_increments = if (length(substeps)) do.call(rbind, substeps) else
      matrix(numeric(), 0, 5, dimnames = list(NULL,
        c("row", "column", "proportions", "gaussian", "total"))),
    membership_delta_trace = membership_deltas,
    converged = converged,
    iterations = length(trace) - 1L,
    row_map_nonempty = row_map_nonempty,
    col_map_nonempty = col_map_nonempty,
    map_nonempty = row_map_nonempty && col_map_nonempty,
    effective_row_sizes = colSums(tau),
    effective_col_sizes = colSums(omega),
    mean_row_entropy = {
      terms <- matrix(0, nrow(tau), ncol(tau))
      active <- tau > 0
      terms[active] <- tau[active] * log(tau[active])
      mean(-rowSums(terms))
    },
    mean_col_entropy = {
      terms <- matrix(0, nrow(omega), ncol(omega))
      active <- omega > 0
      terms[active] <- omega[active] * log(omega[active])
      mean(-rowSums(terms))
    },
    variance_floor_diagnostics = floor_diagnostics()
  )
}

# Low-level configuration-specific fitting wrappers.
fit_classification_start <- function(data, model, start, max_iter, tol,
                                     sigma2_min, monotone_tolerance) {
  fit_membership_start(data, model, start, "HH", max_iter, tol, sigma2_min,
                       effective_weight_min = 0,
                       monotone_tolerance = monotone_tolerance)
}

fit_variational_start <- function(data, model, start, max_iter, tol,
                                  sigma2_min, effective_weight_min,
                                  monotone_tolerance) {
  fit_membership_start(data, model, start, "SS", max_iter, tol, sigma2_min,
                       effective_weight_min = effective_weight_min,
                       monotone_tolerance = monotone_tolerance)
}

fit_hs_start <- function(data, model, start, max_iter, tol,
                         sigma2_min, effective_weight_min,
                         monotone_tolerance) {
  fit_membership_start(data, model, start, "HS", max_iter, tol, sigma2_min,
                       effective_weight_min = effective_weight_min,
                       monotone_tolerance = monotone_tolerance)
}

fit_sh_start <- function(data, model, start, max_iter, tol,
                         sigma2_min, effective_weight_min,
                         monotone_tolerance) {
  fit_membership_start(data, model, start, "SH", max_iter, tol, sigma2_min,
                       effective_weight_min = effective_weight_min,
                       monotone_tolerance = monotone_tolerance)
}

#' Fit a censored Gaussian latent block model
#'
#' Jointly estimate row memberships, column memberships, class proportions,
#' and Gaussian block parameters from exact, censored, and missing reports.
#' The five profiled numerical kernels use the compiled Rcpp backend; starts,
#' validity rules, convergence, and model control remain in R.
#'
#' @param data A [new_cglbm_data()] object or a numeric matrix interpreted as
#'   exact data.
#' @param g,m Positive numbers of row and column classes.
#' @param algorithm Membership configuration: `"HH"`, `"HS"`, `"SH"`, or
#'   `"SS"`; the first letter refers to rows and the second to columns.
#' @param S Variance structure: block-specific (`"Skl"`), shared by row class
#'   (`"Sk"`), shared by column class (`"Sl"`), or globally shared (`"S0"`).
#' @param P Class-proportion structure: estimated freely (`"free"`) or fixed
#'   equally (`"equal"`).
#' @param n_starts Number of valid starts requested.
#' @param starts Optional ordered list of starts, each containing integer
#'   vectors `z` and `w`.
#' @param seed Integer random seed.
#' @param max_iter Maximum GEM iterations per start.
#' @param tol Relative objective convergence tolerance.
#' @param sigma2_min Positive variance floor.
#' @param effective_weight_min Minimum admissible row, column, and observed
#'   block mass.
#' @param monotone_tolerance Allowed numerical decrease in a retained GEM
#'   substep.
#' @param max_start_attempts Maximum attempted starts.
#' @param require_map_nonempty Require every declared class in the final MAP
#'   partitions.
#' @param verbose Emit per-attempt progress messages.
#'
#' @return A `cglbm_fit` object containing memberships, MAP labels, fitted
#'   parameters, objective diagnostics, and Fisher-information summaries.
#' @export
#' @examples
#' y <- matrix(c(-2.1, -1.9, 1.9, 2.1,
#'               -2.0, -1.8, 1.8, 2.0), nrow = 4)
#' dat <- exact_cglbm_data(y)
#' fit <- censored_lbm(
#'   dat, g = 2, m = 1, algorithm = "HH",
#'   n_starts = 1,
#'   starts = list(list(z = c(1L, 1L, 2L, 2L), w = rep(1L, 2))),
#'   max_iter = 30, max_start_attempts = 1, seed = 1
#' )
#' predict(fit, type = "labels")
#' summary(fit)
censored_lbm <- function(data, g, m,
                         algorithm = "SS",
                         S = c("Skl", "Sk", "Sl", "S0"),
                         P = c("free", "equal"),
                         n_starts = 20L, starts = NULL, seed = 1L,
                         max_iter = 400L, tol = 1e-7, sigma2_min = 1e-8,
                         effective_weight_min = 1e-8,
                         monotone_tolerance = 1e-7,
                         max_start_attempts = NULL,
                         require_map_nonempty = TRUE,
                         verbose = FALSE) {
  algorithm <- normalize_algorithm(algorithm)
  S <- match.arg(S)
  P <- match.arg(P)
  if (!inherits(data, "cglbm_data")) data <- exact_cglbm_data(data)
  data <- .cglbm_prepare_cpp_data(data)
  model <- validate_model(list(g = g, m = m, S = S, P = P), data$n, data$d)
  n_starts <- as.integer(n_starts); max_iter <- as.integer(max_iter)
  .cglbm_assert(n_starts >= 1L, "n_starts must be at least 1.")
  .cglbm_assert(max_iter >= 1L, "max_iter must be at least 1.")
  .cglbm_assert(tol > 0, "tol must be positive.")
  .cglbm_assert(sigma2_min > 0, "sigma2_min must be positive.")
  .cglbm_assert(effective_weight_min >= 0,
                 "effective_weight_min cannot be negative.")
  if (is.null(max_start_attempts)) {
    max_start_attempts <- max(10L * n_starts, n_starts)
  }
  max_start_attempts <- as.integer(max_start_attempts)
  .cglbm_assert(max_start_attempts >= n_starts,
                 "max_start_attempts must be at least n_starts.")
  supplied_starts <- .cglbm_or(starts, list())

  cached_complete_censoring_reason <- NULL
  if (algorithm == "SS") {
    cached_complete_censoring_reason <- complete_censoring_in_one_direction(
      data,
      matrix(1, nrow = data$n, ncol = data$d)
    )
  }

  successful <- list(); summaries <- list(); attempt_diagnostics <- list(); attempt <- 0L
  elapsed_begin <- proc.time()[["elapsed"]]
  set.seed(
  as.integer(seed),
  kind = "Mersenne-Twister",
  normal.kind = "Inversion")
  while (length(successful) < n_starts && attempt < max_start_attempts) {
    attempt <- attempt + 1L
    if (attempt <= length(supplied_starts)) {
      start <- supplied_starts[[attempt]]
    } else {
      start <- list(
        z = sample_nonempty_labels(data$n, model$g),
        w = sample_nonempty_labels(data$d, model$m),
        start_id = attempt
      )
    }

    one_begin <- proc.time()[["elapsed"]]
    fit <- tryCatch(
      fit_membership_start(
        data = data,
        model = model,
        start = start,
        algorithm = algorithm,
        max_iter = max_iter,
        tol = tol,
        sigma2_min = sigma2_min,
        effective_weight_min = if (algorithm == "HH") 0 else effective_weight_min,
        monotone_tolerance = monotone_tolerance,
        cached_complete_censoring_reason = cached_complete_censoring_reason
      ),
      error = function(e) list(ok = FALSE, reason = conditionMessage(e))
    )
    one_elapsed <- proc.time()[["elapsed"]] - one_begin

    fit_ok <- isTRUE(fit$ok)
    attempt_diagnostics[[attempt]] <- .cglbm_or(fit$variance_floor_diagnostics, NULL)
    eligible_success <- fit_ok &&
      (!isTRUE(require_map_nonempty) || isTRUE(fit$map_nonempty))
    failure_reason <- if (eligible_success) {
      ""
    } else if (fit_ok && isTRUE(require_map_nonempty) && !isTRUE(fit$map_nonempty)) {
      "incomplete MAP partition"
    } else .cglbm_or(fit$reason, "unknown failure")
    summaries[[attempt]] <- data.frame(
      attempt = attempt,
      fit_ok = fit_ok,
      success = eligible_success,
      objective = if (fit_ok) fit$objective else NA_real_,
      converged = if (fit_ok) fit$converged else FALSE,
      iterations = if (fit_ok) fit$iterations else NA_integer_,
      map_nonempty = if (fit_ok) fit$map_nonempty else FALSE,
      row_map_nonempty = if (fit_ok) fit$row_map_nonempty else FALSE,
      col_map_nonempty = if (fit_ok) fit$col_map_nonempty else FALSE,
      min_effective_row_size = if (fit_ok) min(fit$effective_row_sizes) else NA_real_,
      min_effective_col_size = if (fit_ok) min(fit$effective_col_sizes) else NA_real_,
      elapsed_seconds = one_elapsed,
      reason = failure_reason,
      stringsAsFactors = FALSE
    )
    if (eligible_success) successful[[length(successful) + 1L]] <- fit
    if (verbose) {
      message(sprintf(
        "%s attempt %d: %s%s",
        algorithm,
        attempt,
        if (eligible_success) "success" else "failed",
        if (eligible_success) sprintf(", objective=%.6f", fit$objective) else
          sprintf(", %s", failure_reason)
      ))
    }
  }

  if (!length(successful)) {
    start_summary <- if (length(summaries)) do.call(rbind, summaries) else data.frame()
    tab <- if (nrow(start_summary)) table(start_summary$reason) else integer()
    reason_text <- if (length(tab)) {
      paste(sprintf("%s (n=%d)", names(tab), as.integer(tab)), collapse = "; ")
    } else "no diagnostics"
    condition <- structure(
      list(
        message = sprintf("All starts failed for %s CGLBM: %s", algorithm, reason_text),
        call = NULL,
        start_summary = start_summary,
        n_attempted_starts = attempt,
        variance_floor_diagnostics = attempt_diagnostics
      ),
      class = c("cglbm_all_starts_failed", "error", "condition")
    )
    stop(condition)
  }
  if (length(successful) < n_starts) {
    warning(sprintf(
      "Only %d valid starts were obtained out of requested %d after %d attempts.",
      length(successful), n_starts, attempt
    ), call. = FALSE)
  }

  objectives <- vapply(successful, function(x) x$objective, numeric(1))
  best_index <- which.max(objectives)
  best <- successful[[best_index]]
  best$model <- model
  best$data_dimensions <- c(n = data$n, d = data$d)
  best$n_valid_starts <- length(successful)
  best$n_attempted_starts <- attempt
  best$best_valid_start <- best_index
  best$start_summary <- do.call(rbind, summaries)
  best$variance_floor_attempt_diagnostics <- attempt_diagnostics
  best$control <- list(
    n_starts = n_starts,
    seed = seed,
    max_iter = max_iter,
    tol = tol,
    sigma2_min = sigma2_min,
    effective_weight_min = effective_weight_min,
    monotone_tolerance = monotone_tolerance,
    max_start_attempts = max_start_attempts,
    require_map_nonempty = isTRUE(require_map_nonempty)
  )
  best$call <- match.call()
  class(best) <- c("cglbm_fit", paste0("cglbm_", tolower(algorithm)))
  best <- .cglbm_attach_information(best, data)
  best$elapsed_seconds <- proc.time()[["elapsed"]] - elapsed_begin
  best
}

#' @rdname censored_lbm
#' @param x A `cglbm_fit` object.
#' @param ... Additional arguments, currently unused.
#' @export
print.cglbm_fit <- function(x, ...) {
  cat(sprintf("Censored Gaussian LBM fit (%s: %s)\n",
              x$algorithm, x$algorithm_label))
  cat(sprintf("  model: g=%d, m=%d, S=%s, P=%s\n",
              x$model$g, x$model$m, x$model$S, x$model$P))
  cat(sprintf("  objective: %.6f\n", x$objective))
  cat(sprintf("  iterations: %d; converged: %s\n",
              x$iterations, x$converged))
  cat(sprintf("  valid starts: %d/%d attempts\n",
              x$n_valid_starts, x$n_attempted_starts))
  cat(sprintf("  MAP labels contain all classes: %s\n", x$map_nonempty))
  invisible(x)
}

#' @rdname censored_lbm
#' @param object A `cglbm_fit` object.
#' @return `summary.cglbm_fit()` returns a data frame with one row per block;
#'   `predict.cglbm_fit()` returns labels, memberships, or parameters.
#' @export
summary.cglbm_fit <- function(object, ...) {
  .cglbm_assert(inherits(object, "cglbm_fit"), "object must be a cglbm_fit object.")
  .cglbm_assert(length(object$model$S) == 1L &&
                  object$model$S %in% c("Skl", "Sk", "Sl", "S0"),
                "S must be Skl, Sk, Sl, or S0.")
  # Use exact lookup: an absent cache must not match information_unavailable_reason.
  information <- object[["information"]]
  if (is.null(information)) {
    if (!identical(object$information_unavailable_reason, "mi_pooling")) {
      .cglbm_stop(paste0("This fit has no Stage 3 information. Recompute information ",
                         "from the reported data and its fixed reporting design."))
    }
    # A derived MI-pooled object has no valid single-fit information. Do not
    # borrow reference weights or report counts after changing memberships.
    unavailable <- matrix(NA_real_, object$model$g, object$model$m)
    information <- list(
      mean_standard_error = unavailable, effective_block_weight = unavailable,
      exact_count = unavailable, left_censored_count = unavailable,
      right_censored_count = unavailable, effective_sample_size = unavailable,
      nominal_cell_count = tcrossprod(tabulate(object$z, object$model$g),
                                      tabulate(object$w, object$model$m))
    )
    warning(paste0("Standard errors and effective sample sizes are unavailable for MI-pooled fits. ",
                    "Reference-fit information was not retained; MI uncertainty pooling is not implemented."),
            call. = FALSE)
  } else {
    if (!identical(information$method, "expected_fisher_fixed_reporting_rule")) {
      .cglbm_stop(paste0("This fit has an older information definition. Recompute information ",
                         "from the reported data and its fixed reporting design."))
    }
    reasons <- unique(as.vector(information$information_status))
    reasons <- reasons[reasons != "ok"]
    descriptions <- c(
      finite_interval_unsupported = "finite interval-censored reports (information not implemented)",
      reporting_design_unavailable = "missing or invalid fixed reporting design metadata",
      reporting_design_inconsistent = "reporting design limits inconsistent with reports on the analysis scale",
      singular_information = "zero or numerically singular block information",
      numerical_information_unavailable = "numerically unavailable block information",
      shared_variance_group_unavailable = "unavailable information in another block sharing its variance"
    )
    if (length(reasons)) {
      details <- unname(descriptions[reasons])
      details[is.na(details)] <- reasons[is.na(details)]
      warning(paste0("Standard errors and/or effective sample sizes are NA for affected blocks: ",
                      paste(details, collapse = "; "), "."), call. = FALSE)
    }
  }
  rows <- vector("list", object$model$g * object$model$m)
  index <- 0L
  for (k in seq_len(object$model$g)) for (l in seq_len(object$model$m)) {
    index <- index + 1L
    rows[[index]] <- data.frame(
      K = as.integer(k), L = as.integer(l),
      fitted_mean = as.numeric(object$params$mu[k, l]),
      standard_error = as.numeric(information$mean_standard_error[k, l]),
      fitted_variance = as.numeric(object$params$var[k, l]),
      effective_block_weight = as.numeric(information$effective_block_weight[k, l]),
      exact_count = as.integer(information$exact_count[k, l]),
      left_censored_count = as.integer(information$left_censored_count[k, l]),
      right_censored_count = as.integer(information$right_censored_count[k, l]),
      nominal_cell_count = as.integer(information$nominal_cell_count[k, l]),
      effective_sample_size = as.numeric(information$effective_sample_size[k, l]),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' @rdname censored_lbm
#' @param object A `cglbm_fit` object.
#' @param type One of `"labels"`, `"memberships"`, or `"parameters"`.
#' @export
predict.cglbm_fit <- function(object,
                              type = c("labels", "memberships", "parameters"), ...) {
  type <- match.arg(type)
  if (type == "labels") return(list(z = object$z, w = object$w))
  if (type == "parameters") return(object$params)
  list(tau = object$tau, omega = object$omega)
}

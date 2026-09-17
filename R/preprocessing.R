# Censoring-aware column centering, moved from the application implementation.
# Fit one location per column and either one pooled residual variance or one
# variance per column. Only subtract fitted centers; never divide by a scale.
# Censored conditional moments are internal EM quantities, not completed data.
# Source cglbm_core.R first when using these plain R files directly.

.cglbm_center_log2pi <- log(2 * pi)

.cglbm_center_fsum <- function(values) {
  values <- as.numeric(values)
  if (!length(values)) return(0)
  partials <- numeric()
  for (x0 in values) {
    x <- x0
    old <- partials
    new <- numeric()
    if (length(old)) {
      for (y0 in old) {
        y <- y0
        if (abs(x) < abs(y)) {
          tmp <- x; x <- y; y <- tmp
        }
        hi <- x + y
        lo <- y - (hi - x)
        if (lo != 0) new <- c(new, lo)
        x <- hi
      }
    }
    partials <- c(new, x)
  }
  sum(partials)
}

.cglbm_center_mills_tail <- function(x) {
  # For x >= 8, return r(x)-x and the conditional standard-normal variance.
  h <- 0
  dh <- 0
  for (k in 120:1) {
    den <- x + h
    old_dh <- dh
    h <- k / den
    dh <- -k * (1 + old_dh) / (den * den)
  }
  c(excess = h, variance = -dh)
}

.cglbm_center_logsf <- function(x) {
  if (x >= 8) {
    tail <- .cglbm_center_mills_tail(x)
    return(-0.5 * x * x - 0.5 * .cglbm_center_log2pi - log(x + tail[["excess"]]))
  }
  if (x < 0) {
    probability <- stats::pnorm(x, lower.tail = TRUE, log.p = FALSE)
    return(log1p(-probability))
  }
  stats::pnorm(x, lower.tail = FALSE, log.p = TRUE)
}

.cglbm_center_logcdf <- function(x) .cglbm_center_logsf(-x)

.cglbm_center_right_moments <- function(mu, sigma, lower) {
  z <- (lower - mu) / sigma
  if (z >= 8) {
    tail <- .cglbm_center_mills_tail(z)
    return(c(mean = lower + sigma * tail[["excess"]],
             variance = sigma * sigma * tail[["variance"]]))
  }
  ratio <- exp(stats::dnorm(z, log = TRUE) - .cglbm_center_logsf(z))
  variance <- 1 + z * ratio - ratio * ratio
  if (!is.finite(variance) || variance <= 0) {
    stop("invalid conditional normal variance", call. = FALSE)
  }
  c(mean = mu + sigma * ratio, variance = sigma * sigma * variance)
}

.cglbm_center_left_moments <- function(mu, sigma, upper) {
  moments <- .cglbm_center_right_moments(-mu, sigma, -upper)
  c(mean = -moments[["mean"]], variance = moments[["variance"]])
}

.cglbm_center_validate <- function(data) {
  .cglbm_assert(inherits(data, "cglbm_data"), "data must be a cglbm_data object.")
  n <- data$n; d <- data$d
  .cglbm_assert(length(n) == 1L && length(d) == 1L &&
                  is.finite(n) && is.finite(d) && n >= 2L && d >= 1L,
                "Centering requires at least two rows and one column.")
  for (key in c("value", "status", "lower", "upper")) {
    .cglbm_assert(identical(dim(data[[key]]), c(n, d)), "Invalid centering data dimensions.")
  }
  if (any(data$status == "interval", na.rm = TRUE)) {
    .cglbm_stop("Centering of finite interval reports is not implemented.")
  }
  .cglbm_assert(!anyNA(data$status) && all(data$status %in% c("exact", "left", "right", "missing")),
                "Unsupported centering report status.")
  exact <- data$status == "exact"; left <- data$status == "left"
  right <- data$status == "right"; missing <- data$status == "missing"
  .cglbm_assert(all(is.finite(data$value[exact])) &&
                  all(data$lower[exact] == data$value[exact]) &&
                  all(data$upper[exact] == data$value[exact]), "Invalid exact centering report.")
  .cglbm_assert(all(is.na(data$value[!exact])), "Nonexact reports must have missing values.")
  .cglbm_assert(all(data$lower[left] == -Inf) && all(is.finite(data$upper[left])),
                "Invalid left-censored centering report.")
  .cglbm_assert(all(data$upper[right] == Inf) && all(is.finite(data$lower[right])),
                "Invalid right-censored centering report.")
  .cglbm_assert(all(data$lower[missing] == -Inf) && all(data$upper[missing] == Inf),
                "Invalid missing centering report.")
  .cglbm_assert(is.list(data$metadata), "Centering metadata must be a list.")
  invisible(TRUE)
}

.cglbm_center_estep <- function(data, b, sigma) {
  n <- data$n; d <- data$d
  ey <- matrix(0, n, d); ev <- matrix(0, n, d)
  loglik <- numeric(as.numeric(n) * d)
  index <- 0L
  for (i in seq_len(n)) for (j in seq_len(d)) {
    index <- index + 1L
    status <- data$status[i, j]
    if (status == "missing") next
    mu <- b[[j]]
    s <- if (length(sigma) == 1L) sigma[[1L]] else sigma[[j]]
    if (status == "exact") {
      x <- data$value[i, j]
      ey[i, j] <- x
      loglik[[index]] <- -log(s) - 0.5 * .cglbm_center_log2pi - 0.5 * ((x - mu) / s)^2
    } else if (status == "left") {
      upper <- data$upper[i, j]
      moments <- .cglbm_center_left_moments(mu, s, upper)
      ey[i, j] <- moments[["mean"]]; ev[i, j] <- moments[["variance"]]
      loglik[[index]] <- .cglbm_center_logcdf((upper - mu) / s)
    } else if (status == "right") {
      lower <- data$lower[i, j]
      moments <- .cglbm_center_right_moments(mu, s, lower)
      ey[i, j] <- moments[["mean"]]; ev[i, j] <- moments[["variance"]]
      loglik[[index]] <- .cglbm_center_logsf((lower - mu) / s)
    } else .cglbm_stop("Unsupported centering report status.")
  }
  list(ey = ey, ev = ev, loglik = .cglbm_center_fsum(loglik))
}

.cglbm_center_score <- function(ey, ev, b, sigma, observed, variance_mode) {
  n <- nrow(ey); d <- ncol(ey)
  residual <- matrix(0, n, d)
  for (i in seq_len(n)) for (j in seq_len(d)) {
    if (observed[i, j]) residual[i, j] <- ey[i, j] - b[[j]]
  }
  scales <- rep(sigma, length.out = d)
  location_score <- vapply(seq_len(d), function(j) {
    .cglbm_center_fsum(residual[observed[, j], j]) / scales[[j]]^2
  }, numeric(1))
  if (variance_mode == "pooled") {
    scale_score <- (.cglbm_center_fsum(as.vector(ev + residual^2)) -
                      sum(observed) * sigma * sigma) / sigma^3
  } else {
    scale_score <- vapply(seq_len(d), function(j) {
      use <- observed[, j]
      (.cglbm_center_fsum(ev[use, j] + residual[use, j]^2) -
         sum(use) * scales[[j]] * scales[[j]]) / scales[[j]]^3
    }, numeric(1))
  }
  max(abs(c(location_score, scale_score)))
}

.cglbm_fit_centering <- function(data, variance_mode, max_iter,
                                 tolerance, score_tolerance, parameter_tolerance) {
  n <- data$n; d <- data$d
  observed <- data$status != "missing"
  observed_counts <- colSums(observed)
  b <- numeric(d)
  for (j in seq_len(d)) {
    values <- data$value[data$status[, j] == "exact", j]
    if (!length(values)) {
      .cglbm_stop(paste0("Cannot initialize centering column %d without exact reports; ",
                        "no fallback or hidden imputation."), j)
    }
    b[[j]] <- .cglbm_center_fsum(values) / length(values)
  }
  # Preserve the original pooled initialization and summation order.
  if (variance_mode == "pooled") {
    residual_sq <- numeric()
    for (i in seq_len(n)) for (j in seq_len(d)) {
      if (data$status[i, j] == "exact") {
        residual_sq <- c(residual_sq, (data$value[i, j] - b[[j]])^2)
      }
    }
    sigma <- sqrt(.cglbm_center_fsum(residual_sq) / length(residual_sq))
  } else {
    sigma <- vapply(seq_len(d), function(j) {
      values <- data$value[data$status[, j] == "exact", j]
      sqrt(.cglbm_center_fsum((values - b[[j]])^2) / length(values))
    }, numeric(1))
  }
  if (any(!is.finite(sigma)) || any(sigma < 1e-6)) {
    .cglbm_stop("Degenerate auxiliary residual scale.")
  }
  trace <- numeric()
  for (iteration in seq_len(max_iter)) {
    current <- .cglbm_center_estep(data, b, sigma)
    ll <- current$loglik
    if (length(trace) && ll < tail(trace, 1L) - 1e-8) {
      .cglbm_stop("Centering log likelihood decreased.")
    }
    trace <- c(trace, ll)
    bnew <- vapply(seq_len(d), function(j) {
      .cglbm_center_fsum(current$ey[observed[, j], j]) / observed_counts[[j]]
    }, numeric(1))
    centered_sq <- matrix(0, n, d)
    for (i in seq_len(n)) for (j in seq_len(d)) {
      if (observed[i, j]) centered_sq[i, j] <- (current$ey[i, j] - bnew[[j]])^2
    }
    if (variance_mode == "pooled") {
      vnew <- .cglbm_center_fsum(as.vector(current$ev + centered_sq)) / sum(observed)
    } else {
      vnew <- vapply(seq_len(d), function(j) {
        use <- observed[, j]
        .cglbm_center_fsum(current$ev[use, j] + centered_sq[use, j]) / observed_counts[[j]]
      }, numeric(1))
    }
    if (any(!is.finite(vnew)) || any(vnew < 1e-12)) {
      .cglbm_stop("Auxiliary centering variance degeneracy.")
    }
    snew <- sqrt(vnew)
    updated <- .cglbm_center_estep(data, bnew, snew)
    score <- .cglbm_center_score(updated$ey, updated$ev, bnew, snew, observed, variance_mode)
    if (updated$loglik < ll - 1e-8) .cglbm_stop("Centering M step decreased likelihood.")
    change <- max(c(abs(b - bnew), abs(snew - sigma)))
    b <- bnew; sigma <- snew
    if (abs(updated$loglik - ll) <= tolerance * (1 + abs(ll)) &&
        change <= parameter_tolerance * (1 + max(abs(b))) && score <= score_tolerance) {
      trace <- c(trace, updated$loglik)
      return(list(
        method = "center_only_censored_gaussian", centers = b,
        variance_mode = variance_mode, variance = sigma^2, sigma = sigma,
        converged = TRUE, iterations = as.integer(iteration), trace = trace,
        max_score = score, auxiliary_loglik = updated$loglik
      ))
    }
  }
  .cglbm_stop("Centering did not converge within its fixed budget.")
}

.cglbm_shift_centering_metadata <- function(metadata, centers, n, d) {
  # Update both current and legacy complete-rule fields on the analysis scale.
  # If a rule was absent, do not invent it from realized exact/censored statuses.
  for (field in c("reporting_design", "censoring")) {
    rule <- metadata[[field]]
    if (is.null(rule)) next
    .cglbm_assert(is.list(rule), "Centering metadata '%s' must be a list.", field)
    for (key in c("lower_limit", "upper_limit")) {
      if (is.null(rule[[key]])) next
      limits <- expand_to_matrix(rule[[key]], n, d, paste0(field, "$", key))
      rule[[key]] <- sweep(limits, 2L, centers, "-")
    }
    metadata[[field]] <- rule
  }
  metadata
}

#' Center reported matrix columns under censoring
#'
#' Estimate one center per column from exact and one-sided censored reports,
#' using either a residual variance pooled across columns or one variance per
#' column. The function subtracts centers from exact values, report bounds,
#' latent simulation values, and fixed reporting-design limits; it never
#' divides by a scale.
#'
#' @param data A `cglbm_data` object or numeric matrix interpreted as exact
#'   data.
#' @param variance `"pooled"` for one residual variance or `"column"` for one
#'   variance per column.
#' @param max_iter Maximum EM iterations.
#' @param tolerance Relative log-likelihood convergence tolerance.
#' @param score_tolerance Maximum absolute score used in convergence checks.
#' @param parameter_tolerance Maximum relative parameter change used in
#'   convergence checks.
#'
#' @return A centered `cglbm_data` object. Estimation details are stored in
#'   `metadata$centering`.
#' @export
#' @examples
#' x <- matrix(c(1, 2, 3, 4, 10, 11, 12, 13), nrow = 4)
#' centered <- center_cglbm_data(exact_cglbm_data(x), variance = "pooled")
#' centered$metadata$centering$centers
center_cglbm_data <- function(data, variance = c("pooled", "column"),
                              max_iter = 5000L, tolerance = 1e-11,
                              score_tolerance = 1e-6, parameter_tolerance = 1e-7) {
  variance <- match.arg(variance)
  if (!inherits(data, "cglbm_data")) data <- exact_cglbm_data(data)
  .cglbm_center_validate(data)
  .cglbm_assert(is.numeric(max_iter) && length(max_iter) == 1L &&
                  is.finite(max_iter) && max_iter >= 1 && max_iter <= .Machine$integer.max &&
                  max_iter == floor(max_iter), "max_iter must be a positive integer.")
  for (key in c("tolerance", "score_tolerance", "parameter_tolerance")) {
    value <- get(key)
    .cglbm_assert(is.numeric(value) && length(value) == 1L && is.finite(value) && value > 0,
                  "%s must be a positive finite number.", key)
  }
  fit <- .cglbm_fit_centering(data, variance, as.integer(max_iter), tolerance,
                             score_tolerance, parameter_tolerance)
  out <- data
  shifted_value <- sweep(data$value, 2L, fit$centers, "-")
  exact <- data$status == "exact"
  out$value[exact] <- shifted_value[exact]
  out$lower <- sweep(data$lower, 2L, fit$centers, "-")
  out$upper <- sweep(data$upper, 2L, fit$centers, "-")
  if (!is.null(data$latent_truth)) {
    out$latent_truth <- sweep(data$latent_truth, 2L, fit$centers, "-")
  }
  out$metadata <- .cglbm_shift_centering_metadata(data$metadata, fit$centers, data$n, data$d)
  out$metadata$centering <- fit
  .cglbm_center_validate(out)
  out
}

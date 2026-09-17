# Completion procedures used as comparators in the CGLBM manuscript.
#
# This file consolidates the comparison implementations from the formal
# simulation and CARE-SSc analysis code. It intentionally reproduces only the
# four procedures used in the manuscript: Boundary, Half or double, a single
# pooled truncated-normal imputation, and the columnwise Bayesian multiple-
# imputation sampler. Completed matrices are then fitted as exact Gaussian
# latent block models with the existing CGLBM fitting code.

.completion_need <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
  invisible(TRUE)
}

.completion_validate_data <- function(data) {
  .completion_need(inherits(data, "cglbm_data"),
                   "data must be a cglbm_data object")
  invisible(data)
}

.completion_with_seed <- function(seed, expr) {
  seed <- as.integer(seed)
  .completion_need(length(seed) == 1L && !is.na(seed),
                   "seed must be a single integer")

  old_kind <- RNGkind()
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }

  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(
    seed,
    kind = "Mersenne-Twister",
    normal.kind = "Inversion"
  )
  force(expr)
}

comparison_report_matrix <- function(data) {
  .completion_validate_data(data)
  out <- data$value
  left <- data$status == "left"
  right <- data$status == "right"
  interval <- data$status == "interval"
  out[left] <- data$upper[left]
  out[right] <- data$lower[right]
  out[interval] <- 0.5 * (data$lower[interval] + data$upper[interval])
  out
}

comparison_completed_matrix_valid <- function(completed, data,
                                              allow_boundary = FALSE) {
  .completion_validate_data(data)
  completed <- as.matrix(completed)
  if (!identical(dim(completed), c(data$n, data$d)) ||
      any(!is.finite(completed))) {
    return(FALSE)
  }

  exact <- data$status == "exact"
  left <- data$status == "left"
  right <- data$status == "right"
  interval <- data$status == "interval"

  exact_ok <- identical(completed[exact], data$value[exact])
  if (allow_boundary) {
    left_ok <- !any(left) || all(completed[left] <= data$upper[left])
    right_ok <- !any(right) || all(completed[right] >= data$lower[right])
    interval_ok <- !any(interval) || all(
      completed[interval] >= data$lower[interval] &
        completed[interval] <= data$upper[interval]
    )
  } else {
    left_ok <- !any(left) || all(completed[left] < data$upper[left])
    right_ok <- !any(right) || all(completed[right] > data$lower[right])
    interval_ok <- !any(interval) || all(
      completed[interval] > data$lower[interval] &
        completed[interval] < data$upper[interval]
    )
  }

  isTRUE(exact_ok && left_ok && right_ok && interval_ok)
}

comparison_exact_data <- function(completed, data, method,
                                  seed = NA_integer_,
                                  allow_boundary = FALSE) {
  .completion_need(
    comparison_completed_matrix_valid(completed, data, allow_boundary),
    paste(method, "completion violated an exact value or reporting interval")
  )
  exact_cglbm_data(
    completed,
    metadata = c(data$metadata, list(
      comparison_method = method,
      imputation_seed = seed
    ))
  )
}

# Formal Experiment I implementation for Boundary substitution.
comparison_boundary <- function(data) {
  completed <- comparison_report_matrix(data)
  left <- data$status == "left"
  right <- data$status == "right"
  completed[left] <- data$upper[left]
  completed[right] <- data$lower[right]
  comparison_exact_data(
    completed, data, "Boundary", allow_boundary = TRUE
  )
}

# Formal Half or double substitution. On the manuscript simulation scale the
# offset is log(2); on the CARE-SSc log2 scale the offset is 1.
comparison_halfdouble <- function(data, offset) {
  .completion_validate_data(data)
  offset <- as.numeric(offset)
  .completion_need(length(offset) == 1L && is.finite(offset) && offset > 0,
                   "offset must be a positive finite number")
  completed <- comparison_report_matrix(data)
  left <- data$status == "left"
  right <- data$status == "right"
  completed[left] <- data$upper[left] - offset
  completed[right] <- data$lower[right] + offset
  comparison_exact_data(completed, data, "HalfDouble")
}

comparison_rtruncnorm_inverse <- function(n, lower, upper, mean, sd) {
  n <- as.integer(n)
  lower <- rep_len(as.numeric(lower), n)
  upper <- rep_len(as.numeric(upper), n)
  mean <- rep_len(as.numeric(mean), n)
  sd <- rep_len(as.numeric(sd), n)
  .completion_need(
    n >= 0L && all(is.finite(mean)) && all(is.finite(sd)) && all(sd > 0),
    "invalid truncated-normal parameters"
  )
  if (n == 0L) return(numeric())

  # Keep the ordinary inverse-CDF path for the ranges used in the manuscript
  # experiments.  This preserves the historical seeded draws exactly.
  lo <- stats::pnorm(lower, mean = mean, sd = sd)
  hi <- stats::pnorm(upper, mean = mean, sd = sd)
  lo[!is.finite(lower) & lower < 0] <- 0
  hi[!is.finite(upper) & upper > 0] <- 1
  lo <- pmin(pmax(lo, 0), 1)
  hi <- pmin(pmax(hi, 0), 1)

  out <- rep(NA_real_, n)
  ordinary <- is.finite(lo) & is.finite(hi) & hi > lo
  if (any(ordinary)) {
    u <- stats::runif(sum(ordinary), min = lo[ordinary], max = hi[ordinary])
    out[ordinary] <- stats::qnorm(
      u, mean = mean[ordinary], sd = sd[ordinary]
    )
  }

  # In a deep same-side tail, ordinary probabilities can both round to zero or
  # one.  Draw on the log-probability scale instead of replacing the random
  # draw by a conditional mean.  qnorm(..., log.p = TRUE) then performs the
  # inverse without exponentiating the tiny tail probability.
  tail_index <- which(!ordinary)
  logspace_add <- function(log_x, log_y) {
    if (is.infinite(log_x) && log_x < 0) return(log_y)
    if (is.infinite(log_y) && log_y < 0) return(log_x)
    maximum <- max(log_x, log_y)
    maximum + log(exp(log_x - maximum) + exp(log_y - maximum))
  }
  for (index in tail_index) {
    a <- (lower[[index]] - mean[[index]]) / sd[[index]]
    b <- (upper[[index]] - mean[[index]]) / sd[[index]]
    .completion_need(a < b, "invalid truncated-normal interval")
    u <- stats::runif(1L)
    u <- min(max(u, .Machine$double.eps), 1 - .Machine$double.eps)

    if (b <= 0) {
      log_fa <- if (is.infinite(a) && a < 0) -Inf else
        stats::pnorm(a, log.p = TRUE)
      log_fb <- stats::pnorm(b, log.p = TRUE)
      log_mass <- logdiffexp(log_fb, log_fa)
      .completion_need(is.finite(log_mass),
                       "left-tail truncated-normal interval has zero mass")
      log_p <- logspace_add(log_fa, log(u) + log_mass)
      z <- stats::qnorm(log_p, log.p = TRUE)
    } else if (a >= 0) {
      log_qa <- stats::pnorm(a, lower.tail = FALSE, log.p = TRUE)
      log_qb <- if (is.infinite(b) && b > 0) -Inf else
        stats::pnorm(b, lower.tail = FALSE, log.p = TRUE)
      log_mass <- logdiffexp(log_qa, log_qb)
      .completion_need(is.finite(log_mass),
                       "right-tail truncated-normal interval has zero mass")
      log_q <- logspace_add(log_qb, log(u) + log_mass)
      z <- stats::qnorm(log_q, lower.tail = FALSE, log.p = TRUE)
    } else {
      stop(
        "central truncated-normal interval unexpectedly failed probability-scale inversion",
        call. = FALSE
      )
    }

    draw <- mean[[index]] + sd[[index]] * z
    .completion_need(
      is.finite(draw) && draw > lower[[index]] && draw < upper[[index]],
      "tail-stable truncated-normal draw violated its interval"
    )
    out[[index]] <- draw
  }

  out
}

# Formal pooled single truncated-normal imputation from Experiment I.
comparison_single_truncnorm <- function(data, seed) {
  .completion_validate_data(data)
  completed <- comparison_report_matrix(data)
  exact_values <- data$value[data$status == "exact"]
  exact_values <- exact_values[is.finite(exact_values)]
  if (length(exact_values) < 2L) {
    exact_values <- completed[is.finite(completed)]
  }

  mu <- mean(exact_values)
  standard_deviation <- stats::sd(exact_values)
  if (!is.finite(mu)) mu <- 0
  if (!is.finite(standard_deviation) || standard_deviation <= 0) {
    standard_deviation <- 1
  }

  censored <- data$status != "exact"
  if (any(censored)) {
    completed[censored] <- .completion_with_seed(seed, {
      comparison_rtruncnorm_inverse(
        sum(censored),
        data$lower[censored],
        data$upper[censored],
        mu,
        standard_deviation
      )
    })
  }

  out <- comparison_exact_data(completed, data, "TruncNorm", seed)
  out$metadata$imputation_mean <- mu
  out$metadata$imputation_sd <- standard_deviation
  out
}

comparison_mi_prior <- function(data) {
  .completion_validate_data(data)
  exact <- data$value[data$status == "exact"]
  exact <- exact[is.finite(exact)]
  used_pseudo <- length(exact) < 2L
  values <- if (used_pseudo) {
    as.vector(initialization_pseudo_matrix(data))
  } else {
    exact
  }

  mu0 <- mean(values)
  s20 <- stats::var(values)
  if (!is.finite(mu0)) mu0 <- 0
  if (!is.finite(s20) || s20 <= 0) s20 <- 1

  list(
    mu0 = as.numeric(mu0),
    s20 = as.numeric(s20),
    kappa0 = 0.01,
    alpha0 = 2.01,
    beta0 = as.numeric(1.01 * s20),
    used_pseudo = isTRUE(used_pseudo)
  )
}

# Formal columnwise Bayesian censored-Gaussian multiple-imputation sampler.
comparison_draw_mi <- function(data, m, seed, burnin = 50L, thin = 20L) {
  .completion_validate_data(data)
  m <- as.integer(m)
  burnin <- as.integer(burnin)
  thin <- as.integer(thin)
  seed <- as.integer(seed)
  .completion_need(length(m) == 1L && !is.na(m) && m >= 1L,
                   "m must be a positive integer")
  .completion_need(length(burnin) == 1L && !is.na(burnin) && burnin >= 0L,
                   "burnin must be nonnegative")
  .completion_need(length(thin) == 1L && !is.na(thin) && thin >= 1L,
                   "thin must be a positive integer")
  .completion_need(length(seed) == 1L && !is.na(seed),
                   "seed must be a single integer")

  prior <- comparison_mi_prior(data)
  save_iterations <- burnin + seq_len(m) * thin
  total_iterations <- save_iterations[[m]]
  chain_begin <- proc.time()[["elapsed"]]

  sampled <- .completion_with_seed(seed, {
    current <- initialization_pseudo_matrix(data)
    exact <- data$status == "exact"
    current[exact] <- data$value[exact]
    completed <- vector("list", m)
    parameter_draws <- vector("list", m)
    cumulative_seconds <- numeric(m)
    save_index <- 0L
    mu_draw <- rep(NA_real_, data$d)
    sigma2_draw <- rep(NA_real_, data$d)

    for (iteration in seq_len(total_iterations)) {
      for (column in seq_len(data$d)) {
        values <- current[, column]
        n_column <- length(values)
        y_bar <- mean(values)
        kappa_n <- prior$kappa0 + n_column
        alpha_n <- prior$alpha0 + n_column / 2
        beta_n <- prior$beta0 +
          0.5 * sum((values - y_bar)^2) +
          prior$kappa0 * n_column * (y_bar - prior$mu0)^2 /
            (2 * kappa_n)
        .completion_need(
          is.finite(beta_n) && beta_n > 0,
          "MI posterior scale is nonpositive or nonfinite"
        )

        sigma2_draw[[column]] <- 1 / stats::rgamma(
          1L, shape = alpha_n, rate = beta_n
        )
        mu_n <- (prior$kappa0 * prior$mu0 + n_column * y_bar) / kappa_n
        mu_draw[[column]] <- stats::rnorm(
          1L,
          mean = mu_n,
          sd = sqrt(sigma2_draw[[column]] / kappa_n)
        )

        censored <- data$status[, column] != "exact"
        if (any(censored)) {
          current[censored, column] <- comparison_rtruncnorm_inverse(
            sum(censored),
            data$lower[censored, column],
            data$upper[censored, column],
            mu_draw[[column]],
            sqrt(sigma2_draw[[column]])
          )
        }
        current[exact[, column], column] <- data$value[exact[, column], column]
      }

      if (iteration > burnin && (iteration - burnin) %% thin == 0L) {
        save_index <- save_index + 1L
        .completion_need(
          comparison_completed_matrix_valid(current, data),
          "a saved MI matrix violated an exact value or reporting interval"
        )
        completed[[save_index]] <- current
        parameter_draws[[save_index]] <- list(
          iteration = iteration,
          mu = mu_draw,
          sigma2 = sigma2_draw
        )
        cumulative_seconds[[save_index]] <-
          proc.time()[["elapsed"]] - chain_begin
      }
    }

    list(
      completed = completed,
      draws = parameter_draws,
      cumulative_seconds = cumulative_seconds
    )
  })

  list(
    completed = sampled$completed,
    draws = sampled$draws,
    diagnostics = list(
      prior = prior,
      chain_seed = seed,
      burnin = burnin,
      thin = thin,
      save_iterations = save_iterations,
      saved_draw_cumulative_wall_seconds = sampled$cumulative_seconds,
      exact_preserved = all(vapply(sampled$completed, function(x) {
        exact <- data$status == "exact"
        identical(x[exact], data$value[exact])
      }, logical(1))),
      interval_valid = all(vapply(
        sampled$completed,
        comparison_completed_matrix_valid,
        logical(1),
        data = data
      ))
    )
  )
}

comparison_align_fit <- function(fit, reference, K, L) {
  .completion_need(
    length(fit$z) == length(reference$z) &&
      length(fit$w) == length(reference$w),
    "fit and reference partitions have incompatible dimensions"
  )

  choose_mapping <- function(observed, target, classes) {
    contingency <- unclass(table(
      factor(as.integer(observed), levels = seq_len(classes)),
      factor(as.integer(target), levels = seq_len(classes))
    ))
    cost <- max(contingency) - contingency
    as.integer(clue::solve_LSAP(cost))
  }

  row_mapping <- choose_mapping(fit$z, reference$z, K)
  column_mapping <- choose_mapping(fit$w, reference$w, L)
  row_old_order <- match(seq_len(K), row_mapping)
  column_old_order <- match(seq_len(L), column_mapping)

  aligned <- fit
  aligned$z <- row_mapping[as.integer(fit$z)]
  aligned$w <- column_mapping[as.integer(fit$w)]
  if (!is.null(fit$tau)) {
    aligned$tau <- fit$tau[, row_old_order, drop = FALSE]
  }
  if (!is.null(fit$omega)) {
    aligned$omega <- fit$omega[, column_old_order, drop = FALSE]
  }
  if (!is.null(fit$params$pi)) {
    aligned$params$pi <- fit$params$pi[row_old_order]
  }
  if (!is.null(fit$params$rho)) {
    aligned$params$rho <- fit$params$rho[column_old_order]
  }
  if (!is.null(fit$params$mu)) {
    aligned$params$mu <- fit$params$mu[
      row_old_order, column_old_order, drop = FALSE
    ]
  }
  if (!is.null(fit$params$var)) {
    aligned$params$var <- fit$params$var[
      row_old_order, column_old_order, drop = FALSE
    ]
  }

  if (!is.null(fit[["information"]])) {
    aligned$information <- .cglbm_reorder_information(
      fit[["information"]], row_old_order, column_old_order
    )
  }

  list(
    fit = aligned,
    row_mapping = row_mapping,
    column_mapping = column_mapping
  )
}

comparison_adjusted_rand_from_table <- function(tab) {
  .completion_need(
    length(dim(tab)) == 2L && length(tab) > 0L,
    "invalid ARI contingency table"
  )

  choose2 <- function(x) x * (x - 1) / 2
  n <- sum(tab)

  if (n < 2) return(1)

  A <- sum(choose2(rowSums(tab)))
  B <- sum(choose2(colSums(tab)))
  expected <- A * B / choose2(n)
  denominator <- (A + B) / 2 - expected

  if (abs(denominator) < 1e-12) return(1)

  as.numeric(
    (sum(choose2(tab)) - expected) / denominator
  )
}

comparison_adjusted_rand <- function(a, b) {
  .completion_need(
    length(a) == length(b) && length(a) > 0L,
    "invalid ARI labels"
  )

  comparison_adjusted_rand_from_table(
    table(a, b, useNA = "ifany")
  )
}

#' Co-clustering adjusted Rand index
#'
#' Compare two row/column partitions without constructing cellwise character
#' labels.
#'
#' @param z_a,w_a Row and column labels for the first co-clustering.
#' @param z_b,w_b Row and column labels for the second co-clustering.
#'
#' @return A numeric co-clustering adjusted Rand index.
#' @export
#' @examples
#' comparison_cari(c(1, 1, 2), c(1, 2), c(2, 2, 1), c(2, 1))
comparison_cari <- function(z_a, w_a, z_b, w_b) {
  .completion_need(
    length(z_a) == length(z_b) && length(w_a) == length(w_b),
    "CARI partitions have incompatible dimensions"
  )

  row_tab <- unclass(
    table(z_a, z_b, useNA = "ifany")
  )
  column_tab <- unclass(
    table(w_a, w_b, useNA = "ifany")
  )

  storage.mode(row_tab) <- "double"
  storage.mode(column_tab) <- "double"

  cell_tab <- kronecker(row_tab, column_tab)

  comparison_adjusted_rand_from_table(cell_tab)
}

# MI pooling used in the manuscript: choose the fitted partition with the
# largest mean pairwise CARI, align every fit to that medoid, then average the
# aligned Gaussian block parameters.
comparison_pool_mi <- function(fits, K, L) {
  .completion_need(is.list(fits) && length(fits) >= 1L,
                   "MI pooling requires completed-matrix fits")
  K <- as.integer(K)
  L <- as.integer(L)
  .completion_need(K >= 1L && L >= 1L,
                   "K and L must be positive")

  pairwise <- diag(1, length(fits))
  if (length(fits) > 1L) {
    for (left in seq_len(length(fits) - 1L)) {
      for (right in (left + 1L):length(fits)) {
        value <- comparison_cari(
          fits[[left]]$z,
          fits[[left]]$w,
          fits[[right]]$z,
          fits[[right]]$w
        )
        pairwise[left, right] <- value
        pairwise[right, left] <- value
      }
    }
    mean_pairwise <- (rowSums(pairwise) - 1) / (length(fits) - 1L)
  } else {
    mean_pairwise <- 1
  }

  medoid_index <- as.integer(which.max(mean_pairwise))
  reference <- fits[[medoid_index]]
  aligned <- lapply(fits, function(fit) {
    comparison_align_fit(fit, reference, K, L)$fit
  })

  pooled <- reference
  pooled$z <- as.integer(reference$z)
  pooled$w <- as.integer(reference$w)
  pooled$tau <- one_hot(pooled$z, K)
  pooled$omega <- one_hot(pooled$w, L)
  if (identical(pooled$model$P, "equal")) {
    pooled$params$pi <- rep(1 / K, K)
    pooled$params$rho <- rep(1 / L, L)
  } else {
    pooled$params$pi <- tabulate(pooled$z, nbins = K) / length(pooled$z)
    pooled$params$rho <- tabulate(pooled$w, nbins = L) / length(pooled$w)
  }
  pooled$params$mu <- Reduce(
    `+`, lapply(aligned, function(fit) fit$params$mu)
  ) / length(aligned)
  pooled$params$var <- Reduce(
    `+`, lapply(aligned, function(fit) fit$params$var)
  ) / length(aligned)
  pooled$objective <- NA_real_
  # Pooled parameters are not the reference fit. Stage 3 does not implement
  # MI uncertainty pooling; never expose its copied information as pooled SE.
  # Keep a named NULL so $information cannot partially match the reason field.
  pooled["information"] <- list(NULL)
  pooled$information_unavailable_reason <- "mi_pooling"
  pooled$mi_diagnostics <- list(
    medoid_index = medoid_index,
    mean_pairwise_cari = as.numeric(mean_pairwise),
    minimum_pairwise_cari = if (length(fits) > 1L) {
      min(pairwise[lower.tri(pairwise)])
    } else {
      1
    },
    pairwise_cari = pairwise
  )
  pooled
}

# Public manuscript-level entry points. These are deliberately thin wrappers
# around the extracted formal-study implementations above.
#' Complete censored reports before co-clustering
#'
#' Completion comparators used in the accompanying manuscript. Boundary and
#' half/double substitution are deterministic. Truncated-normal and multiple
#' imputation require explicit seeds. Completed matrices are returned as exact
#' `cglbm_data` objects.
#'
#' @param data A `cglbm_data` object.
#' @param analysis_scale Scale used for half/double substitution.
#' @param seed Integer random seed.
#' @param m Number of multiple-imputation matrices.
#' @param burnin Number of initial sampler iterations discarded.
#' @param thin Iterations between retained imputations.
#' @param fits List of fitted CGLBM objects from multiple imputations.
#'
#' @return `complete_boundary()`, `complete_halfdouble()`, and
#'   `complete_truncnorm()` return an exact `cglbm_data` object.
#'   `complete_mi()` returns a list with `completed` numeric matrices,
#'   parameter `draws`, and sampler `diagnostics`. `pool_mi_fits()` returns a
#'   pooled `cglbm_fit` after fitting the completed matrices separately.
#' @export
#' @examples
#' x <- matrix(c(-2, -0.5, 0.5, 2, -1.5, 0, 1, 2.5), nrow = 4)
#' dat <- apply_detection_limits(x, lower_limit = -1, upper_limit = 2)
#' complete_boundary(dat)
#' complete_halfdouble(dat, analysis_scale = "log2")
complete_boundary <- function(data) {
  comparison_boundary(data)
}

#' @rdname complete_boundary
#' @export
complete_halfdouble <- function(data, analysis_scale = c("log", "log2")) {
  analysis_scale <- match.arg(analysis_scale)
  offset <- if (analysis_scale == "log") log(2) else 1
  comparison_halfdouble(data, offset = offset)
}

#' @rdname complete_boundary
#' @export
complete_truncnorm <- function(data, seed) {
  comparison_single_truncnorm(data, seed = seed)
}

#' @rdname complete_boundary
#' @export
complete_mi <- function(data, seed, m = 10L, burnin = 50L, thin = 20L) {
  comparison_draw_mi(
    data = data,
    m = m,
    seed = seed,
    burnin = burnin,
    thin = thin
  )
}

#' @rdname complete_boundary
#' @export
pool_mi_fits <- function(fits) {
  .completion_need(is.list(fits) && length(fits) >= 1L,
                   "MI pooling requires completed-matrix fits")
  first <- fits[[1L]]
  K <- as.integer(.cglbm_or(first$model$g, nrow(first$params$mu)))
  L <- as.integer(.cglbm_or(first$model$m, ncol(first$params$mu)))
  comparison_pool_mi(fits, K = K, L = L)
}


# Simulation-facing aliases retained for the formal Experiment I code.  They
# delegate to the same shared implementations used by the CARE-SSc analysis.
substitute_censored_values <- function(data, method = c("boundary", "halfdouble")) {
  method <- match.arg(method)
  if (identical(method, "boundary")) return(comparison_boundary(data))
  comparison_halfdouble(data, offset = log(2))
}

single_truncnorm_imputation <- function(data, seed = 1L) {
  comparison_single_truncnorm(data, seed = seed)
}

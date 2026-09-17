# Package-native Rcpp dispatch for the five profiled numerical kernels.
# The pure-R reference implementations remain available internally under
# .cglbm_r_* names for regression and parity testing.

.cglbm_status_code <- function(status) {
  status <- as.matrix(status)
  codes <- match(
    as.character(status),
    c("exact", "left", "right", "interval", "missing")
  ) - 1L
  if (anyNA(codes)) {
    .cglbm_stop("Unknown CGLBM report status while preparing the compiled backend.")
  }
  matrix(as.integer(codes), nrow = nrow(status), ncol = ncol(status))
}

.cglbm_prepare_cpp_data <- function(data) {
  if (!inherits(data, "cglbm_data")) data <- exact_cglbm_data(data)
  data$.cglbm_status_code <- .cglbm_status_code(data$status)
  data$.cglbm_status_snapshot <- data$status
  data
}

.cglbm_status_from_data <- function(data) {
  cached <- data$.cglbm_status_code
  if (is.integer(cached) && identical(dim(cached), c(data$n, data$d)) &&
      identical(data$.cglbm_status_snapshot, data$status)) {
    return(cached)
  }
  .cglbm_status_code(data$status)
}

.cglbm_safe_callback <- function(fun) {
  force(fun)
  function(...) {
    tryCatch(
      list(ok = TRUE, value = fun(...)),
      error = function(e) list(ok = FALSE, message = conditionMessage(e))
    )
  }
}

all_cell_loglik <- function(data, params, model) {
  cglbm_cpp_all_cell_loglik_raw(
    value = data$value,
    status = .cglbm_status_from_data(data),
    lower = data$lower,
    upper = data$upper,
    block_mean = params$mu,
    block_variance = params$var,
    log_interval_fallback = .cglbm_safe_callback(
      .log_standard_normal_interval_prob_numerical
    )
  )
}

row_membership_scores <- function(omega, params, logL, model) {
  cglbm_cpp_row_membership_scores_raw(
    omega = omega,
    pi = params$pi,
    log_likelihood = logL
  )
}

col_membership_scores <- function(tau, params, logL, model) {
  cglbm_cpp_col_membership_scores_raw(
    tau = tau,
    rho = params$rho,
    log_likelihood = logL
  )
}

membership_objective <- function(tau, omega, params, logL, model) {
  as.numeric(
    cglbm_cpp_membership_objective_raw(
      tau = tau,
      omega = omega,
      pi = params$pi,
      rho = params$rho,
      log_likelihood = logL
    )
  )
}

weighted_gaussian_em_step <- function(
    data, tau, omega, old_params, model,
    sigma2_min, effective_weight_min = 0,
    iteration = NA_integer_,
    cached_complete_censoring_reason = NULL) {
  group_id <- variance_group_id(model$g, model$m, model$S)
  cached <- if (is.null(cached_complete_censoring_reason)) {
    character()
  } else {
    as.character(cached_complete_censoring_reason)
  }

  result <- cglbm_cpp_weighted_gaussian_em_step_raw(
    value = data$value,
    status = .cglbm_status_from_data(data),
    lower = data$lower,
    upper = data$upper,
    tau = tau,
    omega = omega,
    old_mean = old_params$mu,
    old_variance = old_params$var,
    variance_group = group_id,
    sigma2_min = sigma2_min,
    effective_weight_min = effective_weight_min,
    cached_complete_censoring_reason = cached,
    numerical_moment_fallback = .cglbm_safe_callback(
      .truncated_normal_moments_numerical
    ),
    log_interval_fallback = .cglbm_safe_callback(
      .log_standard_normal_interval_prob_numerical
    )
  )

  if (!isTRUE(result$ok)) {
    return(list(ok = FALSE, reason = as.character(result$reason)))
  }

  floor_history <- lapply(seq_along(result$group_id), function(index) {
    variance_floor_history_row(
      phase = "update",
      iteration = iteration,
      group_id = result$group_id[[index]],
      raw_variance = result$raw_variance[[index]],
      sigma2_min = sigma2_min
    )
  })

  list(
    ok = TRUE,
    mu = result$mu,
    var = result$var,
    block_weight = result$block_weight,
    block_ss = result$block_ss,
    variance_floor_history = do.call(rbind, floor_history)
  )
}

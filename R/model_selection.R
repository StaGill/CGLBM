# Model selection for the Censored Gaussian Latent Block Model (CGLBM).
# The implemented score is the BIC approximation to the integrated
# classification likelihood used in the manuscript. The requested S/P extension
# uses the same Laplace rates with its distinct parameter counts (R/README.md).
# The classification likelihood is evaluated at the retained GEM parameters
# and MAP labels, without a MAP-parameter refit.

.cglbm_candidate_values <- function(x, name, upper) {
  raw <- as.numeric(x)
  if (!length(raw) || any(!is.finite(raw)) || any(raw != floor(raw)) ||
      any(raw < 1) || any(raw > upper)) {
    .cglbm_stop("%s must contain integers between 1 and %d.", name, upper)
  }
  sort(unique(as.integer(raw)))
}

# One fixed structure per call. Vectors are not an implicit S/P search grid.
.cglbm_selection_structures <- function(S, P) {
  if (!is.character(S) || length(S) != 1L || is.na(S) ||
      !S %in% c("Skl", "Sk", "Sl", "S0")) {
    .cglbm_stop("S must be one fixed structure: Skl, Sk, Sl, or S0.")
  }
  if (!is.character(P) || length(P) != 1L || is.na(P) || !P %in% c("free", "equal")) {
    .cglbm_stop("P must be one fixed structure: free or equal.")
  }
  invisible(TRUE)
}

.cglbm_gaussian_parameter_count <- function(K, L, S) {
  K <- as.double(K); L <- as.double(L)
  K * L + switch(S, Skl = K * L, Sk = K, Sl = L, S0 = 1)
}

cglbm_candidate_grid <- function(K, L, n, m, S = "Skl", P = "free") {
  .cglbm_selection_structures(S, P)
  K <- .cglbm_candidate_values(K, "K", as.integer(n))
  L <- .cglbm_candidate_values(L, "L", as.integer(m))
  out <- expand.grid(
    K = K,
    L = L,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  out <- out[order(out$K, out$L), , drop = FALSE]
  out$candidate <- sprintf("K%d_L%d", out$K, out$L)
  K_numeric <- as.double(out$K)
  L_numeric <- as.double(out$L)
  proportion_count <- if (P == "free") (K_numeric - 1) + (L_numeric - 1) else 0
  out$parameter_count <- proportion_count + .cglbm_gaussian_parameter_count(K_numeric, L_numeric, S)
  rownames(out) <- NULL
  out
}

.cglbm_classification_labels <- function(
    labels, size, n_classes, name, require_all_classes = TRUE) {
  raw <- as.numeric(labels)
  if (length(raw) != size || any(!is.finite(raw)) ||
      any(raw != floor(raw)) || any(raw < 1) || any(raw > n_classes)) {
    .cglbm_stop(
      "%s must contain %d integer labels between 1 and %d.",
      name, size, n_classes
    )
  }
  labels <- as.integer(raw)
  if (require_all_classes &&
      !identical(sort(unique(labels)), seq_len(n_classes))) {
    .cglbm_stop(
      "The MAP %s partition must contain all %d declared classes.",
      name, n_classes
    )
  }
  labels
}

cglbm_classification_loglik_components <- function(
    data, z, w, params, model, require_all_classes = TRUE) {
  if (!is.logical(require_all_classes) || length(require_all_classes) != 1L ||
      is.na(require_all_classes)) {
    .cglbm_stop("require_all_classes must be TRUE or FALSE.")
  }
  if (!inherits(data, "cglbm_data")) data <- exact_cglbm_data(data)
  model <- validate_model(model, data$n, data$d)
  z <- .cglbm_classification_labels(
    z, data$n, model$g, "row", require_all_classes
  )
  w <- .cglbm_classification_labels(
    w, data$d, model$m, "column", require_all_classes
  )

  logL <- all_cell_loglik(data, params, model)
  value <- classification_objective(z, w, params, logL, model)
  if (!is.finite(value)) {
    .cglbm_stop("The classification log likelihood is not finite.")
  }
  as.numeric(value)
}

cglbm_classification_loglik <- function(fit, data) {
  .cglbm_assert(inherits(fit, "cglbm_fit"), "fit must be a cglbm_fit object.")
  if (!inherits(data, "cglbm_data")) data <- exact_cglbm_data(data)
  expected <- as.integer(c(n = data$n, d = data$d))
  fitted_dimensions <- as.integer(fit$data_dimensions[c("n", "d")])
  if (length(fitted_dimensions) != 2L || anyNA(fitted_dimensions) ||
      !identical(fitted_dimensions, expected)) {
    .cglbm_stop("data dimensions do not match the fitted model.")
  }
  if (!isTRUE(fit$map_nonempty)) {
    .cglbm_stop("The fitted MAP row and column partitions must contain all declared classes.")
  }

  cglbm_classification_loglik_components(
    data = data,
    z = fit$z,
    w = fit$w,
    params = fit$params,
    model = fit$model
  )
}

cglbm_icl_penalty <- function(K, L, n, m, S = "Skl", P = "free") {
  .cglbm_selection_structures(S, P)
  values <- as.numeric(c(K = K, L = L, n = n, m = m))
  if (length(values) != 4L || any(!is.finite(values)) ||
      any(values != floor(values)) || any(values < 1)) {
    .cglbm_stop("K, L, n, and m must be positive integers.")
  }
  K <- as.double(K)
  L <- as.double(L)
  n <- as.double(n)
  m <- as.double(m)
  # Same Laplace rates as the manuscript: n and m for free proportions,
  # nm for the KL means plus the distinct variance parameters. Fixed equal
  # proportions contribute no estimated parameters, but keep their likelihood.
  proportion_penalty <- if (P == "free") {
    0.5 * (K - 1) * log(n) + 0.5 * (L - 1) * log(m)
  } else 0
  gaussian_coefficient <- if (S == "Skl") K * L else
    0.5 * .cglbm_gaussian_parameter_count(K, L, S)
  proportion_penalty + gaussian_coefficient * log(n * m)
}

#' @rdname select_cglbm_model
#' @param fit A fitted `cglbm_fit` object.
#' @param data The reported data used for fitting.
#' @return `cglbm_icl_bic()` returns one numeric ICL score.
#' @export
cglbm_icl_bic <- function(fit, data) {
  .cglbm_assert(inherits(fit, "cglbm_fit"), "fit must be a cglbm_fit object.")
  if (!inherits(data, "cglbm_data")) data <- exact_cglbm_data(data)
  .cglbm_selection_structures(fit$model$S, fit$model$P)
  model <- validate_model(fit$model, data$n, data$d)
  # Do not score a relabelled model description with incompatible parameters.
  if (model$P == "equal") {
    expected_pi <- rep(1 / model$g, model$g)
    expected_rho <- rep(1 / model$m, model$m)
    if (length(fit$params$pi) != model$g || length(fit$params$rho) != model$m ||
        any(!is.finite(c(fit$params$pi, fit$params$rho))) ||
        any(abs(fit$params$pi - expected_pi) > 64 * .Machine$double.eps * expected_pi) ||
        any(abs(fit$params$rho - expected_rho) > 64 * .Machine$double.eps * expected_rho)) {
      .cglbm_stop("P = equal requires fixed uniform row and column class probabilities.")
    }
  }
  if (model$S != "Skl") {
    .cglbm_assert(identical(dim(fit$params$var), c(model$g, model$m)) &&
                    all(is.finite(fit$params$var)) && all(fit$params$var > 0),
                  "The fitted variance matrix is invalid.")
    gid <- variance_group_id(model$g, model$m, model$S)
    for (gr in unique(as.vector(gid))) {
      if (length(unique(fit$params$var[gid == gr])) != 1L) {
        .cglbm_stop("Fitted block variances violate the declared shared variance structure.")
      }
    }
  }
  loglik <- cglbm_classification_loglik(fit, data)
  penalty <- cglbm_icl_penalty(
    fit$model$g, fit$model$m, data$n, data$d, S = fit$model$S, P = fit$model$P
  )
  as.numeric(loglik - penalty)
}

.cglbm_model_selection_seed <- function(seed, candidate_index, max_start_attempts) {
  seed <- as.numeric(seed)
  candidate_index <- as.integer(candidate_index)
  max_start_attempts <- as.integer(max_start_attempts)
  if (length(seed) != 1L || !is.finite(seed) || seed != floor(seed) ||
      seed < 1 || seed > .Machine$integer.max) {
    .cglbm_stop("seed must be an integer between 1 and .Machine$integer.max.")
  }
  .cglbm_assert(candidate_index >= 1L, "candidate_index must be positive.")
  .cglbm_assert(max_start_attempts >= 1L, "max_start_attempts must be positive.")

  seed_space <- as.double(.Machine$integer.max) - max_start_attempts - 1
  if (seed_space < 1) .cglbm_stop("max_start_attempts is too large for integer seeds.")
  stride <- as.double(max_start_attempts) + 1
  value <- ((seed - 1 + (candidate_index - 1) * stride) %% seed_space) + 1
  as.integer(value)
}

#' Select row and column class numbers by ICL
#'
#' Fit a grid of candidate class-number pairs under one fixed membership,
#' variance, and class-proportion structure, then maximize the implemented
#' BIC approximation to the integrated classification likelihood.
#'
#' @param data A `cglbm_data` object or numeric matrix.
#' @param K,L Candidate row and column class numbers.
#' @param algorithm,S,P Fixed model structures passed to [censored_lbm()].
#' @param n_starts Number of valid starts requested for each candidate.
#' @param seed Base integer seed used to construct deterministic candidate seeds.
#' @param max_iter Maximum GEM iterations per start.
#' @param tol Relative objective convergence tolerance.
#' @param sigma2_min Positive variance floor.
#' @param effective_weight_min Minimum admissible row, column, and observed
#'   block mass.
#' @param monotone_tolerance Allowed numerical decrease in a retained GEM
#'   substep.
#' @param max_start_attempts Maximum attempted starts for each candidate.
#' @param verbose Emit per-candidate progress messages.
#'
#' @return A `cglbm_model_selection` object containing the candidate table,
#'   retained fits, and selected fit.
#' @export
#' @examples
#' y <- matrix(seq_len(24), nrow = 6)
#' sel <- select_cglbm_model(
#'   exact_cglbm_data(y), K = 1, L = 1,
#'   n_starts = 1, max_start_attempts = 1, max_iter = 10, seed = 1
#' )
#' sel
select_cglbm_model <- function(data,
                               K = 1:5,
                               L = 1:5,
                               algorithm = "SS",
                               n_starts = 20L,
                               seed = 1L,
                               max_iter = 400L,
                               tol = 1e-7,
                               sigma2_min = 1e-8,
                               effective_weight_min = 1e-8,
                               monotone_tolerance = 1e-7,
                               max_start_attempts = 200L,
                               verbose = FALSE,
                               S = "Skl",
                               P = "free") {
  if (!inherits(data, "cglbm_data")) data <- exact_cglbm_data(data)
  data <- .cglbm_prepare_cpp_data(data)
  algorithm <- normalize_algorithm(algorithm)
  .cglbm_selection_structures(S, P)
  grid <- cglbm_candidate_grid(K, L, data$n, data$d, S = S, P = P)

  n_starts <- as.integer(n_starts)
  max_iter <- as.integer(max_iter)
  max_start_attempts <- as.integer(max_start_attempts)
  .cglbm_assert(n_starts >= 1L, "n_starts must be at least 1.")
  .cglbm_assert(max_iter >= 1L, "max_iter must be at least 1.")
  .cglbm_assert(max_start_attempts >= n_starts,
                "max_start_attempts must be at least n_starts.")
  .cglbm_assert(is.finite(tol) && tol > 0, "tol must be positive and finite.")
  .cglbm_assert(is.finite(sigma2_min) && sigma2_min > 0,
                "sigma2_min must be positive and finite.")
  .cglbm_assert(is.finite(effective_weight_min) && effective_weight_min >= 0,
                "effective_weight_min must be nonnegative and finite.")
  .cglbm_assert(is.finite(monotone_tolerance) && monotone_tolerance >= 0,
                "monotone_tolerance must be nonnegative and finite.")

  fits <- vector("list", nrow(grid))
  names(fits) <- grid$candidate
  rows <- vector("list", nrow(grid))

  for (index in seq_len(nrow(grid))) {
    K_now <- grid$K[index]
    L_now <- grid$L[index]
    candidate_seed <- .cglbm_model_selection_seed(
      seed, index, max_start_attempts
    )
    if (verbose) {
      message(sprintf("Fitting candidate K=%d, L=%d", K_now, L_now))
    }

    fit <- tryCatch(
      censored_lbm(
        data = data,
        g = K_now,
        m = L_now,
        algorithm = algorithm,
        S = S,
        P = P,
        n_starts = n_starts,
        seed = candidate_seed,
        max_iter = max_iter,
        tol = tol,
        sigma2_min = sigma2_min,
        effective_weight_min = effective_weight_min,
        monotone_tolerance = monotone_tolerance,
        max_start_attempts = max_start_attempts,
        require_map_nonempty = TRUE,
        verbose = FALSE
      ),
      cglbm_all_starts_failed = function(e) e
    )

    eligible <- inherits(fit, "cglbm_fit")
    if (eligible) {
      classification_loglik <- cglbm_classification_loglik(fit, data)
      penalty <- cglbm_icl_penalty(K_now, L_now, data$n, data$d, S = S, P = P)
      icl <- classification_loglik - penalty
      eligible <- is.finite(classification_loglik) && is.finite(icl)
    } else {
      classification_loglik <- NA_real_
      penalty <- cglbm_icl_penalty(K_now, L_now, data$n, data$d, S = S, P = P)
      icl <- NA_real_
    }

    if (eligible) {
      fits[[index]] <- fit
      reason <- ""
    } else {
      reason <- if (inherits(fit, "error")) {
        conditionMessage(fit)
      } else {
        "non-finite classification likelihood or ICL score"
      }
    }

    rows[[index]] <- data.frame(
      K = K_now,
      L = L_now,
      candidate = grid$candidate[index],
      parameter_count = grid$parameter_count[index],
      seed = candidate_seed,
      eligible = eligible,
      n_valid_starts = if (inherits(fit, "cglbm_fit")) fit$n_valid_starts else 0L,
      n_attempted_starts = if (inherits(fit, "cglbm_fit")) fit$n_attempted_starts else
        if (inherits(fit, "cglbm_all_starts_failed")) .cglbm_or(fit$n_attempted_starts, NA_integer_) else
          NA_integer_,
      objective = if (inherits(fit, "cglbm_fit")) fit$objective else NA_real_,
      classification_loglik = classification_loglik,
      penalty = penalty,
      icl = if (eligible) icl else NA_real_,
      elapsed_seconds = if (inherits(fit, "cglbm_fit")) fit$elapsed_seconds else NA_real_,
      reason = reason,
      stringsAsFactors = FALSE
    )
  }

  table <- do.call(rbind, rows)
  eligible_indices <- which(table$eligible & is.finite(table$icl))
  success <- length(eligible_indices) > 0L
  selected_index <- if (success) {
    ranked <- eligible_indices[order(
      -table$icl[eligible_indices],
      table$parameter_count[eligible_indices],
      table$K[eligible_indices] + table$L[eligible_indices],
      table$K[eligible_indices],
      table$L[eligible_indices]
    )]
    ranked[[1L]]
  } else {
    NA_integer_
  }

  out <- list(
    success = success,
    selected_K = if (success) table$K[selected_index] else NA_integer_,
    selected_L = if (success) table$L[selected_index] else NA_integer_,
    selected_icl = if (success) table$icl[selected_index] else NA_real_,
    selected_fit = if (success) fits[[selected_index]] else NULL,
    table = table,
    fits = fits,
    control = list(
      algorithm = algorithm,
      n_starts = n_starts,
      seed = as.integer(seed),
      max_iter = max_iter,
      tol = tol,
      sigma2_min = sigma2_min,
      effective_weight_min = effective_weight_min,
      monotone_tolerance = monotone_tolerance,
      max_start_attempts = max_start_attempts,
      variance_structure = S,
      proportion_structure = P,
      require_map_nonempty = TRUE
    ),
    call = match.call()
  )
  class(out) <- "cglbm_model_selection"
  if (!success) {
    warning("No candidate produced a finite eligible ICL score.", call. = FALSE)
  }
  out
}

#' @rdname select_cglbm_model
#' @param x A `cglbm_model_selection` object.
#' @param ... Additional arguments, currently unused.
#' @export
print.cglbm_model_selection <- function(x, ...) {
  cat("CGLBM model selection by ICL\n")
  if (!isTRUE(x$success)) {
    cat("  no eligible candidate\n")
    return(invisible(x))
  }
  cat(sprintf("  selected model: K=%d, L=%d\n", x$selected_K, x$selected_L))
  cat(sprintf("  ICL: %.6f\n", x$selected_icl))
  cat(sprintf("  eligible candidates: %d/%d\n",
              sum(x$table$eligible), nrow(x$table)))
  invisible(x)
}

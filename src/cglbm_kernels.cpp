#include <Rcpp.h>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

// Match the retained R reference by preventing implicit fused multiply-add.
// The Stage 1 sourceCpp prototype used -ffp-contract=off; source pragmas keep
// the same arithmetic contract without a non-portable package compiler flag.
#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize ("fp-contract=off")
#endif

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

namespace {

inline bool finite_number(const double x) {
  return R_finite(x);
}

inline R_xlen_t array_index(const int i, const int j, const int k, const int l,
                            const int n, const int d, const int g) {
  return static_cast<R_xlen_t>(i) +
    static_cast<R_xlen_t>(n) *
    (static_cast<R_xlen_t>(j) +
     static_cast<R_xlen_t>(d) *
     (static_cast<R_xlen_t>(k) + static_cast<R_xlen_t>(g) * l));
}

inline double logdiffexp_scalar(const double log_x, const double log_y) {
  if (ISNAN(log_x) || ISNAN(log_y)) return NA_REAL;
  if (log_y > log_x) {
    stop("logdiffexp requires log_x >= log_y.");
  }
  if (log_x == log_y) return R_NegInf;
  if (log_x == R_NegInf) return R_NegInf;
  return log_x + std::log(-std::expm1(log_y - log_x));
}

struct LogProbability {
  double value;
  bool needs_fallback;
};

LogProbability standard_normal_interval_logprob(const double a,
                                                const double b) {
  // Validate raw bounds at the caller. Distinct raw endpoints may round to
  // equal standardized endpoints; the reference routes those to quadrature.
  if (ISNAN(a) || ISNAN(b)) return {NA_REAL, false};
  const bool left_inf = !finite_number(a) && a < 0.0;
  const bool right_inf = !finite_number(b) && b > 0.0;

  if (left_inf && right_inf) return {0.0, false};
  if (left_inf) {
    const double value = R::pnorm(b, 0.0, 1.0, 1, 1);
    // The reference R implementation only invokes numerical quadrature for
    // finite intervals. A non-finite one-sided result is rejected later with
    // the same interval-probability error rather than being sent to a fallback
    // that requires two finite endpoints.
    return {value, false};
  }
  if (right_inf) {
    const double value = R::pnorm(a, 0.0, 1.0, 0, 1);
    return {value, false};
  }
  if (!finite_number(a) || !finite_number(b)) {
    return {NA_REAL, false};
  }

  const double width = b - a;
  const double scale = std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
  if (width <= 1e-5 * scale) return {NA_REAL, true};

  double value = NA_REAL;
  if (b <= 0.0) {
    const double log_b = R::pnorm(b, 0.0, 1.0, 1, 1);
    const double log_a = R::pnorm(a, 0.0, 1.0, 1, 1);
    value = logdiffexp_scalar(log_b, log_a);
  } else if (a >= 0.0) {
    const double log_a = R::pnorm(a, 0.0, 1.0, 0, 1);
    const double log_b = R::pnorm(b, 0.0, 1.0, 0, 1);
    value = logdiffexp_scalar(log_a, log_b);
  } else {
    const double probability =
      R::pnorm(b, 0.0, 1.0, 1, 0) -
      R::pnorm(a, 0.0, 1.0, 1, 0);
    value = probability > 0.0 ? std::log(probability) : R_NegInf;
  }

  return {value, !finite_number(value)};
}

// The R reference resolves all narrow finite intervals first, then retries
// the remaining non-finite finite intervals. Keep that vector-wide order
// rather than interleaving quadrature with the ordinary per-cell path.
struct IntervalFallback {
  R_xlen_t index;
  double a;
  double b;
  bool narrow;
};

inline bool is_narrow_standard_interval(const double a, const double b) {
  return finite_number(a) && finite_number(b) &&
    b - a <= 1e-5 * std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
}

void resolve_interval_fallbacks(double* values,
                               const std::vector<IntervalFallback>& pending,
                               Function fallback) {
  for (const auto& cell : pending) {
    if (!cell.narrow) continue;
    List result = fallback(cell.a, cell.b, 1e-12);
    if (!as<bool>(result["ok"])) stop(as<std::string>(result["message"]));
    values[cell.index] = as<double>(result["value"]);
  }
  // This also retries a narrow interval whose first fallback stayed non-finite,
  // exactly as bad_finite does after narrow_finite in log_normal_interval_prob.
  for (const auto& cell : pending) {
    if (finite_number(values[cell.index])) continue;
    List result = fallback(cell.a, cell.b, 1e-12);
    if (!as<bool>(result["ok"])) stop(as<std::string>(result["message"]));
    values[cell.index] = as<double>(result["value"]);
  }
  for (const auto& cell : pending) {
    if (!finite_number(values[cell.index])) {
      stop("Failed to compute one or more positive Gaussian interval probabilities.");
    }
  }
}

struct ConditionalMoment {
  double mean;
  double variance;
};

// R callbacks return a status envelope. This preserves the original
// conditionMessage without allowing Rcpp's callback machinery to reword it.
double interval_logprob_with_fallback(const double a, const double b,
                                      Function fallback) {
  LogProbability probability = standard_normal_interval_logprob(a, b);
  if (probability.needs_fallback) {
    List result = fallback(a, b, 1e-12);
    if (!as<bool>(result["ok"])) stop(as<std::string>(result["message"]));
    probability.value = as<double>(result["value"]);
  }
  if (!finite_number(probability.value)) {
    stop("Failed to compute one or more positive Gaussian interval probabilities.");
  }
  return probability.value;
}

ConditionalMoment numerical_moment_fallback(Function fallback,
                                             const double a,
                                             const double b) {
  List result = fallback(a, b, 1e-10);
  if (!as<bool>(result["ok"])) stop(as<std::string>(result["message"]));
  List stable = as<List>(result["value"]);
  const double standardized_mean = as<double>(stable["mean"]);
  const double standardized_variance = as<double>(stable["variance"]);
  if (!finite_number(standardized_mean) ||
      !finite_number(standardized_variance) ||
      standardized_variance < 0.0) {
    stop("Failed to compute valid truncated-normal moments.");
  }
  return {standardized_mean, standardized_variance};
}

ConditionalMoment conditional_moment(const double lower,
                                     const double upper,
                                     const double mean,
                                     const double variance,
                                     Function log_interval_fallback,
                                     Function numerical_fallback,
                                     bool& invalid_truncated_moments,
                                     const double validated_log_delta = NA_REAL) {
  if (!finite_number(variance) || variance <= 0.0) {
    stop("All variances must be positive and finite.");
  }
  if (lower >= upper) {
    stop("Every censored interval must satisfy lower < upper.");
  }

  const double sd = std::sqrt(variance);
  const double a = (lower - mean) / sd;
  const double b = (upper - mean) / sd;
  // The R reference always validates this probability BEFORE considering
  // a numerical-moment shortcut. Do not bypass that failure boundary.
  const double log_delta = finite_number(validated_log_delta) ? validated_log_delta :
    interval_logprob_with_fallback(a, b, log_interval_fallback);

  double ratio_a = 0.0;
  double ratio_b = 0.0;
  if (finite_number(a)) {
    ratio_a = std::exp(R::dnorm(a, 0.0, 1.0, 1) - log_delta);
  }
  if (finite_number(b)) {
    ratio_b = std::exp(R::dnorm(b, 0.0, 1.0, 1) - log_delta);
  }

  const double lambda = ratio_a - ratio_b;
  const double a_term = finite_number(a) ? a * ratio_a : 0.0;
  const double b_term = finite_number(b) ? b * ratio_b : 0.0;
  const double chi = a_term - b_term;
  const double standardized_second = 1.0 + chi;
  const double standardized_variance = standardized_second - lambda * lambda;
  double first = mean + sd * lambda;
  double conditional_variance = variance * standardized_variance;
  double second = first * first + conditional_variance;

  const bool finite_interval = finite_number(lower) && finite_number(upper);
  const double standardized_width = b - a;
  const double finite_scale =
    std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
  const bool narrow_finite =
    finite_interval && standardized_width <= 1e-3 * finite_scale;
  const bool deep_finite_tail = finite_interval && (a >= 8.0 || b <= -8.0);
  const bool extreme_one_sided =
    ((!finite_number(a) && a < 0.0 && finite_number(b) && b <= -30.0) ||
     (finite_number(a) && a >= 30.0 && !finite_number(b) && b > 0.0));

  bool support_violation = false;
  bool mean_outside = false;
  if (finite_interval) {
    const double width = upper - lower;
    const double width_squared = width * width;
    const double variance_bound = width_squared / 4.0;
    const double variance_slack =
      10.0 * DBL_EPSILON * std::max(1.0, width_squared);
    support_violation =
      conditional_variance > variance_bound + variance_slack;
    const double endpoint_scale =
      std::max(1.0, std::max(std::fabs(lower), std::fabs(upper)));
    const double mean_slack = 10.0 * DBL_EPSILON * endpoint_scale;
    mean_outside = first < lower - mean_slack || first > upper + mean_slack;
  }

  const bool bad =
    !finite_number(first) || !finite_number(second) ||
    !finite_number(conditional_variance) || conditional_variance < 0.0 ||
    narrow_finite || deep_finite_tail || extreme_one_sided ||
    support_violation || mean_outside;
  if (bad) {
    const ConditionalMoment stable = numerical_moment_fallback(numerical_fallback, a, b);
    first = mean + sd * stable.mean;
    conditional_variance = variance * stable.variance;
    second = first * first + conditional_variance;
  }
  // The vectorized reference checks final moments only after ALL scalar
  // fallback calls. Record failure here, but let later callbacks run first.
  invalid_truncated_moments = invalid_truncated_moments ||
    !finite_number(first) || !finite_number(second) || conditional_variance < 0.0;
  return {first, conditional_variance};
}

IntegerVector checked_loglik_dimensions(const NumericVector& log_likelihood) {
  SEXP dim_sexp = log_likelihood.attr("dim");
  if (Rf_isNull(dim_sexp)) stop("logL must be a four-dimensional array.");
  IntegerVector dimensions(dim_sexp);
  if (dimensions.size() != 4) {
    stop("logL must be a four-dimensional array.");
  }
  R_xlen_t expected = 1;
  for (int axis = 0; axis < 4; ++axis) {
    const int size = dimensions[axis];
    if (size < 1 || expected > std::numeric_limits<R_xlen_t>::max() / size) {
      stop("logL dimensions must be positive and non-overflowing.");
    }
    expected *= size;
  }
  if (expected != log_likelihood.size()) {
    stop("logL dimensions do not match its length.");
  }
  return dimensions;
}

}  // namespace

// [[Rcpp::export(rng = false)]]
NumericVector cglbm_cpp_all_cell_loglik_raw(
    const NumericMatrix& value,
    const IntegerMatrix& status,
    const NumericMatrix& lower,
    const NumericMatrix& upper,
    const NumericMatrix& block_mean,
    const NumericMatrix& block_variance,
    Function log_interval_fallback) {
  const int n = value.nrow();
  const int d = value.ncol();
  const int g = block_mean.nrow();
  const int m = block_mean.ncol();

  if (status.nrow() != n || status.ncol() != d ||
      lower.nrow() != n || lower.ncol() != d ||
      upper.nrow() != n || upper.ncol() != d) {
    stop("Reported-data matrices must have identical dimensions.");
  }
  if (block_variance.nrow() != g || block_variance.ncol() != m) {
    stop("Block mean and variance matrices must have identical dimensions.");
  }
  if (n < 1 || d < 1 || g < 1 || m < 1) {
    stop("Likelihood inputs must have positive dimensions.");
  }

  R_xlen_t total = n;
  for (const int size : {d, g, m}) {
    if (total > std::numeric_limits<R_xlen_t>::max() / size) {
      stop("Likelihood array dimensions overflow the supported length.");
    }
    total *= size;
  }
  NumericVector output(total, NA_REAL);

  // Preserve the reference R block traversal order: k is outer and l is
  // inner. The order matters for the first returned block-specific failure
  // reason when more than one block becomes invalid in the same iteration.
  for (int k = 0; k < g; ++k) {
    for (int l = 0; l < m; ++l) {
      checkUserInterrupt();
      const double mu = block_mean(k, l);
      const double variance = block_variance(k, l);
      if (!finite_number(variance) || variance <= 0.0) {
        stop("Variance must be positive and finite.");
      }
      const double sd = std::sqrt(variance);
      bool invalid_exact = false;
      bool invalid_probability = false;
      std::vector<IntervalFallback> interval_fallbacks;
      // The vectorized R function validates all raw intervals first.
      for (int j = 0; j < d; ++j) for (int i = 0; i < n; ++i) {
        if (status(i, j) >= 1 && status(i, j) <= 4 && lower(i, j) >= upper(i, j)) {
          stop("Every censored interval must satisfy lower < upper.");
        }
      }

      for (int j = 0; j < d; ++j) {
        for (int i = 0; i < n; ++i) {
          const int report_status = status(i, j);
          double log_likelihood = NA_REAL;

          if (report_status == 0) {
            const double observed = value(i, j);
            log_likelihood = R::dnorm(observed, mu, sd, 1);
            invalid_exact = invalid_exact || !finite_number(log_likelihood);
          } else if (report_status >= 1 && report_status <= 4) {
            // Status 4 is missing: canonical bounds (-Inf, Inf) give log(1).
            const double a = (lower(i, j) - mu) / sd;
            const double b = (upper(i, j) - mu) / sd;
            if (lower(i, j) >= upper(i, j)) {
              stop("Every censored interval must satisfy lower < upper.");
            }
            LogProbability probability = standard_normal_interval_logprob(a, b);
            if (probability.needs_fallback) {
              interval_fallbacks.push_back({
                array_index(i, j, k, l, n, d, g), a, b,
                is_narrow_standard_interval(a, b)
              });
            } else {
              invalid_probability = invalid_probability || !finite_number(probability.value);
            }
            log_likelihood = probability.value;
          } else {
            stop("Unknown CGLBM report status code.");
          }

          output[array_index(i, j, k, l, n, d, g)] = log_likelihood;
        }
      }
      resolve_interval_fallbacks(output.begin(), interval_fallbacks,
                                 log_interval_fallback);
      if (invalid_probability) {
        stop("Failed to compute one or more positive Gaussian interval probabilities.");
      }
      if (invalid_exact) stop("A non-finite cell log-likelihood was produced.");
    }
  }

  output.attr("dim") = IntegerVector::create(n, d, g, m);
  return output;
}

// [[Rcpp::export(rng = false)]]
NumericMatrix cglbm_cpp_row_membership_scores_raw(
    const NumericMatrix& omega,
    const NumericVector& pi,
    const NumericVector& log_likelihood) {
  const IntegerVector dimensions = checked_loglik_dimensions(log_likelihood);
  const int n = dimensions[0];
  const int d = dimensions[1];
  const int g = dimensions[2];
  const int m = dimensions[3];
  if (omega.nrow() != d || omega.ncol() != m || pi.size() != g) {
    stop("Inputs have incompatible dimensions for row membership scores.");
  }

  NumericMatrix scores(n, g);
  for (int k = 0; k < g; ++k) {
    const double prior = std::log(pi[k]);
    for (int i = 0; i < n; ++i) scores(i, k) = prior;
    std::vector<double> block_dot(n, 0.0);
    for (int l = 0; l < m; ++l) {
      checkUserInterrupt();
      std::fill(block_dot.begin(), block_dot.end(), 0.0);
      for (int j = 0; j < d; ++j) {
        const double membership = omega(j, l);
        for (int i = 0; i < n; ++i) {
          block_dot[i] +=
            log_likelihood[array_index(i, j, k, l, n, d, g)] * membership;
        }
      }
      // Match the R expression: scores + (ll_kl %*% omega), rather than
      // putting the prior into the inner dot-product accumulator.
      for (int i = 0; i < n; ++i) scores(i, k) += block_dot[i];
    }
  }
  return scores;
}

// [[Rcpp::export(rng = false)]]
NumericMatrix cglbm_cpp_col_membership_scores_raw(
    const NumericMatrix& tau,
    const NumericVector& rho,
    const NumericVector& log_likelihood) {
  const IntegerVector dimensions = checked_loglik_dimensions(log_likelihood);
  const int n = dimensions[0];
  const int d = dimensions[1];
  const int g = dimensions[2];
  const int m = dimensions[3];
  if (tau.nrow() != n || tau.ncol() != g || rho.size() != m) {
    stop("Inputs have incompatible dimensions for column membership scores.");
  }

  NumericMatrix scores(d, m);
  for (int l = 0; l < m; ++l) {
    const double prior = std::log(rho[l]);
    for (int j = 0; j < d; ++j) scores(j, l) = prior;
    for (int k = 0; k < g; ++k) {
      checkUserInterrupt();
      for (int j = 0; j < d; ++j) {
        double contribution = 0.0;
        for (int i = 0; i < n; ++i) {
          contribution +=
            log_likelihood[array_index(i, j, k, l, n, d, g)] * tau(i, k);
        }
        scores(j, l) += contribution;
      }
    }
  }
  return scores;
}

// [[Rcpp::export(rng = false)]]
double cglbm_cpp_membership_objective_raw(
    const NumericMatrix& tau,
    const NumericMatrix& omega,
    const NumericVector& pi,
    const NumericVector& rho,
    const NumericVector& log_likelihood) {
  const IntegerVector dimensions = checked_loglik_dimensions(log_likelihood);
  const int n = dimensions[0];
  const int d = dimensions[1];
  const int g = dimensions[2];
  const int m = dimensions[3];
  if (tau.nrow() != n || tau.ncol() != g ||
      omega.nrow() != d || omega.ncol() != m ||
      pi.size() != g || rho.size() != m) {
    stop("Inputs have incompatible dimensions for the membership objective.");
  }

  double row_term = 0.0;
  for (int k = 0; k < g; ++k) {
    long double mass_sum = 0.0L;
    for (int i = 0; i < n; ++i) mass_sum += tau(i, k);
    const double mass = static_cast<double>(mass_sum);
    if (pi[k] < 0.0) return R_NegInf;
    if (pi[k] == 0.0 && mass > 0.0) return R_NegInf;
    if (pi[k] > 0.0) row_term += mass * std::log(pi[k]);
  }

  double column_term = 0.0;
  for (int l = 0; l < m; ++l) {
    long double mass_sum = 0.0L;
    for (int j = 0; j < d; ++j) mass_sum += omega(j, l);
    const double mass = static_cast<double>(mass_sum);
    if (rho[l] < 0.0) return R_NegInf;
    if (rho[l] == 0.0 && mass > 0.0) return R_NegInf;
    if (rho[l] > 0.0) column_term += mass * std::log(rho[l]);
  }

  double data_term = 0.0;
  std::vector<double> row_dot(n, 0.0);
  for (int k = 0; k < g; ++k) {
    for (int l = 0; l < m; ++l) {
      checkUserInterrupt();
      std::fill(row_dot.begin(), row_dot.end(), 0.0);
      for (int j = 0; j < d; ++j) {
        const double membership = omega(j, l);
        for (int i = 0; i < n; ++i) {
          row_dot[i] +=
            log_likelihood[array_index(i, j, k, l, n, d, g)] * membership;
        }
      }
      long double block_sum = 0.0L;
      for (int i = 0; i < n; ++i) {
        const double term = tau(i, k) * row_dot[i];
        block_sum += term;
      }
      data_term += static_cast<double>(block_sum);
    }
  }

  long double row_entropy_sum = 0.0L;
  for (int k = 0; k < g; ++k) {
    for (int i = 0; i < n; ++i) {
      const double membership = tau(i, k);
      if (membership > 0.0) row_entropy_sum -= membership * std::log(membership);
    }
  }

  long double column_entropy_sum = 0.0L;
  for (int l = 0; l < m; ++l) {
    for (int j = 0; j < d; ++j) {
      const double membership = omega(j, l);
      if (membership > 0.0) {
        column_entropy_sum -= membership * std::log(membership);
      }
    }
  }

  const double row_entropy = static_cast<double>(row_entropy_sum);
  const double column_entropy = static_cast<double>(column_entropy_sum);
  return row_term + column_term + data_term + row_entropy + column_entropy;
}

// [[Rcpp::export(rng = false)]]
List cglbm_cpp_weighted_gaussian_em_step_raw(
    const NumericMatrix& value,
    const IntegerMatrix& status,
    const NumericMatrix& lower,
    const NumericMatrix& upper,
    const NumericMatrix& tau,
    const NumericMatrix& omega,
    const NumericMatrix& old_mean,
    const NumericMatrix& old_variance,
    const IntegerMatrix& variance_group,
    const double sigma2_min,
    const double effective_weight_min,
    const CharacterVector& cached_complete_censoring_reason,
    Function numerical_moment_fallback,
    Function log_interval_fallback) {
  const int n = value.nrow();
  const int d = value.ncol();
  const int g = tau.ncol();
  const int m = omega.ncol();

  if (status.nrow() != n || status.ncol() != d ||
      lower.nrow() != n || lower.ncol() != d ||
      upper.nrow() != n || upper.ncol() != d ||
      tau.nrow() != n || omega.nrow() != d ||
      old_mean.nrow() != g || old_mean.ncol() != m ||
      old_variance.nrow() != g || old_variance.ncol() != m ||
      variance_group.nrow() != g || variance_group.ncol() != m) {
    stop("Inputs have incompatible dimensions for the Gaussian update.");
  }
  if (!finite_number(sigma2_min) || sigma2_min <= 0.0 ||
      !finite_number(effective_weight_min) || effective_weight_min < 0.0) {
    stop("Gaussian update thresholds must be finite and valid.");
  }

  const bool has_cached_reason = cached_complete_censoring_reason.size() == 1;
  const std::string cached_reason = has_cached_reason ?
    as<std::string>(cached_complete_censoring_reason[0]) : std::string();

  NumericMatrix block_weight(g, m);
  NumericMatrix new_mean(g, m);
  NumericMatrix block_ss(g, m);
  const R_xlen_t cell_count = static_cast<R_xlen_t>(n) * d;
  std::vector<double> conditional_mean(cell_count, 0.0);
  std::vector<double> conditional_variance(cell_count, 0.0);

  // weighted_gaussian_em_step() in the R reference visits blocks with k
  // outermost and l innermost. Preserve that order so the first reported
  // block-specific failure remains identical.
  for (int k = 0; k < g; ++k) {
    for (int l = 0; l < m; ++l) {
      checkUserInterrupt();
      long double weight_sum = 0.0L;
      double maximum_weight = R_NegInf;
      int anchor_i = -1;
      int anchor_j = -1;
      bool all_left = true;
      bool all_right = true;
      bool has_active = false;

      for (int j = 0; j < d; ++j) {
        for (int i = 0; i < n; ++i) {
          // Missing cells have no block mass and cannot be the anchor or
          // participate in the complete-one-sided-censoring check.
          if (status(i, j) == 4) continue;
          const double weight = tau(i, k) * omega(j, l);
          weight_sum += weight;
          if (finite_number(weight) && weight > 0.0) {
            has_active = true;
            if (weight > maximum_weight) {
              maximum_weight = weight;
              anchor_i = i;
              anchor_j = j;
            }
            const int report_status = status(i, j);
            all_left = all_left && report_status == 1 && finite_number(upper(i, j));
            all_right = all_right && report_status == 2 && finite_number(lower(i, j));
          }
        }
      }

      const double total_weight = static_cast<double>(weight_sum);
      block_weight(k, l) = total_weight;
      if (!finite_number(total_weight) ||
          total_weight <= effective_weight_min) {
        return List::create(
          _["ok"] = false,
          _["reason"] = "effective block weight too small at (" +
            std::to_string(k + 1) + "," + std::to_string(l + 1) + ")"
        );
      }

      if (has_cached_reason && !cached_reason.empty()) {
        return List::create(
          _["ok"] = false,
          _["reason"] = cached_reason + " in block (" +
            std::to_string(k + 1) + "," + std::to_string(l + 1) + ")"
        );
      }
      if (!has_cached_reason && has_active && all_left) {
        return List::create(
          _["ok"] = false,
          _["reason"] = "complete left censoring in one direction in block (" +
            std::to_string(k + 1) + "," + std::to_string(l + 1) + ")"
        );
      }
      if (!has_cached_reason && has_active && all_right) {
        return List::create(
          _["ok"] = false,
          _["reason"] = "complete right censoring in one direction in block (" +
            std::to_string(k + 1) + "," + std::to_string(l + 1) + ")"
        );
      }

      if (anchor_i < 0 || anchor_j < 0) {
        return List::create(
          _["ok"] = false,
          _["reason"] = "invalid conditional moments"
        );
      }
      // Validate all active censored probabilities before ANY conditional
      // moment fallback, matching the reference's vectorized call. Reuse the
      // variance buffer for log probabilities, then overwrite it with moments.
      bool has_censored = false;
      for (int j = 0; j < d; ++j) for (int i = 0; i < n; ++i) {
        const double weight = tau(i, k) * omega(j, l);
        if (!(finite_number(weight) && weight > 0.0) || status(i, j) == 0 ||
            status(i, j) == 4) continue;
        has_censored = true;
        if (status(i, j) < 1 || status(i, j) > 3) stop("Unknown CGLBM report status code.");
      }
      if (has_censored) {
        const double variance = old_variance(k, l);
        if (!finite_number(variance) || variance <= 0.0) {
          stop("All variances must be positive and finite.");
        }
        const double sd = std::sqrt(variance);
        for (int j = 0; j < d; ++j) for (int i = 0; i < n; ++i) {
          const double weight = tau(i, k) * omega(j, l);
          if (!(finite_number(weight) && weight > 0.0) || status(i, j) == 0 ||
              status(i, j) == 4) continue;
          if (lower(i, j) >= upper(i, j)) {
            stop("Every censored interval must satisfy lower < upper.");
          }
        }
        bool invalid_probability = false;
        std::vector<IntervalFallback> interval_fallbacks;
        for (int j = 0; j < d; ++j) for (int i = 0; i < n; ++i) {
          const double weight = tau(i, k) * omega(j, l);
          if (!(finite_number(weight) && weight > 0.0) || status(i, j) == 0 ||
              status(i, j) == 4) continue;
          const R_xlen_t index = static_cast<R_xlen_t>(i) + static_cast<R_xlen_t>(n) * j;
          const double a = (lower(i, j) - old_mean(k, l)) / sd;
          const double b = (upper(i, j) - old_mean(k, l)) / sd;
          LogProbability probability = standard_normal_interval_logprob(a, b);
          if (probability.needs_fallback) {
            interval_fallbacks.push_back({index, a, b, is_narrow_standard_interval(a, b)});
          } else {
            invalid_probability = invalid_probability || !finite_number(probability.value);
          }
          conditional_variance[index] = probability.value;
        }
        resolve_interval_fallbacks(conditional_variance.data(), interval_fallbacks,
                                   log_interval_fallback);
        if (invalid_probability) {
          stop("Failed to compute one or more positive Gaussian interval probabilities.");
        }
      }
      // Compute every active moment once before reading the anchor. This
      // avoids a duplicate quadrature call and preserves cell traversal.
      bool invalid_truncated_moments = false;
      bool invalid_conditional_moments = false;
      for (int j = 0; j < d; ++j) {
        for (int i = 0; i < n; ++i) {
          const R_xlen_t index = static_cast<R_xlen_t>(i) +
            static_cast<R_xlen_t>(n) * j;
          const double weight = tau(i, k) * omega(j, l);
          if (!(finite_number(weight) && weight > 0.0) || status(i, j) == 4) continue;

          if (status(i, j) == 0) {
            conditional_mean[index] = value(i, j);
            conditional_variance[index] = 0.0;
          } else if (status(i, j) >= 1 && status(i, j) <= 3) {
            const ConditionalMoment moment = conditional_moment(
              lower(i, j), upper(i, j), old_mean(k, l), old_variance(k, l),
              log_interval_fallback, numerical_moment_fallback,
              invalid_truncated_moments, conditional_variance[index]
            );
            conditional_mean[index] = moment.mean;
            conditional_variance[index] = moment.variance;
          } else {
            stop("Unknown CGLBM report status code.");
          }

          invalid_conditional_moments = invalid_conditional_moments ||
            !finite_number(conditional_mean[index]) ||
            !finite_number(conditional_variance[index]) ||
            conditional_variance[index] < 0.0;
        }
      }
      if (invalid_truncated_moments) {
        stop("Failed to compute valid truncated-normal moments.");
      }
      if (invalid_conditional_moments) {
        return List::create(_["ok"] = false, _["reason"] = "invalid conditional moments");
      }

      const R_xlen_t anchor_index = static_cast<R_xlen_t>(anchor_i) +
        static_cast<R_xlen_t>(n) * anchor_j;
      const double anchor = conditional_mean[anchor_index];
      long double centered_sum = 0.0L;
      for (int j = 0; j < d; ++j) {
        for (int i = 0; i < n; ++i) {
          const double weight = tau(i, k) * omega(j, l);
          if (!(finite_number(weight) && weight > 0.0) || status(i, j) == 4) continue;
          const R_xlen_t index = static_cast<R_xlen_t>(i) + static_cast<R_xlen_t>(n) * j;
          const double term = weight * (conditional_mean[index] - anchor);
          centered_sum += term;
        }
      }
      new_mean(k, l) = anchor + static_cast<double>(centered_sum) / total_weight;
      long double residual_sum = 0.0L;
      for (int j = 0; j < d; ++j) {
        for (int i = 0; i < n; ++i) {
          const double weight = tau(i, k) * omega(j, l);
          if (!(finite_number(weight) && weight > 0.0) || status(i, j) == 4) continue;
          const R_xlen_t index = static_cast<R_xlen_t>(i) +
            static_cast<R_xlen_t>(n) * j;
          const double centered = conditional_mean[index] - new_mean(k, l);
          const double term = weight *
            (conditional_variance[index] + centered * centered);
          residual_sum += term;
        }
      }
      block_ss(k, l) = static_cast<double>(residual_sum);
    }
  }

  for (int l = 0; l < m; ++l) {
    for (int k = 0; k < g; ++k) {
      if (!finite_number(block_ss(k, l)) || block_ss(k, l) < 0.0) {
        return List::create(
          _["ok"] = false,
          _["reason"] = "invalid expected residual sum of squares"
        );
      }
    }
  }

  int maximum_group = 0;
  for (int l = 0; l < m; ++l) {
    for (int k = 0; k < g; ++k) {
      maximum_group = std::max(maximum_group, variance_group(k, l));
    }
  }
  NumericMatrix new_variance(g, m);
  NumericVector raw_variance(maximum_group, NA_REAL);
  IntegerVector group_id(maximum_group);
  for (int group = 1; group <= maximum_group; ++group) {
    long double denominator_sum = 0.0L;
    long double numerator_sum = 0.0L;
    for (int l = 0; l < m; ++l) {
      for (int k = 0; k < g; ++k) {
        if (variance_group(k, l) == group) {
          denominator_sum += block_weight(k, l);
          numerator_sum += block_ss(k, l);
        }
      }
    }
    const double denominator = static_cast<double>(denominator_sum);
    const double numerator = static_cast<double>(numerator_sum);
    if (!finite_number(denominator) ||
        denominator <= effective_weight_min) {
      return List::create(
        _["ok"] = false,
        _["reason"] = "variance-group weight too small for group " +
          std::to_string(group)
      );
    }
    const double raw = numerator / denominator;
    raw_variance[group - 1] = raw;
    group_id[group - 1] = group;
    const double fitted = std::max(raw, sigma2_min);
    for (int l = 0; l < m; ++l) {
      for (int k = 0; k < g; ++k) {
        if (variance_group(k, l) == group) new_variance(k, l) = fitted;
      }
    }
  }

  return List::create(
    _["ok"] = true,
    _["mu"] = new_mean,
    _["var"] = new_variance,
    _["block_weight"] = block_weight,
    _["block_ss"] = block_ss,
    _["group_id"] = group_id,
    _["raw_variance"] = raw_variance
  );
}

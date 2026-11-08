#include "ChemicalSolver.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace ChemicalSolver {
namespace {

bool FiniteNonNegative(double value) {
  return std::isfinite(value) && value >= 0.0;
}

Status EvaluateRates(const std::vector<double>& state, double time,
                     RateFunction rates, void* user_data,
                     std::vector<double>& production,
                     std::vector<double>& loss_frequency,
                     Diagnostics& diagnostics) {
  ++diagnostics.rate_evaluations;
  std::fill(production.begin(), production.end(), 0.0);
  std::fill(loss_frequency.begin(), loss_frequency.end(), 0.0);
  const Status callback_status = rates(
      state.data(), static_cast<int>(state.size()), time, production.data(),
      loss_frequency.data(), user_data);
  if (callback_status != Status::Ok) return Status::CallbackFailure;

  for (std::size_t i = 0; i < state.size(); ++i) {
    if (!std::isfinite(production[i]) ||
        !std::isfinite(loss_frequency[i]))
      return Status::NonFiniteRate;
    if (production[i] < 0.0 || loss_frequency[i] < 0.0)
      return Status::NegativeRate;
    if (!std::isfinite(loss_frequency[i] * state[i]))
      return Status::NonFiniteRate;
  }
  return Status::Ok;
}

/* Form the backward-Euler residual
 *
 *   F_i(z) = z_i-y_i-h[P_i(z,t+h)-L_i(z,t+h)z_i].
 *
 * A zero residual is the fully implicit update. Rate validation is repeated
 * at every Newton/finite-difference point, so a callback cannot introduce a
 * NaN or negative kinetic coefficient through a perturbed state. */
Status Residual(const std::vector<double>& candidate,
                const std::vector<double>& old_state, double end_time,
                double step, RateFunction rates, void* user_data,
                std::vector<double>& residual, Diagnostics& diagnostics) {
  const std::size_t n = candidate.size();
  std::vector<double> production(n), loss_frequency(n);
  const Status status = EvaluateRates(candidate, end_time, rates, user_data,
                                      production, loss_frequency, diagnostics);
  if (status != Status::Ok) return status;
  for (std::size_t i = 0; i < n; ++i) {
    residual[i] = candidate[i] - old_state[i] - step *
        (production[i] - loss_frequency[i] * candidate[i]);
    if (!std::isfinite(residual[i])) return Status::NonFiniteRate;
  }
  return Status::Ok;
}

double ScaledInfinityNorm(const std::vector<double>& value,
                          const std::vector<double>& first,
                          const std::vector<double>& second,
                          const Options& options) {
  double result = 0.0;
  for (std::size_t i = 0; i < value.size(); ++i) {
    const double scale = options.absolute_tolerance +
        options.relative_tolerance *
            std::max(std::fabs(first[i]), std::fabs(second[i]));
    result = std::max(result, std::fabs(value[i]) / scale);
  }
  return result;
}

/* Dense Gaussian elimination with partial pivoting. Chemistry systems used
 * through this API are cell-local and typically small; a dense solve keeps
 * the dependency-light core deterministic. A pivot below machine epsilon
 * times the largest row coefficient is a typed singular-Jacobian result,
 * never silently replaced by an explicit step. */
bool SolveLinear(std::vector<double>& matrix, std::vector<double>& rhs,
                 int n) {
  for (int column = 0; column < n; ++column) {
    int pivot = column;
    double pivot_magnitude = std::fabs(matrix[column * n + column]);
    double row_scale = pivot_magnitude;
    for (int row = column; row < n; ++row) {
      for (int j = column; j < n; ++j)
        row_scale = std::max(row_scale, std::fabs(matrix[row * n + j]));
      const double magnitude = std::fabs(matrix[row * n + column]);
      if (magnitude > pivot_magnitude) {
        pivot = row;
        pivot_magnitude = magnitude;
      }
    }
    if (!(pivot_magnitude >
          64.0 * std::numeric_limits<double>::epsilon() *
              std::max(1.0, row_scale)))
      return false;
    if (pivot != column) {
      for (int j = column; j < n; ++j)
        std::swap(matrix[column * n + j], matrix[pivot * n + j]);
      std::swap(rhs[column], rhs[pivot]);
    }
    for (int row = column + 1; row < n; ++row) {
      const double factor = matrix[row * n + column] /
          matrix[column * n + column];
      matrix[row * n + column] = 0.0;
      for (int j = column + 1; j < n; ++j)
        matrix[row * n + j] -= factor * matrix[column * n + j];
      rhs[row] -= factor * rhs[column];
    }
  }
  for (int row = n - 1; row >= 0; --row) {
    double sum = rhs[row];
    for (int j = row + 1; j < n; ++j)
      sum -= matrix[row * n + j] * rhs[j];
    rhs[row] = sum / matrix[row * n + row];
    if (!std::isfinite(rhs[row])) return false;
  }
  return true;
}

Status ImplicitEuler(const std::vector<double>& old_state, double start_time,
                     double step, RateFunction rates, void* user_data,
                     const Options& options, Diagnostics& diagnostics,
                     std::vector<double>& result) {
  const int n = static_cast<int>(old_state.size());
  std::vector<double> production(n), loss_frequency(n);
  Status status = EvaluateRates(old_state, start_time, rates, user_data,
                                production, loss_frequency, diagnostics);
  if (status != Status::Ok) return status;

  // Positive Patankar-Euler prediction provides a robust Newton initial state
  // for stiff destruction: (y+hP)/(1+hL) >= 0 for valid rates.
  result.resize(n);
  for (int i = 0; i < n; ++i) {
    result[i] = (old_state[i] + step * production[i]) /
        (1.0 + step * loss_frequency[i]);
    if (!FiniteNonNegative(result[i])) return Status::NonFiniteRate;
  }

  std::vector<double> residual(n), perturbed_residual(n), trial_residual(n);
  status = Residual(result, old_state, start_time + step, step, rates,
                    user_data, residual, diagnostics);
  if (status != Status::Ok) return status;

  for (int iteration = 0; iteration < options.maximum_newton_iterations;
       ++iteration) {
    ++diagnostics.newton_iterations;
    const double norm =
        ScaledInfinityNorm(residual, old_state, result, options);
    /* The nonlinear algebraic solve must be substantially more accurate than
     * the requested time-integration error.  Stopping at O(1) in this scaled
     * norm lets a small residual in every accepted step accumulate as false
     * creation or destruction of a linear invariant.  1e-6 keeps that drift
     * below the temporal truncation error while remaining above the
     * finite-difference/roundoff floor for the supported double-precision
     * density range. */
    if (norm <= 1.0e-6) return Status::Ok;

    std::vector<double> jacobian(n * n, 0.0);
    for (int column = 0; column < n; ++column) {
      std::vector<double> perturbed = result;
      const double delta = std::sqrt(
          std::numeric_limits<double>::epsilon()) *
          std::max(1.0, std::fabs(result[column]));
      perturbed[column] += delta;
      status = Residual(perturbed, old_state, start_time + step, step, rates,
                        user_data, perturbed_residual, diagnostics);
      if (status != Status::Ok) return status;
      for (int row = 0; row < n; ++row)
        jacobian[row * n + column] =
            (perturbed_residual[row] - residual[row]) / delta;
    }

    std::vector<double> correction(n);
    for (int i = 0; i < n; ++i) correction[i] = -residual[i];
    if (!SolveLinear(jacobian, correction, n))
      return Status::SingularJacobian;

    // Backtracking prevents Newton from crossing the non-negative density
    // domain and requires actual residual reduction. No clipping is used:
    // clipping would violate a conservative reaction budget silently.
    bool accepted = false;
    double fraction = 1.0;
    std::vector<double> trial(n);
    for (int line_search = 0; line_search < 24; ++line_search) {
      bool admissible = true;
      for (int i = 0; i < n; ++i) {
        trial[i] = result[i] + fraction * correction[i];
        admissible = admissible && FiniteNonNegative(trial[i]);
      }
      if (admissible) {
        status = Residual(trial, old_state, start_time + step, step, rates,
                          user_data, trial_residual, diagnostics);
        if (status != Status::Ok) return status;
        const double trial_norm =
            ScaledInfinityNorm(trial_residual, old_state, trial, options);
        if (trial_norm < norm) {
          result.swap(trial);
          residual.swap(trial_residual);
          accepted = true;
          break;
        }
      }
      fraction *= 0.5;
    }
    if (!accepted) return Status::NewtonFailure;
  }
  return Status::NewtonFailure;
}

bool ValidOptions(const Options& options) {
  return std::isfinite(options.relative_tolerance) &&
         options.relative_tolerance > 0.0 &&
         std::isfinite(options.absolute_tolerance) &&
         options.absolute_tolerance > 0.0 &&
         std::isfinite(options.initial_step) && options.initial_step >= 0.0 &&
         std::isfinite(options.minimum_step) && options.minimum_step >= 0.0 &&
         std::isfinite(options.maximum_step) && options.maximum_step >= 0.0 &&
         options.maximum_steps > 0 && options.maximum_newton_iterations > 0;
}

}  // namespace

const char* StatusName(Status status) {
  switch (status) {
    case Status::Ok: return "ok";
    case Status::InvalidArgument: return "invalid-argument";
    case Status::InvalidState: return "invalid-state";
    case Status::InvalidNetwork: return "invalid-network";
    case Status::DisconnectedNetwork: return "disconnected-network";
    case Status::CallbackFailure: return "callback-failure";
    case Status::NegativeRate: return "negative-rate";
    case Status::NonFiniteRate: return "non-finite-rate";
    case Status::SingularJacobian: return "singular-jacobian";
    case Status::NewtonFailure: return "newton-failure";
    case Status::StepUnderflow: return "step-underflow";
    case Status::MaximumStepsExceeded: return "maximum-steps-exceeded";
  }
  return "unknown-status";
}

Status ValidateNetwork(int species_count, const ReactionDefinition* reactions,
                       int reaction_count) {
  if (species_count <= 0 || reactions == nullptr || reaction_count <= 0)
    return Status::InvalidArgument;

  /* Each reaction forms a hyperedge: all reactants and products participating
   * in that event are mutually connected.  A compact union-find computes the
   * undirected connectivity without assigning kinetics or assuming reaction
   * reversibility.  Connectivity is a structural input property and has no
   * units; it does not establish elemental/charge conservation, which must be
   * checked against application-owned species compositions. */
  std::vector<int> parent(species_count);
  std::vector<bool> present(species_count, false);
  for (int i = 0; i < species_count; ++i) parent[i] = i;
  const auto root = [&parent](int value) {
    int result = value;
    while (parent[result] != result) result = parent[result];
    while (parent[value] != value) {
      const int next = parent[value];
      parent[value] = result;
      value = next;
    }
    return result;
  };

  for (int reaction_index = 0; reaction_index < reaction_count;
       ++reaction_index) {
    const ReactionDefinition& reaction = reactions[reaction_index];
    if (reaction.reactants == nullptr || reaction.products == nullptr ||
        reaction.reactant_count <= 0 || reaction.product_count <= 0)
      return Status::InvalidNetwork;

    int anchor = -1;
    for (int side = 0; side < 2; ++side) {
      const ReactionTerm* terms = side == 0 ? reaction.reactants
                                            : reaction.products;
      const int count = side == 0 ? reaction.reactant_count
                                  : reaction.product_count;
      for (int term_index = 0; term_index < count; ++term_index) {
        const ReactionTerm& term = terms[term_index];
        if (term.species < 0 || term.species >= species_count ||
            !std::isfinite(term.coefficient) || term.coefficient <= 0.0)
          return Status::InvalidNetwork;
        for (int previous = 0; previous < term_index; ++previous)
          if (terms[previous].species == term.species)
            return Status::InvalidNetwork;

        present[term.species] = true;
        if (anchor < 0) {
          anchor = term.species;
        }
        else {
          const int anchor_root = root(anchor);
          const int term_root = root(term.species);
          if (anchor_root != term_root) parent[term_root] = anchor_root;
        }
      }
    }
  }

  for (int species = 0; species < species_count; ++species)
    if (!present[species]) return Status::DisconnectedNetwork;
  const int connected_root = root(0);
  for (int species = 1; species < species_count; ++species)
    if (root(species) != connected_root)
      return Status::DisconnectedNetwork;
  return Status::Ok;
}

Status Integrate(double* state, int n, double start_time, double duration,
                 RateFunction rates, void* user_data, const Options& options,
                 Diagnostics* diagnostics) {
  Diagnostics local_diagnostics;
  Diagnostics& report = diagnostics == nullptr ? local_diagnostics : *diagnostics;
  report = Diagnostics();
  report.final_time = start_time;

  if (state == nullptr || rates == nullptr || n <= 0 ||
      !std::isfinite(start_time) || !std::isfinite(duration) ||
      duration < 0.0 || !ValidOptions(options))
    return Status::InvalidArgument;
  for (int i = 0; i < n; ++i)
    if (!FiniteNonNegative(state[i])) return Status::InvalidState;
  if (duration == 0.0) return Status::Ok;

  const double end_time = start_time + duration;
  if (!std::isfinite(end_time) || !(end_time > start_time))
    return Status::InvalidArgument;
  const double minimum_step = options.minimum_step > 0.0
      ? options.minimum_step
      : 64.0 * std::numeric_limits<double>::epsilon() *
            std::max(1.0, std::fabs(end_time));
  const double maximum_step = options.maximum_step > 0.0
      ? std::min(options.maximum_step, duration) : duration;
  double step = options.initial_step > 0.0
      ? std::min(options.initial_step, maximum_step) : maximum_step;
  if (!(maximum_step >= minimum_step) || !(step >= minimum_step))
    return Status::InvalidArgument;

  std::vector<double> current(state, state + n);
  double time = start_time;
  int attempts = 0;
  while (time < end_time) {
    if (++attempts > options.maximum_steps)
      return Status::MaximumStepsExceeded;
    step = std::min(step, end_time - time);
    if (!(step >= minimum_step) || time + step == time)
      return Status::StepUnderflow;

    std::vector<double> full, first_half, second_half;
    Status status = ImplicitEuler(current, time, step, rates, user_data,
                                  options, report, full);
    if (status == Status::Ok)
      status = ImplicitEuler(current, time, 0.5 * step, rates, user_data,
                             options, report, first_half);
    if (status == Status::Ok)
      status = ImplicitEuler(first_half, time + 0.5 * step, 0.5 * step,
                             rates, user_data, options, report, second_half);

    if (status != Status::Ok) {
      if (status == Status::NegativeRate ||
          status == Status::NonFiniteRate ||
          status == Status::CallbackFailure)
        return status;
      ++report.rejected_steps;
      step *= 0.5;
      continue;
    }

    std::vector<double> difference(n);
    for (int i = 0; i < n; ++i) difference[i] = second_half[i] - full[i];
    const double error =
        ScaledInfinityNorm(difference, full, second_half, options);
    if (!std::isfinite(error)) return Status::NonFiniteRate;

    if (error <= 1.0) {
      /* For backward Euler (order p=1), the leading local error of the two
       * half steps is one half that of the full step.  Richardson's formula
       * therefore gives
       *
       *   y_R = y_(h/2,h/2) + [y_(h/2,h/2)-y_h]/(2^p-1)
       *       = 2 y_(h/2,h/2) - y_h,
       *
       * which cancels the O(h^2) local term and makes accepted states
       * globally second order.  Extrapolation is never clipped: if roundoff
       * or a nonsmooth rate law makes any component negative, the step is
       * rejected and retried at half size so positivity and reaction budgets
       * are not silently altered. */
      std::vector<double> extrapolated(n);
      bool admissible = true;
      for (int i = 0; i < n; ++i) {
        extrapolated[i] = 2.0 * second_half[i] - full[i];
        admissible = admissible && FiniteNonNegative(extrapolated[i]);
      }
      if (!admissible) {
        ++report.rejected_steps;
        step *= 0.5;
        continue;
      }
      current.swap(extrapolated);
      time += step;
      ++report.accepted_steps;
      report.final_time = time;
      report.minimum_accepted_step = report.accepted_steps == 1
          ? step : std::min(report.minimum_accepted_step, step);
      report.maximum_accepted_step =
          std::max(report.maximum_accepted_step, step);
      const double factor = error > 0.0
          ? std::max(0.25, std::min(2.0, 0.9 / std::sqrt(error))) : 2.0;
      step = std::min(maximum_step, step * factor);
      std::copy(current.begin(), current.end(), state);
    }
    else {
      ++report.rejected_steps;
      step *= std::max(0.1, 0.9 / std::sqrt(error));
    }
  }
  report.final_time = end_time;
  return Status::Ok;
}

}  // namespace ChemicalSolver

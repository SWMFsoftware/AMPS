/*
 * Dependency-light surface-source normalization and accounting primitives.
 *
 * Full exosphere injection owns particles through PIC/AMR.  The calculations
 * below are the deterministic physical layer underneath that ownership:
 * integration of a surface shape, emission-angle sampling, and duplicate
 * physical/sampled ledgers.  Keeping these operations independent of global
 * PIC state permits exact tests and input validation before particle creation.
 */
#ifndef AMPS_MODELS_EXOSPHERE_SURFACE_SOURCE_H
#define AMPS_MODELS_EXOSPHERE_SURFACE_SOURCE_H

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace Exosphere {
namespace SurfaceSource {

/**
 * Convert a nonnegative surface shape to per-element rates.
 *
 * `shape[i]` is a relative source density and `area_m2[i]` is the curved cell
 * area in m^2.  If Q is the requested total physical source in s^-1, the
 * returned element rate is
 *
 *   Q_i = Q shape_i A_i / sum_j(shape_j A_j).
 *
 * Thus a constant shape is uniform per unit area, not uniform per mesh cell.
 * Negative/nonfinite shapes, nonpositive/nonfinite areas, a negative or
 * nonfinite Q, and an all-zero shape are rejected with typed exceptions.
 * Long-double accumulation reduces loss of small cells beside dominant cells;
 * it does not repair an under-resolved physical shape.
 */
inline std::vector<double> NormalizeShape(
    const std::vector<double>& shape,
    const std::vector<double>& area_m2,
    double requested_rate_per_s) {
  if (shape.empty() || shape.size() != area_m2.size())
    throw std::invalid_argument("shape and area arrays must have equal nonzero size");
  if (!(requested_rate_per_s >= 0.0) || !std::isfinite(requested_rate_per_s))
    throw std::invalid_argument("requested source rate must be finite and nonnegative");

  long double integral = 0.0L;
  for (std::size_t i = 0; i < shape.size(); ++i) {
    if (!(shape[i] >= 0.0) || !std::isfinite(shape[i]))
      throw std::invalid_argument("source shape must be finite and nonnegative");
    if (!(area_m2[i] > 0.0) || !std::isfinite(area_m2[i]))
      throw std::invalid_argument("surface area must be finite and positive");
    integral += static_cast<long double>(shape[i]) * area_m2[i];
  }
  if (!(integral > 0.0L) || !std::isfinite(integral))
    throw std::invalid_argument("surface source shape has zero or nonfinite integral");

  std::vector<double> result(shape.size(), 0.0);
  for (std::size_t i = 0; i < shape.size(); ++i) {
    const long double fraction =
        static_cast<long double>(shape[i]) * area_m2[i] / integral;
    result[i] = static_cast<double>(requested_rate_per_s * fraction);
  }
  return result;
}

enum class AngularLaw {
  UniformSolidAngle,
  KnudsenCosine
};

/**
 * Map two independent U[0,1) variates to an outward unit direction.
 *
 * Let mu = cos(theta) relative to the outward normal.  Uniform solid angle on
 * a hemisphere has p(mu)=1 and therefore mu=u.  A Knudsen effusive flux has
 * p(mu)=2mu and therefore mu=sqrt(u).  Azimuth is phi=2*pi*v.  The resulting
 * vector is rotated through an orthonormal tangent frame, which remains valid
 * at both poles (including an exactly anti-z normal).
 *
 * Inputs are dimensionless.  A nonfinite/zero normal or variate outside
 * [0,1) is a typed invalid_argument.  The routine samples an angular law only;
 * it makes no claim about speed, temperature, AMR ownership, or particle
 * statistical weight.
 */
inline std::array<double, 3> DirectionFromUnitSamples(
    const std::array<double, 3>& outward_normal,
    double u, double v, AngularLaw law) {
  if (!(u >= 0.0 && u < 1.0) || !(v >= 0.0 && v < 1.0) ||
      !std::isfinite(u) || !std::isfinite(v))
    throw std::invalid_argument("angular variates must be finite in [0,1)");
  if (law != AngularLaw::UniformSolidAngle &&
      law != AngularLaw::KnudsenCosine)
    throw std::invalid_argument("unknown surface angular law");

  double n2 = 0.0;
  for (int d = 0; d < 3; ++d) {
    if (!std::isfinite(outward_normal[d]))
      throw std::invalid_argument("surface normal must be finite");
    n2 += outward_normal[d] * outward_normal[d];
  }
  if (!(n2 > 0.0) || !std::isfinite(n2))
    throw std::invalid_argument("surface normal must be nonzero");

  const double inverse_n = 1.0 / std::sqrt(n2);
  std::array<double, 3> n = {{outward_normal[0] * inverse_n,
                              outward_normal[1] * inverse_n,
                              outward_normal[2] * inverse_n}};

  // Choose an axis safely separated from n, then build a right-handed local
  // frame.  This avoids the singular 1/sin(theta) rotation used by several
  // older surface generators at a pole.
  const std::array<double, 3> reference =
      std::fabs(n[2]) < 0.9 ? std::array<double, 3>{{0.0, 0.0, 1.0}}
                            : std::array<double, 3>{{1.0, 0.0, 0.0}};
  std::array<double, 3> tangent1 = {{
      reference[1] * n[2] - reference[2] * n[1],
      reference[2] * n[0] - reference[0] * n[2],
      reference[0] * n[1] - reference[1] * n[0]}};
  const double tangent_norm = std::sqrt(
      tangent1[0] * tangent1[0] + tangent1[1] * tangent1[1] +
      tangent1[2] * tangent1[2]);
  for (int d = 0; d < 3; ++d) tangent1[d] /= tangent_norm;
  const std::array<double, 3> tangent2 = {{
      n[1] * tangent1[2] - n[2] * tangent1[1],
      n[2] * tangent1[0] - n[0] * tangent1[2],
      n[0] * tangent1[1] - n[1] * tangent1[0]}};

  const double mu = law == AngularLaw::UniformSolidAngle ? u : std::sqrt(u);
  const double transverse = std::sqrt(std::fmax(0.0, 1.0 - mu * mu));
  const double pi = 3.141592653589793238462643383279502884;
  const double phi = 2.0 * pi * v;
  const double c = transverse * std::cos(phi);
  const double s = transverse * std::sin(phi);
  return {{mu * n[0] + c * tangent1[0] + s * tangent2[0],
           mu * n[1] + c * tangent1[1] + s * tangent2[1],
           mu * n[2] + c * tangent1[2] + s * tangent2[2]}};
}

/**
 * Invert the exponential waiting-time CDF for a stationary Poisson source.
 *
 * `rate_per_s` is the model-particle event intensity in s^-1 and `u` is an
 * independent U[0,1) variate.  The returned interval
 *
 *   dt = -log(1-u) / rate
 *
 * is in seconds and has E[dt]=1/rate.  Using log1p preserves relative
 * accuracy when u is close to zero.  A zero rate means that no event can
 * occur and returns +infinity; negative/nonfinite rates and invalid variates
 * are configuration errors and therefore raise invalid_argument.  This
 * function performs no RNG access, which lets production and tests share the
 * same event law without sharing a random stream.
 */
inline double PoissonWaitingTime(double rate_per_s, double u) {
  if (!(u >= 0.0 && u < 1.0) || !std::isfinite(u))
    throw std::invalid_argument("Poisson variate must be finite in [0,1)");
  if (rate_per_s == 0.0)
    return std::numeric_limits<double>::infinity();
  if (!(rate_per_s > 0.0) || !std::isfinite(rate_per_s))
    throw std::invalid_argument("Poisson event rate must be finite and nonnegative");
  return -std::log1p(-u) / rate_per_s;
}

/**
 * Select one compiled source process from nonnegative physical rates.
 *
 * For rates Q_i [s^-1], the half-open cumulative intervals have lengths
 * Q_i/Q, Q=sum(Q_i), so this inverse-CDF selector returns i with exactly that
 * probability.  `u` is U[0,1).  The routine rejects an inconsistent supplied
 * total rather than silently renormalizing it: the tolerance is 64 machine
 * epsilons times the larger magnitude of the supplied and recomputed totals.
 * This relative scale remains meaningful for physically valid sub-unit rates;
 * an absolute max(1,Q) scale would accept order-one relative disagreement
 * when Q is much smaller than 1 s^-1.  Zero-rate entries are never selected.
 */
inline std::size_t SelectSourceProcess(
    const double* rates_per_s, std::size_t count,
    double total_rate_per_s, double u) {
  if (rates_per_s == NULL || count == 0)
    throw std::invalid_argument("source-rate table must be nonempty");
  if (!(u >= 0.0 && u < 1.0) || !std::isfinite(u))
    throw std::invalid_argument("source selector variate must be finite in [0,1)");
  if (!(total_rate_per_s > 0.0) || !std::isfinite(total_rate_per_s))
    throw std::invalid_argument("total source rate must be finite and positive");

  long double sum = 0.0L;
  for (std::size_t i = 0; i < count; ++i) {
    if (!(rates_per_s[i] >= 0.0) || !std::isfinite(rates_per_s[i]))
      throw std::invalid_argument("source rates must be finite and nonnegative");
    sum += rates_per_s[i];
  }
  const long double rate_scale = std::fmax(
      std::fabs(sum),
      std::fabs(static_cast<long double>(total_rate_per_s)));
  const long double tolerance =
      64.0L * std::numeric_limits<double>::epsilon() * rate_scale;
  if (std::fabs(sum - static_cast<long double>(total_rate_per_s)) > tolerance)
    throw std::invalid_argument("source-process rates do not close to total rate");

  const long double target =
      static_cast<long double>(u) * static_cast<long double>(total_rate_per_s);
  long double cumulative = 0.0L;
  std::size_t last_positive = count;
  for (std::size_t i = 0; i < count; ++i) {
    if (rates_per_s[i] == 0.0) continue;
    last_positive = i;
    cumulative += rates_per_s[i];
    if (target < cumulative) return i;
  }

  // Roundoff can place a target infinitesimally above the accumulated double
  // sum even though u<1.  Returning the last positive bin preserves the
  // inverse-CDF endpoint without ever selecting an unavailable source.
  if (last_positive != count) return last_positive;
  throw std::invalid_argument("source-rate table contains no positive process");
}

inline std::size_t SelectSourceProcess(
    const std::vector<double>& rates_per_s,
    double total_rate_per_s, double u) {
  return SelectSourceProcess(rates_per_s.data(),rates_per_s.size(),
      total_rate_per_s,u);
}

/**
 * Minimal physical/sampling ledger for one surface element and species.
 *
 * Particle weights are numbers of represented physical particles.  An
 * injection with local statistical time step dt contributes W/dt [s^-1] to
 * both the physical source diagnostic and its sampled duplicate, W to the
 * sampled injected count, and (only for an inventory-depleting source) W to
 * desorption.  An impact contributes W to both global and sampled return
 * counts; the caller supplies the subset A in [0,W] actually adsorbed.
 * Inventory change over the recorded events is A_adsorbed_down - W_up.
 *
 * This ledger intentionally excludes stochastic event scheduling and surface
 * chemistry.  It checks accounting after those choices have been made.
 */
class FluxLedger {
 public:
  FluxLedger()
      : physical_source_rate_(0.0), sampled_source_rate_(0.0),
        sampled_injected_(0.0), desorbed_(0.0), physical_returned_(0.0),
        sampled_returned_(0.0), adsorbed_(0.0) {}

  void RecordInjection(double physical_weight, double time_step_s,
                       bool depletes_surface) {
    ValidateWeight(physical_weight);
    if (!(time_step_s > 0.0) || !std::isfinite(time_step_s))
      throw std::invalid_argument("source time step must be finite and positive");
    const double rate = physical_weight / time_step_s;
    physical_source_rate_ += rate;
    sampled_source_rate_ += rate;
    sampled_injected_ += physical_weight;
    if (depletes_surface) desorbed_ += physical_weight;
  }

  void RecordImpact(double returned_weight, double adsorbed_weight) {
    ValidateWeight(returned_weight);
    ValidateWeight(adsorbed_weight);
    if (adsorbed_weight > returned_weight)
      throw std::invalid_argument("adsorbed weight cannot exceed returned weight");
    physical_returned_ += returned_weight;
    sampled_returned_ += returned_weight;
    adsorbed_ += adsorbed_weight;
  }

  double physical_source_rate_per_s() const { return physical_source_rate_; }
  double sampled_source_rate_per_s() const { return sampled_source_rate_; }
  double sampled_injected_particles() const { return sampled_injected_; }
  double desorbed_particles() const { return desorbed_; }
  double physical_returned_particles() const { return physical_returned_; }
  double sampled_returned_particles() const { return sampled_returned_; }
  double adsorbed_particles() const { return adsorbed_; }
  double inventory_change_particles() const { return adsorbed_ - desorbed_; }

  bool DuplicateLedgersClose(double absolute_tolerance) const {
    if (!(absolute_tolerance >= 0.0) || !std::isfinite(absolute_tolerance))
      throw std::invalid_argument("closure tolerance must be finite and nonnegative");
    return std::fabs(physical_source_rate_ - sampled_source_rate_) <=
               absolute_tolerance &&
           std::fabs(physical_returned_ - sampled_returned_) <=
               absolute_tolerance;
  }

 private:
  static void ValidateWeight(double weight) {
    if (!(weight >= 0.0) || !std::isfinite(weight))
      throw std::invalid_argument("particle weight must be finite and nonnegative");
  }

  double physical_source_rate_;
  double sampled_source_rate_;
  double sampled_injected_;
  double desorbed_;
  double physical_returned_;
  double sampled_returned_;
  double adsorbed_;
};

}  // namespace SurfaceSource
}  // namespace Exosphere

#endif

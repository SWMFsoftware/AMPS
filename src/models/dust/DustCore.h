/*
 * Dependency-light physical kernels shared by the AMPS dust model and its
 * focused tests.  This header deliberately has no PIC, MPI, mesh, or
 * application dependency: the formulas below can therefore be checked
 * against independent SI oracles without constructing a configured AMPS
 * executable.
 *
 * These routines are algebraic contracts, not a second dust model.  The
 * PIC-facing code remains responsible for particle ownership, random-number
 * generation, current collection, time integration, and sampling-buffer
 * layout.  Invalid physical domains are reported with typed C++ exceptions;
 * callers in the legacy PIC path supply already validated configuration.
 */

#ifndef AMPS_SRC_MODELS_DUST_DUST_CORE_H
#define AMPS_SRC_MODELS_DUST_DUST_CORE_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace ElectricallyChargedDust {
namespace Core {

const double kPi = 3.141592653589793238462643383279502884;

// Map a nonnegative radius/speed to equal-width logarithmic bins.  Bin i is
// [xmin exp(i Delta), xmin exp((i+1) Delta)); the final bin also owns xmax.
// Values outside the configured physical interval are clamped, matching the
// dust velocity-group convention and preventing an endpoint from indexing one
// element beyond a sampling array.
inline int LogarithmicBin(double value, double minimum, double maximum,
                          int bin_count) {
  if (!std::isfinite(value) || !std::isfinite(minimum) ||
      !std::isfinite(maximum) || value < 0.0 || minimum <= 0.0 ||
      maximum <= minimum || bin_count <= 0) {
    throw std::invalid_argument(
        "dust logarithmic bin requires 0 < minimum < maximum, a nonnegative "
        "finite value, and at least one bin");
  }
  if (value <= minimum) return 0;
  if (value >= maximum) return bin_count - 1;

  const double coordinate =
      std::log(value / minimum) / std::log(maximum / minimum);
  int index = static_cast<int>(coordinate * bin_count);
  // Roundoff at a representable edge must not escape the allocated interval.
  if (index < 0) index = 0;
  if (index >= bin_count) index = bin_count - 1;
  return index;
}

// A compact spherical grain of material density rho [kg m^-3] and radius a
// [m] has mass m = (4 pi / 3) rho a^3 [kg].  Porosity or nonspherical shape
// factors are outside this kernel and must be supplied by a model-specific
// grain-property callback.
inline double SphericalGrainMass(double radius, double density) {
  if (!std::isfinite(radius) || !std::isfinite(density) || radius <= 0.0 ||
      density <= 0.0) {
    throw std::invalid_argument(
        "spherical dust mass requires positive finite radius and density");
  }
  return (4.0 * kPi / 3.0) * density * radius * radius * radius;
}

// Isolated-sphere capacitance in vacuum is C = 4 pi eps0 a, hence the grain
// potential relative to infinity is phi = q/C [V].  radius is in metres,
// charge in coulombs, and vacuum_permittivity in F m^-1.  Plasma shielding is
// intentionally absent: it belongs to the current/charging closure.
inline double IsolatedSpherePotential(double radius, double charge,
                                      double vacuum_permittivity) {
  if (!std::isfinite(radius) || !std::isfinite(charge) ||
      !std::isfinite(vacuum_permittivity) || radius <= 0.0 ||
      vacuum_permittivity <= 0.0) {
    throw std::invalid_argument(
        "dust potential requires finite charge and positive finite radius "
        "and permittivity");
  }
  return charge / (4.0 * kPi * vacuum_permittivity * radius);
}

// Evaluate a = (q/m) [E + v x B] in SI units.  E is V m^-1, B is tesla,
// v is m s^-1, q is coulombs, and m is kg, so the result is m s^-2.  This is
// the nonrelativistic point-grain Lorentz law used by TotalGrainAcceleration;
// q=0 and E=B=0 are exact force-free limits.
inline void LorentzAcceleration(double charge, double mass,
                                const double electric_field[3],
                                const double magnetic_field[3],
                                const double velocity[3],
                                double acceleration[3]) {
  if (!std::isfinite(charge) || !std::isfinite(mass) || mass <= 0.0) {
    throw std::invalid_argument(
        "dust Lorentz acceleration requires finite charge and positive mass");
  }
  for (int i = 0; i != 3; ++i) {
    if (!std::isfinite(electric_field[i]) ||
        !std::isfinite(magnetic_field[i]) || !std::isfinite(velocity[i])) {
      throw std::invalid_argument("dust Lorentz vectors must be finite");
    }
  }

  const double scale = charge / mass;
  acceleration[0] = scale * (electric_field[0] +
      velocity[1] * magnetic_field[2] - velocity[2] * magnetic_field[1]);
  acceleration[1] = scale * (electric_field[1] -
      velocity[0] * magnetic_field[2] + velocity[2] * magnetic_field[0]);
  acceleration[2] = scale * (electric_field[2] +
      velocity[0] * magnetic_field[1] - velocity[1] * magnetic_field[0]);
}

// Normalize a nonzero three-component Gaussian variate.  Independent
// zero-mean, equal-variance Gaussian components are rotationally invariant;
// dividing by their norm therefore produces a uniform direction on S^2.
// Position a*n and velocity u*n implement the existing isotropic spherical
// surface source.  The caller owns the random Gaussian draw.
inline void NormalizeInjectionDirection(double direction[3]) {
  if (!std::isfinite(direction[0]) || !std::isfinite(direction[1]) ||
      !std::isfinite(direction[2])) {
    throw std::invalid_argument("dust injection direction must be finite");
  }
  const double norm = std::sqrt(direction[0] * direction[0] +
      direction[1] * direction[1] + direction[2] * direction[2]);
  if (!(norm > 0.0) || !std::isfinite(norm)) {
    throw std::invalid_argument(
        "dust injection direction requires a finite nonzero norm");
  }
  for (int i = 0; i != 3; ++i) direction[i] /= norm;
}

// Convert a unit outward normal into the absolute state of a spherical source.
// center and position are metres, sphere_radius is metres, radial_speed is
// m s^-1, and velocity is m s^-1.  The nonzero center term is essential for
// internal spheres not located at the simulation origin.  The direction is
// validated as unit length instead of silently renormalized so source geometry
// and its random angular law remain separate, auditable operations.
inline void SphericalSourceState(const double center[3], double sphere_radius,
                                 double radial_speed,
                                 const double unit_direction[3],
                                 double position[3], double velocity[3]) {
  if (!std::isfinite(sphere_radius) || sphere_radius <= 0.0 ||
      !std::isfinite(radial_speed) || radial_speed < 0.0) {
    throw std::invalid_argument(
        "dust spherical source requires positive radius and nonnegative speed");
  }
  double norm_squared = 0.0;
  for (int i = 0; i != 3; ++i) {
    if (!std::isfinite(center[i]) || !std::isfinite(unit_direction[i]))
      throw std::invalid_argument("dust spherical source vectors must be finite");
    norm_squared += unit_direction[i] * unit_direction[i];
  }
  if (std::fabs(norm_squared - 1.0) >
      64.0 * std::numeric_limits<double>::epsilon()) {
    throw std::invalid_argument(
        "dust spherical source direction must have unit length");
  }
  for (int i = 0; i != 3; ++i) {
    position[i] = center[i] + sphere_radius * unit_direction[i];
    velocity[i] = radial_speed * unit_direction[i];
  }
}

// The physical differential size law is f(a) proportional to a^-p da.  A
// logarithmic numerical draw is uniform in d(ln a), so its physical density
// in that coordinate is proportional to a^(1-p).  This normalization is
// integral_a_min^a_max a^(1-p) d(ln a); the p=1 branch is the analytic limit
// and avoids the former 0/0 singularity.
inline double PowerLawLogNormalization(double minimum_radius,
                                       double maximum_radius,
                                       double power_index) {
  if (!std::isfinite(minimum_radius) || !std::isfinite(maximum_radius) ||
      !std::isfinite(power_index) || minimum_radius <= 0.0 ||
      maximum_radius <= minimum_radius) {
    throw std::invalid_argument(
        "dust power law requires 0 < minimum radius < maximum radius");
  }
  const double exponent = 1.0 - power_index;
  if (std::fabs(exponent) < 1.0e-12)
    return std::log(maximum_radius / minimum_radius);
  return (std::pow(maximum_radius, exponent) -
          std::pow(minimum_radius, exponent)) / exponent;
}

inline double PowerLawLogDensity(double radius, double minimum_radius,
                                 double maximum_radius,
                                 double power_index) {
  if (!std::isfinite(radius) || radius < minimum_radius ||
      radius > maximum_radius) {
    throw std::invalid_argument("dust radius is outside the power-law domain");
  }
  return std::pow(radius, 1.0 - power_index) /
      PowerLawLogNormalization(minimum_radius, maximum_radius, power_index);
}

inline double PowerLawLogMaximum(double minimum_radius, double maximum_radius,
                                 double power_index) {
  (void)PowerLawLogNormalization(minimum_radius, maximum_radius, power_index);
  const double exponent = 1.0 - power_index;
  return std::max(std::pow(minimum_radius, exponent),
                  std::pow(maximum_radius, exponent));
}

// Raw weighted moments accumulated by the PIC sampling buffers.  weight is
// the represented physical grain count, velocity is m s^-1, and charge is C.
// The fields consequently carry N, N m s^-1, N m s^-1, C N, C N m s^-1,
// and |C| N m s^-1 for absolute_current.
// Division by number is intentionally left to the output stage so empty bins
// remain distinguishable from a physical zero mean.
struct WeightedMoments {
  double number;
  double velocity[3];
  double speed;
  double charge;
  double current[3];
  double absolute_current;

  WeightedMoments()
      : number(0.0), speed(0.0), charge(0.0), absolute_current(0.0) {
    velocity[0] = velocity[1] = velocity[2] = 0.0;
    current[0] = current[1] = current[2] = 0.0;
  }
};

inline void AccumulateWeightedMoments(double weight, const double velocity[3],
                                      double charge,
                                      WeightedMoments& moments) {
  if (!std::isfinite(weight) || !std::isfinite(charge) || weight < 0.0) {
    throw std::invalid_argument(
        "dust moments require nonnegative finite weight and finite charge");
  }
  double speed_squared = 0.0;
  for (int i = 0; i != 3; ++i) {
    if (!std::isfinite(velocity[i]))
      throw std::invalid_argument("dust moment velocity must be finite");
    speed_squared += velocity[i] * velocity[i];
    moments.velocity[i] += weight * velocity[i];
    moments.current[i] += weight * charge * velocity[i];
  }
  const double speed = std::sqrt(speed_squared);
  moments.number += weight;
  moments.speed += weight * speed;
  moments.charge += weight * charge;
  moments.absolute_current += weight * std::fabs(charge) * speed;
}

}  // namespace Core
}  // namespace ElectricallyChargedDust

#endif  // AMPS_SRC_MODELS_DUST_DUST_CORE_H

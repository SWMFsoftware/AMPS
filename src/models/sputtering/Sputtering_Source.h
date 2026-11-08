/*
 * Dependency-light sputtering emission and physical-accounting layer.
 *
 * Yield tables answer "products per incident projectile" while PIC owns mesh
 * placement and numerical particles.  This header supplies the deterministic
 * physics between those layers: an explicit emission-angle selector, source
 * rate conversion, statistical event planning, and product conservation
 * ledgers.  It deliberately has no PIC, MPI, AMR, or application dependency.
 */
#ifndef AMPS_MODELS_SPUTTERING_SOURCE_H
#define AMPS_MODELS_SPUTTERING_SOURCE_H

#include "../exosphere/Exosphere_SurfaceSource.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>

namespace Sputtering {
namespace Source {

/** Runtime-selectable angular distributions about the outward surface normal. */
enum class EmissionAngleLaw {
  // Knudsen/Lambert emission: p(Omega)=mu/pi, or p(mu)=2 mu on 0<=mu<=1.
  // This is the default because particles crossing a plane from an isotropic
  // velocity distribution carry the projected-flux factor mu=cos(theta).
  KnudsenCosine,

  // Uniform probability per unit solid angle on the outward hemisphere:
  // p(Omega)=1/(2 pi), or p(mu)=1.  This contains no fitted constants and is
  // retained as an explicit comparison/control law.
  UniformHemisphere
};

constexpr EmissionAngleLaw kDefaultEmissionAngleLaw =
    EmissionAngleLaw::KnudsenCosine;

inline const char* EmissionAngleLawName(EmissionAngleLaw law) {
  switch (law) {
    case EmissionAngleLaw::KnudsenCosine: return "cosine-knudsen";
    case EmissionAngleLaw::UniformHemisphere: return "uniform-hemisphere";
  }
  throw std::invalid_argument("unsupported sputtering emission-angle law");
}

/**
 * Parse the stable, input-facing selector name.
 *
 * Exact lower-case names keep input receipts unambiguous.  Unknown values are
 * a typed error; the model never silently substitutes the default law.
 */
inline EmissionAngleLaw ParseEmissionAngleLaw(const std::string& name) {
  if (name == "cosine-knudsen") return EmissionAngleLaw::KnudsenCosine;
  if (name == "uniform-hemisphere")
    return EmissionAngleLaw::UniformHemisphere;
  throw std::invalid_argument("unsupported sputtering emission-angle selector: " +
                              name);
}

/**
 * Convert two independent U[0,1) samples into an outward unit direction.
 *
 * The output and `outward_normal` use the caller's Cartesian frame.  The
 * normal may have any positive magnitude and is normalized internally.  For
 * mu=cos(theta), inverse CDF sampling gives mu=sqrt(u) for Knudsen cosine and
 * mu=u for a uniform hemisphere; phi=2*pi*v.  The common AMPS surface helper
 * owns construction of the nonsingular right-handed tangent frame, including
 * normals parallel or antiparallel to a coordinate axis.
 *
 * Inputs are dimensionless.  A zero/nonfinite normal, a variate outside
 * [0,1), or an unknown enum value raises invalid_argument.  This operation
 * supplies direction only and adds no empirical energy-angle correlation.
 */
inline std::array<double, 3> DirectionFromUnitSamples(
    const std::array<double, 3>& outward_normal,
    double u, double v,
    EmissionAngleLaw law = kDefaultEmissionAngleLaw) {
  Exosphere::SurfaceSource::AngularLaw common_law;
  switch (law) {
    case EmissionAngleLaw::KnudsenCosine:
      common_law = Exosphere::SurfaceSource::AngularLaw::KnudsenCosine;
      break;
    case EmissionAngleLaw::UniformHemisphere:
      common_law = Exosphere::SurfaceSource::AngularLaw::UniformSolidAngle;
      break;
    default:
      throw std::invalid_argument("unsupported sputtering emission-angle law");
  }
  return Exosphere::SurfaceSource::DirectionFromUnitSamples(
      outward_normal, u, v, common_law);
}

/** Analytic first two moments of mu=cos(theta) for validation/closures. */
struct AngularMoments {
  double mean_mu;
  double mean_mu_squared;
};

inline AngularMoments GetAngularMoments(EmissionAngleLaw law) {
  switch (law) {
    case EmissionAngleLaw::KnudsenCosine:
      // integral_0^1 mu^n (2 mu) dmu = 2/(n+2).
      return AngularMoments{2.0 / 3.0, 1.0 / 2.0};
    case EmissionAngleLaw::UniformHemisphere:
      // integral_0^1 mu^n dmu = 1/(n+1).
      return AngularMoments{1.0 / 2.0, 1.0 / 3.0};
  }
  throw std::invalid_argument("unsupported sputtering emission-angle law");
}

namespace Detail {
inline void RequireFiniteNonnegative(double value, const char* name) {
  if (!(value >= 0.0) || !std::isfinite(value))
    throw std::invalid_argument(std::string(name) +
                                " must be finite and nonnegative");
}

inline void RequireFinitePositive(double value, const char* name) {
  if (!(value > 0.0) || !std::isfinite(value))
    throw std::invalid_argument(std::string(name) +
                                " must be finite and positive");
}

inline double CheckedProduct(double left, double right, const char* name) {
  const long double value = static_cast<long double>(left) * right;
  if (!std::isfinite(value) ||
      value > static_cast<long double>(std::numeric_limits<double>::max()) ||
      value < -static_cast<long double>(std::numeric_limits<double>::max()))
    throw std::overflow_error(std::string(name) + " is not representable");
  return static_cast<double>(value);
}

inline double CheckedSum(double left, double right, const char* name) {
  const long double value = static_cast<long double>(left) + right;
  if (!std::isfinite(value) ||
      value > static_cast<long double>(std::numeric_limits<double>::max()) ||
      value < -static_cast<long double>(std::numeric_limits<double>::max()))
    throw std::overflow_error(std::string(name) + " is not representable");
  return static_cast<double>(value);
}

inline double CheckedDivide(double numerator, double denominator,
                            const char* name) {
  const long double value =
      static_cast<long double>(numerator) / denominator;
  if (!std::isfinite(value) ||
      value > static_cast<long double>(std::numeric_limits<double>::max()) ||
      value < -static_cast<long double>(std::numeric_limits<double>::max()))
    throw std::overflow_error(std::string(name) + " is not representable");
  return static_cast<double>(value);
}
}  // namespace Detail

/** Deterministic source and numerical-event rates for one product species. */
struct SourceRateAccounting {
  double incident_projectiles_per_s;
  double physical_products_per_s;
  double physical_products_per_numerical_event;
  double numerical_events_per_s;
  double expected_numerical_events;
};

/**
 * Convert incident flux and yield to a physical rate and an event expectation.
 *
 * With incident number flux Gamma [projectiles m^-2 s^-1], physical area A
 * [m^2], and yield Y [products/projectile],
 *
 *   R_incident = Gamma A                         [projectiles s^-1]
 *   R_product  = Gamma A Y                       [products s^-1].
 *
 * If a numerical emission event represents W physical products, its mean
 * event rate is lambda=R_product/W [events s^-1] and an interval dt [s] has
 * expectation lambda*dt.  W must be the *effective* physical weight of that
 * event (base statistical weight times any individual weight correction).
 * Variable-weight scheduling cannot replace W by an arbitrary nominal value;
 * callers must use a qualified expectation or close the realized weights with
 * ProductLedger below.
 */
inline SourceRateAccounting ComputeSourceRateAccounting(
    double incident_flux_per_m2_s,
    double surface_area_m2,
    double yield_products_per_projectile,
    double physical_products_per_numerical_event,
    double interval_s) {
  Detail::RequireFiniteNonnegative(incident_flux_per_m2_s, "incident flux");
  Detail::RequireFiniteNonnegative(surface_area_m2, "surface area");
  Detail::RequireFiniteNonnegative(yield_products_per_projectile, "yield");
  Detail::RequireFinitePositive(physical_products_per_numerical_event,
                                "physical products per numerical event");
  Detail::RequireFiniteNonnegative(interval_s, "accounting interval");

  const double incident_rate = Detail::CheckedProduct(
      incident_flux_per_m2_s, surface_area_m2, "incident rate");
  const double product_rate = Detail::CheckedProduct(
      incident_rate, yield_products_per_projectile, "product rate");
  const double event_rate = Detail::CheckedDivide(
      product_rate, physical_products_per_numerical_event,
      "numerical event rate");
  const double expected_events = Detail::CheckedProduct(
      event_rate, interval_s, "expected numerical event count");
  return SourceRateAccounting{incident_rate, product_rate,
                              physical_products_per_numerical_event,
                              event_rate, expected_events};
}

/** Convert base statistical weight and an individual correction to W. */
inline double EffectivePhysicalWeight(double base_physical_weight,
                                      double individual_weight_correction) {
  Detail::RequireFinitePositive(base_physical_weight,
                                "base statistical particle weight");
  Detail::RequireFinitePositive(individual_weight_correction,
                                "individual weight correction");
  return Detail::CheckedProduct(base_physical_weight,
                                individual_weight_correction,
                                "effective physical particle weight");
}

/** Time-normalized physical and numerical ledger for one product identity. */
struct ProductRates {
  int product_species;
  double physical_products_per_s;
  double numerical_events_per_s;
  double mass_kg_per_s;
  std::array<double, 3> momentum_newton;
  double kinetic_power_watt;
};

/**
 * Conservative event ledger for one sputtered product species.
 *
 * RecordEvent consumes product mass m [kg], velocity v [m s^-1] in one stated
 * inertial/body-fixed Cartesian frame, base statistical weight W0 [physical
 * products/event], and individual correction c.  It accumulates
 *
 *   N += W0 c;  M += m W0 c;
 *   P += m W0 c v                              [kg m s^-1];
 *   K += (1/2) m W0 c |v|^2                    [J].
 *
 * Rates(duration) divides those extensive quantities by elapsed physical time
 * [s], yielding product rate [s^-1], mass rate [kg s^-1], momentum transfer
 * rate [N], and kinetic power [W].  Positive momentum is outward when the
 * caller's velocity and surface normal use that convention.  No binding,
 * recoil, surface-energy, or incident-projectile budget is silently added.
 */
class ProductLedger {
 public:
  explicit ProductLedger(int product_species)
      : product_species_(product_species), numerical_events_(0),
        physical_products_(0.0), mass_kg_(0.0),
        momentum_kg_m_per_s_{{0.0, 0.0, 0.0}}, kinetic_energy_joule_(0.0) {
    if (product_species < 0)
      throw std::invalid_argument("product species identity must be nonnegative");
  }

  void RecordEvent(int product_species, double product_mass_kg,
                   const std::array<double, 3>& velocity_m_per_s,
                   double base_physical_weight,
                   double individual_weight_correction = 1.0) {
    if (product_species != product_species_)
      throw std::invalid_argument("sputtering product species mismatch");
    Detail::RequireFinitePositive(product_mass_kg, "product mass");
    const double weight = EffectivePhysicalWeight(
        base_physical_weight, individual_weight_correction);

    long double speed_squared = 0.0L;
    for (int d = 0; d < 3; ++d) {
      if (!std::isfinite(velocity_m_per_s[d]))
        throw std::invalid_argument("product velocity must be finite");
      speed_squared += static_cast<long double>(velocity_m_per_s[d]) *
                       velocity_m_per_s[d];
    }
    if (!std::isfinite(speed_squared))
      throw std::overflow_error("product speed squared is not representable");

    const double event_mass = Detail::CheckedProduct(
        product_mass_kg, weight, "represented event mass");
    std::array<double, 3> event_momentum;
    std::array<double, 3> new_momentum;
    for (int d = 0; d < 3; ++d) {
      event_momentum[d] = Detail::CheckedProduct(
          event_mass, velocity_m_per_s[d], "represented event momentum");
      new_momentum[d] = Detail::CheckedSum(
          momentum_kg_m_per_s_[d], event_momentum[d],
          "cumulative represented momentum");
    }
    const long double event_energy_long =
        0.5L * static_cast<long double>(event_mass) * speed_squared;
    if (!std::isfinite(event_energy_long) ||
        event_energy_long > std::numeric_limits<double>::max())
      throw std::overflow_error("represented event kinetic energy is not representable");
    const double new_energy = Detail::CheckedSum(
        kinetic_energy_joule_, static_cast<double>(event_energy_long),
        "cumulative represented kinetic energy");
    const double new_products = Detail::CheckedSum(
        physical_products_, weight, "cumulative represented products");
    const double new_mass = Detail::CheckedSum(
        mass_kg_, event_mass, "cumulative represented mass");
    if (numerical_events_ == std::numeric_limits<std::uint64_t>::max())
      throw std::overflow_error("cumulative numerical event count is not representable");
    const std::uint64_t new_events = numerical_events_ + 1;

    // Commit only after every candidate quantity has passed validation so a
    // rejected event cannot leave a partially updated conservation ledger.
    physical_products_ = new_products;
    mass_kg_ = new_mass;
    momentum_kg_m_per_s_ = new_momentum;
    kinetic_energy_joule_ = new_energy;
    numerical_events_ = new_events;
  }

  ProductRates Rates(double duration_s) const {
    Detail::RequireFinitePositive(duration_s, "ledger duration");
    ProductRates rates;
    rates.product_species = product_species_;
    rates.physical_products_per_s = Detail::CheckedDivide(
        physical_products_, duration_s, "physical product rate");
    rates.numerical_events_per_s = Detail::CheckedDivide(
        static_cast<double>(numerical_events_), duration_s,
        "numerical event rate");
    rates.mass_kg_per_s = Detail::CheckedDivide(
        mass_kg_, duration_s, "product mass rate");
    for (int d = 0; d < 3; ++d)
      rates.momentum_newton[d] = Detail::CheckedDivide(
          momentum_kg_m_per_s_[d], duration_s, "product momentum rate");
    rates.kinetic_power_watt = Detail::CheckedDivide(
        kinetic_energy_joule_, duration_s, "product kinetic power");
    return rates;
  }

  int product_species() const { return product_species_; }
  std::uint64_t numerical_events() const { return numerical_events_; }
  double physical_products() const { return physical_products_; }
  double mass_kg() const { return mass_kg_; }
  const std::array<double, 3>& momentum_kg_m_per_s() const {
    return momentum_kg_m_per_s_;
  }
  double kinetic_energy_joule() const { return kinetic_energy_joule_; }

 private:
  int product_species_;
  std::uint64_t numerical_events_;
  double physical_products_;
  double mass_kg_;
  std::array<double, 3> momentum_kg_m_per_s_;
  double kinetic_energy_joule_;
};

}  // namespace Source
}  // namespace Sputtering

#endif

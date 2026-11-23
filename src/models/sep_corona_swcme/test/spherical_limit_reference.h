#ifndef SEP_CORONA_SWCME_SPHERICAL_LIMIT_REFERENCE_H
#define SEP_CORONA_SWCME_SPHERICAL_LIMIT_REFERENCE_H

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

// Independent TEST-ONLY one-dimensional reference for the exact controlled
// Stage13 isothermal/PFSS-polar/Parker laws. No production wind, PFSS, EOS,
// shock, gradient, trajectory or critical-point-selection routine is called.
// It is not the distinct Stage12 adiabatic constant-speed outer-wind profile.
namespace spherical_limit_reference {
using Real = long double;
using Vector = std::array<Real, 3>;
using Matrix = std::array<Real, 9>;
constexpr Real Gm = 1.32712440018e20L, Boltzmann = 1.380649e-23L;
constexpr Real ProtonMass = 1.67262192369e-27L, AlphaMass = 6.6446573357e-27L;
constexpr Real ElectronMass = 9.1093837015e-31L, Mu0 = 1.25663706212e-6L;
// Frozen BEFORE this added comparison runs. Existing Stage13 gates retain
// their original, independent limits. Scale zeros with the full physical
// vector/tensor scale, never max(1, a tiny dimensionful quantity).
constexpr Real PrimitiveTolerance = 3e-7L, GradientTolerance = 3e-6L;
constexpr Real PositionToleranceM = .2L, DeformationTolerance = 2e-6L;
struct Parameters {
  Real solarRadius = 0, sourceRadius = 0, referenceRadius = 0, maximumRadius = 0;
  Real electronDensityAtReference = 0, electronTemperature = 0, protonTemperature = 0;
  Real alphaTemperature = 0, alphaRatio = 0, rotation = 0, gamma = 0,
       dipoleCoefficient = 0;
  bool includeElectronMass = false;
};
struct State {
  Real radius = 0, speed = 0, speedDerivative = 0, density = 0, pressure = 0,
       electronDensity = 0;
  Real protonDensity = 0, magneticZ = 0, densityDerivative = 0, pressureDerivative = 0;
  Vector velocity{}, magnetic{}, densityGradient{}, pressureGradient{};
  Matrix velocityGradient{}, magneticGradient{};
};
struct Flow {
  Real radius = 0, jacobian = 0;
  Matrix deformation{};
};
struct Jump {
  Real compression = 0, density = 0, pressure = 0, velocityZ = 0, magneticZ = 0,
       fastMach = 0;
};
class Profile {
 public:
  explicit Profile(Parameters inputs) : p_(inputs) {
    if (!(p_.solarRadius > 0 && p_.sourceRadius > p_.solarRadius &&
          p_.maximumRadius > p_.sourceRadius && p_.referenceRadius >= p_.sourceRadius &&
          p_.referenceRadius <= p_.maximumRadius && p_.electronDensityAtReference > 0 &&
          p_.electronTemperature > 0 && p_.protonTemperature > 0 &&
          p_.alphaTemperature > 0 && p_.alphaRatio >= 0 && p_.gamma > 1 &&
          p_.dipoleCoefficient != 0))
      throw std::runtime_error(
          "independent1D reference requires complete explicit SI inputs");
    // Quasineutral e/p/alpha EOS, independently written from species sums.
    const Real charge = 1 + 2 * p_.alphaRatio;
    densityPerElectron_ = (ProtonMass + p_.alphaRatio * AlphaMass) / charge +
                          (p_.includeElectronMass ? ElectronMass : 0);
    pressurePerElectron_ =
        Boltzmann *
        (p_.electronTemperature +
         (p_.protonTemperature + p_.alphaRatio * p_.alphaTemperature) / charge);
    cs2_ = pressurePerElectron_ / densityPerElectron_;
    rc_ = Gm / (2 * cs2_);
    massFlux_ = densityPerElectron_ * p_.electronDensityAtReference *
                Speed(p_.referenceRadius) * p_.referenceRadius * p_.referenceRadius;
    // Y_10=sqrt(3/(4pi))*cos(theta). This is a flux-balanced dipole,
    // not a fictitious magnetic monopole. The polar ray is a valid open tube.
    coefficient_ = p_.dipoleCoefficient * std::sqrt(3 / (4 * std::acos(-1.L)));
    denominator_ = 1 + 2 * std::pow(p_.sourceRadius / p_.solarRadius, 3);
  }
  Real CriticalRadius() const {
    return rc_;
  }
  Real Speed(Real radius) const {
    if (!(radius > 0))
      throw std::runtime_error("invalid independent1D Parker radius");
    // Solve the unique transonic Parker branch in log(M_iso^2), by bisection
    // of exp(q)-q=A(r). No sampled critical candidate or table supplies it.
    const Real delta = (radius - rc_) / rc_;
    const Real rhs = 1 + 4 * (std::log1p(delta) - delta / (1 + delta));
    if (rhs < 1)
      throw std::runtime_error("independent Parker invariant is numerically invalid");
    if (rhs == 1)
      return std::sqrt(cs2_);
    Real lo = radius < rc_ ? -rhs - 2 : 0;
    Real hi = radius < rc_ ? 0 : std::log(rhs) + 2;
    for (unsigned i = 0; i < 180; ++i) {
      const Real middle = (lo + hi) / 2;
      if (middle == lo || middle == hi)
        break;
      const Real residual = std::exp(middle) - middle - rhs;
      if (radius < rc_ ? (residual > 0) : (residual < 0))
        lo = middle;
      else
        hi = middle;
    }
    return std::sqrt(cs2_) * std::exp((lo + hi) / 4);
  }
  State At(Real radius, int hemisphere) const {
    if ((hemisphere != 1 && hemisphere != -1) || radius < p_.solarRadius ||
        radius > p_.maximumRadius)
      throw std::runtime_error(
          "independent1D polar sample is outside declared support");
    State s;
    s.radius = radius;
    s.speed = Speed(radius);
    // A separate fourth-order scalar difference avoids the singular 0/0
    // implicit derivative at the regular critical point. Long-double h/r=1e-4
    // resolves roundoff while making its smooth truncation O(1e-16).
    const Real h = radius * 1e-4L;
    s.speedDerivative = (Speed(radius - 2 * h) - 8 * Speed(radius - h) +
                         8 * Speed(radius + h) - Speed(radius + 2 * h)) /
                        (12 * h);
    Real dlogB = 0, transverseStrain = 0, angularRate = 0, magneticDerivative = 0;
    const Real sourceField = 3 * coefficient_ / denominator_;
    if (radius < p_.sourceRadius) {
      const Real q = std::pow(p_.sourceRadius / radius, 3);
      s.magneticZ = coefficient_ * (1 + 2 * q) / denominator_;
      magneticDerivative = -6 * coefficient_ * q / (denominator_ * radius);
      dlogB = magneticDerivative / s.magneticZ;
      s.density = massFlux_ * std::abs(s.magneticZ) /
                  (p_.sourceRadius * p_.sourceRadius * std::abs(sourceField) * s.speed);
      transverseStrain = -s.speed * dlogB / 2;
      angularRate = p_.rotation;
    } else {
      s.magneticZ = sourceField * std::pow(p_.sourceRadius / radius, 2);
      magneticDerivative = -2 * s.magneticZ / radius;
      dlogB = -2 / radius;
      s.density = massFlux_ / (s.speed * radius * radius);
      transverseStrain = s.speed / radius;
      angularRate = p_.rotation * std::pow(p_.sourceRadius / radius, 2);
    }
    s.electronDensity = s.density / densityPerElectron_;
    s.protonDensity = s.electronDensity / (1 + 2 * p_.alphaRatio);
    s.pressure = s.electronDensity * pressurePerElectron_;
    s.densityDerivative = s.density * (dlogB - s.speedDerivative / s.speed);
    s.pressureDerivative = cs2_ * s.densityDerivative;
    s.velocity = {{0, 0, hemisphere * s.speed}};
    s.magnetic = {{0, 0, s.magneticZ}};
    s.densityGradient = {{0, 0, hemisphere * s.densityDerivative}};
    s.pressureGradient = {{0, 0, hemisphere * s.pressureDerivative}};
    s.velocityGradient = {{transverseStrain, -angularRate, 0, angularRate,
                           transverseStrain, 0, 0, 0, s.speedDerivative}};
    // Transverse derivatives follow the independent analytic dipole potential
    // below Rss, or the signed radial-flux/Parker field above it. Thus a 1D
    // tube also predicts the complete local Cartesian tensors at its axis.
    const Real diagonal = -hemisphere * magneticDerivative / 2;
    const Real winding =
        radius < p_.sourceRadius
            ? 0
            : p_.rotation * (1 - std::pow(p_.sourceRadius / radius, 2)) / s.speed;
    s.magneticGradient = {{diagonal, hemisphere * s.magneticZ * winding, 0,
                           -hemisphere * s.magneticZ * winding, diagonal, 0, 0, 0,
                           hemisphere * magneticDerivative}};
    return s;
  }
  Flow Advance(Real initial, Real duration) const {
    if (duration < 0 || initial < p_.solarRadius || initial >= p_.maximumRadius)
      throw std::runtime_error("independent1D forward flow has invalid support/time");
    Flow out;
    out.radius = initial;
    if (duration > 0) {
      Real lo = initial,
           hi = std::min(p_.maximumRadius, initial + 2 * Speed(initial) * duration);
      while (Travel(initial, hi) < duration) {
        if (hi == p_.maximumRadius)
          throw std::runtime_error("independent1D flow exceeds actual guard");
        hi = std::min(p_.maximumRadius, initial + 2 * (hi - initial));
      }
      for (unsigned i = 0; i < 100; ++i) {
        const Real middle = (lo + hi) / 2;
        if (middle == lo || middle == hi)
          break;
        if (Travel(initial, middle) < duration)
          lo = middle;
        else
          hi = middle;
      }
      out.radius = (lo + hi) / 2;
    }
    const auto a = At(initial, 1), b = At(out.radius, 1);
    const Real transverse = std::sqrt(std::abs(a.magneticZ / b.magneticZ));
    // On the polar material tube the transverse 2x2 generators commute:
    // isotropic area stretch times a genuine corotation/angular-momentum
    // rotation, with radial stretch U_end/U_begin from the travel-time map.
    const Real angle = Rotation(initial, out.radius);
    out.deformation = {{transverse * std::cos(angle), -transverse * std::sin(angle), 0,
                        transverse * std::sin(angle), transverse * std::cos(angle), 0,
                        0, 0, b.speed / a.speed}};
    out.jacobian = transverse * transverse * b.speed / a.speed;
    return out;
  }
  Jump ParallelJump(Real radius, Real frontSpeed, int hemisphere) const {
    const auto s = At(radius, hemisphere);
    const Real incoming = frontSpeed - s.speed;
    const Real sound = std::sqrt(p_.gamma * s.pressure / s.density);
    const Real alfven = std::abs(s.magneticZ) / std::sqrt(Mu0 * s.density);
    if (!(incoming > std::max(sound, alfven)))
      throw std::runtime_error("independent parallel front is not fast");
    const Real mach2 = incoming * incoming / (sound * sound);
    Jump out;
    out.compression = (p_.gamma + 1) * mach2 / ((p_.gamma - 1) * mach2 + 2);
    if (incoming * incoming / (alfven * alfven) <= out.compression)
      throw std::runtime_error(
          "independent gas-parallel branch is not a regular fast MHD shock");
    out.density = s.density * out.compression;
    out.pressure =
        s.pressure * (2 * p_.gamma * mach2 - (p_.gamma - 1)) / (p_.gamma + 1);
    out.velocityZ = hemisphere * (frontSpeed - incoming / out.compression);
    out.magneticZ = s.magneticZ;
    out.fastMach = incoming / std::max(sound, alfven);
    return out;
  }

 private:
  template <class F>
  Real Integrate(const F& function, Real a, Real b) const {
    // Composite 32-panel Simpson quadrature in the scalar radial coordinate,
    // independent of the production Cartesian RK4/variational integration.
    // The material steps used here span <1% in r; a 64-panel duplicate checks
    // that reference quadrature itself is resolved before it can grade code.
    const auto sum = [&](unsigned n) {
      const Real h = (b - a) / n;
      Real value = function(a) + function(b);
      for (unsigned i = 1; i < n; ++i)
        value += (i % 2 ? 4 : 2) * function(a + i * h);
      return value * h / 3;
    };
    const Real coarse = sum(32), fine = sum(64);
    if (std::abs(fine - coarse) > 1e-12L * std::max(std::abs(fine), 1e-30L))
      throw std::runtime_error("independent1D reference quadrature is unresolved");
    return fine;
  }
  Real Travel(Real a, Real b) const {
    return Integrate([&](Real r) { return 1 / Speed(r); }, a, b);
  }
  Real Rotation(Real a, Real b) const {
    if (a == b)
      return 0;
    const auto integrate = [&](Real lo, Real hi, bool exterior) {
      return Integrate(
          [&](Real r) {
            return p_.rotation * (exterior ? std::pow(p_.sourceRadius / r, 2) : 1) /
                   Speed(r);
          },
          lo, hi);
    };
    if (a < p_.sourceRadius && b > p_.sourceRadius)
      return integrate(a, p_.sourceRadius, false) + integrate(p_.sourceRadius, b, true);
    return integrate(a, b, a >= p_.sourceRadius);
  }
  Parameters p_;
  Real densityPerElectron_ = 0, pressurePerElectron_ = 0, cs2_ = 0, rc_ = 0;
  Real massFlux_ = 0, coefficient_ = 0, denominator_ = 0;
};
}  // namespace spherical_limit_reference
#endif

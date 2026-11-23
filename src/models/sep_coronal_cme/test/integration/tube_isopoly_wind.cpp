#include "sep_coronal_cme/tube_isopoly_wind.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>

namespace C = SEP::CoronalCME;
namespace {
using L = long double;
constexpr L Kb = 1.380649e-23L;
// Frozen before execution. The references use long-double algebra/root solves
// or independent fourth-order differences, never a production solver helper.
constexpr L ValueTolerance = 2e-5L, DerivativeTolerance = 4e-5L;
void Require(bool value, const std::string& message) {
  if (!value)
    throw std::runtime_error(message);
}
void Match(L actual, L expected, L scale, L tolerance, const std::string& label) {
  Require(std::isfinite(actual) && std::isfinite(expected) && scale > 0 &&
              std::abs(actual - expected) <= tolerance * scale,
          label);
}
C::IsoPolySpecies Species(double gamma = 1) {
  // Manufactured unit-sound-speed gas. This is a numerical fixture, not a
  // proton/electron observational profile or a new scientific default.
  return {static_cast<double>(1 / Kb), 1, 1, gamma, .1, .2, 5. / 3.};
}
C::IsoPolyTubeInput Radial(double gamma = 1) {
  C::IsoPolyTubeInput in;
  // Dimensionful manufactured GM=2 m^3/s^2 and A=(r/1m)^2 m^2.
  // In the isothermal limit c_iso^2=1 m^2/s^2, giving r_c=1 m.
  in.geometry = [](double r) {
    return SEP::Core::Result<C::TubeGeometryPoint>::Success(
        {r, r, r * r, -2 / r, 2 / r, 2 / (r * r)});
  };
  auto species = Species(gamma);
  species.referenceTemperatureK = 1 / gamma;
  in.species = {species};
  in.beginM = .5;
  in.endM = 8;
  in.massFluxKgPerS = 1;
  in.criticalScanIntervals = 384;
  in.maximumIntegrationAttempts = 65536;
  in.relativeTolerance = 1e-10;
  return in;
}
C::IsoPolyTubeWind Solve(const C::IsoPolyTubeInput& input) {
  const auto result = C::IsoPolyTubeWind::Solve(input);
  Require(result.ok(), result.status.message);
  return result.value;
}
C::IsoPolyTubeSample At(const C::IsoPolyTubeWind& wind, double x) {
  const auto sample = wind.Evaluate(x);
  Require(sample.ok(), sample.status.message);
  return sample.value;
}
L IsothermalSpeed(L radius, L soundSquared, L critical) {
  const L epsilon = (radius - critical) / critical;
  const L rhs = 1 + 4 * (std::log1p(epsilon) - epsilon / (1 + epsilon));
  Require(rhs >= 1, "independent Parker invariant is invalid");
  if (rhs == 1)
    return std::sqrt(soundSquared);
  L lo = radius < critical ? -rhs - 2 : 0,
    hi = radius < critical ? 0 : std::log(rhs) + 2;
  for (unsigned i = 0; i < 180; ++i) {
    const L mid = (lo + hi) / 2;
    if (mid == lo || mid == hi)
      break;
    const L f = std::exp(mid) - mid - rhs;
    if (radius < critical ? (f > 0) : (f < 0))
      lo = mid;
    else
      hi = mid;
  }
  return std::sqrt(soundSquared) * std::exp((lo + hi) / 4);
}
void IPW01() {
  const auto input = Radial();
  const auto wind = Solve(input);
  const L a2 =
      Kb * input.species[0].particlesPerKg * input.species[0].referenceTemperatureK;
  const L critical = 1 / a2;
  Match(wind.CriticalCoordinateM(), critical, critical, ValueTolerance,
        "radial isothermal critical point");
  L maximum = 0;
  for (double r : {.5, .7, .9, 1., 1.2, 2., 4., 8.}) {
    const auto s = At(wind, r);
    const L speed = IsothermalSpeed(r, a2, critical);
    const L density = 1 / (speed * r * r);
    Match(s.wind.speedMPerS, speed, speed, ValueTolerance,
          "independent isothermal Parker speed");
    Match(s.wind.densityKgM3, density, density, ValueTolerance,
          "independent isothermal density normalization");
    Match(s.wind.pressurePa, a2 * density, a2 * density, ValueTolerance,
          "independent isothermal pressure");
    Match(s.wind.densityKgM3 * s.wind.speedMPerS * s.wind.areaM2, 1, 1, 2e-13,
          "mass flux invariant");
    maximum = std::max(maximum, std::abs(s.wind.speedMPerS / speed - 1));
  }
  const auto initial = At(wind, input.beginM), final = At(wind, input.endM);
  Require(initial.wind.speedMPerS < std::sqrt(initial.wind.soundSpeedSquaredM2S2) &&
              final.wind.speedMPerS > std::sqrt(final.wind.soundSpeedSquaredM2S2),
          "solution did not span the declared subsonic/supersonic branches");
  std::cout << "max_radial_isothermal_speed_relative=" << maximum << '\n';
}
struct PolytropicReference {
  L gamma, K, critical, speedCritical, energy;
  explicit PolytropicReference(const C::IsoPolyTubeInput& input) {
    const auto& s = input.species[0];
    gamma = s.polytropicIndex;
    K = Kb * s.particlesPerKg * s.referenceTemperatureK /
        std::pow(static_cast<L>(s.referenceDensityKgM3), gamma - 1);
    // Mass flux q=1 and GM=2 give rho_c=r_c^(-3/2), a_c^2=1/r_c.
    critical = std::pow(1 / (gamma * K), 1 / (1 - 1.5L * (gamma - 1)));
    speedCritical = 1 / std::sqrt(critical);
    energy = speedCritical * speedCritical * (1 / (gamma - 1) - 1.5L);
  }
  L Speed(L r) const {
    if (std::abs(r / critical - 1) < 1e-7L) {
      // Analytic regular sonic slope; this local expansion avoids subtracting
      // two nearly equal Bernoulli minima at a represented critical radius.
      const L g = gamma - 1;
      const L slope = (-2 * g + std::sqrt(2 * (5 - 3 * gamma))) / (gamma + 1) *
                      speedCritical / critical;
      return speedCritical + slope * (r - critical);
    }
    const L area = r * r;
    const L sonic =
        std::pow(gamma * K * std::pow(1 / area, gamma - 1), 1 / (gamma + 1));
    const auto residual = [&](L logu) {
      const L u = std::exp(logu), rho = 1 / (u * area);
      return u * u / 2 + gamma * K * std::pow(rho, gamma - 1) / (gamma - 1) - 2 / r -
             energy;
    };
    L lo = std::log(sonic) + (r < critical ? -40 : 0),
      hi = std::log(sonic) + (r < critical ? 0 : 40);
    Require(residual(std::log(sonic)) < 0,
            "independent polytropic branch has no Bernoulli roots");
    for (unsigned i = 0; i < 180; ++i) {
      const L middle = (lo + hi) / 2;
      if (middle == lo || middle == hi)
        break;
      if (r < critical ? (residual(middle) > 0) : (residual(middle) < 0))
        lo = middle;
      else
        hi = middle;
    }
    return std::exp((lo + hi) / 2);
  }
};
void IPW02() {
  const auto input = Radial(1.2);
  const auto wind = Solve(input);
  const PolytropicReference reference(input);
  Match(wind.CriticalCoordinateM(), reference.critical, reference.critical,
        ValueTolerance, "independent polytropic sonic radius");
  L maximum = 0;
  for (double r : {.5, .7, .9, 1., 1.2, 2., 4., 8.}) {
    const auto s = At(wind, r);
    const L u = reference.Speed(r), rho = 1 / (u * r * r);
    const L pressure = reference.K * std::pow(rho, reference.gamma);
    Match(s.wind.speedMPerS, u, u, ValueTolerance,
          "independent variable-temperature Bernoulli speed");
    Match(s.wind.densityKgM3, rho, rho, ValueTolerance,
          "independent polytropic density");
    Match(s.wind.pressurePa, pressure, pressure, ValueTolerance,
          "independent polytropic pressure");
    const L energy = s.wind.speedMPerS * s.wind.speedMPerS / 2 +
                     reference.gamma * s.wind.pressurePa /
                         ((reference.gamma - 1) * s.wind.densityKgM3) -
                     2 / r;
    Match(energy, reference.energy, std::abs(reference.energy), ValueTolerance,
          "independent polytropic Bernoulli invariant");
    maximum = std::max(maximum, std::abs(s.wind.speedMPerS / u - 1));
  }
  Require(At(wind, .5).temperatureK[0] > At(wind, 8).temperatureK[0],
          "nonisothermal limit was replaced by an isothermal wind");
  std::cout << "max_radial_polytropic_speed_relative=" << maximum << '\n';
}
struct ThermalReference {
  L temperature, pressure, pressureDerivative, heat, explicitTerm;
};
ThermalReference Thermal(const C::IsoPolySpecies& species, L x, L rho, L rhoPrime,
                         L u) {
  L weight = 0, derivative = 0;
  if (species.polytropicIndex != 1 && x > species.blendBeginM) {
    const L width = species.blendEndM - species.blendBeginM;
    const L z = (x - species.blendBeginM) / width;
    if (z >= 1)
      weight = 1;
    else {
      weight = 10 * std::pow(z, 3) - 15 * std::pow(z, 4) + 6 * std::pow(z, 5);
      derivative = (30 * z * z - 60 * z * z * z + 30 * z * z * z * z) / width;
    }
  }
  const L q = species.polytropicIndex - 1;
  const L temperature = species.referenceTemperatureK *
                        std::pow(rho / species.referenceDensityKgM3, q * weight);
  const L pressure = Kb * species.particlesPerKg * rho * temperature;
  const L explicitTerm = q * derivative * std::log(rho / species.referenceDensityKgM3);
  const L dp = pressure * ((1 + q * weight) * rhoPrime / rho + explicitTerm);
  return {temperature, pressure, dp,
          u / (species.adiabaticIndex - 1) *
              (dp - species.adiabaticIndex * pressure * rhoPrime / rho),
          pressure * explicitTerm};
}
void IPW03() {
  auto input = Radial();
  input.species = {Species(1.1), Species(1.3)};
  for (auto& s : input.species) {
    s.particlesPerKg *= .5;
    s.blendBeginM = 1.5;
    s.blendEndM = 2.5;
  }
  const auto wind = Solve(input);
  L maxMomentum = 0, maxExplicit = 0;
  for (double x : {.7, 1.1, 1.5, 1.7, 2., 2.3, 2.5, 3., 5., 7.}) {
    const auto sample = At(wind, x);
    Require(sample.temperatureK.size() == input.species.size() &&
                sample.partialPressurePa.size() == input.species.size() &&
                sample.netHeatWPerM3.size() == input.species.size(),
            "thermal owner omitted a configured species output slot");
    const L u = sample.wind.speedMPerS, rho = sample.wind.densityKgM3;
    // Differentiate actual DENSE values with an independent fourth-order
    // stencil; never use the solver's ODE evaluator to grade its momentum.
    const double h = 1e-4 * x;
    const auto m2 = At(wind, x - 2 * h), m1 = At(wind, x - h), p1 = At(wind, x + h),
               p2 = At(wind, x + 2 * h);
    const L du = (static_cast<L>(m2.wind.speedMPerS) - 8 * m1.wind.speedMPerS +
                  8 * p1.wind.speedMPerS - p2.wind.speedMPerS) /
                 (12 * h);
    const L drho = -rho * (du / u + 2 / x);
    L pressure = 0, dp = 0;
    for (std::size_t i = 0; i < input.species.size(); ++i) {
      const auto t = Thermal(input.species[i], x, rho, drho, u);
      Match(sample.temperatureK[i], t.temperature, t.temperature, ValueTolerance,
            "independent log-temperature blend");
      Match(sample.partialPressurePa[i], t.pressure, t.pressure, ValueTolerance,
            "independent species pressure");
      Match(sample.netHeatWPerM3[i], t.heat, u * t.pressure / x, DerivativeTolerance,
            "independent nonadiabatic species energy budget");
      pressure += t.pressure;
      dp += t.pressureDerivative;
      maxExplicit = std::max(maxExplicit, std::abs(t.explicitTerm));
    }
    Match(sample.wind.pressurePa, pressure, pressure, ValueTolerance,
          "sum of species thermal pressures");
    Match(sample.speedDerivativePerS, du, u / x, DerivativeTolerance,
          "dense speed derivative");
    Match(sample.densityDerivativeKgM4, drho, rho / x, DerivativeTolerance,
          "density derivative from independent mass balance");
    Match(sample.pressureDerivativePaPerM, dp, pressure / x, DerivativeTolerance,
          "explicit thermal-gradient pressure derivative");
    const L residual = u * du + dp / rho + 2 / (x * x);
    const L scale = std::abs(u * du) + std::abs(dp / rho) + 2 / (x * x);
    maxMomentum = std::max(maxMomentum, std::abs(residual) / scale);
    Require(std::abs(residual) <= DerivativeTolerance * scale,
            "nonisothermal momentum residual omits explicit spatial thermodynamics");
  }
  Require(
      maxExplicit > .01L,
      "iso-poly fixture did not exercise a resolved explicit thermal-gradient term");
  std::cout << "max_isopoly_momentum_relative=" << maxMomentum
            << " max_explicit_pressure_gradient=" << maxExplicit << '\n';
}
L Height(L x, L epsilon) {
  return .1L * ((x - 1) * (x - 1) * (x - 3) * (x - 3) + epsilon * (x - 1) * (x - 1));
}
L HeightDerivative(L x, L epsilon) {
  return .2L * (x - 1) * ((x - 3) * (2 * x - 4) + epsilon);
}
C::IsoPolyTubeInput Multiple(double epsilon) {
  // For a constant-temperature nozzle, 0.5*u^2-c^2*log(u) equals
  // H=c^2*log(A)-Phi plus its Bernoulli constant. The left side has its
  // minimum at the sonic speed, so a global sonic branch must use the GLOBAL
  // minimum of H. For epsilon=0.1 the stationary points are exactly
  // 1 and (5 +/- sqrt(0.8))/2: the inner minimum H(1)=0 is unique, while
  // the outer positive minimum is only a local candidate. Selecting that
  // outer root would require a negative Bernoulli inversion at s=1.
  // These are manufactured smooth geometry/potential data, not a solar
  // gravitational profile or a newly selected physical coronal closure.
  auto input = Radial();
  input.beginM = 0;
  input.endM = 4;
  input.massFluxKgPerS = std::exp(1.);
  input.geometry = [epsilon](double x) {
    return SEP::Core::Result<C::TubeGeometryPoint>::Success(
        {x, 10 + x, std::exp(x), x - static_cast<double>(Height(x, epsilon)), 1,
         1 - static_cast<double>(HeightDerivative(x, epsilon))});
  };
  return input;
}
L NozzleSpeed(L height, bool subsonic) {
  // 0.5(u^2-1)-log(u)=H-H_critical. Both roots are independently bracketed.
  if (height == 0)
    return 1;
  const auto f = [&](L q) { return (std::exp(2 * q) - 1) / 2 - q - height; };
  L lo = subsonic ? -height - 2 : 0, hi = subsonic ? 0 : std::log(2 * height + 2) + 2;
  for (unsigned i = 0; i < 180; ++i) {
    const L mid = (lo + hi) / 2;
    if (mid == lo || mid == hi)
      break;
    if (subsonic ? (f(mid) > 0) : (f(mid) < 0))
      lo = mid;
    else
      hi = mid;
  }
  return std::exp((lo + hi) / 2);
}
void IPW04() {
  const auto input = Multiple(.1);
  const auto wind = Solve(input);
  const L outer = (5 + std::sqrt(.8L)) / 2;
  Require(Height(1, .1L) == 0 && Height(outer, .1L) > 0,
          "independent global critical-height certificate invalid");
  Match(wind.CriticalCoordinateM(), 1, 1, ValueTolerance,
        "outermost positive critical candidate was selected instead of the unique "
        "global branch");
  for (double x : {0., .5, 1., 1.5, 2., static_cast<double>(outer), 3.5, 4.}) {
    const auto p = At(wind, x);
    const L expected = NozzleSpeed(Height(x, .1L), x < 1);
    Match(p.wind.speedMPerS, expected, expected, ValueTolerance,
          "independent multiple-critical Bernoulli branch");
  }
}
void IPW05() {
  auto input = Multiple(0);
  // Constant area makes the two equal potential barriers EXACTLY equal for
  // any represented kB*particlesPerKg value. It removes an unintended EOS
  // roundoff tilt of the earlier exp(s) area from this degeneracy control.
  input.geometry = [](double x) {
    return SEP::Core::Result<C::TubeGeometryPoint>::Success(
        {x, 10 + x, 1, -static_cast<double>(Height(x, 0)), 0,
         -static_cast<double>(HeightDerivative(x, 0))});
  };
  input.massFluxKgPerS = 1;
  const auto result = C::IsoPolyTubeWind::Solve(input);
  if (result.ok()) {
    const auto a = At(result.value, 1), b = At(result.value, 3);
    std::cout << "unexpected_degenerate_critical=" << result.value.CriticalCoordinateM()
              << " mach_at_equal_barriers="
              << a.wind.speedMPerS / std::sqrt(a.wind.soundSpeedSquaredM2S2) << ','
              << b.wind.speedMPerS / std::sqrt(b.wind.soundSpeedSquaredM2S2) << '\n';
  }
  Require(!result.ok(),
          "two equal global sonic minima were admitted as one unique full-support "
          "transonic branch");
  std::cout << "degenerate_global_branch_rejected=" << result.status.message << '\n';
  // The physical state is still isothermal when a non-unit outer index is
  // configured but its entire transition lies beyond the solved support.
  // Admission must use the actual closure on this interval, not just gamma.
  input.species[0].polytropicIndex = 1.2;
  input.species[0].blendBeginM = 10;
  input.species[0].blendEndM = 11;
  Require(!C::IsoPolyTubeWind::Solve(input).ok(),
          "dormant polytropic band bypassed the equal-isothermal-minimum rejection");
  // An equal minimum at the declared inner endpoint would make that endpoint
  // sonic, contradicting the required strictly subsonic initial branch. It
  // must not be hidden because the critical scan stores interior roots only.
  input.endM = 3;
  input.geometry = [](double x) {
    const double h = .1 * x * x * (x - 1) * (x - 1);
    const double derivative = .2 * x * (x - 1) * (2 * x - 1);
    return SEP::Core::Result<C::TubeGeometryPoint>::Success(
        {x, 10 + x, 1, -h, 0, -derivative});
  };
  Require(!C::IsoPolyTubeWind::Solve(input).ok(),
          "tied endpoint sonic minimum was admitted as a strictly subsonic boundary");
}
void IPW06() {
  const auto original = Radial();
  auto invalid = original;
  invalid.geometry = {};
  Require(!C::IsoPolyTubeWind::Solve(invalid).ok(), "missing geometry accepted");
  invalid = original;
  invalid.species.clear();
  Require(!C::IsoPolyTubeWind::Solve(invalid).ok(),
          "missing thermal composition accepted");
  invalid = original;
  invalid.massFluxKgPerS = 0;
  Require(!C::IsoPolyTubeWind::Solve(invalid).ok(),
          "missing mass normalization accepted");
  invalid = original;
  invalid.species[0].referenceTemperatureK = 0;
  Require(!C::IsoPolyTubeWind::Solve(invalid).ok(), "missing temperature accepted");
  invalid = original;
  invalid.species[0].polytropicIndex = 1.2;
  invalid.species[0].blendEndM = invalid.species[0].blendBeginM;
  Require(!C::IsoPolyTubeWind::Solve(invalid).ok(),
          "unresolved thermal transition accepted");
  invalid = original;
  invalid.geometry = [](double x) {
    return SEP::Core::Result<C::TubeGeometryPoint>::Success(
        {x, x, -1, -2 / x, 2 / x, 2 / (x * x)});
  };
  Require(!C::IsoPolyTubeWind::Solve(invalid).ok(), "nonpositive tube area repaired");
  invalid = original;
  invalid.geometry = [](double x) {
    return SEP::Core::Result<C::TubeGeometryPoint>::Success(
        {x + .01, x, x * x, -2 / x, 2 / x, 2 / (x * x)});
  };
  Require(!C::IsoPolyTubeWind::Solve(invalid).ok(),
          "geometry coordinate mismatch accepted");
  invalid = original;
  invalid.maximumIntegrationAttempts = 1;
  Require(!C::IsoPolyTubeWind::Solve(invalid).ok(),
          "exhausted numerical work silently returned a profile");
}
void IPW07() {
  auto input = Radial();
  const auto geometry = input.geometry;
  input.geometry = [geometry](double x) {
    if (x > 2 && x < 3)
      return SEP::Core::Result<C::TubeGeometryPoint>::Failure(
          SEP::Core::StatusCode::DataIntegrityFailure,
          "explicit missing open-tube support");
    return geometry(x);
  };
  Require(!C::IsoPolyTubeWind::Solve(input).ok(),
          "missing interior physical support was filled/extrapolated");
  const auto wind = Solve(Radial());
  Require(!wind.Evaluate(.5 - 1e-8).ok() && !wind.Evaluate(8 + 1e-8).ok(),
          "solved profile silently extrapolated outside its support");
  Require(!wind.Evaluate(std::numeric_limits<double>::quiet_NaN()).ok(),
          "NaN profile query accepted");
  input = Multiple(.1);
  input.geometry = [](double x) {
    return SEP::Core::Result<C::TubeGeometryPoint>::Success(
        {x, 10 + x, std::exp(x), 0, 1, 0});
  };
  Require(
      !C::IsoPolyTubeWind::Solve(input).ok(),
      "geometry without a regular sonic point acquired a fabricated transonic wind");
}
}  // namespace
int main(int argc, char** argv) {
  const std::map<std::string, std::function<void()>> tests = {
      {"IPW01", IPW01}, {"IPW02", IPW02}, {"IPW03", IPW03}, {"IPW04", IPW04},
      {"IPW05", IPW05}, {"IPW06", IPW06}, {"IPW07", IPW07}};
  if (argc == 2 && std::string(argv[1]) == "--list") {
    for (const auto& test : tests)
      std::cout << test.first << '\n';
    return 0;
  }
  const bool all = argc == 2 && std::string(argv[1]) == "--all";
  if (!all &&
      (argc != 3 || std::string(argv[1]) != "--test" || !tests.count(argv[2]))) {
    std::cerr << "--all | --list | --test IPWxx required\n";
    return 2;
  }
  std::cout << std::setprecision(17);
  bool failed = false;
  for (const auto& test : tests)
    if (all || test.first == argv[2])
      try {
        test.second();
        std::cout << '[' << test.first << "] PASS\n";
      } catch (const std::exception& e) {
        failed = true;
        std::cerr << '[' << test.first << "] FAIL " << e.what() << '\n';
      }
  return failed ? 1 : 0;
}

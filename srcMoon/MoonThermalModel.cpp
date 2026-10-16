#include "MoonThermalModel.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <utility>

namespace Moon {
namespace Thermal {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kStefanBoltzmannWm2K4 = 5.670374419e-8;

bool PositiveFinite(double value) {
  return std::isfinite(value) && value > 0.0;
}

// Return the exact integral of the configured heat-capacity polynomial from
// referenceTemperatureK to temperatureK.  The result has units J kg^-1.  The
// implicit finite-volume equation below uses the corresponding secant heat
// capacity, so its storage term is exactly the enthalpy change even when a
// time step spans a large lunar day/night temperature excursion.
double HeatCapacityEnthalpyDifference(
    const std::vector<double>& coefficients, double referenceTemperatureK,
    double temperatureK) {
  double difference = 0.0;
  for (std::size_t power = 0; power < coefficients.size(); ++power) {
    difference += coefficients[power] / static_cast<double>(power + 1) *
        (std::pow(temperatureK, power + 1) -
         std::pow(referenceTemperatureK, power + 1));
  }
  return difference;
}

// Solve a tridiagonal system in O(N).  The arrays are passed by value because
// the Thomas elimination overwrites them; the caller retains the assembled
// coefficients for diagnostics.  A singular/non-finite pivot is an explicit
// solver failure, not a reason to clip a physically invalid temperature.
bool SolveTridiagonal(std::vector<double> lower, std::vector<double> diagonal,
                      std::vector<double> upper, std::vector<double> rhs,
                      std::vector<double>* solution) {
  const std::size_t count = diagonal.size();
  if (count == 0 || lower.size() != count || upper.size() != count ||
      rhs.size() != count || solution == nullptr) return false;
  for (std::size_t row = 1; row < count; ++row) {
    if (!std::isfinite(diagonal[row - 1]) || diagonal[row - 1] == 0.0)
      return false;
    const double multiplier = lower[row] / diagonal[row - 1];
    diagonal[row] -= multiplier * upper[row - 1];
    rhs[row] -= multiplier * rhs[row - 1];
  }
  solution->assign(count, 0.0);
  if (!std::isfinite(diagonal.back()) || diagonal.back() == 0.0) return false;
  solution->back() = rhs.back() / diagonal.back();
  for (std::size_t reverse = 1; reverse < count; ++reverse) {
    const std::size_t row = count - reverse - 1;
    if (!std::isfinite(diagonal[row]) || diagonal[row] == 0.0) return false;
    (*solution)[row] =
        (rhs[row] - upper[row] * (*solution)[row + 1]) / diagonal[row];
  }
  return std::all_of(solution->begin(), solution->end(),
      [](double value) { return PositiveFinite(value); });
}

}  // namespace

bool Column::Initialize(const ColumnConfiguration& configuration,
                        double initialTemperatureK, std::string* error) {
  if (error == nullptr) return false;
  error->clear();
  if (!PositiveFinite(configuration.solarConstantWm2) ||
      !PositiveFinite(configuration.infraredEmissivity) ||
      configuration.infraredEmissivity > 1.0 ||
      !std::isfinite(configuration.normalBondAlbedo) ||
      configuration.normalBondAlbedo < 0.0 ||
      configuration.normalBondAlbedo >= 1.0 ||
      !std::isfinite(configuration.albedoCoefficientA) ||
      !std::isfinite(configuration.albedoCoefficientB) ||
      !PositiveFinite(configuration.surfaceDensityKgM3) ||
      configuration.deepDensityKgM3 <= configuration.surfaceDensityKgM3 ||
      !PositiveFinite(configuration.densityScaleHeightM) ||
      !PositiveFinite(configuration.surfaceConductivityWmK) ||
      configuration.deepConductivityWmK <
          configuration.surfaceConductivityWmK ||
      configuration.radiativeConductivityParameter < 0.0 ||
      configuration.heatCapacityPolynomial.size() != 5 ||
      configuration.geothermalFluxWm2 < 0.0 ||
      !PositiveFinite(configuration.bottomDepthM) ||
      configuration.layerCount < 2 ||
      !PositiveFinite(configuration.nonlinearToleranceK) ||
      configuration.nonlinearMaximumIterations < 1 ||
      !PositiveFinite(initialTemperatureK)) {
    *error = "invalid or incomplete lunar thermal-column configuration";
    return false;
  }

  configuration_ = configuration;
  layerThicknessM_ = configuration.bottomDepthM / configuration.layerCount;
  layerDepthM_.resize(configuration.layerCount);
  layerTemperatureK_.assign(configuration.layerCount, initialTemperatureK);
  for (int layer = 0; layer < configuration.layerCount; ++layer)
    layerDepthM_[layer] = (layer + 0.5) * layerThicknessM_;
  surfaceTemperatureK_ = initialTemperatureK;
  lastEnergyResidualJm2_ = 0.0;
  lastAbsorbedSolarFluxWm2_ = 0.0;
  lastNonlinearIterations_ = 0;
  initialized_ = true;
  return true;
}

double Column::DensityKgM3(double depthM) const {
  if (!initialized_ || !std::isfinite(depthM) || depthM < 0.0)
    return std::numeric_limits<double>::quiet_NaN();
  const double surface = configuration_.surfaceDensityKgM3;
  const double deep = configuration_.deepDensityKgM3;
  return deep - (deep - surface) *
      std::exp(-depthM / configuration_.densityScaleHeightM);
}

double Column::HeatCapacityJkgK(double temperatureK) const {
  if (!initialized_ || !PositiveFinite(temperatureK))
    return std::numeric_limits<double>::quiet_NaN();
  // Horner evaluation minimizes operations and preserves the coefficient order
  // c0 + c1*T + ... + c4*T^4 stated in the input receipt and Hayne et al.
  double result = configuration_.heatCapacityPolynomial.back();
  for (int index = 3; index >= 0; --index)
    result = result * temperatureK +
        configuration_.heatCapacityPolynomial[index];
  return PositiveFinite(result) ? result :
      std::numeric_limits<double>::quiet_NaN();
}

double Column::ConductivityWmK(double depthM, double temperatureK) const {
  const double density = DensityKgM3(depthM);
  if (!PositiveFinite(density) || !PositiveFinite(temperatureK))
    return std::numeric_limits<double>::quiet_NaN();
  const double fraction = (density - configuration_.surfaceDensityKgM3) /
      (configuration_.deepDensityKgM3 -
       configuration_.surfaceDensityKgM3);
  const double contact = configuration_.surfaceConductivityWmK + fraction *
      (configuration_.deepConductivityWmK -
       configuration_.surfaceConductivityWmK);
  return contact * (1.0 + configuration_.radiativeConductivityParameter *
      std::pow(temperatureK / 350.0, 3));
}

double Column::DirectionalBondAlbedo(double incidenceAngleRad) const {
  if (!initialized_ || !std::isfinite(incidenceAngleRad) ||
      incidenceAngleRad < 0.0 || incidenceAngleRad > 0.5 * kPi)
    return std::numeric_limits<double>::quiet_NaN();
  const double angleDeg = incidenceAngleRad * 180.0 / kPi;
  const double albedo = configuration_.normalBondAlbedo +
      configuration_.albedoCoefficientA * std::pow(angleDeg / 45.0, 3) +
      configuration_.albedoCoefficientB * std::pow(angleDeg / 90.0, 8);
  // The published directional law is only physical while A<1.  A clipped
  // value would hide an unsupported incidence/parameter combination.
  return (std::isfinite(albedo) && albedo >= 0.0 && albedo < 1.0) ? albedo :
      std::numeric_limits<double>::quiet_NaN();
}

bool Column::Advance(double timeStepS, double solarIncidenceCosine,
                     double visibleSolarDiscFraction,
                     double heliocentricDistanceAu,
                     std::string* error) {
  if (error == nullptr) return false;
  error->clear();
  if (!initialized_ || !PositiveFinite(timeStepS) ||
      !std::isfinite(solarIncidenceCosine) ||
      solarIncidenceCosine < -1.0 || solarIncidenceCosine > 1.0 ||
      !std::isfinite(visibleSolarDiscFraction) ||
      visibleSolarDiscFraction < 0.0 || visibleSolarDiscFraction > 1.0 ||
      !PositiveFinite(heliocentricDistanceAu)) {
    *error = "invalid thermal-column step input";
    return false;
  }

  double absorbed = 0.0;
  if (visibleSolarDiscFraction > 0.0 && solarIncidenceCosine > 0.0) {
    const double incidence = std::acos(solarIncidenceCosine);
    const double albedo = DirectionalBondAlbedo(incidence);
    if (!std::isfinite(albedo)) {
      *error = "directional Bond albedo is outside its physical domain";
      return false;
    }
    absorbed = configuration_.solarConstantWm2 /
        (heliocentricDistanceAu * heliocentricDistanceAu) *
        (1.0 - albedo) * solarIncidenceCosine * visibleSolarDiscFraction;
  }
  lastAbsorbedSolarFluxWm2_ = absorbed;

  const std::vector<double> oldTemperature = layerTemperatureK_;
  std::vector<double> guess(layerTemperatureK_.size() + 1);
  guess[0] = surfaceTemperatureK_;
  std::copy(layerTemperatureK_.begin(), layerTemperatureK_.end(),
            guess.begin() + 1);

  bool converged = false;
  for (int iteration = 1;
       iteration <= configuration_.nonlinearMaximumIterations; ++iteration) {
    const std::size_t count = guess.size();
    std::vector<double> lower(count, 0.0), diagonal(count, 0.0),
        upper(count, 0.0), rhs(count, 0.0), solution;

    std::vector<double> capacity(configuration_.layerCount);
    std::vector<double> conductivity(configuration_.layerCount);
    for (int layer = 0; layer < configuration_.layerCount; ++layer) {
      const double temperature = guess[layer + 1];
      const double oldLayerTemperature = oldTemperature[layer];
      const double temperatureChange = temperature - oldLayerTemperature;
      // Backward Euler written with cp(T_new) is not energy conservative for
      // the strongly temperature-dependent lunar-regolith heat capacity.  A
      // secant of the exact polynomial enthalpy makes
      // C_secant*(T_new-T_old) == h(T_new)-h(T_old).  At an indistinguishable
      // temperature change the analytic derivative cp(T) avoids cancellation.
      const double secantHeatCapacity =
          std::abs(temperatureChange) >
                  16.0 * std::numeric_limits<double>::epsilon() *
                      std::max(std::abs(temperature),
                               std::abs(oldLayerTemperature))
              ? HeatCapacityEnthalpyDifference(
                    configuration_.heatCapacityPolynomial,
                    oldLayerTemperature, temperature) /
                    temperatureChange
              : HeatCapacityJkgK(temperature);
      capacity[layer] = DensityKgM3(layerDepthM_[layer]) *
          secantHeatCapacity * layerThicknessM_;
      conductivity[layer] = ConductivityWmK(layerDepthM_[layer], temperature);
      if (!PositiveFinite(capacity[layer]) ||
          !PositiveFinite(conductivity[layer])) {
        *error = "non-physical thermal material property during iteration";
        return false;
      }
    }

    // The surface-to-first-cell distance is dz/2.  Linearize epsilon*sigma*T^4
    // about the current surface guess.  This preserves the exact nonlinear
    // balance on convergence without assigning heat capacity to the boundary.
    const double topConductance = conductivity[0] / (0.5 * layerThicknessM_);
    const double surfaceGuess = guess[0];
    const double radiationDerivative = 4.0 *
        configuration_.infraredEmissivity * kStefanBoltzmannWm2K4 *
        std::pow(surfaceGuess, 3);
    diagonal[0] = radiationDerivative + topConductance;
    upper[0] = -topConductance;
    rhs[0] = absorbed + 3.0 * configuration_.infraredEmissivity *
        kStefanBoltzmannWm2K4 * std::pow(surfaceGuess, 4);

    for (int layer = 0; layer < configuration_.layerCount; ++layer) {
      const std::size_t row = static_cast<std::size_t>(layer + 1);
      const double storage = capacity[layer] / timeStepS;
      diagonal[row] = storage;
      rhs[row] = storage * oldTemperature[layer];
      if (layer == 0) {
        lower[row] = -topConductance;
        diagonal[row] += topConductance;
      } else {
        const double conductance = 2.0 * conductivity[layer - 1] *
            conductivity[layer] /
            ((conductivity[layer - 1] + conductivity[layer]) *
             layerThicknessM_);
        lower[row] = -conductance;
        diagonal[row] += conductance;
      }
      if (layer + 1 < configuration_.layerCount) {
        const double conductance = 2.0 * conductivity[layer] *
            conductivity[layer + 1] /
            ((conductivity[layer] + conductivity[layer + 1]) *
             layerThicknessM_);
        upper[row] = -conductance;
        diagonal[row] += conductance;
      } else {
        // A positive lunar geothermal flux enters the deepest control volume
        // from below; hence it appears as a positive source on the RHS.
        rhs[row] += configuration_.geothermalFluxWm2;
      }
    }

    if (!SolveTridiagonal(lower, diagonal, upper, rhs, &solution)) {
      *error = "thermal-column tridiagonal solve failed";
      return false;
    }
    double maximumChange = 0.0;
    for (std::size_t index = 0; index < guess.size(); ++index)
      maximumChange = std::max(maximumChange,
                               std::abs(solution[index] - guess[index]));
    guess = std::move(solution);
    lastNonlinearIterations_ = iteration;
    if (maximumChange <= configuration_.nonlinearToleranceK) {
      converged = true;
      break;
    }
  }
  if (!converged) {
    std::ostringstream message;
    message << "thermal-column nonlinear iteration did not converge within "
            << configuration_.nonlinearMaximumIterations
            << " iterations; reduce thermal_time_step_s, relax only after "
               "a documented convergence study, or increase the configured "
               "iteration limit";
    *error = message.str();
    return false;
  }

  surfaceTemperatureK_ = guess[0];
  std::copy(guess.begin() + 1, guess.end(), layerTemperatureK_.begin());

  // Independently recompute the column-integrated energy receipt.  Heat
  // capacity is temperature-dependent, so integrate its polynomial exactly
  // between old and new temperatures rather than using cp at one endpoint.
  double storedEnergy = 0.0;
  for (int layer = 0; layer < configuration_.layerCount; ++layer) {
    const double oldT = oldTemperature[layer];
    const double newT = layerTemperatureK_[layer];
    const double enthalpyDifference = HeatCapacityEnthalpyDifference(
        configuration_.heatCapacityPolynomial, oldT, newT);
    storedEnergy += DensityKgM3(layerDepthM_[layer]) * layerThicknessM_ *
        enthalpyDifference;
  }
  const double emitted = configuration_.infraredEmissivity *
      kStefanBoltzmannWm2K4 * std::pow(surfaceTemperatureK_, 4);
  const double boundaryEnergy = timeStepS *
      (absorbed + configuration_.geothermalFluxWm2 - emitted);
  lastEnergyResidualJm2_ = storedEnergy - boundaryEnergy;
  return true;
}

double Column::TotalSensibleEnthalpyJm2() const {
  if (!initialized_) return std::numeric_limits<double>::quiet_NaN();
  double result=0.0;
  for (int layer=0;layer<configuration_.layerCount;++layer) {
    result+=DensityKgM3(layerDepthM_[layer])*layerThicknessM_*
        HeatCapacityEnthalpyDifference(configuration_.heatCapacityPolynomial,
                                       0.0,layerTemperatureK_[layer]);
  }
  return result;
}

}  // namespace Thermal
}  // namespace Moon

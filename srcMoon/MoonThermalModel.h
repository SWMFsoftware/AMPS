#ifndef SRCMOON_MOON_THERMAL_MODEL_H
#define SRCMOON_MOON_THERMAL_MODEL_H

#include <string>
#include <vector>

namespace Moon {
namespace Thermal {

// Physical and numerical receipt for one vertical regolith column.  Every
// dimensional field uses SI units.  The production input parser, rather than
// this structure, owns policy about which values are mandatory.
struct ColumnConfiguration {
  double solarConstantWm2 = 0.0;
  double infraredEmissivity = 0.0;
  double normalBondAlbedo = 0.0;
  double albedoCoefficientA = 0.0;
  double albedoCoefficientB = 0.0;
  double surfaceDensityKgM3 = 0.0;
  double deepDensityKgM3 = 0.0;
  double densityScaleHeightM = 0.0;
  double surfaceConductivityWmK = 0.0;
  double deepConductivityWmK = 0.0;
  double radiativeConductivityParameter = 0.0;
  std::vector<double> heatCapacityPolynomial;
  double geothermalFluxWm2 = 0.0;
  double bottomDepthM = 0.0;
  int layerCount = 0;
  double nonlinearToleranceK = 0.0;
  int nonlinearMaximumIterations = 0;
};

// Depth-resolved, conservative finite-volume solution of the one-dimensional
// lunar regolith heat equation.  Layer temperatures are cell averages; the
// separately solved surface temperature has no artificial heat capacity.
// Positive conductive flux is downward, whereas geothermalFluxWm2 is defined
// positive upward at the bottom boundary.
class Column {
 public:
  bool Initialize(const ColumnConfiguration& configuration,
                  double initialTemperatureK, std::string* error);

  // Advance one backward-Euler step. ``visibleSolarDiscFraction`` is the
  // unocculted fraction of the uniformly bright projected solar disc in
  // [0,1]. It is zero for terrain shadow or totality, between zero and one in
  // penumbra, and one in unobstructed sunlight. ``solarIncidenceCosine`` is
  // the cosine with the outward terrain normal. ``heliocentricDistanceAu``
  // supplies the inverse-square solar-flux correction.  A scalar fraction is
  // valid here because the approved first I17 model assumes uniform solar
  // radiance and does not model limb darkening or atmospheric transmission.
  bool Advance(double timeStepS, double solarIncidenceCosine,
               double visibleSolarDiscFraction,
               double heliocentricDistanceAu,
               std::string* error);

  double SurfaceTemperatureK() const { return surfaceTemperatureK_; }
  const std::vector<double>& LayerTemperatureK() const {
    return layerTemperatureK_;
  }
  double LastEnergyResidualJm2() const { return lastEnergyResidualJm2_; }
  double LastAbsorbedSolarFluxWm2() const {
    return lastAbsorbedSolarFluxWm2_;
  }
  int LastNonlinearIterations() const { return lastNonlinearIterations_; }

  // Whole-column sensible enthalpy relative to the exact polynomial's zero-K
  // reference [J m^-2].  Only differences are physically used.  Exposing the
  // conservative state measure lets periodic spin-up test the deep profile,
  // rather than declaring convergence from surface temperature alone.
  double TotalSensibleEnthalpyJm2() const;

  // Public material kernels are used by U11's independent limiting-case
  // checks.  Invalid input returns NaN and never supplies a plausible default.
  double DensityKgM3(double depthM) const;
  double HeatCapacityJkgK(double temperatureK) const;
  double ConductivityWmK(double depthM, double temperatureK) const;
  double DirectionalBondAlbedo(double incidenceAngleRad) const;

 private:
  ColumnConfiguration configuration_;
  std::vector<double> layerDepthM_;
  std::vector<double> layerTemperatureK_;
  double layerThicknessM_ = 0.0;
  double surfaceTemperatureK_ = 0.0;
  double lastEnergyResidualJm2_ = 0.0;
  double lastAbsorbedSolarFluxWm2_ = 0.0;
  int lastNonlinearIterations_ = 0;
  bool initialized_ = false;
};

}  // namespace Thermal
}  // namespace Moon

#endif  // SRCMOON_MOON_THERMAL_MODEL_H

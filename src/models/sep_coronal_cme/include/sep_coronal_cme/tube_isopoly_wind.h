#ifndef SEP_CORONAL_CME_TUBE_ISOPOLY_WIND_H
#define SEP_CORONAL_CME_TUBE_ISOPOLY_WIND_H

#include "sep_coronal_cme/flux_tube_wind.h"

#include <cstddef>
#include <functional>
#include <vector>

namespace SEP::CoronalCME {

struct IsoPolySpecies {
  // Constant composition: species number density = particlesPerKg * rho.
  // The outer law is T_ref * (rho/rho_ref)^(polytropicIndex - 1); at
  // the stated reference mass density it agrees with the inner temperature.
  double particlesPerKg=0, referenceTemperatureK=0, referenceDensityKgM3=0;
  double polytropicIndex=1, blendBeginM=0, blendEndM=0, adiabaticIndex=0;
};

struct IsoPolyTubeInput {
  // Coordinate is outward distance [m]. A, Phi, and their coordinate
  // derivatives [1/m, m/s^2] must come from the same smooth geometry. The
  // callback may fail (e.g. closed field, magnetic null, missing coverage).
  // It owns its captures and remains immutable for the resulting profile.
  std::function<Core::Result<TubeGeometryPoint>(double)> geometry;
  std::vector<IsoPolySpecies> species;
  double beginM=0, endM=0, massFluxKgPerS=0;
  unsigned criticalScanIntervals=256, maximumIntegrationAttempts=65536;
  double relativeTolerance=1e-9;
};

struct IsoPolyTubeSample {
  WindPoint wind;
  double coordinateM=0, speedDerivativePerS=0, densityDerivativeKgM4=0;
  double pressureDerivativePaPerM=0;
  std::vector<double> temperatureK, partialPressurePa, netHeatWPerM3;
};

class IsoPolyTubeWind {
 public:
  // Thermal-only steady flow. The unique regular sonic branch must span the
  // declared interval from strictly subsonic to strictly supersonic flow.
  static Core::Result<IsoPolyTubeWind> Solve(const IsoPolyTubeInput& input);
  Core::Result<IsoPolyTubeSample> Evaluate(double coordinateM) const;
  double CriticalCoordinateM() const { return criticalCoordinateM_; }

 private:
  struct Node { double x=0, logSpeed=0, slope=0; };
  IsoPolyTubeInput input_;
  double coordinateScaleM_=0, coordinateExtent_=0, criticalCoordinateM_=0;
  std::vector<Node> nodes_;
};

}  // namespace SEP::CoronalCME

#endif

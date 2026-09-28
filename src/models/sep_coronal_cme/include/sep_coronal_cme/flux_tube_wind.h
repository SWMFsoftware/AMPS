#ifndef SEP_CORONAL_CME_FLUX_TUBE_WIND_H
#define SEP_CORONAL_CME_FLUX_TUBE_WIND_H

#include "sep_status.h"

#include <functional>
#include <vector>

namespace SEP { namespace CoronalCME {

struct WindPoint {
  double radiusM=0.0, areaM2=0.0, speedMPerS=0.0;
  double densityKgM3=0.0, pressurePa=0.0, soundSpeedSquaredM2S2=0.0;
  double effectivePotentialM2S2=0.0;
};

struct CriticalCandidate {
  double coordinateM=0.0, soundSpeedSquaredM2S2=0.0;
  bool globallyAdmissible=false;
};

struct TubeGeometryPoint {
  double coordinateM=0.0;
  double radiusM=0.0;
  double areaM2=0.0;
  double effectivePotentialM2S2=0.0;
  double dLogAreaDsPerM=0.0;
  double dPotentialDsMPerS2=0.0;
};

struct PolytropicWindSolution {
  std::vector<WindPoint> points;
  std::vector<CriticalCandidate> criticalCandidates;
  std::size_t selectedCriticalIndex=0;
  double massFluxKgPerS=0.0;
  double bernoulliM2S2=0.0;
  double polytropicConstantSI=0.0;
};

double FluxTubeArea(double magneticFluxWb,double magneticFieldMagnitudeT);
double EffectivePotential(double radiusM,double cylindricalRadiusM,
                          double rotationRateRadPerS);

// Exact radial isothermal Parker relation, solved on the subsonic branch below
// r_c and supersonic branch above it. density follows one independent mass flux.
Core::Result<std::vector<WindPoint>> SolveRadialIsothermalParker(
    const std::vector<double>& radiiM,double soundSpeedMPerS,double criticalRadiusM,
    double massFluxKgPerS,double referenceAreaAtOneM2=1.0);

Core::Status CheckWindInvariants(const std::vector<WindPoint>& points,
    double gammaWind,double relativeTolerance,double* maximumMassResidual=nullptr,
    double* maximumBernoulliResidual=nullptr,double* maximumEntropyResidual=nullptr);

Core::Result<std::vector<CriticalCandidate>> FindCriticalCandidates(
    const std::vector<double>& coordinateM,
    const std::vector<double>& dLogAreaDsPerM,
    const std::vector<double>& dPotentialDsMPerS2);

// Solves the algebraic steady nozzle invariants on both sides of a declared
// regular critical point. The caller supplies completed smooth tube geometry;
// Stage 3 is responsible for tracing and for splitting physical interfaces.
Core::Result<PolytropicWindSolution> SolvePolytropicTube(
    const std::vector<TubeGeometryPoint>& geometry,double gammaWind,
    double polytropicConstantSI,double criticalDensityKgM3,
    std::size_t selectedCriticalIndex,double criticalRelativeTolerance=1e-8);

Core::Result<double> InvertTargetSpeed(
    const std::function<Core::Result<double>(double)>& speedFromTemperature,
    double targetSpeedMPerS,double minimumTemperatureK,double maximumTemperatureK,
    double relativeTolerance=1.0e-10);

Core::Result<double> ResolveMassLoading(double densityKgM3,double speedMPerS,
    double magneticFieldT);
Core::Result<double> ResolveMassLoadingFromRadialFlux(
    double radialMassFluxKgM2S,double radialMagneticFieldT);

} }  // namespace SEP::CoronalCME
#endif

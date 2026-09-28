#ifndef SEP_CORONAL_CME_TURBULENCE_TRANSPORT_H
#define SEP_CORONAL_CME_TURBULENCE_TRANSPORT_H

#include "sep_coronal_cme/source_surface_coupling.h"
#include "sep_status.h"

#include <functional>

namespace SEP { namespace CoronalCME {

// Outward/inward are physical radial labels and are authoritative.  The
// parallel/antiparallel labels are derived only after the categorical magnetic
// sector is known, which makes an ideal HCS reversal a pure label exchange.
struct DirectionalWaveState {
  bool applicable = true;
  double outwardJPerM3 = 0.0;
  double inwardJPerM3 = 0.0;
  double parallelJPerM3 = 0.0;
  double antiparallelJPerM3 = 0.0;
  double totalJPerM3 = 0.0;
  double outwardCrossHelicity = 0.0;
  double fieldCrossHelicity = 0.0;
  double deltaBSquaredT2 = 0.0;
};
Core::Result<DirectionalWaveState> BuildDirectionalWaveState(
    double outwardJPerM3, double inwardJPerM3, MagneticSector sector);
DirectionalWaveState InapplicableWaveState();

struct PrescribedWaveParameters {
  double referenceRadiusM = 0.0;
  double referenceEnergyJPerM3 = 0.0;
  double energyRadialExponent = 0.0;
  double outwardCrossHelicity = 1.0;
  double referenceCorrelationLengthM = 0.0;
  double correlationLengthExponent = 0.0;
};
struct PrescribedWaveEvaluation {
  DirectionalWaveState wave;
  double correlationLengthM = 0.0;
};
Core::Result<PrescribedWaveEvaluation> EvaluatePrescribedWave(
    const PrescribedWaveParameters& parameters, double radiusM,
    MagneticSector sector);

struct WkbReferenceState {
  double areaM2 = 0.0;
  double fieldAlignedSpeedMPerS = 0.0;
  double alfvenSpeedMPerS = 0.0;
  double outwardEnergyJPerM3 = 0.0;
};
Core::Result<DirectionalWaveState> EvaluateWkbOutwardWave(
    const WkbReferenceState& reference, double areaM2,
    double fieldAlignedSpeedMPerS, double alfvenSpeedMPerS,
    MagneticSector sector);
Core::Result<double> WkbWaveAction(const DirectionalWaveState& wave,
                                   double areaM2,
                                   double fieldAlignedSpeedMPerS,
                                   double alfvenSpeedMPerS);

// Symmetric closed-loop baseline.  Two equal footpoint sources decay over a
// correlation length and are normalized so the total equals the prescribed
// footpoint value at both ends; it is deliberately distinct from outward WKB.
Core::Result<DirectionalWaveState> EvaluateClosedLoopWave(
    double distanceFromFirstFootpointM, double loopLengthM,
    double footpointTotalEnergyJPerM3, double decayLengthM,
    MagneticSector sector);

struct WaveValidityInput {
  double magneticMagnitudeT = 0.0;
  double thermalPressurePa = 0.0;
  double massDensityKgM3 = 0.0;
  double fieldAlignedSpeedMPerS = 0.0;
  double wavePressureGradientPaPerM = 0.0;
  double inertialAccelerationMPerS2 = 0.0;
  double pressureAccelerationMPerS2 = 0.0;
  double potentialAccelerationMPerS2 = 0.0;
  double absoluteWaveAccelerationUncertaintyMPerS2 = 0.0;
  double maximumAbsoluteUncertaintyMPerS2 = 0.0;
  double maximumForceFraction = 0.0;
  double maximumDeltaBOverB = 0.0;
  bool enforceSmallAmplitude = true;
};
struct WaveValidityDiagnostics {
  double deltaBOverB = 0.0;
  double waveToThermalPressure = 0.0;
  double waveToRamPlusThermalPressure = 0.0;
  double waveToMagneticPressure = 0.0;
  double signedWaveAccelerationMPerS2 = 0.0;
  double retainedAccelerationMPerS2 = 0.0;
  double forceExcessFraction = 0.0;
  bool passed = false;
};
Core::Result<WaveValidityDiagnostics> EvaluateWaveValidity(
    const DirectionalWaveState& wave, const WaveValidityInput& input);

Core::Result<double> SinglePowerLawMeanFreePath(
    double referenceMeanFreePathM, double radiusM, double referenceRadiusM,
    double rigidityV, double referenceRigidityV, double radialExponent,
    double rigidityExponent);
Core::Result<double> SmoothBrokenMeanFreePath(
    double breakMeanFreePathM, double radiusM, double breakRadiusM,
    double rigidityV, double referenceRigidityV, double innerRadialExponent,
    double outerRadialExponent, double smoothness,
    double rigidityExponent);

// Implements lambda_parallel=(3v/8) integral[(1-mu^2)^2/D_mumu]dmu.
// Endpoints are excluded from the midpoint quadrature, so a physical
// D_mumu proportional to (1-mu^2) is integrable without 0/0 evaluation.
Core::Result<double> IntegrateParallelMeanFreePath(
    double speedMPerS, int intervals,
    const std::function<double(double)>& pitchAngleDiffusionPerS);

// One Euler--Maruyama step for D_mumu=D0(1-mu^2), including the Ito drift
// dD/dmu and reflecting endpoints.  The supplied normal deviate keeps this
// dependency-light kernel deterministic and independently testable.
Core::Result<double> StepIsotropicPitchAngleDiffusion(
    double mu, double diffusionScalePerS, double timeStepS,
    double normalDeviate);

struct DiscreteScatterResult {
  double mu = 0.0;
  bool scattered = false;
  double scatterProbability = 0.0;
};
Core::Result<DiscreteScatterResult> StepIsotropicPoissonScattering(
    double mu, double speedMPerS, double meanFreePathM, double timeStepS,
    double eventUniform01, double directionUniform01);

enum class MissingCoefficientPolicy { Fail, Ballistic };
enum class CollisionModel { None, PitchAngleDiffusion,
                            DiscreteIsotropicPoisson };
Core::Status ValidateCollisionModel(CollisionModel model,
                                    bool finiteMeanFreePathAvailable,
                                    MissingCoefficientPolicy missingPolicy);

} }  // namespace SEP::CoronalCME

#endif  // SEP_CORONAL_CME_TURBULENCE_TRANSPORT_H

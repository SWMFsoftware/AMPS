#include "sep_coronal_cme/turbulence_transport.h"

#include "sep_coronal_cme/constants.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace SEP { namespace CoronalCME { namespace {

bool Finite(double value) { return std::isfinite(value); }

double SectorSign(MagneticSector sector) {
  return sector == MagneticSector::Positive ? 1.0 : -1.0;
}

}  // namespace

Core::Result<DirectionalWaveState> BuildDirectionalWaveState(
    double outwardJPerM3, double inwardJPerM3, MagneticSector sector) {
  if (!(Finite(outwardJPerM3) && Finite(inwardJPerM3) &&
        outwardJPerM3 >= 0.0 && inwardJPerM3 >= 0.0)) {
    return Core::Result<DirectionalWaveState>::Failure(
        Core::StatusCode::InvalidState,
        "directional wave energies must be finite and nonnegative");
  }
  const double total = outwardJPerM3 + inwardJPerM3;
  DirectionalWaveState result;
  result.outwardJPerM3 = outwardJPerM3;
  result.inwardJPerM3 = inwardJPerM3;
  result.totalJPerM3 = total;
  if (sector == MagneticSector::Positive) {
    result.parallelJPerM3 = outwardJPerM3;
    result.antiparallelJPerM3 = inwardJPerM3;
  } else {
    result.parallelJPerM3 = inwardJPerM3;
    result.antiparallelJPerM3 = outwardJPerM3;
  }
  result.outwardCrossHelicity = total > 0.0
      ? (outwardJPerM3 - inwardJPerM3) / total : 0.0;
  result.fieldCrossHelicity = SectorSign(sector) *
      result.outwardCrossHelicity;
  // The model uses total Alfvén-wave energy, so ideal equipartition gives
  // delta(B)^2=mu0*(w_out+w_in), not twice the energy of one direction.
  result.deltaBSquaredT2 = Constants::kVacuumPermeabilityHPerM * total;
  return Core::Result<DirectionalWaveState>::Success(result);
}

DirectionalWaveState InapplicableWaveState() {
  DirectionalWaveState result;
  result.applicable = false;
  return result;
}

Core::Result<PrescribedWaveEvaluation> EvaluatePrescribedWave(
    const PrescribedWaveParameters& parameters, double radiusM,
    MagneticSector sector) {
  if (!(parameters.referenceRadiusM > 0.0 &&
        parameters.referenceEnergyJPerM3 > 0.0 &&
        parameters.referenceCorrelationLengthM > 0.0 && radiusM > 0.0 &&
        Finite(parameters.energyRadialExponent) &&
        Finite(parameters.correlationLengthExponent) &&
        parameters.outwardCrossHelicity >= -1.0 &&
        parameters.outwardCrossHelicity <= 1.0)) {
    return Core::Result<PrescribedWaveEvaluation>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "prescribed waves require positive references and |sigma_c|<=1");
  }
  const double radialRatio = radiusM / parameters.referenceRadiusM;
  const double total = parameters.referenceEnergyJPerM3 *
      std::pow(radialRatio, -parameters.energyRadialExponent);
  const double outward = 0.5 *
      (1.0 + parameters.outwardCrossHelicity) * total;
  const double inward = total - outward;
  auto wave = BuildDirectionalWaveState(outward, inward, sector);
  if (!wave.ok()) {
    return Core::Result<PrescribedWaveEvaluation>::Failure(
        wave.status.code, wave.status.message);
  }
  PrescribedWaveEvaluation result;
  result.wave = wave.value;
  result.correlationLengthM = parameters.referenceCorrelationLengthM *
      std::pow(radialRatio, parameters.correlationLengthExponent);
  if (!(Finite(result.correlationLengthM) && result.correlationLengthM > 0.0)) {
    return Core::Result<PrescribedWaveEvaluation>::Failure(
        Core::StatusCode::InvalidState,
        "prescribed correlation length is nonpositive or nonfinite");
  }
  return Core::Result<PrescribedWaveEvaluation>::Success(result);
}

Core::Result<DirectionalWaveState> EvaluateWkbOutwardWave(
    const WkbReferenceState& reference, double areaM2,
    double fieldAlignedSpeedMPerS, double alfvenSpeedMPerS,
    MagneticSector sector) {
  if (!(reference.areaM2 > 0.0 &&
        reference.fieldAlignedSpeedMPerS >= 0.0 &&
        reference.alfvenSpeedMPerS > 0.0 &&
        reference.outwardEnergyJPerM3 > 0.0 && areaM2 > 0.0 &&
        fieldAlignedSpeedMPerS >= 0.0 && alfvenSpeedMPerS > 0.0)) {
    return Core::Result<DirectionalWaveState>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "WKB wave action requires positive area, v_A, and reference energy");
  }
  const double referenceGroupSpeed = reference.fieldAlignedSpeedMPerS +
      reference.alfvenSpeedMPerS;
  const double groupSpeed = fieldAlignedSpeedMPerS + alfvenSpeedMPerS;
  const double outward = reference.outwardEnergyJPerM3 *
      reference.areaM2 / areaM2 *
      alfvenSpeedMPerS / reference.alfvenSpeedMPerS *
      std::pow(referenceGroupSpeed / groupSpeed, 2);
  return BuildDirectionalWaveState(outward, 0.0, sector);
}

Core::Result<double> WkbWaveAction(
    const DirectionalWaveState& wave, double areaM2,
    double fieldAlignedSpeedMPerS, double alfvenSpeedMPerS) {
  if (!wave.applicable || wave.outwardJPerM3 < 0.0 || areaM2 <= 0.0 ||
      fieldAlignedSpeedMPerS < 0.0 || alfvenSpeedMPerS <= 0.0) {
    return Core::Result<double>::Failure(Core::StatusCode::InvalidState,
        "wave action requested for an inapplicable/nonphysical state");
  }
  return Core::Result<double>::Success(areaM2 * wave.outwardJPerM3 *
      std::pow(fieldAlignedSpeedMPerS + alfvenSpeedMPerS, 2) /
      alfvenSpeedMPerS);
}

Core::Result<DirectionalWaveState> EvaluateClosedLoopWave(
    double distanceFromFirstFootpointM, double loopLengthM,
    double footpointTotalEnergyJPerM3, double decayLengthM,
    MagneticSector sector) {
  if (!(loopLengthM > 0.0 && distanceFromFirstFootpointM >= 0.0 &&
        distanceFromFirstFootpointM <= loopLengthM &&
        footpointTotalEnergyJPerM3 > 0.0 && decayLengthM > 0.0)) {
    return Core::Result<DirectionalWaveState>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "closed-loop wave requires 0<=s<=L and positive energy/decay length");
  }
  const double endCoupling = std::exp(-loopLengthM / decayLengthM);
  const double scale = footpointTotalEnergyJPerM3 / (1.0 + endCoupling);
  const double fromFirst = scale *
      std::exp(-distanceFromFirstFootpointM / decayLengthM);
  const double fromSecond = scale *
      std::exp(-(loopLengthM - distanceFromFirstFootpointM) / decayLengthM);
  // On a closed loop there is no globally unique outward direction.  The two
  // footpoint directions are stored symmetrically as the two local channels.
  return BuildDirectionalWaveState(fromFirst, fromSecond, sector);
}

Core::Result<WaveValidityDiagnostics> EvaluateWaveValidity(
    const DirectionalWaveState& wave, const WaveValidityInput& input) {
  if (!wave.applicable) {
    return Core::Result<WaveValidityDiagnostics>::Failure(
        Core::StatusCode::OutOfDomain,
        "wave validity is inapplicable for a direct-mean-free-path state");
  }
  if (!(input.magneticMagnitudeT > 0.0 && input.thermalPressurePa > 0.0 &&
        input.massDensityKgM3 > 0.0 &&
        input.absoluteWaveAccelerationUncertaintyMPerS2 >= 0.0 &&
        input.maximumAbsoluteUncertaintyMPerS2 >= 0.0 &&
        input.maximumForceFraction >= 0.0 &&
        input.maximumDeltaBOverB > 0.0)) {
    return Core::Result<WaveValidityDiagnostics>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "wave validity requires positive B, p, rho and nonnegative bounds");
  }
  if (input.absoluteWaveAccelerationUncertaintyMPerS2 >
      input.maximumAbsoluteUncertaintyMPerS2) {
    return Core::Result<WaveValidityDiagnostics>::Failure(
        Core::StatusCode::DataIntegrityFailure,
        "derivative uncertainty exceeds its preregistered absolute cap");
  }
  WaveValidityDiagnostics result;
  const double wavePressure = 0.5 * wave.totalJPerM3;
  const double magneticPressure = input.magneticMagnitudeT *
      input.magneticMagnitudeT /
      (2.0 * Constants::kVacuumPermeabilityHPerM);
  result.deltaBOverB = std::sqrt(std::max(0.0, wave.deltaBSquaredT2)) /
      input.magneticMagnitudeT;
  result.waveToThermalPressure = wavePressure / input.thermalPressurePa;
  result.waveToRamPlusThermalPressure = wavePressure /
      (input.thermalPressurePa + input.massDensityKgM3 *
       input.fieldAlignedSpeedMPerS * input.fieldAlignedSpeedMPerS);
  result.waveToMagneticPressure = wavePressure / magneticPressure;
  result.signedWaveAccelerationMPerS2 =
      -input.wavePressureGradientPaPerM / input.massDensityKgM3;
  result.retainedAccelerationMPerS2 =
      std::abs(input.inertialAccelerationMPerS2) +
      std::abs(input.pressureAccelerationMPerS2) +
      std::abs(input.potentialAccelerationMPerS2);
  const double excess = std::max(0.0,
      std::abs(result.signedWaveAccelerationMPerS2) -
      input.absoluteWaveAccelerationUncertaintyMPerS2);
  if (result.retainedAccelerationMPerS2 > 0.0) {
    result.forceExcessFraction = excess /
        result.retainedAccelerationMPerS2;
  } else {
    result.forceExcessFraction = excess == 0.0 ? 0.0
        : std::numeric_limits<double>::infinity();
  }
  const bool forcePassed = excess <= input.maximumForceFraction *
      result.retainedAccelerationMPerS2;
  const bool amplitudePassed = !input.enforceSmallAmplitude ||
      result.deltaBOverB <= input.maximumDeltaBOverB;
  result.passed = forcePassed && amplitudePassed;
  return Core::Result<WaveValidityDiagnostics>::Success(result);
}

Core::Result<double> SinglePowerLawMeanFreePath(
    double referenceMeanFreePathM, double radiusM, double referenceRadiusM,
    double rigidityV, double referenceRigidityV, double radialExponent,
    double rigidityExponent) {
  if (!(referenceMeanFreePathM > 0.0 && radiusM > 0.0 &&
        referenceRadiusM > 0.0 && rigidityV > 0.0 &&
        referenceRigidityV > 0.0 && Finite(radialExponent) &&
        Finite(rigidityExponent))) {
    return Core::Result<double>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "single-power-law mean free path requires positive references");
  }
  const double value = referenceMeanFreePathM *
      std::pow(radiusM / referenceRadiusM, radialExponent) *
      std::pow(rigidityV / referenceRigidityV, rigidityExponent);
  if (!(Finite(value) && value > 0.0)) {
    return Core::Result<double>::Failure(Core::StatusCode::InvalidState,
        "single-power-law mean free path is nonpositive or nonfinite");
  }
  return Core::Result<double>::Success(value);
}

Core::Result<double> SmoothBrokenMeanFreePath(
    double breakMeanFreePathM, double radiusM, double breakRadiusM,
    double rigidityV, double referenceRigidityV, double innerRadialExponent,
    double outerRadialExponent, double smoothness,
    double rigidityExponent) {
  if (!(breakMeanFreePathM > 0.0 && radiusM > 0.0 && breakRadiusM > 0.0 &&
        rigidityV > 0.0 && referenceRigidityV > 0.0 && smoothness > 0.0)) {
    return Core::Result<double>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "smooth-broken mean free path requires positive scales and nu");
  }
  const double x = radiusM / breakRadiusM;
  const double transition = std::pow((1.0 + std::pow(x, smoothness)) / 2.0,
      (outerRadialExponent - innerRadialExponent) / smoothness);
  const double value = breakMeanFreePathM *
      std::pow(x, innerRadialExponent) * transition *
      std::pow(rigidityV / referenceRigidityV, rigidityExponent);
  if (!(Finite(value) && value > 0.0)) {
    return Core::Result<double>::Failure(Core::StatusCode::InvalidState,
        "smooth-broken mean free path is nonpositive or nonfinite");
  }
  return Core::Result<double>::Success(value);
}

Core::Result<double> IntegrateParallelMeanFreePath(
    double speedMPerS, int intervals,
    const std::function<double(double)>& pitchAngleDiffusionPerS) {
  if (!(speedMPerS > 0.0 && intervals >= 16 && pitchAngleDiffusionPerS)) {
    return Core::Result<double>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "mean-free-path quadrature requires speed, callback, and >=16 bins");
  }
  const double width = 2.0 / intervals;
  double integral = 0.0;
  for (int index = 0; index < intervals; ++index) {
    const double mu = -1.0 + (index + 0.5) * width;
    const double diffusion = pitchAngleDiffusionPerS(mu);
    if (!(Finite(diffusion) && diffusion > 0.0)) {
      return Core::Result<double>::Failure(Core::StatusCode::InvalidState,
          "D_mumu must be finite and positive inside (-1,1)");
    }
    integral += std::pow(1.0 - mu * mu, 2) / diffusion * width;
  }
  return Core::Result<double>::Success(3.0 * speedMPerS * integral / 8.0);
}

Core::Result<double> StepIsotropicPitchAngleDiffusion(
    double mu, double diffusionScalePerS, double timeStepS,
    double normalDeviate) {
  if (!(mu >= -1.0 && mu <= 1.0 && diffusionScalePerS >= 0.0 &&
        timeStepS >= 0.0 && Finite(normalDeviate))) {
    return Core::Result<double>::Failure(Core::StatusCode::InvalidConfiguration,
        "pitch-angle step received an invalid state/time/deviate");
  }
  const double diffusion = diffusionScalePerS * (1.0 - mu * mu);
  const double drift = -2.0 * diffusionScalePerS * mu;
  double next = mu + drift * timeStepS +
      std::sqrt(2.0 * diffusion * timeStepS) * normalDeviate;
  // Repeated mirror reflection conserves probability for a finite numerical
  // overshoot without pinning a stochastic trajectory at mu=+/-1.
  while (next > 1.0 || next < -1.0) {
    if (next > 1.0) next = 2.0 - next;
    if (next < -1.0) next = -2.0 - next;
  }
  return Core::Result<double>::Success(next);
}

Core::Result<DiscreteScatterResult> StepIsotropicPoissonScattering(
    double mu, double speedMPerS, double meanFreePathM, double timeStepS,
    double eventUniform01, double directionUniform01) {
  if (!(mu >= -1.0 && mu <= 1.0 && speedMPerS > 0.0 &&
        meanFreePathM > 0.0 && timeStepS >= 0.0 &&
        eventUniform01 >= 0.0 && eventUniform01 < 1.0 &&
        directionUniform01 >= 0.0 && directionUniform01 < 1.0)) {
    return Core::Result<DiscreteScatterResult>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "discrete scattering requires physical state and [0,1) deviates");
  }
  DiscreteScatterResult result;
  result.scatterProbability = 1.0 -
      std::exp(-speedMPerS * timeStepS / meanFreePathM);
  result.scattered = eventUniform01 < result.scatterProbability;
  result.mu = result.scattered ? 2.0 * directionUniform01 - 1.0 : mu;
  return Core::Result<DiscreteScatterResult>::Success(result);
}

Core::Status ValidateCollisionModel(
    CollisionModel model, bool finiteMeanFreePathAvailable,
    MissingCoefficientPolicy missingPolicy) {
  if (model == CollisionModel::None) return Core::Status::Success();
  if (!finiteMeanFreePathAvailable &&
      missingPolicy == MissingCoefficientPolicy::Fail) {
    return Core::Status::Failure(Core::StatusCode::OutOfDomain,
        "focused collision model has no finite mean-free-path authority");
  }
  if (!finiteMeanFreePathAvailable &&
      missingPolicy == MissingCoefficientPolicy::Ballistic) {
    return Core::Status::Success();
  }
  return Core::Status::Success();
}

} }  // namespace SEP::CoronalCME

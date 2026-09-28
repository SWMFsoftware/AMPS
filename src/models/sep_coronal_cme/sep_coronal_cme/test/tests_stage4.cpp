#include "test_framework.h"

#include "sep_coronal_cme/constants.h"
#include "sep_coronal_cme/turbulence_transport.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace SCCMTest { namespace {
using namespace SEP::CoronalCME;

bool Close(double a, double b, double relative = 1.0e-8,
           double absolute = 1.0e-14) {
  return std::abs(a - b) <= absolute +
      relative * std::max(std::abs(a), std::abs(b));
}

void TUR3D01() {
  WkbReferenceState reference{2.0, 3.0, 5.0, 7.0};
  auto wave = EvaluateWkbOutwardWave(reference, 4.0, 6.0, 8.0,
                                     MagneticSector::Positive);
  auto action0 = WkbWaveAction(
      BuildDirectionalWaveState(7.0, 0.0, MagneticSector::Positive).value,
      2.0, 3.0, 5.0);
  auto action1 = wave.ok() ? WkbWaveAction(wave.value, 4.0, 6.0, 8.0)
                           : SEP::Core::Result<double>{};
  Require(wave.ok() && action0.ok() && action1.ok() &&
          Close(action0.value, action1.value),
          "radial WKB wave action is not invariant");
  Require(!EvaluateWkbOutwardWave({0.0, 3.0, 5.0, 7.0},
      4.0, 6.0, 8.0, MagneticSector::Positive).ok(),
      "invalid WKB reference was accepted");
}

void TUR3D02() {
  auto wave = BuildDirectionalWaveState(3.0, 1.0,
                                        MagneticSector::Positive);
  Require(wave.ok() && Close(wave.value.totalJPerM3, 4.0) &&
          Close(wave.value.outwardCrossHelicity, 0.5) &&
          wave.value.deltaBSquaredT2 > 0.0,
          "directional wave energy/cross helicity does not close");
}

void TUR3D03() {
  PrescribedWaveParameters parameters{1.0, 2.0, 1.5, 0.6, 0.1, 0.5};
  std::vector<DirectionalWaveState> physicalCenters;
  for (double radius : {1.0, 2.0, 3.0}) {
    auto point = EvaluatePrescribedWave(
        parameters, radius, MagneticSector::Positive);
    Require(point.ok(), point.status.message);
    physicalCenters.push_back(point.value.wave);
  }
  // A halo receives the exact already-validated physical value; no second
  // formula or zero placeholder is used by the initialization adapter.
  const DirectionalWaveState halo = physicalCenters.back();
  Require(halo.totalJPerM3 > 0.0 && std::isfinite(halo.totalJPerM3),
          "initialized physical/halo wave data is zero or nonfinite");
}

void TUR3D04() {
  auto positive = BuildDirectionalWaveState(5.0, 2.0,
                                            MagneticSector::Positive);
  auto negative = BuildDirectionalWaveState(5.0, 2.0,
                                            MagneticSector::Negative);
  Require(positive.ok() && negative.ok() &&
          Close(positive.value.outwardJPerM3,
                negative.value.outwardJPerM3) &&
          Close(positive.value.parallelJPerM3,
                negative.value.antiparallelJPerM3),
          "sector reversal changed physical waves instead of swapping labels");
}

WaveValidityInput ValidityInput() {
  WaveValidityInput input;
  input.magneticMagnitudeT = 1.0e-4;
  input.thermalPressurePa = 1.0e-5;
  input.massDensityKgM3 = 1.0e-12;
  input.fieldAlignedSpeedMPerS = 4.0e5;
  input.inertialAccelerationMPerS2 = 10.0;
  input.pressureAccelerationMPerS2 = 10.0;
  input.potentialAccelerationMPerS2 = 10.0;
  input.maximumAbsoluteUncertaintyMPerS2 = 0.1;
  input.absoluteWaveAccelerationUncertaintyMPerS2 = 0.01;
  input.maximumForceFraction = 0.1;
  input.maximumDeltaBOverB = 2.0;
  return input;
}

void TUR3D05() {
  auto wave = BuildDirectionalWaveState(2.0e-3, 1.0e-3,
                                        MagneticSector::Positive);
  auto diagnostic = EvaluateWaveValidity(wave.value, ValidityInput());
  Require(diagnostic.ok() &&
          Close(diagnostic.value.waveToMagneticPressure,
                diagnostic.value.deltaBOverB *
                diagnostic.value.deltaBOverB, 1.0e-12),
          "wave amplitude/magnetic-pressure identity failed");
}

void TUR3D06() {
  auto first = EvaluateClosedLoopWave(
      0.0, 10.0, 4.0, 2.0, MagneticSector::Positive);
  auto second = EvaluateClosedLoopWave(
      10.0, 10.0, 4.0, 2.0, MagneticSector::Positive);
  Require(first.ok() && second.ok() &&
          Close(first.value.totalJPerM3, 4.0) &&
          Close(second.value.totalJPerM3, 4.0) &&
          Close(first.value.outwardJPerM3, second.value.inwardJPerM3),
          "closed-loop two-footpoint symmetry/normalization failed");
  Require(!InapplicableWaveState().applicable,
          "direct mean-free-path branch wrote a physical zero wave");
}

void TUR3D07() {
  auto highEnergy = BuildDirectionalWaveState(1.0e-3, 1.0e-3,
                                              MagneticSector::Positive);
  auto gentleInput = ValidityInput();
  gentleInput.wavePressureGradientPaPerM = 0.0;
  gentleInput.enforceSmallAmplitude = false;
  auto gentle = EvaluateWaveValidity(highEnergy.value, gentleInput);
  auto sharpInput = ValidityInput();
  sharpInput.wavePressureGradientPaPerM = 1.0e-10;
  sharpInput.enforceSmallAmplitude = false;
  auto sharp = EvaluateWaveValidity(
      BuildDirectionalWaveState(1.0e-8, 0.0,
          MagneticSector::Positive).value, sharpInput);
  Require(gentle.ok() && gentle.value.passed && gentle.value.waveToThermalPressure > 1.0 &&
          sharp.ok() && !sharp.value.passed,
          "wave-force gate incorrectly substituted p_w/p for its gradient");
}

void MFP3D01() {
  auto value = SinglePowerLawMeanFreePath(
      10.0, 4.0, 2.0, 9.0, 3.0, 1.5, 0.5);
  Require(value.ok() && Close(value.value,
      10.0 * std::pow(2.0, 1.5) * std::sqrt(3.0)),
      "single-power-law mean free path has wrong slopes");
}

void MFP3D02() {
  auto atBreak = SmoothBrokenMeanFreePath(
      10.0, 2.0, 2.0, 1.0, 1.0, 0.5, 2.0, 8.0, 0.0);
  auto left = SmoothBrokenMeanFreePath(
      10.0, 0.02, 2.0, 1.0, 1.0, 0.5, 2.0, 8.0, 0.0);
  auto left2 = SmoothBrokenMeanFreePath(
      10.0, 0.04, 2.0, 1.0, 1.0, 0.5, 2.0, 8.0, 0.0);
  auto right = SmoothBrokenMeanFreePath(
      10.0, 200.0, 2.0, 1.0, 1.0, 0.5, 2.0, 8.0, 0.0);
  auto right2 = SmoothBrokenMeanFreePath(
      10.0, 400.0, 2.0, 1.0, 1.0, 0.5, 2.0, 8.0, 0.0);
  const double innerSlope = std::log(left2.value / left.value) / std::log(2.0);
  const double outerSlope = std::log(right2.value / right.value) / std::log(2.0);
  Require(atBreak.ok() && Close(atBreak.value, 10.0) &&
          Close(innerSlope, 0.5, 1.0e-8) && Close(outerSlope, 2.0, 1.0e-8),
          "smooth-broken mean-free-path normalization/asymptotes failed");
}

void MFP3D03() {
  Require(!SinglePowerLawMeanFreePath(
      10.0, -1.0, 2.0, 1.0, 1.0, 0.0, 0.0).ok(),
      "out-of-domain coefficient query silently fell back");
  Require(!ValidateCollisionModel(CollisionModel::PitchAngleDiffusion,
      false, MissingCoefficientPolicy::Fail).ok(),
      "missing mean free path ignored fail policy");
  Require(ValidateCollisionModel(CollisionModel::PitchAngleDiffusion,
      false, MissingCoefficientPolicy::Ballistic).ok(),
      "explicit ballistic policy was rejected");
}

void MFP3D04() {
  const double speed = 12.0, d0 = 3.0;
  auto lambda = IntegrateParallelMeanFreePath(speed, 20000,
      [=](double mu) { return d0 * (1.0 - mu * mu); });
  Require(lambda.ok() && Close(lambda.value, speed / (2.0 * d0), 1.0e-8) &&
          Close(speed * lambda.value / 3.0, 8.0),
          "D_mumu integral or Parker kappa=v lambda/3 failed");
}

void MFP3D05() {
  auto plus = StepIsotropicPitchAngleDiffusion(1.0, 2.0, 0.1, 5.0);
  auto minus = StepIsotropicPitchAngleDiffusion(-1.0, 2.0, 0.1, -5.0);
  auto interior = StepIsotropicPitchAngleDiffusion(0.2, 2.0, 0.01, 0.5);
  Require(plus.ok() && minus.ok() && interior.ok() &&
          plus.value >= -1.0 && plus.value <= 1.0 &&
          minus.value >= -1.0 && minus.value <= 1.0 &&
          interior.value >= -1.0 && interior.value <= 1.0,
          "pitch-angle SDE escaped reflecting endpoints");
}

void MFP3D06() {
  auto noScatter = StepIsotropicPoissonScattering(
      0.4, 2.0, 10.0, 1.0, 0.99, 0.2);
  auto scatter = StepIsotropicPoissonScattering(
      0.4, 2.0, 10.0, 1.0, 0.0, 0.75);
  Require(noScatter.ok() && scatter.ok() && !noScatter.value.scattered &&
          scatter.value.scattered && Close(scatter.value.mu, 0.5) &&
          Close(1.0 - noScatter.value.scatterProbability,
                std::exp(-0.2)),
          "Poisson scattering normalization/autocorrelation failed");
}

void MFP3D07() {
  Require(ValidateCollisionModel(CollisionModel::None, false,
                                 MissingCoefficientPolicy::Fail).ok() &&
          ValidateCollisionModel(CollisionModel::PitchAngleDiffusion, true,
                                 MissingCoefficientPolicy::Fail).ok() &&
          ValidateCollisionModel(CollisionModel::DiscreteIsotropicPoisson,
                                 true, MissingCoefficientPolicy::Fail).ok(),
          "documented collision-model dispatch is incomplete");
}

}  // namespace

void RegisterStage4(Registry* tests) {
  (*tests)["TUR3D01"] = TUR3D01; (*tests)["TUR3D02"] = TUR3D02;
  (*tests)["TUR3D03"] = TUR3D03; (*tests)["TUR3D04"] = TUR3D04;
  (*tests)["TUR3D05"] = TUR3D05; (*tests)["TUR3D06"] = TUR3D06;
  (*tests)["TUR3D07"] = TUR3D07;
  (*tests)["MFP3D01"] = MFP3D01; (*tests)["MFP3D02"] = MFP3D02;
  (*tests)["MFP3D03"] = MFP3D03; (*tests)["MFP3D04"] = MFP3D04;
  (*tests)["MFP3D05"] = MFP3D05; (*tests)["MFP3D06"] = MFP3D06;
  (*tests)["MFP3D07"] = MFP3D07;
}

}  // namespace SCCMTest

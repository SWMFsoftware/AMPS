#ifndef SEP_UTIL_SEP_COEFFICIENT_PHYSICS_H
#define SEP_UTIL_SEP_COEFFICIENT_PHYSICS_H

#include "sep_transport_common.h"

#include <cstdint>
#include <functional>
#include <string>

namespace SEP {
namespace Transport {
namespace CoefficientPhysics {

// A floating-point sentinel is not sufficient to distinguish a deliberate
// ballistic limit from absent data or numerical failure.  This enum accompanies
// every repaired coefficient result; value is finite only for Finite and is
// positive infinity only for Ballistic.
enum class ValueState { Finite, Ballistic, Unavailable, Invalid };

// Species data are explicit inputs to every resonant calculation.  Charge is
// signed [C], rest mass is [kg], and nucleonCount is dimensionless.  Electrons
// use nucleonCount=0 because energy-per-nucleon is not meaningful for them.
struct SpeciesProperties {
  int modelSpecies = -1;
  std::string name;
  double signedChargeC = 0.0;
  double restMassKg = 0.0;
  double nucleonCount = 0.0;
};

Status ValidateSpecies(const SpeciesProperties& species);
ScalarResult GyrofrequencyRadPerS(double magneticFieldT,
                                  const SpeciesProperties& species);
ScalarResult LarmorRadiusM(double perpendicularMomentumKgMPerS,
                           double magneticFieldT,
                           const SpeciesProperties& species);
ScalarResult RigidityVolt(double momentumKgMPerS,
                          const SpeciesProperties& species,
                          double speedOfLightMPerS);

// Spectrum bounds are specified at referenceRadiusM and scaled independently.
// Keeping both exponents visible avoids the former mixture of 1/r and 1/r^2
// rules hidden in otherwise similarly named QLT helpers.
struct SpectrumParameters {
  double referenceRadiusM = 1.495978707e11;
  double kMinAtReferencePerM = 1.0e-10;
  double kMaxAtReferencePerM = 1.0e-7;
  double kMinRadialExponent = 2.0;
  double kMaxRadialExponent = 2.0;
  double spectralIndex = 5.0 / 3.0;
};

Status ValidateSpectrum(const SpectrumParameters& spectrum);

// This is the immutable, source-bound input consumed by pure coefficient
// kernels.  Magnetic variances are [T^2], wave energy densities are converted
// by the PIC adapter before this boundary, and all identity fields belong to
// the same background/turbulence generation.
struct LocalInputView {
  std::string source;
  std::string representation;
  std::uint64_t generation = 0;
  std::uint64_t checksum = 0;
  double heliocentricRadiusM = 0.0;
  double magneticFieldT = 0.0;
  double deltaB2T2 = 0.0;
  double deltaBPlus2T2 = 0.0;
  double deltaBMinus2T2 = 0.0;
  double alfvenSpeedMPerS = 0.0;
};

Status ValidateLocalInput(const LocalInputView& input);

enum class DerivativeState { Smooth, OneSidedLimit, Nondifferentiable };

struct PitchAngleResult {
  Status status;
  ValueState valueState = ValueState::Invalid;
  DerivativeState derivativeState = DerivativeState::Smooth;
  double dMuMuPerS = 0.0;
  double dDmuMuDmuPerS = 0.0;
};

// The constant provider is intentionally a pure kernel.  It validates the
// configured SI rate before assigning either output, so negative/NaN values
// cannot masquerade as ballistic transport through a zeroed legacy output.
PitchAngleResult EvaluateConstantDmumu(double configuredDmumuPerS,
                                       double mu);

// Magnetostatic slab QLT with a normalized power-law spectrum.  D and dD/dmu
// are evaluated from one piecewise expression, including resonant-k dependence
// and finite one-sided endpoint limits.  Exact band transitions are returned
// as Nondifferentiable rather than hidden by division through 1-mu^2.
PitchAngleResult EvaluateJokipiiSlab(
    const LocalInputView& input, const SpectrumParameters& spectrum,
    const SpeciesProperties& species, double speedMPerS, double mu);

struct FlorinskiyParameters {
  double spectralIndex = 5.0 / 3.0;
  double parallelCorrelationLengthM = 0.03 * 1.495978707e11;
};

// Dynamic-slab branch resonance used by the historical Florinskiy provider.
// The plus and minus contributions use their own |v*mu -/+ v_A| denominators;
// this removes the previous copy/paste t0 denominator from the minus branch.
// A bounded derivative stencil never evaluates outside physical pitch angle.
PitchAngleResult EvaluateFlorinskiySlab(
    const LocalInputView& input, const FlorinskiyParameters& parameters,
    const SpeciesProperties& species, double speedMPerS, double mu);

// Species-aware QLT1 closure.  momentum is total relativistic momentum [kg
// m/s], not energy per nucleon; signed charge enters through |q| in r_L.
struct MeanFreePathResult {
  Status status;
  ValueState valueState = ValueState::Invalid;
  double lambdaParallelM = 0.0;
};

MeanFreePathResult EvaluateCorrelationMeanFreePath(
    const LocalInputView& input, const SpeciesProperties& species,
    double momentumKgMPerS, double correlationLengthAtReferenceM,
    double referenceRadiusM);

enum class ResonanceGapPolicy { Reject, Ballistic };

struct SpatialQuadratureConfiguration {
  double absoluteToleranceM2PerS = 1.0;
  double relativeTolerance = 1.0e-6;
  int maximumRecursion = 20;
  ResonanceGapPolicy gapPolicy = ResonanceGapPolicy::Reject;
};

struct SpatialDiffusionResult {
  Status status;
  ValueState valueState = ValueState::Invalid;
  double kappaParallelM2PerS = 0.0;
  double estimatedAbsoluteErrorM2PerS = 0.0;
  std::size_t evaluations = 0;
};

typedef std::function<PitchAngleResult(double)> PitchAngleFunction;

// Error-controlled adaptive Simpson quadrature for
// kappa_parallel=(v^2/8) integral_-1^1 (1-mu^2)^2/D_mumu dmu.
// The interval is split at mu=0 and every reported zero inside the physical
// domain is treated according to gapPolicy.  No division by zero or implicit
// floor is permitted.
SpatialDiffusionResult IntegrateSpatialDiffusion(
    double speedMPerS, const PitchAngleFunction& dmumu,
    const SpatialQuadratureConfiguration& configuration);

}  // namespace CoefficientPhysics
}  // namespace Transport
}  // namespace SEP

#endif  // SEP_UTIL_SEP_COEFFICIENT_PHYSICS_H

#ifndef SEP_COMMON_SEP_COEFFICIENT_PHYSICS_H
#define SEP_COMMON_SEP_COEFFICIENT_PHYSICS_H

#include "sep_transport_common.h"

#include <cstdint>
#include <cmath>
#include <functional>
#include <string>

namespace SEP {
namespace Transport {
namespace CoefficientPhysics {

// All SEP branch energies at this boundary are magnetic fluctuation energy,
// as specified by the model owner. This value matches AMPS constants.h; the
// PIC bridge checks the equality at compile time rather than allowing a second
// slightly different permeability to alter an amplitude during conversion.
constexpr double MagneticEnergyPermeabilitySI = 1.256637062e-6;

struct MagneticVariance {
  Status status;
  double plusT2 = 0.0;
  double minusT2 = 0.0;
  double totalT2 = 0.0;
};

// W_B is magnetic energy density [J/m^3], not magnetic-plus-kinetic wave
// energy. Consequently deltaB_branch^2=2*mu0*W_B,branch [T^2]. Zero is a real
// empty branch; negative, nonfinite, overflowing, or underflowed positive input
// is reported explicitly. This inline kernel is the single production formula
// used by the common turbulence core and the PIC coefficient-source adapters.
inline MagneticVariance MagneticVarianceFromEnergyDensity(
    double plusJPerM3, double minusJPerM3) {
  MagneticVariance result;
  if (!std::isfinite(plusJPerM3) || !std::isfinite(minusJPerM3) ||
      plusJPerM3 < 0.0 || minusJPerM3 < 0.0) {
    result.status = Status::Error(StatusCode::InvalidCoefficient,
        "magnetic branch energy densities must be finite and nonnegative");
    return result;
  }
  result.plusT2 = 2.0 * MagneticEnergyPermeabilitySI * plusJPerM3;
  result.minusT2 = 2.0 * MagneticEnergyPermeabilitySI * minusJPerM3;
  result.totalT2 = result.plusT2 + result.minusT2;
  if (!std::isfinite(result.totalT2)) {
    result.status = Status::Error(StatusCode::InvalidCoefficient,
        "magnetic energy-to-variance conversion overflowed");
    return result;
  }
  if ((plusJPerM3 > 0.0 && result.plusT2 == 0.0) ||
      (minusJPerM3 > 0.0 && result.minusT2 == 0.0)) {
    result.status = Status::Error(StatusCode::StepUnderflow,
        "positive magnetic branch energy is below represented variance range");
    return result;
  }
  result.status = Status::Ok();
  return result;
}

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
// Rest/nonrelativistic cyclotron angular frequency Omega0=|q|B/m_rest.
ScalarResult GyrofrequencyRadPerS(double magneticFieldT,
                                  const SpeciesProperties& species);
// Lab/plasma-frame gyration frequency Omega=Omega0/gamma [rad/s].
// v is total particle speed [m/s]; luminal/superluminal values are invalid.
ScalarResult RelativisticGyrofrequencyRadPerS(double magneticFieldT,
    const SpeciesProperties& species,double speedMPerS,
    double speedOfLightMPerS=299792458.0);
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
  // PIC native pitch uses geometric coordinate; pure QLT pitch uses B.
  int magneticPolarityAlongLine = 1;
  bool branchDirectionsKnown = true;
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

// Magnetostatic slab QLT restricted to equal left/right magnetic helicities.
// The positive-k power-law integral is the total real transverse variance;
// each resonant helicity therefore contains half that power. This historical
// scalar API cannot represent independently polarized spectra. D and dD/dmu
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
// Explicit coordinate boundary: mu is relative to increasing field-line
// coordinate. E+ remains along B; the derivative follows the chain rule.
PitchAngleResult EvaluateFlorinskiyAlongFieldLine(
    const LocalInputView& input,const FlorinskiyParameters& parameters,
    const SpeciesProperties& species,double speedMPerS,double muAlongLine);

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

#endif  // SEP_COMMON_SEP_COEFFICIENT_PHYSICS_H

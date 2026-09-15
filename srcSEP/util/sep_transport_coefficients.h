#ifndef SEP_UTIL_SEP_TRANSPORT_COEFFICIENTS_H
#define SEP_UTIL_SEP_TRANSPORT_COEFFICIENTS_H

#include "sep_transport_common.h"
#include "sep_coefficient_physics.h"

#include <cstdint>
#include <string>

namespace SEP {
namespace Transport {

struct SpatialDiffusionSample {
  Status status;
  CoefficientPhysics::ValueState valueState =
      CoefficientPhysics::ValueState::Invalid;
  double kappaParallelM2PerS = 0.0;
  double dKappaParallelDsMPerS = 0.0;
  double estimatedAbsoluteErrorM2PerS = 0.0;
  std::string provenance;
};

class SpatialDiffusionProvider {
 public:
  virtual ~SpatialDiffusionProvider() {}

  // sM is physical arc length [m].  The gradient therefore has units m/s,
  // which makes d(kappa)/ds*dt a physical displacement in the Ito SDE.
  virtual SpatialDiffusionSample Evaluate(double sM,
                                          double speedMPerS) const = 0;
};

struct PitchAngleDiffusionSample {
  Status status;
  CoefficientPhysics::ValueState valueState =
      CoefficientPhysics::ValueState::Invalid;
  CoefficientPhysics::DerivativeState derivativeState =
      CoefficientPhysics::DerivativeState::Smooth;
  double dMuMuPerS = 0.0;
  double dDmuMuDmuPerS = 0.0;
  std::string provenance;
  std::string turbulenceStateIdentity;
  std::string source;
  std::string representation;
  std::uint64_t generation = 0;
  std::uint64_t sourceChecksum = 0;
  double spectrumKMinPerM = 0.0;
  double spectrumKMaxPerM = 0.0;
};

class PitchAngleDiffusionProvider {
 public:
  virtual ~PitchAngleDiffusionProvider() {}

  // Both D_mumu and its mu derivative carry s^-1 because mu is dimensionless.
  // Implementations must identify the turbulence generation used so the mover
  // can prove that scattering and wave-feedback deposition share one state.
  virtual PitchAngleDiffusionSample Evaluate(double sM, double momentumKgMPerS,
                                              double mu) const = 0;
};

struct MeanFreePathSample {
  Status status;
  CoefficientPhysics::ValueState valueState =
      CoefficientPhysics::ValueState::Invalid;
  double lambdaParallelM = 0.0;
  // Branch-resolved event rates are measured in s^-1 in the local plasma
  // frame.  ``nuPlus`` samples waves propagating along the positive magnetic-
  // field direction and ``nuMinus`` samples the opposite branch.  Keeping the
  // two rates in the provider result prevents the mover from inferring a wave
  // population from sign(mu), which is wrong for imbalanced turbulence.
  double nuPlusPerS = 0.0;
  double nuMinusPerS = 0.0;
  bool hasBranchResolvedRates = false;
  double minimumMomentumKgMPerS = 0.0;
  double maximumMomentumKgMPerS = 0.0;
  std::string provenance;
  std::string turbulenceStateIdentity;
  std::string source;
  std::string representation;
  std::uint64_t generation = 0;
  std::uint64_t sourceChecksum = 0;
};

class MeanFreePathProvider {
 public:
  virtual ~MeanFreePathProvider() {}

  // lambda_parallel is measured in metres. Positive infinity is the explicit
  // ballistic limit; finite non-positive values and NaN are invalid. The
  // optional momentum interval is closed and uses SI kg m/s; two zero bounds
  // mean that the provider declares no additional momentum restriction.
  virtual MeanFreePathSample Evaluate(double sM, double momentumKgMPerS,
                                      double mu) const = 0;
};

}  // namespace Transport
}  // namespace SEP

#endif  // SEP_UTIL_SEP_TRANSPORT_COEFFICIENTS_H

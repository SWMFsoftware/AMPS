#ifndef SEP_COMMON_COHERENT_TRANSPORT_H
#define SEP_COMMON_COHERENT_TRANSPORT_H
#include "sep_status.h"
#include <array>
#include <cstdint>
#include <functional>
#include <string>

namespace SEP { namespace Coherent {
enum class Equation { ParkerIsotropic, FocusedGyrotropic };
enum class Ownership { AdditiveParkerAntisymmetricOperator, FullFocusedDeterministicCharacteristic };
enum InvalidRegion : std::uint32_t { Hcs=1, Separatrix=2, Null=4, Transition=8, Shock=16 };
struct Phase {
  Equation equation = Equation::FocusedGyrotropic;
  std::array<double,3> positionM{};
  double momentumSi=0, pitchCosine=0, massKg=0, chargeC=0, epochS=0;
  bool pitchApplicable=true;
  std::string speciesId, equationFrameId;
};
struct SmoothSnapshot {
  std::array<double,3> magneticFieldT{}, electricFieldVPerM{}, gradientMagnitudeTPerM{}, curlUnitFieldPerM{};
  double gradientScaleM=0, curvatureScaleM=0, stepS=0, maximumOrderingError=0;
  std::uint32_t invalidRegion=0;
  bool stationaryFields=true, stationaryInertialEquationFrame=true;
  std::uint64_t generation=0;
  std::string frameId;
};
struct ParkerOperator {
  Ownership ownership=Ownership::AdditiveParkerAntisymmetricOperator;
  std::array<double,9> antisymmetricTensorM2PerS{};
  std::array<double,3> diagnosticDriftVelocityMPerS{};
  double estimatedOrderingError=0;
};
struct FocusedCharacteristic {
  Ownership ownership=Ownership::FullFocusedDeterministicCharacteristic;
  std::array<double,3> positionRateMPerS{};
  double momentumRateSiPerS=0, pitchRatePerS=0, estimatedOrderingError=0;
};
// This qualified mathematical domain is a stationary inertial Hamiltonian
// guiding-center reduction. Moving/plasma-frame fields require a separately
// derived provider; rejecting them prevents adiabatic/electric double work.
Core::Result<ParkerOperator> EvaluateParker(const SmoothSnapshot&, const Phase&);
Core::Result<FocusedCharacteristic> EvaluateFocused(const SmoothSnapshot&, const Phase&);
Core::Result<FocusedCharacteristic> SelectFocused(bool enabled, const SmoothSnapshot&, const Phase&,
    const std::function<Core::Result<FocusedCharacteristic>()>& baseline);
} }
#endif

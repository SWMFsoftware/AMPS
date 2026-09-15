#ifndef SEP_UTIL_SEP_FOCUSED_TRANSPORT_CORE_H
#define SEP_UTIL_SEP_FOCUSED_TRANSPORT_CORE_H

#include "sep_transport_coefficients.h"

#include <string>
#include <vector>

namespace SEP {
namespace Transport {

// FullGyrotropic implements the frozen-background Skilling/Ruffolo focused
// transport coefficients.  ReducedFieldAligned1D is retained only as an
// explicit compatibility gate for campaigns whose documented equation omits
// transverse flow expansion.  Production setup must choose one deliberately.
enum class FocusedEquationMode { FullGyrotropic, ReducedFieldAligned1D };

// The pitch-angle Fokker-Planck operator is in conservative form,
//   partial_t f = partial_mu(D_mumu partial_mu f),
// with zero probability flux at mu=+-1.  Its Ito SDE has drift dD/dmu and
// noise sqrt(2D).  ReflectingMilstein is the production default because the
// derivative-dependent correction improves weak behavior for state-dependent
// D while mirror folding enforces the same no-flux boundary after a finite
// stochastic increment.  ReflectingEulerMaruyama remains only as an explicit
// comparison mode for historical validation.
enum class BoundedPitchAngleScheme {
  ReflectingMilstein,
  ReflectingEulerMaruyama
};

struct FocusedTransportState {
  FocusedTransportState() {}
  FocusedTransportState(double arcLength, double momentum, double pitchMu)
      : arcLengthM(arcLength), momentumKgMPerS(momentum), mu(pitchMu) {}

  double arcLengthM = 0.0;
  double momentumKgMPerS = 0.0;
  double mu = 0.0;
};

struct FocusedTransportBackground {
  FocusedTransportBackground() {}
  FocusedTransportBackground(double dLnBds, double advection,
                             double parallelGradient, double divergence)
      : dLnAbsBdsPerM(dLnBds), plasmaAdvectionMPerS(advection),
        parallelVelocityGradientPerS(parallelGradient),
        velocityDivergencePerS(divergence) {}

  double dLnAbsBdsPerM = 0.0;
  double plasmaAdvectionMPerS = 0.0;
  double parallelVelocityGradientPerS = 0.0;
  double velocityDivergencePerS = 0.0;
  double fieldAlignedStrainPerS = 0.0;
  FocusedEquationMode equationMode =
      FocusedEquationMode::ReducedFieldAligned1D;
  BoundedPitchAngleScheme pitchAngleScheme =
      BoundedPitchAngleScheme::ReflectingMilstein;
};

struct WaveContribution {
  std::string turbulenceStateIdentity;
  double signedStreaming = 0.0;
  double intervalS = 0.0;
  double displacementM = 0.0;
  std::uint64_t eventIndex = 0;
  bool scatteringEventAtEnd = false;
};

// Each worker owns one accumulator.  Sorting on reduction makes the numerical
// result independent of thread completion order for a fixed set of keyed
// particle streams; no mover writes the shared turbulence state directly.
class ThreadLocalWaveAccumulator {
 public:
  void Deposit(const WaveContribution& contribution);
  const std::vector<WaveContribution>& Contributions() const {
    return contributions_;
  }
  static std::vector<WaveContribution> DeterministicReduce(
      const std::vector<const ThreadLocalWaveAccumulator*>& locals);

 private:
  std::vector<WaveContribution> contributions_;
};

struct FocusedTransportIncrement {
  Status status;
  FocusedTransportState state;
  double displacementM = 0.0;
  double deterministicMuIncrement = 0.0;
  double stochasticMuIncrement = 0.0;
  unsigned pitchAngleReflections = 0;
  std::string coefficientProvenance;
};

double ReflectPitchAngle(double mu, unsigned* reflections);

// Public coefficient helpers freeze the governing SDE contract in one tested
// location.  dMu/dt is dimensionless per second; dLnP/dt is per second.
double FocusedPitchDriftPerS(double mu, double speedMPerS,
                            const FocusedTransportBackground& background);
double FocusedLogMomentumRatePerS(
    double mu, const FocusedTransportBackground& background);

// Symmetric deterministic/stochastic/deterministic splitting is declared here
// so every caller uses the same order: half focusing/velocity-gradient drift
// and half cooling, one Ito D_mumu kick, then the remaining deterministic half;
// streaming uses the midpoint pitch angle.  D_mumu is evaluated from one
// provider sample carrying the turbulence-state identity used for deposition.
FocusedTransportIncrement AdvanceFocusedTransportDmumu(
    const FocusedTransportState& initial,
    const FocusedTransportBackground& background,
    double massKg,
    double speedOfLightMPerS,
    double dtS,
    const PitchAngleDiffusionProvider& provider,
    RandomStream& random,
    ThreadLocalWaveAccumulator* waveAccumulator);

}  // namespace Transport
}  // namespace SEP

#endif  // SEP_UTIL_SEP_FOCUSED_TRANSPORT_CORE_H

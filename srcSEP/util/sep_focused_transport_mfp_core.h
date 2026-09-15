#ifndef SEP_UTIL_SEP_FOCUSED_TRANSPORT_MFP_CORE_H
#define SEP_UTIL_SEP_FOCUSED_TRANSPORT_MFP_CORE_H

#include "sep_focused_transport_core.h"

#include <cstddef>
#include <cstdint>
#include <limits>

namespace SEP {
namespace Transport {

struct FocusedTransportMfpState {
  FocusedTransportMfpState() {}
  FocusedTransportMfpState(double arcLength, double momentum, double pitchMu)
      : arcLengthM(arcLength), momentumKgMPerS(momentum), mu(pitchMu) {}

  double arcLengthM = 0.0;
  double momentumKgMPerS = 0.0;
  double mu = 0.0;
  // Unit-exponential optical depth remaining to the next event.  NaN means
  // "not drawn yet" (new particle or backward-compatible restart).  Carrying
  // this state makes segment/snapshot subcycling exactly partition invariant.
  double remainingOpticalDepth = std::numeric_limits<double>::quiet_NaN();
  std::uint64_t nextEventIndex = 0;
};

struct FocusedTransportMfpBackground {
  FocusedTransportMfpBackground() {}
  FocusedTransportMfpBackground(double dLnBds, double advection,
                                double parallelGradient, double divergence,
                                double alfvenSpeed)
      : dLnAbsBdsPerM(dLnBds), plasmaAdvectionMPerS(advection),
        parallelVelocityGradientPerS(parallelGradient),
        velocityDivergencePerS(divergence), alfvenSpeedMPerS(alfvenSpeed) {}

  // Particle velocity and pitch angle are in the local plasma frame. The
  // signed Alfven speed is not stored: the resonance rule selects the branch
  // opposite to the particle's parallel direction at each event.
  double dLnAbsBdsPerM = 0.0;
  double plasmaAdvectionMPerS = 0.0;
  double parallelVelocityGradientPerS = 0.0;
  double velocityDivergencePerS = 0.0;
  double alfvenSpeedMPerS = 0.0;
  double fieldAlignedStrainPerS = 0.0;
  FocusedEquationMode equationMode =
      FocusedEquationMode::ReducedFieldAligned1D;
};

struct MfpEventDiagnostics {
  std::size_t deterministicIntervals = 0;
  std::size_t scatteringEvents = 0;
  std::size_t ballisticIntervals = 0;
  std::size_t providerEvaluations = 0;
  std::size_t plusBranchEvents = 0;
  std::size_t minusBranchEvents = 0;
  std::size_t opticalDepthRootIterations = 0;
};

struct WaveFrameScatterResult {
  Status status;
  double speedMPerS = 0.0;
  double mu = 0.0;
  double waveFrameSpeedBeforeMPerS = 0.0;
  double waveFrameSpeedAfterMPerS = 0.0;
};

struct FocusedTransportMfpIncrement {
  Status status;
  FocusedTransportMfpState state;
  double displacementM = 0.0;
  MfpEventDiagnostics diagnostics;
  std::string coefficientProvenance;
};

// Draw t=-ln(xi)/nu from an open-interval uniform stream. A zero rate returns
// positive infinity, which is the declared no-event/ballistic limit.
ScalarResult SampleExponentialWaitingTime(double ratePerS,
                                          RandomStream& random);

// Redistribute pitch angle isotropically in the selected Alfven-wave frame.
// The exact one-dimensional Lorentz velocity transformation preserves the
// particle speed in that wave frame and returns a strictly subluminal plasma-
// frame state. waveSpeedMPerS is signed in the positive field-line direction.
WaveFrameScatterResult ScatterIsotropicallyInWaveFrame(
    double speedMPerS, double mu, double waveSpeedMPerS,
    double speedOfLightMPerS, RandomStream& random);

// Event-driven focused transport with the explicit closure nu=v/lambda. The
// maximum deterministic interval is supplied by the adapter after composing
// segment, focusing, cooling, snapshot-validity, and shock-crossing limits.
// A unit-exponential optical depth is carried in FocusedTransportMfpState and
// consumed by the integrated, location-dependent total hazard.  Deterministic
// splits therefore never consume a new random variate or change event timing.
FocusedTransportMfpIncrement AdvanceFocusedTransportMfp(
    const FocusedTransportMfpState& initial,
    const FocusedTransportMfpBackground& background,
    double massKg,
    double speedOfLightMPerS,
    double dtS,
    double maximumDeterministicIntervalS,
    const MeanFreePathProvider& provider,
    RandomStream& random,
    ThreadLocalWaveAccumulator* waveAccumulator);

}  // namespace Transport
}  // namespace SEP

#endif  // SEP_UTIL_SEP_FOCUSED_TRANSPORT_MFP_CORE_H

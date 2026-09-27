// ============================================================================
// AMPS-independent gyrotropic focused-transport core.
//
// The frozen-local-background SDE is advanced with a symmetric split:
// deterministic focusing/flow/cooling half-step, one Ito D_mumu kick, then a
// second deterministic half-step.  The production stochastic kick is scalar
// Milstein.  Mirror folding at mu=+-1 implements the no-probability-flux
// boundary without clipping particles onto the endpoints.
//
// In the local plasma frame the deterministic coefficients are
//
// dmu/dt = (1-mu^2)/2 [v div(b) + mu(divU-3 bb:gradU)],
// dlnp/dt = -1/2 [(1-mu^2)divU +(3mu^2-1)bb:gradU].
//
// Isotropic pitch averaging of the second equation recovers Parker cooling,
// dlnp/dt=-divU/3. Spatial motion is U + mu v b plus the controlled V01
// guiding-centre drift and two separately keyed perpendicular Wiener terms.
// ============================================================================

#ifndef SEP3D_TRANSPORT_FOCUSED_TRANSPORT_H
#define SEP3D_TRANSPORT_FOCUSED_TRANSPORT_H

#include "../core/sep3d_types.h"
#include "keyed_random.h"
#include "perpendicular_transport.h"

#include <cstdint>
#include <limits>

namespace SEP3D {
namespace Transport {

enum class PitchAngleScheme { ReflectingMilstein, ReflectingEulerMaruyama };

struct FocusedParticleState {
  Core::Vec3 positionM;
  double momentumKgMPerS = 0.0;
  double mu = 0.0;
};

struct FocusedLocalState {
  Core::Vec3 bulkVelocityMPerS;
  Core::Vec3 bHat;
  double divBhatPerM = 0.0;
  double divUPerS = 0.0;
  double fieldAlignedStrainPerS = 0.0;
  double dMuMuPerS = 0.0;
  double dDmuMuDmuPerS = 0.0;
  PitchAngleScheme scheme = PitchAngleScheme::ReflectingMilstein;
  double kappaPerpendicularM2PerS = 0.0;
  Core::Vec3 driftVelocityMPerS;
};

struct FocusedRandomStreams {
  KeyedRandomStream* pitch = nullptr;
  KeyedRandomStream* perpendicularFirst = nullptr;
  KeyedRandomStream* perpendicularSecond = nullptr;
};

struct FocusedStepResult {
  Core::Status status;
  FocusedParticleState state;
  Core::Vec3 displacementM;
  double deterministicMuIncrement = 0.0;
  double stochasticMuIncrement = 0.0;
  double logMomentumIncrement = 0.0;
  unsigned pitchReflections = 0;
};

enum class FocusedScatteringFrame {
  PlasmaFrameIsotropic,
  AlfvenWaveFrameIsotropic
};

// Persistent state for the event-driven mean-free-path formulation.  Optical
// depth, rather than a sampled waiting time, is carried across AMR and output
// boundaries so repartitioning a requested interval cannot redraw an event.
struct FocusedScatteringState {
  FocusedParticleState particle;
  double remainingOpticalDepth =
      std::numeric_limits<double>::quiet_NaN();
  std::uint64_t nextEventIndex = 0;
};

struct FocusedScatteringLocalState {
  FocusedLocalState deterministic;
  double meanFreePathM = 0.0;
  double alfvenSpeedMPerS = 0.0;
  // Fractions are the directional magnetic-variance partition and must sum
  // to one.  They select the resonant wave frame without assuming balance.
  double plusWaveFraction = 0.5;
  double minusWaveFraction = 0.5;
  FocusedScatteringFrame frame =
      FocusedScatteringFrame::PlasmaFrameIsotropic;
};

struct FocusedScatteringStepResult {
  Core::Status status;
  FocusedScatteringState state;
  Core::Vec3 displacementM;
  std::uint64_t scatteringEvents = 0;
};

double RelativisticSpeed(double momentumKgMPerS, double massKg);
double FocusedPitchDriftPerS(double mu, double speedMPerS,
                             const FocusedLocalState& local);
double FocusedLogMomentumRatePerS(double mu,
                                  const FocusedLocalState& local);
double ReflectPitchAngle(double mu, unsigned* reflections = nullptr);

FocusedStepResult AdvanceFocused(const FocusedParticleState& initial,
                                  const FocusedLocalState& local,
                                  double massKg,
                                  double dtS,
                                  const FocusedRandomStreams& random);
// Compatibility overload for callers selecting no perpendicular diffusion.
FocusedStepResult AdvanceFocused(const FocusedParticleState& initial,
                                  const FocusedLocalState& local,
                                  double massKg,
                                  double dtS,
                                  KeyedRandomStream* random);

// Advance a piecewise-constant local background with the explicit event rate
// nu=v/lambda_parallel.  Deterministic focusing, streaming, and cooling use
// the same reviewed focused core as the D_mumu mover.  Each event receives
// independent semantic random keys based on nextEventIndex.
FocusedScatteringStepResult AdvanceFocusedScattering(
    const FocusedScatteringState& initial,
    const FocusedScatteringLocalState& local,
    double massKg, double dtS, std::uint64_t campaignSeed,
    std::uint64_t stableParticleId, std::uint64_t completedStep,
    std::uint64_t maximumEvents);

}  // namespace Transport
}  // namespace SEP3D

#endif  // SEP3D_TRANSPORT_FOCUSED_TRANSPORT_H

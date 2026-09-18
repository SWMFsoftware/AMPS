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

}  // namespace Transport
}  // namespace SEP3D

#endif  // SEP3D_TRANSPORT_FOCUSED_TRANSPORT_H

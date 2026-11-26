// ============================================================================
// Phase-A host-neutral particle mover adapter.
//
// This layer is deliberately free of pic.h.  It validates the complete record
// produced by an AMPS (or test) buffer translator, selects exactly one Phase-P
// core, classifies physical boundaries, and records a possible expanding-shock
// intersection.  It contains no duplicate Parker or focused equations.
// ============================================================================

#ifndef SEP3D_ADAPTERS_TRANSPORT_ADAPTER_H
#define SEP3D_ADAPTERS_TRANSPORT_ADAPTER_H

#include "../background/bg_provider.h"
#include "../runtime/run_configuration.h"
#include "../transport/focused_transport.h"
#include "../transport/parker_transport.h"
#include "../transport/time_step.h"
#include "shock_geometry.h"

#include <cstdint>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace SEP3D {
namespace Adapters {

struct ParticleRecord {
  std::uint64_t stableId = 0;
  int species = -1;
  Core::Vec3 positionM;
  double momentumKgMPerS = 0.0;
  double mu = 0.0;
  double gyrophaseRad = 0.0;
  double statisticalWeight = 0.0;
  std::uint64_t completedStep = 0;
  std::uint64_t substep = 0;
  std::uint64_t lastShockGeneration = 0;
  // Event-driven focused-scattering state.  NaN denotes a new particle whose
  // first unit-exponential optical depth has not yet been drawn.
  double remainingScatteringOpticalDepth =
      std::numeric_limits<double>::quiet_NaN();
  std::uint64_t nextScatteringEvent = 0;
};

struct LocalTransportRecord {
  Background::BackgroundSample background;
  double cellSizeM = 0.0;
  double kappaParallelM2PerS = 0.0;
  double meanFreePathM = 0.0;
  double dKappaParallelDsMPerS = 0.0;
  double dMuMuPerS = 0.0;
  double dDmuMuDmuPerS = 0.0;
  // Scattering step limit [s] supplied by the pitch-angle library's step rule
  // (focused-diffusion with pitch_angle_diffusion_model = pitch-angle-library);
  // +inf otherwise, which leaves the substep selection of the other D_mumu
  // models unchanged.
  double pitchAngleScatteringLimitS = std::numeric_limits<double>::infinity();
  std::uint32_t parallelDiffusionDiagnosticMask = 0;
  std::string parallelDiffusionModelId;
  std::string parallelDiffusionConfigurationFingerprint;
  double plusWaveFraction = 0.5;
  double minusWaveFraction = 0.5;
  double fractionalFieldVariationPerS = 0.0;
  double timeToSnapshotBoundaryS = 0.0;
};

struct ShockIntersection {
  Core::Status status;
  bool crossed = false;
  double stepFraction = 0.0;
  Core::Vec3 positionM;
  // Geometry at the accepted crossing, including SSE's local normal speed.
  // These values do not assert that a geometric flank is super-fast.
  Core::Vec3 outwardNormal;
  double normalSpeedMPerS = 0.0;
  std::uint64_t generation = 0;
};

enum class ParticleDisposition { Active, Escaped, Absorbed, Failed };

const char* Name(ParticleDisposition disposition);

struct MoverInput {
  RuntimeModel::TransportModel model = RuntimeModel::TransportModel::Parker3D;
  ParticleRecord particle;
  LocalTransportRecord local;
  ExpandingShock shock;
  double speciesMassKg = 0.0;
  double speciesChargeC = 0.0;
  double requestedDtS = 0.0;
  double innerRadiusM = 0.0;
  double outerRadiusM = 0.0;
  std::uint64_t campaignSeed = 0;
  Transport::TimeStepControls timeStepControls;
  RuntimeModel::PitchAngleSchemeMode pitchScheme =
      RuntimeModel::PitchAngleSchemeMode::ReflectingMilstein;
  RuntimeModel::PerpendicularDiffusionMode perpendicularDiffusion =
      RuntimeModel::PerpendicularDiffusionMode::None;
  double constantKappaPerpendicularM2PerS = 0.0;
  double kappaPerpendicularToParallelRatio = 0.0;
  RuntimeModel::DriftMode drift = RuntimeModel::DriftMode::None;
  RuntimeModel::FocusedScatteringFrame focusedScatteringFrame =
      RuntimeModel::FocusedScatteringFrame::PlasmaFrameIsotropic;
  std::uint64_t maximumScatteringEventsPerSubstep = 100000;
  // Absorbing shock (source.shock.returning_particles = absorb, 2026-11-26).
  // A particle that crosses the front into its downstream side, where the
  // front is an accepted shock, is removed: the 1-D net source already counts
  // only escaping particles.  moverEpochS is the model time at the start of
  // this (sub)step; shockAccepted answers whether the front point is an
  // accepted shock at an epoch.
  //   absorbAtShock  enables the test after every substep (AdvanceParticle;
  //                  geometry and statistics in shock_return.h).  The front
  //                  geometry is `shock`, whose radiusAtStepStartM the caller
  //                  sets to the radius at moverEpochS.
  //   moverEpochS    model time [s] at the start of this (sub)step; the hit
  //                  epoch passed to shockAccepted is moverEpochS + tau dt.
  //   shockAccepted  optional query (frame point [m], epoch [s]) ->
  //                  accepted; main_lib.cpp binds it to the ShockSurface
  //                  front state (supported && accepted).  nullptr accepts
  //                  every point of the geometric front.  A failed query
  //                  fails the particle.
  bool absorbAtShock = false;
  double moverEpochS = 0.0;
  // Rotating-frame inertial forces of the focused movers (corotating frame).
  //   inertialForces  application.sep.transport.frame_inertial_forces:
  //                   None, Coriolis, or CoriolisAndCentrifugal.  Only the
  //                   two focused movers evaluate it; Parker ignores it.
  //   frameAngularVelocityRadPerS  Omega of the computation frame [rad/s]
  //                   (0 in a non-rotating frame, where the terms vanish).
  // The terms are described at EvaluateFrameInertialTerms
  // (transport_adapter.cpp) and in transport/focused_transport.h.
  RuntimeModel::FrameInertialForces inertialForces =
      RuntimeModel::FrameInertialForces::None;
  Core::Vec3 frameAngularVelocityRadPerS;
  Core::Status (*shockAccepted)(const Core::Vec3& frontPointM, double epochS,
                                bool* accepted) = nullptr;
};

struct MoverResult {
  Core::Status status;
  ParticleRecord particle;
  ParticleDisposition disposition = ParticleDisposition::Failed;
  Transport::TimeStepSelection selectedStep;
  ShockIntersection shockIntersection;
  // R02 accounting fields. AdvanceParticle() reports one accepted substep;
  // AdvanceParticleRequestedTime() accumulates these fields over the complete
  // AMPS request.  Keeping the consumed time explicit makes a terminal
  // boundary crossing distinguishable from an accidental short advance.
  double consumedTimeS = 0.0;
  std::uint64_t acceptedSubsteps = 0;
  LocalTransportRecord finalLocal;
};

// Resolve the immutable cell-local state at the particle's *current* position.
// elapsedTimeS is measured from the beginning of this AMPS mover call.  The
// callback may use it to select the correct staged snapshot, but it must never
// mutate the particle.  A void context keeps this AMPS-independent while
// allowing the production boundary to carry its current AMR-node hint.
using LocalTransportResolver = Core::Status (*)(
    const ParticleRecord& particle, double elapsedTimeS,
    void* resolverContext, LocalTransportRecord* local);

struct RequestedTimeAdvance {
  MoverInput input;
  LocalTransportResolver resolveLocal = nullptr;
  void* resolverContext = nullptr;
  std::uint64_t maximumSubsteps = 100000;
};

// Solve the first intersection between a straight particle substep and a
// Sphere or finite SSE front with linearly evolving apex over the interval.
// Both quadratic roots are checked: an SSE rear-sphere hit is rejected before
// selecting the first outward-cap root. A generation recorded by the particle
// is suppressed so a trajectory cannot inject twice on one shock surface.
ShockIntersection FirstShockIntersection(
    const Core::Vec3& initialPositionM,
    const Core::Vec3& finalPositionM,
    double dtS,
    const ExpandingShock& shock,
    std::uint64_t lastShockGeneration);

class ProductionMoverRegistry final {
 public:
  static const std::vector<std::string>& CanonicalNames();
  static Core::Status Validate(RuntimeModel::TransportModel selected);
};

MoverResult AdvanceParticle(const MoverInput& input);

// Consume one complete host-requested interval through accepted local
// substeps.  The resolver is called before every substep, so a particle that
// enters another AMR cell cannot continue with stale B, U, gradients, or
// scattering coefficients.  Active results consume requestedDtS to a stable
// floating-point tolerance; terminal results report the physically consumed
// prefix, and failures report the prefix completed before the error.
MoverResult AdvanceParticleRequestedTime(const RequestedTimeAdvance& request);

}  // namespace Adapters
}  // namespace SEP3D

#endif  // SEP3D_ADAPTERS_TRANSPORT_ADAPTER_H

#ifndef SEP_TRANSPORT_COMMON_H
#define SEP_TRANSPORT_COMMON_H

#include "sep.h"

#include <cstdint>
#include <string>

namespace SEP {
namespace Transport {
namespace PICAdapter {

struct ParticleContext {
  long int pointer = -1;
  PIC::ParticleBuffer::byte* data = NULL;
  PIC::FieldLine::cFieldLineSegment* segment = NULL;
  ParticleState state;
  std::uint64_t stableParticleId = 0;
  std::uint64_t snapshotGeneration = 0;
  double particleStepEpochS = 0.0;
  double statisticalWeight = 0.0;
};

// All three production movers enter through these helpers.  Validation and
// attachment therefore cannot diverge as the mover implementations evolve.
Status LoadParticle(long int pointer, ParticleContext* context);
// Injection sites call this once with a physical source/sequence key and the
// generated three-component momentum [kg m/s].  The resulting ID is persisted
// in the AMPS particle record and remains independent of later buffer moves.
void InitializeParticleTransportState(long int pointer,
                                      std::uint64_t sourceKey,
                                      std::uint64_t sequenceKey,
                                      const double momentumKgMPerS[3]);
Status AdvanceAlongFieldLine(ParticleContext* context,
                             double displacementM);
Status CommitAndAttach(const ParticleContext& context);

// Clip an attempted physical displacement to the first absorbing field-line
// boundary.  The returned coordinate is an exact segment endpoint and the
// signed path contains only in-domain work; callers scale interval time by the
// path ratio before enqueuing a BoundaryExit coupling record.
Status ClipPathToFieldLineBoundary(int fieldLineId, double startCoordinate,
                                   double attemptedDisplacementM,
                                   double* boundaryCoordinate,
                                   double* inDomainDisplacementM);

struct LocalBackgroundView {
  // All vectors and scalars are SI and belong to one immutable snapshot.
  // ``tangent`` and ``bUnit`` are kept separately because field-line geometry
  // and magnetic polarity need not share the same orientation.
  double magneticFieldT[3] = {0.0, 0.0, 0.0};
  double magneticFieldMagnitudeT = 0.0;
  double tangent[3] = {0.0, 0.0, 0.0};
  double bUnit[3] = {0.0, 0.0, 0.0};
  double plasmaVelocityMPerS[3] = {0.0, 0.0, 0.0};
  double numberDensityPerM3 = 0.0;
  double massDensityKgPerM3 = 0.0;
  double dLnAbsBdsPerM = 0.0;
  double plasmaAdvectionMPerS = 0.0;
  double parallelVelocityGradientPerS = 0.0;
  double velocityDivergencePerS = 0.0;
  double fieldAlignedStrainPerS = 0.0;
  double alfvenSpeedMPerS = 0.0;
  double snapshotSecondsRemaining = 0.0;
  double sampleCoordinate = 0.0;
  std::uint64_t generation = 0;
  std::string identity;
  std::string provenance;
};

typedef LocalBackgroundView LocalBackground;

// Reads geometry, flow, and density from one segment while the immutable
// BackgroundSnapshot read phase is active.  The derivative convention is the
// same for Parker and both focused movers.
Status EvaluateLocalBackgroundAt(const ParticleContext& context,
                                 double relativeArcLengthM,
                                 LocalBackgroundView* background);
Status EvaluateLocalBackground(const ParticleContext& context,
                               LocalBackgroundView* background);

extern std::uint64_t CampaignRandomSeed;

enum class CouplingEventType {
  DeterministicInterval = 0,
  ScatteringEvent = 1,
  BoundaryExit = 2
};

struct CouplingRecord {
  // This payload is self-contained by design.  FlushWaveContributions never
  // dereferences a particle-buffer handle, so an absorbing-boundary exit cannot
  // discard already completed in-domain work or leave a dangling pointer.
  std::string turbulenceStateIdentity;
  int fieldLineId = -1;
  int species = -1;
  double statisticalWeight = 0.0;
  std::uint64_t stableParticleId = 0;
  std::uint64_t snapshotGeneration = 0;
  std::uint64_t eventIndex = 0;
  std::uint64_t intervalIndex = 0;
  double dtS = 0.0;
  double preMomentumKgMPerS = 0.0;
  double postMomentumKgMPerS = 0.0;
  double midpointParallelVelocityMPerS = 0.0;
  double midpointNormalVelocityMPerS = 0.0;
  double startCoordinate = 0.0;
  double finishCoordinate = 0.0;
  double signedPathM = 0.0;
  int resonantBranch = 0;
  CouplingEventType eventType = CouplingEventType::DeterministicInterval;
  bool pitchAngleResolved = false;
  bool crossedBoundary = false;
};

// Movers enqueue coupling records in worker-local storage.  The driver calls
// FlushWaveContributions() after PIC::TimeStep(), when no particle thread is
// writing, to sort and apply them deterministically to the legacy G+/G- arrays.
void QueueWaveContribution(const CouplingRecord& record);
void FlushWaveContributions();

}  // namespace PICAdapter
}  // namespace Transport
}  // namespace SEP

#endif  // SEP_TRANSPORT_COMMON_H

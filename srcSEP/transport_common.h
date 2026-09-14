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
};

// All three production movers enter through these helpers.  Validation and
// attachment therefore cannot diverge as the mover implementations evolve.
Status LoadParticle(long int pointer, ParticleContext* context);
Status AdvanceAlongFieldLine(ParticleContext* context,
                             double displacementM);
Status CommitAndAttach(const ParticleContext& context);

struct LocalBackground {
  double dLnAbsBdsPerM = 0.0;
  double plasmaAdvectionMPerS = 0.0;
  double parallelVelocityGradientPerS = 0.0;
  double velocityDivergencePerS = 0.0;
  double alfvenSpeedMPerS = 0.0;
  double snapshotSecondsRemaining = 0.0;
  std::string identity;
};

// Reads geometry, flow, and density from one segment while the immutable
// BackgroundSnapshot read phase is active.  The derivative convention is the
// same for Parker and both focused movers.
Status EvaluateLocalBackground(const ParticleContext& context,
                               double densityEvolutionIntervalS,
                               LocalBackground* background);

extern std::uint64_t CampaignRandomSeed;

// Movers enqueue coupling records in worker-local storage.  The driver calls
// FlushWaveContributions() after PIC::TimeStep(), when no particle thread is
// writing, to sort and apply them deterministically to the legacy G+/G- arrays.
void QueueAveragedWaveContribution(int fieldLineId, long int particlePointer,
                                   double dtS, double speedMPerS,
                                   double startCoordinate,
                                   double finishCoordinate,
                                   double signedPathM,
                                   std::uint64_t eventIndex);
void QueueFocusedWaveContribution(const std::string& turbulenceStateIdentity,
                                  int fieldLineId, long int particlePointer,
                                  double dtS, double vParallelMPerS,
                                  double vNormalMPerS,
                                  double startCoordinate,
                                  double finishCoordinate,
                                  double signedPathM,
                                  std::uint64_t eventIndex);
void FlushWaveContributions();

}  // namespace PICAdapter
}  // namespace Transport
}  // namespace SEP

#endif  // SEP_TRANSPORT_COMMON_H

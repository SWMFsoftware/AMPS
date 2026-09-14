#include "transport_common.h"

#include "amps2swmf.h"

#include <cmath>
#include <algorithm>
#include <mutex>
#include <sstream>
#include <vector>

namespace SEP {
namespace Transport {
namespace PICAdapter {

std::uint64_t CampaignRandomSeed = 0;

namespace {

struct PendingWaveContribution {
  std::string turbulenceStateIdentity;
  int fieldLineId;
  long int particlePointer;
  double dtS;
  double speedOrParallelMPerS;
  double normalMPerS;
  double startCoordinate;
  double finishCoordinate;
  double signedPathM;
  std::uint64_t eventIndex;
  bool pitchAngleResolved;
};

struct PendingQueue {
  std::vector<PendingWaveContribution> values;
};

std::mutex gQueueRegistrationMutex;
std::vector<PendingQueue*> gRegisteredQueues;

PendingQueue* RegisterThreadQueue() {
  // Queues intentionally have process lifetime.  This avoids dangling registry
  // entries when an OpenMP worker is retired and is bounded by the number of
  // worker threads, not by the number of particles or timesteps.
  PendingQueue* queue = new PendingQueue;
  std::lock_guard<std::mutex> lock(gQueueRegistrationMutex);
  gRegisteredQueues.push_back(queue);
  return queue;
}

PendingQueue& LocalQueue() {
  thread_local PendingQueue* queue = RegisterThreadQueue();
  return *queue;
}

bool PendingLess(const PendingWaveContribution& left,
                 const PendingWaveContribution& right) {
  if (left.turbulenceStateIdentity != right.turbulenceStateIdentity)
    return left.turbulenceStateIdentity < right.turbulenceStateIdentity;
  if (left.fieldLineId != right.fieldLineId)
    return left.fieldLineId < right.fieldLineId;
  if (left.particlePointer != right.particlePointer)
    return left.particlePointer < right.particlePointer;
  if (left.eventIndex != right.eventIndex)
    return left.eventIndex < right.eventIndex;
  return left.startCoordinate < right.startCoordinate;
}

}  // namespace

void QueueAveragedWaveContribution(int fieldLineId, long int particlePointer,
                                   double dtS, double speedMPerS,
                                   double startCoordinate,
                                   double finishCoordinate,
                                   double signedPathM,
                                   std::uint64_t eventIndex) {
  PendingWaveContribution value;
  value.turbulenceStateIdentity = "pitch-angle-averaged";
  value.fieldLineId = fieldLineId;
  value.particlePointer = particlePointer;
  value.dtS = dtS;
  value.speedOrParallelMPerS = speedMPerS;
  value.normalMPerS = 0.0;
  value.startCoordinate = startCoordinate;
  value.finishCoordinate = finishCoordinate;
  value.signedPathM = signedPathM;
  value.eventIndex = eventIndex;
  value.pitchAngleResolved = false;
  LocalQueue().values.push_back(value);
}

void QueueFocusedWaveContribution(const std::string& turbulenceStateIdentity,
                                  int fieldLineId, long int particlePointer,
                                  double dtS, double vParallelMPerS,
                                  double vNormalMPerS,
                                  double startCoordinate,
                                  double finishCoordinate,
                                  double signedPathM,
                                  std::uint64_t eventIndex) {
  PendingWaveContribution value;
  // Preserve the provider's immutable identity through the worker queue.  The
  // physical legacy accumulator cannot consume this label yet, but ordering by
  // it prevents contributions from distinct published states from interleaving
  // if a future driver intentionally batches more than one snapshot.
  value.turbulenceStateIdentity = turbulenceStateIdentity;
  value.fieldLineId = fieldLineId;
  value.particlePointer = particlePointer;
  value.dtS = dtS;
  value.speedOrParallelMPerS = vParallelMPerS;
  value.normalMPerS = vNormalMPerS;
  value.startCoordinate = startCoordinate;
  value.finishCoordinate = finishCoordinate;
  value.signedPathM = signedPathM;
  value.eventIndex = eventIndex;
  value.pitchAngleResolved = true;
  LocalQueue().values.push_back(value);
}

void FlushWaveContributions() {
  std::vector<PendingWaveContribution> merged;
  {
    // PIC::TimeStep() has joined its worker region before this call, so the
    // queues are quiescent.  The mutex protects only registry traversal and
    // future worker registration; particle stepping never contends on it.
    std::lock_guard<std::mutex> lock(gQueueRegistrationMutex);
    for (PendingQueue* queue : gRegisteredQueues) {
      if (!queue) continue;
      merged.insert(merged.end(), queue->values.begin(), queue->values.end());
      queue->values.clear();
    }
  }
  std::sort(merged.begin(), merged.end(), PendingLess);
  for (const PendingWaveContribution& value : merged) {
    if (value.pitchAngleResolved) {
      SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::
          AccumulateParticleFluxForWaveCoupling(
              value.fieldLineId, value.particlePointer, value.dtS,
              value.speedOrParallelMPerS, value.normalMPerS,
              value.startCoordinate, value.finishCoordinate,
              value.signedPathM);
    }
    else {
      SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::
          AccumulateParticleFluxForWaveCoupling(
              value.fieldLineId, value.particlePointer, value.dtS,
              value.speedOrParallelMPerS, value.startCoordinate,
              value.finishCoordinate, value.signedPathM);
    }
  }
}

Status LoadParticle(long int pointer, ParticleContext* context) {
  if (!context || pointer < 0) {
    return Status::Error(StatusCode::InvalidArgument,
                         "particle context and handle must be valid");
  }
  context->pointer = pointer;
  context->data = PIC::ParticleBuffer::GetParticleDataPointer(pointer);
  if (!context->data) {
    return Status::Error(StatusCode::InvalidParticleState,
                         "particle buffer returned a null record");
  }

  context->state.species = PIC::ParticleBuffer::GetI(context->data);
  context->state.fieldLineId = PIC::ParticleBuffer::GetFieldLineId(context->data);
  context->state.coordinate =
      PIC::ParticleBuffer::GetFieldLineCoord(context->data);
  context->state.vParallelMPerS =
      PIC::ParticleBuffer::GetVParallel(context->data);
  context->state.vNormalMPerS =
      PIC::ParticleBuffer::GetVNormal(context->data);

  if (context->state.species < 0 ||
      context->state.species >= PIC::nTotalSpecies) {
    return Status::Error(StatusCode::InvalidParticleState,
                         "particle species is outside the configured table");
  }
  context->state.massKg =
      PIC::MolecularData::GetMass(context->state.species);
  const Status validation = ValidateParticleState(context->state, SpeedOfLight);
  if (!validation.ok()) return validation;

  if (!PIC::FieldLine::FieldLinesAll || context->state.fieldLineId < 0 ||
      context->state.fieldLineId >= PIC::FieldLine::nFieldLine) {
    return Status::Error(StatusCode::InvalidParticleState,
                         "particle field-line identifier is invalid");
  }
  context->segment =
      PIC::FieldLine::FieldLinesAll[context->state.fieldLineId].GetSegment(
          context->state.coordinate);
  if (!context->segment) {
    return Status::Error(StatusCode::OutOfDomain,
                         "particle coordinate is outside its field line");
  }
  return Status::Ok();
}

Status AdvanceAlongFieldLine(ParticleContext* context,
                             double displacementM) {
  if (!context || !context->segment || !std::isfinite(displacementM)) {
    return Status::Error(StatusCode::InvalidArgument,
                         "field-line advance requires a finite displacement");
  }

  // PIC::FieldLine::move consumes physical distance [m] and applies the line's
  // segment metric internally.  No mover may add a raw metre displacement to
  // the dimensionless segment-plus-fraction coordinate.
  context->state.coordinate =
      PIC::FieldLine::FieldLinesAll[context->state.fieldLineId].move(
          context->state.coordinate, displacementM, context->segment);
  context->segment =
      PIC::FieldLine::FieldLinesAll[context->state.fieldLineId].GetSegment(
          context->state.coordinate);
  if (!context->segment) {
    return Status::Error(StatusCode::OutOfDomain,
                         "particle crossed an absorbing field-line boundary");
  }
  return Status::Ok();
}

Status CommitAndAttach(const ParticleContext& context) {
  if (!context.data || !context.segment ||
      !std::isfinite(context.state.coordinate) ||
      !std::isfinite(context.state.vParallelMPerS) ||
      !std::isfinite(context.state.vNormalMPerS)) {
    return Status::Error(StatusCode::InvalidParticleState,
                         "cannot commit an invalid field-line particle");
  }

  PIC::ParticleBuffer::SetVParallel(context.state.vParallelMPerS, context.data);
  PIC::ParticleBuffer::SetVNormal(context.state.vNormalMPerS, context.data);
  PIC::ParticleBuffer::SetFieldLineCoord(context.state.coordinate, context.data);

  // The temporary list is the only production attachment destination in
  // srcSEP.  Atomic exchange preserves the existing OpenMP insertion contract.
  const long int previous =
      context.segment->tempFirstParticleIndex.exchange(context.pointer);
  PIC::ParticleBuffer::SetNext(previous, context.data);
  PIC::ParticleBuffer::SetPrev(-1, context.data);
  if (previous != -1) PIC::ParticleBuffer::SetPrev(context.pointer, previous);
  return Status::Ok();
}

Status EvaluateLocalBackground(const ParticleContext& context,
                               double densityEvolutionIntervalS,
                               LocalBackground* background) {
  namespace FL = PIC::FieldLine;
  if (!background || !context.segment ||
      !std::isfinite(densityEvolutionIntervalS) ||
      densityEvolutionIntervalS <= 0.0) {
    return Status::Error(StatusCode::InvalidArgument,
                         "local background requires a segment and positive interval");
  }

  FL::cFieldLineVertex* begin = context.segment->GetBegin();
  FL::cFieldLineVertex* end = context.segment->GetEnd();
  if (!begin || !end) {
    return Status::Error(StatusCode::InvalidParticleState,
                         "field-line segment has incomplete vertices");
  }
  const double lengthM = context.segment->GetLength();
  if (!std::isfinite(lengthM) || lengthM <= 0.0) {
    return Status::Error(StatusCode::InvalidParticleState,
                         "field-line segment length is invalid");
  }

  double* b0 = begin->GetDatum_ptr(FL::DatumAtVertexMagneticField);
  double* b1 = end->GetDatum_ptr(FL::DatumAtVertexMagneticField);
  if (!b0 || !b1) {
    return Status::Error(StatusCode::InvalidParticleState,
                         "magnetic field is absent from a segment vertex");
  }
  const double absB0 = Vector3D::Length(b0);
  const double absB1 = Vector3D::Length(b1);
  if (!std::isfinite(absB0) || !std::isfinite(absB1) ||
      absB0 <= 0.0 || absB1 <= 0.0) {
    return Status::Error(StatusCode::InvalidParticleState,
                         "magnetic-field magnitude must be finite and positive");
  }
  background->dLnAbsBdsPerM = (std::log(absB1) - std::log(absB0)) / lengthM;

  double u0[3] = {0.0, 0.0, 0.0};
  double u1[3] = {0.0, 0.0, 0.0};
  begin->GetPlasmaVelocity(u0);
  end->GetPlasmaVelocity(u1);
  const double uParallel0 = Vector3D::DotProduct(u0, b0) / absB0;
  const double uParallel1 = Vector3D::DotProduct(u1, b1) / absB1;
  const double fraction = context.state.coordinate -
                          std::floor(context.state.coordinate);
  background->plasmaAdvectionMPerS =
      (1.0 - fraction) * uParallel0 + fraction * uParallel1;
  background->parallelVelocityGradientPerS =
      (uParallel1 - uParallel0) / lengthM;

  double densityCurrent0 = 0.0;
  double densityCurrent1 = 0.0;
  double densityPrevious0 = 0.0;
  double densityPrevious1 = 0.0;
  begin->GetDatum(FL::DatumAtVertexPlasmaDensity, &densityCurrent0);
  end->GetDatum(FL::DatumAtVertexPlasmaDensity, &densityCurrent1);
  begin->GetDatum(FL::DatumAtVertexPrevious::DatumAtVertexPlasmaDensity,
                  &densityPrevious0);
  end->GetDatum(FL::DatumAtVertexPrevious::DatumAtVertexPlasmaDensity,
                &densityPrevious1);
  const double current = (1.0 - fraction) * densityCurrent0 +
                         fraction * densityCurrent1;
  const double previous = (1.0 - fraction) * densityPrevious0 +
                          fraction * densityPrevious1;
  if (!std::isfinite(current) || !std::isfinite(previous) ||
      current <= 0.0 || previous <= 0.0) {
    return Status::Error(StatusCode::InvalidParticleState,
                         "current and previous plasma density must be positive");
  }

  double intervalS = densityEvolutionIntervalS;
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
  if (AMPS2SWMF::MagneticFieldLineUpdate::SecondCouplingFlag) {
    intervalS = AMPS2SWMF::MagneticFieldLineUpdate::LastCouplingTime -
                AMPS2SWMF::MagneticFieldLineUpdate::LastLastCouplingTime;
  }
#endif
  if (!std::isfinite(intervalS) || intervalS <= 0.0) {
    return Status::Error(StatusCode::InvalidParticleState,
                         "background density epochs are not time ordered");
  }
  background->velocityDivergencePerS = -std::log(current / previous) / intervalS;
  if (!std::isfinite(background->velocityDivergencePerS)) {
    return Status::Error(StatusCode::InvalidParticleState,
                         "plasma velocity divergence is not finite");
  }

  const Background::BackgroundSnapshot& snapshot =
      Background::SnapshotStore::Instance().AcquireForMover();

  // The vertex density is a number density [m^-3].  Convert it with the same
  // mean-ion-mass convention used by the SWMF coupling before evaluating
  // v_A=|B|/sqrt(mu0*rho).  Keeping this calculation in the shared adapter
  // guarantees that both focused movers use the same plasma/wave frames.
  const double absB = SEP::FieldLineData::GetAbsB(
      context.state.coordinate, context.segment, context.state.fieldLineId);
  const double massDensity =
      current * PIC::CPLR::SWMF::MeanPlasmaAtomicMass;
  if (!std::isfinite(absB) || absB <= 0.0 ||
      !std::isfinite(massDensity) || massDensity <= 0.0) {
    return Status::Error(StatusCode::InvalidParticleState,
        "interpolated magnetic field and mass density must be positive");
  }
  background->alfvenSpeedMPerS =
      absB / std::sqrt(VacuumPermeability * massDensity);
  if (!std::isfinite(background->alfvenSpeedMPerS) ||
      background->alfvenSpeedMPerS < 0.0 ||
      background->alfvenSpeedMPerS >= SpeedOfLight) {
    return Status::Error(StatusCode::InvalidParticleState,
                         "local Alfven speed is invalid or superluminal");
  }

  // A mover may never integrate beyond the immutable background generation it
  // acquired.  The driver opens the read phase at the authoritative PIC time;
  // subtracting that time from valid_until gives the exact remaining interval.
  background->snapshotSecondsRemaining =
      snapshot.valid_until_seconds() - Background::SimulationTimeSeconds();
  if (!std::isfinite(background->snapshotSecondsRemaining) ||
      background->snapshotSecondsRemaining < 0.0) {
    return Status::Error(StatusCode::InvalidParticleState,
                         "background snapshot validity interval has expired");
  }
  std::ostringstream identity;
  identity << snapshot.configuration_fingerprint() << ":g"
           << snapshot.field_line_generation() << ":t"
           << snapshot.epoch_seconds();
  background->identity = identity.str();
  return Status::Ok();
}

}  // namespace PICAdapter
}  // namespace Transport
}  // namespace SEP

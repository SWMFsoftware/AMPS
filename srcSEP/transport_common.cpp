#include "transport_common.h"
#include "util/sep_reproducible_reduction.h"
#include "util/sep_runtime_contracts.h"

#include "amps2swmf.h"

#include <cmath>
#include <algorithm>
#include <mutex>
#include <sstream>
#include <vector>
#include <cstring>
#include <limits>

namespace SEP {
namespace Transport {
namespace PICAdapter {

std::uint64_t CampaignRandomSeed = 0;

namespace {

std::uint64_t HashParticleBirthState(const ParticleState& state) {
  // FNV-1a hashes an immutable physical tuple.  Particle-buffer addresses,
  // worker indices, and MPI ranks are excluded from the persistent identity.
  std::uint64_t hash = UINT64_C(1469598103934665603);
  const auto append = [&hash](const void* value, std::size_t size) {
    const unsigned char* bytes =
        reinterpret_cast<const unsigned char*>(value);
    for (std::size_t i = 0; i < size; ++i) {
      hash ^= bytes[i];
      hash *= UINT64_C(1099511628211);
    }
  };
  append(&state.species, sizeof(state.species));
  append(&state.fieldLineId, sizeof(state.fieldLineId));
  append(&state.coordinate, sizeof(state.coordinate));
  append(&state.vParallelMPerS, sizeof(state.vParallelMPerS));
  append(&state.vNormalMPerS, sizeof(state.vNormalMPerS));
  append(&state.massKg, sizeof(state.massKg));
  return hash == 0 ? UINT64_C(1) : hash;
}

Status ContextAtRelativeArcLength(const ParticleContext& base,
                                  double relativeArcLengthM,
                                  ParticleContext* sampled) {
  if (!sampled || !std::isfinite(relativeArcLengthM))
    return Status::Error(StatusCode::InvalidArgument,
                         "background sample displacement must be finite");
  *sampled = base;
  sampled->state.coordinate =
      PIC::FieldLine::FieldLinesAll[base.state.fieldLineId].move(
          base.state.coordinate, relativeArcLengthM, sampled->segment);
  sampled->segment =
      PIC::FieldLine::FieldLinesAll[base.state.fieldLineId].GetSegment(
          sampled->state.coordinate);
  return sampled->segment
      ? Status::Ok()
      : Status::Error(StatusCode::OutOfDomain,
                      "requested sample is outside the field line");
}

struct PendingQueue {
  std::vector<CouplingRecord> values;
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

bool PendingLess(const CouplingRecord& left,
                 const CouplingRecord& right) {
  if (left.turbulenceStateIdentity != right.turbulenceStateIdentity)
    return left.turbulenceStateIdentity < right.turbulenceStateIdentity;
  if (left.fieldLineId != right.fieldLineId)
    return left.fieldLineId < right.fieldLineId;
  if (left.stableParticleId != right.stableParticleId)
    return left.stableParticleId < right.stableParticleId;
  if (left.eventIndex != right.eventIndex)
    return left.eventIndex < right.eventIndex;
  if (left.intervalIndex != right.intervalIndex)
    return left.intervalIndex < right.intervalIndex;
  return static_cast<int>(left.eventType) < static_cast<int>(right.eventType);
}

}  // namespace

void InitializeParticleTransportState(long int pointer,
                                      std::uint64_t sourceKey,
                                      std::uint64_t sequenceKey,
                                      const double momentumKgMPerS[3]) {
  if (!momentumKgMPerS) return;
  std::uint64_t hash = UINT64_C(1469598103934665603);
  const auto append = [&hash](const void* value, std::size_t size) {
    const unsigned char* bytes =
        reinterpret_cast<const unsigned char*>(value);
    for (std::size_t i = 0; i < size; ++i) {
      hash ^= bytes[i];
      hash *= UINT64_C(1099511628211);
    }
  };
  append(&sourceKey, sizeof(sourceKey));
  append(&sequenceKey, sizeof(sequenceKey));
  append(momentumKgMPerS, 3 * sizeof(double));
  if (hash == 0) hash = UINT64_C(1);
  InitializeParticleTransportStateWithIdentity(pointer, hash);
}

void InitializeParticleTransportStateWithIdentity(
    long int pointer, std::uint64_t stableParticleId) {
  // Particle-buffer addresses are deliberately excluded.  AMPS may compact or
  // migrate the buffer, whereas the caller-provided source-event identity is a
  // persistent physical label and therefore survives restart/decomposition.
  PIC::ParticleBuffer::byte* data =
      PIC::ParticleBuffer::GetParticleDataPointer(pointer);
  if (!data || stableParticleId == 0) return;
  *reinterpret_cast<std::uint64_t*>(data + SEP::Offset::TransportSchema) =
      UINT64_C(0x5352435345500001);
  *reinterpret_cast<std::uint64_t*>(data + SEP::Offset::StableParticleId) =
      stableParticleId;
  *reinterpret_cast<double*>(data + SEP::Offset::MfpOpticalDepth) =
      std::numeric_limits<double>::quiet_NaN();
  *reinterpret_cast<std::uint64_t*>(data + SEP::Offset::MfpEventIndex) = 0;
}

Status QueueWaveContribution(const CouplingRecord& value) {
  if (value.fieldLineId < 0 || value.species < 0 ||
      !(value.statisticalWeight > 0.0) ||
      value.stableParticleId == 0 || value.snapshotGeneration == 0 ||
      !(value.dtS > 0.0) || !std::isfinite(value.dtS) ||
      !std::isfinite(value.signedPathM)) {
    return Status::Error(StatusCode::InvalidArgument,
        "invalid self-contained particle-wave coupling record");
  }
  LocalQueue().values.push_back(value);
  return Status::Ok();
}

Status DrainWaveContributions(std::vector<CouplingRecord>* records) {
  if (!records)
    return Status::Error(StatusCode::InvalidArgument,
                         "coupling drain destination is null");
  std::vector<CouplingRecord> merged;
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

  // Exercise the same physical-key reducer used by decomposition tests before
  // touching legacy G arrays.  The reduced values are an audit ledger here;
  // path-resolved deposition below still decodes exact segment overlaps.
  std::vector<std::vector<Reproducibility::Contribution> > partitions(1);
  partitions[0].reserve(merged.size());
  for (const CouplingRecord& value : merged) {
    Reproducibility::Contribution audit;
    audit.key.schema = Reproducibility::ProductionContributionKeySchema;
    audit.key.source = static_cast<std::uint64_t>(value.eventType);
    audit.key.fieldLine = static_cast<std::uint64_t>(value.fieldLineId);
    audit.key.segment = static_cast<std::uint64_t>(
        std::max(0.0, std::floor(value.startCoordinate)));
    audit.key.branch = static_cast<std::uint64_t>(value.resonantBranch + 1);
    audit.key.species = static_cast<std::uint64_t>(value.species);
    audit.key.particle = value.stableParticleId;
    audit.key.step = value.snapshotGeneration;
    audit.key.event = value.eventIndex;
    audit.key.interval = value.intervalIndex;
    audit.key.purpose = value.pitchAngleResolved ? 2 : 1;
    audit.streaming = value.signedPathM * value.statisticalWeight;
    audit.resonantCount = value.eventType == CouplingEventType::ScatteringEvent;
    partitions[0].push_back(audit);
  }
  std::vector<Reproducibility::SegmentAccumulator> auditLedger;
  const Status reductionStatus = Reproducibility::CanonicalPartitionReduction(
      partitions, &auditLedger);
  if (!reductionStatus.ok()) return reductionStatus;

  std::sort(merged.begin(), merged.end(), PendingLess);
  records->swap(merged);
  return Status::Ok();
}

void FlushWaveContributions() {
  // Kept only for source compatibility with out-of-tree callers during the
  // WP43 migration.  The production drivers use DrainWaveContributions and
  // propagate its Status.  Discarding here is safer than resurrecting direct
  // mutation of legacy G+/G- arrays behind an untyped void interface.
  std::vector<CouplingRecord> discarded;
  (void)DrainWaveContributions(&discarded);
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
  const Background::BackgroundSnapshot& snapshot =
      Background::SnapshotStore::Instance().AcquireForMover();
  context->snapshotGeneration = snapshot.field_line_generation();
  context->particleStepEpochS = Background::SimulationTimeSeconds();
  const std::uint64_t transportSchema = UINT64_C(0x5352435345500001);
  std::uint64_t* storedSchema = reinterpret_cast<std::uint64_t*>(
      context->data + SEP::Offset::TransportSchema);
  std::uint64_t* storedIdentity = reinterpret_cast<std::uint64_t*>(
      context->data + SEP::Offset::StableParticleId);
  if (*storedSchema != transportSchema || *storedIdentity == 0) {
    // Newly injected AMPS extension bytes are normally zero.  The explicit
    // schema marker also detects pre-WP09 restart records and initializes their
    // event state without interpreting arbitrary legacy bytes as optical depth.
    *storedSchema = transportSchema;
    *storedIdentity = HashParticleBirthState(context->state);
    *reinterpret_cast<double*>(
        context->data + SEP::Offset::MfpOpticalDepth) =
            std::numeric_limits<double>::quiet_NaN();
    *reinterpret_cast<std::uint64_t*>(
        context->data + SEP::Offset::MfpEventIndex) = 0;
  }
  context->stableParticleId = *storedIdentity;
  context->statisticalWeight =
      PIC::ParticleWeightTimeStep::GlobalParticleWeight[context->state.species] *
      PIC::ParticleBuffer::GetIndividualStatWeightCorrection(pointer);
  if (!(context->statisticalWeight > 0.0) ||
      !std::isfinite(context->statisticalWeight)) {
    return Status::Error(StatusCode::InvalidParticleState,
                         "particle statistical weight is invalid");
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

Status ClipPathToFieldLineBoundary(int fieldLineId, double startCoordinate,
                                   double attemptedDisplacementM,
                                   double* boundaryCoordinate,
                                   double* inDomainDisplacementM) {
  if (!boundaryCoordinate || !inDomainDisplacementM ||
      fieldLineId < 0 || fieldLineId >= PIC::FieldLine::nFieldLine ||
      !std::isfinite(startCoordinate) ||
      !std::isfinite(attemptedDisplacementM) ||
      attemptedDisplacementM == 0.0) {
    return Status::Error(StatusCode::InvalidArgument,
                         "boundary clipping input is invalid");
  }
  PIC::FieldLine::cFieldLine* line =
      &PIC::FieldLine::FieldLinesAll[fieldLineId];
  const int count = line->GetTotalSegmentNumber();
  const int startSegment = static_cast<int>(std::floor(startCoordinate));
  const double fraction = startCoordinate - std::floor(startCoordinate);
  if (count <= 0 || startSegment < 0 || startSegment >= count)
    return Status::Error(StatusCode::OutOfDomain,
                         "boundary clipping start is outside the field line");

  double distanceM = 0.0;
  if (attemptedDisplacementM > 0.0) {
    distanceM += (1.0 - fraction) * line->GetSegment(startSegment)->GetLength();
    for (int i = startSegment + 1; i < count; ++i)
      distanceM += line->GetSegment(i)->GetLength();
    *boundaryCoordinate = static_cast<double>(count);
    *inDomainDisplacementM = distanceM;
  }
  else {
    distanceM += fraction * line->GetSegment(startSegment)->GetLength();
    for (int i = startSegment - 1; i >= 0; --i)
      distanceM += line->GetSegment(i)->GetLength();
    *boundaryCoordinate = 0.0;
    *inDomainDisplacementM = -distanceM;
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

Status EvaluateLocalBackgroundAt(const ParticleContext& context,
                                 double relativeArcLengthM,
                                 LocalBackgroundView* background) {
  namespace FL = PIC::FieldLine;
  if (!background || !context.segment || !std::isfinite(relativeArcLengthM)) {
    return Status::Error(StatusCode::InvalidArgument,
                         "local background requires a segment and finite location");
  }

  ParticleContext sampled;
  Status sampledStatus = ContextAtRelativeArcLength(
      context, relativeArcLengthM, &sampled);
  if (!sampledStatus.ok()) return sampledStatus;

  FL::cFieldLineVertex* begin = sampled.segment->GetBegin();
  FL::cFieldLineVertex* end = sampled.segment->GetEnd();
  if (!begin || !end) {
    return Status::Error(StatusCode::InvalidParticleState,
                         "field-line segment has incomplete vertices");
  }
  const double lengthM = sampled.segment->GetLength();
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
  const double fraction = sampled.state.coordinate -
                          std::floor(sampled.state.coordinate);
  for (int d = 0; d < 3; ++d) {
    background->magneticFieldT[d] =
        (1.0 - fraction) * b0[d] + fraction * b1[d];
    background->plasmaVelocityMPerS[d] =
        (1.0 - fraction) * u0[d] + fraction * u1[d];
  }
  sampled.segment->GetDir(background->tangent);
  const double tangentNorm = Vector3D::Length(background->tangent);
  background->magneticFieldMagnitudeT =
      Vector3D::Length(background->magneticFieldT);
  if (!(tangentNorm > 0.0) ||
      !(background->magneticFieldMagnitudeT > 0.0)) {
    return Status::Error(StatusCode::InvalidParticleState,
        "interpolated tangent and magnetic field must be nonzero");
  }
  for (int d = 0; d < 3; ++d) {
    background->tangent[d] /= tangentNorm;
    background->bUnit[d] = background->magneticFieldT[d] /
        background->magneticFieldMagnitudeT;
  }
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

  const Background::BackgroundSnapshot& snapshot =
      Background::SnapshotStore::Instance().AcquireForMover();
  const double intervalS = snapshot.physical_epoch_interval_seconds();
  if (intervalS > 0.0) {
    const double materialLogDensityDerivativePerS =
        std::log(current / previous) / intervalS;
    background->velocityDivergencePerS =
        -materialLogDensityDerivativePerS;
    const ScalarResult residual = RuntimeContracts::ContinuityResidualPerS(
        materialLogDensityDerivativePerS,
        background->velocityDivergencePerS);
    if (!residual.status.ok()) return residual.status;
    background->continuityResidualPerS = residual.value;
  }
  else if (std::fabs(current - previous) <=
           32.0 * std::numeric_limits<double>::epsilon() * current) {
    // A static provider may publish a single epoch.  Equal density states then
    // represent exactly zero temporal compression/expansion.
    background->velocityDivergencePerS = 0.0;
  }
  else {
    return Status::Error(StatusCode::InvalidParticleState,
        "distinct density states require ordered physical background epochs");
  }
  if (!std::isfinite(background->velocityDivergencePerS)) {
    return Status::Error(StatusCode::InvalidParticleState,
                         "plasma velocity divergence is not finite");
  }

  // The vertex density is a number density [m^-3].  Convert it with the same
  // mean-ion-mass convention used by the SWMF coupling before evaluating
  // v_A=|B|/sqrt(mu0*rho).  Keeping this calculation in the shared adapter
  // guarantees that both focused movers use the same plasma/wave frames.
  const double absB = SEP::FieldLineData::GetAbsB(
      sampled.state.coordinate, sampled.segment, sampled.state.fieldLineId);
  const double massDensity =
      current * PIC::CPLR::SWMF::MeanPlasmaAtomicMass;
  if (!std::isfinite(absB) || absB <= 0.0 ||
      !std::isfinite(massDensity) || massDensity <= 0.0) {
    return Status::Error(StatusCode::InvalidParticleState,
        "interpolated magnetic field and mass density must be positive");
  }
  background->alfvenSpeedMPerS =
      absB / std::sqrt(VacuumPermeability * massDensity);
  background->numberDensityPerM3 = current;
  background->massDensityKgPerM3 = massDensity;
  background->sampleCoordinate = sampled.state.coordinate;
  // Resolve bb:grad(U) separately from d(U.b)/ds.  The rank-one tensor below is
  // the derivative information available from a field-line chord:
  // grad(U) ~= (dU/ds) b.  Curvature db/ds is evaluated from the endpoint B
  // directions, so a rigid vector translation on a curved line yields zero
  // strain while d(U.b)/ds retains its geometric U.kappa contribution.
  RuntimeContracts::VelocityGradientInput derivativeInput;
  for (int i = 0; i < 3; ++i) {
    derivativeInput.velocityMPerS[i] = background->plasmaVelocityMPerS[i];
    derivativeInput.bUnit[i] = background->bUnit[i];
    derivativeInput.curvaturePerM[i] =
        (b1[i] / absB1 - b0[i] / absB0) / lengthM;
    const double dUds = (u1[i] - u0[i]) / lengthM;
    for (int j = 0; j < 3; ++j)
      derivativeInput.gradientPerS[i][j] =
          dUds * background->bUnit[j];
  }
  const RuntimeContracts::VelocityDerivatives velocityDerivatives =
      RuntimeContracts::ComputeVelocityDerivatives(derivativeInput);
  if (!velocityDerivatives.status.ok()) return velocityDerivatives.status;
  background->fieldAlignedStrainPerS =
      velocityDerivatives.fieldAlignedStrainPerS;
  background->parallelVelocityGradientPerS =
      velocityDerivatives.parallelVelocityGradientPerS;
  background->velocityDerivativeMethod = velocityDerivatives.method;
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
  background->generation = snapshot.field_line_generation();
  background->provenance = snapshot.provenance();
  return Status::Ok();
}

Status EvaluateLocalBackground(const ParticleContext& context,
                               LocalBackgroundView* background) {
  return EvaluateLocalBackgroundAt(context, 0.0, background);
}

}  // namespace PICAdapter
}  // namespace Transport
}  // namespace SEP

#include "transport_common.h"
#include "util/sep_provider_radius_frame.h"
#include "util/sep_reproducible_reduction.h"
#include "util/sep_background_derivatives.h"
#include "util/sep_analytic_radial_background.h"
#include "util/sep_population_control.h"
#include "util/sep_resonant_collision.h"
#include "util/sep_turbulence_core.h"
#include "adapters/swcme1d_adapter.h"
#include "adapters/reduced_shock_background_adapter.h"

#include "amps2swmf.h"

#include <cmath>
#include <algorithm>
#include <mutex>
#include <sstream>
#include <vector>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace SEP {
namespace Transport {
namespace PICAdapter {

std::uint64_t CampaignRandomSeed = 0;

namespace {RunPurpose gRunPurpose=RunPurpose::Production;}
void SetStandaloneRunPurpose(RunPurpose purpose) {gRunPurpose=purpose;}
RunPurpose CurrentRunPurpose() {return gRunPurpose;}
Status ProductionClosurePreflight() {
  if(gRunPurpose==RunPurpose::ControlledDiagnostic)return Status::Ok();
  return ResonantCollision::RequireImplementedProductionClosure(
      SEP::AlfvenTurbulence_Kolmogorov::ActiveFlag&&
          SEP::AlfvenTurbulence_Kolmogorov::ParticleCouplingMode,
      SEP::AlfvenTurbulence_Kolmogorov::ActiveFlag&&
          Turbulence::EvolvesLocally(Turbulence::ActiveConfiguration().source));
}

namespace {
bool gVelocityGradientFixtureActive=false;
std::array<double,9> gVelocityGradientFixture{};
bool gVelocityGradientFixtureTransverse=false;
double gVelocityGradientFixtureDivergence=0.0;
void AssertDerivativeFixtureMayWrite() {
  const auto current=Background::SnapshotStore::Instance().Current();
  if(current)Background::SnapshotStore::Instance().AssertProviderMayWrite(
      current->provider());
}
}

ScopedVelocityGradientFixture::ScopedVelocityGradientFixture(
    const std::array<double,9>& gradient)
    :previousActive_(gVelocityGradientFixtureActive),
     previousGradient_(gVelocityGradientFixture),
     previousTransverse_(gVelocityGradientFixtureTransverse),
     previousDivergence_(gVelocityGradientFixtureDivergence) {
  AssertDerivativeFixtureMayWrite();
  for(double value:gradient)if(!std::isfinite(value))
    throw std::invalid_argument("controlled velocity gradient must be finite");
  gVelocityGradientFixture=gradient;
  gVelocityGradientFixtureActive=true;
  gVelocityGradientFixtureTransverse=false;
}
ScopedVelocityGradientFixture::ScopedVelocityGradientFixture(double divergence)
    :ScopedVelocityGradientFixture(std::array<double,9>{}) {
  if(!std::isfinite(divergence))
    throw std::invalid_argument("transverse fixture divergence must be finite");
  gVelocityGradientFixtureTransverse=true;
  gVelocityGradientFixtureDivergence=divergence;
}
ScopedVelocityGradientFixture::~ScopedVelocityGradientFixture() {
  AssertDerivativeFixtureMayWrite();
  gVelocityGradientFixture=previousGradient_;
  gVelocityGradientFixtureActive=previousActive_;
  gVelocityGradientFixtureTransverse=previousTransverse_;
  gVelocityGradientFixtureDivergence=previousDivergence_;
}

namespace {

const std::uint64_t kOperatorIds[3] = {7, 8, 9};

int OperatorIndex(StochasticOperator operation) {
  const int index = static_cast<int>(operation);
  return index >= 0 && index < 3 ? index : -1;
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
  PIC::ParticleBuffer::byte* data =
      PIC::ParticleBuffer::GetParticleDataPointer(pointer);
  if (!data || !momentumKgMPerS) return;
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
  InitializeParticleTransportStateWithIdentity(pointer,hash);
}

void InitializeParticleTransportStateWithIdentity(long int pointer,
                                                  std::uint64_t stableId) {
  if(NativeDriverWindowActive())exit(__LINE__,__FILE__,
      "fresh births inside an active source window must declare their source origin and remaining duration explicitly");
  Birth::Record birth;
  Status status=Birth::InitialCondition(Background::SimulationTimeSeconds(),&birth);
  if(status.ok())status=InitializeParticleTransportStateWithIdentity(pointer,stableId,birth);
  if(!status.ok())exit(__LINE__,__FILE__,status.message.c_str());
}

Status InitializeParticleTransportStateWithIdentity(long int pointer,
    std::uint64_t stableId,const Birth::Record& birth) {
  const Status valid=Birth::ValidateRecord(birth);if(!valid.ok())return valid;
  if(birth.origin==Birth::Origin::InitialCondition) {
    if(NativeDriverWindowActive())return Status::Error(StatusCode::InvalidParticleState,
        "a fresh initial condition cannot replace a scheduled birth inside an active driver window");
  } else {
    if(!birth.pendingFirstMove)return Status::Error(StatusCode::InvalidParticleState,
        "fresh source initialization cannot invent an already consumed movement history");
    Birth::DriverWindow window;const auto windowStatus=GetNativeDriverWindow(&window);
    if(!windowStatus.ok())return windowStatus;
    Birth::FirstMovePlan plan;const auto planned=Birth::PlanFirstMove(birth,window,window.durationS,&plan);
    if(!planned.ok())return planned;
  }
  if(pointer<0||stableId==0)return Status::Error(StatusCode::InvalidArgument,
      "fresh transport birth requires a valid particle and nonzero source identity");
  PIC::ParticleBuffer::byte* data=
      PIC::ParticleBuffer::GetParticleDataPointer(pointer);
  if(!data)return Status::Error(StatusCode::InvalidParticleState,"fresh transport birth storage is absent");
  // Validate all birth input before writing any extension bytes. Allocation is
  // owned by the caller; this function never draws source phase/velocity or
  // converts an undeclared frame on the caller's behalf.
  *reinterpret_cast<std::uint64_t*>(data + SEP::Offset::TransportSchema) =
      ParticleTransportSchema;
  *reinterpret_cast<std::uint64_t*>(data + SEP::Offset::StableParticleId) =
      stableId;
  *reinterpret_cast<double*>(data + SEP::Offset::MfpOpticalDepth) =
      std::numeric_limits<double>::quiet_NaN();
  *reinterpret_cast<std::uint64_t*>(data + SEP::Offset::MfpEventIndex) = 0;
  *reinterpret_cast<std::uint64_t*>(data + SEP::Offset::SplitEventCounter) = 0;
  // This is the only fresh-birth boundary. Every operator receives one keyed
  // stream now, before motion; loading a later particle never guesses birth
  // state from its current coordinate or silently reconstructs lost history.
  OperatorRandomState randomStates[3];
  for (int operation = 0; operation < 3; ++operation) {
    KeyedRandomStream random(CampaignRandomSeed, stableId, kOperatorIds[operation], 0);
    randomStates[operation].stream = random.SaveState();
  }
  std::memcpy(data + SEP::Offset::OperatorRandomStates,
              randomStates, sizeof(randomStates));
  std::memcpy(data + SEP::Offset::BirthRecord,&birth,sizeof(birth));
  return Status::Ok();
}

Status ResumeRandomStream(const ParticleContext& context,
                          StochasticOperator operation,
                          KeyedRandomStream* random,
                          std::uint64_t* acceptedIntervals) {
  const int index = OperatorIndex(operation);
  if (index < 0)
    return Status::Error(StatusCode::InvalidArgument,
                          "unknown production stochastic operator");
  return ResumeOperatorRandomStream(context.randomStates[index], random,
                                     acceptedIntervals);
}

Status StageRandomStream(ParticleContext* context,
                         StochasticOperator operation,
                         const KeyedRandomStream& random,
                         std::uint64_t acceptedIntervals) {
  const int index = OperatorIndex(operation);
  if (!context || index < 0)
    return Status::Error(StatusCode::InvalidArgument,
                          "invalid production RNG stage target");
  return StageOperatorRandomStream(&context->randomStates[index], random,
                                    acceptedIntervals);
}

namespace {
// Population control is serial after the particle worker join. Keep its
// reservation state separate from stepping queues. The collective preflight
// detects hash collisions rather than silently assigning a rank-dependent ID.
std::unordered_set<std::uint64_t> gPopulationIdentities;
bool gPopulationPrepared=false;
long int gPreparedSplitParent=-1;
std::uint64_t gPreparedSplitChild=0;

std::uint64_t& ParticleWord(long int particle,int offset) {
  return *reinterpret_cast<std::uint64_t*>(
      PIC::ParticleBuffer::GetParticleDataPointer(particle)+offset);
}
std::uint64_t CheckedPopulationIdentity(long int particle) {
  if(ParticleWord(particle,SEP::Offset::TransportSchema)!=ParticleTransportSchema ||
     ParticleWord(particle,SEP::Offset::StableParticleId)==0)
    exit(__LINE__,__FILE__,"population control requires initialized corrected-model particle history");
  Birth::Record birth;std::memcpy(&birth,PIC::ParticleBuffer::GetParticleDataPointer(particle)+
      Offset::BirthRecord,sizeof(birth));
  const auto status=Birth::ValidateRecord(birth);
  if(!status.ok())exit(__LINE__,__FILE__,status.message.c_str());
  return ParticleWord(particle,SEP::Offset::StableParticleId);
}
std::uint64_t NextSplitIdentity(long int parent) {
  const std::uint64_t counter=ParticleWord(parent,SEP::Offset::SplitEventCounter);
  if(counter==std::numeric_limits<std::uint64_t>::max())
    exit(__LINE__,__FILE__,"exact split lineage counter exhausted before mutation");
  return PopulationControl::DeriveChildId(CheckedPopulationIdentity(parent),
      counter+1,1,ParticleTransportSchema);
}
void PrepareExactSplit(long int parent) {
  if(!gPopulationPrepared || gPreparedSplitParent!=-1)
    exit(__LINE__,__FILE__,"exact split lacks a collective identity preflight");
  const std::uint64_t child=NextSplitIdentity(parent);
  if(gPopulationIdentities.count(child)!=1)
    exit(__LINE__,__FILE__,"exact split child identity was not reserved");
  // Validate every stream before the core changes weights or allocates a child.
  // This prevents a corrupt checkpoint from being hidden by fresh child state.
  OperatorRandomState states[3];
  std::memcpy(states,PIC::ParticleBuffer::GetParticleDataPointer(parent)+
      SEP::Offset::OperatorRandomStates,sizeof(states));
  for(int operation=0;operation<3;++operation) {
    KeyedRandomStream random(CampaignRandomSeed,CheckedPopulationIdentity(parent),
                             kOperatorIds[operation],0);
    std::uint64_t acceptedIntervals=0;
    const Status status=ResumeOperatorRandomStream(states[operation],&random,
                                                    &acceptedIntervals);
    if(!status.ok())exit(__LINE__,__FILE__,status.message.c_str());
  }
  gPreparedSplitParent=parent;
  gPreparedSplitChild=child;
}
void FinishExactSplit(long int parent,long int child) {
  // Prepare has already validated this no-fail committed path. The original
  // parent's phase space, hazard, event count, RNG and cached normal are retained.
  // Only its split serial increments. The child starts an independent keyed
  // future: conditional survival of a Poisson optical-depth clock is exponential
  // (memorylessness), so drawing a child residual preserves the one-particle law.
  // Cloned phase space remains correlated at birth by the existing exact splitter.
  if(gPreparedSplitParent!=parent || gPreparedSplitChild==0)
    exit(__LINE__,__FILE__,"exact split lifecycle callbacks are out of order");
  ++ParticleWord(parent,SEP::Offset::SplitEventCounter);
  ParticleWord(child,SEP::Offset::StableParticleId)=gPreparedSplitChild;
  ParticleWord(child,SEP::Offset::SplitEventCounter)=0;
  ParticleWord(child,SEP::Offset::MfpEventIndex)=0;
  PIC::ParticleBuffer::byte* data=PIC::ParticleBuffer::GetParticleDataPointer(child);
  *reinterpret_cast<double*>(data+SEP::Offset::MfpOpticalDepth)=
      std::numeric_limits<double>::quiet_NaN();
  OperatorRandomState states[3];
  for(int operation=0;operation<3;++operation) {
    KeyedRandomStream random(CampaignRandomSeed,gPreparedSplitChild,kOperatorIds[operation],0);
    states[operation].stream=random.SaveState();
  }
  std::memcpy(data+SEP::Offset::OperatorRandomStates,states,sizeof(states));
  gPreparedSplitParent=-1;
  gPreparedSplitChild=0;
}
} // namespace

void RegisterPopulationControlHooks() {
  namespace Split=PIC::ParticleSplitting::FledLine;
  Split::ParticleIdentity=CheckedPopulationIdentity;
  Split::PrepareExactSplit=PrepareExactSplit;
  Split::PostExactSplit=FinishExactSplit;
}

void BeginPopulationControl() {
  namespace FL=PIC::FieldLine;
  namespace PB=PIC::ParticleBuffer;
  gPopulationPrepared=false;
  gPreparedSplitParent=-1;
  gPreparedSplitChild=0;
  std::vector<std::uint64_t> local;
  for(int line=0;line<FL::nFieldLine;++line) {
    if(!FL::FieldLinesAll[line].IsInitialized())continue;
    for(FL::cFieldLineSegment* segment=FL::FieldLinesAll[line].GetFirstSegment();
        segment;segment=segment->GetNext()) {
      for(long int particle=segment->FirstParticleIndex;particle!=-1;particle=PB::GetNext(particle)) {
        local.push_back(CheckedPopulationIdentity(particle));
        // Each initial candidate is split at most once in the core controller.
        // Reserve all potential children, including on other MPI owners, before
        // mutation. Exhausted counters are rejected only if selected for split.
        if(ParticleWord(particle,SEP::Offset::SplitEventCounter)!=
           std::numeric_limits<std::uint64_t>::max())
          local.push_back(NextSplitIdentity(particle));
      }
    }
  }
  if(local.size()>static_cast<std::size_t>(std::numeric_limits<int>::max()))
    exit(__LINE__,__FILE__,"population identity preflight exceeds MPI count range");
  const int localCount=static_cast<int>(local.size());
  int ranks=0;
  MPI_Comm_size(MPI_GLOBAL_COMMUNICATOR,&ranks);
  std::vector<int> counts(ranks),displacements(ranks);
  MPI_Allgather(&localCount,1,MPI_INT,counts.data(),1,MPI_INT,MPI_GLOBAL_COMMUNICATOR);
  std::size_t total=0;
  for(int rank=0;rank<ranks;++rank) {
    if(counts[rank]<0 || total>static_cast<std::size_t>(std::numeric_limits<int>::max()-counts[rank]))
      exit(__LINE__,__FILE__,"global population identity preflight exceeds MPI count range");
    displacements[rank]=static_cast<int>(total);
    total+=static_cast<std::size_t>(counts[rank]);
  }
  std::vector<std::uint64_t> all(total);
  MPI_Allgatherv(local.data(),localCount,MPI_UINT64_T,all.data(),counts.data(),
                 displacements.data(),MPI_UINT64_T,MPI_GLOBAL_COMMUNICATOR);
  gPopulationIdentities.clear();
  for(std::uint64_t identity:all)
    if(identity==0 || !gPopulationIdentities.insert(identity).second)
      exit(__LINE__,__FILE__,"duplicate current or prospective split identity before population mutation");
  gPopulationPrepared=true;
}

void QueueWaveContribution(const CouplingRecord& value) {
  if (value.fieldLineId < 0 || value.species < 0 ||
      !(value.statisticalWeight > 0.0) ||
      value.stableParticleId == 0 || value.snapshotGeneration == 0 ||
      !(value.dtS > 0.0) || !std::isfinite(value.dtS) ||
      !std::isfinite(value.signedPathM)) {
    exit(__LINE__, __FILE__,
         "invalid self-contained particle-wave coupling record");
  }
  LocalQueue().values.push_back(value);
}

bool PendingWaveContributionsEmpty() {
  std::lock_guard<std::mutex> lock(gQueueRegistrationMutex);
  for(const PendingQueue* queue:gRegisteredQueues)
    if(queue&&!queue->values.empty())return false;
  return true;
}

void FlushWaveContributions() {
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
  if (!reductionStatus.ok())
    exit(__LINE__, __FILE__, reductionStatus.message.c_str());

  std::sort(merged.begin(), merged.end(), PendingLess);
  for (const CouplingRecord& value : merged) {
    if (value.pitchAngleResolved) {
      SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::
          AccumulateParticleFluxForWaveCoupling(
              value.fieldLineId, value.species, value.statisticalWeight,
              value.dtS, value.midpointParallelVelocityMPerS,
              value.midpointNormalVelocityMPerS,
              value.startCoordinate, value.finishCoordinate,
              value.signedPathM);
    }
    else {
      SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::
          AccumulateParticleFluxForWaveCoupling(
              value.fieldLineId, value.species, value.statisticalWeight,
              value.dtS, std::hypot(value.midpointParallelVelocityMPerS,
                                    value.midpointNormalVelocityMPerS),
              value.startCoordinate,
              value.finishCoordinate, value.signedPathM);
    }
  }
}

Status LoadParticle(long int pointer, ParticleContext* context) {
  if (!context || pointer < 0) {
    return Status::Error(StatusCode::InvalidArgument,
                         "particle context and handle must be valid");
  }
  const Status closure=ProductionClosurePreflight();
  if(!closure.ok())return closure;
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
  const auto driverSnapshot=ValidateNativeDriverSnapshot(snapshot);
  if(!driverSnapshot.ok())return driverSnapshot;
  context->snapshotGeneration = snapshot.field_line_generation();
  context->particleStepEpochS = Background::SimulationTimeSeconds();
  context->snapshotFromS=snapshot.valid_from_seconds();
  context->snapshotUntilS=snapshot.valid_until_seconds();
  context->stepPrepared=false;context->consumesPendingBirth=false;context->usesDriverWindow=false;
  std::uint64_t* storedSchema = reinterpret_cast<std::uint64_t*>(
      context->data + SEP::Offset::TransportSchema);
  std::uint64_t* storedIdentity = reinterpret_cast<std::uint64_t*>(
      context->data + SEP::Offset::StableParticleId);
  if (*storedSchema != ParticleTransportSchema || *storedIdentity == 0)
    return Status::Error(StatusCode::InvalidParticleState,
        "SEP particle has absent/incompatible transport schema or identity: "
        "fresh sources must call InitializeParticleTransportState; legacy "
        "records omit complete RNG history and cannot be resumed implicitly");
  context->stableParticleId = *storedIdentity;
  std::memcpy(&context->birth,context->data+SEP::Offset::BirthRecord,sizeof(context->birth));
  const Status birthStatus=Birth::ValidateRecord(context->birth);
  if(!birthStatus.ok())return birthStatus;
  std::memcpy(context->randomStates,
      context->data + SEP::Offset::OperatorRandomStates,
      sizeof(context->randomStates));
  for (int operation = 0; operation < 3; ++operation) {
    KeyedRandomStream expected(CampaignRandomSeed, context->stableParticleId,
                               kOperatorIds[operation], 0);
    const Status randomStatus = expected.RestoreState(
        context->randomStates[operation].stream);
    if (!randomStatus.ok()) return randomStatus;
  }
  context->mfpOpticalDepth = *reinterpret_cast<double*>(
      context->data + SEP::Offset::MfpOpticalDepth);
  context->mfpEventIndex = *reinterpret_cast<std::uint64_t*>(
      context->data + SEP::Offset::MfpEventIndex);
  if (!std::isnan(context->mfpOpticalDepth) &&
      (!(context->mfpOpticalDepth > 0.0) ||
       !std::isfinite(context->mfpOpticalDepth)))
    return Status::Error(StatusCode::InvalidParticleState,
                          "stored MFP optical depth is invalid");
  if (std::isnan(context->mfpOpticalDepth) &&
      (context->mfpEventIndex != 0 ||
       context->randomStates[2].acceptedIntervals != 0))
    return Status::Error(StatusCode::InvalidParticleState,
        "MFP history is missing its residual optical depth; only an unused "
        "fresh scattering process may have an undrawn hazard");
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

Status PrepareParticleStep(ParticleContext* context,double requestedHostDtS,
                            double* effectiveDtS) {
  if(!context||!effectiveDtS||!std::isfinite(requestedHostDtS)||requestedHostDtS<0.0)
    return Status::Error(StatusCode::InvalidArgument,"particle step requires a finite nonnegative host duration and output");
  const auto valid=Birth::ValidateRecord(context->birth);if(!valid.ok())return valid;
  ParticleContext candidate=*context;
  candidate.consumesPendingBirth=context->birth.pendingFirstMove!=0;
  candidate.usesDriverWindow=NativeDriverWindowActive();
  if(candidate.consumesPendingBirth||candidate.usesDriverWindow) {
    Birth::DriverWindow window;auto status=GetNativeDriverWindow(&window);if(!status.ok())return status;
    if(requestedHostDtS!=window.durationS)
      return Status::Error(StatusCode::UnsupportedConfiguration,
          "mover host duration differs from the frozen common source clock; unequal species timesteps are unsupported");
    if(candidate.consumesPendingBirth) {
      status=Birth::PlanFirstMove(candidate.birth,window,requestedHostDtS,&candidate.firstMove);
      if(!status.ok())return status;
    } else {candidate.firstMove.parent=window;candidate.firstMove.durationS=requestedHostDtS;}
    const auto coverage=Birth::SnapshotRemaining(candidate.firstMove,0.0,
        candidate.snapshotFromS,candidate.snapshotUntilS);
    if(!coverage.status.ok())return coverage.status;
    candidate.moveDurationS=candidate.firstMove.durationS;
    // Display/legacy metadata only; validity arithmetic below keeps the exact
    // remainder separately and never recovers it from this rounded epoch.
    candidate.particleStepEpochS=window.endEpochS-candidate.moveDurationS;
  } else {
    if(CurrentRunPurpose()!=RunPurpose::ControlledDiagnostic)
      return Status::Error(StatusCode::UnsupportedConfiguration,
          "production particle step requires the explicit driver window before PIC::TimeStep");
    candidate.moveDurationS=requestedHostDtS;
    if(!std::isfinite(candidate.particleStepEpochS)||candidate.snapshotFromS>candidate.particleStepEpochS||
        std::isnan(candidate.snapshotUntilS)||candidate.snapshotUntilS<candidate.particleStepEpochS||
        candidate.snapshotUntilS-candidate.particleStepEpochS<requestedHostDtS)
      return Status::Error(StatusCode::InvalidParticleState,"direct diagnostic interval is outside its declared immutable snapshot");
  }
  candidate.stepPrepared=true;*effectiveDtS=candidate.moveDurationS;*context=candidate;
  return Status::Ok();
}

ScalarResult ParticleSnapshotRemaining(const ParticleContext& context,
    double elapsedS,std::uint64_t sampledGeneration) {
  ScalarResult result;
  if(!context.stepPrepared||sampledGeneration!=context.snapshotGeneration||
      !std::isfinite(elapsedS)||elapsedS<0.0||elapsedS>context.moveDurationS) {
    result.status=Status::Error(StatusCode::InvalidParticleState,
        "particle time plan is absent, stale, or outside its declared duration");return result;
  }
  if(context.usesDriverWindow||context.consumesPendingBirth)
    return Birth::SnapshotRemaining(context.firstMove,elapsedS,context.snapshotFromS,context.snapshotUntilS);
  result.value=(context.snapshotUntilS-context.particleStepEpochS)-elapsedS;
  if(std::isnan(result.value)||result.value<0.0) {
    result.status=Status::Error(StatusCode::StepUnderflow,"particle reached the immutable snapshot boundary");return result;
  }
  result.status=Status::Ok();return result;
}

ScalarResult SelectParticleSubstep(const ParticleContext& context,double elapsedS,
    const std::vector<StepLimit>& limits,double minimumAdaptiveStepS,StepDiagnostics* diagnostics) {
  if(context.consumesPendingBirth)
    return Birth::SelectFirstMoveSubstep(context.firstMove,elapsedS,limits,minimumAdaptiveStepS,diagnostics);
  return SelectSubstep(context.moveDurationS-elapsedS,limits,minimumAdaptiveStepS,diagnostics);
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
  if(!context.stepPrepared)return Status::Error(StatusCode::InvalidParticleState,
      "cannot commit particle motion without an explicit time plan");
  // The field-line coordinate is the sole position authority. Generic PIC
  // consumers also read PB::X, so derive its Cartesian view from that exact
  // coordinate before publishing ANY phase, random, birth or list state.
  // Validate bounds before GetCartesian: the core accessor otherwise invokes
  // the host fatal-error path for a missing segment instead of returning a
  // typed, atomic rejection. Geometry is immutable during particle reads.
  if(!PIC::FieldLine::FieldLinesAll||context.state.fieldLineId<0||
      context.state.fieldLineId>=PIC::FieldLine::nFieldLine)
    return Status::Error(StatusCode::InvalidParticleState,"cannot commit an invalid field-line position owner");
  auto& line=PIC::FieldLine::FieldLinesAll[context.state.fieldLineId];
  const int count=line.GetTotalSegmentNumber();
  if(context.state.coordinate<0.0||context.state.coordinate>=count||
      line.GetSegment(context.state.coordinate)!=context.segment||
      !context.segment->GetBegin()||!context.segment->GetEnd())
    return Status::Error(StatusCode::InvalidParticleState,"cannot derive Cartesian position from an incompatible field-line coordinate/segment");
  double begin[3],end[3],position[3];
  context.segment->GetBegin()->GetX(begin);context.segment->GetEnd()->GetX(end);
  for(unsigned component=0;component<3;++component)
    if(!std::isfinite(begin[component])||!std::isfinite(end[component]))
      return Status::Error(StatusCode::InvalidParticleState,"cannot commit a nonfinite field-line Cartesian geometry view");
  line.GetCartesian(position,context.state.coordinate);
  for(double component:position)if(!std::isfinite(component))
    return Status::Error(StatusCode::InvalidParticleState,"derived field-line Cartesian position is nonfinite");
  Birth::Record committedBirth=context.birth;
  if(context.consumesPendingBirth) {
    const auto consumed=Birth::StageConsumedFirstMove(context.birth,context.firstMove,&committedBirth);
    if(!consumed.ok())return consumed;
  }

  PIC::ParticleBuffer::SetVParallel(context.state.vParallelMPerS, context.data);
  PIC::ParticleBuffer::SetVNormal(context.state.vNormalMPerS, context.data);
  PIC::ParticleBuffer::SetFieldLineCoord(context.state.coordinate, context.data);
  PIC::ParticleBuffer::SetX(position,context.data);
  // The context contains private accepted stochastic history. Publish it only
  // with the matching phase-space commit, never partway through a mover call.
  std::memcpy(context.data + SEP::Offset::OperatorRandomStates,
              context.randomStates, sizeof(context.randomStates));
  *reinterpret_cast<double*>(context.data + SEP::Offset::MfpOpticalDepth) =
      context.mfpOpticalDepth;
  *reinterpret_cast<std::uint64_t*>(context.data + SEP::Offset::MfpEventIndex) =
      context.mfpEventIndex;
  std::memcpy(context.data+SEP::Offset::BirthRecord,&committedBirth,sizeof(committedBirth));

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
  sampled.segment->GetDir(background->tangent);
  const auto advection0=FieldLineAdvectionMPerS(u0,background->tangent);
  const auto advection1=FieldLineAdvectionMPerS(u1,background->tangent);
  if(!advection0.status.ok())return advection0.status;
  if(!advection1.status.ok())return advection1.status;
  const double uParallel0 = advection0.value;
  const double uParallel1 = advection1.value;
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
  const auto polarity=MagneticPolarityAlongFieldLine(
      background->magneticFieldT,background->tangent);
  if(!polarity.status.ok())return polarity.status;
  background->magneticPolarityAlongLine=polarity.value>0?1:-1;
  background->plasmaAdvectionMPerS =
      (1.0 - fraction) * uParallel0 + fraction * uParallel1;
  background->parallelVelocityGradientPerS =
      (uParallel1 - uParallel0) / lengthM;

  double densityCurrent0 = 0.0;
  double densityCurrent1 = 0.0;
  begin->GetDatum(FL::DatumAtVertexPlasmaDensity, &densityCurrent0);
  end->GetDatum(FL::DatumAtVertexPlasmaDensity, &densityCurrent1);
  const double current = (1.0 - fraction) * densityCurrent0 +
                         fraction * densityCurrent1;
  if (!std::isfinite(current) || current <= 0.0) {
    return Status::Error(StatusCode::InvalidParticleState,
                         "current plasma density must be positive");
  }

  const Background::BackgroundSnapshot& snapshot =
      Background::SnapshotStore::Instance().AcquireForMover();
  double positionM[3];
  const double* x0=begin->GetX();
  const double* x1=end->GetX();
  for(int d=0;d<3;++d)positionM[d]=(1.0-fraction)*x0[d]+fraction*x1[d];
  const double radiusM=std::hypot(positionM[0],std::hypot(positionM[1],positionM[2]));
  // Establish the provider's radial frame BEFORE selecting a derivative
  // implementation. Controlled tensor fixtures replace grad(U), not the
  // radial origin; they must neither lose a known analytic radius nor acquire
  // new endpoint-U validation. Conversely, importing coordinates/grad(U)
  // alone does not declare the Sun/origin for radius-dependent coefficients.
  background->radiusFrameIdentity=Background::RadiusFrame::ContractIdentity(snapshot.provider());
  background->heliocentricRadiusKnown=false;background->heliocentricRadiusM=0.0;
  for(unsigned d=0;d<3;++d) {
    background->radiusOriginM[d]=0.0;background->heliocentricPositionM[d]=0.0;
  }
  if(snapshot.provider()==Background::Provider::Analytic) {
    const auto declaration=Background::AnalyticRadial::DeclaredNativeWind();
    double relative[3];bool finite=true;
    for(unsigned d=0;d<3;++d) {
      background->radiusOriginM[d]=declaration.originM[d];
      relative[d]=positionM[d]-declaration.originM[d];
      finite=finite&&std::isfinite(relative[d])&&std::isfinite(declaration.originM[d]);
    }
    const double declaredRadius=std::hypot(relative[0],std::hypot(relative[1],relative[2]));
    if(finite&&std::isfinite(declaredRadius)&&declaredRadius>0.0) {
      background->heliocentricRadiusM=declaredRadius;background->heliocentricRadiusKnown=true;
      for(unsigned d=0;d<3;++d)background->heliocentricPositionM[d]=relative[d];
    }
  } else if(snapshot.provider()==Background::Provider::Swcme||snapshot.provider()==Background::Provider::ReducedShock) {
    if(std::isfinite(radiusM)&&radiusM>0.0) {
      background->heliocentricRadiusM=radiusM;background->heliocentricRadiusKnown=true;
      for(unsigned d=0;d<3;++d)background->heliocentricPositionM[d]=positionM[d];
    }
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
  } else if(snapshot.provider()==Background::Provider::Swmf) {
    // The PT producer explicitly declares HGR/HGI, whose origin is the Sun,
    // and its actual SI coordinate unit. This is provider metadata, not an
    // origin inferred from the numerical coordinate array or a radial wind.
    const auto& frame=AMPS2SWMF::FieldLineData::CoordinateFrame();
    background->radiusFrameIdentity="provider-radius-v2:SWMF-declared-heliocentric-"+frame;
    if(std::isfinite(radiusM)&&radiusM>0.0) {
      background->heliocentricRadiusKnown=true;background->heliocentricRadiusM=radiusM;
      for(unsigned d=0;d<3;++d)background->heliocentricPositionM[d]=positionM[d];
    }
#endif
  }
  std::array<double,9> gradient;
  double providerMassDensity=0.0;
  Status derivativeStatus=Status::Ok();
  if(gVelocityGradientFixtureActive) {
    gradient=gVelocityGradientFixture;
    if(gVelocityGradientFixtureTransverse) {
      if(Vector3D::Length(u0)!=0.0||Vector3D::Length(u1)!=0.0)
        return Status::Error(StatusCode::UnsupportedConfiguration,
            "transverse tube fixture requires zero centerline velocity");
      for(int i=0;i<3;++i)for(int j=0;j<3;++j)
        gradient[3*i+j]=0.5*gVelocityGradientFixtureDivergence*
            ((i==j?1.0:0.0)-background->tangent[i]*background->tangent[j]);
    }
  } else if(snapshot.provider()==Background::Provider::Analytic) {
    // Use the same initialization-owned radial declaration as generated
    // native endpoints. A translated Parker origin or configured speed must
    // affect both U and its full tensor; global-zero/400 km/s is not an
    // alternative gradient model. Arbitrary imported/native test velocities
    // require their explicit provider or the controlled derivative fixture.
    if(SEP::DomainType!=SEP::DomainType_ParkerSpiral&&SEP::DomainType!=SEP::DomainType_StraitLine)
      return Status::Error(StatusCode::UnsupportedConfiguration,
          "analytic radial derivative requires a generated analytic domain; imported velocity needs its own derivative contract");
    const auto declaration=Background::AnalyticRadial::DeclaredNativeWind();
    const double* positions[2]={x0,x1};
    const double* velocities[2]={u0,u1};
    for(int endpoint=0;endpoint<2;++endpoint) {
      const auto endpointStatus=Background::AnalyticRadial::ValidateVelocity(
          declaration,positions[endpoint],velocities[endpoint],2e-12);
      if(!endpointStatus.ok())return endpointStatus;
    }
    Background::AnalyticRadial::Sample radial;
    derivativeStatus=Background::AnalyticRadial::Evaluate(declaration,positionM,&radial);
    if(derivativeStatus.ok())gradient=radial.gradientPerS;
  } else if(snapshot.provider()==Background::Provider::Swcme) {
    const auto query=SW1DAdapter::QueryAtRadius(radiusM,false);
    if(query.code!=SW1DAdapter::StatusCode::Success||
        query.epoch_seconds!=snapshot.epoch_seconds())
      return Status::Error(StatusCode::UnresolvedCoefficient,
          "SWCME derivative query is unavailable or differs from frozen snapshot epoch: "+query.detail);
    const double radialDerivative=query.sample.divergence_s_inv()-
        2.0*query.sample.speed_m_s()/radiusM;
    derivativeStatus=RadialVelocityGradient(positionM,
        query.sample.speed_m_s(),radialDerivative,&gradient);
    providerMassDensity=current*query.sample.mass_density_kg_m3()/
        query.sample.number_density_m3();
    if(!std::isfinite(providerMassDensity)||providerMassDensity<=0.0)
      return Status::Error(StatusCode::InvalidParticleState,
          "SWCME query lacks canonical composition mass density");
  } else if(snapshot.provider()==Background::Provider::ReducedShock) {
    ReducedShock::AmbientSample sample;
    std::string error;
    if(ReducedShock::Metadata().epochS!=snapshot.epoch_seconds()||
        ReducedShock::Metadata().generation!=snapshot.field_line_generation()||
        !ReducedShock::EvaluateAmbient({{positionM[0],positionM[1],positionM[2]}},
            &sample,&error)||!sample.derivativesValid)
      return Status::Error(StatusCode::UnresolvedCoefficient,
          "reduced provider gradient is unavailable or differs from frozen snapshot: "+error);
    gradient=sample.velocityGradientPerS;
    providerMassDensity=current*sample.massDensityKgM3/sample.numberDensityM3;
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
  } else if(snapshot.provider()==Background::Provider::Swmf) {
    // Received three-dimensional U(x) samples, not a temporal density proxy,
    // determine these numerical derivatives. The receive phase fits an affine
    // Cartesian field on full-rank transverse stencils using pivoted QR and
    // rejects underdetermined data. Interpolate only that frozen epoch here.
    std::string error;
    if(!AMPS2SWMF::FieldLineData::VelocityGradientAt(sampled.state.fieldLineId,
          static_cast<int>(std::floor(sampled.state.coordinate)),fraction,
          snapshot.epoch_seconds(),&gradient,&error))
      return Status::Error(StatusCode::UnresolvedCoefficient,error);
    providerMassDensity=current*AMPS2SWMF::FieldLineData::DensityMassUnitKg();
#endif
  } else {
    // Fixed/Lagrangian vertex labels, full grad U and density source terms
    // are not specified by the legacy SWMF line import. Density history alone
    // cannot determine div U or the independent contraction bb:grad U.
    return Status::Error(StatusCode::UnsupportedConfiguration,
        "imported background requires an explicit Cartesian velocity-gradient contract; density history does not determine div U and bb:grad U");
  }
  if(!derivativeStatus.ok())return derivativeStatus;
  VelocityDerivatives derivatives;
  derivativeStatus=ContractVelocityGradient(gradient,background->bUnit,&derivatives);
  if(!derivativeStatus.ok())return derivativeStatus;
  background->velocityGradientPerS=derivatives.gradientPerS;
  background->velocityDivergencePerS=derivatives.divergencePerS;
  background->fieldAlignedStrainPerS=derivatives.fieldAlignedStrainPerS;

  // Convert native number density with the authority's frozen composition.
  // SWCME/reduced providers expose rho/ne explicitly; the analytic legacy
  // authority retains its declared mean-plasma-atomic-mass convention.
  const double absB = SEP::FieldLineData::GetAbsB(
      sampled.state.coordinate, sampled.segment, sampled.state.fieldLineId);
  const double massDensity = providerMassDensity>0.0?providerMassDensity:
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
  if (!std::isfinite(background->alfvenSpeedMPerS) ||
      background->alfvenSpeedMPerS < 0.0 ||
      background->alfvenSpeedMPerS >= SpeedOfLight) {
    return Status::Error(StatusCode::InvalidParticleState,
                         "local Alfven speed is invalid or superluminal");
  }

  // A mover may never integrate beyond its immutable background generation.
  // The explicit host plan supplies departure time independently of whether
  // this PIC build advances the mutable clock at step start or step end. Birth
  // remainders are retained separately from rounded absolute display epochs.
  // A static (analytic) or open-ended (SWMF) snapshot is valid until +inf
  // (SEP::Background::PrepareSnapshotForParticleStep()); the remaining time is
  // then +inf, which the movers and ComposeSnapshotValidityLimit() accept.
  // Rejecting every non-finite value stopped each analytic-background run at
  // its first particle.  Only NaN and an expired interval are errors.
  const auto remaining=ParticleSnapshotRemaining(context,0.0,snapshot.field_line_generation());
  if(!remaining.status.ok())return remaining.status;
  background->snapshotSecondsRemaining=remaining.value;
  if (std::isnan(background->snapshotSecondsRemaining) ||
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

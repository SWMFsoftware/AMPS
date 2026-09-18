#include "amps_particle_adapter.h"

#include "../SEP3D.h"
#include "../adapters/source_runtime.h"

#include <algorithm>
#include <cmath>
#include <cstring>

namespace SEP3D {
namespace AMPS {
namespace Movers {
namespace {

// Stored as raw bytes through memcpy: AMPS does not promise that an extension
// offset is naturally aligned. The schema tag makes stale checkpoints fail
// visibly instead of interpreting an older byte layout as valid state.
constexpr std::uint64_t kParticleSchema = UINT64_C(0x5345503344413031);
struct PersistentState {
  std::uint64_t schema = kParticleSchema;
  std::uint64_t stableId = 0;
  std::uint64_t completedStep = 0;
  std::uint64_t substep = 0;
  std::uint64_t lastShockGeneration = 0;
  double momentumKgMPerS = 0.0;
  double mu = 0.0;
  double gyrophaseRad = 0.0;
};

long int gParticleStateOffset = -1;
bool gStorageRequested = false;
Context gContext;
bool gContextInstalled = false;

Core::Status Invalid(const char* message) {
  return Core::Status(Core::StatusCode::InvalidInput, message);
}

void LoadPersistent(const PIC::ParticleBuffer::byte* data,
                    PersistentState* state) {
  std::memcpy(state, data + gParticleStateOffset, sizeof(*state));
}

void StorePersistent(PIC::ParticleBuffer::byte* data,
                     const PersistentState& state) {
  std::memcpy(data + gParticleStateOffset, &state, sizeof(state));
}

// Reconstruct a physical Cartesian velocity from the gyrotropic variables.
// The perpendicular basis is deterministic: choose the Cartesian axis least
// aligned with b, then apply cross products. This avoids pole singularities
// and makes restart output independent of platform-specific branch noise.
Core::Vec3 GyrotropicVelocity(double speedMPerS, double mu,
                              double gyrophaseRad,
                              const Core::Vec3& bHat) {
  const double ax = std::fabs(bHat.x);
  const double ay = std::fabs(bHat.y);
  const double az = std::fabs(bHat.z);
  const Core::Vec3 reference = ax <= ay && ax <= az ? Core::Vec3(1.0, 0.0, 0.0)
      : (ay <= az ? Core::Vec3(0.0, 1.0, 0.0)
                  : Core::Vec3(0.0, 0.0, 1.0));
  const Core::Vec3 e1 = bHat.Cross(reference).Normalized();
  const Core::Vec3 e2 = bHat.Cross(e1);
  const double perpendicular = std::sqrt(std::max(0.0, 1.0 - mu * mu));
  return speedMPerS * (mu * bHat + perpendicular *
      (std::cos(gyrophaseRad) * e1 + std::sin(gyrophaseRad) * e2));
}

void RecordOutcome(const Adapters::MoverResult& moved, int species,
                   std::uint64_t step) {
  if (gContext.ledger == nullptr) return;
  // Begin/Close are host orchestration operations. The mover records only one
  // outcome, so a failed row lookup cannot recursively change the disposition.
  (void)gContext.ledger->RecordMover(
      step, species, moved.disposition, moved.shockIntersection.crossed);
}

}  // namespace

Core::Status InstallContext(const Context& context) {
  if (context.resolveLocal == nullptr)
    return Invalid("AMPS mover context requires a local-state resolver");
  if (context.maximumSubsteps == 0)
    return Invalid("AMPS mover context requires a positive substep cap");
  if (SEP3D::ApplicationRuntime().state() ==
      RuntimeModel::LifecycleState::Running)
    return Core::Status(Core::StatusCode::InvalidTransition,
                        "cannot replace mover context during PIC::TimeStep");
  gContext = context;
  gContextInstalled = true;
  return Core::Status::OK();
}

bool ContextInstalled() { return gContextInstalled; }

Core::Status UpdateShock(const Adapters::ExpandingSphericalShock& shock) {
  if (!gContextInstalled)
    return Invalid("AMPS mover context is not installed");
  if (SEP3D::ApplicationRuntime().state() ==
      RuntimeModel::LifecycleState::Running)
    return Core::Status(Core::StatusCode::InvalidTransition,
                        "cannot replace shock state during particle motion");
  if (shock.active && (shock.generation == 0 ||
      !std::isfinite(shock.radiusAtStepStartM) ||
      shock.radiusAtStepStartM <= 0.0 ||
      !std::isfinite(shock.radialSpeedMPerS)))
    return Invalid("active AMPS shock state is invalid");
  gContext.shock = shock;
  return Core::Status::OK();
}

Core::Status RequestParticleStorage() {
  if (gStorageRequested) return Core::Status::OK();
  long int offset = -1;
  PIC::ParticleBuffer::RequestDataStorage(offset,
                                          sizeof(PersistentState));
  if (offset < 0)
    return Core::Status(Core::StatusCode::LayoutMismatch,
                        "AMPS rejected srcSEP3D particle storage");
  gParticleStateOffset = offset;
  gStorageRequested = true;
  return Core::Status::OK();
}

long int ParticleStateOffset() { return gParticleStateOffset; }

Core::Status InitializeParticle(long int ptr,
                                const Adapters::ParticleRecord& particle) {
  if (!gStorageRequested || ptr < 0 || particle.stableId == 0 ||
      particle.species < 0 || !std::isfinite(particle.momentumKgMPerS) ||
      particle.momentumKgMPerS < 0.0 || !std::isfinite(particle.mu) ||
      particle.mu < -1.0 || particle.mu > 1.0)
    return Invalid("new AMPS particle or persistent state is invalid");
  PIC::ParticleBuffer::byte* data =
      PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  if (data == nullptr) return Invalid("AMPS particle pointer is null");
  PersistentState state;
  state.stableId = particle.stableId;
  state.completedStep = particle.completedStep;
  state.substep = particle.substep;
  state.lastShockGeneration = particle.lastShockGeneration;
  state.momentumKgMPerS = particle.momentumKgMPerS;
  state.mu = particle.mu;
  state.gyrophaseRad = particle.gyrophaseRad;
  StorePersistent(data, state);
  return Core::Status::OK();
}

Core::Status ReadParticle(long int ptr, Adapters::ParticleRecord* particle) {
  if (!gStorageRequested || ptr < 0 || particle == nullptr)
    return Invalid("AMPS particle read request is invalid");
  PIC::ParticleBuffer::byte* data =
      PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  if (data == nullptr) return Invalid("AMPS particle pointer is null");
  PersistentState state;
  LoadPersistent(data, &state);
  if (state.schema != kParticleSchema || state.stableId == 0)
    return Invalid("AMPS particle has no valid srcSEP3D persistent state");
  double x[3]; PIC::ParticleBuffer::GetX(x, data);
  particle->stableId = state.stableId;
  particle->species = PIC::ParticleBuffer::GetI(data);
  particle->positionM = Core::Vec3(x);
  particle->momentumKgMPerS = state.momentumKgMPerS;
  particle->mu = state.mu;
  particle->gyrophaseRad = state.gyrophaseRad;
  particle->statisticalWeight =
      PIC::ParticleWeightTimeStep::GlobalParticleWeight[particle->species] *
      PIC::ParticleBuffer::GetIndividualStatWeightCorrection(data);
  particle->completedStep = state.completedStep;
  particle->substep = state.substep;
  particle->lastShockGeneration = state.lastShockGeneration;
  return Core::Status::OK();
}

InjectionOutcome InjectParticles(const Adapters::InjectionPlan& plan) {
  InjectionOutcome outcome;
  if (!gContextInstalled || !plan.status.ok()) {
    outcome.status = !plan.status.ok()
        ? plan.status
        : Invalid("AMPS mover context must be installed before injection");
    outcome.rejected = plan.particles.size();
    return outcome;
  }
  for (const Adapters::InjectedParticle& injected : plan.particles) {
    const Adapters::ParticleRecord& particle = injected.particle;
    double x[3]; particle.positionM.CopyTo(x);
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
        PIC::Mesh::mesh->findTreeNode(x);
    if (node != nullptr && node->Thread != PIC::ThisThread) continue;
    if (!injected.status.ok() || node == nullptr || node->block == nullptr ||
        particle.species < 0 || particle.species >= PIC::nTotalSpecies) {
      ++outcome.rejected;
      continue;
    }
    Adapters::LocalTransportRecord local;
    const Core::Status localStatus = gContext.resolveLocal(
        particle.positionM, particle.species, particle.momentumKgMPerS,
        particle.mu, node, &local);
    if (!localStatus.ok()) { ++outcome.rejected; continue; }
    const double speed = Transport::RelativisticSpeed(
        particle.momentumKgMPerS,
        PIC::MolecularData::GetMass(particle.species));
    const Core::Vec3 velocity = GyrotropicVelocity(
        speed, particle.mu, particle.gyrophaseRad, local.background.bHat);
    double v[3]; velocity.CopyTo(v);
    int species = particle.species;
    double correction = particle.statisticalWeight /
        PIC::ParticleWeightTimeStep::GlobalParticleWeight[species];
    const long int ptr = PIC::ParticleBuffer::InitiateParticle(
        x, v, &correction, &species, nullptr,
        _PIC_INIT_PARTICLE_MODE__ADD2LIST_, static_cast<void*>(node));
    if (ptr < 0) { ++outcome.rejected; continue; }
    const Core::Status initialized = InitializeParticle(ptr, particle);
    if (!initialized.ok()) {
      PIC::ParticleBuffer::DeleteParticle(ptr);
      ++outcome.rejected;
      continue;
    }
    ++outcome.allocated;
  }
  outcome.status = outcome.rejected == 0
      ? Core::Status::OK()
      : Core::Status(Core::StatusCode::Error,
                     "one or more planned source particles were rejected by AMPS");
  return outcome;
}

int MoveParticle(long int ptr, double dtTotal,
                 cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  if (!gStorageRequested || !gContextInstalled || ptr < 0 ||
      startNode == nullptr || !std::isfinite(dtTotal) || dtTotal <= 0.0) {
    if (ptr >= 0) PIC::ParticleBuffer::DeleteParticle(ptr);
    return _PARTICLE_LEFT_THE_DOMAIN_;
  }
  PIC::ParticleBuffer::byte* data =
      PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  if (data == nullptr) return _PARTICLE_LEFT_THE_DOMAIN_;

  PersistentState persistent;
  LoadPersistent(data, &persistent);
  const int species = PIC::ParticleBuffer::GetI(data);
  const auto& runtime = SEP3D::ApplicationRuntime();
  const auto& configuration = runtime.configuration();
  if (persistent.schema != kParticleSchema || persistent.stableId == 0 ||
      species < 0 || species >= PIC::nTotalSpecies || !configuration ||
      runtime.state() != RuntimeModel::LifecycleState::Running) {
    PIC::ParticleBuffer::DeleteParticle(ptr);
    return _PARTICLE_LEFT_THE_DOMAIN_;
  }

  double position[3];
  PIC::ParticleBuffer::GetX(position, data);
  Adapters::MoverInput input;
  input.model = configuration->options().transport;
  input.particle.stableId = persistent.stableId;
  input.particle.species = species;
  input.particle.positionM = Core::Vec3(position);
  input.particle.momentumKgMPerS = persistent.momentumKgMPerS;
  input.particle.mu = persistent.mu;
  input.particle.gyrophaseRad = persistent.gyrophaseRad;
  input.particle.completedStep = persistent.completedStep;
  input.particle.substep = persistent.substep;
  input.particle.lastShockGeneration = persistent.lastShockGeneration;
  input.particle.statisticalWeight =
      PIC::ParticleWeightTimeStep::GlobalParticleWeight[species] *
      PIC::ParticleBuffer::GetIndividualStatWeightCorrection(data);
  input.shock = gContext.shock;
  input.speciesMassKg = PIC::MolecularData::GetMass(species);
  input.speciesChargeC = PIC::MolecularData::GetElectricCharge(species);
  input.requestedDtS = dtTotal;
  input.innerRadiusM = configuration->options().innerRadiusM;
  input.outerRadiusM = configuration->options().outerRadiusM;
  input.campaignSeed = configuration->options().campaignSeed;
  input.pitchScheme = configuration->options().pitchAngleScheme;
  input.perpendicularDiffusion =
      configuration->options().perpendicularDiffusion;
  input.constantKappaPerpendicularM2PerS =
      configuration->options().constantKappaPerpendicularM2PerS;
  input.kappaPerpendicularToParallelRatio =
      configuration->options().kappaPerpendicularToParallelRatio;
  input.drift = configuration->options().drift;
  input.timeStepControls.cellCrossingFraction =
      configuration->options().cellCrossingFraction;
  input.timeStepControls.diffusionFraction =
      configuration->options().diffusionFraction;
  input.timeStepControls.focusingFraction =
      configuration->options().focusingFraction;
  input.timeStepControls.coolingFraction =
      configuration->options().coolingFraction;
  input.timeStepControls.fieldVariationFraction =
      configuration->options().fieldVariationFraction;
  input.timeStepControls.shockCrossingFraction =
      configuration->options().shockCrossingFraction;
  input.timeStepControls.minimumSubstepS =
      configuration->options().minimumTransportSubstepS;

  struct ResolverBridge {
    LocalRecordResolver resolver = nullptr;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node = nullptr;
  } bridge{gContext.resolveLocal, startNode};
  auto resolveEverySubstep = [](
      const Adapters::ParticleRecord& particle, double elapsedTimeS,
      void* opaque, Adapters::LocalTransportRecord* local) {
    ResolverBridge* state = static_cast<ResolverBridge*>(opaque);
    if (state == nullptr || state->resolver == nullptr || local == nullptr)
      return Invalid("AMPS local-state resolver bridge is invalid");
    double x[3];
    particle.positionM.CopyTo(x);
    state->node = PIC::Mesh::mesh->findTreeNode(x, state->node);
    if (state->node == nullptr || state->node->block == nullptr)
      return Core::Status(Core::StatusCode::DomainExit,
                          "particle left the allocated AMR tree during subcycling");
    Core::Status resolved = state->resolver(
        particle.positionM, particle.species, particle.momentumKgMPerS,
        particle.mu, state->node, local);
    if (resolved.ok() && local->timeToSnapshotBoundaryS > 0.0)
      local->timeToSnapshotBoundaryS = std::max(
          0.0, local->timeToSnapshotBoundaryS - elapsedTimeS);
    return resolved;
  };
  Adapters::RequestedTimeAdvance request;
  request.input = input;
  request.resolveLocal = resolveEverySubstep;
  request.resolverContext = &bridge;
  request.maximumSubsteps = gContext.maximumSubsteps;
  Adapters::MoverResult moved =
      Adapters::AdvanceParticleRequestedTime(request);
  const std::uint64_t step = runtime.counters().completedSteps;
  RecordOutcome(moved, species, step);

  if (moved.disposition != Adapters::ParticleDisposition::Active ||
      !moved.status.ok()) {
    PIC::ParticleBuffer::DeleteParticle(ptr);
    return _PARTICLE_LEFT_THE_DOMAIN_;
  }

  // completedSteps is the start tick while Runtime is Running.  The accepted
  // particle has consumed the complete requested interval and therefore uses
  // the next integer tick for all future semantic random keys.
  persistent.completedStep = step + 1;
  persistent.substep = moved.particle.substep;
  persistent.lastShockGeneration = moved.particle.lastShockGeneration;
  persistent.momentumKgMPerS = moved.particle.momentumKgMPerS;
  persistent.mu = moved.particle.mu;
  persistent.gyrophaseRad = moved.particle.gyrophaseRad;
  StorePersistent(data, persistent);
  double finalPosition[3];
  moved.particle.positionM.CopyTo(finalPosition);
  PIC::ParticleBuffer::SetX(finalPosition, data);
  const double speed = Transport::RelativisticSpeed(
      persistent.momentumKgMPerS, input.speciesMassKg);
  Core::Vec3 velocity = GyrotropicVelocity(
      speed, persistent.mu, persistent.gyrophaseRad,
      moved.finalLocal.background.bHat);
  double finalVelocity[3];
  velocity.CopyTo(finalVelocity);
  PIC::ParticleBuffer::SetV(finalVelocity, data);

  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* finalNode =
      PIC::Mesh::mesh->findTreeNode(finalPosition, bridge.node);
  int i = 0, j = 0, k = 0;
  if (finalNode == nullptr || finalNode->block == nullptr ||
      PIC::Mesh::mesh->FindCellIndex(
          finalPosition, i, j, k, finalNode, false) == -1) {
    PIC::ParticleBuffer::DeleteParticle(ptr);
    return _PARTICLE_LEFT_THE_DOMAIN_;
  }
#if _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
  long int* head = finalNode->block->tempParticleMovingListTable + i +
      _BLOCK_CELLS_X_ * (j + _BLOCK_CELLS_Y_ * k);
  const long int previousHead = *head;
  PIC::ParticleBuffer::SetNext(previousHead, data);
  PIC::ParticleBuffer::SetPrev(-1, data);
  if (previousHead != -1) PIC::ParticleBuffer::SetPrev(ptr, previousHead);
  *head = ptr;
#else
  PIC::ParticleBuffer::DeleteParticle(ptr);
  return _PARTICLE_LEFT_THE_DOMAIN_;
#endif
  return _PARTICLE_MOTION_FINISHED_;
}

}  // namespace Movers
}  // namespace AMPS
}  // namespace SEP3D

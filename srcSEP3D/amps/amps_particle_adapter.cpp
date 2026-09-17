#include "amps_particle_adapter.h"

#include "../SEP3D.h"

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
  if (SEP3D::ApplicationRuntime().state() ==
      RuntimeModel::LifecycleState::Running)
    return Core::Status(Core::StatusCode::InvalidTransition,
                        "cannot replace mover context during PIC::TimeStep");
  gContext = context;
  gContextInstalled = true;
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
  input.requestedDtS = dtTotal;
  input.innerRadiusM = configuration->options().innerRadiusM;
  input.outerRadiusM = configuration->options().outerRadiusM;
  input.campaignSeed = configuration->options().campaignSeed;
  input.pitchScheme = configuration->options().pitchAngleScheme;
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

  Core::Status localStatus = gContext.resolveLocal(
      input.particle.positionM, species, input.particle.momentumKgMPerS,
      input.particle.mu, startNode, &input.local);
  Adapters::MoverResult moved;
  if (!localStatus.ok()) {
    moved.status = localStatus;
    moved.particle = input.particle;
    moved.disposition = Adapters::ParticleDisposition::Failed;
  } else {
    moved = Adapters::AdvanceParticle(input);
  }
  const std::uint64_t step = runtime.counters().completedSteps;
  RecordOutcome(moved, species, step);

  if (moved.disposition != Adapters::ParticleDisposition::Active ||
      !moved.status.ok()) {
    PIC::ParticleBuffer::DeleteParticle(ptr);
    return _PARTICLE_LEFT_THE_DOMAIN_;
  }

  persistent.completedStep = step;
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
      input.local.background.bHat);
  double finalVelocity[3];
  velocity.CopyTo(finalVelocity);
  PIC::ParticleBuffer::SetV(finalVelocity, data);

  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* finalNode =
      PIC::Mesh::mesh->findTreeNode(finalPosition, startNode);
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

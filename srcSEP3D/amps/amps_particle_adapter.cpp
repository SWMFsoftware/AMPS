#include "amps_particle_adapter.h"

#include "../SEP3D.h"
#include "../adapters/source_runtime.h"
#include "../transport/population_control.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <unordered_set>
#include <vector>

namespace SEP3D {
namespace AMPS {
namespace Movers {
namespace {

// Stored as raw bytes through memcpy: AMPS does not promise that an extension
// offset is naturally aligned. The schema tag makes stale checkpoints fail
// visibly instead of interpreting an older byte layout as valid state.
constexpr std::uint64_t kParticleSchema = UINT64_C(0x5345503344413032);
struct PersistentState {
  std::uint64_t schema = kParticleSchema;
  std::uint64_t stableId = 0;
  std::uint64_t completedStep = 0;
  std::uint64_t substep = 0;
  std::uint64_t lastShockGeneration = 0;
  double momentumKgMPerS = 0.0;
  double mu = 0.0;
  double gyrophaseRad = 0.0;
  double remainingScatteringOpticalDepth =
      std::numeric_limits<double>::quiet_NaN();
  std::uint64_t nextScatteringEvent = 0;
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

struct PopulationParticle {
  long int ptr = -1;
  PersistentState persistent;
  Core::Vec3 positionM;
  Core::Vec3 momentumKgMPerS;
  double weightCorrection = 0.0;
};

Core::Vec3 VelocityFromMomentum(const Core::Vec3& momentumKgMPerS,
                                double massKg) {
  const double magnitude = momentumKgMPerS.Norm();
  if (magnitude == 0.0) return Core::Vec3();
  const double speed = Transport::RelativisticSpeed(magnitude, massKg);
  return momentumKgMPerS * (speed / magnitude);
}

Core::Status ReadPopulationParticle(long int ptr, double massKg,
                                    PopulationParticle* output) {
  if (output == nullptr || ptr < 0 || !std::isfinite(massKg) || massKg <= 0.0)
    return Invalid("population particle read request is invalid");
  PIC::ParticleBuffer::byte* data =
      PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  if (data == nullptr) return Invalid("population particle buffer is null");
  PopulationParticle candidate;
  candidate.ptr = ptr;
  LoadPersistent(data, &candidate.persistent);
  if (candidate.persistent.schema != kParticleSchema ||
      candidate.persistent.stableId == 0 ||
      !std::isfinite(candidate.persistent.momentumKgMPerS) ||
      candidate.persistent.momentumKgMPerS < 0.0) {
    return Invalid("population particle has invalid persistent state");
  }
  double x[3], v[3];
  PIC::ParticleBuffer::GetX(x, data);
  PIC::ParticleBuffer::GetV(v, data);
  candidate.positionM = Core::Vec3(x);
  const Core::Vec3 velocity(v);
  const double speed = velocity.Norm();
  candidate.momentumKgMPerS = speed > 0.0
      ? velocity * (candidate.persistent.momentumKgMPerS / speed)
      : Core::Vec3();
  candidate.weightCorrection =
      PIC::ParticleBuffer::GetIndividualStatWeightCorrection(data);
  if (!std::isfinite(candidate.weightCorrection) ||
      candidate.weightCorrection <= 0.0) {
    return Invalid("population particle weight correction is invalid");
  }
  *output = candidate;
  return Core::Status::OK();
}

Core::Status ResolvePopulationMagneticDirection(
    const PopulationParticle& particle, int species,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node, Core::Vec3* bHat) {
  (void)species;
  if (bHat == nullptr || gContext.resolveMagneticDirection == nullptr)
    return Invalid("population magnetic-direction request is invalid");
  const Core::Status resolved = gContext.resolveMagneticDirection(
      particle.positionM, node, bHat);
  if (!resolved.ok()) return resolved;
  if (!std::isfinite(bHat->x) || !std::isfinite(bHat->y) ||
      !std::isfinite(bHat->z) ||
      std::fabs(bHat->Norm() - 1.0) > 1.0e-12) {
    return Invalid("population control resolved an invalid magnetic direction");
  }
  return Core::Status::OK();
}

double GyrophaseFromMomentum(const Core::Vec3& momentum,
                             const Core::Vec3& bHat) {
  const double magnitude = momentum.Norm();
  if (magnitude == 0.0) return 0.0;
  const double ax = std::fabs(bHat.x);
  const double ay = std::fabs(bHat.y);
  const double az = std::fabs(bHat.z);
  const Core::Vec3 reference = ax <= ay && ax <= az
      ? Core::Vec3(1.0, 0.0, 0.0)
      : (ay <= az ? Core::Vec3(0.0, 1.0, 0.0)
                  : Core::Vec3(0.0, 0.0, 1.0));
  const Core::Vec3 e1 = bHat.Cross(reference).Normalized();
  const Core::Vec3 e2 = bHat.Cross(e1);
  const Core::Vec3 direction = momentum / magnitude;
  double phase = std::atan2(direction.Dot(e2), direction.Dot(e1));
  if (phase < 0.0) phase += 2.0 * Core::Const::kPi;
  return phase;
}

Core::Status WritePopulationParticle(
    const PopulationParticle& original, const Core::Vec3& positionM,
    const Core::Vec3& momentumKgMPerS, const Core::Vec3& bHat,
    double weightCorrection, double massKg, std::uint64_t stableId,
    std::uint64_t completedStep, std::uint64_t lastShockGeneration) {
  if (original.ptr < 0 || stableId == 0 ||
      !std::isfinite(weightCorrection) || weightCorrection <= 0.0)
    return Invalid("population particle write request is invalid");
  PIC::ParticleBuffer::byte* data =
      PIC::ParticleBuffer::GetParticleDataPointer(original.ptr);
  if (data == nullptr) return Invalid("population output buffer is null");
  PersistentState state;
  state.stableId = stableId;
  state.completedStep = completedStep;
  state.substep = 0;
  state.lastShockGeneration = lastShockGeneration;
  state.momentumKgMPerS = momentumKgMPerS.Norm();
  state.mu = state.momentumKgMPerS > 0.0
      ? std::max(-1.0, std::min(1.0,
          momentumKgMPerS.Dot(bHat) / state.momentumKgMPerS))
      : 0.0;
  state.gyrophaseRad = GyrophaseFromMomentum(momentumKgMPerS, bHat);
  // Population resampling creates a new stochastic history.  Drawing a fresh
  // optical depth under the new stable ID is unbiased; copying the parent's
  // residual would make split descendants collide at the same event time.
  state.remainingScatteringOpticalDepth =
      std::numeric_limits<double>::quiet_NaN();
  state.nextScatteringEvent = 0;
  StorePersistent(data, state);
  double x[3], v[3];
  positionM.CopyTo(x);
  VelocityFromMomentum(momentumKgMPerS, massKg).CopyTo(v);
  PIC::ParticleBuffer::SetX(x, data);
  PIC::ParticleBuffer::SetV(v, data);
  PIC::ParticleBuffer::SetIndividualStatWeightCorrection(
      weightCorrection, data);
  return Core::Status::OK();
}

std::uint64_t NewPopulationStableId(
    const PopulationControlRequest& request, std::uint64_t parentIdentity,
    std::uint64_t operationIndex, std::uint64_t childIndex,
    std::unordered_set<std::uint64_t>* used) {
  Transport::RandomKey key;
  key.campaignSeed = request.campaignSeed;
  key.particleId = parentIdentity;
  key.step = request.step;
  key.substep = operationIndex;
  key.purpose = Transport::RandomPurpose::PopulationSplitIdentity;
  std::uint64_t draw = childIndex;
  for (;;) {
    const std::uint64_t candidate =
        Transport::KeyedRandomStream::Hash(key, draw++);
    if (candidate != 0 && used->insert(candidate).second) return candidate;
  }
}

std::vector<PopulationParticle> CollectPopulationSpecies(
    long int firstParticle, int species, double massKg,
    Core::Status* status) {
  std::vector<PopulationParticle> particles;
  for (long int ptr = firstParticle; ptr != -1;
       ptr = PIC::ParticleBuffer::GetNext(ptr)) {
    PIC::ParticleBuffer::byte* data =
        PIC::ParticleBuffer::GetParticleDataPointer(ptr);
    if (PIC::ParticleBuffer::GetI(data) != species) continue;
    PopulationParticle particle;
    *status = ReadPopulationParticle(ptr, massKg, &particle);
    if (!status->ok()) {
      particles.clear();
      return particles;
    }
    particles.push_back(particle);
  }
  *status = Core::Status::OK();
  return particles;
}

double RelativeResidual(double observed, double reference) {
  return std::fabs(observed - reference) /
      std::max(std::fabs(reference), std::numeric_limits<double>::min());
}

}  // namespace

Core::Status InstallContext(const Context& context) {
  if (context.resolveLocal == nullptr ||
      context.resolveMagneticDirection == nullptr)
    return Invalid("AMPS mover context requires local-state and magnetic-direction resolvers");
  if (context.maximumSubsteps == 0)
    return Invalid("AMPS mover context requires a positive substep cap");
  if (SEP3D::ApplicationRuntime().state() ==
      RuntimeModel::LifecycleState::Running)
    return Core::Status(Core::StatusCode::InvalidTransition,
                        "cannot replace mover context during PIC::TimeStep");
  if (context.shock.active) {
    const auto geometry=Adapters::ValidateShockGeometry(context.shock);
    if (!geometry.ok()) return geometry;
  }
  gContext = context;
  gContextInstalled = true;
  return Core::Status::OK();
}

bool ContextInstalled() { return gContextInstalled; }

Core::Status UpdateShock(const Adapters::ExpandingShock& shock) {
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
  // Reject an unknown shape or malformed SSE axis before replacing live
  // context; failure leaves the previous published epoch intact.
  if (shock.active) {
    const auto geometry=Adapters::ValidateShockGeometry(shock);
    if (!geometry.ok()) return geometry;
  }
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
  state.remainingScatteringOpticalDepth =
      particle.remainingScatteringOpticalDepth;
  state.nextScatteringEvent = particle.nextScatteringEvent;
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
  particle->remainingScatteringOpticalDepth =
      state.remainingScatteringOpticalDepth;
  particle->nextScatteringEvent = state.nextScatteringEvent;
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

PopulationControlReport ApplyPopulationControl(
    const PopulationControlRequest& request) {
  PopulationControlReport report;
  if (!gStorageRequested || !gContextInstalled || request.step == 0 ||
      request.campaignSeed == 0 ||
      request.minimumParticlesPerCellPerSpecies < 2 ||
      request.targetParticlesPerCellPerSpecies <
          request.minimumParticlesPerCellPerSpecies ||
      request.maximumParticlesPerCellPerSpecies <
          request.targetParticlesPerCellPerSpecies) {
    report.status = Invalid("population-control request is invalid");
    return report;
  }

  std::unordered_set<std::uint64_t> usedStableIds;
  // Seed the collision guard with every owner-local active ID.  Source IDs are
  // already semantic hashes; retaining deleted IDs in this set also prevents
  // an operation later in the same boundary from recycling an identity.
  for (unsigned int blockIndex = 0;
       blockIndex < PIC::DomainBlockDecomposition::nLocalBlocks;
       ++blockIndex) {
    auto* node = PIC::DomainBlockDecomposition::BlockTable[blockIndex];
    if (node == nullptr || node->block == nullptr) continue;
    for (int cell = 0;
         cell < _BLOCK_CELLS_X_ * _BLOCK_CELLS_Y_ * _BLOCK_CELLS_Z_;
         ++cell) {
      for (long int ptr = node->block->FirstCellParticleTable[cell];
           ptr != -1; ptr = PIC::ParticleBuffer::GetNext(ptr)) {
        PersistentState state;
        LoadPersistent(PIC::ParticleBuffer::GetParticleDataPointer(ptr),
                       &state);
        if (state.schema != kParticleSchema || state.stableId == 0 ||
            !usedStableIds.insert(state.stableId).second) {
          report.status = Invalid(
              "population control found a missing or duplicate stable ID");
          return report;
        }
        ++report.particlesBefore;
      }
    }
  }

  for (unsigned int blockIndex = 0;
       blockIndex < PIC::DomainBlockDecomposition::nLocalBlocks;
       ++blockIndex) {
    auto* node = PIC::DomainBlockDecomposition::BlockTable[blockIndex];
    if (node == nullptr || node->block == nullptr) continue;
    for (int cell = 0;
         cell < _BLOCK_CELLS_X_ * _BLOCK_CELLS_Y_ * _BLOCK_CELLS_Z_;
         ++cell) {
      long int& first = node->block->FirstCellParticleTable[cell];
      for (int species = 0; species < PIC::nTotalSpecies; ++species) {
        const double massKg = PIC::MolecularData::GetMass(species);
        Core::Status status;
        std::vector<PopulationParticle> particles =
            CollectPopulationSpecies(first, species, massKg, &status);
        if (!status.ok()) {
          report.status = status;
          return report;
        }
        if (particles.empty()) continue;
        ++report.occupiedCellSpecies;

        // This ordinal is local to one physical cell/species population.  It
        // must not count operations in preceding owner-local blocks: those
        // blocks change when the same mesh is repartitioned, which would make
        // post-resampling IDs and merge directions depend on MPI ownership.
        // Within one population, candidate ordering is stable-ID sorted and
        // each count transition is deterministic, so this local ordinal is
        // invariant to rank count, block traversal, and unrelated cells.
        std::uint64_t operationIndex = 0;

        const bool mergeTriggered = particles.size() >
            request.maximumParticlesPerCellPerSpecies;
        while (mergeTriggered && particles.size() >
               request.targetParticlesPerCellPerSpecies) {
          if (particles.size() < 3) {
            report.status = Invalid(
                "relativistic 3-to-2 merge cannot meet the requested target");
            return report;
          }
          std::sort(particles.begin(), particles.end(),
              [](const PopulationParticle& left,
                 const PopulationParticle& right) {
                if (left.weightCorrection != right.weightCorrection)
                  return left.weightCorrection < right.weightCorrection;
                return left.persistent.stableId < right.persistent.stableId;
              });
          const PopulationParticle a = particles[0];
          const PopulationParticle b = particles[1];
          const PopulationParticle c = particles[2];
          const std::array<Transport::WeightedPhasePoint, 3> mergeInput = {{
              {a.weightCorrection, a.positionM, a.momentumKgMPerS},
              {b.weightCorrection, b.positionM, b.momentumKgMPerS},
              {c.weightCorrection, c.positionM, c.momentumKgMPerS}}};

          std::uint64_t ids[3] = {a.persistent.stableId,
                                  b.persistent.stableId,
                                  c.persistent.stableId};
          std::sort(ids, ids + 3);
          Transport::RandomKey identityKey;
          identityKey.campaignSeed = request.campaignSeed;
          identityKey.particleId = ids[0];
          identityKey.step = ids[1];
          identityKey.substep = ids[2];
          identityKey.purpose =
              Transport::RandomPurpose::PopulationMergeDirection;
          const std::uint64_t parentIdentity =
              Transport::KeyedRandomStream::Hash(identityKey, request.step);
          Transport::RandomKey directionKey;
          directionKey.campaignSeed = request.campaignSeed;
          directionKey.particleId = parentIdentity;
          directionKey.step = request.step;
          directionKey.substep = operationIndex;
          directionKey.purpose =
              Transport::RandomPurpose::PopulationMergeDirection;
          Transport::KeyedRandomStream directionRandom(directionKey);
          const double cosine = 1.0 - 2.0 * directionRandom.UniformOpen01();
          const double phi = 2.0 * Core::Const::kPi *
              directionRandom.UniformOpen01();
          const double sine = std::sqrt(std::max(0.0, 1.0 - cosine*cosine));
          const Core::Vec3 direction(sine * std::cos(phi),
                                     sine * std::sin(phi), cosine);

          const Transport::RelativisticMergeResult merged =
              Transport::MergeRelativisticThreeToTwo(
                  mergeInput, massKg, direction);
          if (!merged.status.ok()) {
            report.status = merged.status;
            return report;
          }
          // Both output particles are placed at the conserved weighted
          // centroid.  Resolve B at that actual output point before converting
          // their Cartesian momenta back to (mu, gyrophase); using B at one of
          // the three input positions would be inconsistent on an AMR cell
          // spanning a curved or strongly focusing field.
          PopulationParticle outputLocation = a;
          outputLocation.positionM = merged.outputPositionM;
          outputLocation.momentumKgMPerS = merged.firstMomentumKgMPerS;
          outputLocation.persistent.momentumKgMPerS =
              merged.firstMomentumKgMPerS.Norm();
          Core::Vec3 bHat;
          status = ResolvePopulationMagneticDirection(
              outputLocation, species, node, &bHat);
          if (!status.ok()) {
            report.status = status;
            return report;
          }
          const std::uint64_t lastShock = std::max(
              a.persistent.lastShockGeneration,
              std::max(b.persistent.lastShockGeneration,
                       c.persistent.lastShockGeneration));
          const std::uint64_t idA = NewPopulationStableId(
              request, parentIdentity, operationIndex, 0, &usedStableIds);
          const std::uint64_t idB = NewPopulationStableId(
              request, parentIdentity, operationIndex, 1, &usedStableIds);
          status = WritePopulationParticle(
              a, merged.outputPositionM, merged.firstMomentumKgMPerS, bHat,
              merged.outputWeight, massKg,
              idA, request.step, lastShock);
          if (status.ok()) status = WritePopulationParticle(
              b, merged.outputPositionM, merged.secondMomentumKgMPerS, bHat,
              merged.outputWeight, massKg,
              idB, request.step, lastShock);
          if (!status.ok()) {
            report.status = status;
            return report;
          }
          PIC::ParticleBuffer::DeleteParticle(c.ptr, first);

          report.maximumRelativeWeightResidual = std::max(
              report.maximumRelativeWeightResidual,
              merged.relativeWeightResidual);
          report.maximumRelativeMomentumResidual = std::max(
              report.maximumRelativeMomentumResidual,
              merged.relativeMomentumResidual);
          report.maximumRelativeEnergyResidual = std::max(
              report.maximumRelativeEnergyResidual,
              merged.relativeEnergyResidual);
          ++report.mergeOperations;
          ++operationIndex;
          particles = CollectPopulationSpecies(
              first, species, massKg, &status);
          if (!status.ok()) {
            report.status = status;
            return report;
          }
        }

        const bool splitTriggered = !particles.empty() &&
            particles.size() < request.minimumParticlesPerCellPerSpecies;
        while (splitTriggered && particles.size() <
               request.targetParticlesPerCellPerSpecies) {
          auto heaviest = std::max_element(
              particles.begin(), particles.end(),
              [](const PopulationParticle& left,
                 const PopulationParticle& right) {
                if (left.weightCorrection != right.weightCorrection)
                  return left.weightCorrection < right.weightCorrection;
                return left.persistent.stableId > right.persistent.stableId;
              });
          const PopulationParticle parent = *heaviest;
          Core::Vec3 bHat;
          status = ResolvePopulationMagneticDirection(
              parent, species, node, &bHat);
          if (!status.ok()) {
            report.status = status;
            return report;
          }
          const long int newPtr = PIC::ParticleBuffer::GetNewParticle(first);
          if (newPtr < 0) {
            report.status = Invalid(
                "AMPS particle buffer could not allocate a split child");
            return report;
          }
          PIC::ParticleBuffer::byte* childData =
              PIC::ParticleBuffer::GetParticleDataPointer(newPtr);
          PIC::ParticleBuffer::byte* parentData =
              PIC::ParticleBuffer::GetParticleDataPointer(parent.ptr);
          PIC::ParticleBuffer::CloneParticle(childData, parentData);
          PopulationParticle child = parent;
          child.ptr = newPtr;
          const double childWeight = 0.5 * parent.weightCorrection;
          const std::uint64_t childId = NewPopulationStableId(
              request, parent.persistent.stableId, operationIndex, 0,
              &usedStableIds);
          status = WritePopulationParticle(
              parent, parent.positionM, parent.momentumKgMPerS, bHat,
              childWeight, massKg, parent.persistent.stableId,
              request.step, parent.persistent.lastShockGeneration);
          if (status.ok()) status = WritePopulationParticle(
              child, child.positionM, child.momentumKgMPerS, bHat,
              childWeight, massKg, childId, request.step,
              parent.persistent.lastShockGeneration);
          if (!status.ok()) {
            PIC::ParticleBuffer::DeleteParticle(newPtr, first);
            report.status = status;
            return report;
          }
          report.maximumRelativeWeightResidual = std::max(
              report.maximumRelativeWeightResidual,
              RelativeResidual(2.0 * childWeight,
                               parent.weightCorrection));
          ++report.splitOperations;
          ++operationIndex;
          particles = CollectPopulationSpecies(
              first, species, massKg, &status);
          if (!status.ok()) {
            report.status = status;
            return report;
          }
        }
        report.particlesAfter += particles.size();
      }
    }
  }
  report.status = Core::Status::OK();
  return report;
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
  input.particle.remainingScatteringOpticalDepth =
      persistent.remainingScatteringOpticalDepth;
  input.particle.nextScatteringEvent = persistent.nextScatteringEvent;
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
  input.focusedScatteringFrame =
      configuration->options().focusedScatteringFrame;
  input.maximumScatteringEventsPerSubstep =
      configuration->options().maximumScatteringEventsPerSubstep;
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
  persistent.remainingScatteringOpticalDepth =
      moved.particle.remainingScatteringOpticalDepth;
  persistent.nextScatteringEvent = moved.particle.nextScatteringEvent;
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

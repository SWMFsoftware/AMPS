// ============================================================================
// srcSEP3D/main_lib.cpp
//
// AMPS application boundary through Phases R2, M, B, T, P, A, and O.
//
// The production boundary now owns the typed Runtime introduced in R2.  Both
// standalone and coupled hosts install a validated immutable configuration and
// use the same Runtime transitions; this file does not parse process arguments
// or parameter files.  Phase M builds the Cartesian AMR mesh and freezes AMPS
// storage offsets.  Phase B publishes a complete immutable ambient snapshot,
// and Phase T validates/fills the selected turbulence input.  Particle motion
// is dispatched through the single Phase-A AMPS adapter. Phase-O coordinators
// are available to the host only at joined step/checkpoint boundaries.
// ============================================================================

#include "SEP3D.h"

#include "background/bg_parker.h"
#include "background/background_snapshot.h"
#include "adapters/source_runtime.h"
#include "amps/amps_mover_status.h"
#include "amps/amps_particle_adapter.h"
#include "mesh/mesh_model.h"
#include "output/observer_runtime.h"
#include "output/publication.h"
#include "output/restart.h"
#include "runtime/runtime_adapters.h"
#include "turbulence/turbulence_models.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

[[noreturn]] void StopWithStatus(const char* operation,
                                 const SEP3D::Core::Status& status) {
  std::cerr << "[srcSEP3D] " << operation << " failed: "
            << status.message << '\n';
  std::abort();
}

std::shared_ptr<const SEP3D::Background::BackgroundSnapshot>
    gInstalledBackground;
std::shared_ptr<const SEP3D::Background::BackgroundSnapshot>
    gStagedBackground;
std::shared_ptr<SEP3D::Background::AnalyticParkerProvider>
    gAnalyticBackgroundProvider;
std::shared_ptr<SEP3D::Turbulence::TurbulenceProvider>
    gInstalledTurbulence;
std::shared_ptr<SEP3D::Turbulence::TurbulenceProvider>
    gStagedTurbulence;
std::shared_ptr<SEP3D::Adapters::ShockProvider> gInstalledShock;
SEP3D::Adapters::ParticleLedger gParticleLedger;
// Closed rows are global (MPI-reduced) evidence.  The mover-facing ledger
// above is cleared and reopened every particle phase; retaining the reduced
// history separately prevents rank migration from appearing as creation or
// loss and makes the same rows available to sampling and restart.
std::vector<SEP3D::Adapters::LedgerRow> gClosedParticleLedger;
std::vector<SEP3D::Adapters::SourceLedgerRow> gSourceLedger;
SEP3D::Output::ObserverRuntime gObserverRuntime;
std::unique_ptr<SEP3D::Output::RestartState> gPendingRestart;

int gStaticCellDataOffset = -1;
int gSamplingDataOffset = -1;
bool gStorageCallbacksRegistered = false;
std::unordered_map<PIC::Mesh::cDataCenterNode*, std::size_t> gCellSampleIndex;

const SEP3D::RuntimeModel::RunConfiguration3D& Configuration() {
  const auto& configuration = SEP3D::ApplicationRuntime().configuration();
  if (!configuration) {
    SEP3D::Core::Status status(
        SEP3D::Core::StatusCode::InvalidTransition,
        "the host must install RunConfiguration3D before mesh setup");
    StopWithStatus("configuration lookup", status);
  }
  return *configuration;
}

SEP3D::Mesh::ResolutionConfiguration ResolutionConfiguration() {
  const auto& options = Configuration().options();
  SEP3D::Mesh::ResolutionConfiguration result;
  result.innerRadiusM = options.innerRadiusM;
  result.outerRadiusM = options.outerRadiusM;
  result.minimumCellSizeM = options.minimumCellSizeM;
  result.backgroundCellSizeM = options.backgroundCellSizeM;
  result.enableRadialRefinement = options.enableRadialRefinement;
  result.solarSurfaceCellSizeM = options.solarSurfaceCellSizeM;
  result.solarRefinementOuterRadiusM =
      options.solarRefinementOuterRadiusM;
  result.solarRefinementProfile = options.solarRefinementProfile;
  result.solarRefinementExponent = options.solarRefinementExponent;
  result.enableTubeRefinement = options.enableTubeRefinement;
  result.tubeLongitudeRad = options.tubeLongitudeRad;
  result.tubeColatitudeRad = options.tubeColatitudeRad;
  result.tubeReferenceRadiusM = options.tubeReferenceRadiusM;
  result.tubeRadiusAtReferenceM = options.tubeRadiusAtReferenceM;
  result.tubeRadiusMode = options.tubeRadiusMode;
  result.tubeCellSizeM = options.tubeCellSizeM;
  result.tubeTransverseProfile = options.tubeTransverseProfile;
  result.tubeTransverseExponent = options.tubeTransverseExponent;
  result.solarWindSpeedMPerS = options.parker.solarWindSpeedMPerS;
  result.solarRotationRateRadPerS = options.parker.solarRotationRateRadPerS;
  result.cellsPerBlockEdge = options.meshCellsPerBlockEdge;
  result.maximumLevel = options.maximumMeshLevel;
  result.blockOverheadBytes = options.meshBlockOverheadBytes;
  result.memoryBudgetBytes = options.meshMemoryBudgetBytes;
  result.memoryModel = options.memoryModel;
  return result;
}

int RequestStaticCellData(int offset) {
  if (gStaticCellDataOffset >= 0) {
    SEP3D::Core::Status status(
        SEP3D::Core::StatusCode::LayoutMismatch,
        "AMPS requested srcSEP3D static cell storage more than once");
    StopWithStatus("static cell-data allocation", status);
  }
  const std::size_t bytes = Configuration().storage_layout().cellAssociatedBytes;
  if (bytes > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    SEP3D::Core::Status status(
        SEP3D::Core::StatusCode::LayoutMismatch,
        "srcSEP3D static cell-data length exceeds the AMPS integer ABI");
    StopWithStatus("static cell-data allocation", status);
  }
  gStaticCellDataOffset = offset;
  return static_cast<int>(bytes);
}

int RequestSamplingData(int offset) {
  if (gSamplingDataOffset >= 0) {
    SEP3D::Core::Status status(
        SEP3D::Core::StatusCode::LayoutMismatch,
        "AMPS requested srcSEP3D sampling storage more than once");
    StopWithStatus("sampling-data allocation", status);
  }
  const std::size_t bytes = Configuration().storage_layout().samplingBytesPerCell;
  if (bytes > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    SEP3D::Core::Status status(
        SEP3D::Core::StatusCode::LayoutMismatch,
        "srcSEP3D sampling length exceeds the AMPS integer ABI");
    StopWithStatus("sampling-data allocation", status);
  }
  gSamplingDataOffset = offset;
  return static_cast<int>(bytes);
}

struct AmpsCellReference {
  PIC::Mesh::cDataCenterNode* cell = nullptr;
  SEP3D::Core::Vec3 positionM;
  std::uint64_t stableId = 0;
  double volumeM3 = 0.0;
};

std::uint64_t StableCellId(const SEP3D::Core::Vec3& positionM) {
  std::uint64_t hash = UINT64_C(14695981039346656037);
  const double values[3] = {positionM.x, positionM.y, positionM.z};
  for (double value : values) {
    std::uint64_t bits = 0; std::memcpy(&bits, &value, sizeof(bits));
    for (unsigned shift = 0; shift < 64; shift += 8) {
      hash ^= (bits >> shift) & 0xffU;
      hash *= UINT64_C(1099511628211);
    }
  }
  return hash == 0 ? 1 : hash;
}

std::vector<AmpsCellReference> CollectOwnedPhysicalCells() {
  std::vector<AmpsCellReference> result;
  const double innerRadiusM = Configuration().options().innerRadiusM;
  for (unsigned int blockIndex = 0;
       blockIndex < PIC::DomainBlockDecomposition::nLocalBlocks;
       ++blockIndex) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
        PIC::DomainBlockDecomposition::BlockTable[blockIndex];
    if (node == nullptr || node->block == nullptr) continue;
    const double dx[3] = {
        (node->xmax[0] - node->xmin[0]) / _BLOCK_CELLS_X_,
        (node->xmax[1] - node->xmin[1]) / _BLOCK_CELLS_Y_,
        (node->xmax[2] - node->xmin[2]) / _BLOCK_CELLS_Z_};
    for (int k = 0; k < _BLOCK_CELLS_Z_; ++k) {
      for (int j = 0; j < _BLOCK_CELLS_Y_; ++j) {
        for (int i = 0; i < _BLOCK_CELLS_X_; ++i) {
          PIC::Mesh::cDataCenterNode* cell = node->block->GetCenterNode(
              PIC::Mesh::mesh->getCenterNodeLocalNumber(i, j, k));
          if (cell == nullptr) continue;
          SEP3D::Core::Vec3 position(
              node->xmin[0] + (i + 0.5) * dx[0],
              node->xmin[1] + (j + 0.5) * dx[1],
              node->xmin[2] + (k + 0.5) * dx[2]);
          // The inner sphere is a physical boundary, not a second background
          // domain.  Its allocated Cartesian cells remain zero-initialized and
          // are excluded from the complete ambient snapshot.
          if (position.Norm() >= innerRadiusM)
            result.push_back({cell, position, StableCellId(position),
                              dx[0] * dx[1] * dx[2]});
        }
      }
    }
  }
  return result;
}

std::vector<std::uint64_t> CountLocalParticlesBySpecies() {
  std::vector<std::uint64_t> counts(
      static_cast<std::size_t>(PIC::nTotalSpecies), 0);
  for (unsigned int blockIndex = 0;
       blockIndex < PIC::DomainBlockDecomposition::nLocalBlocks;
       ++blockIndex) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
        PIC::DomainBlockDecomposition::BlockTable[blockIndex];
    if (node == nullptr || node->block == nullptr) continue;
    for (int k = 0; k < _BLOCK_CELLS_Z_; ++k)
      for (int j = 0; j < _BLOCK_CELLS_Y_; ++j)
        for (int i = 0; i < _BLOCK_CELLS_X_; ++i) {
          long int ptr = node->block->FirstCellParticleTable[
              i + _BLOCK_CELLS_X_ * (j + _BLOCK_CELLS_Y_ * k)];
          while (ptr != -1) {
            PIC::ParticleBuffer::byte* data =
                PIC::ParticleBuffer::GetParticleDataPointer(ptr);
            const int species = PIC::ParticleBuffer::GetI(data);
            if (species < 0 || species >= PIC::nTotalSpecies)
              StopWithStatus("particle ledger count", SEP3D::Core::Status(
                  SEP3D::Core::StatusCode::LayoutMismatch,
                  "AMPS particle contains an out-of-range species"));
            ++counts[static_cast<std::size_t>(species)];
            ptr = PIC::ParticleBuffer::GetNext(ptr);
          }
        }
  }
  return counts;
}

void BeginParticleLedger(std::uint64_t step) {
  const std::vector<std::uint64_t> local = CountLocalParticlesBySpecies();
  gParticleLedger.Clear();
  for (int species = 0; species < PIC::nTotalSpecies; ++species) {
    const SEP3D::Core::Status status = gParticleLedger.Begin(
        step, species, local[static_cast<std::size_t>(species)]);
    if (!status.ok()) StopWithStatus("particle ledger begin", status);
  }
}

void CloseParticleLedgerCollectively(std::uint64_t step) {
  using namespace SEP3D;
  const std::vector<std::uint64_t> localEnd =
      CountLocalParticlesBySpecies();
  std::vector<Adapters::LedgerRow> globallyClosed;
  globallyClosed.reserve(static_cast<std::size_t>(PIC::nTotalSpecies));
  for (int species = 0; species < PIC::nTotalSpecies; ++species) {
    const Adapters::LedgerRow* local = gParticleLedger.Find(step, species);
    if (local == nullptr)
      StopWithStatus("particle ledger reduction", Core::Status(
          Core::StatusCode::NotFound,
          "rank-local mover ledger row is absent"));
    unsigned long long send[8] = {
        local->activeStart, local->injected, local->advanced,
        local->escaped, local->absorbed, local->failed,
        local->shockCrossings,
        localEnd[static_cast<std::size_t>(species)]};
    unsigned long long global[8] = {};
    MPI_Allreduce(send, global, 8, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                  MPI_GLOBAL_COMMUNICATOR);
    Adapters::LedgerRow row;
    row.key = {step, species};
    row.activeStart = global[0]; row.injected = global[1];
    row.advanced = global[2]; row.escaped = global[3];
    row.absorbed = global[4]; row.failed = global[5];
    row.shockCrossings = global[6]; row.activeEnd = global[7];
    row.closed = true;
    globallyClosed.push_back(row);
  }

  // Replace rank-local working rows only after every collective has finished.
  // ImportClosed performs overflow and exact conservation checks, so a lost,
  // duplicated, or post-dispatch deleted AMPS particle is a hard run failure.
  gParticleLedger.Clear();
  for (const Adapters::LedgerRow& row : globallyClosed) {
    const Core::Status status = gParticleLedger.ImportClosed(row);
    if (!status.ok()) StopWithStatus("global particle ledger closure", status);
    gClosedParticleLedger.push_back(row);
  }
}

bool SamePosition(const SEP3D::Core::Vec3& left,
                  const SEP3D::Core::Vec3& right) {
  const double scale = std::max(1.0, std::max(left.Norm(), right.Norm()));
  return (left - right).Norm() <= 1.0e-12 * scale;
}

void StoreBytes(PIC::Mesh::cDataCenterNode* cell, std::size_t offset,
                const void* source, std::size_t bytes) {
  std::memcpy(cell->GetAssociatedDataBufferPointer() + gStaticCellDataOffset +
                  offset,
              source, bytes);
}

void LoadBytes(PIC::Mesh::cDataCenterNode* cell, std::size_t offset,
               void* destination, std::size_t bytes) {
  std::memcpy(destination,
              cell->GetAssociatedDataBufferPointer() + gStaticCellDataOffset +
                  offset,
              bytes);
}

SEP3D::Core::Status ResolveLocalTransport(
    const SEP3D::Core::Vec3& positionM, int species,
    double momentumKgMPerS, double mu,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node,
    SEP3D::Adapters::LocalTransportRecord* local) {
  using namespace SEP3D;
  if (node == nullptr || node->block == nullptr || local == nullptr)
    return Core::Status(Core::StatusCode::NotFound,
                        "local transport lookup has no allocated AMR block");
  double position[3]; positionM.CopyTo(position);
  int i = 0, j = 0, k = 0;
  const int localCell = PIC::Mesh::mesh->FindCellIndex(
      position, i, j, k, node, false);
  if (localCell < 0)
    return Core::Status(Core::StatusCode::NotFound,
                        "particle position has no AMR center cell");
  PIC::Mesh::cDataCenterNode* cell = node->block->GetCenterNode(localCell);
  if (cell == nullptr)
    return Core::Status(Core::StatusCode::BackgroundInvalid,
                        "particle entered an uninitialized AMR center cell");

  const RuntimeModel::StorageLayout& layout = Configuration().storage_layout();
  Background::BackgroundSample background;
  const auto mapped = gCellSampleIndex.find(cell);
  if (gInstalledBackground && mapped != gCellSampleIndex.end() &&
      mapped->second < gInstalledBackground->samples().size()) {
    // R03 double-buffer read path: Runtime swaps the immutable snapshot
    // pointer only at a joined boundary. Movers never observe a cell array
    // while it is being filled for the next coupled generation.
    background = gInstalledBackground->samples()[mapped->second];
  } else {
  double magnetic[3], velocity[3], curvature[3];
  LoadBytes(cell, layout.magneticFieldOffset, magnetic, sizeof(magnetic));
  LoadBytes(cell, layout.bulkVelocityOffset, velocity, sizeof(velocity));
  LoadBytes(cell, layout.curvatureOffset, curvature, sizeof(curvature));
  background.B = Core::Vec3(magnetic);
  background.U = Core::Vec3(velocity);
  background.curvature = Core::Vec3(curvature);
  background.absB = background.B.Norm();
  background.bHat = background.B.Normalized();
  LoadBytes(cell, layout.numberDensityOffset, &background.numberDensityM3,
            sizeof(double));
  LoadBytes(cell, layout.velocityDivergenceOffset, &background.divU,
            sizeof(double));
  LoadBytes(cell, layout.temperatureOffset, &background.temperatureK,
            sizeof(double));
  LoadBytes(cell, layout.pressureOffset, &background.pressurePa,
            sizeof(double));
  LoadBytes(cell, layout.alfvenSpeedOffset, &background.alfvenSpeedMpS,
            sizeof(double));
  LoadBytes(cell, layout.divBhatOffset, &background.divBhat, sizeof(double));
  LoadBytes(cell, layout.focusingLengthOffset, &background.focusingLenM,
            sizeof(double));
  LoadBytes(cell, layout.fieldAlignedStrainOffset,
            &background.fieldAlignedStrain, sizeof(double));
  background.valid = std::isfinite(background.absB) && background.absB > 0.0;
  background.status = background.valid
      ? Core::Status::OK()
      : Core::Status(Core::StatusCode::BackgroundInvalid,
                     "stored AMR magnetic field is invalid");
  }
  if (!background.status.ok()) return background.status;
  local->background = background;
  const double dx = (node->xmax[0] - node->xmin[0]) / _BLOCK_CELLS_X_;
  const double dy = (node->xmax[1] - node->xmin[1]) / _BLOCK_CELLS_Y_;
  const double dz = (node->xmax[2] - node->xmin[2]) / _BLOCK_CELLS_Z_;
  local->cellSizeM = std::min(dx, std::min(dy, dz));

  if (!gInstalledTurbulence)
    return Core::Status(Core::StatusCode::SnapshotUnavailable,
                        "no prepared turbulence provider is installed");
  Turbulence::TurbulenceSample waves =
      gInstalledTurbulence->Evaluate(positionM, background);
  if (layout.waveEnergyOffset != RuntimeModel::kNoOffset) {
    double variance[2];
    LoadBytes(cell, layout.waveEnergyOffset, variance, sizeof(variance));
    waves.deltaBPlus2T2 = variance[0];
    waves.deltaBMinus2T2 = variance[1];
    waves.deltaB2T2 = variance[0] + variance[1];
  }
  if (!waves.status.usable() || !waves.valid) return waves.status;
  const Turbulence::LocalScatteringCoefficients coefficients =
      Turbulence::EvaluateLocalScattering(
          waves, background, positionM, species,
          PIC::MolecularData::GetMass(species),
          PIC::MolecularData::GetElectricCharge(species),
          momentumKgMPerS, mu);
  if (!coefficients.status.ok()) return coefficients.status;
  local->kappaParallelM2PerS = coefficients.kappaParallelM2PerS;
  local->dKappaParallelDsMPerS = 0.0;
  local->dMuMuPerS = coefficients.dMuMuPerS;
  local->dDmuMuDmuPerS = coefficients.dDmuMuDmuPerS;
  local->fractionalFieldVariationPerS = std::fabs(background.divU);
  const RuntimeModel::SnapshotDescriptor* active =
      ApplicationRuntime().active_snapshot();
  local->timeToSnapshotBoundaryS = active == nullptr ? 0.0 :
      std::max(0.0, active->validUntilS - ApplicationRuntime().CurrentTimeS());
  return Core::Status::OK();
}

void StoreBackground(PIC::Mesh::cDataCenterNode* cell,
                     const SEP3D::Background::BackgroundSample& sample) {
  const auto& layout = Configuration().storage_layout();
  const double magnetic[3] = {sample.B.x, sample.B.y, sample.B.z};
  const double velocity[3] = {sample.U.x, sample.U.y, sample.U.z};
  const double curvature[3] = {
      sample.curvature.x, sample.curvature.y, sample.curvature.z};
  StoreBytes(cell, layout.magneticFieldOffset, magnetic, sizeof(magnetic));
  StoreBytes(cell, layout.bulkVelocityOffset, velocity, sizeof(velocity));
  StoreBytes(cell, layout.numberDensityOffset, &sample.numberDensityM3,
             sizeof(double));
  StoreBytes(cell, layout.velocityDivergenceOffset, &sample.divU,
             sizeof(double));
  StoreBytes(cell, layout.temperatureOffset, &sample.temperatureK,
             sizeof(double));
  StoreBytes(cell, layout.pressureOffset, &sample.pressurePa, sizeof(double));
  StoreBytes(cell, layout.alfvenSpeedOffset, &sample.alfvenSpeedMpS,
             sizeof(double));
  StoreBytes(cell, layout.divBhatOffset, &sample.divBhat, sizeof(double));
  StoreBytes(cell, layout.focusingLengthOffset, &sample.focusingLenM,
             sizeof(double));
  StoreBytes(cell, layout.curvatureOffset, curvature, sizeof(curvature));
  StoreBytes(cell, layout.fieldAlignedStrainOffset,
             &sample.fieldAlignedStrain, sizeof(double));
  if (layout.magneticGradientOffset != SEP3D::RuntimeModel::kNoOffset)
    StoreBytes(cell, layout.magneticGradientOffset, sample.gradB.m,
               sizeof(sample.gradB.m));
  if (layout.velocityGradientOffset != SEP3D::RuntimeModel::kNoOffset)
    StoreBytes(cell, layout.velocityGradientOffset, sample.gradU.m,
               sizeof(sample.gradU.m));
}

std::shared_ptr<SEP3D::Turbulence::TurbulenceProvider>
MakePrescribedTurbulence() {
  const auto& options = Configuration().options();
  SEP3D::Turbulence::PrescribedKolmogorovConfiguration model;
  model.deltaBOverB = options.prescribedDeltaBOverB;
  model.kMinAtReferencePerM = options.turbulenceKMinPerM;
  model.kMaxAtReferencePerM = options.turbulenceKMaxPerM;
  model.spectralIndex = options.turbulenceSpectralIndex;
  model.parallelCorrelationLengthAtReferenceM =
      options.turbulenceCorrelationLengthM;
  return std::shared_ptr<SEP3D::Turbulence::TurbulenceProvider>(
      new SEP3D::Turbulence::PrescribedKolmogorovProvider(model));
}

void FillAndPublishBackground() {
  using namespace SEP3D;
  const std::vector<AmpsCellReference> cells = CollectOwnedPhysicalCells();
  std::vector<Core::Vec3> positions;
  positions.reserve(cells.size());
  for (const AmpsCellReference& cell : cells) positions.push_back(cell.positionM);

  std::shared_ptr<const Background::BackgroundSnapshot> snapshot =
      gInstalledBackground;
  if (!snapshot) {
    if (Configuration().options().background !=
        RuntimeModel::BackgroundAuthority::AnalyticParker) {
      StopWithStatus("background initialization", Core::Status(
          Core::StatusCode::SnapshotUnavailable,
          "SWMF authority requires InstallBackgroundSnapshot before amps_init"));
    }
    Background::ParkerConfiguration parker;
    const RuntimeModel::ParkerPhysicsOptions& configured =
        Configuration().options().parker;
    parker.sourceRadiusM = configured.sourceRadiusM;
    parker.sourceLongitudeRad = configured.sourceLongitudeRad;
    parker.sourceColatitudeRad = configured.sourceColatitudeRad;
    parker.referenceRadiusM = configured.referenceRadiusM;
    parker.radialFieldAtReferenceT = configured.radialFieldAtReferenceT;
    parker.numberDensityAtReferenceM3 = configured.numberDensityAtReferenceM3;
    parker.temperatureK = configured.temperatureK;
    parker.solarWindSpeedMPerS = configured.solarWindSpeedMPerS;
    parker.solarRotationRateRadPerS = configured.solarRotationRateRadPerS;
    parker.magneticPolarity = configured.magneticPolarity;
    parker.validityCadenceS = configured.validityCadenceS;
    parker.coordinateFrame = configured.coordinateFrame;
    gAnalyticBackgroundProvider.reset(
        new Background::AnalyticParkerProvider(parker));
    const double initialEpochS = ApplicationRuntime().CurrentTimeS();
    Core::Status status = gAnalyticBackgroundProvider->Prepare(initialEpochS);
    if (!status.ok()) StopWithStatus("Parker preparation", status);
    Background::BackgroundSnapshotBuilder builder;
    status = builder.Build(*gAnalyticBackgroundProvider, positions, &snapshot);
    if (!status.ok()) StopWithStatus("Parker snapshot build", status);
    if (gPendingRestart) {
      Background::SnapshotMetadata restored = snapshot->metadata();
      restored.epochS = gPendingRestart->activeSnapshot.epochS;
      restored.validFromS = gPendingRestart->activeSnapshot.validFromS;
      restored.validUntilS = gPendingRestart->activeSnapshot.validUntilS;
      restored.generation = gPendingRestart->backgroundGeneration;
      std::vector<Background::BackgroundSample> samples = snapshot->samples();
      for (Background::BackgroundSample& sample : samples)
        sample.generation = restored.generation;
      snapshot.reset(new Background::BackgroundSnapshot(
          restored, snapshot->capabilities(), snapshot->positions(), samples));
    }
    gInstalledBackground = snapshot;
  }

  if (snapshot->positions().size() != cells.size() ||
      snapshot->samples().size() != cells.size()) {
    StopWithStatus("background cell mapping", Core::Status(
        Core::StatusCode::LayoutMismatch,
        "snapshot size does not match the owner-local physical-cell list"));
  }
  if (gPendingRestart &&
      Configuration().options().background ==
          RuntimeModel::BackgroundAuthority::Swmf &&
      snapshot->metadata().generation != gPendingRestart->backgroundGeneration) {
    StopWithStatus("coupled restart snapshot", Core::Status(
        Core::StatusCode::SnapshotUnavailable,
        "installed SWMF snapshot generation differs from the checkpoint"));
  }
  for (std::size_t i = 0; i < cells.size(); ++i) {
    if (!SamePosition(snapshot->positions()[i], cells[i].positionM)) {
      StopWithStatus("background cell mapping", Core::Status(
          Core::StatusCode::LayoutMismatch,
          "snapshot positions are not in deterministic owner-cell order"));
    }
    StoreBackground(cells[i].cell, snapshot->samples()[i]);
    gCellSampleIndex[cells[i].cell] = i;
  }

  RuntimeModel::Runtime& runtime = ApplicationRuntime();
  Core::Status status;
  if (Configuration().options().background ==
      RuntimeModel::BackgroundAuthority::AnalyticParker) {
    RuntimeModel::StandaloneAdapter adapter;
    status = adapter.Initialize(&runtime);
    if (status.ok()) status = adapter.PublishSnapshot(&runtime, *snapshot);
  } else {
    RuntimeModel::SwmfAdapter adapter;
    status = adapter.Initialize(&runtime);
    if (status.ok()) status = adapter.PublishSnapshot(&runtime, *snapshot);
  }
  if (!status.ok()) StopWithStatus("background snapshot publication", status);

  std::shared_ptr<Turbulence::TurbulenceProvider> turbulence =
      gInstalledTurbulence;
  if (!turbulence) {
    if (Configuration().options().turbulence ==
        RuntimeModel::TurbulenceAuthority::Swmf) {
      StopWithStatus("turbulence initialization", Core::Status(
          Core::StatusCode::SnapshotUnavailable,
          "SWMF turbulence authority requires InstallTurbulenceProvider"));
    }
    turbulence = MakePrescribedTurbulence();
    gInstalledTurbulence = turbulence;
  }
  if (gPendingRestart) {
    std::shared_ptr<Turbulence::PrescribedKolmogorovProvider> prescribed =
        std::dynamic_pointer_cast<Turbulence::PrescribedKolmogorovProvider>(
            turbulence);
    status = prescribed
        ? prescribed->PrepareGeneration(snapshot->metadata().epochS,
              gPendingRestart->turbulenceGeneration)
        : turbulence->Prepare(snapshot->metadata().epochS);
  } else {
    status = turbulence->Prepare(snapshot->metadata().epochS);
  }
  if (!status.ok()) StopWithStatus("turbulence preparation", status);
  for (std::size_t i = 0; i < cells.size(); ++i) {
    const Turbulence::TurbulenceSample waves = turbulence->Evaluate(
        cells[i].positionM, snapshot->samples()[i]);
    if (!waves.status.usable() || !waves.valid)
      StopWithStatus("turbulence cell evaluation", waves.status);
    const std::size_t offset =
        Configuration().storage_layout().waveEnergyOffset;
    if (offset != RuntimeModel::kNoOffset) {
      const double variance[2] = {
          waves.deltaBPlus2T2, waves.deltaBMinus2T2};
      StoreBytes(cells[i].cell, offset, variance, sizeof(variance));
    }
  }
}

void RefreshBackgroundAtBoundary() {
  using namespace SEP3D;
  RuntimeModel::Runtime& runtime = ApplicationRuntime();
  if (!runtime.EventDue(RuntimeModel::ScheduledEvent::Background)) return;

  const std::vector<AmpsCellReference> cells = CollectOwnedPhysicalCells();
  std::vector<Core::Vec3> positions;
  positions.reserve(cells.size());
  for (const AmpsCellReference& cell : cells) positions.push_back(cell.positionM);
  std::shared_ptr<const Background::BackgroundSnapshot> candidate;
  std::shared_ptr<Turbulence::TurbulenceProvider> turbulence;
  Core::Status status;
  const double epochS = runtime.CurrentTimeS();

  if (Configuration().options().background ==
      RuntimeModel::BackgroundAuthority::AnalyticParker) {
    if (!gAnalyticBackgroundProvider)
      StopWithStatus("analytic snapshot update", Core::Status(
          Core::StatusCode::SnapshotUnavailable,
          "analytic provider ownership was lost after initialization"));
    status = gAnalyticBackgroundProvider->Prepare(epochS);
    if (status.ok()) {
      Background::BackgroundSnapshotBuilder builder;
      status = builder.Build(*gAnalyticBackgroundProvider, positions, &candidate);
    }
    turbulence = gInstalledTurbulence;
  } else {
    candidate = gStagedBackground;
    turbulence = gStagedTurbulence;
    if (!candidate || !turbulence)
      status = Core::Status(
          Core::StatusCode::SnapshotUnavailable,
          "SWMF cadence reached without staged background and turbulence");
  }
  if (!status.ok()) StopWithStatus("background update fill", status);
  if (!candidate || candidate->positions().size() != cells.size() ||
      candidate->samples().size() != cells.size())
    StopWithStatus("background update grid", Core::Status(
        Core::StatusCode::LayoutMismatch,
        "staged snapshot does not match deterministic owner-cell order"));

  const Background::SnapshotMetadata& metadata = candidate->metadata();
  status = runtime.RequestSnapshotUpdate(metadata.epochS, metadata.generation);
  if (status.ok()) status = runtime.BeginSnapshotFill();
  if (!status.ok()) StopWithStatus("Runtime snapshot update begin", status);

  RuntimeModel::SnapshotDescriptor descriptor;
  descriptor.authority = Configuration().options().background;
  descriptor.epochS = metadata.epochS;
  descriptor.validFromS = metadata.validFromS;
  descriptor.validUntilS = metadata.validUntilS;
  descriptor.generation = metadata.generation;
  descriptor.complete = true;
  descriptor.coordinateFrame = metadata.coordinateFrame;
  descriptor.providerIdentity = metadata.providerIdentity;
  descriptor.configurationFingerprint = Configuration().physics_fingerprint();
  status = runtime.StageSnapshot(descriptor);
  if (!status.ok()) {
    (void)runtime.FailSnapshotUpdate(status.message);
    StopWithStatus("Runtime snapshot staging", status);
  }

  status = turbulence->Prepare(epochS);
  if (!status.ok()) {
    (void)runtime.FailSnapshotUpdate(status.message);
    StopWithStatus("turbulence update preparation", status);
  }
  for (std::size_t index = 0; index < cells.size(); ++index) {
    if (!SamePosition(candidate->positions()[index], cells[index].positionM)) {
      status = Core::Status(Core::StatusCode::LayoutMismatch,
                            "staged snapshot position order changed");
      break;
    }
    const Turbulence::TurbulenceSample waves = turbulence->Evaluate(
        cells[index].positionM, candidate->samples()[index]);
    if (!waves.status.usable() || !waves.valid) { status = waves.status; break; }
    // Cell storage is a diagnostic/cache copy. The immutable shared_ptr above
    // remains the mover authority and is not swapped until collective commit.
    StoreBackground(cells[index].cell, candidate->samples()[index]);
    const std::size_t waveOffset =
        Configuration().storage_layout().waveEnergyOffset;
    if (waveOffset != RuntimeModel::kNoOffset) {
      const double variance[2] = {
          waves.deltaBPlus2T2, waves.deltaBMinus2T2};
      StoreBytes(cells[index].cell, waveOffset, variance, sizeof(variance));
    }
  }
  int localReady = status.ok() ? 1 : 0;
  int globallyReady = 0;
  MPI_Allreduce(&localReady, &globallyReady, 1, MPI_INT, MPI_MIN,
                MPI_GLOBAL_COMMUNICATOR);
  if (!status.ok() || globallyReady == 0) {
    (void)runtime.FailSnapshotUpdate(
        status.ok() ? "another rank rejected the staged snapshot" : status.message);
    StopWithStatus("collective snapshot validation",
                   status.ok() ? Core::Status(
                       Core::StatusCode::SnapshotUnavailable,
                       "another rank rejected the staged snapshot") : status);
  }
  status = runtime.PublishStagedSnapshot(true);
  if (!status.ok()) StopWithStatus("collective snapshot publication", status);
  gInstalledBackground = candidate;
  gInstalledTurbulence = turbulence;
  gStagedBackground.reset();
  gStagedTurbulence.reset();
}

struct PackedObservation {
  std::uint64_t stableId, cellId;
  std::int32_t species;
  double x, y, z, momentum, mass, mu, gyrophase, weight;
  std::uint64_t completedStep, substep, lastShockGeneration;
};

struct PackedCell {
  std::uint64_t cellId;
  double x, y, z, volume;
};

template <typename Record>
std::vector<Record> GatherRecordsToRoot(const std::vector<Record>& local) {
  const std::size_t localBytesSize = local.size() * sizeof(Record);
  if (localBytesSize > static_cast<std::size_t>(std::numeric_limits<int>::max()))
    StopWithStatus("observer MPI gather", SEP3D::Core::Status(
        SEP3D::Core::StatusCode::LayoutMismatch,
        "one rank's observer record buffer exceeds MPI int count"));
  const int localBytes = static_cast<int>(localBytesSize);
  std::vector<int> counts(PIC::nTotalThreads), displacements(PIC::nTotalThreads);
  MPI_Gather(&localBytes, 1, MPI_INT, counts.data(), 1, MPI_INT, 0,
             MPI_GLOBAL_COMMUNICATOR);
  int totalBytes = 0;
  if (PIC::ThisThread == 0) {
    for (int rank = 0; rank < PIC::nTotalThreads; ++rank) {
      displacements[rank] = totalBytes;
      if (counts[rank] < 0 || totalBytes >
          std::numeric_limits<int>::max() - counts[rank])
        StopWithStatus("observer MPI gather", SEP3D::Core::Status(
            SEP3D::Core::StatusCode::LayoutMismatch,
            "global observer record buffer exceeds MPI int displacement"));
      totalBytes += counts[rank];
    }
  }
  std::vector<Record> gathered(PIC::ThisThread == 0
      ? static_cast<std::size_t>(totalBytes) / sizeof(Record) : 0);
  MPI_Gatherv(local.empty() ? nullptr : local.data(), localBytes, MPI_BYTE,
              gathered.empty() ? nullptr : gathered.data(), counts.data(),
              displacements.data(), MPI_BYTE, 0, MPI_GLOBAL_COMMUNICATOR);
  return gathered;
}

void PublishObserversAtBoundary() {
  using namespace SEP3D;
  RuntimeModel::Runtime& runtime = ApplicationRuntime();
  if (!runtime.EventDue(RuntimeModel::ScheduledEvent::Sampling)) return;

  std::vector<PackedCell> localCells;
  std::vector<PackedObservation> localParticles;
  const std::vector<AmpsCellReference> cells = CollectOwnedPhysicalCells();
  localCells.reserve(cells.size());
  for (const AmpsCellReference& cell : cells) {
    localCells.push_back({cell.stableId, cell.positionM.x, cell.positionM.y,
                          cell.positionM.z, cell.volumeM3});
  }
  for (unsigned int blockIndex = 0;
       blockIndex < PIC::DomainBlockDecomposition::nLocalBlocks;
       ++blockIndex) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
        PIC::DomainBlockDecomposition::BlockTable[blockIndex];
    if (node == nullptr || node->block == nullptr) continue;
    for (int k = 0; k < _BLOCK_CELLS_Z_; ++k)
      for (int j = 0; j < _BLOCK_CELLS_Y_; ++j)
        for (int i = 0; i < _BLOCK_CELLS_X_; ++i) {
          PIC::Mesh::cDataCenterNode* center = node->block->GetCenterNode(
              PIC::Mesh::mesh->getCenterNodeLocalNumber(i, j, k));
          const auto found = gCellSampleIndex.find(center);
          if (found == gCellSampleIndex.end() || found->second >= cells.size())
            continue;
          long int ptr = node->block->FirstCellParticleTable[
              i + _BLOCK_CELLS_X_ * (j + _BLOCK_CELLS_Y_ * k)];
          while (ptr != -1) {
            Adapters::ParticleRecord particle;
            const Core::Status read = AMPS::Movers::ReadParticle(ptr, &particle);
            if (!read.ok()) StopWithStatus("observer particle translation", read);
            localParticles.push_back({
                particle.stableId, cells[found->second].stableId,
                static_cast<std::int32_t>(particle.species),
                particle.positionM.x, particle.positionM.y, particle.positionM.z,
                particle.momentumKgMPerS,
                PIC::MolecularData::GetMass(particle.species), particle.mu,
                particle.gyrophaseRad, particle.statisticalWeight,
                particle.completedStep, particle.substep,
                particle.lastShockGeneration});
            ptr = PIC::ParticleBuffer::GetNext(ptr);
          }
        }
  }

  const std::vector<PackedCell> globalCells = GatherRecordsToRoot(localCells);
  const std::vector<PackedObservation> globalParticles =
      GatherRecordsToRoot(localParticles);
  int publicationOK = 1;
  std::string failure;
  if (PIC::ThisThread == 0) {
    Output::SamplingRequest request;
    for (const PackedCell& cell : globalCells)
      request.cells.push_back({cell.cellId, Core::Vec3(cell.x, cell.y, cell.z),
                               cell.volume});
    for (const PackedObservation& packed : globalParticles) {
      Output::ParticleObservation particle;
      particle.stableId = packed.stableId; particle.cellId = packed.cellId;
      particle.species = packed.species;
      particle.positionM = Core::Vec3(packed.x, packed.y, packed.z);
      particle.momentumKgMPerS = packed.momentum;
      particle.restMassKg = packed.mass; particle.mu = packed.mu;
      particle.statisticalWeight = packed.weight;
      request.particles.push_back(particle);
    }
    request.ledgerRows = gClosedParticleLedger;
    Core::Status status = Output::BuildObserverDefinitions(
        Configuration(), runtime.CurrentTimeS(), &request.spacecraft);
    if (status.ok()) status = gObserverRuntime.Capture(std::move(request));
    const Output::SamplingSnapshot prepared = status.ok()
        ? gObserverRuntime.PreparePublication() : Output::SamplingSnapshot{};
    if (status.ok()) status = prepared.status;
    Output::PublicationResult published;
    if (status.ok()) {
      Output::PublicationMetadata metadata;
      metadata.sequence = runtime.counters().outputSequence;
      metadata.simulationTimeS = runtime.CurrentTimeS();
      metadata.snapshotGeneration = runtime.active_snapshot()->generation;
      metadata.configurationFingerprint = Configuration().physics_fingerprint();
      metadata.codeIdentity = "srcSEP3D-R01-R07";
      metadata.snapshotFingerprint =
          runtime.active_snapshot()->providerIdentity + ":" +
          std::to_string(runtime.active_snapshot()->generation);
      published = Output::Publish(Configuration().options().outputDirectory,
                                  Configuration().options().outputPrefix,
                                  metadata, prepared);
      status = published.status;
    }
    if (status.ok()) status = gObserverRuntime.CommitPublication(prepared);
    if (!status.ok()) { publicationOK = 0; failure = status.message; }
  }
  MPI_Bcast(&publicationOK, 1, MPI_INT, 0, MPI_GLOBAL_COMMUNICATOR);
  if (publicationOK == 0)
    StopWithStatus("observer publication", Core::Status(
        Core::StatusCode::Error,
        PIC::ThisThread == 0 ? failure :
            "root rank rejected observer publication"));
}

void WriteCheckpointAtBoundary() {
  using namespace SEP3D;
  RuntimeModel::Runtime& runtime = ApplicationRuntime();
  if (!runtime.EventDue(RuntimeModel::ScheduledEvent::Checkpoint)) return;
  Core::Status status = runtime.BeginCheckpoint();
  if (!status.ok()) StopWithStatus("checkpoint begin", status);

  std::vector<PackedObservation> localParticles;
  const std::vector<AmpsCellReference> cells = CollectOwnedPhysicalCells();
  for (unsigned int blockIndex = 0;
       blockIndex < PIC::DomainBlockDecomposition::nLocalBlocks;
       ++blockIndex) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
        PIC::DomainBlockDecomposition::BlockTable[blockIndex];
    if (node == nullptr || node->block == nullptr) continue;
    for (int k = 0; k < _BLOCK_CELLS_Z_; ++k)
      for (int j = 0; j < _BLOCK_CELLS_Y_; ++j)
        for (int i = 0; i < _BLOCK_CELLS_X_; ++i) {
          PIC::Mesh::cDataCenterNode* center = node->block->GetCenterNode(
              PIC::Mesh::mesh->getCenterNodeLocalNumber(i, j, k));
          const auto found = gCellSampleIndex.find(center);
          if (found == gCellSampleIndex.end() || found->second >= cells.size())
            continue;
          long int ptr = node->block->FirstCellParticleTable[
              i + _BLOCK_CELLS_X_ * (j + _BLOCK_CELLS_Y_ * k)];
          while (ptr != -1) {
            Adapters::ParticleRecord particle;
            status = AMPS::Movers::ReadParticle(ptr, &particle);
            if (!status.ok()) StopWithStatus("checkpoint particle translation", status);
            localParticles.push_back({
                particle.stableId, cells[found->second].stableId,
                static_cast<std::int32_t>(particle.species),
                particle.positionM.x, particle.positionM.y, particle.positionM.z,
                particle.momentumKgMPerS,
                PIC::MolecularData::GetMass(particle.species), particle.mu,
                particle.gyrophaseRad, particle.statisticalWeight,
                particle.completedStep, particle.substep,
                particle.lastShockGeneration});
            ptr = PIC::ParticleBuffer::GetNext(ptr);
          }
        }
  }
  const std::vector<PackedObservation> globalParticles =
      GatherRecordsToRoot(localParticles);
  const std::vector<Adapters::SourceLedgerRow> globalSourceLedger =
      GatherRecordsToRoot(gSourceLedger);

  int writeOK = 1;
  std::string failure;
  if (PIC::ThisThread == 0) {
    Output::RestartState checkpoint;
    checkpoint.configurationFingerprint = Configuration().physics_fingerprint();
    checkpoint.resolvedConfigurationManifest = Configuration().resolved_manifest();
    checkpoint.storageLayoutFingerprint =
        Configuration().storage_layout().fingerprint;
    checkpoint.codeIdentity = "srcSEP3D-R01-R07";
    const RuntimeModel::SnapshotDescriptor* active = runtime.active_snapshot();
    if (active != nullptr) {
      checkpoint.activeSnapshot = *active;
      checkpoint.backgroundGeneration = active->generation;
      checkpoint.snapshotFingerprint = active->providerIdentity + ":" +
          std::to_string(active->generation);
    }
    checkpoint.runtimeCounters = runtime.counters();
    ++checkpoint.runtimeCounters.checkpointSequence;
    checkpoint.eventSchedule = runtime.event_schedule();
    checkpoint.baseTimeStepS = Configuration().options().requestedTimeStepS;
    const Turbulence::TurbulenceMetadata* waves =
        gInstalledTurbulence ? gInstalledTurbulence->PreparedMetadata() : nullptr;
    checkpoint.turbulenceGeneration = waves ? waves->generation : 0;
    if (gInstalledShock) {
      checkpoint.shockState = gInstalledShock->Evaluate(runtime.CurrentTimeS());
      checkpoint.sourceGeneration = checkpoint.shockState.generation;
    } else {
      checkpoint.shockState.status = Core::Status::OK();
    }
    checkpoint.campaignSeed = Configuration().options().campaignSeed;
    checkpoint.savedRankCount = PIC::nTotalThreads;
    checkpoint.samplingState = gObserverRuntime.state();
    checkpoint.ledgerRows = gClosedParticleLedger;
    checkpoint.sourceLedgerRows = globalSourceLedger;
    std::uint64_t maximumId = 0;
    for (const PackedObservation& packed : globalParticles) {
      Adapters::ParticleRecord particle;
      particle.stableId = packed.stableId; particle.species = packed.species;
      particle.positionM = Core::Vec3(packed.x, packed.y, packed.z);
      particle.momentumKgMPerS = packed.momentum; particle.mu = packed.mu;
      particle.gyrophaseRad = packed.gyrophase;
      particle.statisticalWeight = packed.weight;
      particle.completedStep = packed.completedStep;
      particle.substep = packed.substep;
      particle.lastShockGeneration = packed.lastShockGeneration;
      checkpoint.particles.push_back(particle);
      maximumId = std::max(maximumId, particle.stableId);
    }
    checkpoint.nextStableParticleId = maximumId == UINT64_MAX
        ? 0 : std::max<std::uint64_t>(1, maximumId + 1);
    status = Output::WriteRestart(
        Configuration().options().restartOutputPath, checkpoint);
    if (!status.ok()) { writeOK = 0; failure = status.message; }
  }
  MPI_Bcast(&writeOK, 1, MPI_INT, 0, MPI_GLOBAL_COMMUNICATOR);
  if (writeOK != 0)
    status = runtime.CompleteCheckpoint();
  else
    status = runtime.AbortCheckpoint();
  if (writeOK == 0 || !status.ok())
    StopWithStatus("checkpoint commit", !status.ok() ? status : Core::Status(
        Core::StatusCode::Error,
        PIC::ThisThread == 0 ? failure : "root rank checkpoint write failed"));
}

} // namespace

// ---------------------------------------------------------------------------
// Required Exosphere hooks
//
// AMPS links these symbols for applications based on the Exosphere module.
// They are inert compatibility hooks, not srcSEP3D physics.  Named casts make
// the intentionally unused inputs explicit and keep strict-warning builds
// clean.
// ---------------------------------------------------------------------------
double Exosphere::OrbitalMotion::GetTAA(SpiceDouble et) {
  (void)et;
  return 0.0;
}

int Exosphere::ColumnIntegral::GetVariableList(char* variableList) {
  (void)variableList;
  return 0;
}

void Exosphere::ColumnIntegral::ProcessColumnIntegrationVector(
    double* result, int resultLength) {
  (void)result;
  (void)resultLength;
}

double Exosphere::GetSurfaceTemperature(double cosSubsolarAngle,
                                         double* position) {
  (void)cosSubsolarAngle;
  (void)position;
  return 0.0;
}

char Exosphere::SO_FRAME[_MAX_STRING_LENGTH_PIC_] = "HCI_like_inertial";
char Exosphere::ObjectName[_MAX_STRING_LENGTH_PIC_] = "Sun";

void Exosphere::ColumnIntegral::CoulumnDensityIntegrant(
    double* result, int resultLength, double* position,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  (void)result;
  (void)resultLength;
  (void)position;
  (void)node;
}

double Exosphere::SurfaceInteraction::StickingProbability(
    int spec, double& reemissionParticleFraction, double temperature) {
  (void)spec;
  (void)temperature;
  reemissionParticleFraction = 0.0;
  return 0.0;
}

SEP3D::RuntimeModel::Runtime& SEP3D::ApplicationRuntime() {
  // AMPS exposes process-level application callbacks, so one process-owned
  // Runtime is the matching ownership scope.  The object's state and counters
  // are encapsulated and change only through typed, transactional methods.
  static RuntimeModel::Runtime runtime;
  return runtime;
}

SEP3D::Core::Status SEP3D::ConfigureApplication(
    const std::shared_ptr<const RuntimeModel::RunConfiguration3D>& configuration) {
  return ApplicationRuntime().Configure(configuration);
}

SEP3D::Core::Status SEP3D::InstallBackgroundSnapshot(
    const std::shared_ptr<const Background::BackgroundSnapshot>& snapshot) {
  if (!snapshot) {
    return Core::Status(Core::StatusCode::InvalidInput,
                        "installed background snapshot is null");
  }
  if (!ApplicationRuntime().configuration()) {
    return Core::Status(Core::StatusCode::InvalidTransition,
                        "configure the application before installing a snapshot");
  }
  const RuntimeModel::LifecycleState state = ApplicationRuntime().state();
  if (state != RuntimeModel::LifecycleState::Configured &&
      state != RuntimeModel::LifecycleState::MeshReady &&
      state != RuntimeModel::LifecycleState::SnapshotReady) {
    return Core::Status(
        Core::StatusCode::InvalidTransition,
        "background installation is legal before acquisition or at a joined snapshot boundary");
  }
  const bool analytic =
      snapshot->metadata().provider == Background::ProviderKind::AnalyticParker;
  const bool expectedAnalytic =
      ApplicationRuntime().configuration()->options().background ==
      RuntimeModel::BackgroundAuthority::AnalyticParker;
  if (analytic != expectedAnalytic) {
    return Core::Status(Core::StatusCode::ConfigurationConflict,
                        "installed snapshot authority differs from configuration");
  }
  const Background::SnapshotMetadata& metadata = snapshot->metadata();
  if (snapshot->positions().empty() ||
      snapshot->positions().size() != snapshot->samples().size() ||
      metadata.coordinateFrame != "HCI-like-inertial" ||
      !std::isfinite(metadata.epochS) ||
      !std::isfinite(metadata.validFromS) ||
      !std::isfinite(metadata.validUntilS) ||
      metadata.validUntilS < metadata.validFromS ||
      metadata.epochS < metadata.validFromS ||
      metadata.epochS > metadata.validUntilS || metadata.generation == 0 ||
      metadata.providerIdentity.empty() ||
      metadata.configurationFingerprint.empty()) {
    return Core::Status(Core::StatusCode::SnapshotUnavailable,
                        "installed snapshot grid or coordinate frame is invalid");
  }
  for (const Background::BackgroundSample& sample : snapshot->samples()) {
    const Core::Status complete = Background::ValidateCompleteSample(
        sample, snapshot->capabilities());
    if (!complete.ok()) return complete;
  }
  if (state == RuntimeModel::LifecycleState::SnapshotReady)
    gStagedBackground = snapshot;
  else
    gInstalledBackground = snapshot;
  return Core::Status::OK();
}

SEP3D::Core::Status SEP3D::InstallTurbulenceProvider(
    const std::shared_ptr<Turbulence::TurbulenceProvider>& provider) {
  if (!provider) {
    return Core::Status(Core::StatusCode::InvalidInput,
                        "installed turbulence provider is null");
  }
  if (!ApplicationRuntime().configuration()) {
    return Core::Status(
        Core::StatusCode::InvalidTransition,
        "configure the application before installing turbulence");
  }
  const RuntimeModel::LifecycleState state = ApplicationRuntime().state();
  if (state != RuntimeModel::LifecycleState::Configured &&
      state != RuntimeModel::LifecycleState::MeshReady &&
      state != RuntimeModel::LifecycleState::SnapshotReady) {
    return Core::Status(
        Core::StatusCode::InvalidTransition,
        "turbulence installation is legal before acquisition or at a joined snapshot boundary");
  }
  const Core::Status status = provider->Validate();
  if (!status.ok()) return status;
  const bool imported =
      provider->Source() == Turbulence::TurbulenceSource::SwmfAwsom;
  const bool expectedImported =
      ApplicationRuntime().configuration()->options().turbulence ==
          RuntimeModel::TurbulenceAuthority::Swmf;
  if (imported != expectedImported) {
    return Core::Status(
        Core::StatusCode::ConfigurationConflict,
        "installed turbulence authority differs from configuration");
  }
  if (state == RuntimeModel::LifecycleState::SnapshotReady)
    gStagedTurbulence = provider;
  else
    gInstalledTurbulence = provider;
  return Core::Status::OK();
}

SEP3D::Core::Status SEP3D::InstallShockProvider(
    const std::shared_ptr<Adapters::ShockProvider>& provider) {
  if (!provider)
    return Core::Status(Core::StatusCode::InvalidInput,
                        "installed shock provider is null");
  if (!ApplicationRuntime().configuration() ||
      ApplicationRuntime().state() != RuntimeModel::LifecycleState::Configured)
    return Core::Status(Core::StatusCode::InvalidTransition,
                        "install the shock provider after configuration and before mesh binding");
  gInstalledShock = provider;
  return Core::Status::OK();
}

SEP3D::Core::Status SEP3D::InstallRestartState(
    const Output::RestartState& state) {
  if (!ApplicationRuntime().configuration() ||
      ApplicationRuntime().state() != RuntimeModel::LifecycleState::Configured)
    return Core::Status(Core::StatusCode::InvalidTransition,
                        "restart state must be installed while Runtime is Configured");
  if (state.configurationFingerprint !=
          ApplicationRuntime().configuration()->physics_fingerprint() ||
      state.runtimeCounters.currentTick !=
          ApplicationRuntime().counters().currentTick)
    return Core::Status(Core::StatusCode::ConfigurationConflict,
                        "restart state differs from restored Runtime identity or clock");
  gPendingRestart.reset(new Output::RestartState(state));
  if (state.shockState.active && !gInstalledShock) {
    std::shared_ptr<Adapters::PublishedShockProvider> restored(
        new Adapters::PublishedShockProvider("restart-shock"));
    const Core::Status published = restored->Publish(state.shockState);
    if (!published.ok()) { gPendingRestart.reset(); return published; }
    gInstalledShock = restored;
  }
  return Core::Status::OK();
}

void SEP3D::Init_BeforeParser() {
  // Register byte requests after PIC has created its request registries but
  // before initCellSamplingDataBuffer freezes the center-node layout.  The
  // callbacks return the exact sizes already fingerprinted by RunConfiguration.
  // Coupled entry points still do not inspect argc/argv or AMPS_PARAM.in.
  (void)Configuration();
  const Core::Status particleStorage = AMPS::Movers::RequestParticleStorage();
  if (!particleStorage.ok())
    StopWithStatus("particle-buffer storage request", particleStorage);
  if (!gStorageCallbacksRegistered) {
    PIC::IndividualModelSampling::RequestStaticCellData.push_back(
        RequestStaticCellData);
    if (Configuration().storage_layout().samplingBytesPerCell != 0) {
      PIC::IndividualModelSampling::RequestSamplingData.push_back(
          RequestSamplingData);
    }
    gStorageCallbacksRegistered = true;
  }
}

double localResolution(double* position) {
  if (position == nullptr) {
    StopWithStatus("localResolution", SEP3D::Core::Status(
        SEP3D::Core::StatusCode::InvalidInput,
        "AMPS supplied a null mesh position"));
  }
  const SEP3D::Mesh::ResolutionConfiguration resolution =
      ResolutionConfiguration();
  const SEP3D::Core::Status valid = SEP3D::Mesh::Validate(resolution);
  if (!valid.ok()) StopWithStatus("localResolution", valid);
  return SEP3D::Mesh::RequestedCellSizeM(
      SEP3D::Core::Vec3(position), resolution);
}

double InitLoadMeasure(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  return node != nullptr && node->IsUsedInCalculationFlag ? 1.0 : 0.0;
}

bool TrajectoryTrackingCondition(double* position, double* velocity, int spec,
                                 void* particleData) {
  (void)position;
  (void)velocity;
  (void)spec;
  (void)particleData;
  return false;
}

void amps_init_mesh() {
  if (SEP3D::ApplicationRuntime().state() ==
      SEP3D::RuntimeModel::LifecycleState::Created) {
    std::cerr
        << "[srcSEP3D] amps_init_mesh requires the host to install an "
        << "immutable RunConfiguration3D before AMPS allocates mesh storage.\n";
    std::abort();
  }
  const SEP3D::Mesh::ResolutionConfiguration resolution =
      ResolutionConfiguration();
  const SEP3D::Core::Status resolutionStatus =
      SEP3D::Mesh::Validate(resolution);
  if (!resolutionStatus.ok())
    StopWithStatus("mesh configuration", resolutionStatus);

  // PIC initializes MPI and its allocation registries.  srcSEP3D then adds
  // the frozen static/sampling requests before AMPS calculates final offsets.
  PIC::Init_BeforeParser();
  SEP3D::Init_BeforeParser();
  PIC::Mesh::initCellSamplingDataBuffer();
  if (gStaticCellDataOffset < 0) {
    StopWithStatus("mesh storage freeze", SEP3D::Core::Status(
        SEP3D::Core::StatusCode::LayoutMismatch,
        "AMPS did not invoke the srcSEP3D static-data request"));
  }

  const auto& options = Configuration().options();
  const SEP3D::Mesh::DomainBounds domain =
      SEP3D::Mesh::MakeDomain(options);
  double minimum[3] = {
      domain.minimumM.x, domain.minimumM.y, domain.minimumM.z};
  double maximum[3] = {
      domain.maximumM.x, domain.maximumM.y, domain.maximumM.z};

  // Build and partition the AMR tree before block allocation.  This is the
  // same ordering used by mature AMPS applications and guarantees that each
  // MPI rank fills only the blocks it owns after decomposition.
  PIC::Mesh::mesh->AllowBlockAllocation = false;
  PIC::Mesh::mesh->init(minimum, maximum, localResolution);
  PIC::Mesh::mesh->buildMesh();
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  PIC::Mesh::mesh->SetParallelLoadMeasure(InitLoadMeasure);
  PIC::Mesh::mesh->CreateNewParallelDistributionLists();
  PIC::Mesh::mesh->AllowBlockAllocation = true;
  PIC::Mesh::mesh->AllocateTreeBlocks();
  PIC::Mesh::mesh->InitCellMeasure();
  PIC::Mesh::mesh->memoryAllocationReport();
  PIC::Mesh::mesh->GetMeshTreeStatistics();

  SEP3D::RuntimeModel::MeshBinding binding;
  binding.layout = Configuration().storage_layout();
  const SEP3D::Core::Status bound =
      SEP3D::ApplicationRuntime().BindMesh(binding);
  if (!bound.ok()) StopWithStatus("Runtime mesh binding", bound);
}

void amps_init() {
  PIC::Init_AfterParser();
  // R04: one configured base step is authoritative for Runtime and every AMPS
  // species/block.  Assigning this after PIC initialization prevents legacy
  // default setup from silently replacing the user-resolved cadence.
  const double configuredDt = Configuration().options().requestedTimeStepS;
  for (int species = 0; species < PIC::nTotalSpecies; ++species)
    PIC::ParticleWeightTimeStep::GlobalTimeStep[species] = configuredDt;
  for (int species = 0; species < PIC::nTotalSpecies; ++species)
    PIC::ParticleWeightTimeStep::GlobalParticleWeight[species] =
        Configuration().options().species.macroparticleWeight;
  PIC::ParticleWeightTimeStep::GlobalTimeStepInitialized = true;
  for (unsigned int blockIndex = 0;
       blockIndex < PIC::DomainBlockDecomposition::nLocalBlocks;
       ++blockIndex) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
        PIC::DomainBlockDecomposition::BlockTable[blockIndex];
    if (node == nullptr || node->block == nullptr) continue;
    for (int species = 0; species < PIC::nTotalSpecies; ++species)
      node->block->SetLocalTimeStep(configuredDt, species);
  }
  FillAndPublishBackground();
  SEP3D::AMPS::Movers::Context mover;
  mover.resolveLocal = ResolveLocalTransport;
  mover.ledger = &gParticleLedger;
  mover.maximumSubsteps = Configuration().options().maximumTransportSubsteps;
  if (gInstalledShock) {
    const SEP3D::Adapters::ShockState shock =
        gInstalledShock->Evaluate(SEP3D::ApplicationRuntime().CurrentTimeS());
    if (!shock.status.ok()) StopWithStatus("initial shock state", shock.status);
    mover.shock.active = shock.active;
    mover.shock.centerM = shock.centerM;
    mover.shock.radiusAtStepStartM = shock.radiusM;
    mover.shock.radialSpeedMPerS = shock.radialSpeedMPerS;
    mover.shock.generation = shock.generation;
  } else if (Configuration().options().source.enabled) {
    StopWithStatus("source initialization", SEP3D::Core::Status(
        SEP3D::Core::StatusCode::SnapshotUnavailable,
        "enabled shock source requires InstallShockProvider before mesh setup"));
  }
  const SEP3D::Core::Status installed =
      SEP3D::AMPS::Movers::InstallContext(mover);
  if (!installed.ok()) StopWithStatus("AMPS mover context installation", installed);
  if (gPendingRestart) {
    PIC::SimulationTime::SetInitialValue(
        SEP3D::ApplicationRuntime().CurrentTimeS());
    SEP3D::Adapters::InjectionPlan restored;
    restored.status = SEP3D::Core::Status::OK();
    for (const SEP3D::Adapters::ParticleRecord& particle :
         gPendingRestart->particles) {
      SEP3D::Adapters::InjectedParticle wrapped;
      wrapped.status = SEP3D::Core::Status::OK();
      wrapped.particle = particle;
      restored.particles.push_back(wrapped);
    }
    const SEP3D::AMPS::Movers::InjectionOutcome loaded =
        SEP3D::AMPS::Movers::InjectParticles(restored);
    if (!loaded.status.ok()) StopWithStatus("restart particle restore", loaded.status);
    gObserverRuntime.RestoreState(gPendingRestart->samplingState);
    gClosedParticleLedger = gPendingRestart->ledgerRows;
    // Source rows are rank-local until the next checkpoint gather.  Restoring
    // the canonical historical table on root only prevents an N-rank restart
    // from multiplying all pre-restart source evidence by N.
    if (PIC::ThisThread == 0)
      gSourceLedger = gPendingRestart->sourceLedgerRows;
    else
      gSourceLedger.clear();
  }
  if (PIC::ThisThread == 0) {
    const SEP3D::RuntimeModel::EventSchedule& events =
        SEP3D::ApplicationRuntime().event_schedule();
    std::cout << "[srcSEP3D] runtime clock: tick="
              << SEP3D::ApplicationRuntime().counters().currentTick
              << " time_s=" << SEP3D::ApplicationRuntime().CurrentTimeS()
              << " dt_s=" << configuredDt
              << " next_background_tick=" << events.nextBackgroundTick
              << " next_injection_tick=" << events.nextInjectionTick
              << " next_sampling_tick=" << events.nextSamplingTick
              << " next_checkpoint_tick=" << events.nextCheckpointTick
              << '\n';
  }
}

int amps_time_step() {
  // Pin one immutable background generation for the whole AMPS particle
  // phase. All workers therefore observe identical coefficients even if a
  // coupled host stages the next SWMF snapshot concurrently.
  SEP3D::RuntimeModel::Runtime& runtime = SEP3D::ApplicationRuntime();
  SEP3D::RuntimeModel::ClockObservation clocks;
  clocks.picTimeS = PIC::SimulationTime::Get();
  clocks.picTimeStepS = PIC::ParticleWeightTimeStep::GlobalTimeStep[0];
  clocks.snapshotTimeS = runtime.CurrentTimeS();
  clocks.shockTimeS = runtime.CurrentTimeS();
  SEP3D::Core::Status status = runtime.VerifyClockAgreement(clocks);
  if (!status.ok()) StopWithStatus("R04 clock agreement", status);
  status = runtime.BeginStep(PIC::SimulationTime::Get());
  if (!status.ok()) StopWithStatus("Runtime BeginStep", status);
  const std::uint64_t particleStep = runtime.counters().completedSteps;
  BeginParticleLedger(particleStep);
  const int returnCode = PIC::TimeStep();
  CloseParticleLedgerCollectively(particleStep);
  status = runtime.CompleteStep();
  if (!status.ok()) StopWithStatus("Runtime CompleteStep", status);
  RefreshBackgroundAtBoundary();

  // R05 source creation is a joined-boundary operation.  The provider state,
  // stochastic rounding, AMPS allocation, and conservation ledger therefore
  // all refer to the same authoritative integer tick.
  if (gInstalledShock) {
    const SEP3D::Adapters::ShockState shock =
        gInstalledShock->Evaluate(runtime.CurrentTimeS());
    if (!shock.status.ok()) StopWithStatus("shock update", shock.status);
    SEP3D::Adapters::ExpandingSphericalShock moverShock;
    moverShock.active = shock.active;
    moverShock.centerM = shock.centerM;
    moverShock.radiusAtStepStartM = shock.radiusM;
    moverShock.radialSpeedMPerS = shock.radialSpeedMPerS;
    moverShock.generation = shock.generation;
    status = SEP3D::AMPS::Movers::UpdateShock(moverShock);
    if (!status.ok()) StopWithStatus("mover shock update", status);

    if (Configuration().options().source.enabled &&
        runtime.EventDue(SEP3D::RuntimeModel::ScheduledEvent::Injection)) {
      const auto& options = Configuration().options();
      for (const SEP3D::Adapters::ShockSourceRecord& patch : shock.patches) {
        SEP3D::Adapters::SourceRequest source;
        source.patch = patch;
        source.step = runtime.counters().currentTick;
        source.species = 0;
        source.speciesMassKg = PIC::MolecularData::GetMass(0);
        source.intervalS = options.requestedTimeStepS *
            options.injectionCadenceSteps;
        source.physicalParticleRatePerS =
            options.source.physicalParticleRatePerS *
            patch.relativePatchWeight;
        source.macroparticleWeight = options.species.macroparticleWeight;
        source.maximumMacroparticles = options.source.samplesPerStep;
        SEP3D::Adapters::InjectionPlan plan =
            SEP3D::Adapters::BuildInjectionPlan(source);
        if (!plan.status.ok()) StopWithStatus("shock source plan", plan.status);
        const SEP3D::AMPS::Movers::InjectionOutcome injected =
            SEP3D::AMPS::Movers::InjectParticles(plan);
        if (!injected.status.ok())
          StopWithStatus("AMPS source allocation", injected.status);
        // All ranks construct the same deterministic plan, but AMPS allocates
        // it only on the rank owning the shock-patch position.  Record the
        // physical source row there so the checkpoint gather contains exactly
        // one copy rather than one replicated row per MPI rank.
        if (injected.allocated != 0 ||
            (plan.particles.empty() && PIC::ThisThread == 0))
          gSourceLedger.push_back(plan.ledger);
      }
    }
  }
  // R06 ordering is intentional: observers see particles after both motion
  // and this boundary's shock injection.  Publication is globally gathered,
  // deterministic by stable ID, and its accumulator clears only on commit.
  PublishObserversAtBoundary();
  WriteCheckpointAtBoundary();
  return returnCode;
}

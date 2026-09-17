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
#include "mesh/mesh_model.h"
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
std::shared_ptr<SEP3D::Turbulence::TurbulenceProvider>
    gInstalledTurbulence;

int gStaticCellDataOffset = -1;
int gSamplingDataOffset = -1;
bool gStorageCallbacksRegistered = false;

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
  result.enableTubeRefinement = options.enableTubeRefinement;
  result.tubeLongitudeRad = options.tubeLongitudeRad;
  result.tubeColatitudeRad = options.tubeColatitudeRad;
  result.tubePolarity = options.tubePolarity;
  result.tubeCoreRadiusM = options.tubeCoreRadiusM;
  result.tubeShoulderRadiusM = options.tubeShoulderRadiusM;
  result.tubeCellSizeM = options.tubeCellSizeM;
  result.cellsPerBlockEdge = options.meshCellsPerBlockEdge;
  result.maximumLevel = options.maximumMeshLevel;
  result.blockOverheadBytes = options.meshBlockOverheadBytes;
  result.memoryBudgetBytes = options.meshMemoryBudgetBytes;
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
};

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
            result.push_back({cell, position});
        }
      }
    }
  }
  return result;
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
    parker.sourceRadiusM = Configuration().options().innerRadiusM;
    parker.validityCadenceS =
        Configuration().options().requestedTimeStepS *
        Configuration().options().backgroundCadenceSteps;
    Background::AnalyticParkerProvider provider(parker);
    Core::Status status = provider.Prepare(0.0);
    if (!status.ok()) StopWithStatus("Parker preparation", status);
    Background::BackgroundSnapshotBuilder builder;
    status = builder.Build(provider, positions, &snapshot);
    if (!status.ok()) StopWithStatus("Parker snapshot build", status);
  }

  if (snapshot->positions().size() != cells.size() ||
      snapshot->samples().size() != cells.size()) {
    StopWithStatus("background cell mapping", Core::Status(
        Core::StatusCode::LayoutMismatch,
        "snapshot size does not match the owner-local physical-cell list"));
  }
  for (std::size_t i = 0; i < cells.size(); ++i) {
    if (!SamePosition(snapshot->positions()[i], cells[i].positionM)) {
      StopWithStatus("background cell mapping", Core::Status(
          Core::StatusCode::LayoutMismatch,
          "snapshot positions are not in deterministic owner-cell order"));
    }
    StoreBackground(cells[i].cell, snapshot->samples()[i]);
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
  }
  status = turbulence->Prepare(snapshot->metadata().epochS);
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
      state != RuntimeModel::LifecycleState::MeshReady) {
    return Core::Status(
        Core::StatusCode::InvalidTransition,
        "background installation is legal only before acquisition begins");
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
      state != RuntimeModel::LifecycleState::MeshReady) {
    return Core::Status(
        Core::StatusCode::InvalidTransition,
        "turbulence installation is legal only before acquisition begins");
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
  gInstalledTurbulence = provider;
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
  const SEP3D::Mesh::DomainBounds domain = SEP3D::Mesh::MakeDomain(
      options.domain, options.innerRadiusM, options.outerRadiusM);
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
  FillAndPublishBackground();
}

int amps_time_step() {
  // Pin one immutable background generation for the whole AMPS particle
  // phase. All workers therefore observe identical coefficients even if a
  // coupled host stages the next SWMF snapshot concurrently.
  SEP3D::RuntimeModel::Runtime& runtime = SEP3D::ApplicationRuntime();
  SEP3D::Core::Status status = runtime.BeginStep(PIC::SimulationTime::Get());
  if (!status.ok()) StopWithStatus("Runtime BeginStep", status);
  const int returnCode = PIC::TimeStep();
  status = runtime.CompleteStep();
  if (!status.ok()) StopWithStatus("Runtime CompleteStep", status);
  return returnCode;
}

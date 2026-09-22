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
#include "output/sampling.h"
#include "runtime/runtime_adapters.h"
#include "turbulence/turbulence_models.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

namespace fs = std::filesystem;
constexpr double kMagneticPermeabilityVacuum =
    4.0e-7 * SEP3D::Core::Const::kPi;

[[noreturn]] void StopWithStatus(const char* operation,
                                 const SEP3D::Core::Status& status) {
  std::cerr << "[srcSEP3D] " << operation << " failed: "
            << status.message << '\n';
  std::abort();
}

// Create all initialization-product parents once on rank zero, then publish the
// result to every rank before AMPS opens its distributed Tecplot mesh.  Doing
// this after PIC has initialized MPI avoids a many-rank mkdir race on shared
// filesystems.  It also covers paths declared directly in the input deck, not
// only paths rewritten by --initialization-output-dir.
void EnsureInitializationOutputParents(
    const SEP3D::RuntimeModel::RunConfiguration3DOptions& options) {
  int directoryStatus = 1;
  std::string rootMessage;
  if (PIC::ThisThread == 0) {
    const std::string paths[] = {
        options.initializationMeshTecplotFile,
        options.initializationParkerLineTecplotFile,
        options.initializationDataTecplotFile};
    for (const std::string& pathText : paths) {
      const fs::path parent = fs::path(pathText).parent_path();
      if (parent.empty()) continue;
      std::error_code error;
      fs::create_directories(parent, error);
      const bool isDirectory = fs::is_directory(parent, error);
      if (error || !isDirectory) {
        directoryStatus = 0;
        rootMessage = "cannot create initialization output directory '" +
            parent.string() + "'";
        if (error) rootMessage += ": " + error.message();
        break;
      }
    }
  }
  MPI_Bcast(&directoryStatus, 1, MPI_INT, 0, MPI_GLOBAL_COMMUNICATOR);
  if (directoryStatus == 0) {
    StopWithStatus("initialization output directory",
        SEP3D::Core::Status(SEP3D::Core::StatusCode::Error,
            PIC::ThisThread == 0 ? rootMessage
                                 : "rank zero could not create the directory"));
  }
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
}

// AMPS' native data writer accepts one DataSetNumber, which is the compiled
// species index used for block-local time step and weight. Preserve the exact
// configured name for a one-species build; for mixed SpeciesList builds,
// create deterministic sibling files so no species silently overwrites or
// masquerades as another.
std::string InitializationDataPath(const std::string& base, int species,
                                   int speciesCount) {
  if (speciesCount == 1) return base;
  const fs::path path(base);
  const std::string name = path.stem().string() + ".species-" +
      std::to_string(species) + path.extension().string();
  return (path.parent_path() / name).string();
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

// The AMPS configurator rewrites both constants from the same SpeciesList.
// Keeping this compile-time assertion at the application boundary catches a
// partially regenerated tree before either array can be indexed.
static_assert(_TOTAL_SPECIES_NUMBER_ == PIC::nTotalSpecies,
              "AMPS species count constants disagree; regenerate the application");

// Read-only mirror of the generated AMPS molecular table.  It is populated
// once after PIC::Init_BeforeParser and is then used as the authoritative
// iteration order for numerical initialization and injection.  No entry is
// synthesized from a species macro and no AMPS mass/charge is overwritten.
std::vector<SEP3D::RuntimeModel::CompiledSpeciesRecord> gCompiledSpecies;

int gStaticCellDataOffset = -1;
int gSamplingDataOffset = -1;
bool gStorageCallbacksRegistered = false;
// The initialization Tecplot product is not permitted to run until the same
// validated background generation has been copied to both srcSEP3D's frozen
// application storage and the native AMPS coupler storage used by AMPS' own
// output callback.  This flag records completion of that collective boundary;
// it is deliberately reset during each background refresh.
bool gNativeAmpsBackgroundReady = false;
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

void BindCompiledSpeciesTable() {
  gCompiledSpecies.clear();
  gCompiledSpecies.reserve(static_cast<std::size_t>(PIC::nTotalSpecies));

  // PIC::nTotalSpecies and ChemTable are generated by the AMPS preprocessing
  // of SpeciesList.  Enumerating the table is the only general solution: an
  // ELECTRON-only executable has no active _H_PLUS_SPEC_, and mixed tables do
  // not promise that any particular chemical macro occupies slot zero.
  for (int speciesIndex = 0; speciesIndex < PIC::nTotalSpecies;
       ++speciesIndex) {
    SEP3D::RuntimeModel::CompiledSpeciesRecord record;
    record.ampsIndex = speciesIndex;
    const char* symbol = PIC::MolecularData::GetChemSymbol(speciesIndex);
    if (symbol != nullptr) record.symbol = symbol;
    record.massKg = PIC::MolecularData::GetMass(speciesIndex);
    record.chargeC = PIC::MolecularData::GetElectricCharge(speciesIndex);
    gCompiledSpecies.push_back(std::move(record));
  }

  const SEP3D::Core::Status binding =
      SEP3D::RuntimeModel::ValidateCompiledSpeciesBinding(
          Configuration().options(), static_cast<int>(PIC::nTotalSpecies),
          gCompiledSpecies);
  if (!binding.ok()) StopWithStatus("AMPS compiled species binding", binding);

  if (PIC::ThisThread == 0) {
    std::ostringstream message;
    message << std::scientific << std::setprecision(17)
            << "[srcSEP3D] bound " << gCompiledSpecies.size()
            << " immutable AMPS species from SpeciesList\n";
    for (const auto& species : gCompiledSpecies) {
      message << "  index=" << species.ampsIndex
              << " symbol='" << species.symbol << "'"
              << " mass_kg=" << species.massKg
              << " charge_C=" << species.chargeC << '\n';
    }
    std::cout << message.str();
  }
}

SEP3D::Mesh::ResolutionConfiguration ResolutionConfiguration() {
  const auto& options = Configuration().options();
  SEP3D::Mesh::ResolutionConfiguration result;
  result.originM = options.coordinateOriginM;
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
  result.parkerInitialPointM = options.parkerSpiralInitialPointM;
  result.parkerLengthM = options.parkerSpiralLengthM;
  result.parkerPointCount = options.parkerSpiralPointCount;
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
  const double outerRadiusM = Configuration().options().outerRadiusM;
  const SEP3D::Core::Vec3 originM =
      Configuration().options().coordinateOriginM;
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
          // Both spherical boundaries delimit the physical background domain.
          // The enclosing Cartesian cube also allocates cells inside the Sun
          // and outside the requested heliocentric radius; those padding cells
          // stay zero-initialized and are explicitly marked background_valid=0
          // in Tecplot output.  Measure radius from the configured heliocentric
          // origin rather than silently assuming an origin at (0,0,0).
          const double radiusM = (position - originM).Norm();
          if (radiusM >= innerRadiusM && radiusM <= outerRadiusM)
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

// AMPS writes a FEBRICK zone at mesh vertices.  Its writer obtains every
// vertex value by creating a temporary cDataCenterNode and calling
// cDataCenterNode::Interpolate() with the surrounding physical centre nodes.
// That AMPS method knows how to interpolate built-in particle sampling and
// DATAFILE fields, but deliberately knows nothing about static bytes requested
// by an application.  Without this callback, the real cells contain the
// initialized srcSEP3D turbulence state while the temporary node printed by
// Tecplot retains allocation-time zeros.
//
// Registering this callback before AMPS freezes the centre-node layout makes
// the complete srcSEP3D static slice participate in the same interpolation as
// the native IMF/plasma fields.  The frozen slice consists only of doubles,
// so background primitives, optional gradients, and deltaB_+/-^2 are
// interpolated together and remain registered at every output vertex.
void InterpolateInitializationCellData(
    PIC::Mesh::cDataCenterNode** interpolationList,
    double* interpolationCoefficients, int interpolationCount,
    PIC::Mesh::cDataCenterNode* destinationNode) {
  const std::size_t bytes =
      Configuration().storage_layout().cellAssociatedBytes;
  if (gStaticCellDataOffset < 0 || bytes == 0 ||
      bytes % sizeof(double) != 0 || interpolationList == nullptr ||
      interpolationCoefficients == nullptr || interpolationCount <= 0 ||
      destinationNode == nullptr) {
    StopWithStatus("srcSEP3D center-node interpolation",
        SEP3D::Core::Status(SEP3D::Core::StatusCode::LayoutMismatch,
            "AMPS requested interpolation before the frozen srcSEP3D "
            "static center-node layout was available"));
  }

  const std::size_t valueCount = bytes / sizeof(double);
  // The byte offset is owned by AMPS and is not assumed to satisfy C++ double
  // alignment.  Copy each slice into aligned vector storage before doing
  // floating-point arithmetic; StoreBytes/LoadBytes use memcpy for the same
  // reason everywhere else at this boundary.
  std::vector<double> alignedSources(
      static_cast<std::size_t>(interpolationCount) * valueCount, 0.0);
  std::vector<const double*> sourceValues(
      static_cast<std::size_t>(interpolationCount), nullptr);
  for (int index = 0; index < interpolationCount; ++index) {
    if (interpolationList[index] == nullptr) {
      StopWithStatus("srcSEP3D center-node interpolation",
          SEP3D::Core::Status(SEP3D::Core::StatusCode::InvalidInput,
              "AMPS supplied a null source center node"));
    }
    double* aligned = alignedSources.data() +
        static_cast<std::size_t>(index) * valueCount;
    std::memcpy(aligned,
                interpolationList[index]->GetAssociatedDataBufferPointer() +
                    gStaticCellDataOffset,
                bytes);
    sourceValues[static_cast<std::size_t>(index)] = aligned;
  }
  std::vector<double> interpolated(valueCount, 0.0);
  const SEP3D::Core::Status status =
      SEP3D::Output::InterpolateStaticCenterState(
          sourceValues.data(), interpolationCoefficients,
          sourceValues.size(), valueCount, interpolated.data());
  if (!status.ok())
    StopWithStatus("srcSEP3D center-node interpolation", status);
  std::memcpy(destinationNode->GetAssociatedDataBufferPointer() +
                  gStaticCellDataOffset,
              interpolated.data(), bytes);
}

// Append names for every initialized srcSEP3D macroscopic field stored in the
// AMPS center-node buffer. The block object itself appends the selected
// species' local time step and particle weight, so the resulting file is a
// native data-bearing AMPS Tecplot product rather than a geometry-only mesh.
void PrintInitializationVariableList(FILE* output, int dataSetNumber) {
  (void)dataSetNumber;
  std::fprintf(output,
      ", \"B_x_T\", \"B_y_T\", \"B_z_T\""
      ", \"U_x_m_per_s\", \"U_y_m_per_s\", \"U_z_m_per_s\""
      ", \"number_density_m-3\", \"div_U_s-1\", \"temperature_K\""
      ", \"pressure_Pa\", \"alfven_speed_m_per_s\", \"div_bhat_m-1\""
      ", \"focusing_length_m\""
      ", \"curvature_x_m-1\", \"curvature_y_m-1\", \"curvature_z_m-1\""
      ", \"field_aligned_strain_s-1\"");
  const auto& layout = Configuration().storage_layout();
  if (layout.magneticGradientOffset != SEP3D::RuntimeModel::kNoOffset)
    for (int row = 0; row < 3; ++row)
      for (int column = 0; column < 3; ++column)
        std::fprintf(output, ", \"gradB_%d%d_T_per_m\"", row, column);
  if (layout.velocityGradientOffset != SEP3D::RuntimeModel::kNoOffset)
    for (int row = 0; row < 3; ++row)
      for (int column = 0; column < 3; ++column)
        std::fprintf(output, ", \"gradU_%d%d_s-1\"", row, column);
  // BuildLayout reserves directional wave variance for every authority.  The
  // six public columns are therefore mandatory and include both the requested
  // total turbulence wave-energy density and its field-aligned partition.
  // Keeping the variable fragment in the AMPS-independent output module lets
  // the standalone suite verify the exact output contract.
  if (layout.waveEnergyOffset == SEP3D::RuntimeModel::kNoOffset) {
    StopWithStatus("initialization turbulence output layout",
        SEP3D::Core::Status(SEP3D::Core::StatusCode::LayoutMismatch,
            "mandatory directional turbulence storage is absent"));
  }
  std::fprintf(output, "%s",
               SEP3D::Output::TurbulenceTecplotVariableList());
  // These flags make the two independent empty-data cases machine-readable:
  // background_valid=0 identifies padding cells outside the heliocentric
  // shell, while particle_sample_present=0 identifies a cell/species with no
  // sampled macroparticles.  particle_sampling_window_valid=0 additionally
  // identifies initialization output written before the first sample window.
  std::fprintf(output,
      ", \"background_valid\""
      ", \"particle_sampling_window_valid\""
      ", \"particle_sample_present\"");
}

void PrintInitializationCellData(
    FILE* output, int dataSetNumber, CMPI_channel* pipe,
    int centerNodeThread, PIC::Mesh::cDataCenterNode* centerNode) {
  (void)dataSetNumber;
  const auto& layout = Configuration().storage_layout();
  // Seventeen mandatory background values, six mandatory turbulence values,
  // and three final validity flags.
  // Cells outside the declared heliocentric shell are allocated by the
  // enclosing Cartesian AMR cube but do not represent physical background
  // samples. Empty particle samples are valid and are marked independently.
  std::size_t valueCount = 26;
  if (layout.magneticGradientOffset != SEP3D::RuntimeModel::kNoOffset)
    valueCount += 9;
  if (layout.velocityGradientOffset != SEP3D::RuntimeModel::kNoOffset)
    valueCount += 9;
  std::vector<double> values(valueCount, 0.0);

  const bool ownsNode = pipe == nullptr || pipe->ThisThread == centerNodeThread;
  if (ownsNode) {
    std::size_t cursor = 0;
    auto append = [&](std::size_t offset, std::size_t count) {
      LoadBytes(centerNode, offset, values.data() + cursor,
                count * sizeof(double));
      cursor += count;
    };
    append(layout.magneticFieldOffset, 3);
    append(layout.bulkVelocityOffset, 3);
    append(layout.numberDensityOffset, 1);
    append(layout.velocityDivergenceOffset, 1);
    append(layout.temperatureOffset, 1);
    append(layout.pressureOffset, 1);
    append(layout.alfvenSpeedOffset, 1);
    append(layout.divBhatOffset, 1);
    append(layout.focusingLengthOffset, 1);
    append(layout.curvatureOffset, 3);
    append(layout.fieldAlignedStrainOffset, 1);
    if (layout.magneticGradientOffset != SEP3D::RuntimeModel::kNoOffset)
      append(layout.magneticGradientOffset, 9);
    if (layout.velocityGradientOffset != SEP3D::RuntimeModel::kNoOffset)
      append(layout.velocityGradientOffset, 9);
    if (layout.waveEnergyOffset == SEP3D::RuntimeModel::kNoOffset) {
      StopWithStatus("initialization turbulence cell output",
          SEP3D::Core::Status(SEP3D::Core::StatusCode::LayoutMismatch,
              "mandatory directional turbulence storage is absent"));
    }
    double variance[2] = {};
    LoadBytes(centerNode, layout.waveEnergyOffset, variance, sizeof(variance));
    SEP3D::Output::TurbulenceTecplotPresentation turbulence;
    const SEP3D::Core::Status turbulenceStatus =
        SEP3D::Output::PrepareTurbulenceTecplotPresentation(
            variance[0], variance[1], &turbulence);
    if (!turbulenceStatus.ok())
      StopWithStatus("initialization turbulence cell output",
                     turbulenceStatus);
    values[cursor++] = turbulence.deltaB2T2;
    values[cursor++] = turbulence.deltaBPlus2T2;
    values[cursor++] = turbulence.deltaBMinus2T2;
    values[cursor++] = turbulence.waveEnergyJPerM3;
    values[cursor++] = turbulence.waveEnergyPlusJPerM3;
    values[cursor++] = turbulence.waveEnergyMinusJPerM3;

    double position[3] = {};
    centerNode->GetX(position);
    const SEP3D::Core::Vec3 relative(
        position[0] - Configuration().options().coordinateOriginM.x,
        position[1] - Configuration().options().coordinateOriginM.y,
        position[2] - Configuration().options().coordinateOriginM.z);
    const double radius = relative.Norm();
    const bool insidePhysicalShell =
        radius >= Configuration().options().innerRadiusM &&
        radius <= Configuration().options().outerRadiusM;

    // AMPS' native particle sampler already returns finite zeros for weighted
    // moments when the total sampled weight is zero.  Read its independently
    // sampled particle count only to publish explicit availability flags; the
    // background fields never depend on particle occupancy.
    const double sampledParticleNumber = centerNode->GetDatumAverage(
        PIC::Mesh::DatumParticleNumber, dataSetNumber);
    std::vector<double> storedBackground(values.begin(),
                                         values.begin() + cursor);
    const SEP3D::Output::TecplotCellPresentation presentation =
        SEP3D::Output::PrepareTecplotCellPresentation(
            storedBackground, insidePhysicalShell, PIC::LastSampleLength,
            sampledParticleNumber);
    std::copy(presentation.backgroundValues.begin(),
              presentation.backgroundValues.end(), values.begin());
    values[cursor++] = presentation.backgroundValid;
    values[cursor++] = presentation.particleSamplingWindowValid;
    values[cursor++] = presentation.particleSamplePresent;
  }

  if (PIC::ThisThread == 0 || pipe == nullptr) {
    if (centerNodeThread != 0 && pipe != nullptr)
      pipe->recv(values.data(), static_cast<int>(values.size()),
                 centerNodeThread);
    for (double value : values) std::fprintf(output, "%e ", value);
  } else {
    pipe->send(values.data(), static_cast<int>(values.size()));
  }
}

SEP3D::Core::Status ResolveLocalTransportImpl(
    const SEP3D::Core::Vec3& positionM, int species,
    double momentumKgMPerS, double mu,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node,
    SEP3D::Adapters::LocalTransportRecord* local,
    bool evaluateParallelGradient) {
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
  if (layout.magneticGradientOffset != RuntimeModel::kNoOffset)
    LoadBytes(cell, layout.magneticGradientOffset, background.gradB.m,
              sizeof(background.gradB.m));
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
  if (layout.waveEnergyOffset == RuntimeModel::kNoOffset)
    return Core::Status(Core::StatusCode::LayoutMismatch,
                        "mandatory directional turbulence storage is absent");
  double variance[2];
  LoadBytes(cell, layout.waveEnergyOffset, variance, sizeof(variance));
  waves.deltaBPlus2T2 = variance[0];
  waves.deltaBMinus2T2 = variance[1];
  waves.deltaB2T2 = variance[0] + variance[1];
  if (!waves.status.usable() || !waves.valid) return waves.status;
  const Turbulence::LocalScatteringCoefficients coefficients =
      Turbulence::EvaluateLocalScattering(
          waves, background, positionM, species,
          PIC::MolecularData::GetMass(species),
          PIC::MolecularData::GetElectricCharge(species),
          momentumKgMPerS, mu);
  if (!coefficients.status.ok()) return coefficients.status;
  local->kappaParallelM2PerS = coefficients.kappaParallelM2PerS;
  local->dMuMuPerS = coefficients.dMuMuPerS;
  local->dDmuMuDmuPerS = coefficients.dDmuMuDmuPerS;

  // The Parker Ito drift requires b-hat dot grad(kappa_parallel), not merely
  // kappa itself.  Evaluate the same provider/coefficient chain one local-cell
  // spacing in both field-aligned directions.  The spacing is enlarged by the
  // largest b component so at least one Cartesian coordinate crosses a cell
  // centre spacing even when the field is oblique to the AMR axes.
  if (evaluateParallelGradient) {
    // Do not prefill the production derivative with zero.  A successful
    // top-level resolution must write a value through the validated stencil;
    // every failure returns before the partially filled record can reach a
    // mover.  Recursive neighbour records need only kappa_parallel itself and
    // therefore deliberately skip this field.
    const double maximumDirectionComponent = std::max(
        std::fabs(background.bHat.x),
        std::max(std::fabs(background.bHat.y),
                 std::fabs(background.bHat.z)));
    if (!(maximumDirectionComponent > 0.0) ||
        !std::isfinite(maximumDirectionComponent)) {
      return Core::Status(Core::StatusCode::BackgroundInvalid,
                          "parallel-kappa stencil has invalid magnetic direction");
    }
    const double stencilStepM = local->cellSizeM / maximumDirectionComponent;
    const Core::Vec3 displacement = background.bHat * stencilStepM;
    const Core::Vec3 minusPosition = positionM - displacement;
    const Core::Vec3 plusPosition = positionM + displacement;

    auto nodeAt = [node](const Core::Vec3& position) {
      double coordinates[3];
      position.CopyTo(coordinates);
      return PIC::Mesh::mesh->findTreeNode(coordinates, node);
    };
    Adapters::LocalTransportRecord minus;
    Adapters::LocalTransportRecord plus;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* minusNode = nodeAt(minusPosition);
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* plusNode = nodeAt(plusPosition);
    const Core::Status minusStatus = ResolveLocalTransportImpl(
        minusPosition, species, momentumKgMPerS, mu, minusNode, &minus, false);
    const Core::Status plusStatus = ResolveLocalTransportImpl(
        plusPosition, species, momentumKgMPerS, mu, plusNode, &plus, false);

    Turbulence::ParallelKappaGradientStencil stencil;
    stencil.centerKappaM2PerS = local->kappaParallelM2PerS;
    stencil.stepM = stencilStepM;
    stencil.hasMinus = minusStatus.ok();
    stencil.minusKappaM2PerS = minus.kappaParallelM2PerS;
    stencil.hasPlus = plusStatus.ok();
    stencil.plusKappaM2PerS = plus.kappaParallelM2PerS;
    const Core::Status gradientStatus =
        Turbulence::EvaluateParallelKappaGradient(
            stencil, &local->dKappaParallelDsMPerS);
    if (!gradientStatus.ok()) {
      // A particle for which neither neighbour is locally resolvable cannot be
      // advanced with a fabricated zero drift.  Return the explicit stencil
      // failure; AMPS will classify the particle through its normal fail-closed
      // mover path and retain the diagnostic in the particle ledger.
      return gradientStatus;
    }
  }
  local->fractionalFieldVariationPerS = std::fabs(background.divU);
  const RuntimeModel::SnapshotDescriptor* active =
      ApplicationRuntime().active_snapshot();
  local->timeToSnapshotBoundaryS = active == nullptr ? 0.0 :
      std::max(0.0, active->validUntilS - ApplicationRuntime().CurrentTimeS());
  return Core::Status::OK();
}

SEP3D::Core::Status ResolveLocalTransport(
    const SEP3D::Core::Vec3& positionM, int species,
    double momentumKgMPerS, double mu,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node,
    SEP3D::Adapters::LocalTransportRecord* local) {
  // Only the top-level resolver constructs a derivative stencil.  Its
  // neighbour evaluations call the implementation with the flag cleared, so
  // the work is bounded to two additional coefficient evaluations rather than
  // recursively expanding across the mesh.
  return ResolveLocalTransportImpl(positionM, species, momentumKgMPerS, mu,
                                   node, local, true);
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

// Give every owner-local interior cell a deterministic application state
// before filling the physical heliocentric shell.  AMPS allocates a Cartesian
// cube, while srcSEP3D's physical background occupies only the configured
// spherical shell.  Explicitly zeroing the complete frozen slice makes the
// nonphysical padding well-defined and, critically, gives the Tecplot
// interpolation callback finite source values when a vertex stencil straddles
// either spherical boundary.
void ZeroApplicationStateOnOwnedCells() {
  const std::size_t bytes =
      Configuration().storage_layout().cellAssociatedBytes;
  if (gStaticCellDataOffset < 0 || bytes == 0) {
    StopWithStatus("srcSEP3D static center-node initialization",
        SEP3D::Core::Status(SEP3D::Core::StatusCode::LayoutMismatch,
            "the frozen application cell layout is not allocated"));
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
          PIC::Mesh::cDataCenterNode* cell = node->block->GetCenterNode(
              PIC::Mesh::mesh->getCenterNodeLocalNumber(i, j, k));
          if (cell == nullptr) continue;
          std::memset(cell->GetAssociatedDataBufferPointer() +
                          gStaticCellDataOffset,
                      0, bytes);
        }
  }
}

// Store the directional magnetic variances used by the scattering kernel in
// the application-owned center-node slice and immediately read them back.
// This function is intentionally called only after the selected turbulence
// provider has been prepared and before the associated-data halo exchange.
// Thus one initialized value is authoritative for movers, interpolation, and
// the initialization Tecplot product.
double StoreTurbulenceAtCellCenter(
    PIC::Mesh::cDataCenterNode* cell,
    const SEP3D::Turbulence::TurbulenceSample& waves) {
  using namespace SEP3D;
  const std::size_t offset = Configuration().storage_layout().waveEnergyOffset;
  if (cell == nullptr || offset == RuntimeModel::kNoOffset) {
    StopWithStatus("turbulence cell storage",
        Core::Status(Core::StatusCode::LayoutMismatch,
            "mandatory directional turbulence center-node storage is absent"));
  }
  const double variance[2] = {
      waves.deltaBPlus2T2, waves.deltaBMinus2T2};
  const double total = variance[0] + variance[1];
  if (!std::isfinite(variance[0]) || !std::isfinite(variance[1]) ||
      !std::isfinite(total) || variance[0] < 0.0 || variance[1] < 0.0) {
    StopWithStatus("turbulence cell storage",
        Core::Status(Core::StatusCode::InvalidInput,
            "provider returned negative or non-finite directional variance"));
  }
  const double comparisonScale = std::max(
      std::numeric_limits<double>::min(),
      std::max(std::fabs(total), std::fabs(waves.deltaB2T2)));
  if (!std::isfinite(waves.deltaB2T2) ||
      std::fabs(total - waves.deltaB2T2) >
          128.0 * std::numeric_limits<double>::epsilon() * comparisonScale) {
    StopWithStatus("turbulence cell storage",
        Core::Status(Core::StatusCode::ConfigurationConflict,
            "provider total variance disagrees with its directional partition"));
  }
  // Every accepted prescribed amplitude law has a strictly positive
  // normalization and Parker/SWMF background validation requires |B|>0.
  // Therefore an ordinary prescribed record cannot physically be zero.  A
  // zero state remains legal only for an explicitly selected coupled
  // ballistic/missing-data record, which carries waves.ballistic=true.
  if (Configuration().options().turbulence ==
          RuntimeModel::TurbulenceAuthority::Prescribed &&
      !waves.ballistic && !(total > 0.0)) {
    StopWithStatus("turbulence cell storage",
        Core::Status(Core::StatusCode::BackgroundInvalid,
            "prescribed turbulence evaluated to zero magnetic variance"));
  }

  StoreBytes(cell, offset, variance, sizeof(variance));
  double readback[2] = {};
  LoadBytes(cell, offset, readback, sizeof(readback));
  if (readback[0] != variance[0] || readback[1] != variance[1]) {
    StopWithStatus("turbulence cell storage",
        Core::Status(Core::StatusCode::LayoutMismatch,
            "directional turbulence variance failed center-node readback"));
  }
  return total;
}

// ---------------------------------------------------------------------------
// Native AMPS background bridge
//
// srcSEP3D keeps a complete, versioned BackgroundSample in application-owned
// associated data because that is the mover/restart contract.  AMPS' native
// Tecplot callback, however, does not read those offsets: in DATAFILE coupler
// builds it reads PIC::CPLR::DATAFILE's independent center-node region.  A
// Parker snapshot can therefore be fully initialized while the native Bx/By/
// Bz columns still contain allocation-time zeros.  The bridge below copies the
// *same validated sample* into that native region; it does not evaluate a
// second model or invent a second set of plasma parameters.
// ---------------------------------------------------------------------------

#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__DATAFILE_

void ValidateNativeAmpsBackgroundLayout() {
  using namespace PIC::CPLR::DATAFILE;
  if (CenterNodeAssociatedDataOffsetBegin < 0 ||
      MULTIFILE::CurrDataFileOffset < 0 ||
      nTotalBackgroundVariables <= 0) {
    StopWithStatus("native AMPS background layout",
        SEP3D::Core::Status(SEP3D::Core::StatusCode::LayoutMismatch,
            "DATAFILE center-node storage was not allocated before the "
            "srcSEP3D background fill"));
  }

  // BackgroundSample represents one canonical solar-wind state.  Replicating
  // it into an arbitrary number of AMPS ion-fluid slots would silently assign
  // identities/compositions that are not present in the input contract.  The
  // current standalone interface therefore supports the one-fluid DATAFILE
  // layout used by srcSEP3D and fails closed for a different AMPS build.
  if (nIonFluids != 1) {
    StopWithStatus("native AMPS background layout",
        SEP3D::Core::Status(
            SEP3D::Core::StatusCode::ConfigurationConflict,
            "srcSEP3D has one canonical solar-wind state but the AMPS "
            "DATAFILE layout declares " + std::to_string(nIonFluids) +
            " ion fluids; an explicit fluid-to-physics mapping is required"));
  }
}

// Return the current DATAFILE storage slot and, when AMPS has initialized a
// distinct time-interpolation slot, that slot as well.  Initializing both
// avoids a later interpolation between a valid Parker state and uninitialized
// bytes.  DATAFILE leaves NextDataFileOffset negative when interpolation is
// not in use, so no speculative schedule or offset is constructed here.
int NativeAmpsDataSlots(int slots[2]) {
  slots[0] = PIC::CPLR::DATAFILE::MULTIFILE::CurrDataFileOffset;
  const int next = PIC::CPLR::DATAFILE::MULTIFILE::NextDataFileOffset;
  if (next >= 0 && next != slots[0]) {
    slots[1] = next;
    return 2;
  }
  return 1;
}

void StoreNativeAmpsField(
    PIC::Mesh::cDataCenterNode* cell,
    const PIC::CPLR::DATAFILE::cOffsetElement& field,
    const double* values, int valueCount) {
  // An unallocated optional field has no storage and is intentionally skipped.
  // An allocated field with an invalid offset is a frozen-layout corruption,
  // not a condition that may be hidden by omitting a Tecplot column.
  if (!field.allocate) return;
  if (field.RelativeOffset < 0 || field.nVars != valueCount) {
    StopWithStatus("native AMPS background field",
        SEP3D::Core::Status(SEP3D::Core::StatusCode::LayoutMismatch,
            "allocated DATAFILE field '" + std::string(field.VarList) +
            "' has an invalid offset or component count"));
  }

  int slots[2] = {};
  const int slotCount = NativeAmpsDataSlots(slots);
  for (int slotIndex = 0; slotIndex < slotCount; ++slotIndex) {
    char* destination = cell->GetAssociatedDataBufferPointer() +
        PIC::CPLR::DATAFILE::CenterNodeAssociatedDataOffsetBegin +
        slots[slotIndex] + field.RelativeOffset;
    std::memcpy(destination, values,
                static_cast<std::size_t>(valueCount) * sizeof(double));
  }
}

void ZeroNativeAmpsBackgroundOnOwnedCells() {
  ValidateNativeAmpsBackgroundLayout();
  int slots[2] = {};
  const int slotCount = NativeAmpsDataSlots(slots);
  const std::size_t bytes = static_cast<std::size_t>(
      PIC::CPLR::DATAFILE::nTotalBackgroundVariables) * sizeof(double);

  // AMPS allocates a Cartesian cube around the physical heliocentric shell.
  // Zero every owner-local interior cell first so padding cells have a finite,
  // deterministic placeholder.  Physical cells are overwritten below and are
  // distinguished from padding by srcSEP3D's background_valid column.
  for (unsigned int blockIndex = 0;
       blockIndex < PIC::DomainBlockDecomposition::nLocalBlocks;
       ++blockIndex) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
        PIC::DomainBlockDecomposition::BlockTable[blockIndex];
    if (node == nullptr || node->block == nullptr) continue;
    for (int k = 0; k < _BLOCK_CELLS_Z_; ++k)
      for (int j = 0; j < _BLOCK_CELLS_Y_; ++j)
        for (int i = 0; i < _BLOCK_CELLS_X_; ++i) {
          PIC::Mesh::cDataCenterNode* cell = node->block->GetCenterNode(
              PIC::Mesh::mesh->getCenterNodeLocalNumber(i, j, k));
          if (cell == nullptr) continue;
          for (int slotIndex = 0; slotIndex < slotCount; ++slotIndex) {
            char* destination = cell->GetAssociatedDataBufferPointer() +
                PIC::CPLR::DATAFILE::CenterNodeAssociatedDataOffsetBegin +
                slots[slotIndex];
            std::memset(destination, 0, bytes);
          }
        }
  }
}

void StoreNativeAmpsBackground(
    PIC::Mesh::cDataCenterNode* cell,
    const SEP3D::Background::BackgroundSample& sample) {
  namespace Offset = PIC::CPLR::DATAFILE::Offset;

  const double magnetic[3] = {sample.B.x, sample.B.y, sample.B.z};
  const double velocity[3] = {sample.U.x, sample.U.y, sample.U.z};

  // Both currently implemented background authorities describe an ideal-MHD
  // solar wind.  The electric field is therefore the SI motional field
  // E=-U x B.  This is derived from the already validated U and B rather than
  // being supplied as an independent, potentially inconsistent input.
  const double electric[3] = {
      -(sample.U.y * sample.B.z - sample.U.z * sample.B.y),
      -(sample.U.z * sample.B.x - sample.U.x * sample.B.z),
      -(sample.U.x * sample.B.y - sample.U.y * sample.B.x)};

  // The native field named PlasmaNumberDensity mirrors BackgroundSample's
  // documented electron density; PlasmaTemperature mirrors proton
  // temperature; and PlasmaIonPressure mirrors the canonical total thermal
  // pressure.  These exact meanings are stated in BACKGROUND_FIELD.md and are
  // also exposed by the unit-bearing srcSEP3D columns in the same file.
  StoreNativeAmpsField(cell, Offset::PlasmaNumberDensity,
                       &sample.numberDensityM3, 1);
  StoreNativeAmpsField(cell, Offset::PlasmaBulkVelocity, velocity, 3);
  StoreNativeAmpsField(cell, Offset::PlasmaTemperature,
                       &sample.temperatureK, 1);
  StoreNativeAmpsField(cell, Offset::PlasmaIonPressure,
                       &sample.pressurePa, 1);
  StoreNativeAmpsField(cell, Offset::PlasmaDivU, &sample.divU, 1);
  StoreNativeAmpsField(cell, Offset::MagneticField, magnetic, 3);
  StoreNativeAmpsField(cell, Offset::ElectricField, electric, 3);
  StoreNativeAmpsField(cell, Offset::MagneticFieldGradient,
                       &sample.gradB.m[0][0], 9);

  // Relativistic-GCA builds allocate current explicitly.  It is determined
  // without finite differencing from the provider's analytic gradient via
  // Ampere's law (displacement current is absent in the stationary Parker
  // initialization): J = curl(B)/mu0.
  const double current[3] = {
      (sample.gradB.m[2][1] - sample.gradB.m[1][2]) /
          kMagneticPermeabilityVacuum,
      (sample.gradB.m[0][2] - sample.gradB.m[2][0]) /
          kMagneticPermeabilityVacuum,
      (sample.gradB.m[1][0] - sample.gradB.m[0][1]) /
          kMagneticPermeabilityVacuum};
  StoreNativeAmpsField(cell, Offset::Current, current, 3);

  // Electron pressure has a separate native slot only in selected AMPS
  // readers.  Populate it from the resolved SWCME closure when present.  In
  // proton-only closure the canonical pressure deliberately excludes an
  // electron component; in multi-species closure ne*kB*Te is exact.
  double electronPressurePa = 0.0;
  const SEP3D::RuntimeModel::ParkerPhysicsOptions& parker =
      Configuration().options().parker;
  if (parker.thermodynamicClosure ==
      SEP3D::RuntimeModel::SolarWindThermodynamicClosure::MultiSpecies) {
    electronPressurePa = sample.numberDensityM3 * SEP3D::Core::Const::k_B *
        parker.electronTemperatureK;
  }
  StoreNativeAmpsField(cell, Offset::PlasmaElectronPressure,
                       &electronPressurePa, 1);
}

void CompleteNativeAmpsBackgroundInstallation() {
  // AMPS' center-to-corner interpolation can consume neighboring ghost cells
  // while producing output.  Exchange the complete associated-data buffer
  // only after all owner cells contain background and turbulence values, so a
  // rank boundary cannot introduce zero seams into the Tecplot product.
  PIC::Mesh::mesh->ParallelBlockDataExchange();

#if _PIC_MOVER_INTEGRATOR_MODE_ == \
    _PIC_MOVER_INTEGRATOR_MODE__RELATIVISTIC_GCA_
  // These higher-order relativistic-GCA quantities depend on neighboring B/E
  // values.  Generate them only after the first halo exchange, then publish
  // the derived values to ghosts with a second exchange.
  for (unsigned int blockIndex = 0;
       blockIndex < PIC::DomainBlockDecomposition::nLocalBlocks;
       ++blockIndex) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
        PIC::DomainBlockDecomposition::BlockTable[blockIndex];
    if (node != nullptr && node->block != nullptr)
      PIC::CPLR::DATAFILE::GenerateVarForRelativisticGCA(node);
  }
  PIC::Mesh::mesh->ParallelBlockDataExchange();
#endif

  gNativeAmpsBackgroundReady = true;
}

#else

// SWMF-coupler builds own their native AMPS field buffers outside DATAFILE.
// srcSEP3D still fills its application cache and publishes the immutable
// snapshot, while this no-op keeps the initialization boundary uniform.
void ZeroNativeAmpsBackgroundOnOwnedCells() {}
void StoreNativeAmpsBackground(
    PIC::Mesh::cDataCenterNode*,
    const SEP3D::Background::BackgroundSample&) {}
void CompleteNativeAmpsBackgroundInstallation() {
  PIC::Mesh::mesh->ParallelBlockDataExchange();
  gNativeAmpsBackgroundReady = true;
}

#endif

std::shared_ptr<SEP3D::Turbulence::TurbulenceProvider>
MakePrescribedTurbulence() {
  const auto& options = Configuration().options();
  SEP3D::Turbulence::PrescribedKolmogorovConfiguration model;
  switch (options.prescribedTurbulenceModel) {
    case SEP3D::RuntimeModel::PrescribedTurbulenceModel::PowerLaw:
      model.spectrumModel =
          SEP3D::Turbulence::PrescribedSpectrumModel::PowerLaw;
      break;
    case SEP3D::RuntimeModel::PrescribedTurbulenceModel::Kolmogorov:
      model.spectrumModel =
          SEP3D::Turbulence::PrescribedSpectrumModel::Kolmogorov;
      break;
    case SEP3D::RuntimeModel::PrescribedTurbulenceModel::Kraichnan:
      model.spectrumModel =
          SEP3D::Turbulence::PrescribedSpectrumModel::Kraichnan;
      break;
  }
  switch (options.prescribedTurbulenceAmplitudeModel) {
    case SEP3D::RuntimeModel::PrescribedTurbulenceAmplitudeModel::
        ConstantDeltaBOverB:
      model.amplitudeModel =
          SEP3D::Turbulence::PrescribedAmplitudeModel::ConstantDeltaBOverB;
      break;
    case SEP3D::RuntimeModel::PrescribedTurbulenceAmplitudeModel::
        WaveEnergyPowerLaw:
      model.amplitudeModel =
          SEP3D::Turbulence::PrescribedAmplitudeModel::WaveEnergyPowerLaw;
      break;
  }
  model.deltaBOverB = options.prescribedDeltaBOverB;
  model.waveEnergyAtReferenceJPerM3 =
      options.turbulenceWaveEnergyAtReferenceJPerM3;
  model.waveEnergyRadialExponent =
      options.turbulenceWaveEnergyRadialExponent;
  model.normalizedCrossHelicity =
      options.turbulenceNormalizedCrossHelicity;
  model.referenceRadiusM = options.turbulenceReferenceRadiusM;
  model.kMinAtReferencePerM = options.turbulenceKMinPerM;
  model.kMaxAtReferencePerM = options.turbulenceKMaxPerM;
  model.kMinRadialExponent = options.turbulenceKMinRadialExponent;
  model.kMaxRadialExponent = options.turbulenceKMaxRadialExponent;
  model.spectralIndex = options.turbulenceSpectralIndex;
  model.parallelCorrelationLengthAtReferenceM =
      options.turbulenceCorrelationLengthM;
  model.correlationLengthRadialExponent =
      options.turbulenceCorrelationLengthRadialExponent;
  model.validityCadenceS = options.turbulenceValidityCadenceS;
  model.coordinateFrame = options.coordinateFrame;
  return std::shared_ptr<SEP3D::Turbulence::TurbulenceProvider>(
      new SEP3D::Turbulence::PrescribedKolmogorovProvider(model));
}

void FillAndPublishBackground() {
  using namespace SEP3D;
  gNativeAmpsBackgroundReady = false;
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
    parker.densityReferenceRadiusM = configured.densityReferenceRadiusM;
    parker.temperatureK = configured.temperatureK;
    parker.adiabaticIndex = configured.adiabaticIndex;
    parker.thermodynamicClosure =
        configured.thermodynamicClosure ==
                RuntimeModel::SolarWindThermodynamicClosure::MultiSpecies
            ? swcme::solarwind::ThermodynamicClosure::MultiSpecies
            : swcme::solarwind::ThermodynamicClosure::ProtonOnly;
    parker.alphaToProtonRatio = configured.alphaToProtonRatio;
    parker.electronTemperatureK = configured.electronTemperatureK;
    parker.alphaTemperatureK = configured.alphaTemperatureK;
    parker.referenceSinColatitude = configured.referenceSinColatitude;
    parker.solarWindSpeedMPerS = configured.solarWindSpeedMPerS;
    parker.solarRotationRateRadPerS = configured.solarRotationRateRadPerS;
    parker.rotationAxis = configured.rotationAxis;
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

  Core::Status status;
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

  // The application slice and AMPS' DATAFILE slice are independent regions of
  // each center node.  Clear both before one joined physical-cell pass.  That
  // pass is the correct initialization boundary: both providers are prepared,
  // but neither the Runtime snapshot nor the halo-visible AMPS state has yet
  // been advertised as complete.
  ZeroApplicationStateOnOwnedCells();
  ZeroNativeAmpsBackgroundOnOwnedCells();
  gCellSampleIndex.clear();

  unsigned long long localCellCount = 0;
  unsigned long long localPositiveTurbulenceCount = 0;
  double localMinimumVariance = std::numeric_limits<double>::infinity();
  double localMaximumVariance = 0.0;
  for (std::size_t i = 0; i < cells.size(); ++i) {
    if (!SamePosition(snapshot->positions()[i], cells[i].positionM)) {
      StopWithStatus("background cell mapping", Core::Status(
          Core::StatusCode::LayoutMismatch,
          "snapshot positions are not in deterministic owner-cell order"));
    }
    const Turbulence::TurbulenceSample waves = turbulence->Evaluate(
        cells[i].positionM, snapshot->samples()[i]);
    if (!waves.status.usable() || !waves.valid)
      StopWithStatus("turbulence cell evaluation", waves.status);

    // Prescribe background and turbulence to the same physical center-node
    // object before it can be consumed by either movers or Tecplot.  Keeping
    // these writes adjacent eliminates the former state in which the native
    // plasma/IMF slice was populated but the application turbulence slice was
    // still absent from AMPS' output interpolation path.
    StoreBackground(cells[i].cell, snapshot->samples()[i]);
    StoreNativeAmpsBackground(cells[i].cell, snapshot->samples()[i]);
    const double totalVariance =
        StoreTurbulenceAtCellCenter(cells[i].cell, waves);
    gCellSampleIndex[cells[i].cell] = i;
    ++localCellCount;
    if (totalVariance > 0.0) {
      ++localPositiveTurbulenceCount;
      localMinimumVariance = std::min(localMinimumVariance, totalVariance);
      localMaximumVariance = std::max(localMaximumVariance, totalVariance);
    }
  }

  // Make a zero-filled prescribed mesh impossible to misreport as a completed
  // initialization.  The count and range also provide a concise run-time
  // diagnostic that can be compared with the Tecplot columns.  Coupled AWSoM
  // may explicitly select ballistic missing-data handling; those zero cells
  // remain visible through the positive/total count instead of being replaced
  // by an invented amplitude.
  unsigned long long globalCellCount = 0;
  unsigned long long globalPositiveTurbulenceCount = 0;
  double globalMinimumVariance = 0.0;
  double globalMaximumVariance = 0.0;
  MPI_Allreduce(&localCellCount, &globalCellCount, 1,
                MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&localPositiveTurbulenceCount,
                &globalPositiveTurbulenceCount, 1,
                MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&localMinimumVariance, &globalMinimumVariance, 1, MPI_DOUBLE,
                MPI_MIN, MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&localMaximumVariance, &globalMaximumVariance, 1, MPI_DOUBLE,
                MPI_MAX, MPI_GLOBAL_COMMUNICATOR);
  if (Configuration().options().turbulence ==
          RuntimeModel::TurbulenceAuthority::Prescribed &&
      globalPositiveTurbulenceCount != globalCellCount) {
    StopWithStatus("turbulence center-node initialization",
        Core::Status(Core::StatusCode::BackgroundInvalid,
            "not every physical center node received positive prescribed "
            "turbulence variance"));
  }
  if (PIC::ThisThread == 0) {
    std::cout << "[srcSEP3D] turbulence center-node initialization: physical_cells="
              << globalCellCount << " positive_variance_cells="
              << globalPositiveTurbulenceCount;
    if (globalPositiveTurbulenceCount != 0)
      std::cout << " deltaB2_min_T2=" << std::scientific
                << globalMinimumVariance << " deltaB2_max_T2="
                << globalMaximumVariance << std::defaultfloat;
    std::cout << '\n';
  }

  // Publish only after the complete background/turbulence state has survived
  // provider validation, center-node write/readback, and collective coverage
  // checks.  The subsequent halo exchange then makes precisely this published
  // generation available to neighboring AMPS blocks and output interpolation.
  RuntimeModel::Runtime& runtime = ApplicationRuntime();
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

  // This collective halo exchange is the final background-installation
  // boundary.  The data-bearing initialization file is deliberately emitted
  // only after it returns on every rank.
  CompleteNativeAmpsBackgroundInstallation();
}

void WriteInitializationDataTecplotAfterBackground() {
  if (Configuration().options().inputSchemaVersion < 3) return;

  // Treat the file location in the initialization sequence as a contract:
  // geometry may be written after AMR construction, but this data-bearing
  // product may be written only after the immutable snapshot, turbulence,
  // srcSEP3D cell cache, native AMPS DATAFILE cache, species weights, and time
  // steps all describe the same completed initialization state.
  const SEP3D::RuntimeModel::SnapshotDescriptor* active =
      SEP3D::ApplicationRuntime().active_snapshot();
  if (!gInstalledBackground || !gInstalledTurbulence || active == nullptr ||
      !active->complete || !gNativeAmpsBackgroundReady) {
    StopWithStatus("initialization data output",
        SEP3D::Core::Status(
            SEP3D::Core::StatusCode::SnapshotUnavailable,
            "sep3d-initialization-data.dat was requested before the complete "
            "background/turbulence/native-AMPS installation boundary"));
  }
  if (!PIC::ParticleWeightTimeStep::GlobalTimeStepInitialized) {
    StopWithStatus("initialization data output",
        SEP3D::Core::Status(SEP3D::Core::StatusCode::InvalidTransition,
            "AMPS particle time steps are not initialized"));
  }
  for (const auto& species : gCompiledSpecies) {
    const double timeStep =
        PIC::ParticleWeightTimeStep::GlobalTimeStep[species.ampsIndex];
    const double weight =
        PIC::ParticleWeightTimeStep::GlobalParticleWeight[species.ampsIndex];
    if (!std::isfinite(timeStep) || timeStep <= 0.0 ||
        !std::isfinite(weight) || weight <= 0.0) {
      StopWithStatus("initialization data output",
          SEP3D::Core::Status(SEP3D::Core::StatusCode::InvalidTransition,
              "compiled species " + std::to_string(species.ampsIndex) +
              " has no positive finite AMPS time step/particle weight"));
    }
  }

  // CompleteNativeAmpsBackgroundInstallation() is already collective.  This
  // extra barrier makes the output boundary explicit and prevents a fast rank
  // from entering the distributed AMPS writer while another rank is still
  // checking species numerics.
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  for (const auto& species : gCompiledSpecies) {
    const std::string path = InitializationDataPath(
        Configuration().options().initializationDataTecplotFile,
        species.ampsIndex, PIC::nTotalSpecies);
    PIC::Mesh::mesh->outputMeshDataTECPLOT(
        path.c_str(), species.ampsIndex);
  }
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
}

void RefreshBackgroundAtBoundary() {
  using namespace SEP3D;
  RuntimeModel::Runtime& runtime = ApplicationRuntime();
  if (!runtime.EventDue(RuntimeModel::ScheduledEvent::Background)) return;

  gNativeAmpsBackgroundReady = false;

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
    StoreNativeAmpsBackground(cells[index].cell,
                              candidate->samples()[index]);
    // Use the identical validated center-node write/readback path as initial
    // publication so a later analytic/SWMF generation cannot reintroduce the
    // zero-turbulence output defect.
    (void)StoreTurbulenceAtCellCenter(cells[index].cell, waves);
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
  // Publish the committed generation to AMPS ghost cells as one final
  // collective operation.  The native buffer can now be consumed by the next
  // particle phase and by any later AMPS data output.
  CompleteNativeAmpsBackgroundInstallation();
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
    // Register before AMPS freezes its output callback lists. The callbacks
    // read only the immutable storage layout and the completed center-node
    // buffers; they do not introduce a second background representation.
    PIC::Mesh::PrintVariableListCenterNode.push_back(
        PrintInitializationVariableList);
    PIC::Mesh::PrintDataCenterNode.push_back(PrintInitializationCellData);
    // outputMeshDataTECPLOT prints vertex records through temporary
    // center-node objects.  This hook is therefore as essential as the print
    // callback: it copies the initialized application-owned background and
    // turbulence slice into those objects before PrintInitializationCellData
    // reads it.
    PIC::Mesh::InterpolateCenterNode.push_back(
        InterpolateInitializationCellData);
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
  // Capture and validate the complete generated species table immediately
  // after AMPS base initialization.  The table is read-only: its count,
  // symbols, masses, charges, and indices were fixed by SpeciesList when this
  // executable was built and cannot be redefined by the runtime SEP deck.
  BindCompiledSpeciesTable();
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
  if (options.inputSchemaVersion >= 3) {
    // The AMPS writer emits the final distributed octree after refinement and
    // decomposition.  The finite Parker centreline is a deterministic ordered
    // zone and is written once by rank zero from the same ResolutionConfiguration
    // consumed by localResolution().  Failure is fatal because both files are
    // declared initialization products, not optional diagnostics.
    EnsureInitializationOutputParents(options);
    PIC::Mesh::mesh->outputMeshTECPLOT(
        options.initializationMeshTecplotFile.c_str());
    if (PIC::ThisThread == 0) {
      const SEP3D::Core::Status lineOutput =
          SEP3D::Mesh::WriteParkerCenterlineTecplot(
              resolution, options.initializationParkerLineTecplotFile);
      if (!lineOutput.ok())
        StopWithStatus("Parker initialization Tecplot", lineOutput);
    }
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  }
  PIC::Mesh::mesh->AllowBlockAllocation = true;
  PIC::Mesh::mesh->AllocateTreeBlocks();

  // AllocateTreeBlocks() materializes the blocks referenced by the final
  // parallel distribution list, but it does not populate AMPS' cached
  // DomainBlockDecomposition::BlockTable.  Every srcSEP3D owner-local pass
  // below (block particle numerics, background fill, observers, checkpoints,
  // and ledgers) deliberately iterates that cache.  Refresh it here, after
  // allocation and before amps_init(), so a valid nonempty mesh cannot appear
  // as an empty background grid.  PIC::TimeStep() also refreshes this table,
  // but waiting until the first step is too late because initialization must
  // publish a complete Parker snapshot before particle motion is permitted.
  PIC::DomainBlockDecomposition::UpdateBlockTable();

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
  // BindCompiledSpeciesTable() already captured and validated the generated
  // AMPS table in amps_init_mesh().  No molecular-data setter is called: AMPS
  // remains the sole authority for the immutable species identity and physics.
  const SEP3D::RuntimeModel::SpeciesOptions& configuredSpecies =
      Configuration().options().species;

  // R04: one explicitly configured base step and base statistical weight are
  // installed on every compiled species.  Source sampling later creates an
  // individual correction when exact conservation requires it; that does not
  // change the global/block base values initialized here.
  const double configuredDt = Configuration().options().requestedTimeStepS;
  for (const auto& species : gCompiledSpecies) {
    PIC::ParticleWeightTimeStep::GlobalTimeStep[species.ampsIndex] =
        configuredDt;
    PIC::ParticleWeightTimeStep::GlobalParticleWeight[species.ampsIndex] =
        configuredSpecies.macroparticleWeight;
  }
  PIC::ParticleWeightTimeStep::GlobalTimeStepInitialized = true;
  for (unsigned int blockIndex = 0;
       blockIndex < PIC::DomainBlockDecomposition::nLocalBlocks;
       ++blockIndex) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
        PIC::DomainBlockDecomposition::BlockTable[blockIndex];
    if (node == nullptr || node->block == nullptr) continue;
    for (const auto& species : gCompiledSpecies) {
      node->block->SetLocalTimeStep(configuredDt, species.ampsIndex);
      // Keep AMPS' block-local normalization identical to the global base
      // weight for every compiled species.  No slot may retain -1 or a legacy
      // value from a different initialization path.
      node->block->SetLocalParticleWeight(
          configuredSpecies.macroparticleWeight, species.ampsIndex);
    }
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

  // Keep the data-bearing output as the final operation in amps_init().  At
  // this point background/turbulence publication, native AMPS synchronization,
  // mover installation, optional restart restoration, and all species
  // numerical initialization have completed.
  WriteInitializationDataTecplotAfterBackground();
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

    // A valid CME may deliberately start after simulation time zero.  Its
    // provider publishes active=false and no surface patches before
    // event.valid_from; that interval represents zero physical source, not an
    // allocation failure and not N synthetic particles.  Once active, an
    // empty/malformed surface still reaches AllocateExactPatchMacroparticles
    // and fails closed.
    if (shock.active && Configuration().options().source.enabled &&
        runtime.EventDue(SEP3D::RuntimeModel::ScheduledEvent::Injection)) {
      const auto& options = Configuration().options();
      // `samples_per_step` is an exact count per compiled species.  Allocate
      // it independently across the full shock surface so every SpeciesList
      // entry is injected; sharing one count among species would invent an
      // undeclared composition and could omit a low-index or low-mass species.
      for (const auto& compiled : gCompiledSpecies) {
        std::vector<std::uint64_t> exactPatchCounts;
        if (options.inputSchemaVersion >= 3) {
          status = SEP3D::Adapters::AllocateExactPatchMacroparticles(
              shock.patches, options.source.samplesPerStep,
              &exactPatchCounts);
          if (!status.ok())
            StopWithStatus("exact per-species shock source allocation", status);
        }

        for (std::size_t patchIndex = 0;
             patchIndex < shock.patches.size(); ++patchIndex) {
          SEP3D::Adapters::ShockSourceRecord patch =
              shock.patches[patchIndex];
          // SWCME supplies shock geometry and compression.  Convert the input
          // deck's total kinetic-energy interval with this compiled species'
          // immutable AMPS mass; reusing a proton momentum interval for an
          // electron or heavy ion would inject the wrong energies.
          status = SEP3D::Adapters::ConfigureSpeciesSpectrum(
              &patch, compiled.massKg, options.source.minimumEnergyJ,
              options.source.maximumEnergyJ);
          if (!status.ok())
            StopWithStatus("per-species source spectrum", status);

          SEP3D::Adapters::SourceRequest source;
          source.patch = patch;
          source.step = runtime.counters().currentTick;
          source.species = compiled.ampsIndex;
          source.speciesMassKg = compiled.massKg;
          source.intervalS = options.requestedTimeStepS *
              options.injectionCadenceSteps;
          // The declared physical rate is per compiled species.  Apply the
          // common efficiency once, then partition that species' represented
          // population by the normalized physical shock-patch weight.
          source.physicalParticleRatePerS =
              options.source.physicalParticleRatePerS *
              options.source.injectionEfficiency *
              patch.relativePatchWeight;
          source.macroparticleWeight = options.species.macroparticleWeight;
          if (options.inputSchemaVersion >= 3) {
            source.prescribedMacroparticles = exactPatchCounts[patchIndex];
            source.maximumMacroparticles = source.prescribedMacroparticles;
          } else {
            source.maximumMacroparticles = options.source.samplesPerStep;
          }
          SEP3D::Adapters::InjectionPlan plan =
              SEP3D::Adapters::BuildInjectionPlan(source);
          if (!plan.status.ok())
            StopWithStatus("shock source plan", plan.status);
          const SEP3D::AMPS::Movers::InjectionOutcome injected =
              SEP3D::AMPS::Movers::InjectParticles(plan);
          if (!injected.status.ok())
            StopWithStatus("AMPS source allocation", injected.status);
          // All ranks construct the same deterministic plan, but AMPS allocates
          // it only on the rank owning the patch position.  Record the row
          // there so each (step,species,patch) appears exactly once globally.
          if (injected.allocated != 0 ||
              (plan.particles.empty() && PIC::ThisThread == 0))
            gSourceLedger.push_back(plan.ledger);
        }
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

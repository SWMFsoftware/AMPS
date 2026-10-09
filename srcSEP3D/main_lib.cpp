// ============================================================================
// srcSEP3D/main_lib.cpp
//
// AMPS application boundary through Phases R2, M, B, T, P, A, and O.
//
// The production boundary now owns the typed Runtime introduced in R2.  Both
// standalone and coupled hosts install a validated immutable configuration and
// use the same Runtime transitions. This file never parses process arguments;
// in shared-file mode it does own the collective srcSEP3D section parse at the
// required post-Init_BeforeParser boundary. Phase M builds the Cartesian AMR mesh and freezes AMPS
// storage offsets.  Phase B publishes a complete immutable ambient snapshot,
// and Phase T validates/fills the selected turbulence input.  Particle motion
// is dispatched through the single Phase-A AMPS adapter. Phase-O coordinators
// are available to the host only at joined step/checkpoint boundaries.
// ============================================================================

#include "SEP3D.h"

#include "background/bg_parker.h"
#include "background/background_snapshot.h"
#include "adapters/source_runtime.h"
#include "adapters/shock_front_background_adapter.h"
#include "amps/amps_mover_status.h"
#include "amps/amps_particle_adapter.h"
#include "mesh/mesh_model.h"
#include "output/observer_runtime.h"
#include "output/publication.h"
#include "output/restart.h"
#include "output/reduced_front_output.h"
#include "output/sampling.h"
#include "output/shock_history.h"
#include "runtime/runtime_adapters.h"
#include "runtime/background_factory.h"
#include "runtime/application_input.h"
#include "runtime/particle_normalization.h"
#include "turbulence/turbulence_models.h"
#include "validation/coronal_cme_application_test.h"
#include "diagnostics.h"
#include "sep_coronal_cme/constants.h"

#include <algorithm>
#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <list>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// srcSEP3D registers the photosphere through AMPS' internal-boundary manager;
// compiling that manager out would otherwise leave a source-compatible call
// that aborts only at run time.  Fail at compilation with the actual contract.
#if _INTERNAL_BOUNDARY_MODE_ != _INTERNAL_BOUNDARY_MODE_ON_
#error "srcSEP3D requires AMPS internal-boundary support for the solar surface"
#endif
#if _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ != \
    _USER_DEFINED_INTERNAL_BOUNDARY_SPHERE_MODE_ON_
#error "srcSEP3D requires AMPS user-defined sphere callbacks for solar absorption"
#endif

namespace {

namespace fs = std::filesystem;
constexpr double kMagneticPermeabilityVacuum =
    4.0e-7 * SEP3D::Core::Const::kPi;

// Empty for coupled hosts and maintained schema-4 ``--input`` runs.  In the
// new shared-file mode main.cpp installs exactly one path before native
// initialization; the file is deliberately not opened until both AMPS and
// srcSEP3D have completed their Init_BeforeParser hooks.
std::string gApplicationInputPath;
SEP3D::RuntimeModel::Sep3dApplicationInput gParsedApplicationInput;
bool gHasParsedApplicationInput = false;

// Defined below after the process-owned state declarations.  The early input
// transaction uses the same fail-closed accessor as all later native paths.
const SEP3D::RuntimeModel::RunConfiguration3D& Configuration();

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

// Scan a published text product token by token and reject only values that
// parse as floating-point numbers but are non-finite.  Header words containing
// character sequences such as "inf" are not numbers and therefore cannot
// produce a false failure.  The native test calls this only after every rank
// has returned from the collective AMPS writer.
bool FileHasOnlyFiniteNumericTokens(const std::string& path) {
  std::FILE* input = std::fopen(path.c_str(), "r");
  if (input == nullptr) return false;
  char rawToken[1024];
  bool finite = true;
  while (std::fscanf(input, "%1023s", rawToken) == 1) {
    std::string token(rawToken);
    while (!token.empty() &&
           (token.back() == ',' || token.back() == ';')) token.pop_back();
    if (token.empty()) continue;
    char* end = nullptr;
    errno = 0;
    const double value = std::strtod(token.c_str(), &end);
    if (end != token.c_str() && *end == '\0' &&
        (errno == ERANGE || !std::isfinite(value))) {
      finite = false;
      break;
    }
  }
  const bool readSucceeded = std::ferror(input) == 0;
  std::fclose(input);
  return finite && readSucceeded;
}

// Parse the shared-file srcSEP3D section on rank zero and distribute the one
// resolved setting, rather than allowing ranks to see different filesystem
// contents during startup.  The provisional configuration exists only so
// Init_BeforeParser can register its storage callbacks.  The replacement is
// accepted only when its byte layout is identical, and it commits before the
// first AMPS buffer offset is frozen or provider snapshot is constructed.
void ParseInstalledApplicationInput() {
  if (gApplicationInputPath.empty()) return;

  SEP3D::RuntimeModel::Sep3dApplicationInput parsed;
  SEP3D::Core::Status parseStatus = SEP3D::Core::Status::OK();
  if (PIC::ThisThread == 0) {
    parseStatus = SEP3D::RuntimeModel::ParseSep3dApplicationInput(
        gApplicationInputPath, &parsed);
  }

  int parseOk = parseStatus.ok() ? 1 : 0;
  MPI_Bcast(&parseOk, 1, MPI_INT, 0, MPI_GLOBAL_COMMUNICATOR);
  if (parseOk == 0) {
    std::string message = PIC::ThisThread == 0
        ? parseStatus.message : std::string();
    unsigned long long length =
        static_cast<unsigned long long>(message.size());
    MPI_Bcast(&length, 1, MPI_UNSIGNED_LONG_LONG, 0,
              MPI_GLOBAL_COMMUNICATOR);
    if (PIC::ThisThread != 0) message.resize(static_cast<std::size_t>(length));
    if (length != 0) {
      MPI_Bcast(message.data(), static_cast<int>(length), MPI_CHAR, 0,
                MPI_GLOBAL_COMMUNICATOR);
    }
    StopWithStatus("application input parser", SEP3D::Core::Status(
        SEP3D::Core::StatusCode::InvalidInput, message));
  }

  auto broadcastString=[](std::string* value) {
    unsigned long long length=PIC::ThisThread==0?
        static_cast<unsigned long long>(value->size()):0ULL;
    MPI_Bcast(&length,1,MPI_UNSIGNED_LONG_LONG,0,MPI_GLOBAL_COMMUNICATOR);
    if(length>static_cast<unsigned long long>(
          std::numeric_limits<int>::max()))
      StopWithStatus("application input broadcast",SEP3D::Core::Status(
          SEP3D::Core::StatusCode::InvalidInput,
          "application input string exceeds the MPI count range"));
    if(PIC::ThisThread!=0)value->resize(static_cast<std::size_t>(length));
    if(length!=0)MPI_Bcast(value->data(),static_cast<int>(length),MPI_CHAR,0,
        MPI_GLOBAL_COMMUNICATOR);
  };

  unsigned long long integerInputs[3] = {
      PIC::ThisThread == 0
          ? static_cast<unsigned long long>(parsed.particlesPerIteration)
          : 0ULL,
      PIC::ThisThread == 0
          ? static_cast<unsigned long long>(
                parsed.maximumInjectionEventsPerSpeciesPerStep)
          : 0ULL,
      PIC::ThisThread == 0
          ? static_cast<unsigned long long>(parsed.maximumTimeSteps)
          : 0ULL};
  MPI_Bcast(integerInputs, 3, MPI_UNSIGNED_LONG_LONG, 0,
            MPI_GLOBAL_COMMUNICATOR);
  double numericalInputs[6]={parsed.maximumParticleSpeedMPerS,
      parsed.timeStepMarginFactor,parsed.sourceNormalizationRadiusM,
      parsed.minimumInjectionEnergyJ,parsed.maximumInjectionEnergyJ,
      parsed.fixedPhaseSpacePowerIndex};
  MPI_Bcast(numericalInputs,6,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);
  broadcastString(&parsed.shockModel);
  broadcastString(&parsed.backgroundPlasmaModel);
  broadcastString(&parsed.sourceModel);
  broadcastString(&parsed.particleWeightingModel);
  broadcastString(&parsed.momentumPowerLawModel);
  broadcastString(&parsed.reducedShockConfiguration);
  broadcastString(&parsed.reducedShockAssetDirectory);
  if(PIC::ThisThread!=0) {
    parsed.particlesPerIteration=
        static_cast<std::uint64_t>(integerInputs[0]);
    parsed.maximumInjectionEventsPerSpeciesPerStep=
        static_cast<std::uint64_t>(integerInputs[1]);
    parsed.maximumTimeSteps=static_cast<std::uint64_t>(integerInputs[2]);
    parsed.maximumParticleSpeedMPerS=numericalInputs[0];
    parsed.timeStepMarginFactor=numericalInputs[1];
    parsed.sourceNormalizationRadiusM=numericalInputs[2];
    parsed.minimumInjectionEnergyJ=numericalInputs[3];
    parsed.maximumInjectionEnergyJ=numericalInputs[4];
    parsed.fixedPhaseSpacePowerIndex=numericalInputs[5];
  }

  SEP3D::RuntimeModel::RunConfiguration3DOptions options =
      Configuration().options();
  // The shared section owns the selected provider pair. The reduced provider
  // remains the sole front/ambient authority; `shock=none` deliberately keeps
  // the legacy SWCME adapter inactive.  ShockInjection is nevertheless a
  // truthful intent because the new callback samples the same reduced
  // provider's accepted triangular faces rather than constructing a second
  // geometry authority.
  options.inputSchemaVersion=4;
  options.maximumTimeSteps=parsed.maximumTimeSteps;
  options.background=
      SEP3D::RuntimeModel::BackgroundAuthority::RuntimeModel;
  options.backgroundModelId="sep-corona-swcme-shock-front-v1";
  options.backgroundModelAssetPath.clear();
  options.backgroundModelInlineConfiguration=
      parsed.reducedShockConfiguration;
  options.backgroundModelAssetDirectory=parsed.reducedShockAssetDirectory;
  options.coordinateFrame="HCI";
  // The generic configuration keeps an analytic Parker record even when it
  // is not authoritative, because mesh and legacy diagnostics share that
  // typed structure. Keep its declared frame consistent; this does not make
  // Parker the selected background or copy any of its plasma values.
  options.parker.coordinateFrame=options.coordinateFrame;
  options.shock=SEP3D::RuntimeModel::ShockAuthority::None;
  options.intent=SEP3D::RuntimeModel::RunIntent::ShockInjection;
  options.source.enabled=true;
  options.source.samplesPerStep =
      static_cast<std::uint64_t>(integerInputs[0]);
  options.source.minimumEnergyJ=parsed.minimumInjectionEnergyJ;
  options.source.maximumEnergyJ=parsed.maximumInjectionEnergyJ;
  options.source.maximumMacroparticlesPerSpeciesPerStep=
      parsed.maximumInjectionEventsPerSpeciesPerStep;
  if(parsed.momentumPowerLawModel=="compression-ratio") {
    options.source.spectrumModel=
        SEP3D::RuntimeModel::SourceSpectrumModel::LocalCompressionDsa;
    options.source.fixedPhaseSpacePowerIndex=0.0;
  } else {
    options.source.spectrumModel=
        SEP3D::RuntimeModel::SourceSpectrumModel::FixedPhaseSpacePowerLaw;
    options.source.fixedPhaseSpacePowerIndex=
        parsed.fixedPhaseSpacePowerIndex;
  }
  options.source.weightingModel=
      parsed.particleWeightingModel=="constant-statistical-weight"
          ? SEP3D::RuntimeModel::SourceWeightingModel::
                ConstantStatisticalWeight
          : SEP3D::RuntimeModel::SourceWeightingModel::
                LogUniformMomentumImportance;
  // AcceptedShockIncidentFlux already is the declared physical seed rate.
  // No hidden efficiency is applied to the live per-face sum.
  options.source.injectionEfficiency=1.0;
  options.particleNumerics.deriveFromMeshAndShock=true;
  options.particleNumerics.sourceRateModel=
      SEP3D::RuntimeModel::SourceRateNormalizationModel::
          AcceptedShockIncidentFlux;
  options.particleNumerics.maximumParticleSpeedMPerS=numericalInputs[0];
  options.particleNumerics.timeStepMarginFactor=numericalInputs[1];
  options.particleNumerics.sourceNormalizationRadiusM=numericalInputs[2];
  std::shared_ptr<const SEP3D::RuntimeModel::RunConfiguration3D> resolved;
  SEP3D::Core::Status status =
      SEP3D::RuntimeModel::RunConfiguration3D::Create(options, &resolved);
  if (status.ok())
    status = SEP3D::ApplicationRuntime().ReplaceConfigurationBeforeMesh(
        resolved);
  if (!status.ok()) StopWithStatus("application input commit", status);
  gParsedApplicationInput=parsed;
  gHasParsedApplicationInput=true;

  if (PIC::ThisThread == 0) {
    std::cout << SEP3D::RuntimeModel::Sep3dApplicationInputSummary(parsed)
              << "  pre_mesh_configuration_fingerprint="
              << resolved->physics_fingerprint() << '\n';
    if (!resolved->options().parallelDiffusionModelId.empty())
      std::cout << "  parallel_diffusion_model="
                << resolved->options().parallelDiffusionModelId << '\n'
                << "  parallel_diffusion_configuration_fingerprint="
                << resolved->options()
                       .parallelDiffusionConfigurationFingerprint << '\n';
  }
}

std::uint64_t Fnv1a64(const std::string& text) {
  std::uint64_t value = UINT64_C(1469598103934665603);
  for (unsigned char c : text) {
    value ^= c;
    value *= UINT64_C(1099511628211);
  }
  return value;
}

std::shared_ptr<const SEP3D::Background::BackgroundSnapshot>
    gInstalledBackground;
std::shared_ptr<const SEP3D::Background::BackgroundSnapshot>
    gStagedBackground;
// Model ownership is generic; no SWCME-specific switch belongs in cell writes.
// The live provider may already have prepared a candidate while movers still
// resolve the previous immutable gInstalledBackground. Update boundaries are
// joined, so no particle phase overlaps preparation/publication.
std::shared_ptr<SEP3D::Background::BackgroundProvider>
    gRuntimeBackgroundProvider;
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
// Counts completed refreshes only: initial fill is generation setup, not an
// update. Increment after halo completion AND Runtime commit for native evidence.
std::uint64_t gBackgroundPublishedUpdates = 0;
// The native propagation driver may ask for a cadence snapshot and a nearby
// physics landmark on the same tick.  One process-wide tick guard prevents a
// second collective AMPS write while preserving restart-independent naming.
std::uint64_t gLastReducedOutputTick = UINT64_MAX;
std::unordered_map<PIC::Mesh::cDataCenterNode*, std::size_t> gCellSampleIndex;
// The registered AMPS object is owned by Sphere::InternalSpheres for the
// process lifetime.  This non-owning pointer is retained only to prove that
// registration occurred exactly once before mesh construction.
cInternalSphericalData* gSolarSurfaceBoundary = nullptr;

// Static leaf-mask expectations captured before AMPS allocates blocks.  The
// mask is the union of the optional Parker transport corridor and leaves that
// lie wholly inside the solid solar photosphere.  It is checked immediately
// after allocation so an API or ordering change cannot silently turn a
// correct replicated flag plan into a partially resident mesh.
bool gStaticLeafMaskInstalled = false;
// Installation describes a verified plan, including a full-domain identity
// plan.  Pruning and subsequent block-allocation verification are independent
// facts; zero removed leaves must not make a valid plan look uninstalled.
bool gStaticLeafMaskPruningApplied = false;
bool gActiveRegionAllocationVerified = false;
std::size_t gPlannedActiveLeafCount = 0;
std::size_t gPlannedInactiveLeafCount = 0;
std::size_t gPlannedSolarInteriorLeafCount = 0;

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

void InitializeSelectedModelsAfterParser() {
  if(!gHasParsedApplicationInput)return;
  if(gRuntimeBackgroundProvider)
    StopWithStatus("parsed model initialization",SEP3D::Core::Status(
        SEP3D::Core::StatusCode::InvalidTransition,
        "runtime background provider was initialized more than once"));
  SEP3D::Core::Status status=
      SEP3D::RuntimeModel::CreateBackgroundProvider(
          Configuration(),&gRuntimeBackgroundProvider);
  int local=status.ok()?1:0,global=0;
  MPI_Allreduce(&local,&global,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  if(!global)StopWithStatus("parsed model initialization",status.ok()?
      SEP3D::Core::Status(SEP3D::Core::StatusCode::BackgroundInvalid,
          "another MPI rank rejected the parsed reduced model"):status);
  const auto reduced=std::dynamic_pointer_cast<
      SEP3D::Adapters::ShockFrontBackgroundAdapter>(
          gRuntimeBackgroundProvider);
  if(!reduced||!reduced->SharedProvider())
    StopWithStatus("parsed model initialization",SEP3D::Core::Status(
        SEP3D::Core::StatusCode::ConfigurationConflict,
        "selected shock/background models did not construct the reduced provider"));
  if(PIC::ThisThread==0)
    std::cout << "[srcSEP3D] parsed reduced shock/background initialized"
              << " event_fingerprint="
              << reduced->SharedProvider()->Event().physicsFingerprint
              << " handoff_radius_m="
              << reduced->SharedProvider()->Event().handoffApexRadiusM
              << '\n';
}

void ValidateParsedRuntimeMeshBackgroundABI() {
  // Call this function unconditionally on every rank after parsing.  Keeping
  // the compile-time DATAFILE check behind one all-rank function boundary is
  // more than source organization: solar-boundary registration follows it and
  // must never become rank-local because every rank constructs the same AMR
  // cut surface.  Unsupported coupler builds still fail before any mesh or
  // internal-boundary state is allocated.
#if _PIC_COUPLER_MODE_ != _PIC_COUPLER_MODE__DATAFILE_
  if (Configuration().options().background ==
      SEP3D::RuntimeModel::BackgroundAuthority::RuntimeModel)
    StopWithStatus("parsed runtime mesh background", SEP3D::Core::Status(
        SEP3D::Core::StatusCode::ConfigurationConflict,
        "the reduced runtime provider requires the one-fluid AMPS DATAFILE buffer layout"));
#endif
}

double ConfiguredParticleWeight(int ampsIndex) {
  const auto& options=Configuration().options();
  if(!options.speciesParticleNormalizations.empty()) {
    for(const auto& item:options.speciesParticleNormalizations)
      if(item.ampsIndex==ampsIndex)return item.macroparticleWeight;
    StopWithStatus("particle weight lookup",SEP3D::Core::Status(
        SEP3D::Core::StatusCode::ConfigurationConflict,
        "no derived particle normalization exists for compiled species "+
        std::to_string(ampsIndex)));
  }
  return options.species.macroparticleWeight;
}

void FinalizeParticleNumericsAfterMeshAllocation() {
  if(!gHasParsedApplicationInput)return;

  // AMPS' characteristic cell size is the same length used by mature
  // application LocalTimeStep callbacks. Reduce only allocated owner blocks;
  // ghost blocks do not create an independent stability restriction. A rank
  // with no blocks contributes +infinity and the collective minimum remains
  // the smallest physical scale represented anywhere in the domain.
  double localMinimum=std::numeric_limits<double>::infinity();
  for(unsigned int blockIndex=0;
      blockIndex<PIC::DomainBlockDecomposition::nLocalBlocks;++blockIndex) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=
        PIC::DomainBlockDecomposition::BlockTable[blockIndex];
    if(node==nullptr||node->block==nullptr||!node->IsUsedInCalculationFlag)
      continue;
    localMinimum=std::min(localMinimum,node->GetCharacteristicCellSize());
  }
  double globalMinimum=0;
  MPI_Allreduce(&localMinimum,&globalMinimum,1,MPI_DOUBLE,MPI_MIN,
      MPI_GLOBAL_COMMUNICATOR);
  double timeStep=0;
  SEP3D::Core::Status status=
      SEP3D::RuntimeModel::CalculateMeshGlobalTimeStep(globalMinimum,
          gParsedApplicationInput.maximumParticleSpeedMPerS,
          gParsedApplicationInput.timeStepMarginFactor,&timeStep);
  if(!status.ok())StopWithStatus("global particle time step",status);

  const auto reduced=std::dynamic_pointer_cast<
      SEP3D::Adapters::ShockFrontBackgroundAdapter>(
          gRuntimeBackgroundProvider);
  if(!reduced||!reduced->SharedProvider())
    StopWithStatus("particle source normalization",SEP3D::Core::Status(
        SEP3D::Core::StatusCode::SnapshotUnavailable,
        "reduced provider is unavailable after mesh allocation"));
  const auto flux=
      SEP::CoronaSwcme::ShockFront::EvaluateIncidentParticleFluxAtApexRadius(
          *reduced->SharedProvider(),
          gParsedApplicationInput.sourceNormalizationRadiusM);
  if(!flux.ok())StopWithStatus("particle source normalization",
      SEP3D::Core::Status(SEP3D::Core::StatusCode::InvalidInput,
          flux.status.message));

  std::vector<SEP3D::RuntimeModel::SpeciesParticleNormalization>
      normalizations;
  status=SEP3D::RuntimeModel::CalculateSpeciesParticleNormalizations(
      gCompiledSpecies,flux.value,timeStep,
      gParsedApplicationInput.particlesPerIteration,&normalizations);
  if(!status.ok())StopWithStatus("per-species particle weight",status);

  // Deterministic shared physics should produce bit-identical values, but
  // explicitly checking the cross-rank range turns an asset/filesystem or
  // floating-environment disagreement into a pre-transport failure.
  std::vector<double> rankValues={globalMinimum,timeStep,flux.value.epochS,
      flux.value.apexRadiusM,flux.value.acceptedAreaM2,
      flux.value.excludedPhysicalAreaM2,flux.value.protonRatePerS,
      flux.value.electronRatePerS,flux.value.alphaRatePerS};
  for(const auto& item:normalizations) {
    rankValues.push_back(item.physicalSourceRatePerS);
    rankValues.push_back(item.macroparticleWeight);
  }
  std::vector<double> minimum(rankValues.size()),maximum(rankValues.size());
  MPI_Allreduce(rankValues.data(),minimum.data(),
      static_cast<int>(rankValues.size()),MPI_DOUBLE,MPI_MIN,
      MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(rankValues.data(),maximum.data(),
      static_cast<int>(rankValues.size()),MPI_DOUBLE,MPI_MAX,
      MPI_GLOBAL_COMMUNICATOR);
  for(std::size_t index=0;index<rankValues.size();++index)
    if(minimum[index]!=maximum[index])
      StopWithStatus("particle numerics rank agreement",SEP3D::Core::Status(
          SEP3D::Core::StatusCode::ConfigurationConflict,
          "MPI ranks derived different time-step/source/weight values"));

  // This is the second and final immutable replacement. AMPS has allocated
  // blocks, but Runtime has not yet bound the mesh or acquired a background.
  // Storage layout is unchanged; only deterministic numerics derived from the
  // now-observable mesh and already-resolved provider are committed.
  SEP3D::RuntimeModel::RunConfiguration3DOptions options=
      Configuration().options();
  std::vector<double> requestedObserverCadences;
  requestedObserverCadences.reserve(options.observers.size());
  for(const auto& observer:options.observers)
    requestedObserverCadences.push_back(observer.cadenceS);
  options.requestedTimeStepS=timeStep;
  options.particleNumerics.resolvedMinimumCellSizeM=globalMinimum;
  options.speciesParticleNormalizations=normalizations;
  // Retain the first species in the legacy scalar for old diagnostics; every
  // AMPS assignment below uses the explicit per-species table.
  options.species.macroparticleWeight=normalizations.front().macroparticleWeight;
  std::shared_ptr<const SEP3D::RuntimeModel::RunConfiguration3D> resolved;
  status=SEP3D::RuntimeModel::AlignObserverCadencesToGlobalStep(
      timeStep,&options.observers);
  if(status.ok())
    status=SEP3D::RuntimeModel::RunConfiguration3D::Create(options,&resolved);
  if(status.ok())status=SEP3D::RuntimeModel::ValidateCompiledSpeciesBinding(
      resolved->options(),PIC::nTotalSpecies,gCompiledSpecies);
  if(status.ok())status=SEP3D::ApplicationRuntime().
      ReplaceConfigurationBeforeMesh(resolved);
  if(!status.ok())StopWithStatus("particle numerics commit",status);

  if(PIC::ThisThread==0) {
    std::cout << std::scientific << std::setprecision(17)
              << "[srcSEP3D] global particle numerics summary\n"
              << "  minimum_allocated_cell_size_m=" << globalMinimum << '\n'
              << "  maximum_particle_speed_m_s="
              << gParsedApplicationInput.maximumParticleSpeedMPerS << '\n'
              << "  time_step_margin_factor="
              << gParsedApplicationInput.timeStepMarginFactor << '\n'
              << "  global_time_step_s=" << timeStep << '\n'
              << "  source_model=" << gParsedApplicationInput.sourceModel
              << '\n' << "  source_normalization_radius_m="
              << flux.value.apexRadiusM << '\n'
              << "  source_normalization_epoch_s=" << flux.value.epochS
              << '\n' << "  accepted_shock_area_m2="
              << flux.value.acceptedAreaM2 << '\n'
              << "  excluded_physical_area_m2="
              << flux.value.excludedPhysicalAreaM2 << '\n';
    for(const auto& item:normalizations)
      std::cout << "  species[" << item.ampsIndex << "].symbol="
                << item.symbol << '\n' << "  species[" << item.ampsIndex
                << "].source_rate_s-1=" << item.physicalSourceRatePerS
                << '\n' << "  species[" << item.ampsIndex
                << "].particle_weight=" << item.macroparticleWeight << '\n';
    for(std::size_t index=0;index<options.observers.size();++index)
      std::cout << "  observer[" << index << "].id="
                << options.observers[index].id << '\n'
                << "  observer[" << index << "].requested_cadence_s="
                << requestedObserverCadences[index] << '\n'
                << "  observer[" << index << "].resolved_cadence_s="
                << options.observers[index].cadenceS << '\n';
    std::cout << "  final_configuration_fingerprint="
              << resolved->physics_fingerprint() << std::defaultfloat << '\n';
    if (!resolved->options().parallelDiffusionModelId.empty())
      std::cout << "  parallel_diffusion_model="
                << resolved->options().parallelDiffusionModelId << '\n'
                << "  parallel_diffusion_configuration_fingerprint="
                << resolved->options()
                       .parallelDiffusionConfigurationFingerprint << '\n';
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
  result.solarRefinementAnchor = options.solarRefinementAnchor;
  result.activeSolarSphereRadiusM = options.activeSolarSphereRadiusM;
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
  // The active-domain corridor deliberately has its own width law.  A user
  // may refine a narrow core while retaining a wider particle-transport halo;
  // silently reusing the refinement width here could deactivate blocks that
  // are needed by an observer, shock source, or finite-difference stencil.
  result.activeRegion = options.activeRegion;
  result.activeTubeReferenceRadiusM = options.activeTubeReferenceRadiusM;
  result.activeTubeRadiusAtReferenceM =
      options.activeTubeRadiusAtReferenceM;
  result.activeTubeRadiusMode = options.activeTubeRadiusMode;
  result.activeTubeBufferBlocks = options.activeTubeBufferBlocks;
  result.solarWindSpeedMPerS = options.parker.solarWindSpeedMPerS;
  result.solarRotationRateRadPerS = options.parker.solarRotationRateRadPerS;
  result.rotationAxis = options.parker.rotationAxis;
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
          // Background coverage begins at the configurable Parker/CME source
          // shell (normally 20 R_sun), outside the registered 1-R_sun solid
          // photosphere.  The enclosing Cartesian cube also contains padding
          // inside that source shell and beyond the outer heliocentric sphere;
          // those cells stay zero-initialized and are explicitly marked
          // background_valid=0 in Tecplot output.  Measure radius from the
          // configured heliocentric origin rather than assuming (0,0,0).
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

void ApplyPopulationControlAtBoundary() {
  using namespace SEP3D;
  const RuntimeModel::RunConfiguration3DOptions& options =
      Configuration().options();
  if (options.populationControl !=
      RuntimeModel::PopulationControlMode::SplitMerge) {
    return;
  }
  const std::uint64_t step = ApplicationRuntime().counters().currentTick;
  if (step % options.populationControlCadenceSteps != 0) return;

  AMPS::Movers::PopulationControlRequest request;
  request.step = step;
  request.campaignSeed = options.campaignSeed;
  request.minimumParticlesPerCellPerSpecies =
      options.minimumParticlesPerCellPerSpecies;
  request.targetParticlesPerCellPerSpecies =
      options.targetParticlesPerCellPerSpecies;
  request.maximumParticlesPerCellPerSpecies =
      options.maximumParticlesPerCellPerSpecies;
  const AMPS::Movers::PopulationControlReport local =
      AMPS::Movers::ApplyPopulationControl(request);
  if (!local.status.ok()) {
    std::cerr << "[srcSEP3D] rank " << PIC::ThisThread
              << " population control failed: "
              << local.status.message << '\n';
  }
  int localOK = local.status.ok() ? 1 : 0;
  int globalOK = 0;
  MPI_Allreduce(&localOK, &globalOK, 1, MPI_INT, MPI_MIN,
                MPI_GLOBAL_COMMUNICATOR);
  if (globalOK == 0) {
    StopWithStatus("population control", Core::Status(
        Core::StatusCode::Error,
        "one or more MPI ranks rejected SEP-aware split/merge"));
  }

  const unsigned long long localCounts[5] = {
      static_cast<unsigned long long>(local.occupiedCellSpecies),
      static_cast<unsigned long long>(local.splitOperations),
      static_cast<unsigned long long>(local.mergeOperations),
      static_cast<unsigned long long>(local.particlesBefore),
      static_cast<unsigned long long>(local.particlesAfter)};
  unsigned long long globalCounts[5] = {};
  MPI_Allreduce(localCounts, globalCounts, 5, MPI_UNSIGNED_LONG_LONG,
                MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
  const double localResiduals[3] = {
      local.maximumRelativeWeightResidual,
      local.maximumRelativeMomentumResidual,
      local.maximumRelativeEnergyResidual};
  double globalResiduals[3] = {};
  MPI_Allreduce(localResiduals, globalResiduals, 3, MPI_DOUBLE, MPI_MAX,
                MPI_GLOBAL_COMMUNICATOR);
  if (PIC::ThisThread == 0) {
    std::cout << "[srcSEP3D] population control step=" << step
              << " occupied_cell_species=" << globalCounts[0]
              << " split_operations=" << globalCounts[1]
              << " merge_operations=" << globalCounts[2]
              << " particles_before=" << globalCounts[3]
              << " particles_after=" << globalCounts[4]
              << " max_weight_residual=" << globalResiduals[0]
              << " max_momentum_residual=" << globalResiduals[1]
              << " max_relativistic_energy_residual="
              << globalResiduals[2] << '\n';
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

// AMPS writes whole cells as FEBRICK zones and boundary cut cells as
// tetrahedral zones. Its writers obtain vertex values by creating a temporary
// cDataCenterNode and calling
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
      bytes % sizeof(double) != 0 || destinationNode == nullptr) {
    StopWithStatus("srcSEP3D center-node interpolation",
        SEP3D::Core::Status(SEP3D::Core::StatusCode::LayoutMismatch,
            "AMPS requested interpolation before the frozen srcSEP3D "
            "static center-node layout was available"));
  }
  if (interpolationCount < 0 ||
      (interpolationCount > 0 &&
       (interpolationList == nullptr || interpolationCoefficients == nullptr))) {
    StopWithStatus("srcSEP3D center-node interpolation",
        SEP3D::Core::Status(SEP3D::Core::StatusCode::InvalidInput,
            "AMPS supplied an invalid center-node interpolation stencil"));
  }

  const std::size_t valueCount = bytes / sizeof(double);
  // Internal-sphere/corridor output vertices can have no positive-volume
  // donors. AMPS intentionally permits count==0 with internal boundaries.
  // The shared interpolator then writes a finite zero placeholder into the
  // entire application slice, including optional gradients and turbulence.
  // Do not skip PrintData: every MPI owner must still send the expected row.
  // Positive density and temperature distinguish populated background below;
  // an empty stencil must never be presented as a physical zero-density wind.
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

// Keep reduced-provider discovery and event-composition physics outside the
// generic AMPS output callback.  The callback is also compiled by a portable
// byte-slice/channel probe; embedding provider RTTI there would make a generic
// interpolation safety test depend on the complete native runtime.  These
// helpers form the narrow native seam, while the callback still owns row
// sizing, ordering, and owner/root transport.
bool ReducedProductionColumnsEnabled() {
  const auto reduced=std::dynamic_pointer_cast<
      SEP3D::Adapters::ShockFrontBackgroundAdapter>(gRuntimeBackgroundProvider);
  return reduced&&Configuration().options().intent==
      SEP3D::RuntimeModel::RunIntent::ShockPropagation;
}

void AppendReducedProductionColumns(
    PIC::Mesh::cDataCenterNode* centerNode,
    std::vector<double>* values,std::size_t* cursor) {
  const auto reduced=std::dynamic_pointer_cast<
      SEP3D::Adapters::ShockFrontBackgroundAdapter>(gRuntimeBackgroundProvider);
  if(!centerNode||!values||!cursor||!reduced||
      Configuration().options().intent!=
          SEP3D::RuntimeModel::RunIntent::ShockPropagation||
      *cursor+4>values->size()) {
    StopWithStatus("reduced production output columns",
        SEP3D::Core::Status(SEP3D::Core::StatusCode::LayoutMismatch,
            "reduced output metadata is absent or its row layout is incomplete"));
  }
  double electronDensityM3=0.0;
  LoadBytes(centerNode,Configuration().storage_layout().numberDensityOffset,
      &electronDensityM3,sizeof(electronDensityM3));
  const auto& composition=
      reduced->SharedProvider()->Event().ambient.composition;
  const double alpha=composition.alphaToProtonNumberRatio;
  // Charge neutrality gives ne=np+2*nalpha.  Reconstruct the event's declared
  // mass density from the actually installed electron density and frozen
  // composition; do not substitute the common rho=mp*ne approximation, which
  // is wrong when alpha particles or electron mass are enabled.
  const double protonDensityM3=electronDensityM3/(1.0+2.0*alpha);
  const double massDensityKgM3=protonDensityM3*(
      SEP::CoronalCME::Constants::kProtonMassKg+
      alpha*SEP::CoronalCME::Constants::kAlphaMassKg)+
      (composition.includeElectronMass?electronDensityM3*
       SEP::CoronalCME::Constants::kElectronMassKg:0.0);
  const auto front=reduced->FrontEpoch();
  const auto* active=SEP3D::ApplicationRuntime().active_snapshot();
  (*values)[(*cursor)++]=massDensityKgM3;
  (*values)[(*cursor)++]=SEP3D::ApplicationRuntime().CurrentTimeS();
  (*values)[(*cursor)++]=active?static_cast<double>(active->generation):0.0;
  (*values)[(*cursor)++]=front?static_cast<double>(front->generation):0.0;
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
  if(ReducedProductionColumnsEnabled()) {
    // These values bind every plotted ambient vertex to the committed native
    // epoch.  Mass density is reconstructed from the event's frozen species
    // composition and the actually stored electron density; it is not an
    // independent model or an assumed downstream CME density.
    std::fprintf(output,
        ", \"mass_density_kg_m-3\""
        ", \"simulation_time_s\""
        ", \"background_generation\""
        ", \"front_generation\"");
  }
  // These flags make the two independent empty-data cases machine-readable:
  // background_valid=0 identifies padding or vertices without physical
  // background donors, while particle_sample_present=0 identifies a cell/species with no
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
  const bool reducedPropagation=ReducedProductionColumnsEnabled();
  std::size_t valueCount = 26+(reducedPropagation?4:0);
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

    if(reducedPropagation)
      AppendReducedProductionColumns(centerNode,&values,&cursor);

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
    // Validated background snapshots require positive density and temperature.
    // Their interpolated values therefore identify initialized donors without
    // changing the frozen storage/restart ABI. A no-donor stencil is all zero
    // even when its vertex lies inside the radial shell (e.g. a corridor edge).
    double numberDensityM3 = 0.0, temperatureK = 0.0;
    LoadBytes(centerNode, layout.numberDensityOffset, &numberDensityM3,
              sizeof(numberDensityM3));
    LoadBytes(centerNode, layout.temperatureOffset, &temperatureK,
              sizeof(temperatureK));
    const bool backgroundStateAvailable =
        numberDensityM3 > 0.0 && temperatureK > 0.0;
    const SEP3D::Output::TecplotCellPresentation presentation =
        SEP3D::Output::PrepareTecplotCellPresentation(
            storedBackground, insidePhysicalShell, PIC::LastSampleLength,
            sampledParticleNumber, backgroundStateAvailable);
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

SEP3D::Core::Status ResolveStoredBackground(
    const SEP3D::Core::Vec3& positionM,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node,
    SEP3D::Background::BackgroundSample* resolved,
    PIC::Mesh::cDataCenterNode** resolvedCell) {
  using namespace SEP3D;
  if (node == nullptr || node->block == nullptr || resolved == nullptr)
    return Core::Status(Core::StatusCode::NotFound,
                        "background lookup has no allocated AMR block");
  double position[3]; positionM.CopyTo(position);
  int i = 0, j = 0, k = 0;
  const int localCell = PIC::Mesh::mesh->FindCellIndex(
      position, i, j, k, node, false);
  if (localCell < 0)
    return Core::Status(Core::StatusCode::NotFound,
                        "position has no AMR center cell");
  PIC::Mesh::cDataCenterNode* cell = node->block->GetCenterNode(localCell);
  if (cell == nullptr)
    return Core::Status(Core::StatusCode::BackgroundInvalid,
                        "position has an uninitialized AMR center cell");

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
    // Remote/ghost resolution needs the same strain/compression tensor as
    // owner snapshots when its optional application storage was requested.
    if (layout.velocityGradientOffset != RuntimeModel::kNoOffset)
      LoadBytes(cell, layout.velocityGradientOffset, background.gradU.m,
                sizeof(background.gradU.m));
    if (gInstalledBackground) {
      background.generation=gInstalledBackground->metadata().generation;
      // Ghost data carry the same generation after the joined halo exchange.
      background.configurationDigest=gInstalledBackground->samples().empty()?0:
          gInstalledBackground->samples().front().configurationDigest;
    }
    background.valid =
        std::isfinite(background.absB) && background.absB > 0.0;
    background.status = background.valid
        ? Core::Status::OK()
        : Core::Status(Core::StatusCode::BackgroundInvalid,
                       "stored AMR magnetic field is invalid");
  }
  if (!background.status.ok()) return background.status;
  *resolved = background;
  if (resolvedCell != nullptr) *resolvedCell = cell;
  return Core::Status::OK();
}

SEP3D::Core::Status ResolvePopulationMagneticDirection(
    const SEP3D::Core::Vec3& positionM,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node,
    SEP3D::Core::Vec3* bHat) {
  if (bHat == nullptr)
    return SEP3D::Core::Status(SEP3D::Core::StatusCode::InvalidInput,
                              "magnetic-direction output is null");
  SEP3D::Background::BackgroundSample background;
  const SEP3D::Core::Status status = ResolveStoredBackground(
      positionM, node, &background, nullptr);
  if (!status.ok()) return status;
  if (!std::isfinite(background.bHat.x) ||
      !std::isfinite(background.bHat.y) ||
      !std::isfinite(background.bHat.z) ||
      std::fabs(background.bHat.Norm() - 1.0) > 1.0e-12) {
    return SEP3D::Core::Status(
        SEP3D::Core::StatusCode::BackgroundInvalid,
        "stored background has no unit magnetic direction");
  }
  *bHat = background.bHat;
  return SEP3D::Core::Status::OK();
}

SEP3D::Adapters::ExpandingShock ReducedMoverShock(
    const SEP::CoronaSwcme::ShockFront::Epoch& epoch,
    const SEP::CoronaSwcme::ShockFront::Configuration& event) {
  SEP3D::Adapters::ExpandingShock shock;
  shock.centerM=SEP3D::Core::Vec3();
  shock.radiusAtStepStartM=epoch.trajectory.apexRadiusM;
  shock.radialSpeedMPerS=epoch.trajectory.apexSpeedMPerS;
  shock.generation=epoch.generation;
  shock.active=epoch.apexShockAccepted;
  shock.geometry=SEP3D::Adapters::ShockGeometryKind::FiniteSSE;
  shock.cmeDirection=SEP3D::Core::Vec3(
      event.direction.x,event.direction.y,event.direction.z);
  shock.halfWidthRad=event.halfWidthRad;
  return shock;
}

double RelativisticKineticEnergyJ(double momentumKgMPerS,double massKg) {
  const long double p=momentumKgMPerS;
  const long double m=massKg;
  const long double c=SEP3D::Core::Const::c;
  return static_cast<double>(
      (std::sqrt(p*p*c*c+m*m*c*c*c*c)-m*c*c));
}

// AMPS calls this function once per rank at the injection phase of every
// PIC::TimeStep. All ranks reconstruct the same Poisson candidates from
// (campaign,species,tick,generation) keys. The replicated AMR tree then
// selects exactly one owner for each point, so MPI decomposition changes who
// allocates a particle but not the stochastic physical source.
long int InjectReducedShockSurfaceParticles() {
  using namespace SEP3D;
  const auto reduced=std::dynamic_pointer_cast<
      Adapters::ShockFrontBackgroundAdapter>(gRuntimeBackgroundProvider);
  const auto provider=reduced?reduced->SharedProvider():nullptr;
  const auto epoch=reduced?reduced->FrontEpoch():nullptr;
  if(!provider||!epoch)
    StopWithStatus("reduced-front particle injection",Core::Status(
        Core::StatusCode::SnapshotUnavailable,
        "the committed reduced-front epoch is unavailable"));
  const auto& options=Configuration().options();
  const double dt=options.requestedTimeStepS;
  const std::uint64_t particleStep=
      ApplicationRuntime().counters().completedSteps;
  unsigned long long localTotal=0;

  for(const auto& species:gCompiledSpecies) {
    RuntimeModel::SurfaceInjectionBatch batch;
    Core::Status status;
    if(options.source.weightingModel==
        RuntimeModel::SourceWeightingModel::ConstantStatisticalWeight) {
      status=RuntimeModel::GenerateConstantWeightSurfaceInjectionBatch(
          *provider,*epoch,species,options.source,
          ConfiguredParticleWeight(species.ampsIndex),dt,
          options.campaignSeed,particleStep,&batch);
    } else {
      status=RuntimeModel::GenerateLogUniformMomentumImportanceBatch(
          *provider,*epoch,species,options.source,
          ConfiguredParticleWeight(species.ampsIndex),dt,
          options.campaignSeed,particleStep,&batch);
    }
    if(!status.ok())StopWithStatus("reduced-front source sampling",status);

    Adapters::InjectionPlan localPlan;
    localPlan.status=Core::Status::OK();
    std::uint64_t disconnected=0;
    double localEnergyJ=0;
    Core::Vec3 localMomentum;
    for(const RuntimeModel::SurfaceInjectionEvent& event:batch.events) {
      const double radius=event.positionM.Norm();
      if(!std::isfinite(radius)||radius<=0)
        StopWithStatus("reduced-front source direction",Core::Status(
            Core::StatusCode::InvalidInput,
            "a sampled shock point has no heliocentric direction"));
      double x[3];event.positionM.CopyTo(x);
      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=
          PIC::Mesh::mesh->findTreeNode(x);
      // A Cartesian AMR block can geometrically cover the solar interior or
      // the deliberately excluded r<innerRadius transport region. Native
      // background storage is undefined there even when the tree node exists.
      // Connectivity therefore requires both radial physics support and an
      // active AMR leaf. These predicates are replicated, so every rank adds
      // the same event once to the global disconnected budget; the event is
      // never moved to another face or used to renormalize the physical rate.
      if(radius<options.innerRadiusM||radius>options.outerRadiusM||
          node==nullptr||!node->IsUsedInCalculationFlag) {
        ++disconnected;
        continue;
      }
      if(node->Thread!=PIC::ThisThread)continue;
      if(node->block==nullptr)
        StopWithStatus("reduced-front source ownership",Core::Status(
            Core::StatusCode::LayoutMismatch,
            "the owning active AMR node has no allocated block"));

      const Core::Vec3 antiSunward=event.positionM/radius;
      Core::Vec3 bHat;
      status=ResolvePopulationMagneticDirection(event.positionM,node,&bHat);
      if(!status.ok())StopWithStatus("reduced-front source magnetic basis",status);
      double mu=0,gyrophase=0;
      status=AMPS::Movers::GyrotropicCoordinatesForDirection(
          antiSunward,bHat,&mu,&gyrophase);
      if(!status.ok())StopWithStatus("reduced-front source direction",status);

      Adapters::InjectedParticle injected;
      injected.status=Core::Status::OK();
      injected.remainingFirstStepFraction=event.remainingStepFraction;
      injected.particle.stableId=event.stableId;
      injected.particle.species=species.ampsIndex;
      injected.particle.positionM=event.positionM;
      injected.particle.momentumKgMPerS=event.momentumKgMPerS;
      injected.particle.mu=mu;
      injected.particle.gyrophaseRad=gyrophase;
      injected.particle.statisticalWeight=
          ConfiguredParticleWeight(species.ampsIndex);
      injected.particle.completedStep=particleStep;
      // The particle is born on this generation. Marking it prevents the
      // geometric intersection finder from treating the t=0 birth point as
      // an additional crossing of the same shock.
      injected.particle.lastShockGeneration=epoch->generation;
      localPlan.particles.push_back(std::move(injected));
    }

    const AMPS::Movers::InjectionOutcome outcome=
        AMPS::Movers::InjectParticles(localPlan);
    if(!outcome.status.ok())
      StopWithStatus("reduced-front AMPS particle allocation",outcome.status);
    status=gParticleLedger.RecordInjection(
        particleStep,species.ampsIndex,outcome.allocated);
    if(!status.ok())StopWithStatus("reduced-front particle ledger",status);
    if(outcome.allocated>static_cast<std::uint64_t>(LONG_MAX)||
        PIC::BC::nInjectedParticles[species.ampsIndex]>
            LONG_MAX-static_cast<long int>(outcome.allocated))
      StopWithStatus("reduced-front injection counter",Core::Status(
          Core::StatusCode::Error,"native injection counter overflow"));
    PIC::BC::nInjectedParticles[species.ampsIndex]+=
        static_cast<long int>(outcome.allocated);
    PIC::BC::ParticleProductionRate[species.ampsIndex]+=
        outcome.allocated*ConfiguredParticleWeight(species.ampsIndex)/dt;
    PIC::BC::ParticleMassProductionRate[species.ampsIndex]+=
        outcome.allocated*ConfiguredParticleWeight(species.ampsIndex)*
        species.massKg/dt;
    localTotal+=outcome.allocated;

    for(const auto& injected:localPlan.particles) {
      const double represented=ConfiguredParticleWeight(species.ampsIndex);
      localEnergyJ+=represented*RelativisticKineticEnergyJ(
          injected.particle.momentumKgMPerS,species.massKg);
      localMomentum+=represented*injected.particle.momentumKgMPerS*
          injected.particle.positionM.Normalized();
    }
    unsigned long long localCount=outcome.allocated,globalCount=0;
    MPI_Allreduce(&localCount,&globalCount,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,
        MPI_GLOBAL_COMMUNICATOR);
    double localConserved[4]={localEnergyJ,localMomentum.x,localMomentum.y,
        localMomentum.z},globalConserved[4]={};
    MPI_Allreduce(localConserved,globalConserved,4,MPI_DOUBLE,MPI_SUM,
        MPI_GLOBAL_COMMUNICATOR);
    if(globalCount+disconnected!=batch.events.size())
      StopWithStatus("reduced-front source MPI ownership",Core::Status(
          Core::StatusCode::LayoutMismatch,
          "Poisson candidates were not assigned exactly once to an owner or "
          "the explicit disconnected budget"));

    if(PIC::ThisThread==0) {
      Adapters::SourceLedgerRow row;
      row.step=particleStep;
      row.species=species.ampsIndex;
      row.shockGeneration=epoch->generation;
      row.sourceId=Fnv1a64("reduced-front:"+
          std::to_string(epoch->generation)+":"+
          std::to_string(particleStep)+":"+
          std::to_string(species.ampsIndex));
      if(row.sourceId==0)row.sourceId=1;
      row.representedParticles=globalCount*
          ConfiguredParticleWeight(species.ampsIndex);
      row.injectedEnergyJ=globalConserved[0];
      row.injectedMomentumKgMPerS=Core::Vec3(
          globalConserved[1],globalConserved[2],globalConserved[3]);
      row.macroparticles=globalCount;
      row.rejected=disconnected;
      row.disconnectedPatches=disconnected==0?0:1;
      gSourceLedger.push_back(row);
      std::cout<<std::scientific<<std::setprecision(17)
          <<"[srcSEP3D] reduced-front source step="<<particleStep
          <<" generation="<<epoch->generation
          <<" species="<<species.symbol
          <<" accepted_faces="<<batch.distribution.faces.size()
          <<" total_source_rate_s-1="
          <<batch.distribution.physicalRatePerS
          <<" poisson_candidates="<<batch.events.size()
          <<" injected_all_ranks="<<globalCount
          <<" disconnected="<<disconnected<<'\n';
    }
  }
  if(localTotal>static_cast<unsigned long long>(LONG_MAX))
    StopWithStatus("reduced-front injection return",Core::Status(
        Core::StatusCode::Error,"rank-local injected count exceeds long int"));
  return static_cast<long int>(localTotal);
}

SEP3D::Core::Status ResolveLocalTransportImpl(
    const SEP3D::Core::Vec3& positionM, int species,
    double momentumKgMPerS, double mu,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node,
    SEP3D::Adapters::LocalTransportRecord* local,
    bool evaluateParallelGradient) {
  using namespace SEP3D;
  if (local == nullptr)
    return Core::Status(Core::StatusCode::InvalidInput,
                        "local transport output is null");
  Background::BackgroundSample background;
  PIC::Mesh::cDataCenterNode* cell = nullptr;
  const Core::Status backgroundStatus = ResolveStoredBackground(
      positionM, node, &background, &cell);
  if (!backgroundStatus.ok()) return backgroundStatus;
  local->background = background;
  const RuntimeModel::StorageLayout& layout = Configuration().storage_layout();
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
  if (waves.deltaB2T2 > 0.0) {
    local->plusWaveFraction = waves.deltaBPlus2T2 / waves.deltaB2T2;
    local->minusWaveFraction = waves.deltaBMinus2T2 / waves.deltaB2T2;
  } else {
    // A zero-variance ballistic sample has no physical branch preference.
    // The event rate is also zero/infinite-MFP, so this symmetric partition is
    // an inert but finite sentinel rather than an assumed turbulence model.
    local->plusWaveFraction = 0.5;
    local->minusWaveFraction = 0.5;
  }
  const auto& coefficientOptions = Configuration().options();
  Turbulence::CoefficientSelection coefficientSelection;
  coefficientSelection.spatial = coefficientOptions.spatialDiffusionModel;
  coefficientSelection.pitchAngle =
      coefficientOptions.pitchAngleDiffusionModel;
  coefficientSelection.meanFreePath = coefficientOptions.meanFreePathModel;
  coefficientSelection.constantDmumuPerS =
      coefficientOptions.constantDmumuPerS;
  coefficientSelection.constantMeanFreePathM =
      coefficientOptions.constantMeanFreePathM;
  coefficientSelection.meanFreePathReferenceM =
      coefficientOptions.meanFreePathReferenceM;
  coefficientSelection.meanFreePathReferenceRadiusM =
      coefficientOptions.meanFreePathReferenceRadiusM;
  coefficientSelection.meanFreePathReferenceRigidityV =
      coefficientOptions.meanFreePathReferenceRigidityV;
  coefficientSelection.meanFreePathRadialExponent =
      coefficientOptions.meanFreePathRadialExponent;
  coefficientSelection.meanFreePathRigidityExponent =
      coefficientOptions.meanFreePathRigidityExponent;
  coefficientSelection.quadratureAbsoluteToleranceM2PerS =
      coefficientOptions.spatialQuadratureAbsoluteToleranceM2PerS;
  coefficientSelection.quadratureRelativeTolerance =
      coefficientOptions.spatialQuadratureRelativeTolerance;
  coefficientSelection.quadratureMaximumRecursion =
      coefficientOptions.spatialQuadratureMaximumRecursion;
  coefficientSelection.timeS = ApplicationRuntime().CurrentTimeS();
  coefficientSelection.solarOriginM = coefficientOptions.coordinateOriginM;
  const bool parkerMover = coefficientOptions.transport ==
      RuntimeModel::TransportModel::Parker3D;
  const bool focusedDiffusionMover = coefficientOptions.transport ==
      RuntimeModel::TransportModel::FocusedDiffusion3D;
  const bool focusedScatteringMover = coefficientOptions.transport ==
      RuntimeModel::TransportModel::FocusedScattering3D;
  coefficientSelection.requireSpatialDiffusion = parkerMover ||
      coefficientOptions.perpendicularDiffusion ==
          RuntimeModel::PerpendicularDiffusionMode::ConstantRatio;
  coefficientSelection.requirePitchAngleDiffusion =
      focusedDiffusionMover;
  coefficientSelection.requireMeanFreePath = focusedScatteringMover;
  const Turbulence::LocalScatteringCoefficients coefficients =
      Turbulence::EvaluateLocalScattering(
          waves, background, positionM, species,
          PIC::MolecularData::GetMass(species),
          PIC::MolecularData::GetElectricCharge(species),
          momentumKgMPerS, mu, coefficientSelection);
  if (!coefficients.status.ok()) return coefficients.status;
  local->kappaParallelM2PerS = coefficients.kappaParallelM2PerS;
  local->meanFreePathM = coefficients.meanFreePathM;
  local->dMuMuPerS = coefficients.dMuMuPerS;
  local->dDmuMuDmuPerS = coefficients.dDmuMuDmuPerS;
  local->parallelDiffusionDiagnosticMask =
      coefficients.parallelDiagnosticMask;
  local->parallelDiffusionModelId = coefficients.parallelModelId;
  local->parallelDiffusionConfigurationFingerprint =
      coefficients.parallelConfigurationFingerprint;

  // The Parker Ito drift requires b-hat dot grad(kappa_parallel), not merely
  // kappa itself.  Evaluate the same provider/coefficient chain one local-cell
  // spacing in both field-aligned directions.  The spacing is enlarged by the
  // largest b component so at least one Cartesian coordinate crosses a cell
  // centre spacing even when the field is oblique to the AMR axes.
  if (evaluateParallelGradient && parkerMover) {
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

// Validate offsets without writing or aborting, so every rank can join the
// candidate gate even when its native layout is invalid. One sample represents
// one plasma fluid; a multifluid build needs an explicit future physics mapping.
SEP3D::Core::Status NativeAmpsBackgroundLayoutStatus() {
  using namespace PIC::CPLR::DATAFILE;
  using namespace SEP3D::Core;
  auto invalid=[](const std::string& message) {
    return Status(StatusCode::LayoutMismatch,message);
  };
  if (CenterNodeAssociatedDataOffsetBegin<0 || MULTIFILE::CurrDataFileOffset<0 ||
      nTotalBackgroundVariables<=0 || nIonFluids!=1)
    return invalid("runtime backgrounds require allocated one-fluid DATAFILE storage");
  const cOffsetElement* fields[]={&Offset::PlasmaNumberDensity,&Offset::PlasmaBulkVelocity,
      &Offset::PlasmaTemperature,&Offset::PlasmaIonPressure,&Offset::PlasmaDivU,
      &Offset::MagneticField,&Offset::ElectricField,&Offset::MagneticFieldGradient,
      &Offset::Current,&Offset::PlasmaElectronPressure};
  const int components[]={1,3,1,1,1,3,3,9,3,1};
  // RelativeOffset is bytes; nVars and nTotalBackgroundVariables count doubles.
  // Keep this list in step with StoreNativeAmpsBackground: validating only B/U
  // could permit a later pressure/current write to overrun an allocated slice.
  for (std::size_t i=0;i<sizeof(fields)/sizeof(fields[0]);++i) {
    const auto& field=*fields[i];
    // U/B/n/T are mandatory; remaining fields are checked if allocated. Catch
    // every layout error before a first field can be written on ANY rank.
    const bool mandatory=i==0 || i==1 || i==2 || i==5;
    if ((mandatory && !field.allocate) || (field.allocate &&
        (field.RelativeOffset<0 || field.nVars!=components[i] ||
         static_cast<std::size_t>(field.RelativeOffset)+components[i]*sizeof(double)>
             static_cast<std::size_t>(nTotalBackgroundVariables)*sizeof(double))))
      return invalid("invalid DATAFILE field layout: "+std::string(field.VarList));
  }
  return Status::OK();
}
void ValidateNativeAmpsBackgroundLayout() {
  const auto status=NativeAmpsBackgroundLayoutStatus();
  if (!status.ok()) StopWithStatus("native AMPS background layout",status);
}

// Return the current DATAFILE storage slot and, when AMPS has initialized a
// distinct next slot, that slot as well. Runtime-owned getters deliberately
// bypass FILE time interpolation: both allocated slots mirror the same
// complete epoch for native output/readback compatibility, not a pair of file
// epochs to blend. DATAFILE may leave NextDataFileOffset negative; do not invent
// a second slot or file schedule when the allocator has not supplied one.
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
    // The provider reports proton T after total-pressure heating. Its explicit
    // partition preserves Te/Tp, so electron pressure must use that same scale,
    // not the unheated configured Te. n is electron density in both closures.
    electronPressurePa = sample.numberDensityM3 * SEP3D::Core::Const::k_B *
        parker.electronTemperatureK * (sample.temperatureK / parker.temperatureK);
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
SEP3D::Core::Status NativeAmpsBackgroundLayoutStatus() { return SEP3D::Core::Status::OK(); }
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

SEP3D::Core::Status ValidateBackgroundCandidateCollectively(
    const std::shared_ptr<const SEP3D::Background::BackgroundSnapshot>&,
    const std::vector<AmpsCellReference>&,
    const std::shared_ptr<SEP3D::Turbulence::TurbulenceProvider>&,
    SEP3D::Core::Status,std::vector<SEP3D::Turbulence::TurbulenceSample>*);

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
    // A typed host may install a provider first; otherwise the named factory
    // selects one. The rest of initialization uses only BackgroundProvider.
    Core::Status status;
    if (!gRuntimeBackgroundProvider)
      status = RuntimeModel::CreateBackgroundProvider(Configuration(),
                                                       &gRuntimeBackgroundProvider);
    const double initialEpochS = ApplicationRuntime().CurrentTimeS();
    if (status.ok()) status = gRuntimeBackgroundProvider->Prepare(initialEpochS);
    if (status.ok()) {
      Background::BackgroundSnapshotBuilder builder;
      status = builder.Build(*gRuntimeBackgroundProvider, positions, &snapshot);
    }
    // Join provider/build failures before dereferencing the candidate. The
    // complete turbulence/grid validation follows after both are prepared.
    int local=status.ok()?1:0,global=0;
    MPI_Allreduce(&local,&global,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
    if (!global) StopWithStatus("runtime background initialization",status.ok()?
        Core::Status(Core::StatusCode::SnapshotUnavailable,"another rank rejected initial background"):status);
    if (gPendingRestart) {
      // Reconstruct numerical fields, then synchronize the provider counter
      // with the checkpoint before exposing restored snapshot tags. Unsupported
      // model restart hooks reject collectively rather than relabeling old data.
      status=gRuntimeBackgroundProvider->RestorePreparedGeneration(gPendingRestart->backgroundGeneration);
      int localRestore=status.ok()?1:0,globalRestore=0;
      MPI_Allreduce(&localRestore,&globalRestore,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
      if (!globalRestore) StopWithStatus("background restart generation",status.ok()?
          Core::Status(Core::StatusCode::SnapshotUnavailable,"another rank rejected background restoration"):status);
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
  // The common gate below joins preparation failures on all MPI ranks.

  std::vector<Turbulence::TurbulenceSample> initialWaves;
  status=ValidateBackgroundCandidateCollectively(snapshot,cells,turbulence,status,&initialWaves);
  if (!status.ok()) StopWithStatus("initial background collective validation",status);

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
    const Turbulence::TurbulenceSample& waves = initialWaves[i];

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
  if (Configuration().options().background !=
      RuntimeModel::BackgroundAuthority::Swmf) {
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
  // This initialization product must use AMPS' cut-cell-aware writer: the
  // single-file brick writer emits complete Cartesian cells across the solar
  // surface even when their physical measures exclude the solid interior.
  // All ranks write their owner-local fragments, then AMPS assembles the
  // complete FEBRICK/tetrahedron zones into the requested filename. This
  // explicit call is independent of _PIC_OUTPUT_MODE_, which controls regular
  // sampling output. Allocation, cut-cell measures and halo exchange have
  // already completed; this writer must never run during preallocation mesh
  // output or inside a rank-zero-only branch. PrintMeshData must be true so
  // initialized background, turbulence and species numerics are emitted.
  PIC::Mesh::mesh->SetAssembleDistributedOutputFileFlag(true);
  for (const auto& species : gCompiledSpecies) {
    const std::string path = InitializationDataPath(
        Configuration().options().initializationDataTecplotFile,
        species.ampsIndex, PIC::nTotalSpecies);
    PIC::Mesh::mesh->OutputDistributedDataTECPLOT(
        path.c_str(), true, species.ampsIndex);
  }
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
}

// This gate is shared by initial publication and every subsequent provider.
// All ranks enter it even when their candidate failed or their pruned mesh owns
// no physical cells. Numerical/identity failure occurs BEFORE any live cell is
// written. Thus a late bad cell cannot leave a half-new/half-old background.
SEP3D::Core::Status ValidateBackgroundCandidateCollectively(
    const std::shared_ptr<const SEP3D::Background::BackgroundSnapshot>& candidate,
    const std::vector<AmpsCellReference>& cells,
    const std::shared_ptr<SEP3D::Turbulence::TurbulenceProvider>& turbulence,
    SEP3D::Core::Status status,
    std::vector<SEP3D::Turbulence::TurbulenceSample>* waves) {
  using namespace SEP3D;
  if (status.ok() && (!candidate || !turbulence ||
      candidate->positions().size()!=cells.size() || candidate->samples().size()!=cells.size()))
    status=Core::Status(Core::StatusCode::LayoutMismatch,"candidate owner-cell grid is incomplete");
  if (status.ok() && (!RuntimeModel::BackgroundAuthorityMatches(
      Configuration().options().background,candidate->metadata().provider) ||
      candidate->metadata().coordinateFrame!=Configuration().options().coordinateFrame ||
      !candidate->Covers(ApplicationRuntime().CurrentTimeS())))
    status=Core::Status(Core::StatusCode::ConfigurationConflict,"background authority or time coverage mismatch");
  if (status.ok())status=NativeAmpsBackgroundLayoutStatus();
  if (status.ok()) {
    // Turbulence is part of the same candidate. Evaluate into scratch, prove
    // directional variances are nonnegative and sum to the total, and require
    // positive prescribed support unless ballistic transport was explicit.
    // No Store* call is legal inside this validation pass.
    waves->reserve(cells.size());
    for (std::size_t i=0;i<cells.size();++i) {
      if (!SamePosition(cells[i].positionM,candidate->positions()[i]) || !cells[i].cell) {
        status=Core::Status(Core::StatusCode::LayoutMismatch,"background cell order/pointer changed");break;
      }
      status=Background::ValidateCompleteSample(candidate->samples()[i],candidate->capabilities());
      if (!status.ok())break;
      auto w=turbulence->Evaluate(cells[i].positionM,candidate->samples()[i]);
      const double total=w.deltaBPlus2T2+w.deltaBMinus2T2;
      const double scale=std::max(std::numeric_limits<double>::min(),std::max(std::fabs(total),std::fabs(w.deltaB2T2)));
      if (!w.status.usable() || !w.valid || !std::isfinite(total) ||
          w.deltaBPlus2T2<0 || w.deltaBMinus2T2<0 || !std::isfinite(w.deltaB2T2) ||
          std::fabs(total-w.deltaB2T2)>128*std::numeric_limits<double>::epsilon()*scale ||
          (Configuration().options().turbulence==RuntimeModel::TurbulenceAuthority::Prescribed && !w.ballistic && total<=0)) {
        status=Core::Status(Core::StatusCode::BackgroundInvalid,"candidate directional turbulence is incomplete");break;
      }
      waves->push_back(w);
    }
  }
  // All ranks reach this readiness reduction, including failed providers and
  // ranks with zero owned cells. Stop here collectively before dereferencing
  // metadata or entering the later reductions if any local candidate failed.
  int local=status.ok()?1:0,global=0;
  MPI_Allreduce(&local,&global,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  if (!global) return status.ok()?Core::Status(Core::StatusCode::SnapshotUnavailable,
      "another MPI rank rejected the background candidate"):status;
  const auto& m=candidate->metadata();
  // Compare a compact deterministic metadata hash plus the exact generation.
  // Cell arrays differ by MPI ownership and must not enter this identity hash.
  // Seventeen digits preserve double-valued epochs/intervals in the encoding.
  std::ostringstream identity;identity<<std::setprecision(17)<<m.epochS<<'|'<<m.validFromS<<'|'<<m.validUntilS
      <<'|'<<m.coordinateFrame<<'|'<<m.providerIdentity<<'|'<<m.configurationFingerprint
      <<'|'<<Configuration().physics_fingerprint();
  unsigned long long send[2]={m.generation,Fnv1a64(identity.str())},lo[2]={},hi[2]={};
  MPI_Allreduce(send,lo,2,MPI_UNSIGNED_LONG_LONG,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(send,hi,2,MPI_UNSIGNED_LONG_LONG,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  unsigned long long localCells=cells.size(),globalCells=0;
  MPI_Allreduce(&localCells,&globalCells,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  if (!globalCells || lo[0]!=hi[0] || lo[1]!=hi[1])
    return Core::Status(Core::StatusCode::ConfigurationConflict,"MPI background coverage, epoch or identity disagree");
  return Core::Status::OK();
}

// Native tests read the SAME bytes consumed by AMPS, including both DATAFILE
// slots. They never call Store* or Prepare. Optional primitive/transport/E
// slots listed below are checked when allocated; U/B/density/temperature are
// mandatory. Derived native current/electron-pressure slots are filled by the
// bridge but are not separately compared by this readback helper.
bool MeshBackgroundBytesMatch(PIC::Mesh::cDataCenterNode* cell,
    const SEP3D::Background::BackgroundSample& s) {
  using namespace SEP3D;
  const auto& layout=Configuration().storage_layout();
  bool matches=true;
  auto application=[&](std::size_t offset,const double* expected,std::size_t n) {
    // Optional application tensors may be disabled by [storage]. Check every
    // allocated component; immutable snapshots retain the complete sample.
    if (offset==RuntimeModel::kNoOffset)return;
    std::vector<double> actual(n);LoadBytes(cell,offset,actual.data(),n*sizeof(double));
    matches=matches&&std::memcmp(actual.data(),expected,n*sizeof(double))==0;
  };
  double B[3],U[3],curvature[3];s.B.CopyTo(B);s.U.CopyTo(U);s.curvature.CopyTo(curvature);
  application(layout.magneticFieldOffset,B,3);application(layout.bulkVelocityOffset,U,3);
  application(layout.numberDensityOffset,&s.numberDensityM3,1);
  application(layout.temperatureOffset,&s.temperatureK,1);application(layout.pressureOffset,&s.pressurePa,1);
  application(layout.alfvenSpeedOffset,&s.alfvenSpeedMpS,1);application(layout.velocityDivergenceOffset,&s.divU,1);
  application(layout.magneticGradientOffset,&s.gradB.m[0][0],9);application(layout.velocityGradientOffset,&s.gradU.m[0][0],9);
  application(layout.divBhatOffset,&s.divBhat,1);application(layout.focusingLengthOffset,&s.focusingLenM,1);
  application(layout.curvatureOffset,curvature,3);application(layout.fieldAlignedStrainOffset,&s.fieldAlignedStrain,1);
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__DATAFILE_
  using namespace PIC::CPLR::DATAFILE;
  auto native=[&](const cOffsetElement& field,const double* expected,int n,bool required=false) {
    // Byte equality is appropriate for publication evidence: these are direct
    // copies, not a numerical-fit test. Verify both allocated time slots so
    // stale next-slot data cannot hide behind a correct current-slot value.
    if (!field.allocate) { if(required)matches=false;return; }
    if (field.nVars!=n || field.RelativeOffset<0) {matches=false;return;}
    int slots[2];const int count=NativeAmpsDataSlots(slots);
    for (int k=0;k<count;++k) {
      const char* bytes=cell->GetAssociatedDataBufferPointer()+CenterNodeAssociatedDataOffsetBegin+slots[k]+field.RelativeOffset;
      matches=matches&&std::memcmp(bytes,expected,static_cast<std::size_t>(n)*sizeof(double))==0;
    }
  };
  const auto electric=s.U.Cross(s.B)*(-1);double E[3];electric.CopyTo(E);
  native(Offset::MagneticField,B,3,true);native(Offset::PlasmaBulkVelocity,U,3,true);
  native(Offset::PlasmaNumberDensity,&s.numberDensityM3,1,true);
  native(Offset::PlasmaTemperature,&s.temperatureK,1,true);
  native(Offset::PlasmaIonPressure,&s.pressurePa,1);
  native(Offset::ElectricField,E,3);native(Offset::PlasmaDivU,&s.divU,1);
  native(Offset::MagneticFieldGradient,&s.gradB.m[0][0],9);
#else
  matches=false;
#endif
  return matches;
}

// Produce a decomposition-independent value for one physical owner cell.
// Hexadecimal floating-point text preserves the exact binary value while
// avoiding structure padding and native endianness.  The corresponding native
// center-node bytes are independently compared above; hashing the canonical
// sample here is therefore equivalent to hashing those bytes after a PASS,
// while remaining independent of AMPS storage offsets and MPI ownership.
std::uint64_t BackgroundCellFingerprint(const AmpsCellReference& cell,
    const SEP3D::Background::BackgroundSample& sample) {
  std::ostringstream encoded;
  encoded << std::hexfloat
      << cell.positionM.x << '|' << cell.positionM.y << '|'
      << cell.positionM.z << '|' << sample.B.x << '|' << sample.B.y << '|'
      << sample.B.z << '|' << sample.U.x << '|' << sample.U.y << '|'
      << sample.U.z << '|' << sample.numberDensityM3 << '|'
      << sample.temperatureK << '|' << sample.pressurePa << '|'
      << sample.alfvenSpeedMpS << '|' << sample.divU << '|'
      << sample.divBhat << '|' << sample.focusingLenM << '|'
      << sample.curvature.x << '|' << sample.curvature.y << '|'
      << sample.curvature.z << '|' << sample.fieldAlignedStrain << '|'
      << sample.generation << '|' << sample.configurationDigest;
  for (int component=0;component<3;++component)
    for (int coordinate=0;coordinate<3;++coordinate)
      encoded << '|' << sample.gradB(component,coordinate);
  for (int component=0;component<3;++component)
    for (int coordinate=0;coordinate<3;++coordinate)
      encoded << '|' << sample.gradU(component,coordinate);
  return Fnv1a64(encoded.str());
}
// Read-only local evidence capture; its caller reduces flags/counts over MPI.
// Ghost coverage is one physical center per received active block, not every
// ghost cell. Absence of received blocks remains a test prerequisite SKIP.
void CaptureRuntimeMeshBackground(bool* owned,bool* ghosts,bool* provider,
                                  unsigned long long* ghostCount,
                                  unsigned long long* ownerFingerprintXor,
                                  unsigned long long* ownerFingerprintSum) {
  using namespace SEP3D;
  *owned=*ghosts=*provider=true;*ghostCount=0;
  *ownerFingerprintXor=*ownerFingerprintSum=0;
  if (!gInstalledBackground || !gRuntimeBackgroundProvider) {*owned=*provider=false;return;}
  const auto cells=CollectOwnedPhysicalCells();const auto& samples=gInstalledBackground->samples();
  if(cells.size()!=samples.size()){*owned=*provider=false;return;}
  for (std::size_t i=0;i<cells.size();++i) {
    *owned=*owned&&MeshBackgroundBytesMatch(cells[i].cell,samples[i]);
    const unsigned long long fingerprint=static_cast<unsigned long long>(
        BackgroundCellFingerprint(cells[i],samples[i]));
    *ownerFingerprintXor^=fingerprint;
    *ownerFingerprintSum+=fingerprint;
  }
  const auto* prepared=gRuntimeBackgroundProvider->PreparedMetadata();
  *provider=prepared && prepared->epochS==gInstalledBackground->metadata().epochS &&
      prepared->generation==gInstalledBackground->metadata().generation && gNativeAmpsBackgroundReady;
  // A remote block may allocate interior centers that were NEVER sent to this
  // rank. AMPS InitLayerBlock sends only selected face/edge/corner layers.
  // Walk the actual receive plan and honor its center-node mask; selecting the
  // first allocated center instead would falsely report a stale halo on the
  // positive faces. Do not initiate another exchange to hide a missed epoch.
  // RecvNodeTableLength is a per-round work counter and is reset after unpack;
  // GlobalSendTable retains the complete per-sender receive block count.
  const auto& exchange=PIC::Mesh::mesh->ParallelBlockDataExchangeData;
  bool reportedMismatch=false;
  if (exchange.GlobalSendTable && exchange.RecvNodeTable) {
    for(int from=0;from<PIC::nTotalThreads;++from) {
      if(from==PIC::ThisThread)continue;
      const int count=exchange.GlobalSendTable[PIC::ThisThread+from*PIC::nTotalThreads];
      if(count>0 && !exchange.RecvNodeTable[from]) { *ghosts=false;continue; }
      for(int block=0;block<count;++block) {
        auto* node=exchange.RecvNodeTable[from][block];
        if(!node || !node->block || !node->IsUsedInCalculationFlag)continue;
        unsigned char* mask=nullptr;
        if(exchange.RecvCenterNodePackingTable && exchange.RecvCenterNodePackingTable[from] &&
           exchange.BlockCenterNodeSendMaskLength>0)
          mask=exchange.RecvCenterNodePackingTable[from]+block*exchange.BlockCenterNodeSendMaskLength;
        bool checked=false;
        for(int k=0;k<_BLOCK_CELLS_Z_ && !checked;++k)
          for(int j=0;j<_BLOCK_CELLS_Y_ && !checked;++j)
            for(int i=0;i<_BLOCK_CELLS_X_ && !checked;++i) {
              if(mask && !PIC::Mesh::BlockElementSendMask::CenterNode::Test(i,j,k,mask))continue;
              Core::Vec3 at(node->xmin[0]+(i+0.5)*((node->xmax[0]-node->xmin[0])/_BLOCK_CELLS_X_),
                  node->xmin[1]+(j+0.5)*((node->xmax[1]-node->xmin[1])/_BLOCK_CELLS_Y_),
                  node->xmin[2]+(k+0.5)*((node->xmax[2]-node->xmin[2])/_BLOCK_CELLS_Z_));
              const double r=(at-Configuration().options().coordinateOriginM).Norm();
              if(r<Configuration().options().innerRadiusM || r>Configuration().options().outerRadiusM)continue;
              auto* cell=node->block->GetCenterNode(PIC::Mesh::mesh->getCenterNodeLocalNumber(i,j,k));
              if(!cell)continue;
              const auto expected=gRuntimeBackgroundProvider->Evaluate(at);
              const bool matches=expected.valid&&MeshBackgroundBytesMatch(cell,expected);
              *ghosts=*ghosts&&matches; ++*ghostCount;checked=true;
              if(!matches && !reportedMismatch) {
                reportedMismatch=true;
                std::cerr<<"[srcSEP3D] received background mismatch: rank="<<PIC::ThisThread
                    <<" owner="<<from<<" block="<<block<<" cell=("<<i<<','<<j<<','<<k<<')'
                    <<" epoch_s="<<(prepared?prepared->epochS:-1)<<" generation="<<(prepared?prepared->generation:0)
                    <<" expected_valid="<<expected.valid<<" status="<<expected.status.message<<'\n';
              }
            }
      }
    }
  }
  // A nonempty owner rank also checks a freshly evaluated representative,
  // detecting a snapshot that is merely self-consistent but at a stale epoch.
  if (!cells.empty()) {
    const auto expected=gRuntimeBackgroundProvider->Evaluate(cells.front().positionM);
    *provider=*provider&&expected.valid&&MeshBackgroundBytesMatch(cells.front().cell,expected);
  }
}

// Provider-neutral cadence transaction. The active snapshot remains mover
// authority until validation, descriptor staging, all writes and halo work
// finish. Typed candidate/staging failures abort before live bytes change;
// this is not a recover-and-continue rollback after a write/MPI failure.
void RefreshBackgroundAtBoundary() {
  using namespace SEP3D;
  RuntimeModel::Runtime& runtime=ApplicationRuntime();
  if (!runtime.EventDue(RuntimeModel::ScheduledEvent::Background)) return;
  const auto cells=CollectOwnedPhysicalCells();
  std::vector<Core::Vec3> positions;positions.reserve(cells.size());
  for (const auto& cell:cells)positions.push_back(cell.positionM);
  std::shared_ptr<const Background::BackgroundSnapshot> candidate;
  std::shared_ptr<Turbulence::TurbulenceProvider> turbulence;
  Core::Status status;
  const double epochS=runtime.CurrentTimeS();
  if (Configuration().options().background!=RuntimeModel::BackgroundAuthority::Swmf) {
    // Runtime models own Prepare/Evaluate; a coupled SWMF host instead stages
    // its next complete snapshot/turbulence pair before reaching this boundary.
    if (!gRuntimeBackgroundProvider)
      status=Core::Status(Core::StatusCode::SnapshotUnavailable,"runtime background provider ownership was lost");
    else {
      status=gRuntimeBackgroundProvider->Prepare(epochS);
      if (status.ok()) {
        Background::BackgroundSnapshotBuilder builder;
        status=builder.Build(*gRuntimeBackgroundProvider,positions,&candidate);
      }
    }
    turbulence=gInstalledTurbulence;
  } else {
    candidate=gStagedBackground;turbulence=gStagedTurbulence;
    if (!candidate || !turbulence)
      status=Core::Status(Core::StatusCode::SnapshotUnavailable,"SWMF cadence has no staged complete state");
  }
  if (status.ok())status=turbulence->Prepare(epochS);
  std::vector<Turbulence::TurbulenceSample> waves;
  status=ValidateBackgroundCandidateCollectively(candidate,cells,turbulence,status,&waves);
  if (!status.ok())StopWithStatus("collective background candidate",status);
  const auto& m=candidate->metadata();
  const auto& previous=gInstalledBackground->metadata();
  // A cadence update may advance time/generation but cannot switch physical
  // authority, ownership or frozen configuration halfway through a run.
  if (m.provider!=previous.provider || m.ownership!=previous.ownership ||
      m.providerIdentity!=previous.providerIdentity ||
      m.configurationFingerprint!=previous.configurationFingerprint ||
      m.epochS!=epochS || m.generation<=previous.generation)
    StopWithStatus("background identity continuity",Core::Status(
        Core::StatusCode::ConfigurationConflict,"background changed identity or did not advance epoch/generation"));
  status=runtime.RequestSnapshotUpdate(m.epochS,m.generation);
  if (status.ok())status=runtime.BeginSnapshotFill();
  RuntimeModel::SnapshotDescriptor descriptor;
  descriptor.authority=Configuration().options().background;
  descriptor.epochS=m.epochS;descriptor.validFromS=m.validFromS;descriptor.validUntilS=m.validUntilS;
  descriptor.generation=m.generation;descriptor.complete=true;
  descriptor.coordinateFrame=m.coordinateFrame;descriptor.providerIdentity=m.providerIdentity;
  descriptor.configurationFingerprint=Configuration().physics_fingerprint();
  if (status.ok())status=runtime.StageSnapshot(descriptor);
  // Stage validation is rank-local too; join it before touching the live cache.
  int local=status.ok()?1:0,global=0;
  MPI_Allreduce(&local,&global,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  if (!global) {
    (void)runtime.FailSnapshotUpdate("collective Runtime staging rejected");
    StopWithStatus("background staging",status.ok()?Core::Status(
        Core::StatusCode::SnapshotUnavailable,"another rank rejected Runtime staging"):status);
  }
  gNativeAmpsBackgroundReady=false;
  // All fallible numerical/layout/staging checks have joined. Write the two
  // independent center-node slices and matching waves from the same candidate.
  for (std::size_t i=0;i<cells.size();++i) {
    StoreBackground(cells[i].cell,candidate->samples()[i]);
    StoreNativeAmpsBackground(cells[i].cell,candidate->samples()[i]);
    (void)StoreTurbulenceAtCellCenter(cells[i].cell,waves[i]);
  }
  // Both DATAFILE slots receive the identical frozen epoch. AMPS must not
  // temporally interpolate a new magnetic field with an old plasma velocity.
  // Ghost exchange and any derived guiding-centre fields complete before the
  // Runtime or a mover may advertise/consume the new generation.
  CompleteNativeAmpsBackgroundInstallation();
  status=runtime.PublishStagedSnapshot(true);
  if (!status.ok())StopWithStatus("background commit",status);
  gInstalledBackground=candidate;gInstalledTurbulence=turbulence;
  gStagedBackground.reset();gStagedTurbulence.reset();
  ++gBackgroundPublishedUpdates;
  if (PIC::ThisThread==0)
    std::cout<<"[srcSEP3D] background mesh update: provider="<<m.providerIdentity
             <<" tick="<<runtime.counters().currentTick<<" time_s="<<m.epochS
             <<" generation="<<m.generation<<" halo=ready\n";
}



struct PackedObservation {
  std::uint64_t stableId, cellId;
  std::int32_t species;
  double x, y, z, momentum, mass, mu, gyrophase, weight;
  std::uint64_t completedStep, substep, lastShockGeneration;
  double remainingScatteringOpticalDepth;
  std::uint64_t nextScatteringEvent;
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
  // The source-free control normally has no particle observers. Its native
  // shock telemetry is published by the driver after every completed step;
  // avoid gathering the entire AMR volume for an empty observer request.
  if (Configuration().options().intent == RuntimeModel::RunIntent::ShockPropagation &&
      Configuration().options().observers.empty()) return;

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
                particle.lastShockGeneration,
                particle.remainingScatteringOpticalDepth,
                particle.nextScatteringEvent});
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
      metadata.parallelDiffusionModelId =
          Configuration().options().parallelDiffusionModelId;
      metadata.parallelDiffusionConfigurationFingerprint =
          Configuration().options().parallelDiffusionConfigurationFingerprint;
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
                particle.lastShockGeneration,
                particle.remainingScatteringOpticalDepth,
                particle.nextScatteringEvent});
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
    // The serialized identity keeps resolved event/physics inputs and cadence
    // semantics, but excludes output and checkpoint filenames.  A resumed
    // native run must be able to write to a distinct evidence directory; the
    // full resolved manifest remains available in publication provenance.
    checkpoint.resolvedConfigurationManifest =
        Configuration().restart_compatibility_manifest();
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
      particle.remainingScatteringOpticalDepth =
          packed.remainingScatteringOpticalDepth;
      particle.nextScatteringEvent = packed.nextScatteringEvent;
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

SEP3D::Core::Status WriteReducedProductionOutputAtBoundaryImpl() {
  using namespace SEP3D;
  namespace SF=SEP::CoronaSwcme::ShockFront;
  const auto& options=Configuration().options();
  if(options.intent!=RuntimeModel::RunIntent::ShockPropagation||
      options.background!=RuntimeModel::BackgroundAuthority::RuntimeModel||
      options.backgroundModelId!="sep-corona-swcme-shock-front-v1")
    return Core::Status::OK();

  const auto reduced=std::dynamic_pointer_cast<
      Adapters::ShockFrontBackgroundAdapter>(gRuntimeBackgroundProvider);
  const auto epoch=reduced?reduced->FrontEpoch():nullptr;
  const auto provider=reduced?reduced->SharedProvider():nullptr;
  const auto* active=ApplicationRuntime().active_snapshot();
  if(!reduced||!epoch||!provider||!active||!active->complete)
    return Core::Status(Core::StatusCode::SnapshotUnavailable,
        "reduced output requires committed front, ambient and runtime epochs");

  const std::uint64_t tick=ApplicationRuntime().counters().currentTick;
  const double timeS=ApplicationRuntime().CurrentTimeS();
  const double dt=options.requestedTimeStepS;
  const auto handoff=provider->HandoffTimeS();
  const auto endpoint=provider->EndpointTimeS();
  if(!handoff.ok()||!endpoint.ok())return Core::Status(
      Core::StatusCode::BackgroundInvalid,"reduced output landmarks are unavailable");

  // Regular cadence supplies heliospheric context.  Landmark windows add
  // launch acceleration and one sample immediately before/at/after the exact
  // handoff and observer passage.  This is output scheduling only: the front
  // and ambient are always evaluated at the committed native epoch.
  auto near=[&](double landmark,double radius) {
    // Handoff/endpoint times are roots of an analytical trajectory, whereas
    // native boundaries are integer multiples of dt.  A root can consequently
    // lie a few ulps to one side of the mathematically exact boundary.  Use a
    // relative 1e-9 clock tolerance (10.2 microseconds at this handoff), far
    // below both the 600-s host step and any resolved physical time scale, so
    // the requested symmetric pre/at/post window cannot lose only its earlier
    // member through roundoff.  This tolerance selects output times only; it
    // neither changes the trajectory nor admits a shock state.
    const double clockTolerance=1.0e-9*std::max(1.0,std::fabs(landmark));
    return std::fabs(timeS-landmark)<=radius*dt+clockTolerance;
  };
  const bool due=tick==0||tick%options.outputCadenceSteps==0||
      near(0.5*provider->Event().accelerationDurationS,0.5)||
      near(provider->Event().accelerationDurationS,0.5)||
      near(handoff.value,1.0)||near(endpoint.value,1.0);
  if(!due||gLastReducedOutputTick==tick)return Core::Status::OK();

  // The surface and volume products are useful only when they describe the
  // same immutable publication.  Check clocks/generations before entering the
  // expensive AMPS writer; a mismatch is an error, never a best-effort file.
  bool localIdentity=active->epochS==timeS&&active->generation==epoch->generation&&
      epoch->ambientGeneration==epoch->generation&&
      epoch->eventIdentity==provider->Event().physicsFingerprint;
  int localOK=localIdentity?1:0,allOK=0;
  MPI_Allreduce(&localOK,&allOK,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  if(!allOK)return Core::Status(Core::StatusCode::ConfigurationConflict,
      "reduced front, ambient and native runtime epochs disagree");

  bool owner=false,ghost=false,providerMatch=false;
  unsigned long long localGhosts=0,localOwnerXor=0,localOwnerSum=0;
  CaptureRuntimeMeshBackground(&owner,&ghost,&providerMatch,&localGhosts,
      &localOwnerXor,&localOwnerSum);
  int flags[3]={owner?1:0,ghost?1:0,providerMatch?1:0};
  int globalFlags[3]={};
  MPI_Allreduce(flags,globalFlags,3,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  unsigned long long globalGhosts=0,globalOwnerXor=0,globalOwnerSum=0;
  MPI_Allreduce(&localGhosts,&globalGhosts,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,
      MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&localOwnerXor,&globalOwnerXor,1,MPI_UNSIGNED_LONG_LONG,
      MPI_BXOR,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&localOwnerSum,&globalOwnerSum,1,MPI_UNSIGNED_LONG_LONG,
      MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  if(!globalFlags[0]||!globalFlags[1]||!globalFlags[2])return Core::Status(
      Core::StatusCode::LayoutMismatch,
      "native owner/received-ghost/provider readback failed before output");

  const auto localBySpecies=CountLocalParticlesBySpecies();
  unsigned long long localParticles=0,globalParticles=0;
  for(auto count:localBySpecies)localParticles+=count;
  MPI_Allreduce(&localParticles,&globalParticles,1,MPI_UNSIGNED_LONG_LONG,
      MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  unsigned long long localInjected=0,globalInjected=0;
  for(const auto& row:gSourceLedger)localInjected+=row.macroparticles;
  MPI_Allreduce(&localInjected,&globalInjected,1,MPI_UNSIGNED_LONG_LONG,
      MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  if(globalParticles!=0||globalInjected!=0)return Core::Status(
      Core::StatusCode::ConfigurationConflict,
      "reduced background output observed particles or source allocation");

  std::ostringstream stem;
  stem<<options.outputPrefix<<"-tick-"<<std::setw(8)<<std::setfill('0')<<tick;
  const fs::path directory(options.outputDirectory);
  int ioOK=1;
  std::string failure;
  if(PIC::ThisThread==0) {
    std::error_code error;fs::create_directories(directory,error);
    if(error) {ioOK=0;failure="cannot create reduced output directory: "+error.message();}
  }
  MPI_Bcast(&ioOK,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);
  if(!ioOK)return Core::Status(Core::StatusCode::Error,
      PIC::ThisThread==0?failure:"rank zero could not create output directory");
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

  // This maintained distributed writer traverses actual AMPS cells, invokes
  // the registered application interpolation callback, and assembles owner
  // fragments.  It therefore publishes installed native ambient values, not
  // a second provider sampling performed solely for visualization.
  PIC::Mesh::mesh->SetAssembleDistributedOutputFileFlag(true);
  for(const auto& species:gCompiledSpecies) {
    const std::string base=(directory/(stem.str()+"-ambient.dat")).string();
    const std::string path=InitializationDataPath(base,species.ampsIndex,
        PIC::nTotalSpecies);
    PIC::Mesh::mesh->OutputDistributedDataTECPLOT(path.c_str(),true,
        species.ampsIndex);
  }
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

  if(PIC::ThisThread==0) {
    const std::string surface=Output::SerializeReducedFrontTecplot(
        *epoch,provider->Event());
    if(surface.empty()) {ioOK=0;failure="front surface serializer rejected its quadrature layout";}
    auto write=[&](const fs::path& path,const std::string& bytes) {
      if(!ioOK)return;
      std::ofstream output(path,std::ios::binary);
      output<<bytes;output.close();
      if(!output){ioOK=0;failure="cannot write/close '"+path.string()+"'";}
    };
    write(directory/(stem.str()+"-front.dat"),surface);
    write(directory/(stem.str()+"-front.json"),SF::SerializeEpochJson(*epoch)+"\n");

    std::string observerStatus="not-reached";
    bool observerAccepted=false;
    double observerMach=std::numeric_limits<double>::quiet_NaN();
    if(timeS+1e-8>=endpoint.value) {
      const auto passage=provider->EvaluateFrontPoint(
          provider->Event().observerPositionM,endpoint.value,UINT64_C(1));
      if(passage.ok()) {
        observerStatus=SF::Name(passage.value.status);
        observerAccepted=passage.value.status==SF::FrontStatus::SolvedFastShock;
        observerMach=passage.value.fastMach;
      } else {ioOK=0;failure="cannot evaluate exact observer passage: "+
          passage.status.message;}
    }
    if(ioOK) {
      std::ostringstream receipt;receipt<<std::setprecision(17)
        <<"{\n  \"schema\": \"srcsep3d-reduced-background-output-v1\",\n"
        <<"  \"time_s\": "<<timeS<<",\n  \"tick\": "<<tick
        <<",\n  \"background_generation\": "<<active->generation
        <<",\n  \"front_generation\": "<<epoch->generation
        <<",\n  \"event_identity\": "<<std::quoted(epoch->eventIdentity)
        <<",\n  \"phase\": "<<std::quoted(SF::Name(epoch->trajectory.phase))
        <<",\n  \"apex_radius_m\": "<<epoch->trajectory.apexRadiusM
        <<",\n  \"apex_speed_m_s\": "<<epoch->trajectory.apexSpeedMPerS
        <<",\n  \"mesh_volume_role\": \"ambient-reference-only\",\n"
        <<"  \"surface_downstream_role\": \"immediate-rh-limit-only\",\n"
        <<"  \"owner_readback_match\": true,\n"
        <<"  \"received_ghost_readback_match\": true,\n"
        <<"  \"provider_epoch_match\": true,\n"
        <<"  \"owner_fingerprint_xor\": "<<globalOwnerXor
        <<",\n  \"owner_fingerprint_sum\": "<<globalOwnerSum
        <<",\n  \"received_ghost_cells_checked\": "<<globalGhosts
        <<",\n  \"particle_count\": "<<globalParticles
        <<",\n  \"injected_particle_count\": "<<globalInjected
        <<",\n  \"geometric_observer_arrival\": "
        <<(timeS+1e-8>=endpoint.value?"true":"false")
        <<",\n  \"exact_observer_arrival_time_s\": "<<endpoint.value
        <<",\n  \"observer_shock_status\": "<<std::quoted(observerStatus)
        <<",\n  \"observer_shock_accepted\": "
        <<(observerAccepted?"true":"false")<<",\n  \"observer_fast_mach\": ";
      if(std::isfinite(observerMach))receipt<<observerMach;else receipt<<"null";
      receipt<<"\n}\n";
      write(directory/(stem.str()+"-receipt.json"),receipt.str());
    }
  }
  MPI_Bcast(&ioOK,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);
  if(!ioOK)return Core::Status(Core::StatusCode::Error,
      PIC::ThisThread==0?failure:"rank zero could not publish reduced products");
  gLastReducedOutputTick=tick;
  return Core::Status::OK();
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

SEP3D::Core::Status SEP3D::WriteReducedProductionOutputAtBoundary() {
  return WriteReducedProductionOutputAtBoundaryImpl();
}

SEP3D::Core::Status SEP3D::ConfigureApplication(
    const std::shared_ptr<const RuntimeModel::RunConfiguration3D>& configuration) {
  return ApplicationRuntime().Configure(configuration);
}

SEP3D::Core::Status SEP3D::InstallApplicationInputFile(
    const std::string& path) {
  if (path.empty())
    return Core::Status(Core::StatusCode::InvalidInput,
                        "application input path is empty");
  if (ApplicationRuntime().state() !=
      RuntimeModel::LifecycleState::Configured) {
    return Core::Status(Core::StatusCode::InvalidTransition,
                        "install application input after provisional configuration and before mesh setup");
  }
  if (!gApplicationInputPath.empty())
    return Core::Status(Core::StatusCode::InvalidTransition,
                        "application input path was already installed");
  gApplicationInputPath = path;
  return Core::Status::OK();
}

SEP3D::Core::Status SEP3D::InstallBackgroundProvider(
    const std::shared_ptr<Background::BackgroundProvider>& provider) {
  using namespace SEP3D;
  // Host injection is legal after immutable configuration (and optionally
  // mesh binding), before acquisition. Retain shared ownership for subsequent
  // cadence preparation; an imported SWMF snapshot follows its separate API.
  const auto state=ApplicationRuntime().state();
  if (!provider || !ApplicationRuntime().configuration() ||
      (state!=RuntimeModel::LifecycleState::Configured && state!=RuntimeModel::LifecycleState::MeshReady))
    return Core::Status(Core::StatusCode::InvalidTransition,"install a non-null background provider before acquisition");
  if (Configuration().options().background==RuntimeModel::BackgroundAuthority::Swmf)
    return Core::Status(Core::StatusCode::ConfigurationConflict,"SWMF authority uses host-installed snapshots");
  const auto status=provider->Validate();if(!status.ok())return status;
  if (gInstalledBackground || gRuntimeBackgroundProvider)
    return Core::Status(Core::StatusCode::InvalidTransition,"background ownership was already installed");
  gRuntimeBackgroundProvider=provider;return Core::Status::OK();
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
  if (!RuntimeModel::BackgroundAuthorityMatches(
          Configuration().options().background,snapshot->metadata().provider)) {
    return Core::Status(Core::StatusCode::ConfigurationConflict,
                        "installed snapshot authority differs from configuration");
  }
  const Background::SnapshotMetadata& metadata = snapshot->metadata();
  if (snapshot->positions().size() != snapshot->samples().size() ||
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
  // The hook itself does not inspect argc/argv or open a file.  A standalone
  // shared-file run invokes ParseInstalledApplicationInput only after this
  // hook has completed; coupled entry points install typed options directly.
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
    // OutputDistributedDataTECPLOT prints brick and cut-cell vertex records
    // through temporary center-node objects. This hook is as essential as the print
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

// AMPS owns particle deletion after an internal-boundary callback returns
// _PARTICLE_DELETED_ON_THE_FACE_.  Calling DeleteParticle here as well would
// double-release the particle-buffer slot, so this callback intentionally has
// no side effect and returns only the disposition code.
int AbsorbParticleAtSolarSurface(int species, long int particle,
                                 double* position, double* velocity,
                                 double& remainingTimeS, void* node,
                                 void* sphere) {
  (void)species;
  (void)particle;
  (void)position;
  (void)velocity;
  (void)remainingTimeS;
  (void)node;
  (void)sphere;
  return _PARTICLE_DELETED_ON_THE_FACE_;
}

void RegisterSolarSurfaceBoundary() {
  if (gSolarSurfaceBoundary != nullptr) {
    StopWithStatus("solar internal-boundary registration",
        SEP3D::Core::Status(
            SEP3D::Core::StatusCode::InvalidTransition,
            "the solar photosphere was registered more than once"));
  }

  const SEP3D::Mesh::SolarBoundaryGeometry geometry =
      SEP3D::Mesh::MakeSolarBoundary(Configuration().options());
  const double localGeometry[4] = {
      geometry.centerM.x, geometry.centerM.y, geometry.centerM.z,
      geometry.radiusM};
  double minimumGeometry[4] = {};
  double maximumGeometry[4] = {};
  MPI_Allreduce(localGeometry, minimumGeometry, 4, MPI_DOUBLE, MPI_MIN,
                MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(localGeometry, maximumGeometry, 4, MPI_DOUBLE, MPI_MAX,
                MPI_GLOBAL_COMMUNICATOR);
  for (int component = 0; component < 4; ++component) {
    if (!std::isfinite(localGeometry[component]) ||
        minimumGeometry[component] != maximumGeometry[component]) {
      StopWithStatus("solar internal-boundary registration",
          SEP3D::Core::Status(
              SEP3D::Core::StatusCode::ConfigurationConflict,
              "MPI ranks disagree on the finite solar-boundary geometry"));
    }
  }
  if (!(geometry.radiusM > 0.0)) {
    StopWithStatus("solar internal-boundary registration",
        SEP3D::Core::Status(SEP3D::Core::StatusCode::InvalidInput,
                           "the physical solar radius is not positive"));
  }

  // This is the Venus registration sequence, placed after
  // PIC::Init_BeforeParser() but before the AMPS mesh exists.  The returned
  // descriptor is already inserted into PIC::Mesh::mesh by
  // RegisterInternalSphere(); a second explicit registration would duplicate
  // the surface in every cut-cell query.
  PIC::BC::InternalBoundary::Sphere::Init();
  const cInternalBoundaryConditionsDescriptor descriptor =
      PIC::BC::InternalBoundary::Sphere::RegisterInternalSphere();
  if (descriptor.BondaryType != _INTERNAL_BOUNDARY_TYPE_SPHERE_) {
    StopWithStatus("solar internal-boundary registration",
        SEP3D::Core::Status(SEP3D::Core::StatusCode::LayoutMismatch,
                           "AMPS returned a non-spherical boundary descriptor"));
  }
  gSolarSurfaceBoundary = static_cast<cInternalSphericalData*>(
      descriptor.BoundaryElement);
  if (gSolarSurfaceBoundary == nullptr) {
    StopWithStatus("solar internal-boundary registration",
        SEP3D::Core::Status(SEP3D::Core::StatusCode::LayoutMismatch,
                           "AMPS returned a null spherical boundary"));
  }

  double centerM[3] = {
      geometry.centerM.x, geometry.centerM.y, geometry.centerM.z};
  gSolarSurfaceBoundary->SetSphereGeometricalParameters(
      centerM, geometry.radiusM);
  // Use the exact application resolution callback at the surface.  For the
  // present refinement law r=R_sun is below the configurable source shell and
  // therefore clamps safely to mesh.solar.surface_cell_size_m.
  gSolarSurfaceBoundary->localResolution = localResolution;
  gSolarSurfaceBoundary->InjectionRate = nullptr;
  gSolarSurfaceBoundary->faceat = 0;
  gSolarSurfaceBoundary->ParticleSphereInteraction =
      AbsorbParticleAtSolarSurface;
  gSolarSurfaceBoundary->InjectionBoundaryCondition = nullptr;

  if (PIC::ThisThread == 0) {
    // Format through a temporary stream so this diagnostic cannot change the
    // caller's persistent std::cout precision/floatfield.
    std::ostringstream message;
    message << std::scientific << std::setprecision(17)
            << "[srcSEP3D] registered absorbing AMPS solar photosphere "
            << "center_m=(" << geometry.centerM.x << ','
            << geometry.centerM.y << ',' << geometry.centerM.z << ')'
            << " radius_m=" << geometry.radiusM
            << "; domain.inner_radius_m remains the independent "
               "Parker/CME transport source shell\n";
    std::cout << message.str();
  }
}

double InitLoadMeasure(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  return node != nullptr && node->IsUsedInCalculationFlag ? 1.0 : 0.0;
}

void ApplyActiveRegionMask(
    const SEP3D::Mesh::ResolutionConfiguration& resolution) {
  using Node = cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>;
  gStaticLeafMaskInstalled = false;
  gStaticLeafMaskPruningApplied = false;
  gActiveRegionAllocationVerified = false;
  gPlannedActiveLeafCount = 0;
  gPlannedInactiveLeafCount = 0;
  gPlannedSolarInteriorLeafCount = 0;

  // BranchBottomNodeList and all neighbor links are replicated after
  // buildMesh(). Preserve that deterministic list order in both vectors so
  // the AMPS node and AMPS-independent physical box at index i remain paired.
  std::vector<Node*> nodes;
  std::vector<SEP3D::Mesh::LeafBlock> leaves;
  for (Node* node = PIC::Mesh::mesh->BranchBottomNodeList;
       node != nullptr; node = node->nextBranchBottomNode) {
    nodes.push_back(node);
    SEP3D::Mesh::LeafBlock leaf;
    leaf.minimumM = SEP3D::Core::Vec3(
        node->xmin[0], node->xmin[1], node->xmin[2]);
    leaf.maximumM = SEP3D::Core::Vec3(
        node->xmax[0], node->xmax[1], node->xmax[2]);
    leaf.level = static_cast<unsigned>(node->RefinmentLevel);
    leaf.globalLeaf = leaves.size();
    leaves.push_back(leaf);
  }
  if (nodes.empty()) {
    StopWithStatus("static AMPS leaf mask", SEP3D::Core::Status(
        SEP3D::Core::StatusCode::ConfigurationConflict,
        "AMPS finalized an empty leaf list before boundary/mask planning"));
  }

  // Construct the two graph views from AMPS' own coarse/fine-aware neighbor
  // API. A face has four sub-neighbor slots in 3-D, an edge has two, and a
  // corner has one. Deduplication is essential because same-level neighbors
  // legitimately appear in more than one sub-slot.
  std::unordered_map<Node*, std::size_t> index;
  index.reserve(nodes.size());
  for (std::size_t i = 0; i < nodes.size(); ++i) index[nodes[i]] = i;
  SEP3D::Mesh::LeafNeighbourGraph graph;
  graph.face.resize(nodes.size());
  graph.full.resize(nodes.size());
  auto addNeighbour = [&](std::size_t from, Node* neighbour,
                          bool isFace) {
    if (neighbour == nullptr) return;
    const auto found = index.find(neighbour);
    if (found == index.end()) {
      StopWithStatus("active Parker corridor", SEP3D::Core::Status(
          SEP3D::Core::StatusCode::LayoutMismatch,
          "AMPS neighbor API returned a node outside BranchBottomNodeList"));
    }
    const std::size_t to = found->second;
    if (to == from) return;
    graph.full[from].push_back(to);
    graph.full[to].push_back(from);
    if (isFace) {
      graph.face[from].push_back(to);
      graph.face[to].push_back(from);
    }
  };
  for (std::size_t from = 0; from < nodes.size(); ++from) {
    Node* node = nodes[from];
    for (int face = 0; face < 6; ++face)
      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
          addNeighbour(from,
              node->GetNeibFace(face, i, j, PIC::Mesh::mesh), true);
    for (int edge = 0; edge < 12; ++edge)
      for (int i = 0; i < 2; ++i)
        addNeighbour(from,
            node->GetNeibEdge(edge, i, PIC::Mesh::mesh), false);
    for (int corner = 0; corner < 8; ++corner)
      addNeighbour(from,
          node->GetNeibCorner(corner, PIC::Mesh::mesh), false);
  }
  for (auto* rows : {&graph.face, &graph.full}) {
    for (std::vector<std::size_t>& row : *rows) {
      std::sort(row.begin(), row.end());
      row.erase(std::unique(row.begin(), row.end()), row.end());
    }
  }

  SEP3D::Mesh::ActiveRegionPlan plan;
  const SEP3D::Core::Status planned = SEP3D::Mesh::BuildActiveRegionPlan(
      leaves, graph, resolution, &plan);
  if (!planned.ok()) StopWithStatus("active Parker corridor", planned);

  // AMPS keeps leaves that are wholly inside an internal sphere when its
  // build-time outside-domain policy is KEEP.  They contain no computational
  // plasma and must not consume block storage.  Compose that solid-body mask
  // with (rather than replace) the optional Parker-corridor mask.  In
  // parker-tube mode this never forces a disconnected photospheric island to
  // become active: the union can only deactivate an already selected leaf.
  const SEP3D::Mesh::SolarBoundaryGeometry solarBoundary =
      SEP3D::Mesh::MakeSolarBoundary(Configuration().options());
  std::vector<unsigned char> solarInterior(nodes.size(), 0);
  double activeBlockVolumeM3 = 0.0;
  for (std::size_t ordinal = 0; ordinal < nodes.size(); ++ordinal) {
    solarInterior[ordinal] =
        SEP3D::Mesh::AxisAlignedBoxEntirelyInsideSolarBoundary(
            leaves[ordinal].minimumM, leaves[ordinal].maximumM,
            solarBoundary)
        ? 1
        : 0;
    if (solarInterior[ordinal] != 0)
      ++gPlannedSolarInteriorLeafCount;
    const bool tubeActive =
        plan.leafClass[ordinal] != SEP3D::Mesh::ActiveLeafClass::Inactive;
    if (tubeActive && solarInterior[ordinal] == 0) {
      ++gPlannedActiveLeafCount;
      const SEP3D::Core::Vec3 side =
          leaves[ordinal].maximumM - leaves[ordinal].minimumM;
      activeBlockVolumeM3 += side.x * side.y * side.z;
    } else {
      ++gPlannedInactiveLeafCount;
    }
  }

  // Give each inactive leaf to exactly one rank by deterministic ordinal;
  // SetTreeNodeActiveUseFlag gathers those disjoint ID lists and broadcasts
  // the resulting changes to every replica of the AMR tree.
  std::list<Node*> inactive;
  for (std::size_t ordinal = 0; ordinal < nodes.size(); ++ordinal) {
    if (ordinal % static_cast<std::size_t>(PIC::nTotalThreads) !=
        static_cast<std::size_t>(PIC::ThisThread)) {
      continue;
    }
    if (plan.leafClass[ordinal] == SEP3D::Mesh::ActiveLeafClass::Inactive ||
        solarInterior[ordinal] != 0)
      inactive.push_back(nodes[ordinal]);
  }

  // Match the MPI datatype exactly.  std::uint64_t is not required to be an
  // alias of unsigned long long on every supported compiler, even when both
  // happen to be 64 bits; using the exact C++ type avoids an ABI-dependent
  // collective buffer mismatch.
  const unsigned long long localInactive =
      static_cast<unsigned long long>(inactive.size());
  unsigned long long globalInactive = 0;
  MPI_Allreduce(&localInactive, &globalInactive, 1, MPI_UNSIGNED_LONG_LONG,
                MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
  if (globalInactive !=
      static_cast<unsigned long long>(gPlannedInactiveLeafCount)) {
    StopWithStatus("static AMPS leaf mask", SEP3D::Core::Status(
        SEP3D::Core::StatusCode::LayoutMismatch,
        "MPI inactive-leaf partition disagrees with the replicated mask plan"));
  }
  if (globalInactive >= nodes.size()) {
    StopWithStatus("static AMPS leaf mask", SEP3D::Core::Status(
        SEP3D::Core::StatusCode::ConfigurationConflict,
        "solar/active-region mask would deactivate every AMR leaf block"));
  }

  // All ranks make the same branch decision because globalInactive is an
  // Allreduce result.  SetTreeNodeActiveUseFlag itself is collective, so it
  // must either be entered by every rank or by none.
  if (globalInactive != 0) {
    PIC::Mesh::mesh->SetTreeNodeActiveUseFlag(
        &inactive, nullptr, false, nullptr);
  }
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    const bool expected =
        plan.leafClass[i] != SEP3D::Mesh::ActiveLeafClass::Inactive &&
        solarInterior[i] == 0;
    if (nodes[i]->IsUsedInCalculationFlag != expected) {
      StopWithStatus("static AMPS leaf mask", SEP3D::Core::Status(
          SEP3D::Core::StatusCode::LayoutMismatch,
          "AMPS active-use flags disagree with the installed mask plan"));
    }
  }
  // Every leaf flag has matched the composed corridor/solar-interior plan.
  // Installation is complete even when full-domain mode removed no leaves.
  // The Parker-tube planner, halo/cavity rules and solar exclusion above are
  // unchanged: this flag records successful verification, not a new mask.
  gStaticLeafMaskPruningApplied = globalInactive != 0;
  gStaticLeafMaskInstalled = true;
  if (PIC::ThisThread == 0) {
    const double activeFraction = static_cast<double>(
        gPlannedActiveLeafCount) / static_cast<double>(nodes.size());
    const double volumeFraction = activeBlockVolumeM3 /
        plan.totalBlockVolumeM3;
    std::cout << "[srcSEP3D] static leaf mask algorithm="
              << SEP3D::Mesh::ActiveRegionAlgorithmName(resolution)
              << " core=" << plan.coreLeafCount
              << " halo=" << plan.haloLeafCount
              << " cavities_filled=" << plan.cavityLeafCount
              << " solar_interior=" << gPlannedSolarInteriorLeafCount
              << " active=" << gPlannedActiveLeafCount
              << " inactive=" << gPlannedInactiveLeafCount
              << " total=" << nodes.size()
              << " active_leaf_fraction=" << activeFraction
              << " active_volume_fraction=" << volumeFraction
              << " finite_segments=" << plan.segmentCount
              << " finite_length_m=" << plan.effectiveLineLengthM << '\n';
    if (resolution.activeRegion ==
            SEP3D::RuntimeModel::ActiveRegionMode::ParkerTube &&
        gPlannedInactiveLeafCount == 0) {
      std::cout << "[srcSEP3D] WARNING: parker-tube mode retained every leaf; "
                   "the selected width/halo provides no memory pruning\n";
    }
  }
}

void VerifyActiveRegionAllocation() {
  using Node = cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>;
  gActiveRegionAllocationVerified = false;
  if (!gStaticLeafMaskInstalled) {
    StopWithStatus("static leaf-mask allocation", SEP3D::Core::Status(
        SEP3D::Core::StatusCode::InvalidTransition,
        "block allocation was reached before active-region plan installation"));
  }
  unsigned long long localOwnedActive = 0;
  std::size_t replicatedActive = 0;
  std::size_t replicatedInactive = 0;
  for (Node* node = PIC::Mesh::mesh->BranchBottomNodeList;
       node != nullptr; node = node->nextBranchBottomNode) {
    if (!node->IsUsedInCalculationFlag) {
      ++replicatedInactive;
      if (node->block != nullptr) {
        StopWithStatus("static leaf-mask allocation", SEP3D::Core::Status(
            SEP3D::Core::StatusCode::LayoutMismatch,
            "inactive AMR leaf unexpectedly owns allocated block storage"));
      }
      continue;
    }
    ++replicatedActive;
    if (node->Thread == PIC::ThisThread) {
      if (node->block == nullptr) {
        StopWithStatus("static leaf-mask allocation", SEP3D::Core::Status(
            SEP3D::Core::StatusCode::LayoutMismatch,
            "owner-local active AMR leaf has no allocated block storage"));
      }
      ++localOwnedActive;
    }
  }
  if (replicatedActive != gPlannedActiveLeafCount ||
      replicatedInactive != gPlannedInactiveLeafCount) {
    StopWithStatus("static leaf-mask allocation", SEP3D::Core::Status(
        SEP3D::Core::StatusCode::LayoutMismatch,
        "replicated active/inactive counts changed after load distribution"));
  }
  unsigned long long globalOwnedActive = 0;
  MPI_Allreduce(&localOwnedActive, &globalOwnedActive, 1,
                MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
  if (globalOwnedActive !=
      static_cast<unsigned long long>(gPlannedActiveLeafCount)) {
    StopWithStatus("static leaf-mask allocation", SEP3D::Core::Status(
        SEP3D::Core::StatusCode::LayoutMismatch,
            "allocated owner-block count disagrees with the active-region plan"));
  }
  gActiveRegionAllocationVerified = true;
}

void CorrectSolarInteriorCellMeasures() {
  if (gSolarSurfaceBoundary == nullptr) {
    StopWithStatus("solar cell-measure correction",
        SEP3D::Core::Status(SEP3D::Core::StatusCode::InvalidTransition,
                           "the AMPS solar sphere is not registered"));
  }
  const SEP3D::Mesh::SolarBoundaryGeometry boundary =
      SEP3D::Mesh::MakeSolarBoundary(Configuration().options());
  if (gSolarSurfaceBoundary->Radius != boundary.radiusM ||
      gSolarSurfaceBoundary->OriginPosition[0] != boundary.centerM.x ||
      gSolarSurfaceBoundary->OriginPosition[1] != boundary.centerM.y ||
      gSolarSurfaceBoundary->OriginPosition[2] != boundary.centerM.z) {
    StopWithStatus("solar cell-measure correction",
        SEP3D::Core::Status(
            SEP3D::Core::StatusCode::LayoutMismatch,
            "registered AMPS sphere differs from the authoritative solar geometry"));
  }

  // AMPS' analytic spherical-volume routine correctly evaluates cut cells,
  // but its early _AMR_BLOCK_OUTSIDE_DOMAIN_ branch returns a full Cartesian
  // volume for a cell wholly inside the solid sphere.  Preserve every
  // fractional cut-cell value and repair only boxes proven wholly interior.
  // Ghost measures are repaired too so interpolation/movers cannot recover a
  // positive-volume interior cell through a neighboring block's halo.
  unsigned long long localCorrectedPhysicalCells = 0;
  for (unsigned int blockIndex = 0;
       blockIndex < PIC::DomainBlockDecomposition::nLocalBlocks;
       ++blockIndex) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
        PIC::DomainBlockDecomposition::BlockTable[blockIndex];
    if (node == nullptr || node->block == nullptr) continue;
    const double spacingM[3] = {
        (node->xmax[0] - node->xmin[0]) / _BLOCK_CELLS_X_,
        (node->xmax[1] - node->xmin[1]) / _BLOCK_CELLS_Y_,
        (node->xmax[2] - node->xmin[2]) / _BLOCK_CELLS_Z_};
    for (int k = -_GHOST_CELLS_Z_;
         k < _BLOCK_CELLS_Z_ + _GHOST_CELLS_Z_; ++k) {
      for (int j = -_GHOST_CELLS_Y_;
           j < _BLOCK_CELLS_Y_ + _GHOST_CELLS_Y_; ++j) {
        for (int i = -_GHOST_CELLS_X_;
             i < _BLOCK_CELLS_X_ + _GHOST_CELLS_X_; ++i) {
          PIC::Mesh::cDataCenterNode* cell = node->block->GetCenterNode(
              PIC::Mesh::mesh->getCenterNodeLocalNumber(i, j, k));
          if (cell == nullptr) continue;
          const SEP3D::Core::Vec3 cellMinimumM(
              node->xmin[0] + i * spacingM[0],
              node->xmin[1] + j * spacingM[1],
              node->xmin[2] + k * spacingM[2]);
          const SEP3D::Core::Vec3 cellMaximumM(
              cellMinimumM.x + spacingM[0],
              cellMinimumM.y + spacingM[1],
              cellMinimumM.z + spacingM[2]);
          if (!SEP3D::Mesh::AxisAlignedBoxEntirelyInsideSolarBoundary(
                  cellMinimumM, cellMaximumM, boundary)) {
            continue;
          }
          const bool physical =
              i >= 0 && i < _BLOCK_CELLS_X_ &&
              j >= 0 && j < _BLOCK_CELLS_Y_ &&
              k >= 0 && k < _BLOCK_CELLS_Z_;
          if (physical && cell->Measure != 0.0)
            ++localCorrectedPhysicalCells;
          cell->Measure = 0.0;
        }
      }
    }
  }
  unsigned long long globalCorrectedPhysicalCells = 0;
  MPI_Allreduce(&localCorrectedPhysicalCells,
                &globalCorrectedPhysicalCells, 1,
                MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                MPI_GLOBAL_COMMUNICATOR);
  if (PIC::ThisThread == 0) {
    std::cout << "[srcSEP3D] solar cell measures: corrected_fully_inside="
              << globalCorrectedPhysicalCells
              << " (fractional AMPS cut-cell measures preserved)\n";
  }
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
#if _PIC_COUPLER_MODE_ != _PIC_COUPLER_MODE__DATAFILE_
  // Reject an unsupported native coupler before allocating AMR storage; the
  // runtime bridge has an explicit one-fluid DATAFILE mapping, not a generic
  // assumption about buffers owned by another AMPS coupler.
  if (Configuration().options().background==SEP3D::RuntimeModel::BackgroundAuthority::Swcme ||
      Configuration().options().background==SEP3D::RuntimeModel::BackgroundAuthority::RuntimeModel)
    StopWithStatus("runtime mesh background",SEP3D::Core::Status(
        SEP3D::Core::StatusCode::ConfigurationConflict,"runtime mesh providers require the one-fluid AMPS DATAFILE buffer layout"));
#endif
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
  // The legacy automatic splitter is nonrelativistic and copies application
  // extension bytes verbatim.  Keep it disabled; the configured SEP-aware
  // controller runs explicitly at a joined boundary after the transport
  // ledger closes.
  PIC::ParticleSplitting::SetMode(PIC::ParticleSplitting::_disactivated);
  SEP3D::Init_BeforeParser();
  // The shared-file application parser belongs precisely at this boundary:
  // MPI, AMPS registries, and srcSEP3D's pre-parser hooks exist, while the
  // application species binding, sampling offsets, AMR tree, providers, and
  // particles do not. Parser errors therefore terminate before any partially
  // initialized model state can be mistaken for a valid run.
  ParseInstalledApplicationInput();
  // Shared-section parsing can replace an initially analytic provisional
  // configuration with the runtime reduced provider. Recheck the compile-time
  // native buffer ABI after that transaction; the all-rank helper keeps this
  // guard independent of the likewise all-rank solar-sphere registration.
  ValidateParsedRuntimeMeshBackgroundABI();
  // Resolve the selected reduced-front parameters and magnetic assets at the
  // first legal post-parser boundary. Construction is mesh-independent and
  // fails before AMR allocation; later epoch publication reuses this exact
  // provider rather than reparsing or creating a second model authority.
  InitializeSelectedModelsAfterParser();
  // Capture and validate the complete generated species table only after the
  // application parser has committed its immutable count. The table itself is
  // read-only: count, symbols, masses, charges, and indices were fixed by
  // SpeciesList and cannot be redefined by the runtime input file.
  BindCompiledSpeciesTable();
  // Internal surfaces must be registered after AMPS creates its global
  // registries and before mesh->init() creates the root tree.  The sphere is
  // the physical 1-R_sun photosphere, not the independently configurable
  // Parker/CME source shell at domain.inner_radius_m.
  RegisterSolarSurfaceBoundary();
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
  if (PIC::ThisThread == 0) {
    std::printf("[srcSEP3D] domain geometry=%s minimum_m=(%.9e,%.9e,%.9e) "
                "maximum_m=(%.9e,%.9e,%.9e) outer_radius_m=%.9e "
                "selected_endpoint_radius_m=%.9e "
                "active_solar_sphere_radius_m=%.9e solar_refinement_anchor=%s\n",
        SEP3D::RuntimeModel::Name(options.domainBoxGeometry),
        minimum[0], minimum[1], minimum[2], maximum[0], maximum[1], maximum[2],
        options.outerRadiusM, options.parkerSpiralEndRadiusM,
        options.activeSolarSphereRadiusM,
        SEP3D::RuntimeModel::Name(options.solarRefinementAnchor));
  }

  // Build and partition the AMR tree before block allocation.  This is the
  // same ordering used by mature AMPS applications and guarantees that each
  // MPI rank fills only the blocks it owns after decomposition.
  PIC::Mesh::mesh->AllowBlockAllocation = false;
  PIC::Mesh::mesh->init(minimum, maximum, localResolution);
  PIC::Mesh::mesh->buildMesh();
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  // AMPS' public active-use API must be called after the complete tree exists
  // and before load measurement/distribution/block allocation.  Inactive
  // leaves then consume neither center-cell storage nor particle lists, while
  // movers see them through the normal DomainExit path.
  ApplyActiveRegionMask(resolution);
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
  // Verify the static mask at the first point where allocation is observable.
  // This check is intentionally before any background, weight, or time-step
  // initialization so an inactive resident block cannot receive plausible
  // physics data and hide an activation-order regression.
  VerifyActiveRegionAllocation();

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

  // The global cell-crossing step depends on the minimum scale AMPS actually
  // allocated, not merely a nominal input resolution. Finalize that step and
  // the corresponding per-species statistical weights before Runtime binds
  // the mesh and before any particle/background state is installed.
  FinalizeParticleNumericsAfterMeshAllocation();

  PIC::Mesh::mesh->InitCellMeasure();
  CorrectSolarInteriorCellMeasures();
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
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__DATAFILE_
  // DATAFILE is the native center-node buffer ABI for this application, while
  // SEP3D's BackgroundProvider/Runtime owns the complete field snapshots. Claim
  // that ownership on EVERY rank before FillAndPublishBackground or any native
  // getter can run. Parker, SWCME and registered runtime models all use this
  // same publication path; hard-coding a SWCME-only coupler would duplicate it.
  //
  // Do not call MULTIFILE::Init or invent a Schedule entry: that file loader
  // resets the simulation clock and can replace the provider's validated data.
  // The core ownership policy disables file scheduling, EOF termination and
  // file time interpolation, without changing any allocated offsets. Ordinary
  // DATAFILE applications retain their default FileSchedule policy. The policy
  // remains fixed during transport; RefreshBackgroundAtBoundary publishes the
  // next complete snapshot and exchanges ghosts at the configured cadence.
  PIC::CPLR::DATAFILE::BackgroundUpdatePolicy =
      PIC::CPLR::DATAFILE::BackgroundUpdateMode::RuntimeProvider;
  if (PIC::ThisThread == 0)
    std::cout << "[srcSEP3D] native background updates: runtime-provider"
              << " (DATAFILE storage; file scheduling/interpolation disabled)\n";
#endif
  // BindCompiledSpeciesTable() already captured and validated the generated
  // AMPS table in amps_init_mesh().  No molecular-data setter is called: AMPS
  // remains the sole authority for the immutable species identity and physics.
  // R04: the mesh-derived step is global for every compiled species. The base
  // statistical weight is species-specific because upstream electron, proton
  // and alpha incident rates need not be equal. Source sampling may later add
  // an individual correction for exact patch allocation; that does not change
  // these global/block base values.
  const double configuredDt = Configuration().options().requestedTimeStepS;
  for (const auto& species : gCompiledSpecies) {
    PIC::ParticleWeightTimeStep::GlobalTimeStep[species.ampsIndex] =
        configuredDt;
    PIC::ParticleWeightTimeStep::GlobalParticleWeight[species.ampsIndex] =
        ConfiguredParticleWeight(species.ampsIndex);
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
          ConfiguredParticleWeight(species.ampsIndex), species.ampsIndex);
    }
  }
  FillAndPublishBackground();
  SEP3D::AMPS::Movers::Context mover;
  mover.resolveLocal = ResolveLocalTransport;
  mover.resolveMagneticDirection = ResolvePopulationMagneticDirection;
  mover.ledger = &gParticleLedger;
  mover.maximumSubsteps = Configuration().options().maximumTransportSubsteps;
  if (gInstalledShock) {
    const SEP3D::Adapters::ShockState shock =
        gInstalledShock->Evaluate(SEP3D::ApplicationRuntime().CurrentTimeS());
    if (!shock.status.ok()) StopWithStatus("initial shock state", shock.status);
    // Preserve finite shape/axis/width together with the epoch's apex state.
    // Rebuilding this record from radius alone would silently restore a sphere.
    mover.shock = shock.MoverGeometry();
  } else if (Configuration().options().source.enabled) {
    const auto reduced=std::dynamic_pointer_cast<
        SEP3D::Adapters::ShockFrontBackgroundAdapter>(
            gRuntimeBackgroundProvider);
    const auto epoch=reduced?reduced->FrontEpoch():nullptr;
    const auto provider=reduced?reduced->SharedProvider():nullptr;
    if(!epoch||!provider)
      StopWithStatus("source initialization", SEP3D::Core::Status(
          SEP3D::Core::StatusCode::SnapshotUnavailable,
          "enabled reduced-front source has no committed surface epoch"));
    mover.shock=ReducedMoverShock(*epoch,provider->Event());
  }
  const SEP3D::Core::Status installed =
      SEP3D::AMPS::Movers::InstallContext(mover);
  if (!installed.ok()) StopWithStatus("AMPS mover context installation", installed);
  if(gHasParsedApplicationInput&&Configuration().options().source.enabled) {
    if(Configuration().options().source.weightingModel!=
        SEP3D::RuntimeModel::SourceWeightingModel::ConstantStatisticalWeight)
      StopWithStatus("reduced-front source selection",
          SEP3D::Core::Status::Reserved(
              "log-uniform momentum importance weighting for reduced-front injection"));
    if(PIC::BC::UserDefinedParticleInjectionFunction!=nullptr&&
        PIC::BC::UserDefinedParticleInjectionFunction!=
            InjectReducedShockSurfaceParticles)
      StopWithStatus("reduced-front source callback",SEP3D::Core::Status(
          SEP3D::Core::StatusCode::ConfigurationConflict,
          "another application particle-injection callback is already installed"));
    PIC::BC::UserDefinedParticleInjectionFunction=
        InjectReducedShockSurfaceParticles;
    if(PIC::ThisThread==0)
      std::cout<<"[srcSEP3D] installed reduced-front constant-weight "
          "particle source through PIC::BC::UserDefinedParticleInjectionFunction\n";
  }
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

SEP3D::Core::Status SEP3D::Validation::CaptureNativeApplicationState(
    int expectedMpiRanks, NativeApplicationState* state) {
  using SEP3D::Core::Status;
  using SEP3D::Core::StatusCode;

  if (state == nullptr) {
    return Status(StatusCode::InvalidInput,
                  "native AMPS application-state output is null");
  }

  // This routine is deliberately a read-only, collective observation made at
  // a joined application boundary.  It does not call a provider's Prepare(),
  // rebuild the mesh, alter particles, or advance an RNG.  Consequently the
  // evidence describes the same state that a production time step consumes.
  NativeApplicationState captured;
  MPI_Comm_size(MPI_GLOBAL_COMMUNICATOR, &captured.mpiRankCount);
  captured.expectedMpiRanks = expectedMpiRanks;
  captured.configurationFingerprint = Configuration().physics_fingerprint();
  captured.inputSchemaVersion = Configuration().options().inputSchemaVersion;
  captured.backgroundAuthority = RuntimeModel::Name(Configuration().options().background);
  captured.shockAuthority = RuntimeModel::Name(Configuration().options().shock);
  captured.plannedActiveLeaves = gPlannedActiveLeafCount;
  captured.plannedInactiveLeaves = gPlannedInactiveLeafCount;
  captured.plannedSolarInteriorLeaves = gPlannedSolarInteriorLeafCount;
  captured.activeRegionMode = RuntimeModel::Name(
      Configuration().options().activeRegion);
  captured.solarBoundaryRegistered = gSolarSurfaceBoundary != nullptr;
  captured.activeMaskInstalled = gStaticLeafMaskInstalled;
  captured.activeRegionPruningApplied = gStaticLeafMaskPruningApplied;
  captured.activeRegionAllocationVerified = gActiveRegionAllocationVerified;
  captured.backgroundReady = gInstalledBackground != nullptr;
  captured.turbulenceReady = gInstalledTurbulence != nullptr &&
      gInstalledTurbulence->PreparedMetadata() != nullptr;
  captured.sourceEnabled = Configuration().options().source.enabled;
  captured.shockRequired =
      Configuration().options().shock != RuntimeModel::ShockAuthority::None ||
      captured.sourceEnabled;
  captured.restartConfigured =
      !Configuration().options().restartInputPath.empty();
  captured.checkpointSequence=
      ApplicationRuntime().counters().checkpointSequence;
  if(gPendingRestart) {
    captured.restartInputTick=gPendingRestart->runtimeCounters.currentTick;
    captured.restartInputBackgroundGeneration=
        gPendingRestart->backgroundGeneration;
    captured.restartSourceRankCount=gPendingRestart->savedRankCount;
  }
  captured.completedSteps = ApplicationRuntime().counters().completedSteps;

  const unsigned long long localBlockCount =
      static_cast<unsigned long long>(
          PIC::DomainBlockDecomposition::nLocalBlocks);
  unsigned long long globalBlockCount = 0;
  MPI_Allreduce(&localBlockCount, &globalBlockCount, 1,
                MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
  captured.globalAllocatedBlocks = globalBlockCount;

  unsigned long long localPhysicalCellCount = 0;
  bool localBackgroundFinite = captured.backgroundReady;
  bool localDerivativeFinite = captured.backgroundReady;
  bool localTurbulenceFinite = captured.turbulenceReady;
  if (gInstalledBackground) {
    localPhysicalCellCount = static_cast<unsigned long long>(
        gInstalledBackground->samples().size());
    const auto& positions = gInstalledBackground->positions();
    const auto& samples = gInstalledBackground->samples();
    const auto& capabilities = gInstalledBackground->capabilities();
    if (positions.size() != samples.size()) {
      localBackgroundFinite = false;
      localDerivativeFinite = false;
      localTurbulenceFinite = false;
    }
    for (std::size_t i = 0; i < samples.size(); ++i) {
      const Background::BackgroundSample& sample = samples[i];
      localBackgroundFinite = localBackgroundFinite &&
          Background::ValidateCompleteSample(sample, capabilities).ok();

      // Validate the full derivative record, including fields not advertised
      // as analytic.  Numerical AMR reconstruction may populate those fields;
      // a non-finite value is never acceptable to a focused-transport mover.
      bool derivativesFinite =
          std::isfinite(sample.divBhat) &&
          std::isfinite(sample.focusingLenM) &&
          std::isfinite(sample.curvature.x) &&
          std::isfinite(sample.curvature.y) &&
          std::isfinite(sample.curvature.z) &&
          std::isfinite(sample.divU) &&
          std::isfinite(sample.fieldAlignedStrain);
      for (int row = 0; row < 3; ++row) {
        for (int column = 0; column < 3; ++column) {
          derivativesFinite = derivativesFinite &&
              std::isfinite(sample.gradB(row, column)) &&
              std::isfinite(sample.gradU(row, column));
        }
      }
      localDerivativeFinite = localDerivativeFinite && derivativesFinite;

      if (gInstalledTurbulence && i < positions.size()) {
        const Turbulence::TurbulenceSample waves =
            gInstalledTurbulence->Evaluate(positions[i], sample);
        // Ballistic is a valid explicit missing-data policy.  In that case all
        // numerical fields must still be finite zeros; NaN is never used as a
        // sentinel in either production output or native test evidence.
        const double values[] = {
            waves.deltaB2T2, waves.deltaBPlus2T2,
            waves.deltaBMinus2T2, waves.deltaBOutward2T2,
            waves.deltaBInward2T2, waves.waveEnergyPlusJPerM3,
            waves.waveEnergyMinusJPerM3, waves.kMinPerM, waves.kMaxPerM,
            waves.spectralIndex, waves.parallelCorrelationLengthM};
        bool finite = waves.status.usable() &&
            (waves.valid || waves.ballistic);
        for (double value : values) finite = finite && std::isfinite(value);
        if (!waves.ballistic) {
          finite = finite && waves.valid && waves.deltaB2T2 > 0.0 &&
              waves.kMinPerM > 0.0 &&
              waves.kMaxPerM > waves.kMinPerM &&
              waves.parallelCorrelationLengthM > 0.0 &&
              waves.generation > 0;
        }
        localTurbulenceFinite = localTurbulenceFinite && finite;
      } else {
        localTurbulenceFinite = false;
      }
    }
  }
  unsigned long long globalPhysicalCellCount = 0;
  MPI_Allreduce(&localPhysicalCellCount, &globalPhysicalCellCount, 1,
                MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_GLOBAL_COMMUNICATOR);
  captured.globalPhysicalCells = globalPhysicalCellCount;

  auto CollectiveAnd = [](bool local) {
    int input = local ? 1 : 0;
    int output = 0;
    MPI_Allreduce(&input, &output, 1, MPI_INT, MPI_MIN,
                  MPI_GLOBAL_COMMUNICATOR);
    return output != 0;
  };
  if (Configuration().options().background==RuntimeModel::BackgroundAuthority::Swcme ||
      Configuration().options().background==RuntimeModel::BackgroundAuthority::RuntimeModel) {
    // Capture at the actual production boundary. The denominator below counts
    // due cadence events from tick zero; gBackgroundPublishedUpdates counts
    // completed commits independently, so a missed refresh becomes a failure.
    bool owned=false,ghosts=false,provider=false;
    unsigned long long localGhosts=0,globalGhosts=0;
    unsigned long long localOwnerXor=0,globalOwnerXor=0;
    unsigned long long localOwnerSum=0,globalOwnerSum=0;
    CaptureRuntimeMeshBackground(&owned,&ghosts,&provider,&localGhosts,
        &localOwnerXor,&localOwnerSum);
    captured.runtimeMeshOwnedFieldsMatch=CollectiveAnd(owned);
    captured.runtimeMeshGhostFieldsMatch=CollectiveAnd(ghosts);
    captured.runtimeMeshProviderMatch=CollectiveAnd(provider);
    MPI_Allreduce(&localGhosts,&globalGhosts,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(&localOwnerXor,&globalOwnerXor,1,
        MPI_UNSIGNED_LONG_LONG,MPI_BXOR,MPI_GLOBAL_COMMUNICATOR);
    MPI_Allreduce(&localOwnerSum,&globalOwnerSum,1,
        MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
    captured.runtimeMeshGhostCellsChecked=globalGhosts;
    captured.runtimeMeshOwnerFingerprintXor=globalOwnerXor;
    captured.runtimeMeshOwnerFingerprintSum=globalOwnerSum;
    captured.runtimeMeshPublishedUpdates=gBackgroundPublishedUpdates;
    captured.runtimeMeshExpectedUpdates=ApplicationRuntime().counters().currentTick/
        Configuration().options().backgroundCadenceSteps;

    // Exercise the exact pre-write collective gate with an injected failure
    // on one rank.  The candidate is the already committed immutable epoch;
    // no Prepare/Store call occurs.  A correct gate makes every rank reject
    // while the installed pointer/generation and native bytes remain intact.
    // This supplies RSH24's rank-local rollback subgate without a test-only
    // branch in production stepping.
    const auto committedBefore=gInstalledBackground;
    const auto generationBefore=committedBefore?committedBefore->metadata().generation:0;
    std::vector<Turbulence::TurbulenceSample> injectedScratch;
    Core::Status injectedStatus=PIC::ThisThread==0?
        Core::Status(Core::StatusCode::BackgroundInvalid,
          "intentional native rank-local candidate rejection"):
        Core::Status::OK();
    const auto rejected=ValidateBackgroundCandidateCollectively(
        gInstalledBackground,CollectOwnedPhysicalCells(),gInstalledTurbulence,
        injectedStatus,&injectedScratch);
    captured.runtimeCollectiveRollbackVerified=!rejected.ok()&&
        gInstalledBackground==committedBefore&&gInstalledBackground&&
        gInstalledBackground->metadata().generation==generationBefore;
  }
  captured.finiteBackgroundAndTurbulence =
      CollectiveAnd(localBackgroundFinite && localTurbulenceFinite);
  captured.finiteBackgroundDerivatives =
      CollectiveAnd(localDerivativeFinite);
  captured.backgroundReady = CollectiveAnd(captured.backgroundReady);
  captured.turbulenceReady = CollectiveAnd(captured.turbulenceReady);
  captured.solarBoundaryRegistered =
      CollectiveAnd(captured.solarBoundaryRegistered);
  captured.activeMaskInstalled =
      CollectiveAnd(captured.activeMaskInstalled);
  captured.activeRegionPruningApplied =
      CollectiveAnd(captured.activeRegionPruningApplied);
  captured.activeRegionAllocationVerified =
      CollectiveAnd(captured.activeRegionAllocationVerified);

  bool localSpeciesNumericsReady =
      PIC::ParticleWeightTimeStep::GlobalTimeStepInitialized &&
      gCompiledSpecies.size() == static_cast<std::size_t>(PIC::nTotalSpecies);
  captured.species.reserve(gCompiledSpecies.size());
  for (const auto& species : gCompiledSpecies) {
    NativeSpeciesState item;
    item.compiledSlot = species.ampsIndex;
    item.chemicalSymbol = species.symbol;
    item.massKg = species.massKg;
    item.chargeC = species.chargeC;
    item.timeStepS =
        PIC::ParticleWeightTimeStep::GlobalTimeStep[species.ampsIndex];
    item.particleWeight =
        PIC::ParticleWeightTimeStep::GlobalParticleWeight[species.ampsIndex];
    localSpeciesNumericsReady = localSpeciesNumericsReady &&
        item.compiledSlot >= 0 && !item.chemicalSymbol.empty() &&
        std::isfinite(item.massKg) && item.massKg > 0.0 &&
        std::isfinite(item.chargeC) && item.chargeC != 0.0 &&
        std::isfinite(item.timeStepS) && item.timeStepS > 0.0 &&
        std::isfinite(item.particleWeight) && item.particleWeight > 0.0;
    captured.species.push_back(item);
  }
  for (unsigned int blockIndex = 0;
       blockIndex < PIC::DomainBlockDecomposition::nLocalBlocks;
       ++blockIndex) {
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
        PIC::DomainBlockDecomposition::BlockTable[blockIndex];
    if (node == nullptr || node->block == nullptr) continue;
    for (const auto& species : gCompiledSpecies) {
      const double localStep =
          node->block->GetLocalTimeStep(species.ampsIndex);
      const double localWeight =
          node->block->GetLocalParticleWeight(species.ampsIndex);
      localSpeciesNumericsReady = localSpeciesNumericsReady &&
          std::isfinite(localStep) && localStep > 0.0 &&
          std::isfinite(localWeight) && localWeight > 0.0 &&
          localStep ==
              PIC::ParticleWeightTimeStep::GlobalTimeStep[species.ampsIndex] &&
          localWeight == PIC::ParticleWeightTimeStep::GlobalParticleWeight[
                             species.ampsIndex];
    }
  }
  const bool speciesNumericsReady = CollectiveAnd(localSpeciesNumericsReady);

  captured.backgroundGeneration = gInstalledBackground
      ? gInstalledBackground->metadata().generation : 0;
  bool localShockReady = !captured.shockRequired;
  if (gInstalledShock) {
    const Adapters::ShockState shock =
        gInstalledShock->Evaluate(ApplicationRuntime().CurrentTimeS());
    // Provider readiness is distinct from physical shock activation.  Before
    // event.valid_from the canonical SWCME provider deliberately publishes a
    // valid inactive state with generation=0: there is no shock to generate
    // yet.  That state still closes the Shock initialization stage.  Require
    // identity, time coverage and finite fields in either branch so missing
    // or stale provider data cannot pass merely because active=false.
    localShockReady =
        shock.status.ok() &&
        shock.Covers(ApplicationRuntime().CurrentTimeS()) &&
        !shock.providerIdentity.empty() &&
        !shock.configurationFingerprint.empty() &&
        std::isfinite(shock.epochS) &&
        std::isfinite(shock.validUntilS) &&
        std::isfinite(shock.centerM.x) &&
        std::isfinite(shock.centerM.y) &&
        std::isfinite(shock.centerM.z) &&
        std::isfinite(shock.radiusM) &&
        std::isfinite(shock.radialSpeedMPerS) &&
        (!shock.active ||
         (shock.generation > 0 && shock.radiusM > 0.0));
  }
  captured.shockReady = CollectiveAnd(localShockReady);
  // The application installs/evaluates the shock only after the current
  // background has been published.  Capturing both at this joined boundary is
  // therefore the authoritative generation pairing, including the explicit
  // no-shock transport-only case required by the common initialization ledger.
  captured.shockBackgroundGeneration = captured.shockReady
      ? captured.backgroundGeneration : 0;

  // Count actual AMPS particles after initialization/stepping.  A disabled
  // source and particles_per_cell=0 are configuration claims; the global
  // linked-list count is the independent native evidence that no restart,
  // baseline seed path, population control or source produced a particle.
  const std::vector<std::uint64_t> localParticleCounts=
      CountLocalParticlesBySpecies();
  unsigned long long localParticles=0,globalParticles=0;
  for(std::uint64_t count:localParticleCounts)localParticles+=count;
  MPI_Allreduce(&localParticles,&globalParticles,1,MPI_UNSIGNED_LONG_LONG,
      MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  unsigned long long localInjected=0,globalInjected=0;
  for(const auto& row:gSourceLedger)localInjected+=row.macroparticles;
  MPI_Allreduce(&localInjected,&globalInjected,1,MPI_UNSIGNED_LONG_LONG,
      MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  captured.globalParticleCount=globalParticles;
  captured.globalInjectedParticleCount=globalInjected;
  captured.zeroParticleAllocationRequested=
      Configuration().options().memoryModel.particlesPerCell==0.0;

  if(Configuration().options().background==
      RuntimeModel::BackgroundAuthority::RuntimeModel) {
    const auto reduced=std::dynamic_pointer_cast<
        Adapters::ShockFrontBackgroundAdapter>(gRuntimeBackgroundProvider);
    const auto epoch=reduced?reduced->FrontEpoch():nullptr;
    const auto shared=reduced?reduced->SharedProvider():nullptr;
    captured.reducedProviderSelected=static_cast<bool>(reduced&&epoch&&shared);
    if(captured.reducedProviderSelected) {
      captured.reducedFrontGeneration=epoch->generation;
      captured.reducedAmbientGeneration=epoch->ambientGeneration;
      captured.reducedEpochS=epoch->trajectory.timeS;
      captured.reducedApexRadiusM=epoch->trajectory.apexRadiusM;
      captured.reducedApexSpeedMPerS=epoch->trajectory.apexSpeedMPerS;
      captured.reducedPhase=SEP::CoronaSwcme::ShockFront::Name(
          epoch->trajectory.phase);
      captured.reducedEventIdentity=epoch->eventIdentity;
      // The front epoch is replicated, but equality of its apex alone is not
      // enough for restart equivalence.  Include every stable surface label,
      // position, normal, speed, area and physical/numerical classification.
      // This remains a diagnostic fingerprint: the independent RH and
      // geometry tests establish correctness of the encoded values.
      std::ostringstream frontState;
      frontState << std::hexfloat << epoch->trajectory.apexRadiusM << '|'
          << epoch->trajectory.apexSpeedMPerS << '|'
          << static_cast<int>(epoch->trajectory.phase);
      // Surface topology is part of native restart/rank equivalence.  Hash
      // every shared vertex and facet before its face-centred state so a pole
      // duplication, changed diagonal, winding reversal, or curved/chord-area
      // substitution cannot hide behind identical physical sample values.
      for(const auto& vertex:epoch->vertices)frontState<<"|v:"
          <<vertex.stableId<<':'<<vertex.positionM.x<<':'<<vertex.positionM.y
          <<':'<<vertex.positionM.z<<':'<<vertex.supportEdge<<':'<<vertex.apex;
      for(const auto& triangle:epoch->triangles)frontState<<"|f:"
          <<triangle.stableId<<':'<<triangle.vertex[0]<<':'<<triangle.vertex[1]
          <<':'<<triangle.vertex[2]<<':'<<triangle.curvedAreaM2<<':'
          <<triangle.planarAreaM2;
      for(const auto& record:epoch->records) {
        frontState << '|' << record.geometry.stableId << ':'
            << record.geometry.positionM.x << ':'
            << record.geometry.positionM.y << ':'
            << record.geometry.positionM.z << ':'
            << record.geometry.outwardNormal.x << ':'
            << record.geometry.outwardNormal.y << ':'
            << record.geometry.outwardNormal.z << ':'
            << record.geometry.normalSpeedMPerS << ':'
            << record.geometry.areaM2 << ':'
            << static_cast<int>(record.status) << ':'
            << record.inflowMPerS << ':' << record.fastMach << ':'
            << record.signedMagneticNormalCosine << ':'
            << record.downstreamValid;
      }
      captured.reducedFrontStateFingerprint=Fnv1a64(frontState.str());
      captured.reducedAcceptedAreaM2=epoch->area.acceptedShockM2;
      captured.reducedNumericalFailureAreaM2=epoch->area.numericalFailureM2;
      captured.reducedGeometricEndpointReached=
          epoch->geometricEndpointReached;
      captured.reducedApexShockAccepted=epoch->apexShockAccepted;
      const auto endpoint=shared->EndpointTimeS();
      if(endpoint.ok()) {
        captured.reducedEndpointTimeS=endpoint.value;
        const auto record=shared->EvaluateFrontPoint(
            shared->Event().observerPositionM,endpoint.value,UINT64_C(1));
        if(record.ok()) {
          captured.reducedEndpointObserverGeometricHit=true;
          captured.reducedEndpointObserverStatus=
              SEP::CoronaSwcme::ShockFront::Name(record.value.status);
          captured.reducedEndpointObserverShockAccepted=
              record.value.status==SEP::CoronaSwcme::ShockFront::FrontStatus::SolvedFastShock;
        }
      }
    }
  }

  std::vector<Output::VirtualSpacecraftDefinition> observers;
  const Status observerStatus = Output::BuildObserverDefinitions(
      Configuration(), ApplicationRuntime().CurrentTimeS(), &observers);
  const bool observersReady = CollectiveAnd(observerStatus.ok());
  const bool outputDictionaryReady = CollectiveAnd(
      gStaticCellDataOffset >= 0 &&
      (Configuration().storage_layout().samplingBytesPerCell == 0 ||
       gSamplingDataOffset >= 0) &&
      gStorageCallbacksRegistered);
  const bool haloReady = CollectiveAnd(gNativeAmpsBackgroundReady);

  if (captured.globalAllocatedBlocks > 0) captured.initializationMask |= 1U << 0;
  if (captured.solarBoundaryRegistered) captured.initializationMask |= 1U << 1;
  if (captured.backgroundReady) captured.initializationMask |= 1U << 2;
  if (captured.turbulenceReady) captured.initializationMask |= 1U << 3;
  if (captured.shockReady) captured.initializationMask |= 1U << 4;
  if (haloReady) captured.initializationMask |= 1U << 5;
  if (!captured.species.empty()) captured.initializationMask |= 1U << 6;
  if (speciesNumericsReady) captured.initializationMask |= 1U << 7;
  if (observersReady) captured.initializationMask |= 1U << 8;
  if (outputDictionaryReady) captured.initializationMask |= 1U << 9;

  int productsExist = 1;
  int productsFinite = 1;
  if (PIC::ThisThread == 0) {
    std::vector<std::string> paths = {
        Configuration().options().initializationMeshTecplotFile,
        Configuration().options().initializationParkerLineTecplotFile};
    for (const auto& species : gCompiledSpecies) {
      paths.push_back(InitializationDataPath(
          Configuration().options().initializationDataTecplotFile,
          species.ampsIndex, static_cast<int>(gCompiledSpecies.size())));
    }
    for (const std::string& path : paths) {
      std::error_code error;
      const bool exists = fs::is_regular_file(path, error) && !error &&
          fs::file_size(path, error) > 0 && !error;
      productsExist = productsExist && exists;
      productsFinite = productsFinite && exists &&
          FileHasOnlyFiniteNumericTokens(path);
    }
  }
  MPI_Bcast(&productsExist, 1, MPI_INT, 0, MPI_GLOBAL_COMMUNICATOR);
  MPI_Bcast(&productsFinite, 1, MPI_INT, 0, MPI_GLOBAL_COMMUNICATOR);
  captured.initializationProductsExist = productsExist != 0;
  captured.initializationProductsFinite = productsFinite != 0;

  // Hash only collective or immutable values.  Rank-local ownership is
  // represented by the globally reduced counts, so a valid decomposition can
  // differ between ranks without creating a false mismatch.
  std::ostringstream identity;
  identity << captured.configurationFingerprint << '|'
           << captured.inputSchemaVersion << '|'
           << captured.backgroundAuthority << '|'
           << captured.shockAuthority << '|'
           << captured.runtimeMeshPublishedUpdates << '|'
           << captured.runtimeMeshExpectedUpdates << '|'
           << captured.runtimeMeshOwnerFingerprintXor << '|'
           << captured.runtimeMeshOwnerFingerprintSum << '|'
           << captured.runtimeMeshOwnedFieldsMatch << '|'
           << captured.runtimeMeshGhostFieldsMatch << '|'
           << captured.runtimeCollectiveRollbackVerified << '|'
           << captured.backgroundGeneration << '|' 
           << captured.globalAllocatedBlocks << '|'
           << captured.globalPhysicalCells << '|'
           << captured.plannedActiveLeaves << '|'
           << captured.plannedInactiveLeaves << '|'
           << captured.plannedSolarInteriorLeaves << '|'
           << captured.activeRegionMode << '|'
           << captured.activeMaskInstalled << '|'
           << captured.activeRegionPruningApplied << '|'
           << captured.activeRegionAllocationVerified << '|'
           << captured.initializationMask << '|'
           << captured.globalParticleCount << '|'
           << captured.globalInjectedParticleCount << '|'
           << captured.reducedProviderSelected << '|'
           << captured.reducedFrontGeneration << '|'
           << captured.reducedAmbientGeneration << '|'
           << captured.reducedFrontStateFingerprint << '|'
           << captured.reducedEpochS << '|'
           << captured.reducedApexRadiusM << '|'
           << captured.reducedApexSpeedMPerS << '|'
           << captured.reducedEventIdentity << '|'
           << captured.reducedPhase;
  identity << std::setprecision(17);
  for (const NativeSpeciesState& species : captured.species) {
    identity << '|' << species.compiledSlot << ':' << species.chemicalSymbol
             << ':' << species.massKg << ':' << species.chargeC
             << ':' << species.timeStepS << ':' << species.particleWeight;
  }
  const unsigned long long localFingerprint =
      static_cast<unsigned long long>(Fnv1a64(identity.str()));
  unsigned long long minimumFingerprint = 0;
  unsigned long long maximumFingerprint = 0;
  MPI_Allreduce(&localFingerprint, &minimumFingerprint, 1,
                MPI_UNSIGNED_LONG_LONG, MPI_MIN, MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&localFingerprint, &maximumFingerprint, 1,
                MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_GLOBAL_COMMUNICATOR);
  captured.mpiFingerprintConsistent =
      minimumFingerprint == maximumFingerprint;

  *state = std::move(captured);
  return Status::OK();
}

SEP3D::Core::Status SEP3D::CaptureNativeShockHistorySample(
    Output::ShockHistorySample* sample) {
  using Core::Status;
  using Core::StatusCode;
  const auto& runtime=ApplicationRuntime();
  const auto& options=Configuration().options();
  Adapters::ShockState shock;
  if (gInstalledShock) shock=gInstalledShock->Evaluate(runtime.CurrentTimeS());
  const auto reduced=std::dynamic_pointer_cast<
      Adapters::ShockFrontBackgroundAdapter>(gRuntimeBackgroundProvider);
  const auto reducedEpoch=reduced?reduced->FrontEpoch():nullptr;
  const auto reducedProvider=reduced?reduced->SharedProvider():nullptr;

  // Normalize the two supported propagation authorities into the compact
  // history record without manufacturing a legacy ShockState for the reduced
  // provider.  In reduced mode `active` means the apex patch is an accepted
  // fast shock; radius still records the geometric apex when that flag is
  // false.  This preserves the model's essential geometric/physical-arrival
  // distinction in every native row.
  const bool reducedMode=static_cast<bool>(reduced&&reducedEpoch&&reducedProvider);
  const bool legacyMode=static_cast<bool>(gInstalledShock);
  const bool geometryOK=reducedMode?
      (std::isfinite(reducedEpoch->trajectory.apexRadiusM)&&
       reducedEpoch->trajectory.apexRadiusM>0&&
       std::isfinite(reducedEpoch->trajectory.apexSpeedMPerS)&&
       reducedEpoch->trajectory.apexSpeedMPerS>0):
      (legacyMode&&shock.status.ok()&&shock.Covers(runtime.CurrentTimeS())&&
       Adapters::ValidateShockGeometry(shock.MoverGeometry()).ok());
  int localOK=sample!=nullptr && (reducedMode!=legacyMode) && geometryOK &&
      options.intent==RuntimeModel::RunIntent::ShockPropagation &&
      !options.source.enabled;
  int allOK=0;
  MPI_Allreduce(&localOK,&allOK,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  if (!allOK) return Status(StatusCode::InvalidInput,"a rank could not capture installed propagation state");

  const std::vector<std::uint64_t> owned=CountLocalParticlesBySpecies();
  unsigned long long localCounts[2]={0,0}, globalCounts[2]={0,0};
  for (std::uint64_t count:owned) localCounts[0]+=count;
  // Source ledger rows are retained on the rank that actually allocated the
  // patch. This sum includes particles that have already escaped or been lost.
  for (const auto& row:gSourceLedger) localCounts[1]+=row.macroparticles;
  MPI_Allreduce(localCounts,globalCounts,2,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  const double radiusM=reducedMode?reducedEpoch->trajectory.apexRadiusM:shock.radiusM;
  const double speedMPerS=reducedMode?reducedEpoch->trajectory.apexSpeedMPerS:
      shock.radialSpeedMPerS;
  const Core::Vec3 direction=reducedMode?Core::Vec3(
      reducedProvider->Event().direction.x,reducedProvider->Event().direction.y,
      reducedProvider->Event().direction.z):shock.cmeDirection;
  const double halfWidthRad=reducedMode?reducedProvider->Event().halfWidthRad:
      shock.halfWidthRad;
  const std::uint64_t generation=reducedMode?reducedEpoch->generation:shock.generation;
  const bool active=reducedMode?reducedEpoch->apexShockAccepted:shock.active;
  const std::string providerIdentity=reducedMode?reduced->CanonicalName():
      shock.providerIdentity;
  const std::string providerFingerprint=reducedMode?reducedEpoch->eventIdentity:
      shock.configurationFingerprint;
  double values[7]={runtime.CurrentTimeS(),radiusM,speedMPerS,
      direction.x,direction.y,direction.z,halfWidthRad}, minimum[7],maximum[7];
  MPI_Allreduce(values,minimum,7,MPI_DOUBLE,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(values,maximum,7,MPI_DOUBLE,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  unsigned long long identity[5]={runtime.counters().currentTick,generation,
      static_cast<unsigned long long>(active),Fnv1a64(providerIdentity+"|"+
      providerFingerprint+"|"+Configuration().physics_fingerprint()),
      reducedMode?UINT64_C(2):static_cast<unsigned long long>(shock.geometry)};
  unsigned long long lo[5],hi[5];
  MPI_Allreduce(identity,lo,5,MPI_UNSIGNED_LONG_LONG,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(identity,hi,5,MPI_UNSIGNED_LONG_LONG,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  for (int i=0;i<5;++i) if (lo[i]!=hi[i])
    return Status(StatusCode::ConfigurationConflict,"MPI ranks disagree on propagation clock/generation/active state/provider identity");
  for(int i=3;i<7;++i) if(maximum[i]!=minimum[i])
    return Status(StatusCode::ConfigurationConflict,"MPI ranks disagree on finite shock axis/half width");
  if (maximum[2]-minimum[2]>1e-6 ||
      std::fabs(PIC::SimulationTime::Get()-runtime.CurrentTimeS())>1e-9) localOK=0;
  MPI_Allreduce(&localOK,&allOK,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  if (!allOK) return Status(StatusCode::ConfigurationConflict,"native propagation PIC/runtime clocks or rank speeds disagree");
  Output::ShockHistorySample result;
  // Use identical reduced values on every rank. In particular the radius-stop
  // branch in main.cpp must never diverge for a target within roundoff spread.
  result.timeS=minimum[0]; result.tick=identity[0]; result.radiusM=minimum[1]; result.speedMPerS=minimum[2];
  result.active=active; result.generation=generation;
  result.particles=globalCounts[0]; result.injections=globalCounts[1];
  result.mpiRadiusSpreadM=maximum[1]-minimum[1]; result.mpiClockSpreadS=maximum[0]-minimum[0];
  result.providerIdentity=providerIdentity;
  result.configurationFingerprint=providerFingerprint;
  Status status=Output::CheckShockHistorySample(result,options.requestedTimeStepS);
  if (status.ok()) *sample=std::move(result);
  return status;
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

  // The reduced provider owns both ambient and surface epochs.  After the
  // joined background commit, update the mover's finite-SSE crossing geometry
  // from that exact committed epoch before any source or population-control
  // operation can observe the next step.  No legacy SWCME object is created.
  if(!gInstalledShock&&Configuration().options().source.enabled) {
    const auto reduced=std::dynamic_pointer_cast<
        SEP3D::Adapters::ShockFrontBackgroundAdapter>(
            gRuntimeBackgroundProvider);
    const auto epoch=reduced?reduced->FrontEpoch():nullptr;
    const auto provider=reduced?reduced->SharedProvider():nullptr;
    if(!epoch||!provider)
      StopWithStatus("reduced-front source epoch update",SEP3D::Core::Status(
          SEP3D::Core::StatusCode::SnapshotUnavailable,
          "committed reduced-front epoch disappeared after background refresh"));
    status=SEP3D::AMPS::Movers::UpdateShock(
        ReducedMoverShock(*epoch,provider->Event()));
    if(!status.ok())StopWithStatus("reduced-front mover shock update",status);
  }

  // R05 source creation is a joined-boundary operation.  The provider state,
  // stochastic rounding, AMPS allocation, and conservation ledger therefore
  // all refer to the same authoritative integer tick.
  if (gInstalledShock) {
    const SEP3D::Adapters::ShockState shock =
        gInstalledShock->Evaluate(runtime.CurrentTimeS());
    if (!shock.status.ok()) StopWithStatus("shock update", shock.status);
    const SEP3D::Adapters::ExpandingShock moverShock=shock.MoverGeometry();
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
      // Domain pruning can remove most of a spherical SWCME front.  Determine
      // connectivity from the replicated AMR active-use flag before allocating
      // the exact macro count.  Counts are apportioned only over represented
      // patches, while every omitted physical patch receives an explicit
      // disconnected ledger row instead of failing later in InitiateParticle.
      std::vector<std::size_t> connectedOrdinal(
          shock.patches.size(), std::numeric_limits<std::size_t>::max());
      std::vector<SEP3D::Adapters::ShockSourceRecord> connectedPatches;
      connectedPatches.reserve(shock.patches.size());
      for (std::size_t patchIndex = 0;
           patchIndex < shock.patches.size(); ++patchIndex) {
        double position[3];
        shock.patches[patchIndex].positionM.CopyTo(position);
        cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node =
            PIC::Mesh::mesh->findTreeNode(position);
        if (node != nullptr && node->IsUsedInCalculationFlag) {
          connectedOrdinal[patchIndex] = connectedPatches.size();
          connectedPatches.push_back(shock.patches[patchIndex]);
        }
      }
      if (connectedPatches.empty()) {
        StopWithStatus("shock/corridor connectivity", SEP3D::Core::Status(
            SEP3D::Core::StatusCode::ConfigurationConflict,
            "no active SWCME source patch intersects the allocated Parker "
            "corridor"));
      }
      // `samples_per_step` is an exact count per compiled species.  Allocate
      // it independently across the connected part of the shock surface so
      // every SpeciesList entry is injected; sharing one count among species
      // would invent an undeclared composition.
      for (const auto& compiled : gCompiledSpecies) {
        std::vector<std::uint64_t> exactPatchCounts;
        if (options.inputSchemaVersion >= 3) {
          status = SEP3D::Adapters::AllocateExactPatchMacroparticles(
              connectedPatches, options.source.samplesPerStep,
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
              options.source.maximumEnergyJ, options.source.spectrumModel,
              options.source.fixedPhaseSpacePowerIndex);
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
          source.macroparticleWeight =
              ConfiguredParticleWeight(compiled.ampsIndex);
          const std::size_t allocationIndex =
              connectedOrdinal[patchIndex];
          source.connected = allocationIndex !=
              std::numeric_limits<std::size_t>::max();
          if (options.inputSchemaVersion >= 3) {
            source.prescribedMacroparticles = source.connected
                ? exactPatchCounts[allocationIndex] : 0;
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
  // Apply resampling after this boundary's source injection.  Running it
  // before injection would let a strong shock source immediately exceed the
  // declared maximum and carry that excess through the entire next transport
  // step.  At this location all movers and source allocators are quiescent,
  // while observers/checkpoints below see the controlled, weight-conserving
  // representation that will enter the next iteration.
  ApplyPopulationControlAtBoundary();
  // R06 ordering is intentional: observers see particles after both motion
  // and this boundary's shock injection.  Publication is globally gathered,
  // deterministic by stable ID, and its accumulator clears only on commit.
  PublishObserversAtBoundary();
  WriteCheckpointAtBoundary();
  return returnCode;
}

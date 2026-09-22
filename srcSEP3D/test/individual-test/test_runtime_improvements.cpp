// ============================================================================
// R01-R09 executable acceptance tests.
//
// Each callback targets one improvement and stays AMPS/MPI independent.  R01
// audits the generated-header hook text; R02-R09 exercise the exact portable
// service invoked by the production AMPS boundary.  Keeping these in the
// native registry means ``test/run_tests.py --all`` cannot omit the new work.
// ============================================================================

#include "sep3d_test_registry.h"

#include "observer_runtime.h"
#include "restart.h"
#include "runtime_adapters.h"
#include "source_runtime.h"
#include "configuration_io.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>

namespace {

namespace A = SEP3D::Adapters;
namespace C = SEP3D::Core;
namespace O = SEP3D::Output;
namespace R = SEP3D::RuntimeModel;
namespace fs = std::filesystem;
using Result = SEP3D::Testing::Result;

Result Pass(const std::string& message) {
  Result value; value.status = SEP3D::Testing::Status::Pass;
  value.message = message; return value;
}
Result Fail(const std::string& message) {
  Result value; value.status = SEP3D::Testing::Status::Fail;
  value.message = message; return value;
}

std::shared_ptr<const R::RunConfiguration3D> Configuration(
    std::uint64_t background = 1, std::uint64_t injection = 1,
    std::uint64_t output = 1, std::uint64_t checkpoint = 0) {
  R::RunConfiguration3DOptions options;
  options.backgroundCadenceSteps = background;
  options.injectionCadenceSteps = injection;
  options.outputCadenceSteps = output;
  options.checkpointCadenceSteps = checkpoint;
  std::shared_ptr<const R::RunConfiguration3D> result;
  if (!R::RunConfiguration3D::Create(options, &result).ok()) return {};
  return result;
}

R::SnapshotDescriptor Snapshot(const R::RunConfiguration3D& cfg,
                               std::uint64_t generation,
                               double epoch, double until = 100.0) {
  R::SnapshotDescriptor value;
  value.authority = cfg.options().background;
  value.epochS = epoch; value.validFromS = 0.0; value.validUntilS = until;
  value.generation = generation; value.complete = true;
  value.coordinateFrame = cfg.options().coordinateFrame;
  value.providerIdentity = "R01-R07-test-provider";
  value.configurationFingerprint = cfg.physics_fingerprint();
  return value;
}

bool Ready(R::Runtime* runtime,
           const std::shared_ptr<const R::RunConfiguration3D>& cfg) {
  if (!runtime->Configure(cfg).ok()) return false;
  R::MeshBinding binding; binding.layout = cfg->storage_layout();
  if (!runtime->BindMesh(binding).ok()) return false;
  R::StandaloneAdapter adapter;
  return adapter.Initialize(runtime).ok() &&
      runtime->PublishSnapshot(Snapshot(*cfg, 1, 0.0)).ok();
}

A::MoverInput Mover() {
  A::MoverInput input;
  input.model = R::TransportModel::Parker3D;
  input.particle.stableId = 9001; input.particle.species = 0;
  input.particle.positionM = C::Vec3(10.0, 0.0, 0.0);
  input.particle.statisticalWeight = 1.0;
  input.speciesMassKg = C::Const::m_p;
  input.requestedDtS = 1.0; input.innerRadiusM = 1.0;
  input.outerRadiusM = 100.0; input.campaignSeed = 17;
  input.timeStepControls.cellCrossingFraction = 0.25;
  return input;
}

struct ResolverRecord {
  std::uint64_t calls = 0;
  double lastX = 0.0;
};

C::Status Resolve(const A::ParticleRecord& particle, double elapsed,
                  void* opaque, A::LocalTransportRecord* local) {
  ResolverRecord* record = static_cast<ResolverRecord*>(opaque);
  if (record == nullptr || local == nullptr || elapsed < 0.0 ||
      particle.positionM.x < record->lastX)
    return C::Status(C::StatusCode::InvalidInput, "resolver order is invalid");
  ++record->calls; record->lastX = particle.positionM.x;
  local->background.status = C::Status::OK();
  local->background.valid = true;
  local->background.B = C::Vec3(1.0, 0.0, 0.0);
  local->background.absB = 1.0;
  local->background.bHat = C::Vec3(1.0, 0.0, 0.0);
  local->background.U = C::Vec3(1.0, 0.0, 0.0);
  local->cellSizeM = 0.2;
  local->timeToSnapshotBoundaryS = 10.0 - elapsed;
  return C::Status::OK();
}

swcme::sep::SEPSourceState CommonSource() {
  swcme::sep::SEPSourceState source;
  source.status = swcme::ModelStatus::success(); source.active = true;
  source.source_id = 5; source.position_m = {{10.0, 0.0, 0.0}};
  source.normal = {{1.0, 0.0, 0.0}}; source.relative_patch_weight = 1.0;
  source.compression = 4.0; source.normal_speed_m_s = 8.0e5;
  source.q_phase_space = 4.0;
  source.spectrum.particle_mass_kg = C::Const::m_p;
  source.spectrum.kinetic_energy_min_MeV = 1.0;
  source.spectrum.kinetic_energy_max_MeV = 10.0;
  source.spectrum.reference_energy_MeV = 2.0;
  return source;
}

Result RunR3D01() {
  std::ifstream hook("amps/install_mover_hook.py");
  std::ifstream makefile("makefile");
  std::ostringstream hookText, makeText;
  hookText << hook.rdbuf(); makeText << makefile.rdbuf();
  const std::string h = hookText.str(), m = makeText.str();
  if (!hook || !makefile ||
      h.find("SEP3D::AMPS::Movers::MoveParticle(ptr,LocalTimeStep,node)") ==
          std::string::npos ||
      h.find("int MoveParticle(long int ptr, double dtTotal") ==
          std::string::npos ||
      m.find("strict-production: prepare-production") == std::string::npos)
    return Fail("generated picGlobal hook or signature declaration is absent");
  return Pass("production build installs one declared SEP3D AMPS mover before pic_mover compilation");
}

Result RunR3D02() {
  ResolverRecord record;
  A::RequestedTimeAdvance request;
  request.input = Mover(); request.resolveLocal = Resolve;
  request.resolverContext = &record; request.maximumSubsteps = 1000;
  const A::MoverResult moved = A::AdvanceParticleRequestedTime(request);
  if (!moved.status.ok() || moved.disposition != A::ParticleDisposition::Active ||
      moved.consumedTimeS != request.input.requestedDtS ||
      moved.acceptedSubsteps <= 1 || record.calls != moved.acceptedSubsteps ||
      std::fabs(moved.particle.positionM.x - 11.0) > 1.0e-12)
    return Fail("complete requested-time loop skipped time or reused local state");
  return Pass("requested host time is consumed by re-resolved deterministic accepted substeps");
}

Result RunR3D03() {
  const auto cfg = Configuration(); R::Runtime runtime;
  if (!cfg || !Ready(&runtime, cfg)) return Fail("could not prepare Runtime");
  if (!runtime.RequestSnapshotUpdate(2.0, 2).ok() ||
      !runtime.BeginSnapshotFill().ok()) return Fail("update did not enter Filling");
  R::SnapshotDescriptor wrong = Snapshot(*cfg, 3, 2.0);
  if (runtime.StageSnapshot(wrong).ok() ||
      runtime.active_snapshot()->generation != 1 ||
      !runtime.FailSnapshotUpdate("provider rejected generation").ok() ||
      runtime.active_snapshot()->generation != 1 ||
      !runtime.AcknowledgeSnapshotFailure().ok())
    return Fail("invalid staging mutated the active generation");
  if (!runtime.RequestSnapshotUpdate(2.0, 2).ok() ||
      !runtime.BeginSnapshotFill().ok() ||
      !runtime.StageSnapshot(Snapshot(*cfg, 2, 2.0)).ok() ||
      runtime.PublishStagedSnapshot(false).ok() ||
      runtime.active_snapshot()->generation != 1 ||
      !runtime.PublishStagedSnapshot(true).ok() ||
      runtime.active_snapshot()->generation != 2)
    return Fail("collective staged publication was not atomic");
  return Pass("failed fills preserve the active snapshot and collective publish swaps one complete generation");
}

Result RunR3D04() {
  const auto cfg = Configuration(2, 3, 2, 4); R::Runtime runtime;
  if (!cfg || !Ready(&runtime, cfg)) return Fail("could not prepare clock Runtime");
  R::ClockObservation clocks{0.0, 1.0, 0.0, 0.0};
  if (!runtime.VerifyClockAgreement(clocks).ok() || runtime.BeginStep(1.0).ok())
    return Fail("clock agreement accepted a one-tick host offset");
  for (std::uint64_t tick = 0; tick < 4; ++tick) {
    if (!runtime.BeginStep(static_cast<double>(tick)).ok() ||
        !runtime.CompleteStep().ok()) return Fail("integer clock step failed");
  }
  if (runtime.counters().currentTick != 4 ||
      !runtime.EventDue(R::ScheduledEvent::Background) ||
      runtime.EventDue(R::ScheduledEvent::Injection) ||
      !runtime.EventDue(R::ScheduledEvent::Sampling) ||
      !runtime.EventDue(R::ScheduledEvent::Checkpoint) ||
      runtime.event_schedule().nextBackgroundTick != 6)
    return Fail("integer event schedule drifted from configured cadences");
  R::Runtime restored;
  if (!restored.Configure(cfg).ok() ||
      !restored.RestoreCounters(runtime.counters()).ok() ||
      !restored.RestoreEventSchedule(runtime.event_schedule()).ok() ||
      restored.CurrentTimeS() != runtime.CurrentTimeS())
    return Fail("restart did not restore the authoritative clock and schedule");
  return Pass("PIC-style clock observations, integer events, and restored schedule remain aligned");
}

Result RunR3D05() {
  const A::ShockSourceRecord patch = A::MakeShockSourceRecord(
      CommonSource(), 9, 91, 16, 0.2);
  // Equal total kinetic-energy bounds must map to different momentum bounds
  // for an electron and proton.  This is the regression guard against reusing
  // SWCME's reference-particle interval for every compiled AMPS species.
  A::ShockSourceRecord protonPatch = patch;
  A::ShockSourceRecord electronPatch = patch;
  const double minimumEnergyJ = 1.0e-15;
  const double maximumEnergyJ = 1.0e-12;
  if (!A::ConfigureSpeciesSpectrum(
           &protonPatch, C::Const::m_p, minimumEnergyJ,
           maximumEnergyJ).ok() ||
      !A::ConfigureSpeciesSpectrum(
           &electronPatch, C::Const::m_e, minimumEnergyJ,
           maximumEnergyJ).ok() ||
      !(electronPatch.injection.spectrum.minimum <
        protonPatch.injection.spectrum.minimum) ||
      electronPatch.sourceFingerprint == protonPatch.sourceFingerprint)
    return Fail("species mass did not produce a distinct valid momentum interval");
  A::SourceRequest request;
  request.patch = protonPatch; request.step = 4; request.species = 0;
  request.speciesMassKg = C::Const::m_p; request.intervalS = 2.0;
  request.physicalParticleRatePerS = 10.0;
  request.macroparticleWeight = 3.0; request.maximumMacroparticles = 4;
  const A::InjectionPlan plan = A::BuildInjectionPlan(request);
  double represented = 0.0;
  for (const auto& p : plan.particles)
    represented += p.particle.statisticalWeight;
  A::SourceRequest nextRequest = request;
  ++nextRequest.step;
  const A::InjectionPlan nextPlan = A::BuildInjectionPlan(nextRequest);
  request.connected = false;
  const A::InjectionPlan disconnected = A::BuildInjectionPlan(request);
  if (!plan.status.ok() || plan.particles.size() != 4 ||
      plan.ledger.capped == 0 || std::fabs(represented - 20.0) > 1.0e-12 ||
      !(plan.ledger.injectedEnergyJ > 0.0) ||
      !nextPlan.status.ok() || nextPlan.particles.empty() ||
      plan.particles[0].particle.completedStep != request.step ||
      nextPlan.particles[0].particle.completedStep != nextRequest.step ||
      plan.particles[0].particle.stableId ==
          nextPlan.particles[0].particle.stableId ||
      !disconnected.status.ok() || !disconnected.particles.empty() ||
      disconnected.ledger.disconnectedPatches != 1)
    return Fail("source rate, cadence identity, cap, weight, or disconnected policy is incorrect");
  return Pass("species-specific spectra and shock source planning conserve represented number and cadence identity");
}

Result RunR3D06() {
  R::RunConfiguration3DOptions options;
  options.observers[0].kind = R::ObserverKind::MovingCartesian;
  options.observers[0].velocityMPerS = C::Vec3(2.0, 0.0, 0.0);
  options.observers[0].collectionRadiusM = 5.0;
  options.observers[0].minimumEnergyJ = 1.0e-20;
  options.observers[0].maximumEnergyJ = 1.0e-10;
  options.observers[0].energyBins = 4;
  options.observers[0].energyChannelSpacing = R::EnergyChannelSpacing::Linear;
  std::shared_ptr<const R::RunConfiguration3D> cfg;
  if (!R::RunConfiguration3D::Create(options, &cfg).ok())
    return Fail("observer configuration was rejected");
  std::vector<O::VirtualSpacecraftDefinition> definitions;
  if (!O::BuildObserverDefinitions(*cfg, 10.0, &definitions).ok() ||
      definitions.size() != 1 ||
      definitions[0].positionM.x != options.observers[0].positionM.x + 20.0 ||
      definitions[0].kineticEnergyEdgesJ.size() != 5 ||
      std::fabs(definitions[0].kineticEnergyEdgesJ[1] -
                (1.0e-20 + (1.0e-10 - 1.0e-20) / 4.0)) > 1.0e-25)
    return Fail("moving observer was not resolved at the authoritative time");
  O::SamplingRequest sample;
  sample.cells.push_back({1, definitions[0].positionM, 1.0});
  O::ParticleObservation particle;
  particle.stableId = 1; particle.cellId = 1; particle.species = 0;
  particle.positionM = definitions[0].positionM;
  particle.momentumKgMPerS = 1.0e-19; particle.restMassKg = C::Const::m_p;
  particle.mu = 0.25; particle.statisticalWeight = 4.0;
  sample.particles.push_back(particle); sample.spacecraft = definitions;
  O::ObserverRuntime observer;
  if (!observer.Capture(sample).ok()) return Fail("observer capture failed");
  const O::SamplingSnapshot prepared = observer.PreparePublication();
  observer.AbortPublication();
  if (!prepared.status.ok() || prepared.spacecraft.size() != 1 ||
      prepared.spacecraft[0].standardUncertaintyPerJ.empty() ||
      !observer.has_pending_window() ||
      !observer.CommitPublication(prepared).ok() ||
      observer.has_pending_window() || observer.state().pendingWindows != 0)
    return Fail("observer uncertainty or commit/abort transaction is invalid");
  return Pass("moving observer products carry acceptance/uncertainty metadata and reset only after commit");
}

fs::path TemporaryRestart() {
  return fs::temp_directory_path() /
      ("srcsep3d-r07-" + std::to_string(
          std::chrono::high_resolution_clock::now().time_since_epoch().count()) +
       ".chk");
}

Result RunR3D07() {
  const auto cfg = Configuration(2, 3, 2, 4);
  if (!cfg) return Fail("could not build restart configuration");
  O::RestartState state;
  state.configurationFingerprint = cfg->physics_fingerprint();
  state.resolvedConfigurationManifest = cfg->resolved_manifest();
  state.storageLayoutFingerprint = cfg->storage_layout().fingerprint;
  state.codeIdentity = "R07-test-code"; state.snapshotFingerprint = "snap-4";
  state.runtimeCounters = {4, 0, 2, 1, 4};
  state.eventSchedule = {6, 6, 6, 8};
  state.activeSnapshot = Snapshot(*cfg, 4, 4.0);
  state.baseTimeStepS = 1.0; state.backgroundGeneration = 4;
  state.turbulenceGeneration = 5; state.sourceGeneration = 6;
  state.campaignSeed = 19; state.nextStableParticleId = 9;
  state.savedRankCount = 2;
  state.shockState.status = C::Status::OK(); state.shockState.active = true;
  state.shockState.generation = 6; state.shockState.epochS = 4.0;
  state.shockState.validUntilS = 5.0; state.shockState.radiusM = 20.0;
  state.shockState.radialSpeedMPerS = 1.0;
  state.shockState.compressionRatio = 4.0;
  state.shockState.providerIdentity = "swcme";
  state.shockState.configurationFingerprint = cfg->physics_fingerprint();
  state.samplingState = {2, 10, 1, 3, 7.0};
  A::ParticleRecord particle; particle.stableId = 7; particle.species = 0;
  particle.positionM = C::Vec3(10.0, 0.0, 0.0);
  particle.momentumKgMPerS = 1.0e-19; particle.mu = 0.1;
  particle.statisticalWeight = 1.0; particle.completedStep = 4;
  particle.substep = 12; state.particles.push_back(particle);
  A::SourceLedgerRow source; source.step = 4; source.species = 0;
  source.shockGeneration = 6; source.representedParticles = 7.0;
  source.injectedEnergyJ = 1.0; source.macroparticles = 3;
  state.sourceLedgerRows.push_back(source);
  const fs::path path = TemporaryRestart();
  const C::Status written = O::WriteRestart(path.string(), state);
  O::RestartLoadOptions options;
  options.expectedConfigurationFingerprint = cfg->physics_fingerprint();
  options.expectedResolvedConfigurationManifest = cfg->resolved_manifest();
  options.expectedStorageLayoutFingerprint = cfg->storage_layout().fingerprint;
  options.expectedCodeIdentity = "R07-test-code";
  options.expectedSnapshotFingerprint = "snap-4";
  options.availableBackgroundGeneration = 4;
  options.currentRankCount = 4;
  options.repartition = O::RepartitionPolicy::DeterministicByStableId;
  O::RestartState loaded;
  const C::Status read = O::ReadRestart(path.string(), options, &loaded);
  std::error_code ec; fs::remove(path, ec);
  if (!written.ok() || !read.ok() || loaded.runtimeCounters.currentTick != 4 ||
      loaded.eventSchedule.nextCheckpointTick != 8 ||
      loaded.samplingState.pendingRepresentedParticles != 7.0 ||
      loaded.sourceLedgerRows.size() != 1 ||
      loaded.particles[0].substep != 12)
    return Fail("versioned complete restart state did not round-trip");
  return Pass("restart round-trips clock, events, snapshots, shock, RNG tuple, ledgers, and pending sampling state");
}

Result RunR3D08() {
  R::RunConfiguration3DOptions options;
  if (!R::LoadConfigurationFile(
          "examples/sep3d_analytic_parker.in", &options).ok())
    return Fail("could not resolve the schema-v3 initialization fixture");
  std::shared_ptr<const R::RunConfiguration3D> configuration;
  if (!R::RunConfiguration3D::Create(options, &configuration).ok())
    return Fail("could not freeze the schema-v3 initialization fixture");
  std::shared_ptr<A::ShockProvider> provider;
  const C::Status created = A::CreateStandaloneSwcmeShockProvider(
      *configuration, &provider);
  if (!created.ok() || !provider)
    return Fail("canonical standalone SWCME3D provider was not created: " +
                created.message);
  // The fixture explicitly declares 3600 s as its first physically valid
  // source epoch; before that boundary the provider must remain inactive.
  const A::ShockState before = provider->Evaluate(3599.0);
  const A::ShockState state = provider->Evaluate(3600.0);
  if (!state.status.ok() || !state.active || state.patches.empty() ||
      !before.status.ok() || before.active ||
      state.configurationFingerprint !=
          options.swcmeConfigurationFingerprint)
    return Fail("canonical SWCME3D provider did not publish its initial shock surface");

  std::vector<std::uint64_t> counts;
  if (!A::AllocateExactPatchMacroparticles(
           state.patches, options.source.samplesPerStep, &counts).ok() ||
      counts.size() != state.patches.size() ||
      std::accumulate(counts.begin(), counts.end(), UINT64_C(0)) !=
          options.source.samplesPerStep ||
      std::find(counts.begin(), counts.end(), UINT64_C(0)) != counts.end())
    return Fail("global exact source count was not conserved over all active patches");
  std::vector<std::uint64_t> repeated;
  if (!A::AllocateExactPatchMacroparticles(
           state.patches, options.source.samplesPerStep, &repeated).ok() ||
      repeated != counts)
    return Fail("largest-remainder source allocation is not deterministic");
  if (A::AllocateExactPatchMacroparticles(
          state.patches,
          static_cast<std::uint64_t>(state.patches.size() - 1),
          &repeated).ok())
    return Fail("source allocation silently omitted a non-zero physical patch");

  A::SourceRequest request;
  request.patch = state.patches.front();
  request.step = 0;
  request.species = 0;
  request.speciesMassKg = C::Const::m_p;
  request.intervalS = options.requestedTimeStepS;
  request.physicalParticleRatePerS =
      options.source.physicalParticleRatePerS *
      options.source.injectionEfficiency *
      request.patch.relativePatchWeight;
  request.macroparticleWeight = options.species.macroparticleWeight;
  request.prescribedMacroparticles = counts.front();
  request.maximumMacroparticles = counts.front();
  const A::InjectionPlan plan = A::BuildInjectionPlan(request);
  if (!plan.status.ok() || plan.particles.size() != counts.front() ||
      plan.ledger.macroparticles != counts.front() || plan.ledger.capped != 0)
    return Fail("an exact patch allocation was rounded or capped downstream");
  return Pass("canonical SWCME initialization publishes a physical surface and allocates one exact per-species count per step");
}

Result RunR3D09() {
  const std::vector<double> physicalValues = {1.0, -2.0, 3.0};

  // Before the first completed AMPS sampling window, particle moments are
  // finite zeros but are explicitly marked as not-yet-sampled.
  const O::TecplotCellPresentation initialization =
      O::PrepareTecplotCellPresentation(physicalValues, true, 0, 0.0);
  if (initialization.backgroundValues != physicalValues ||
      initialization.backgroundValid != 1.0 ||
      initialization.particleSamplingWindowValid != 0.0 ||
      initialization.particleSamplePresent != 0.0)
    return Fail("initialization output confused an absent sampling window with invalid background");

  // A completed window with no particles is a valid zero-occupancy sample,
  // not a failed or undefined numerical result.
  const O::TecplotCellPresentation empty =
      O::PrepareTecplotCellPresentation(physicalValues, true, 8, 0.0);
  if (empty.backgroundValid != 1.0 ||
      empty.particleSamplingWindowValid != 1.0 ||
      empty.particleSamplePresent != 0.0)
    return Fail("empty particle cell was not represented as a valid finite sample");

  const O::TecplotCellPresentation occupied =
      O::PrepareTecplotCellPresentation(physicalValues, true, 8, 2.0);
  if (occupied.particleSamplingWindowValid != 1.0 ||
      occupied.particleSamplePresent != 1.0)
    return Fail("occupied particle cell lost its explicit presence flag");

  // Cartesian padding cells and corrupted diagnostic values must remain
  // parseable. Their zeros are placeholders guarded by backgroundValid=0.
  std::vector<double> invalidValues = physicalValues;
  invalidValues[1] = std::numeric_limits<double>::quiet_NaN();
  const O::TecplotCellPresentation invalid =
      O::PrepareTecplotCellPresentation(invalidValues, true, 8, 0.0);
  const O::TecplotCellPresentation padding =
      O::PrepareTecplotCellPresentation(physicalValues, false, 8, 0.0);
  const auto finiteZero = [](const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(),
                       [](double value) {
                         return std::isfinite(value) && value == 0.0;
                       });
  };
  if (invalid.backgroundValid != 0.0 || !finiteZero(invalid.backgroundValues) ||
      padding.backgroundValid != 0.0 || !finiteZero(padding.backgroundValues))
    return Fail("invalid background emitted a non-finite Tecplot placeholder");

  return Pass("Tecplot output separates invalid background, unsampled windows, empty cells, and occupied cells without NaN");
}

}  // namespace

std::vector<SEP3D::Testing::Descriptor> RegisterRuntimeImprovementTests() {
  using D = SEP3D::Testing::Descriptor;
  using I = SEP3D::Testing::InitializationLevel;
  using RC = SEP3D::Testing::RuntimeClass;
  auto make = [](const char* id, const char* name,
                 SEP3D::Testing::TestCallback callback) {
    D d; d.id = id; d.name = name; d.group = "R3D";
    d.description = "R01-R09 production-runtime improvement acceptance";
    d.initialization = I::None; d.supportedBuildModes = "standalone-no-AMPS";
    d.runtime = RC::Routine; d.seedPolicy = "semantic keyed or deterministic";
    d.stateIsolation = "fresh runtime/service and temporary artifacts";
    d.callback = std::move(callback); return d;
  };
  return {
      make("R3D01", "AMPS mover selection hook", RunR3D01),
      make("R3D02", "Complete requested-time advance", RunR3D02),
      make("R3D03", "Transactional snapshot update", RunR3D03),
      make("R3D04", "Authoritative integer clock", RunR3D04),
      make("R3D05", "Shock source lifecycle", RunR3D05),
      make("R3D06", "Observer publication transaction", RunR3D06),
      make("R3D07", "Complete restart contract", RunR3D07),
      make("R3D08", "Canonical initialization source", RunR3D08),
      make("R3D09", "Finite empty-cell Tecplot output", RunR3D09),
  };
}

// ============================================================================
// C01-C05 configuration, geometry, and resource-planning acceptance tests.
//
// These tests deliberately remain AMPS independent.  They exercise the exact
// parser, immutable factory, Parker geometry, analytic background provider,
// composite resolution function, and dry-run estimator used by production.
// A failure therefore identifies a model-contract defect before MPI or the
// AMPS mesh can obscure it with host state.
// ============================================================================

#include "sep3d_test_registry.h"
#include "bg_parker.h"
#include "mesh_model.h"
#include "configuration_io.h"
#include "application_input.h"
#include "particle_normalization.h"
#include "background_factory.h"
#include "shock_front_background_adapter.h"
#include "runtime.h"
#include "sep_coronal_cme/constants.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace {

namespace BG = SEP3D::Background;
namespace M = SEP3D::Mesh;
namespace RM = SEP3D::RuntimeModel;
using Result = SEP3D::Testing::Result;

Result Pass(const std::string& message) {
  Result result;
  result.status = SEP3D::Testing::Status::Pass;
  result.message = message;
  return result;
}

Result Fail(const std::string& message) {
  Result result;
  result.status = SEP3D::Testing::Status::Fail;
  result.message = message;
  return result;
}

// A complete, annotated-file-shaped fixture.  Most values equal the typed C++
// defaults on purpose: CFG3D01 can then prove file and programmatic creation
// produce the same normalized fingerprint rather than merely similar values.
std::string CompleteInput() {
  return R"SEP3D(
[run]
schema_version = 1
intent = transport-only
transport = parker3d
time_step_s = 1
maximum_time_steps = 100000001
campaign_seed = 1
background_cadence_steps = 1

[domain]
preset = one-au
inner_radius_m = 1.3914e10
inner_boundary = absorb
outer_radius_mode = preset
outer_boundary = escape
coordinate_frame = HCI-like-inertial
origin_x_m = 0
origin_y_m = 0
origin_z_m = 0

[mesh]
global_cell_size_m = 3.7399467675e10
minimum_cell_size_m = 1.495978707e9
cells_per_block_edge = 4
maximum_level = 7
memory_budget_bytes = 1000000000000000
block_overhead_bytes = 1024

[mesh.solar]
enabled = true
surface_cell_size_m = 1.495978707e9
transition_outer_radius_m = 3.7399467675e10
profile = smoothstep
exponent = 1

[mesh.tube]
enabled = false
source_longitude_rad = 0
source_colatitude_rad = 1.5707963267948966
reference_radius_m = 1.495978707e11
radius_at_reference_m = 4.487936121e9
radius_mode = constant-angular-width
center_cell_size_m = 1.495978707e9
transverse_profile = smoothstep
transverse_exponent = 1

[memory]
base_cell_bytes = 256
base_node_bytes = 128
block_structure_bytes = 4096
communication_bytes_per_block = 2048
particle_bytes = 160
particles_per_cell = 2
halo_fraction = 0.2
safety_margin_fraction = 0.25

[background]
provider = analytic-parker
external_script = false

[background.parker]
reference_radius_m = 1.495978707e11
radial_field_at_reference_t = 3e-9
solar_rotation_rate_rad_per_s = 2.865e-6
solar_wind_speed_m_per_s = 400000
magnetic_polarity = 1
number_density_at_reference_m3 = 5e6
temperature_k = 1e5
validity_cadence_s = 3600

[turbulence]
authority = prescribed
delta_b_over_b = 0.3
k_min_per_m = 1e-10
k_max_per_m = 1e-7
spectral_index = 1.6666666666666667
correlation_length_m = 4.487936121e9
missing_data = fail
resonance_range = reject
self_consistent_3d = false

[transport]
cell_crossing_fraction = 0.4
diffusion_fraction = 0.2
focusing_fraction = 0.2
cooling_fraction = 0.2
field_variation_fraction = 0.2
shock_crossing_fraction = 0.5
minimum_substep_s = 1e-12
pitch_angle_scheme = reflecting-milstein
perpendicular_diffusion = false
drifts = false

[shock]
authority = none

[source]
enabled = false

[species]
macroparticle_weight = 1

[storage]
magnetic_gradient = false
velocity_gradient = false
sampling_bytes_per_cell = 0

[observer.default]
position_x_m = 3.7399467675e10
position_y_m = 0
position_z_m = 0
follows_trajectory = false
cadence_s = 60
energy_bins = 32
pitch_angle_bins = 24
products = flux,spectrum

[output]
cadence_steps = 1
directory = output
prefix = sep3d

[restart]
input_path = none
output_path = restart/sep3d.chk
)SEP3D";
}

// Schema version 2 deliberately adds only the finite-line section to the
// complete version-1 fixture.  This keeps the test sensitive to accidental
// changes in every pre-existing required section while exercising all eight
// newly required SI fields.
std::string CompleteVersion2Input() {
  std::string input = CompleteInput();
  const std::string oldVersion = "schema_version = 1";
  input.replace(input.find(oldVersion), oldVersion.size(),
                "schema_version = 2");
  const std::string marker = "\n[mesh]\n";
  const std::string section = R"SEP3D(
[parker_spiral]
origin_x_m = 0
origin_y_m = 0
origin_z_m = 0
initial_x_m = 1.3914e10
initial_y_m = 0
initial_z_m = 0
length_m = 2.0e11
point_count = 257
)SEP3D";
  input.insert(input.find(marker), section);
  return input;
}

M::ResolutionConfiguration Resolution(
    const RM::RunConfiguration3DOptions& options) {
  M::ResolutionConfiguration r;
  r.originM = options.coordinateOriginM;
  r.innerRadiusM = options.innerRadiusM;
  r.outerRadiusM = options.outerRadiusM;
  r.minimumCellSizeM = options.minimumCellSizeM;
  r.backgroundCellSizeM = options.backgroundCellSizeM;
  r.enableRadialRefinement = options.enableRadialRefinement;
  r.solarSurfaceCellSizeM = options.solarSurfaceCellSizeM;
  r.solarRefinementOuterRadiusM = options.solarRefinementOuterRadiusM;
  r.solarRefinementProfile = options.solarRefinementProfile;
  r.solarRefinementExponent = options.solarRefinementExponent;
  r.enableTubeRefinement = options.enableTubeRefinement;
  r.tubeLongitudeRad = options.tubeLongitudeRad;
  r.tubeColatitudeRad = options.tubeColatitudeRad;
  r.tubeReferenceRadiusM = options.tubeReferenceRadiusM;
  r.tubeRadiusAtReferenceM = options.tubeRadiusAtReferenceM;
  r.tubeRadiusMode = options.tubeRadiusMode;
  r.tubeCellSizeM = options.tubeCellSizeM;
  r.tubeTransverseProfile = options.tubeTransverseProfile;
  r.tubeTransverseExponent = options.tubeTransverseExponent;
  r.solarWindSpeedMPerS = options.parker.solarWindSpeedMPerS;
  r.solarRotationRateRadPerS = options.parker.solarRotationRateRadPerS;
  r.parkerInitialPointM = options.parkerSpiralInitialPointM;
  r.parkerLengthM = options.parkerSpiralLengthM;
  r.parkerPointCount = options.parkerSpiralPointCount;
  r.cellsPerBlockEdge = options.meshCellsPerBlockEdge;
  r.maximumLevel = options.maximumMeshLevel;
  r.blockOverheadBytes = options.meshBlockOverheadBytes;
  r.memoryBudgetBytes = options.meshMemoryBudgetBytes;
  r.memoryModel = options.memoryModel;
  return r;
}

Result RunCFG3D01() {
  RM::RunConfiguration3DOptions parsed;
  if (!RM::ParseConfigurationText(CompleteInput(), &parsed).ok())
    return Fail("the complete version-1 input fixture did not parse");
  std::shared_ptr<const RM::RunConfiguration3D> fromFile;
  if (!RM::RunConfiguration3D::Create(parsed, &fromFile).ok())
    return Fail("parsed options did not pass immutable construction");

  RM::RunConfiguration3DOptions programmatic;
  programmatic.meshMemoryBudgetBytes = 1000000000000000ULL;
  std::shared_ptr<const RM::RunConfiguration3D> fromTyped;
  if (!RM::RunConfiguration3D::Create(programmatic, &fromTyped).ok() ||
      fromFile->physics_fingerprint() != fromTyped->physics_fingerprint()) {
    return Fail("file and typed construction did not normalize to one physics fingerprint");
  }

  const char* argv[] = {"srcSEP3D", "--input", "run.in",
                        "--initialization-only",
                        "--initialization-output-dir", "preview",
                        "--output-dir", "products", "--log-level", "verbose"};
  RM::StandaloneCommandLine cli;
  if (!RM::ParseStandaloneCommandLine(10, const_cast<char**>(argv), &cli).ok() ||
      !cli.initializationOnly || cli.dryRun || cli.inputPath != "run.in" ||
      cli.initializationOutputDirectory != "preview" ||
      cli.outputDirectoryOverride != "products" ||
      cli.verbosity != RM::LogVerbosity::Verbose) {
    return Fail("documented standalone CLI options did not normalize correctly");
  }
  const char* conflicting[] = {"srcSEP3D", "--input", "run.in",
                               "--dry-run", "--initialization-only"};
  if (RM::ParseStandaloneCommandLine(
          5, const_cast<char**>(conflicting), &cli).ok()) {
    return Fail("dry-run and initialization-only modes were not rejected");
  }

  // Linked tests are a first-class executable mode.  They must carry an
  // explicit immutable deck plus deterministic evidence locations; these
  // switches are parsed here without importing AMPS or starting MPI.
  const char* native[] = {
      "srcSEP3D", "--all-tests", "--test-input", "run.in",
      "--test-json", "evidence.json", "--artifact-directory", "artifacts",
      "--test-steps", "2", "--expect-mpi-ranks", "4"};
  if (!RM::ParseStandaloneCommandLine(
          12, const_cast<char**>(native), &cli).ok() ||
      !cli.allTests || cli.testInputPath != "run.in" ||
      cli.testJsonPath != "evidence.json" ||
      cli.testArtifactDirectory != "artifacts" || cli.testSteps != 2 ||
      cli.expectedMpiRanks != 4) {
    return Fail("native linked-test CLI options did not normalize correctly");
  }
  const char* orphanEvidence[] = {
      "srcSEP3D", "--test-json", "evidence.json"};
  if (RM::ParseStandaloneCommandLine(
          3, const_cast<char**>(orphanEvidence), &cli).ok()) {
    return Fail("native evidence options were accepted without a test selector");
  }
  const char* conflictingDecks[] = {
      "srcSEP3D", "--test", "SCCM3D01", "--input", "a.in",
      "--test-input", "b.in"};
  if (RM::ParseStandaloneCommandLine(
          7, const_cast<char**>(conflictingDecks), &cli).ok()) {
    return Fail("native test accepted two different immutable input decks");
  }

  std::string bad = CompleteInput();
  bad += "\n[output]\nunknown_key = x\n";
  if (RM::ParseConfigurationText(bad, &parsed).ok())
    return Fail("duplicate/unknown input was not rejected before initialization");
  // Molecular identity belongs exclusively to the AMPS build-time
  // SpeciesList.  A post-compile attempt to redefine mass must be an unknown
  // key rather than a silently ignored compatibility spelling.
  std::string retiredSpeciesField = CompleteInput();
  const std::string weightLine = "macroparticle_weight = 1\n";
  retiredSpeciesField.insert(
      retiredSpeciesField.find(weightLine) + weightLine.size(),
      "mass_kg = 1.67262192369e-27\n");
  if (RM::ParseConfigurationText(retiredSpeciesField, &parsed).ok())
    return Fail("runtime input was allowed to redefine compiled species mass");
  const std::size_t observer = bad.find("[observer.default]");
  bad = CompleteInput();
  bad.erase(bad.find("[observer.default]"),
            bad.find("[output]") - bad.find("[observer.default]"));
  if (observer == std::string::npos || RM::ParseConfigurationText(bad, &parsed).ok())
    return Fail("a missing required observer group was not rejected");
  const SEP3D::Core::Status flatStatus =
      RM::ParseConfigurationText("scattering = prescribed\n", &parsed);
  if (flatStatus.ok() ||
      flatStatus.message.find("before any [section]") == std::string::npos) {
    return Fail("a legacy flat deck did not receive the section-aware diagnostic");
  }

  std::string summary;
  if (!RM::BuildDryRunSummary(*fromFile, &summary).ok() ||
      summary.find("physics_fingerprint=") == std::string::npos ||
      summary.find("estimated_blocks_by_level=") == std::string::npos) {
    return Fail("dry-run did not emit fingerprint and resource preflight data");
  }
  // The shipped example is executable input, not documentation-like prose.
  // Keeping it inside the acceptance test prevents a renamed key or tightened
  // validator from leaving users with an example that only looks plausible.
  RM::RunConfiguration3DOptions exampleOptions;
  const SEP3D::Core::Status exampleStatus = RM::LoadConfigurationFile(
      "examples/sep3d_analytic_parker.in", &exampleOptions);
  if (!exampleStatus.ok()) {
    return Fail("the annotated production example is invalid: " +
                exampleStatus.message);
  }
  std::shared_ptr<const RM::RunConfiguration3D> example;
  if (!RM::RunConfiguration3D::Create(exampleOptions, &example).ok() ||
      !RM::BuildDryRunSummary(*example, &summary).ok()) {
    return Fail("the annotated production example failed resource preflight");
  }
  if (exampleOptions.observers.size() != 2 ||
      !std::all_of(exampleOptions.observers.begin(),
                   exampleOptions.observers.end(),
                   [](const RM::ObserverOptions& configured) {
                     return configured.allCompiledSpecies &&
                            configured.species.empty();
                   })) {
    return Fail("the production example does not preserve species=all as a wildcard");
  }
  if (!RM::ApplyInitializationOutputDirectory("preview", &exampleOptions).ok() ||
      exampleOptions.initializationMeshTecplotFile !=
          "preview/sep3d-initialization-mesh.dat" ||
      exampleOptions.initializationParkerLineTecplotFile !=
          "preview/sep3d-initialization-parker-line.dat" ||
      exampleOptions.initializationDataTecplotFile !=
          "preview/sep3d-initialization-data.dat") {
    return Fail("initialization output-directory override changed product names");
  }
  return Pass("versioned input, CLI normalization, early errors, typed parity, and allocation-free dry-run passed");
}

Result RunCFG3D02() {
  RM::RunConfiguration3DOptions baseline;
  std::shared_ptr<const RM::RunConfiguration3D> a, b, c;
  if (!RM::RunConfiguration3D::Create(baseline, &a).ok())
    return Fail("baseline typed contract is invalid");
  RM::RunConfiguration3DOptions outputOnly = baseline;
  outputOnly.outputDirectory = "different-products";
  outputOnly.outputPrefix = "different-prefix";
  outputOnly.restartOutputPath = "different.chk";
  if (!RM::RunConfiguration3D::Create(outputOnly, &b).ok() ||
      a->physics_fingerprint() != b->physics_fingerprint())
    return Fail("output-only fields changed the physics fingerprint");
  RM::RunConfiguration3DOptions physics = baseline;
  physics.species.macroparticleWeight = 2.0;
  if (!RM::RunConfiguration3D::Create(physics, &c).ok() ||
      a->physics_fingerprint() == c->physics_fingerprint())
    return Fail("a particle/source contract change did not change the physics fingerprint");

  RM::RunConfiguration3DOptions incompatible = baseline;
  incompatible.intent = RM::RunIntent::ShockInjection;
  if (RM::RunConfiguration3D::Create(incompatible, &c).ok())
    return Fail("shock-injection intent accepted an absent shock/source");
  RM::RunConfiguration3DOptions coupled = baseline;
  coupled.background = RM::BackgroundAuthority::Swmf;
  coupled.turbulence = RM::TurbulenceAuthority::Swmf;
  if (!RM::RunConfiguration3D::Create(coupled, &c).ok())
    return Fail("typed SWMF authority could not use the shared immutable factory");
  RM::RunConfiguration3DOptions extension = baseline;
  extension.perpendicularDiffusion =
      RM::PerpendicularDiffusionMode::ConstantRatio;
  extension.kappaPerpendicularToParallelRatio = 0.02;
  extension.drift = RM::DriftMode::GradientAndCurvature;
  if (!RM::RunConfiguration3D::Create(extension, &c).ok() ||
      !c->options().storeMagneticGradient ||
      a->physics_fingerprint() == c->physics_fingerprint())
    return Fail("V01 closure did not freeze gradients or enter the physics fingerprint");
  return Pass("typed groups, field classifications, incompatibility checks, and analytic/SWMF factory parity passed");
}

Result RunCFG3D03() {
  RM::RunConfiguration3DOptions options;
  std::shared_ptr<const RM::RunConfiguration3D> solar, earth, mars, explicitRun;
  options.domain = RM::DomainPreset::Solar;
  if (!RM::RunConfiguration3D::Create(options, &solar).ok())
    return Fail("solar preset did not normalize");
  options.domain = RM::DomainPreset::OneAu;
  if (!RM::RunConfiguration3D::Create(options, &earth).ok())
    return Fail("one-AU preset did not normalize");
  options.domain = RM::DomainPreset::Mars;
  options.maximumMeshLevel = 8;
  if (!RM::RunConfiguration3D::Create(options, &mars).ok())
    return Fail("Mars preset did not normalize");
  options.outerRadiusMode = RM::OuterRadiusMode::Explicit;
  options.outerRadiusM = 0.75 * SEP3D::Core::Const::AU;
  options.domain = RM::DomainPreset::OneAu;
  options.maximumMeshLevel = 7;
  if (!RM::RunConfiguration3D::Create(options, &explicitRun).ok() ||
      solar->options().outerRadiusM != 0.30 * SEP3D::Core::Const::AU ||
      earth->options().outerRadiusM != SEP3D::Core::Const::AU ||
      mars->options().outerRadiusM != 1.666 * SEP3D::Core::Const::AU ||
      explicitRun->options().outerRadiusM != options.outerRadiusM) {
    return Fail("preset or explicit outer-radius resolution is incorrect");
  }

  RM::RunConfiguration3DOptions invalid;
  invalid.observers.front().positionM.x = 2.0 * SEP3D::Core::Const::AU;
  if (RM::RunConfiguration3D::Create(invalid, &explicitRun).ok())
    return Fail("observer outside the domain was accepted");

  // The photosphere is a fixed physical AMPS boundary, while innerRadiusM is
  // the independently configured Parker/CME source and transport shell.  The
  // source shell may coincide with the photosphere, but it cannot be placed
  // inside the solid Sun.
  RM::RunConfiguration3DOptions atSolarSurface;
  atSolarSurface.innerRadiusM = SEP3D::Core::Const::R_sun;
  if (!RM::RunConfiguration3D::Create(atSolarSurface, &explicitRun).ok())
    return Fail("a source shell exactly on the photosphere was rejected");
  const std::string& solarManifest = explicitRun->resolved_manifest();
  const std::string radiusKey = ";solar_radius_m=";
  const std::size_t radiusOffset = solarManifest.find(radiusKey);
  if (solarManifest.find(
          ";solar_boundary=amps-absorbing-sphere-v1") ==
          std::string::npos ||
      radiusOffset == std::string::npos) {
    return Fail(
        "resolved physics identity omits the fixed AMPS solar boundary");
  }
  const double manifestSolarRadiusM = std::stod(
      solarManifest.substr(radiusOffset + radiusKey.size()));
  if (manifestSolarRadiusM != SEP3D::Core::Const::R_sun) {
    return Fail("resolved physics identity records the wrong solar radius");
  }
  RM::RunConfiguration3DOptions belowSolarSurface = atSolarSurface;
  belowSolarSurface.innerRadiusM = std::nextafter(
      SEP3D::Core::Const::R_sun, 0.0);
  if (RM::RunConfiguration3D::Create(
          belowSolarSurface, &explicitRun).ok())
    return Fail("a source shell below the physical photosphere was accepted");

  const M::DomainBounds domain = M::MakeDomain(earth->options());
  const double inner = domain.innerRadiusM;
  const double outer = domain.outerRadiusM;
  if (M::ClassifyBoundaryCrossing({1.1 * inner, 0, 0},
                                  {0.9 * inner, 0, 0}, domain).code !=
          SEP3D::Core::StatusCode::InnerBoundary ||
      M::ClassifyBoundaryCrossing({0.9 * outer, 0, 0},
                                  {1.1 * outer, 0, 0}, domain).code !=
          SEP3D::Core::StatusCode::DomainExit ||
      !M::ClassifyBoundaryCrossing({1.1 * outer, 0, 0},
                                   {0.9 * outer, 0, 0}, domain).ok()) {
    return Fail("boundary crossing status ignored surface or travel direction");
  }
  return Pass("solar, one-AU, Mars, explicit-domain, photospheric-source ordering, containment, and directional boundary contracts passed");
}

Result RunCFG3D04() {
  RM::RunConfiguration3DOptions options;
  options.enableTubeRefinement = true;
  M::ResolutionConfiguration resolution = Resolution(options);
  const double radius = 0.6 * SEP3D::Core::Const::AU;
  const SEP3D::Core::Vec3 point = radius * M::ParkerTubeDirection(radius, resolution);
  const SEP3D::Core::Vec3 tangent = M::ParkerTubeTangent(radius, resolution);

  BG::ParkerConfiguration positive;
  positive.sourceRadiusM = options.innerRadiusM;
  positive.sourceLongitudeRad = options.tubeLongitudeRad;
  positive.sourceColatitudeRad = options.tubeColatitudeRad;
  BG::ParkerConfiguration negative = positive;
  positive.magneticPolarity = 1;
  negative.magneticPolarity = -1;
  BG::AnalyticParkerProvider plus(positive), minus(negative);
  if (!plus.Prepare(0.0).ok() || !minus.Prepare(0.0).ok())
    return Fail("Parker providers did not prepare");
  const BG::BackgroundSample bp = plus.Evaluate(point);
  const BG::BackgroundSample bm = minus.Evaluate(point);
  if (!bp.valid || !bm.valid || bp.bHat.Dot(tangent) < 1.0 - 1.0e-12 ||
      bm.bHat.Dot(tangent) > -1.0 + 1.0e-12 ||
      M::TubeDistanceM(point, resolution) > 1.0e-9 * radius) {
    return Fail("mesh tangent and analytic field disagree or polarity moved geometry");
  }
  return Pass("one polarity-free Parker geometry drives mesh position/tangent while field polarity only reverses B");
}

Result RunCFG3D05() {
  RM::RunConfiguration3DOptions options;
  options.enableTubeRefinement = true;
  options.meshMemoryBudgetBytes = 1000000000000000ULL;
  std::shared_ptr<const RM::RunConfiguration3D> configuration;
  if (!RM::RunConfiguration3D::Create(options, &configuration).ok())
    return Fail("composite mesh configuration is invalid");
  M::ResolutionConfiguration r = Resolution(configuration->options());
  const double inner = r.innerRadiusM;
  const double transition = r.solarRefinementOuterRadiusM;
  double previous = 0.0;
  for (int i = 0; i <= 100; ++i) {
    const double radius = inner + (transition - inner) * i / 100.0;
    const double cell = M::RequestedCellSizeM({radius, 0, 0}, r);
    if (i != 0 && cell + 1.0e-12 * r.backgroundCellSizeM < previous)
      return Fail("near-Sun degradation profile is not monotone");
    previous = cell;
  }
  const double referenceTube = M::TubeRadiusM(r.tubeReferenceRadiusM, r);
  const double halfRadiusTube = M::TubeRadiusM(0.5 * r.tubeReferenceRadiusM, r);
  if (referenceTube != r.tubeRadiusAtReferenceM ||
      std::fabs(halfRadiusTube - 0.5 * referenceTube) > 1.0e-12 * referenceTube)
    return Fail("constant-angular-width tube radius is not scaled from its declared reference");

  M::RefinementPreflight preflight;
  if (!M::BuildRefinementPreflight(M::MakeDomain(configuration->options()), r,
                                   configuration->storage_layout(),
                                   &preflight).ok() ||
      preflight.estimatedBlocksByLevel.size() != r.maximumLevel + 1 ||
      preflight.memory.baseCellBytes == 0 || preflight.memory.nodeBytes == 0 ||
      preflight.memory.particleBytes == 0 ||
      preflight.memory.totalBytes <= preflight.memory.subtotalBytes) {
    return Fail("preflight omitted level counts or a whole-run memory category");
  }

  RM::RunConfiguration3DOptions impossible = options;
  impossible.maximumMeshLevel = 1;
  if (RM::RunConfiguration3D::Create(impossible, &configuration).ok())
    return Fail("an unrealizable requested resolution passed maximum-level validation");
  return Pass("composite profiles, tube scaling, level preflight, full memory accounting, and level rejection passed");
}

Result RunCFG3D06() {
  RM::RunConfiguration3DOptions options;
  const SEP3D::Core::Status status =
      RM::ParseConfigurationText(CompleteVersion2Input(), &options);
  if (!status.ok() || options.inputSchemaVersion != 2 ||
      options.parkerSpiralPointCount != 257 ||
      options.parkerSpiralLengthM != 2.0e11 ||
      options.parkerSpiralInitialPointM.x != options.innerRadiusM) {
    return Fail("complete schema-version-2 Parker definition did not parse exactly");
  }

  std::string missing = CompleteVersion2Input();
  const std::string required = "point_count = 257\n";
  missing.erase(missing.find(required), required.size());
  if (RM::ParseConfigurationText(missing, &options).ok())
    return Fail("schema version 2 accepted a missing finite-line point count");

  std::string inconsistent = CompleteVersion2Input();
  const std::string initial = "initial_x_m = 1.3914e10";
  inconsistent.replace(inconsistent.find(initial), initial.size(),
                       "initial_x_m = 1.4e10");
  if (RM::ParseConfigurationText(inconsistent, &options).ok())
    return Fail("schema version 2 accepted a line source outside the inner boundary");
  return Pass("schema version 2 requires and validates the complete finite Parker-line definition");
}

Result RunCFG3D07() {
  // The post-compile runtime never chooses species identity.  Exercise the
  // neutral binding layer with the same records production reads from AMPS'
  // generated ChemTable/MolMass/ElectricChargeTable arrays.
  std::shared_ptr<const RM::RunConfiguration3D> configuration;
  RM::RunConfiguration3DOptions options;
  std::vector<RM::CompiledSpeciesRecord> compiled = {
      {0, "H_PLUS", SEP3D::Core::Const::m_p, SEP3D::Core::Const::e},
      {1, "ELECTRON", SEP3D::Core::Const::m_e, -SEP3D::Core::Const::e}};
  options.observers.front().allCompiledSpecies = true;
  options.observers.front().species.clear();
  if (!RM::ValidateCompiledSpeciesBinding(options, 2, compiled).ok())
    return Fail("an all-species observer rejected a mixed AMPS table");
  if (RM::ValidateCompiledSpeciesBinding(options, 1, compiled).ok())
    return Fail("a compiled species-count mismatch was accepted");
  std::vector<RM::CompiledSpeciesRecord> invalid = compiled;
  invalid[1].ampsIndex = 2;
  if (RM::ValidateCompiledSpeciesBinding(options, 2, invalid).ok())
    return Fail("a non-contiguous compiled species index was accepted");
  invalid = compiled;
  invalid[1].symbol = "h_plus";
  if (RM::ValidateCompiledSpeciesBinding(options, 2, invalid).ok())
    return Fail("a duplicate compiled chemical symbol was accepted");
  invalid = compiled;
  invalid[1].chargeC = 0.0;
  if (RM::ValidateCompiledSpeciesBinding(options, 2, invalid).ok())
    return Fail("a neutral species was accepted by charged SEP transport");
  invalid = compiled;
  invalid[1].massKg = 0.0;
  if (RM::ValidateCompiledSpeciesBinding(options, 2, invalid).ok())
    return Fail("a zero-mass compiled species was accepted");
  options.observers.front().allCompiledSpecies = false;
  options.observers.front().species = {2};
  if (RM::ValidateCompiledSpeciesBinding(options, 2, compiled).ok())
    return Fail("an observer species outside the compiled table was accepted");

  options = RM::RunConfiguration3DOptions();
  std::shared_ptr<const RM::RunConfiguration3D> baseline, changed;
  if (!RM::RunConfiguration3D::Create(options, &baseline).ok())
    return Fail("baseline proton configuration is invalid");
  options.species.macroparticleWeight = 2.0;
  if (!RM::RunConfiguration3D::Create(options, &changed).ok() ||
      baseline->physics_fingerprint() == changed->physics_fingerprint())
    return Fail("species weight is absent from restart compatibility identity");
  return Pass("the complete immutable AMPS species table binds without proton or slot assumptions");
}

Result RunCFG3D08() {
  RM::RunConfiguration3DOptions options;
  const SEP3D::Core::Status loaded = RM::LoadConfigurationFile(
      "examples/sep3d_analytic_parker.in", &options);
  if (!loaded.ok())
    return Fail("complete schema-version-4 initialization input did not resolve: " +
                loaded.message);
  if (
      options.inputSchemaVersion != 4 ||
      options.injectionCadenceSteps != 1 ||
      options.source.samplesPerStep != 1000 ||
      options.swcmeAssignments.empty() ||
      options.swcmeConfigurationFingerprint.empty() ||
      options.swcmeResolvedManifest.empty()) {
    return Fail("schema-version-4 initialization metadata is incomplete");
  }

  std::ifstream input("examples/sep3d_analytic_parker.in");
  std::ostringstream buffer;
  buffer << input.rdbuf();
  const std::string complete = buffer.str();
  if (!input || complete.empty()) return Fail("could not read schema-v3 fixture");

  std::string missing = complete;
  const std::string required = "surface.phi_points = 24\n";
  const std::size_t requiredAt = missing.find(required);
  if (requiredAt == std::string::npos)
    return Fail("schema-v3 fixture lost its surface resolution field");
  missing.erase(requiredAt, required.size());
  if (RM::ParseConfigurationText(missing, &options).ok())
    return Fail("schema version 3 accepted an incomplete canonical SWCME model");

  std::string wrongWeight = complete;
  const std::string weight = "macroparticle_weight = 1e23";
  const std::size_t weightAt = wrongWeight.find(weight);
  if (weightAt == std::string::npos)
    return Fail("schema-v3 fixture lost its derived particle weight");
  wrongWeight.replace(weightAt, weight.size(),
                      "macroparticle_weight = 2e23");
  if (RM::ParseConfigurationText(wrongWeight, &options).ok())
    return Fail("schema version 3 accepted a particle weight inconsistent with its physical source rate");

  std::string skippedStep = complete;
  const std::string cadence = "injection_cadence_steps = 1";
  const std::size_t cadenceAt = skippedStep.find(cadence);
  if (cadenceAt == std::string::npos)
    return Fail("schema-v3 fixture lost its injection cadence");
  skippedStep.replace(cadenceAt, cadence.size(),
                      "injection_cadence_steps = 2");
  if (RM::ParseConfigurationText(skippedStep, &options).ok())
    return Fail("schema version 3 accepted a source that skips simulation steps");

  std::string absoluteNormalization = complete;
  const std::string relative = "source.normalization = relative_only";
  const std::size_t relativeAt = absoluteNormalization.find(relative);
  if (relativeAt == std::string::npos)
    return Fail("schema-v3 fixture lost its source normalization");
  absoluteNormalization.replace(
      relativeAt, relative.size(),
      "source.normalization = reference_differential_intensity");
  if (RM::ParseConfigurationText(absoluteNormalization, &options).ok())
    return Fail("schema version 3 accepted two competing absolute source normalizations");

  std::string wrongAxis = complete;
  const std::string axisX = "geometry.solar_rotation_axis_x = 0";
  const std::string axisZ = "geometry.solar_rotation_axis_z = 1";
  const std::size_t axisXAt = wrongAxis.find(axisX);
  const std::size_t axisZAt = wrongAxis.find(axisZ);
  if (axisXAt == std::string::npos || axisZAt == std::string::npos)
    return Fail("schema-v3 fixture lost its solar rotation axis");
  wrongAxis.replace(axisXAt, axisX.size(),
                    "geometry.solar_rotation_axis_x = 1");
  wrongAxis.replace(wrongAxis.find(axisZ), axisZ.size(),
                    "geometry.solar_rotation_axis_z = 0");
  if (RM::ParseConfigurationText(wrongAxis, &options).ok())
    return Fail("schema version 3 accepted inconsistent Parker rotation axes");

  // SWCME defines total |B| at one AU, while AnalyticParkerProvider accepts
  // radial Br at its independently declared reference radius.  A half-AU
  // reference therefore has four times the one-AU radial component.  This
  // positive case prevents the cross-model consistency gate from confusing
  // those two distinct physical conventions.
  std::string halfAuReference = complete;
  const std::size_t backgroundAt =
      halfAuReference.find("[background.parker]");
  const std::string oneAuReference =
      "reference_radius_m = 1.495978707e11";
  const std::size_t referenceAt =
      halfAuReference.find(oneAuReference, backgroundAt);
  const std::string oneAuRadialField =
      "radial_field_at_reference_t = 3.585667175612118e-9";
  const std::size_t fieldAt = halfAuReference.find(oneAuRadialField);
  if (backgroundAt == std::string::npos ||
      referenceAt == std::string::npos || fieldAt == std::string::npos)
    return Fail("schema-v3 fixture lost its Parker reference normalization");
  halfAuReference.replace(referenceAt, oneAuReference.size(),
                          "reference_radius_m = 7.479893535e10");
  halfAuReference.replace(
      fieldAt, oneAuRadialField.size(),
      "radial_field_at_reference_t = 1.4342668702448471e-8");
  if (!RM::ParseConfigurationText(halfAuReference, &options).ok())
    return Fail("schema version 3 rejected the correctly r^-2-scaled half-AU Parker reference");

  // The source surface is explicit physics, not a hidden 20-R_sun constant.
  // Move the domain, finite line, Parker source, and CME launch surface to
  // 25 R_sun and supply the corresponding one-AU radial component.
  std::string movedSource = complete;
  auto replaceRequired = [&](const std::string& from,
                             const std::string& to) -> bool {
    const std::size_t at = movedSource.find(from);
    if (at == std::string::npos) return false;
    movedSource.replace(at, from.size(), to);
    return true;
  };
  if (!replaceRequired("inner_radius_m = 1.3914e10",
                       "inner_radius_m = 1.73925e10") ||
      !replaceRequired("initial_x_m = 1.3914e10",
                       "initial_x_m = 1.73925e10") ||
      !replaceRequired(oneAuRadialField,
                       "radial_field_at_reference_t = 3.6305743892405979e-9") ||
      !replaceRequired("parker.source_radius = 20 Rs",
                       "parker.source_radius = 25 Rs") ||
      !replaceRequired("cme.launch_radius = 20 Rs",
                       "cme.launch_radius = 25 Rs"))
    return Fail("schema-v3 fixture lost an explicit source-surface field");
  if (!RM::ParseConfigurationText(movedSource, &options).ok())
    return Fail("schema version 3 retained a hidden 20-R_sun source assumption");

  return Pass("schema version 4 retains the complete schema-3 SWCME physics contract, honors explicit Parker references/source radii, and rejects missing or inconsistent physics");
}

Result RunCFG3D09() {
  std::ifstream input("examples/sep3d_analytic_parker.in");
  std::ostringstream buffer;
  buffer << input.rdbuf();
  const std::string complete = buffer.str();
  if (!input || complete.empty())
    return Fail("could not read the selectable-background/turbulence fixture");

  RM::RunConfiguration3DOptions options;
  std::string kraichnan = complete;
  const std::string model = "model = kolmogorov";
  const std::string index = "spectral_index = 1.6666666666666667";
  const std::size_t modelAt = kraichnan.find(model);
  const std::size_t indexAt = kraichnan.find(index);
  if (modelAt == std::string::npos || indexAt == std::string::npos)
    return Fail("example lost its explicit turbulence model or slope");
  kraichnan.replace(modelAt, model.size(), "model = kraichnan");
  kraichnan.replace(kraichnan.find(index), index.size(),
                    "spectral_index = 1.5");
  if (!RM::ParseConfigurationText(kraichnan, &options).ok() ||
      options.prescribedTurbulenceModel !=
          RM::PrescribedTurbulenceModel::Kraichnan) {
    return Fail("a physically consistent Kraichnan selection did not parse");
  }

  std::string contradictory = complete;
  contradictory.replace(contradictory.find(model), model.size(),
                         "model = kraichnan");
  if (RM::ParseConfigurationText(contradictory, &options).ok())
    return Fail("a named Kraichnan model accepted the Kolmogorov index");

  // The amplitude prescription is a separate input choice from the spectral
  // slope.  Select direct total wave-energy normalization and its radial law,
  // while setting the inactive deltaB/B normalization to the required zero.
  std::string waveEnergy = complete;
  auto replaceAmplitude = [&](const std::string& from,
                              const std::string& to) -> bool {
    const std::size_t at = waveEnergy.find(from);
    if (at == std::string::npos) return false;
    waveEnergy.replace(at, from.size(), to);
    return true;
  };
  if (!replaceAmplitude("amplitude_model = constant-delta-b-over-b",
                        "amplitude_model = wave-energy-power-law") ||
      !replaceAmplitude("delta_b_over_b = 0.3", "delta_b_over_b = 0") ||
      !replaceAmplitude("wave_energy_density_at_reference_j_per_m3 = 0",
                        "wave_energy_density_at_reference_j_per_m3 = 4e-12") ||
      !replaceAmplitude("wave_energy_density_radial_exponent = 0",
                        "wave_energy_density_radial_exponent = 2")) {
    return Fail("example lost its explicit turbulence amplitude controls");
  }
  if (!RM::ParseConfigurationText(waveEnergy, &options).ok() ||
      options.prescribedTurbulenceAmplitudeModel !=
          RM::PrescribedTurbulenceAmplitudeModel::WaveEnergyPowerLaw ||
      options.turbulenceWaveEnergyAtReferenceJPerM3 != 4.0e-12 ||
      options.turbulenceWaveEnergyRadialExponent != 2.0) {
    return Fail("wave-energy-power-law amplitude selection did not parse exactly");
  }

  std::string competingAmplitudes = waveEnergy;
  const std::string inactiveRatio = "delta_b_over_b = 0";
  competingAmplitudes.replace(competingAmplitudes.find(inactiveRatio),
                              inactiveRatio.size(),
                              "delta_b_over_b = 0.3");
  if (RM::ParseConfigurationText(competingAmplitudes, &options).ok())
    return Fail("input accepted two competing turbulence amplitude normalizations");

  std::string python = complete;
  const std::string parker = "provider = analytic-parker";
  const std::size_t parkerAt = python.find(parker);
  if (parkerAt == std::string::npos)
    return Fail("example lost its explicit background provider");
  python.replace(parkerAt, parker.size(), "provider = python-interpolator");
  const SEP3D::Core::Status pythonStatus =
      RM::ParseConfigurationText(python, &options);
  if (pythonStatus.code != SEP3D::Core::StatusCode::ReservedFeature ||
      pythonStatus.message.find("Python") == std::string::npos) {
    return Fail("reserved Python background did not fail with its typed status");
  }

  return Pass(
      "input selects validated spectral/amplitude turbulence models and reserves the future Python background explicitly");
}

Result RunCFG3D10() {
  std::ifstream input("examples/sep3d_analytic_parker.in");
  std::ostringstream buffer;
  buffer << input.rdbuf();
  const std::string complete = buffer.str();
  if (!input || complete.empty())
    return Fail("could not read the CME/Parker linkage fixture");

  RM::RunConfiguration3DOptions linked;
  const SEP3D::Core::Status parsed =
      RM::ParseConfigurationText(complete, &linked);
  if (!parsed.ok() ||
      linked.parkerSpiralStartMode !=
          RM::ParkerSpiralStartMode::CmeLaunchPoint ||
      !linked.cmeLaunchPointResolved ||
      (linked.parkerSpiralInitialPointM - linked.cmeLaunchPointM).Norm() !=
          0.0) {
    return Fail("canonical SWCME launch apex was not frozen as the Parker start");
  }

  // A changed CME direction must not leave the refinement tube/line at the
  // former coordinates.  The linked mode rejects this incomplete edit before
  // RunConfiguration reaches AMPS allocation.
  std::string wrongDirection = complete;
  const std::string yDirection = "geometry.cme_direction_y = 0";
  const std::size_t yAt = wrongDirection.find(yDirection);
  if (yAt == std::string::npos)
    return Fail("linkage fixture lost the CME direction");
  wrongDirection.replace(yAt, yDirection.size(),
                         "geometry.cme_direction_y = 1");
  if (RM::ParseConfigurationText(wrongDirection, &linked).ok())
    return Fail("CME-linked Parker start accepted a different CME direction");

  std::string wrongRadius = complete;
  const std::string launch = "cme.launch_radius = 20 Rs";
  const std::size_t launchAt = wrongRadius.find(launch);
  if (launchAt == std::string::npos)
    return Fail("linkage fixture lost the CME launch radius");
  wrongRadius.replace(launchAt, launch.size(),
                      "cme.launch_radius = 25 Rs");
  if (RM::ParseConfigurationText(wrongRadius, &linked).ok())
    return Fail("CME-linked Parker start accepted a launch radius outside the inner sphere");

  // Explicit mode remains available for a deliberately independent field
  // line.  It must not carry a stale derived CME point into the fingerprint.
  std::string explicitLine = complete;
  const std::string linkedMode = "start_mode = cme-launch-point";
  const std::size_t modeAt = explicitLine.find(linkedMode);
  if (modeAt == std::string::npos)
    return Fail("linkage fixture lost its start mode");
  explicitLine.replace(modeAt, linkedMode.size(), "start_mode = explicit");
  if (!RM::ParseConfigurationText(explicitLine, &linked).ok() ||
      linked.cmeLaunchPointResolved ||
      linked.parkerSpiralStartMode != RM::ParkerSpiralStartMode::Explicit) {
    return Fail("explicit Parker start retained an implicit CME linkage");
  }
  return Pass(
      "input explicitly links the Parker start to the canonical SWCME launch apex and rejects radius/direction mismatches");
}

Result RunCFG3D11() {
  std::ifstream input("examples/sep3d_analytic_parker.in");
  std::ostringstream buffer;
  buffer << input.rdbuf();
  const std::string complete = buffer.str();
  if (!input || complete.empty())
    return Fail("could not read the schema-4 transport/control fixture");

  auto replaceOnce = [](std::string* text, const std::string& from,
                        const std::string& to) -> bool {
    if (text == nullptr) return false;
    const std::size_t at = text->find(from);
    if (at == std::string::npos) return false;
    text->replace(at, from.size(), to);
    return true;
  };

  RM::RunConfiguration3DOptions options;
  if (!RM::ParseConfigurationText(complete, &options).ok() ||
      options.inputSchemaVersion != 4 ||
      options.transport != RM::TransportModel::Parker3D ||
      options.activeRegion != RM::ActiveRegionMode::FullDomain ||
      options.populationControl != RM::PopulationControlMode::Off ||
      options.source.spectrumModel !=
          RM::SourceSpectrumModel::LocalCompressionDsa ||
      options.source.fixedPhaseSpacePowerIndex != 0.0 ||
      options.spatialDiffusionModel !=
          RM::SpatialDiffusionModel::CorrelationMeanFreePath ||
      options.pitchAngleDiffusionModel !=
          RM::PitchAngleDiffusionModel::Jokipii1966 ||
      options.meanFreePathModel != RM::MeanFreePathModel::Correlation) {
    return Fail("canonical schema-4 selectors did not parse exactly");
  }

  // Parse the shipped active-tube example itself, rather than proving only
  // that a string synthesized by this test would be accepted. This prevents
  // documentation/example drift from reintroducing the old observer
  // coordinates, which followed a curve different from the initialized IMF.
  std::ifstream activeInput(
      "examples/sep3d_analytic_parker_active_tube.in");
  std::ostringstream activeBuffer;
  activeBuffer << activeInput.rdbuf();
  RM::RunConfiguration3DOptions activeOptions;
  const SEP3D::Core::Status activeStatus = RM::ParseConfigurationText(
      activeBuffer.str(), &activeOptions);
  const double observerToleranceM = 1.0;
  if (!activeInput || !activeStatus.ok() ||
      activeOptions.activeRegion != RM::ActiveRegionMode::ParkerTube ||
      activeOptions.activeTubeBufferBlocks != 1 ||
      activeOptions.observers.size() != 2 ||
      activeOptions.observers[0].id != "earth" ||
      activeOptions.observers[1].id != "inner" ||
      std::fabs(activeOptions.observers[0].positionM.x -
                1.1096225266380298e11) > observerToleranceM ||
      std::fabs(activeOptions.observers[0].positionM.y +
                1.0033394939774008e11) > observerToleranceM ||
      std::fabs(activeOptions.observers[1].positionM.x -
                4.463181234418467e10) > observerToleranceM ||
      std::fabs(activeOptions.observers[1].positionM.y +
                4.70726985535542e9) > observerToleranceM) {
    return Fail("shipped active-tube example is not a valid field-connected schema-4 deck");
  }
  SEP3D::Core::ParkerSpiralGeometry activeGeometry;
  activeGeometry.sourceRadiusM = activeOptions.innerRadiusM;
  activeGeometry.sourceLongitudeRad = activeOptions.tubeLongitudeRad;
  activeGeometry.sourceColatitudeRad = activeOptions.tubeColatitudeRad;
  activeGeometry.solarWindSpeedMPerS =
      activeOptions.parker.solarWindSpeedMPerS;
  activeGeometry.solarRotationRateRadPerS =
      activeOptions.parker.solarRotationRateRadPerS;
  activeGeometry.rotationAxis = activeOptions.parker.rotationAxis;
  for (const RM::ObserverOptions& observer : activeOptions.observers) {
    const SEP3D::Core::Vec3 relative =
        observer.positionM - activeOptions.coordinateOriginM;
    const SEP3D::Core::Vec3 expected = activeOptions.coordinateOriginM +
        SEP3D::Core::ParkerCurvePoint(relative.Norm(), activeGeometry);
    if ((observer.positionM - expected).Norm() > observerToleranceM ||
        std::fabs(relative.Norm() - observer.shellRadiusM) >
            observerToleranceM) {
      return Fail("shipped active-tube observer is not on the exact finite Parker line");
    }
  }

  // A fixed phase-space law names q in f(p) proportional to p^(-q).  The
  // input therefore uses positive q=5 for the published p^-5 seed shape; the
  // sampling adapter, not the user, performs the dN/dp exponent conversion.
  std::string fixedSource = complete;
  if (!replaceOnce(&fixedSource,
                   "spectrum_model = local-compression-dsa",
                   "spectrum_model = fixed-phase-space-power-law") ||
      !replaceOnce(&fixedSource, "phase_space_power_index = 0",
                   "phase_space_power_index = 5") ||
      !RM::ParseConfigurationText(fixedSource, &options).ok() ||
      options.source.spectrumModel !=
          RM::SourceSpectrumModel::FixedPhaseSpacePowerLaw ||
      options.source.fixedPhaseSpacePowerIndex != 5.0) {
    return Fail("fixed p^-5 phase-space source did not parse exactly");
  }
  std::string signedIndex = fixedSource;
  if (!replaceOnce(&signedIndex, "phase_space_power_index = 5",
                   "phase_space_power_index = -5") ||
      RM::ParseConfigurationText(signedIndex, &options).ok()) {
    return Fail("fixed source accepted signed -5 instead of positive q=5");
  }
  std::string dormantIndex = complete;
  if (!replaceOnce(&dormantIndex, "phase_space_power_index = 0",
                   "phase_space_power_index = 5") ||
      RM::ParseConfigurationText(dormantIndex, &options).ok()) {
    return Fail("local DSA source accepted a dormant fixed phase-space index");
  }
  std::string missingSourceModel = complete;
  const std::string modelLine =
      "spectrum_model = local-compression-dsa\n";
  const std::size_t modelLineAt = missingSourceModel.find(modelLine);
  if (modelLineAt == std::string::npos)
    return Fail("schema-4 fixture lost its explicit source spectrum model");
  missingSourceModel.erase(modelLineAt, modelLine.size());
  if (RM::ParseConfigurationText(missingSourceModel, &options).ok())
    return Fail("schema version 4 accepted a missing source spectrum model");

  // A useful active corridor must contain the separately configured refined
  // tube and must retain a complete-block stencil halo.  This positive case
  // uses a wider one-AU radius than the 0.03-AU refinement tube.
  std::string corridor = complete;
  if (!replaceOnce(&corridor, "mode = full-domain", "mode = parker-tube") ||
      !replaceOnce(&corridor, "radius_at_reference_m = 0\nradius_mode = constant-angular-width\nbuffer_blocks = 0",
                   "radius_at_reference_m = 7.479893535e9\nradius_mode = constant-angular-width\nbuffer_blocks = 1") ||
      !replaceOnce(&corridor,
                   "position_x_m = 1.495978707e11\nposition_y_m = 0",
                   "position_x_m = 1.10962252663803e11\nposition_y_m = -1.00333949397740e11") ||
      !replaceOnce(&corridor,
                   "position_x_m = 4.487936121e10\nposition_y_m = 0",
                   "position_x_m = 4.46318123441847e10\nposition_y_m = -4.70726985535542e9") ||
      !RM::ParseConfigurationText(corridor, &options).ok() ||
      options.activeRegion != RM::ActiveRegionMode::ParkerTube ||
      options.activeTubeBufferBlocks != 1) {
    return Fail("valid Parker active corridor did not parse");
  }
  std::string narrow = corridor;
  if (!replaceOnce(&narrow, "radius_at_reference_m = 7.479893535e9",
                   "radius_at_reference_m = 1e9") ||
      RM::ParseConfigurationText(narrow, &options).ok()) {
    return Fail("active corridor narrower than the refinement tube was accepted");
  }
  std::string mixedRadiusModes = corridor;
  if (!replaceOnce(&mixedRadiusModes,
                   "radius_mode = constant-angular-width",
                   "radius_mode = physical-constant") ||
      RM::ParseConfigurationText(mixedRadiusModes, &options).ok()) {
    return Fail("active corridor accepted endpoint under-coverage from mixed radius modes");
  }
  std::string disconnectedObserver = corridor;
  if (!replaceOnce(&disconnectedObserver,
                   "position_x_m = 1.10962252663803e11\nposition_y_m = -1.00333949397740e11",
                   "position_x_m = 1.495978707e11\nposition_y_m = 0") ||
      RM::ParseConfigurationText(disconnectedObserver, &options).ok()) {
    return Fail("observer outside the stationary active corridor was accepted");
  }
  std::string truncatedBeforeObserver = corridor;
  if (!replaceOnce(&truncatedBeforeObserver, "length_m = 2.0e11",
                   "length_m = 1.0e10") ||
      RM::ParseConfigurationText(truncatedBeforeObserver, &options).ok()) {
    return Fail("observer beyond the finite active-corridor end cap was accepted");
  }

  std::string controlled = complete;
  if (!replaceOnce(&controlled, "mode = off\nminimum_particles_per_cell_per_species = 0\ntarget_particles_per_cell_per_species = 0\nmaximum_particles_per_cell_per_species = 0",
                   "mode = split-merge\nminimum_particles_per_cell_per_species = 16\ntarget_particles_per_cell_per_species = 24\nmaximum_particles_per_cell_per_species = 32") ||
      !RM::ParseConfigurationText(controlled, &options).ok() ||
      options.populationControl != RM::PopulationControlMode::SplitMerge ||
      options.minimumParticlesPerCellPerSpecies != 16 ||
      options.targetParticlesPerCellPerSpecies != 24 ||
      options.maximumParticlesPerCellPerSpecies != 32) {
    return Fail("valid split/merge population limits did not parse exactly");
  }
  std::string inverted = controlled;
  if (!replaceOnce(&inverted, "maximum_particles_per_cell_per_species = 32",
                   "maximum_particles_per_cell_per_species = 20") ||
      RM::ParseConfigurationText(inverted, &options).ok()) {
    return Fail("population limits outside minimum <= target <= maximum were accepted");
  }

  std::string events = complete;
  if (!replaceOnce(&events, "transport = parker",
                   "transport = focused-scattering") ||
      !RM::ParseConfigurationText(events, &options).ok() ||
      options.transport != RM::TransportModel::FocusedScattering3D) {
    return Fail("focused-scattering mover with a correlation MFP did not parse");
  }
  std::string invalidTransverseEvents = events;
  if (!replaceOnce(&invalidTransverseEvents, "perpendicular_diffusion = none",
                   "perpendicular_diffusion = constant") ||
      !replaceOnce(&invalidTransverseEvents,
                   "constant_kappa_perpendicular_m2_per_s = 0",
                   "constant_kappa_perpendicular_m2_per_s = 1e16") ||
      RM::ParseConfigurationText(invalidTransverseEvents, &options).ok()) {
    return Fail("unvalidated perpendicular diffusion with discrete scattering was accepted");
  }

  // The event-validation law is species-general in rigidity. All reference
  // values and exponents are explicit SI inputs; no proton-only energy
  // shortcut or hidden 1-AU/1-GV normalization is permitted by the parser.
  std::string powerLaw = complete;
  if (!replaceOnce(&powerLaw, "mean_free_path_model = correlation",
                   "mean_free_path_model = radial-rigidity-power-law") ||
      !replaceOnce(&powerLaw, "mean_free_path_reference_m = 0",
                   "mean_free_path_reference_m = 4.487936121e10") ||
      !replaceOnce(&powerLaw, "mean_free_path_reference_radius_m = 0",
                   "mean_free_path_reference_radius_m = 1.495978707e11") ||
      !replaceOnce(&powerLaw, "mean_free_path_reference_rigidity_v = 0",
                   "mean_free_path_reference_rigidity_v = 1e9") ||
      !replaceOnce(&powerLaw, "mean_free_path_radial_exponent = 0",
                   "mean_free_path_radial_exponent = 1") ||
      !replaceOnce(&powerLaw, "mean_free_path_rigidity_exponent = 0",
                   "mean_free_path_rigidity_exponent = 0.3333333333333333") ||
      !RM::ParseConfigurationText(powerLaw, &options).ok() ||
      options.meanFreePathModel !=
          RM::MeanFreePathModel::RadialRigidityPowerLaw ||
      options.meanFreePathReferenceRigidityV != 1.0e9) {
    return Fail("radial-rigidity power-law MFP did not parse exactly");
  }
  std::string invalidPowerLaw = powerLaw;
  if (!replaceOnce(&invalidPowerLaw,
                   "mean_free_path_reference_rigidity_v = 1e9",
                   "mean_free_path_reference_rigidity_v = 0") ||
      RM::ParseConfigurationText(invalidPowerLaw, &options).ok()) {
    return Fail("power-law MFP accepted a non-positive reference rigidity");
  }

  std::shared_ptr<const RM::RunConfiguration3D> configuration;
  if (!RM::ParseConfigurationText(controlled, &options).ok() ||
      !RM::RunConfiguration3D::Create(options, &configuration).ok()) {
    return Fail("valid schema-4 fixture could not be frozen");
  }
  std::string summary;
  if (!RM::BuildDryRunSummary(*configuration, &summary).ok() ||
      summary.find("active_region_mode=full-domain") == std::string::npos ||
      summary.find("population_control_mode=split-merge") ==
          std::string::npos ||
      summary.find("transport_model=parker") == std::string::npos ||
      summary.find("spatial_diffusion_model=mean-free-path") ==
          std::string::npos ||
      summary.find("source_spectrum_model=local-compression-dsa") ==
          std::string::npos ||
      summary.find("source_phase_space_power_index=0") ==
          std::string::npos) {
    return Fail("dry-run omitted schema-4 active-region or transport controls");
  }

  return Pass(
      "schema-4 corridor, population hysteresis, fixed/DSA source spectra, mover/coefficient compatibility including radial-rigidity MFP, and dry-run contracts passed");
}

Result RunCFG3D12() {
  std::string input = CompleteInput();
  const std::string old = "outer_radius_mode = preset";
  const auto at = input.find(old);
  if (at == std::string::npos) return Fail("missing domain fixture");
  input.replace(at, old.size(),
      "outer_radius_mode = field-line-endpoint\n"
      "box_geometry = field-line-corner-cube\n"
      "corner_direction_x = 0\ncorner_direction_y = 0\ncorner_direction_z = 0\n"
      "corner_margin_m = 1.495978707e9");
  input += "\n[parker_spiral]\nend_radius_m = 1.495978707e11\nlength_m = 0\n"
      "[mesh.active_region]\nmode = parker-tube\nreference_radius_m = 1.495978707e11\n"
      "radius_at_reference_m = 1.1967829656e11\nradius_mode = constant-angular-width\n"
      "buffer_blocks = 1\nsolar_sphere_radius_m = 2.0871e10\n";
  // Add the new anchor inside the existing solar section, without creating
  // a duplicate section that the strict INI grammar must reject.
  const auto solar = input.find("[mesh.solar]");
  if (solar == std::string::npos) return Fail("missing solar fixture");
  input.insert(solar + std::string("[mesh.solar]").size(), "\nanchor = photosphere");
  RM::RunConfiguration3DOptions parsed;
  const auto status = RM::ParseConfigurationText(input, &parsed);
  if (!status.ok()) return Fail("corner input rejected: " + status.message);
  SEP3D::Core::ParkerSpiralGeometry geometry;
  geometry.sourceRadiusM = parsed.innerRadiusM;
  geometry.sourceLongitudeRad = parsed.tubeLongitudeRad;
  geometry.sourceColatitudeRad = parsed.tubeColatitudeRad;
  geometry.solarWindSpeedMPerS = parsed.parker.solarWindSpeedMPerS;
  geometry.solarRotationRateRadPerS = parsed.parker.solarRotationRateRadPerS;
  const double arc = SEP3D::Core::ParkerCurveArcLengthM(parsed.parkerSpiralEndRadiusM, geometry);
  if (parsed.outerRadiusM != parsed.parkerSpiralEndRadiusM ||
      parsed.parkerSpiralLengthM != arc ||
      parsed.solarRefinementAnchor != RM::SolarRefinementAnchor::Photosphere)
    return Fail("endpoint radius did not normalize the physical extent/arc length/anchor");
  std::shared_ptr<const RM::RunConfiguration3D> frozen, changed;
  if (!RM::RunConfiguration3D::Create(parsed, &frozen).ok())
    return Fail("normalized corner options cannot round trip through the factory");
  auto edited = parsed;
  edited.activeSolarSphereRadiusM += SEP3D::Core::Const::R_sun;
  if (!RM::RunConfiguration3D::Create(edited, &changed).ok() ||
      frozen->physics_fingerprint() == changed->physics_fingerprint())
    return Fail("solar sphere is absent from immutable physics identity");
  for (int bad = 0; bad < 4; ++bad) {
    edited = parsed;
    if (bad == 0) edited.parkerSpiralLengthM = 1.0;
    if (bad == 1) edited.activeSolarSphereRadiusM = 0.5 * parsed.innerRadiusM;
    if (bad == 2) edited.domainCornerDirection.x = 2.0;
    if (bad == 3) edited.maximumMeshLevel = 0;
    if (RM::RunConfiguration3D::Create(edited, &changed).ok())
      return Fail("invalid/conflicting corner geometry was accepted");
  }
  std::string report;
  if (!RM::BuildDryRunSummary(*frozen, &report).ok() ||
      report.find("domain_box_geometry=field-line-corner-cube") == std::string::npos ||
      report.find("solar_refinement_anchor=photosphere") == std::string::npos)
    return Fail("dry run does not expose the resolved corner/sphere geometry");
  // The x-y variant centers z on the Sun and changes the immutable geometry
  // identity. Keep the original three-axis mode covered above for old decks.
  std::string xyInput = input;
  const std::string mode = "field-line-corner-cube";
  xyInput.replace(xyInput.find(mode), mode.size(), "field-line-xy-corner-cube");
  RM::RunConfiguration3DOptions xy;
  const auto xyStatus = RM::ParseConfigurationText(xyInput, &xy);
  if (!xyStatus.ok()) return Fail("x-y corner input rejected: " + xyStatus.message);
  if (xy.domainBoxGeometry != RM::DomainBoxGeometry::FieldLineXYCornerCube ||
      !RM::RunConfiguration3D::Create(xy, &changed).ok() ||
      frozen->physics_fingerprint() == changed->physics_fingerprint())
    return Fail("centered-z mode was not parsed or fingerprinted independently");
  SEP3D::Core::Vec3 low, high;
  if (!RM::ResolveDomainBoundsM(xy, &low, &high).ok() ||
      std::fabs(low.z + high.z) > 1e-12 * (high.z - low.z))
    return Fail("x-y input did not produce symmetric z bounds");
  xy.domainCornerDirection.z = 1.0;
  if (RM::RunConfiguration3D::Create(xy, &changed).ok())
    return Fail("x-y corner mode accepted a conflicting z corner direction");
  return Pass("both corner modes parse/fingerprint, endpoint extent normalizes, z is centered in x-y mode, and conflicting controls fail closed");
}

Result RunCFG3D13() {
  namespace fs = std::filesystem;
  const fs::path directory = "test_output/cfg3d13-application-input";
  std::error_code cleanupError;
  fs::remove_all(directory, cleanupError);
  cleanupError.clear();
  fs::create_directories(directory / "parts", cleanupError);
  if (cleanupError) return Fail("could not create application-input fixture directory");

  const auto write = [](const fs::path& path, const std::string& text) {
    std::ofstream stream(path);
    stream << text;
    stream.close();
    return static_cast<bool>(stream);
  };
  if (!write(directory / "amps.in",
             "! a different application is intentionally ignored\n"
             "#section begin: core\n"
             "unrelated = value\n"
             "#section end\n"
             "#include \"parts/sep3d.in\"\n") ||
      !write(directory / "parts/sep3d.in",
             "#section begin: sep3d\n"
             "shock_model = reduced-shock-surface\n"
             "background_plasma_model = corona-swcme-ambient\n"
             "source_model = accepted-shock-incident-flux\n"
             "maximum_time_steps = 9\n"
             "particles_per_iteration = \\\n"
             "  37 ! joined to the prior physical line\n"
             "maximum_particle_speed_m_s = 200000000\n"
             "time_step_margin_factor = 0.25\n"
             "source_normalization_radius_m = 13914000000\n"
             "#subsection begin: shock-particle-injection\n"
             "statistical_weight_model=constant-statistical-weight\n"
             "minimum_energy_j=1.602176634e-15\n"
             "maximum_energy_j=1.602176634e-11\n"
             "phase_space_power_model=compression-ratio\n"
             "maximum_events_per_species_per_step=1000000\n"
             "#subsection end\n"
             "#subsection begin: reduced-shock-surface\n"
             "schema = parser-fixture\n"
             "#subsection end\n"
             "#section end\n")) {
    return Fail("could not write application-input fixtures");
  }

  RM::Sep3dApplicationInput parsed;
  SEP3D::Core::Status status = RM::ParseSep3dApplicationInput(
      (directory / "amps.in").string(), &parsed);
  if (!status.ok() || parsed.maximumTimeSteps != 9 ||
      parsed.particlesPerIteration != 37 ||
      parsed.expandedFiles.size() != 2 || parsed.valueLine != 6 ||
      parsed.shockModel != "reduced-shock-surface" ||
      parsed.backgroundPlasmaModel != "corona-swcme-ambient" ||
      parsed.sourceModel != "accepted-shock-incident-flux" ||
      parsed.maximumParticleSpeedMPerS != 2.0e8 ||
      parsed.timeStepMarginFactor != 0.25 ||
      parsed.sourceNormalizationRadiusM != 13914000000.0 ||
      parsed.particleWeightingModel != "constant-statistical-weight" ||
      parsed.momentumPowerLawModel != "compression-ratio" ||
      parsed.minimumInjectionEnergyJ != 1.602176634e-15 ||
      parsed.maximumInjectionEnergyJ != 1.602176634e-11 ||
      parsed.maximumInjectionEventsPerSpeciesPerStep != 1000000 ||
      parsed.reducedShockConfiguration != "schema=parser-fixture\n" ||
      RM::Sep3dApplicationInputSummary(parsed).find(
          "particles_per_iteration=37") == std::string::npos) {
    return Fail("recursive include/comment/continuation input did not resolve: " +
                status.message);
  }

  // The new CLI spelling and no-argument default select shared-section mode;
  // the maintained double-dash spelling remains the complete schema deck.
  RM::StandaloneCommandLine cli;
  const char* defaultArgv[] = {"amps"};
  if (!RM::ParseStandaloneCommandLine(
          1, const_cast<char**>(defaultArgv), &cli).ok() ||
      cli.inputPath != "amps.in" || !cli.sectionInput) {
    return Fail("no-argument CLI did not select ./amps.in section mode");
  }
  const char* sectionArgv[] = {"amps", "-input", "global.in"};
  if (!RM::ParseStandaloneCommandLine(
          3, const_cast<char**>(sectionArgv), &cli).ok() ||
      cli.inputPath != "global.in" || !cli.sectionInput) {
    return Fail("-input did not select shared-section mode");
  }
  const char* legacyArgv[] = {"amps", "--input", "schema.in"};
  if (!RM::ParseStandaloneCommandLine(
          3, const_cast<char**>(legacyArgv), &cli).ok() || cli.sectionInput) {
    return Fail("maintained --input schema mode regressed");
  }

  // Commit the parsed value through the same immutable pre-mesh replacement
  // used by production.  A changed fingerprint proves the setting is not an
  // unfingerprinted mutable global, and layout equality proves parser timing
  // cannot invalidate already registered AMPS byte requests.
  RM::RunConfiguration3DOptions options;
  options.meshMemoryBudgetBytes = 1000000000000000ULL;
  std::shared_ptr<const RM::RunConfiguration3D> provisional, resolved;
  if (!RM::RunConfiguration3D::Create(options, &provisional).ok())
    return Fail("could not create provisional configuration");
  options.inputSchemaVersion = 4;
  options.maximumTimeSteps = parsed.maximumTimeSteps;
  options.background = RM::BackgroundAuthority::RuntimeModel;
  options.backgroundModelId = "sep-corona-swcme-shock-front-v1";
  options.backgroundModelInlineConfiguration =
      parsed.reducedShockConfiguration;
  options.backgroundModelAssetDirectory = parsed.reducedShockAssetDirectory;
  options.coordinateFrame = "HCI";
  options.parker.coordinateFrame = options.coordinateFrame;
  options.intent=RM::RunIntent::ShockInjection;
  options.source.enabled=true;
  options.source.samplesPerStep = parsed.particlesPerIteration;
  options.source.minimumEnergyJ=parsed.minimumInjectionEnergyJ;
  options.source.maximumEnergyJ=parsed.maximumInjectionEnergyJ;
  options.source.spectrumModel=RM::SourceSpectrumModel::LocalCompressionDsa;
  options.source.fixedPhaseSpacePowerIndex=0;
  options.source.weightingModel=
      RM::SourceWeightingModel::ConstantStatisticalWeight;
  options.source.maximumMacroparticlesPerSpeciesPerStep=
      parsed.maximumInjectionEventsPerSpeciesPerStep;
  options.particleNumerics.deriveFromMeshAndShock = true;
  options.particleNumerics.sourceRateModel =
      RM::SourceRateNormalizationModel::AcceptedShockIncidentFlux;
  options.particleNumerics.maximumParticleSpeedMPerS =
      parsed.maximumParticleSpeedMPerS;
  options.particleNumerics.timeStepMarginFactor =
      parsed.timeStepMarginFactor;
  options.particleNumerics.sourceNormalizationRadiusM =
      parsed.sourceNormalizationRadiusM;
  status=RM::RunConfiguration3D::Create(options, &resolved);
  if (!status.ok())
    return Fail("could not create parsed immutable configuration: "+
        status.message);
  RM::Runtime runtime;
  if (!runtime.Configure(provisional).ok() ||
      !runtime.ReplaceConfigurationBeforeMesh(resolved).ok() ||
      runtime.configuration()->options().source.samplesPerStep != 37 ||
      runtime.configuration()->physics_fingerprint() ==
          provisional->physics_fingerprint() ||
      runtime.configuration()->storage_layout() !=
          provisional->storage_layout()) {
    return Fail("pre-mesh parser transaction did not replace one immutable authority");
  }

  if (!write(directory / "bad.in",
             "#section begin: sep3d\n"
             "unknown_setting = 4\n"
             "#section end\n"))
    return Fail("could not write malformed fixture");
  status = RM::ParseSep3dApplicationInput(
      (directory / "bad.in").string(), &parsed);
  if (status.ok() || status.message.find("bad.in:2") == std::string::npos ||
      status.message.find("unknown_setting = 4") == std::string::npos ||
      status.message.find("unrecognized") == std::string::npos) {
    return Fail("malformed setting did not report source line, text, and cause");
  }

  if (!write(directory / "cycle-a.in", "#include cycle-b.in\n") ||
      !write(directory / "cycle-b.in", "#include cycle-a.in\n"))
    return Fail("could not write recursive-include fixture");
  status = RM::ParseSep3dApplicationInput(
      (directory / "cycle-a.in").string(), &parsed);
  if (status.ok() || status.message.find("cycle") == std::string::npos ||
      status.message.find(":1:") == std::string::npos) {
    return Fail("recursive include did not fail with directive provenance");
  }

  if (!write(directory / "unterminated.in",
             "#section begin: sep3d\nshock_model = reduced-shock-surface\n"))
    return Fail("could not write unterminated-section fixture");
  status = RM::ParseSep3dApplicationInput(
      (directory / "unterminated.in").string(), &parsed);
  if (status.ok() || status.message.find("without '#section end'") ==
                         std::string::npos) {
    return Fail("unterminated section was not rejected");
  }

  if (!write(directory / "missing-radius.in",
             "#section begin: sep3d\n"
             "shock_model=reduced-shock-surface\n"
             "background_plasma_model=corona-swcme-ambient\n"
             "source_model=accepted-shock-incident-flux\n"
             "maximum_time_steps=3\n"
             "particles_per_iteration=10\n"
             "maximum_particle_speed_m_s=1000\n"
             "time_step_margin_factor=0.2\n"
             "#subsection begin: reduced-shock-surface\n"
             "schema=x\n#subsection end\n#section end\n"))
    return Fail("could not write missing-radius fixture");
  status = RM::ParseSep3dApplicationInput(
      (directory / "missing-radius.in").string(), &parsed);
  if (status.ok() ||
      status.message.find("source_normalization_radius_m") ==
          std::string::npos)
    return Fail("missing source-normalization radius acquired a default");

  fs::remove_all(directory, cleanupError);
  return Pass("shared sep3d section resolves includes, comments, continuations, CLI defaults, immutable commit, and provenance-rich failures");
}

Result RunCFG3D14() {
  double step=0;
  if(!RM::CalculateMeshGlobalTimeStep(1200.0,300.0,0.25,&step).ok()||
      step!=1.0)
    return Fail("global dt does not equal margin*h_min/v_max");
  if(RM::CalculateMeshGlobalTimeStep(1200.0,300.0,1.01,&step).ok())
    return Fail("time-step calculation accepted a margin above one");
  std::vector<RM::ObserverOptions> observers(2);
  observers[0].id="unaligned";observers[0].cadenceS=60;
  observers[1].id="already-aligned";observers[1].cadenceS=70;
  if(!RM::AlignObserverCadencesToGlobalStep(7,&observers).ok()||
      observers[0].cadenceS!=63||observers[1].cadenceS!=70)
    return Fail("observer cadence did not align upward to exact global ticks");

  SEP::CoronaSwcme::ShockFront::IncidentParticleFlux flux;
  flux.protonRatePerS=4.0e20;
  flux.electronRatePerS=5.0e20;
  flux.alphaRatePerS=1.0e19;
  std::vector<RM::CompiledSpeciesRecord> species={
      {0,"ELECTRON",SEP3D::Core::Const::m_e,-SEP3D::Core::Const::e},
      {1,"H_PLUS",SEP3D::Core::Const::m_p,SEP3D::Core::Const::e},
      {2,"HE_PLUS_PLUS",SEP::CoronalCME::Constants::kAlphaMassKg,
          2*SEP3D::Core::Const::e}};
  std::vector<RM::SpeciesParticleNormalization> normalized;
  if(!RM::CalculateSpeciesParticleNormalizations(
      species,flux,2.0,100,&normalized).ok()||normalized.size()!=3||
      normalized[0].macroparticleWeight!=1.0e19||
      normalized[1].macroparticleWeight!=8.0e18||
      normalized[2].macroparticleWeight!=2.0e17)
    return Fail("electron/proton/alpha rate-to-weight normalization is wrong");
  species.push_back({3,"O_PLUS",16*SEP3D::Core::Const::m_p,
      SEP3D::Core::Const::e});
  if(RM::CalculateSpeciesParticleNormalizations(
      species,flux,2.0,100,&normalized).ok())
    return Fail("unknown ambient species silently borrowed another rate");
  return Pass("global dt and per-species W=Ndot*dt/N use independent inputs and unsupported composition fails closed");
}

Result RunCFG3D15() {
  namespace fs=std::filesystem;
  fs::path root;
  for(const fs::path& candidate:{fs::path("."),fs::path(".."),
      fs::path("../../")})
    if(fs::exists(candidate/"srcSEP3D/examples/application-input/amps.in")) {
      root=fs::canonical(candidate);break;
    }
  if(root.empty())return Fail("cannot locate the maintained shared-input example");

  RM::Sep3dApplicationInput parsed;
  SEP3D::Core::Status status=RM::ParseSep3dApplicationInput(
      (root/"srcSEP3D/examples/application-input/amps.in").string(),&parsed);
  if(!status.ok())return Fail("maintained shared input does not parse: "+
      status.message);

  // Mirror the production post-parser transaction, then ask the registered
  // factory to resolve the real magnetic asset and construct the shared
  // provider. This exercises the actual parser/configuration/adapter path;
  // it is not a second reference-only event parser.
  RM::RunConfiguration3DOptions options;
  options.meshMemoryBudgetBytes=1000000000000000ULL;
  options.inputSchemaVersion=4;
  options.maximumTimeSteps=parsed.maximumTimeSteps;
  options.background=RM::BackgroundAuthority::RuntimeModel;
  options.backgroundModelId="sep-corona-swcme-shock-front-v1";
  options.backgroundModelInlineConfiguration=parsed.reducedShockConfiguration;
  options.backgroundModelAssetDirectory=parsed.reducedShockAssetDirectory;
  options.coordinateFrame="HCI";
  options.parker.coordinateFrame=options.coordinateFrame;
  options.shock=RM::ShockAuthority::None;
  options.intent=RM::RunIntent::ShockInjection;
  options.source.enabled=true;
  options.source.samplesPerStep=parsed.particlesPerIteration;
  options.source.minimumEnergyJ=parsed.minimumInjectionEnergyJ;
  options.source.maximumEnergyJ=parsed.maximumInjectionEnergyJ;
  options.source.spectrumModel=parsed.momentumPowerLawModel=="compression-ratio"
      ? RM::SourceSpectrumModel::LocalCompressionDsa
      : RM::SourceSpectrumModel::FixedPhaseSpacePowerLaw;
  options.source.fixedPhaseSpacePowerIndex=parsed.fixedPhaseSpacePowerIndex;
  options.source.weightingModel=
      RM::SourceWeightingModel::ConstantStatisticalWeight;
  options.source.maximumMacroparticlesPerSpeciesPerStep=
      parsed.maximumInjectionEventsPerSpeciesPerStep;
  options.particleNumerics.deriveFromMeshAndShock=true;
  options.particleNumerics.sourceRateModel=
      RM::SourceRateNormalizationModel::AcceptedShockIncidentFlux;
  options.particleNumerics.maximumParticleSpeedMPerS=
      parsed.maximumParticleSpeedMPerS;
  options.particleNumerics.timeStepMarginFactor=
      parsed.timeStepMarginFactor;
  options.particleNumerics.sourceNormalizationRadiusM=
      parsed.sourceNormalizationRadiusM;
  std::shared_ptr<const RM::RunConfiguration3D> configuration;
  status=RM::RunConfiguration3D::Create(options,&configuration);
  if(!status.ok())return Fail("maintained parsed configuration is invalid: "+
      status.message);
  std::shared_ptr<BG::BackgroundProvider> background;
  status=RM::CreateBackgroundProvider(*configuration,&background);
  const auto adapter=std::dynamic_pointer_cast<
      SEP3D::Adapters::ShockFrontBackgroundAdapter>(background);
  if(!status.ok()||!adapter||!adapter->SharedProvider())
    return Fail("maintained parsed reduced provider does not initialize: "+
        status.message);
  const auto flux=SEP::CoronaSwcme::ShockFront::
      EvaluateIncidentParticleFluxAtApexRadius(*adapter->SharedProvider(),
          parsed.sourceNormalizationRadiusM);
  if(!flux.ok()||flux.value.acceptedAreaM2<=0||
      flux.value.electronRatePerS<=0||
      adapter->SharedProvider()->Current()!=nullptr)
    return Fail("maintained model cannot derive its non-mutating accepted-shock source rate");
  return Pass("maintained shared input resolves its magnetic asset, initializes one reduced provider, and derives a positive accepted-shock rate without advancing its epoch");
}

Result RunCFG3D16() {
  namespace fs=std::filesystem;
  fs::path root;
  for(const fs::path& candidate:{fs::path("."),fs::path(".."),
      fs::path("../../")})
    if(fs::exists(candidate/"srcSEP3D/examples/application-input/amps.in")) {
      root=fs::canonical(candidate);break;
    }
  if(root.empty())return Fail("cannot locate maintained source fixture");
  RM::Sep3dApplicationInput parsed;
  auto status=RM::ParseSep3dApplicationInput(
      (root/"srcSEP3D/examples/application-input/amps.in").string(),&parsed);
  if(!status.ok())return Fail("maintained source fixture does not parse: "+
      status.message);

  RM::RunConfiguration3DOptions options;
  options.meshMemoryBudgetBytes=1000000000000000ULL;
  options.inputSchemaVersion=4;
  options.maximumTimeSteps=parsed.maximumTimeSteps;
  options.background=RM::BackgroundAuthority::RuntimeModel;
  options.backgroundModelId="sep-corona-swcme-shock-front-v1";
  options.backgroundModelInlineConfiguration=parsed.reducedShockConfiguration;
  options.backgroundModelAssetDirectory=parsed.reducedShockAssetDirectory;
  options.coordinateFrame="HCI";
  options.parker.coordinateFrame="HCI";
  options.shock=RM::ShockAuthority::None;
  options.intent=RM::RunIntent::ShockInjection;
  options.source.enabled=true;
  options.source.minimumEnergyJ=parsed.minimumInjectionEnergyJ;
  options.source.maximumEnergyJ=parsed.maximumInjectionEnergyJ;
  options.source.spectrumModel=RM::SourceSpectrumModel::LocalCompressionDsa;
  options.source.fixedPhaseSpacePowerIndex=0;
  options.source.weightingModel=
      RM::SourceWeightingModel::ConstantStatisticalWeight;
  options.source.maximumMacroparticlesPerSpeciesPerStep=1000000;
  options.particleNumerics.deriveFromMeshAndShock=true;
  options.particleNumerics.sourceRateModel=
      RM::SourceRateNormalizationModel::AcceptedShockIncidentFlux;
  options.particleNumerics.maximumParticleSpeedMPerS=
      parsed.maximumParticleSpeedMPerS;
  options.particleNumerics.timeStepMarginFactor=parsed.timeStepMarginFactor;
  options.particleNumerics.sourceNormalizationRadiusM=
      parsed.sourceNormalizationRadiusM;
  std::shared_ptr<const RM::RunConfiguration3D> configuration;
  status=RM::RunConfiguration3D::Create(options,&configuration);
  if(!status.ok())return Fail("source configuration is invalid: "+status.message);
  std::shared_ptr<BG::BackgroundProvider> background;
  status=RM::CreateBackgroundProvider(*configuration,&background);
  const auto adapter=std::dynamic_pointer_cast<
      SEP3D::Adapters::ShockFrontBackgroundAdapter>(background);
  if(!status.ok()||!adapter||!adapter->SharedProvider())
    return Fail("source provider initialization failed: "+status.message);
  const auto provider=adapter->SharedProvider();
  const auto flux=SEP::CoronaSwcme::ShockFront::
      EvaluateIncidentParticleFluxAtApexRadius(
          *provider,parsed.sourceNormalizationRadiusM);
  if(!flux.ok())return Fail("reference incident flux failed: "+
      flux.status.message);
  const auto epoch=provider->EvaluateEpoch(flux.value.epochS,41);
  if(!epoch.ok())return Fail("source epoch evaluation failed: "+
      epoch.status.message);
  const RM::CompiledSpeciesRecord electron={0,"ELECTRON",
      SEP3D::Core::Const::m_e,-SEP3D::Core::Const::e};
  RM::SurfaceParticleRateDistribution distribution;
  status=RM::BuildSurfaceParticleRateDistribution(
      *provider,*epoch.value,electron,&distribution);
  const double rateScale=std::max(1.0,std::abs(flux.value.electronRatePerS));
  if(!status.ok()||distribution.faces.empty()||
      std::abs(distribution.physicalRatePerS-flux.value.electronRatePerS)>
          2e-14*rateScale||
      std::abs(distribution.acceptedAreaM2-flux.value.acceptedAreaM2)>
          2e-14*flux.value.acceptedAreaM2)
    return Fail("per-face source sum does not reproduce independent incident-flux diagnostic");

  const double interval=2.0;
  const double targetMean=64.0;
  const double weight=distribution.physicalRatePerS*interval/targetMean;
  RM::SurfaceInjectionBatch first,repeated;
  status=RM::GenerateConstantWeightSurfaceInjectionBatch(
      *provider,*epoch.value,electron,options.source,weight,interval,987654,7,
      &first);
  if(status.ok())status=RM::GenerateConstantWeightSurfaceInjectionBatch(
      *provider,*epoch.value,electron,options.source,weight,interval,987654,7,
      &repeated);
  if(!status.ok()||first.events.empty()||
      first.events.size()!=repeated.events.size())
    return Fail("deterministic constant-weight Poisson batch failed: "+
        status.message);
  const long double c=SEP3D::Core::Const::c;
  const long double mass=electron.massKg;
  const auto momentum=[&](double energy) {
    const long double k=energy;
    return static_cast<double>(std::sqrt(k*(k+2*mass*c*c))/c);
  };
  const double pMin=momentum(options.source.minimumEnergyJ);
  const double pMax=momentum(options.source.maximumEnergyJ);
  for(std::size_t index=0;index<first.events.size();++index) {
    const auto& event=first.events[index];
    const auto& again=repeated.events[index];
    if(event.stableId!=again.stableId||event.triangleStableId!=
        again.triangleStableId||event.eventTimeS!=again.eventTimeS||
        event.positionM.x!=again.positionM.x||
        event.positionM.y!=again.positionM.y||
        event.positionM.z!=again.positionM.z||
        !(event.eventTimeS>=0&&event.eventTimeS<interval)||
        !(event.remainingStepFraction>0&&
          event.remainingStepFraction<=1)||
        event.momentumKgMPerS<pMin||event.momentumKgMPerS>pMax)
      return Fail("source event is non-deterministic or outside declared time/momentum bounds");
    const auto& triangle=epoch.value->triangles[event.triangleIndex];
    const SEP3D::Core::Vec3 a(
        epoch.value->vertices[triangle.vertex[0]].positionM.x,
        epoch.value->vertices[triangle.vertex[0]].positionM.y,
        epoch.value->vertices[triangle.vertex[0]].positionM.z);
    const SEP3D::Core::Vec3 b(
        epoch.value->vertices[triangle.vertex[1]].positionM.x,
        epoch.value->vertices[triangle.vertex[1]].positionM.y,
        epoch.value->vertices[triangle.vertex[1]].positionM.z);
    const SEP3D::Core::Vec3 d(
        epoch.value->vertices[triangle.vertex[2]].positionM.x,
        epoch.value->vertices[triangle.vertex[2]].positionM.y,
        epoch.value->vertices[triangle.vertex[2]].positionM.z);
    const double total=(b-a).Cross(d-a).Norm();
    const double pieces=(b-event.positionM).Cross(d-event.positionM).Norm()+
        (d-event.positionM).Cross(a-event.positionM).Norm()+
        (a-event.positionM).Cross(b-event.positionM).Norm();
    if(total<=0||std::abs(pieces-total)>2e-12*total)
      return Fail("sampled source point is not inside its selected triangle");
  }

  std::uint64_t totalEvents=0;
  constexpr std::uint64_t trials=512;
  for(std::uint64_t trial=1;trial<=trials;++trial) {
    RM::SurfaceInjectionBatch sample;
    status=RM::GenerateConstantWeightSurfaceInjectionBatch(
        *provider,*epoch.value,electron,options.source,weight,interval,987654,
        100+trial,&sample);
    if(!status.ok())return Fail("Poisson convergence sample failed: "+
        status.message);
    totalEvents+=sample.events.size();
  }
  const double measured=static_cast<double>(totalEvents)/trials;
  const double standardError=std::sqrt(targetMean/trials);
  if(std::abs(measured-targetMean)>8*standardError)
    return Fail("Poisson sample mean is outside the preregistered eight-sigma bound");

  // The production-like launch is intentionally sub-fast at t=0.  Parse the
  // separate already-formed-front fixture and exercise the same factory and
  // accepted-face sum at its initial epoch; this guards the exact positive
  // condition used by the native one-/four-rank allocation smoke.
  RM::Sep3dApplicationInput smoke;
  status=RM::ParseSep3dApplicationInput((root/
      "srcSEP3D/examples/application-input/amps-injection-smoke.in").string(),
      &smoke);
  if(!status.ok())return Fail("positive native source fixture does not parse: "+
      status.message);
  RM::RunConfiguration3DOptions smokeOptions=options;
  smokeOptions.maximumTimeSteps=smoke.maximumTimeSteps;
  smokeOptions.backgroundModelInlineConfiguration=
      smoke.reducedShockConfiguration;
  smokeOptions.backgroundModelAssetDirectory=smoke.reducedShockAssetDirectory;
  smokeOptions.source.weightingModel=
      RM::SourceWeightingModel::ConstantStatisticalWeight;
  std::shared_ptr<const RM::RunConfiguration3D> smokeConfiguration;
  std::shared_ptr<BG::BackgroundProvider> smokeBackground;
  status=RM::RunConfiguration3D::Create(smokeOptions,&smokeConfiguration);
  if(status.ok())status=RM::CreateBackgroundProvider(
      *smokeConfiguration,&smokeBackground);
  const auto smokeAdapter=std::dynamic_pointer_cast<
      SEP3D::Adapters::ShockFrontBackgroundAdapter>(smokeBackground);
  if(status.ok())status=smokeBackground->Prepare(0.0);
  RM::SurfaceParticleRateDistribution smokeDistribution;
  if(status.ok()&&smokeAdapter&&smokeAdapter->FrontEpoch())
    status=RM::BuildSurfaceParticleRateDistribution(
        *smokeAdapter->SharedProvider(),*smokeAdapter->FrontEpoch(),electron,
        &smokeDistribution);
  if(!status.ok()||!smokeAdapter||smokeDistribution.faces.empty()||
      !(smokeDistribution.physicalRatePerS>0))
    return Fail("positive native source fixture has no accepted tick-zero rate: "+
        status.message);

  RM::SurfaceInjectionBatch reserved;
  options.source.weightingModel=
      RM::SourceWeightingModel::LogUniformMomentumImportance;
  status=RM::GenerateLogUniformMomentumImportanceBatch(
      *provider,*epoch.value,electron,options.source,weight,interval,987654,7,
      &reserved);
  if(status.code!=SEP3D::Core::StatusCode::ReservedFeature)
    return Fail("unimplemented momentum-importance mode did not fail explicitly");
  return Pass("accepted triangular rates reproduce the reference flux; keyed Poisson time/face/barycentric/momentum samples are deterministic and statistically convergent; reserved weighting fails closed");
}

}  // namespace

std::vector<SEP3D::Testing::Descriptor> RegisterConfigurationTests() {
  using D = SEP3D::Testing::Descriptor;
  using I = SEP3D::Testing::InitializationLevel;
  using R = SEP3D::Testing::RuntimeClass;
  auto make = [](const char* id, const char* name, const char* description,
                 SEP3D::Testing::TestCallback callback) {
    D d;
    d.id = id;
    d.name = name;
    d.group = "CFG3D";
    d.description = description;
    d.initialization = I::None;
    d.supportedBuildModes = "standalone-no-AMPS";
    d.runtime = R::Routine;
    d.seedPolicy = "deterministic";
    d.stateIsolation = "fresh immutable configuration per test";
    d.callback = std::move(callback);
    return d;
  };
  return {
      make("CFG3D01", "Input and CLI", "C01 schema, CLI, parity, errors, and dry-run.", RunCFG3D01),
      make("CFG3D02", "Typed contracts", "C02 complete groups and fingerprint classes.", RunCFG3D02),
      make("CFG3D03", "Domain contract", "C03 presets, containment, and boundary direction.", RunCFG3D03),
      make("CFG3D04", "Parker geometry", "C04 shared polarity-independent geometry.", RunCFG3D04),
      make("CFG3D05", "Mesh preflight", "C05 composite refinement and memory planning.", RunCFG3D05),
      make("CFG3D06", "Initialization schema", "Finite Parker-line fields and fail-closed consistency checks.", RunCFG3D06),
      make("CFG3D07", "Compiled-species binding", "Complete generated AMPS table and fingerprint contract.", RunCFG3D07),
      make("CFG3D08", "Complete initialization", "Schema-v3 canonical SWCME and exact per-step source contract.", RunCFG3D08),
      make("CFG3D09", "Background/turbulence selection", "Named prescribed slopes and reserved Python source.", RunCFG3D09),
      make("CFG3D10", "CME/Parker start linkage", "Canonical launch-apex linkage and fail-closed geometry checks.", RunCFG3D10),
      make("CFG3D11", "Transport/control schema", "Active corridor, population limits, mover coefficients, and dry-run output.", RunCFG3D11),
      make("CFG3D12", "Corner/sphere geometry", "Endpoint extent, input controls, identity and corner validation.", RunCFG3D12),
      make("CFG3D13", "Shared application input", "Global section/include grammar, early immutable commit, and diagnostics.", RunCFG3D13),
      make("CFG3D14", "Derived particle numerics", "Mesh CFL step and incident-flux per-species weight equations.", RunCFG3D14),
      make("CFG3D15", "Parsed reduced model", "Maintained shared input initializes the real reduced provider and source normalization.", RunCFG3D15),
      make("CFG3D16", "Reduced-front particle source", "Accepted-face rate sum, Poisson timing, triangular position, local spectrum, determinism, and reserved weighting.", RunCFG3D16),
  };
}

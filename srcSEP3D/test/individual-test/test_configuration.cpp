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

#include <cmath>
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
  return Pass("solar, one-AU, Mars, explicit-domain, containment, and directional boundary contracts passed");
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
    return Fail("complete schema-version-3 initialization input did not resolve: " +
                loaded.message);
  if (
      options.inputSchemaVersion != 3 ||
      options.injectionCadenceSteps != 1 ||
      options.source.samplesPerStep != 1000 ||
      options.swcmeAssignments.empty() ||
      options.swcmeConfigurationFingerprint.empty() ||
      options.swcmeResolvedManifest.empty()) {
    return Fail("schema-version-3 initialization metadata is incomplete");
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

  return Pass("schema version 3 resolves complete SWCME physics, honors explicit Parker references/source radii, and rejects missing or inconsistent physics");
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
  };
}

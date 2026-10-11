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
#include "particle_normalization.h"
#include "background_factory.h"
#include "shock_front_background_adapter.h"
#include "runtime.h"
#include "parallel_diffusion/parallel_diffusion.h"
#include "sep_input_layout.h"
#include "swcme3d_input.hpp"
#include "sep_coronal_cme/constants.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <set>
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

// Text of a shipped example deck (tests run from the srcSEP3D directory).
std::string ReadText(const std::string& path) {
  std::ifstream input(path);
  std::ostringstream buffer;
  buffer << input.rdbuf();
  return input ? buffer.str() : std::string();
}

// The complete schema-6 fixture is the maintained analytic-Parker example:
// every srcSEP3D key is explicit there, so the parser's completeness rules
// are exercised against the deck users actually start from.
std::string CompleteInput() {
  return ReadText("examples/sep3d_analytic_parker.in");
}

// Replace the first occurrence of 'from' (which must exist) by 'to'; an
// absent anchor returns false so a drifting fixture fails loudly.
bool ReplaceFirst(std::string* text, const std::string& from,
                  const std::string& to) {
  const std::size_t at = text == nullptr ? std::string::npos : text->find(from);
  if (at == std::string::npos) return false;
  text->replace(at, from.size(), to);
  return true;
}

// Presentation-only rewrite used by CFG3D01: drop every '!' comment and move
// the output block to the top of the file.  The configuration fingerprint is
// computed from normalized values, so it must not change.
std::string WithoutCommentsOutputFirst(const std::string& text) {
  std::istringstream input(text);
  std::string line, body;
  while (std::getline(input, line)) {
    const std::size_t bang = line.find('!');
    if (bang != std::string::npos) line.erase(bang);
    if (line.find_first_not_of(" \t") == std::string::npos) continue;
    body += line + "\n";
  }
  const std::size_t begin = body.find("#begin output");
  const std::size_t end = body.find("#end output", begin);
  if (begin == std::string::npos || end == std::string::npos) return body;
  const std::size_t stop = body.find('\n', end) + 1;
  const std::string output = body.substr(begin, stop - begin);
  body.erase(begin, stop - begin);
  return output + body;
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
  const std::string complete = CompleteInput();
  const SEP3D::Core::Status completeStatus =
      RM::ParseConfigurationText(complete, &parsed);
  if (complete.empty() || !completeStatus.ok())
    return Fail("the complete schema-6 fixture did not parse: " +
                completeStatus.message);
  std::shared_ptr<const RM::RunConfiguration3D> fromFile;
  if (!RM::RunConfiguration3D::Create(parsed, &fromFile).ok())
    return Fail("parsed options did not pass immutable construction");
  if (parsed.inputSchemaVersion != 6 || !parsed.particleMoverExplicitlySelected ||
      parsed.particleMover != RM::ParticleMoverFamily::Parker)
    return Fail("schema 6 did not install contract level 6 and the mover "
                "family derived from application.sep.transport");

  // The identity is computed from normalized values: comments, layout and
  // block order must not change it.
  RM::RunConfiguration3DOptions reformatted;
  std::shared_ptr<const RM::RunConfiguration3D> fromReformatted;
  if (!RM::ParseConfigurationText(WithoutCommentsOutputFirst(complete),
                                  &reformatted).ok() ||
      !RM::RunConfiguration3D::Create(reformatted, &fromReformatted).ok() ||
      fromFile->physics_fingerprint() != fromReformatted->physics_fingerprint())
    return Fail("comments or block order changed the physics fingerprint");

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
  bad.insert(bad.find("#begin output\n") + std::string("#begin output\n").size(),
             "  unknown_key = x\n");
  if (RM::ParseConfigurationText(bad, &parsed).ok())
    return Fail("unknown input key was not rejected before initialization");
  std::string duplicateBlock = CompleteInput() + "#begin output\n#end output\n";
  if (RM::ParseConfigurationText(duplicateBlock, &parsed).ok())
    return Fail("a duplicate block was not rejected before initialization");
  std::string strayBlock = CompleteInput() + "#begin photolysis\n#end photolysis\n";
  const SEP3D::Core::Status strayStatus =
      RM::ParseConfigurationText(strayBlock, &parsed);
  if (strayStatus.ok() ||
      strayStatus.message.find("'photolysis' is not part of the SEP input") ==
          std::string::npos)
    return Fail("a block of a model absent from srcSEP3D was not rejected: " +
                strayStatus.message);
  // Molecular identity belongs exclusively to the AMPS build-time
  // SpeciesList.  A post-compile attempt to redefine mass must be an unknown
  // key rather than a silently ignored compatibility spelling.
  std::string retiredSpeciesField = CompleteInput();
  if (!ReplaceFirst(&retiredSpeciesField, "macroparticle_weight = 1e23\n",
                    "macroparticle_weight = 1e23\n  mass_kg = 1.67262192369e-27\n") ||
      RM::ParseConfigurationText(retiredSpeciesField, &parsed).ok())
    return Fail("runtime input was allowed to redefine compiled species mass");
  bad = CompleteInput();
  const std::size_t firstObserver = bad.find("#begin observer:");
  const std::size_t outputBlock = bad.find("#begin output");
  if (firstObserver == std::string::npos || outputBlock < firstObserver)
    return Fail("fixture lost its observer blocks");
  bad.erase(firstObserver, outputBlock - firstObserver);
  if (RM::ParseConfigurationText(bad, &parsed).ok())
    return Fail("a missing required observer group was not rejected");

  // Shared SEP schema 6 (src/INPUT_FILE.md section 5): each rejection names
  // its decision.  Expected phrases are written here, not read from code.
  struct Rejection { const char* label; std::string text; const char* expect; };
  std::string retiredSchema = complete, srcSepKey = complete,
      srcSepTurbulence = complete, parkerAmbient = complete,
      swcmeSource = complete, retiredMover = complete, oldBlock = complete,
      fieldLineObserver = complete;
  const bool edited =
      ReplaceFirst(&retiredSchema, "schema_version = 6", "schema_version = 4") &&
      ReplaceFirst(&srcSepKey, "kind = fixed-cartesian",
                   "kind = fixed-cartesian\n  heliocentric_radius_m = 1e11") &&
      ReplaceFirst(&srcSepTurbulence, "authority = prescribed",
                   "authority = prescribed\n    cascade = true") &&
      ReplaceFirst(&parkerAmbient, "validity_cadence_s = 3600",
                   "validity_cadence_s = 3600\n    solar_wind_speed_m_per_s = 400000") &&
      ReplaceFirst(&swcmeSource, "surface.phi_points = 24",
                   "surface.phi_points = 24\n  source.energy_min = 0.01 MeV") &&
      ReplaceFirst(&retiredMover, "transport = parker",
                   "transport = parker\n  particle_mover = parker") &&
      ReplaceFirst(&oldBlock, "#begin application: sep", "#begin application: sep3d") &&
      ReplaceFirst(&fieldLineObserver, "kind = fixed-cartesian",
                   "kind = field-line-radius");
  if (!edited) return Fail("schema-6 fixture lost an anchor for its negative controls");
  const Rejection rejections[] = {
      {"retired schema (D1)", retiredSchema, "retired"},
      {"retired application name (D1)", oldBlock, "is not part of the SEP input"},
      {"srcSEP-only observer key (D3)", srcSepKey, "not supported by srcSEP3D"},
      {"srcSEP-only observer kind (D3)", fieldLineObserver, "not supported by srcSEP3D"},
      {"srcSEP turbulence key (D2/D3)", srcSepTurbulence, "not supported by srcSEP3D"},
      {"Parker ambient beside SWCME (D4)", parkerAmbient, "owned by physics: swcme"},
      {"SWCME source key (D4)", swcmeSource, "written in source: shock"},
      {"retired particle_mover key", retiredMover, "unknown configuration key"}};
  for (const Rejection& rejection : rejections) {
    const SEP3D::Core::Status status =
        RM::ParseConfigurationText(rejection.text, &parsed);
    if (status.ok() || status.message.find(rejection.expect) == std::string::npos)
      return Fail(std::string("negative control '") + rejection.label +
                  "' was not rejected with '" + rejection.expect + "': " +
                  status.message);
  }
  const SEP3D::Core::Status flatStatus =
      RM::ParseConfigurationText("scattering = prescribed\n", &parsed);
  if (flatStatus.ok() ||
      flatStatus.message.find("outside any block") == std::string::npos) {
    return Fail("a flat deck did not receive the block-aware diagnostic");
  }
  // The retired INI section syntax is not part of the generic grammar.
  const SEP3D::Core::Status iniStatus =
      RM::ParseConfigurationText("[run]\ntime_step_s = 1\n", &parsed);
  if (iniStatus.ok() || iniStatus.message.find(":1:") == std::string::npos)
    return Fail("a retired INI deck was not rejected with its line number");

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
  return Pass("schema-6 input, layout-independent identity, CLI normalization, D1-D4 rejections, early errors, and allocation-free dry-run passed");
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
  // The finite Parker-line definition is a required part of schema 6.
  RM::RunConfiguration3DOptions options;
  const std::string complete = CompleteInput();
  const SEP3D::Core::Status status =
      RM::ParseConfigurationText(complete, &options);
  if (!status.ok() || options.parkerSpiralPointCount != 4001 ||
      options.parkerSpiralLengthM != 2.0e11) {
    return Fail("complete finite Parker-line definition did not parse exactly: " +
                status.message);
  }

  std::string missing = complete;
  if (!ReplaceFirst(&missing, "point_count = 4001\n", "") ||
      RM::ParseConfigurationText(missing, &options).ok())
    return Fail("schema 6 accepted a missing finite-line point count");

  std::string inconsistent = complete;
  if (!ReplaceFirst(&inconsistent, "initial_x_m = 1.3914e10",
                    "initial_x_m = 1.4e10") ||
      RM::ParseConfigurationText(inconsistent, &options).ok())
    return Fail("schema 6 accepted a line source outside the inner boundary");
  return Pass("schema 6 requires and validates the complete finite Parker-line definition");
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
    return Fail("complete schema-6 initialization input did not resolve: " +
                loaded.message);
  if (
      options.inputSchemaVersion != 6 ||
      options.injectionCadenceSteps != 1 ||
      options.source.samplesPerStep != 1000 ||
      options.swcmeAssignments.empty() ||
      options.swcmeConfigurationFingerprint.empty() ||
      options.swcmeResolvedManifest.empty()) {
    return Fail("schema-6 initialization metadata is incomplete");
  }

  std::ifstream input("examples/sep3d_analytic_parker.in");
  std::ostringstream buffer;
  buffer << input.rdbuf();
  const std::string complete = buffer.str();
  if (!input || complete.empty()) return Fail("could not read the schema-6 fixture");

  std::string missing = complete;
  const std::string required = "surface.phi_points = 24\n";
  const std::size_t requiredAt = missing.find(required);
  if (requiredAt == std::string::npos)
    return Fail("schema-6 fixture lost its surface resolution field");
  missing.erase(requiredAt, required.size());
  if (RM::ParseConfigurationText(missing, &options).ok())
    return Fail("schema 6 accepted an incomplete canonical SWCME model");

  std::string wrongWeight = complete;
  const std::string weight = "macroparticle_weight = 1e23";
  const std::size_t weightAt = wrongWeight.find(weight);
  if (weightAt == std::string::npos)
    return Fail("schema-6 fixture lost its derived particle weight");
  wrongWeight.replace(weightAt, weight.size(),
                      "macroparticle_weight = 2e23");
  if (RM::ParseConfigurationText(wrongWeight, &options).ok())
    return Fail("schema 6 accepted a particle weight inconsistent with its physical source rate");

  std::string skippedStep = complete;
  const std::string cadence = "injection_cadence_steps = 1";
  const std::size_t cadenceAt = skippedStep.find(cadence);
  if (cadenceAt == std::string::npos)
    return Fail("schema-6 fixture lost its injection cadence");
  skippedStep.replace(cadenceAt, cadence.size(),
                      "injection_cadence_steps = 2");
  if (RM::ParseConfigurationText(skippedStep, &options).ok())
    return Fail("schema 6 accepted a source that skips simulation steps");

  std::string absoluteNormalization = complete;
  const std::string relative = "spectrum_normalization = relative-only";
  const std::size_t relativeAt = absoluteNormalization.find(relative);
  if (relativeAt == std::string::npos)
    return Fail("schema-6 fixture lost its source normalization");
  absoluteNormalization.replace(
      relativeAt, relative.size(),
      "spectrum_normalization = reference-differential-intensity");
  if (RM::ParseConfigurationText(absoluteNormalization, &options).ok())
    return Fail("schema 6 accepted two competing absolute source normalizations");

  std::string wrongAxis = complete;
  const std::string axisX = "geometry.solar_rotation_axis_x = 0";
  const std::string axisZ = "geometry.solar_rotation_axis_z = 1";
  const std::size_t axisXAt = wrongAxis.find(axisX);
  const std::size_t axisZAt = wrongAxis.find(axisZ);
  if (axisXAt == std::string::npos || axisZAt == std::string::npos)
    return Fail("schema-6 fixture lost its solar rotation axis");
  wrongAxis.replace(axisXAt, axisX.size(),
                    "geometry.solar_rotation_axis_x = 1");
  wrongAxis.replace(wrongAxis.find(axisZ), axisZ.size(),
                    "geometry.solar_rotation_axis_z = 0");
  if (RM::ParseConfigurationText(wrongAxis, &options).ok())
    return Fail("schema 6 accepted inconsistent Parker rotation axes");

  // Decision D4: the SWCME ambient is the only written solar wind.  SWCME
  // defines total |B| at one AU, while AnalyticParkerProvider uses the radial
  // component at its reference radius (1 AU when SWCME is present).  The
  // expected B_r values below were the reviewed hand-computed Parker inputs
  // of the schema-4 decks (|B| = 5 nT, V = 400 km/s, Omega = 2.865e-6 rad/s,
  // source radius 20 and 25 R_sun); they are independent of the parser.
  const auto relativeDifference = [](double a, double b) {
    return std::fabs(a - b) / std::max(std::fabs(a), std::fabs(b));
  };
  if (!RM::ParseConfigurationText(complete, &options).ok() ||
      options.parker.referenceRadiusM != SEP3D::Core::Const::AU ||
      relativeDifference(options.parker.radialFieldAtReferenceT,
                         3.585667175612118e-9) > 1.0e-12 ||
      options.parker.solarWindSpeedMPerS != 4.0e5 ||
      options.parker.numberDensityAtReferenceM3 != 5.0e6)
    return Fail("the analytic Parker record was not derived from the SWCME ambient");
  std::string ownedReference = complete;
  if (!ReplaceFirst(&ownedReference, "validity_cadence_s = 3600",
                    "validity_cadence_s = 3600\n    reference_radius_m = 7.479893535e10") ||
      RM::ParseConfigurationText(ownedReference, &options).ok())
    return Fail("a Parker reference radius was accepted beside physics: swcme");

  // The source surface is explicit physics, not a hidden 20-R_sun constant.
  // Move the domain, finite line, Parker source, and CME launch surface to
  // 25 R_sun; the derived one-AU radial component follows.
  std::string movedSource = complete;
  if (!ReplaceFirst(&movedSource, "inner_radius_m = 1.3914e10",
                    "inner_radius_m = 1.73925e10") ||
      !ReplaceFirst(&movedSource, "initial_x_m = 1.3914e10",
                    "initial_x_m = 1.73925e10") ||
      !ReplaceFirst(&movedSource, "parker.source_radius = 20 Rs",
                    "parker.source_radius = 25 Rs") ||
      !ReplaceFirst(&movedSource, "cme.launch_radius = 20 Rs",
                    "cme.launch_radius = 25 Rs"))
    return Fail("schema-6 fixture lost an explicit source-surface field");
  if (!RM::ParseConfigurationText(movedSource, &options).ok() ||
      relativeDifference(options.parker.radialFieldAtReferenceT,
                         3.6305743892405979e-9) > 1.0e-12)
    return Fail("schema 6 retained a hidden 20-R_sun source assumption");

  return Pass("schema 6 keeps the complete SWCME physics contract, derives the Parker record from the SWCME ambient (D4), honors explicit source radii, and rejects missing or inconsistent physics");
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
      options.inputSchemaVersion != 6 ||
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
      !replaceOnce(&corridor, "radius_at_reference_m = 0\n    radius_mode = constant-angular-width\n    buffer_blocks = 0",
                   "radius_at_reference_m = 7.479893535e9\n    radius_mode = constant-angular-width\n    buffer_blocks = 1") ||
      !replaceOnce(&corridor,
                   "position_x_m = 1.495978707e11\n  position_y_m = 0",
                   "position_x_m = 1.10962252663803e11\n  position_y_m = -1.00333949397740e11") ||
      !replaceOnce(&corridor,
                   "position_x_m = 4.487936121e10\n  position_y_m = 0",
                   "position_x_m = 4.46318123441847e10\n  position_y_m = -4.70726985535542e9") ||
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
                   "position_x_m = 1.10962252663803e11\n  position_y_m = -1.00333949397740e11",
                   "position_x_m = 1.495978707e11\n  position_y_m = 0") ||
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
  if (!replaceOnce(&controlled, "mode = off\n    minimum_particles_per_cell_per_species = 0\n    target_particles_per_cell_per_species = 0\n    maximum_particles_per_cell_per_species = 0",
                   "mode = split-merge\n    minimum_particles_per_cell_per_species = 16\n    target_particles_per_cell_per_species = 24\n    maximum_particles_per_cell_per_species = 32") ||
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
      options.transport != RM::TransportModel::FocusedScattering3D ||
      options.particleMover != RM::ParticleMoverFamily::FocusedTransport) {
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
  // The maintained corner/sphere example selects the x-y corner box; the
  // three-axis corner mode is the same deck with the original geometry name.
  const std::string xyDeck =
      ReadText("examples/sep3d_analytic_parker_corner_sphere.in");
  std::string input = xyDeck;
  if (xyDeck.empty() ||
      !ReplaceFirst(&input, "box_geometry = field-line-xy-corner-cube",
                    "box_geometry = field-line-corner-cube"))
    return Fail("missing corner/sphere example");
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
  const std::string xyInput = xyDeck;
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
  // Unified generic runtime input (src/INPUT_FILE.md) for reduced-front
  // injection: includes, '!' comments, '\' continuation, the derived particle
  // numerics keys, the three CLI spellings, and provenance-rich failures.
  // The fixture includes the maintained shared blocks of the application-input
  // example by absolute path and a parser-only physics: shock_front block (the
  // shared resolver is exercised by CFG3D15/16, not here).
  namespace fs = std::filesystem;
  const fs::path directory = "test_output/cfg3d13-application-input";
  std::error_code cleanupError;
  fs::remove_all(directory, cleanupError);
  cleanupError.clear();
  fs::create_directories(directory / "parts", cleanupError);
  if (cleanupError) return Fail("could not create application-input fixture directory");
  const fs::path shared = fs::absolute("examples/application-input/sep3d.in");
  if (!fs::exists(shared)) return Fail("cannot locate examples/application-input/sep3d.in");

  const auto write = [](const fs::path& path, const std::string& text) {
    std::ofstream stream(path);
    stream << text;
    stream.close();
    return static_cast<bool>(stream);
  };
  const std::string runBlock =
      "#begin run\n"
      "  time_step_mode = derived-from-mesh\n"
      "  maximum_particle_speed_m_s = 200000000\n"
      "  time_step_margin_factor = 0.25\n"
      "  maximum_time_steps = 9\n"
      "  campaign_seed = 1\n"
      "  #begin restart\n"
      "    input_path = none\n"
      "    output_path = restart/sep3d.chk\n"
      "  #end restart\n"
      "#end run\n";
  const std::string sourceBlock =
      "#begin source: shock\n"
      "  enabled = true\n"
      "  normalization = accepted-shock-incident-flux\n"
      "  normalization_radius_m = 13914000000\n"
      "  samples_per_step = \\\n"
      "    37 ! joined to the prior physical line\n"
      "  minimum_energy_j = 1.602176634e-15\n"
      "  maximum_energy_j = 1.602176634e-11\n"
      "  spectrum_model = local-compression-dsa\n"
      "  phase_space_power_index = 0\n"
      "  weight_model = constant-statistical-weight\n"
      "  maximum_events_per_species_per_step = 1000000\n"
      "#end source\n";
  if (!write(directory / "amps.in",
             "! unified srcSEP3D input fixture\n" + runBlock + sourceBlock +
             "#include \"" + shared.string() + "\"\n"
             "#include \"parts/model.in\"\n") ||
      !write(directory / "parts/model.in",
             "#begin physics: shock_front\n"
             "  schema = parser-fixture\n"
             "#end physics\n")) {
    return Fail("could not write application-input fixtures");
  }

  RM::RunConfiguration3DOptions parsed;
  SEP3D::Core::Status status =
      RM::LoadConfigurationFile((directory / "amps.in").string(), &parsed);
  const fs::path expectedAssets =
      fs::weakly_canonical(fs::absolute(directory / "parts"));
  if (!status.ok() || parsed.maximumTimeSteps != 9 ||
      parsed.source.samplesPerStep != 37 ||
      !parsed.particleNumerics.deriveFromMeshAndShock ||
      parsed.particleNumerics.sourceRateModel !=
          RM::SourceRateNormalizationModel::AcceptedShockIncidentFlux ||
      parsed.particleNumerics.maximumParticleSpeedMPerS != 2.0e8 ||
      parsed.particleNumerics.timeStepMarginFactor != 0.25 ||
      parsed.particleNumerics.sourceNormalizationRadiusM != 13914000000.0 ||
      parsed.source.weightingModel !=
          RM::SourceWeightingModel::ConstantStatisticalWeight ||
      parsed.source.spectrumModel != RM::SourceSpectrumModel::LocalCompressionDsa ||
      parsed.source.minimumEnergyJ != 1.602176634e-15 ||
      parsed.source.maximumEnergyJ != 1.602176634e-11 ||
      parsed.source.maximumMacroparticlesPerSpeciesPerStep != 1000000 ||
      parsed.source.injectionEfficiency != 1.0 ||
      parsed.intent != RM::RunIntent::ShockInjection ||
      parsed.shock != RM::ShockAuthority::None ||
      parsed.background != RM::BackgroundAuthority::RuntimeModel ||
      parsed.backgroundModelInlineConfiguration !=
          "#begin physics: shock_front\n  schema = parser-fixture\n#end physics\n" ||
      fs::path(parsed.backgroundModelAssetDirectory) != expectedAssets) {
    return Fail("include/comment/continuation input did not resolve: " +
                status.message);
  }

  // One input format, three spellings: -input FILE, its synonym --input FILE,
  // and no argument (./amps.in, also used by a coupled SWMF run).
  RM::StandaloneCommandLine cli;
  const char* defaultArgv[] = {"amps"};
  if (!RM::ParseStandaloneCommandLine(
          1, const_cast<char**>(defaultArgv), &cli).ok() ||
      cli.inputPath != "amps.in")
    return Fail("no-argument CLI did not select ./amps.in");
  const char* oneDash[] = {"amps", "-input", "global.in"};
  if (!RM::ParseStandaloneCommandLine(
          3, const_cast<char**>(oneDash), &cli).ok() ||
      cli.inputPath != "global.in")
    return Fail("-input did not select the named input file");
  const char* twoDash[] = {"amps", "--input", "global.in"};
  if (!RM::ParseStandaloneCommandLine(
          3, const_cast<char**>(twoDash), &cli).ok() ||
      cli.inputPath != "global.in")
    return Fail("--input is not a synonym of -input");

  // The parsed derived-mode options are a valid immutable configuration.  A
  // later commit of the mesh-derived step (FinalizeParticleNumerics) changes
  // the physics identity but never the storage layout registered with AMPS.
  std::shared_ptr<const RM::RunConfiguration3D> initial, derived;
  status = RM::RunConfiguration3D::Create(parsed, &initial);
  if (!status.ok())
    return Fail("parsed derived-mode configuration is invalid: " + status.message);
  RM::RunConfiguration3DOptions afterMesh = initial->options();
  afterMesh.requestedTimeStepS = 1.25;  // = 0.25 * 1e9 m / 2e8 m/s
  afterMesh.particleNumerics.resolvedMinimumCellSizeM = 1.0e9;
  RM::Runtime runtime;
  status = RM::RunConfiguration3D::Create(afterMesh, &derived);
  if (status.ok()) status = runtime.Configure(initial);
  if (status.ok()) status = runtime.ReplaceConfigurationBeforeMesh(derived);
  if (!status.ok() ||
      runtime.configuration()->physics_fingerprint() ==
          initial->physics_fingerprint() ||
      runtime.configuration()->storage_layout() != initial->storage_layout())
    return Fail("derived particle numerics did not replace one immutable "
                "authority: " + status.message);

  // Failures name the file, line, cause and offending text.
  const auto expectFailure = [&](const std::string& name, const std::string& text,
                                 const std::string& phrase) -> std::string {
    if (!write(directory / name, text)) return "could not write " + name;
    RM::RunConfiguration3DOptions unused;
    const SEP3D::Core::Status failed =
        RM::LoadConfigurationFile((directory / name).string(), &unused);
    if (failed.ok()) return name + " was accepted";
    if (failed.message.find(phrase) == std::string::npos)
      return name + " diagnostic lacks '" + phrase + "': " + failed.message;
    return std::string();
  };
  std::string problem = expectFailure("bad.in",
      runBlock + "#begin source: shock\n  unknown_setting = 4\n#end source\n",
      "bad.in:13:");
  if (problem.empty()) problem = expectFailure("bad-text.in",
      runBlock + "#begin source: shock\n  unknown_setting = 4\n#end source\n",
      "unknown_setting = 4");
  if (problem.empty()) {
    write(directory / "cycle-b.in", "#include cycle-a.in\n");
    problem = expectFailure("cycle-a.in", "#include cycle-b.in\n", "cycle");
  }
  if (problem.empty()) problem = expectFailure("unterminated.in",
      "#begin run\n  maximum_time_steps = 1\n", "is not closed");
  std::string missingRadius = sourceBlock;
  missingRadius.replace(missingRadius.find("  normalization_radius_m = 13914000000\n"),
                        std::string("  normalization_radius_m = 13914000000\n").size(), "");
  if (problem.empty()) problem = expectFailure("missing-radius.in",
      runBlock + missingRadius + "#include \"" + shared.string() +
      "\"\n#include \"parts/model.in\"\n", "source.shock.normalization_radius_m");
  std::string fixedStep = runBlock;
  fixedStep.replace(fixedStep.find("time_step_mode = derived-from-mesh"),
                    std::string("time_step_mode = derived-from-mesh").size(),
                    "time_step_mode = derived-from-mesh\n  time_step_s = 1");
  if (problem.empty()) problem = expectFailure("written-step.in",
      fixedStep + sourceBlock + "#include \"" + shared.string() +
      "\"\n#include \"parts/model.in\"\n", "'run.time_step_s' is derived");
  if (!problem.empty()) return Fail(problem);

  fs::remove_all(directory, cleanupError);
  return Pass("unified input resolves includes, comments, continuations, derived particle numerics, CLI spellings, immutable commit, and provenance-rich failures");
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

  // AMPS species tables use the A*u mass convention: H_PLUS = 1 u and
  // HE_2PLUS = 4 u (u = 1.66053906660e-27 kg, CODATA 2018), not the CODATA
  // proton/alpha masses.  Both must map onto the proton and alpha
  // populations; a deuteron (A=2, Z=1) and 3He++ (A=3, Z=2) must not.
  const double u=1.66053906660e-27;
  const std::vector<RM::CompiledSpeciesRecord> ampsConvention={
      {0,"H_PLUS",1*u,SEP3D::Core::Const::e},
      {1,"HE_2PLUS",4*u,2*SEP3D::Core::Const::e}};
  if(!RM::CalculateSpeciesParticleNormalizations(
      ampsConvention,flux,2.0,100,&normalized).ok()||normalized.size()!=2||
      normalized[0].macroparticleWeight!=8.0e18||
      normalized[1].macroparticleWeight!=2.0e17)
    return Fail("AMPS 1-u H_PLUS / 4-u alpha are not mapped to the proton/alpha populations");
  for(const RM::CompiledSpeciesRecord& nuclide:
      {RM::CompiledSpeciesRecord{0,"D_PLUS",2*u,SEP3D::Core::Const::e},
       RM::CompiledSpeciesRecord{0,"HE3_2PLUS",3*u,2*SEP3D::Core::Const::e},
       RM::CompiledSpeciesRecord{0,"H_2PLUS_BOGUS",1*u,2*SEP3D::Core::Const::e}})
    if(RM::CalculateSpeciesParticleNormalizations(
        {nuclide},flux,2.0,100,&normalized).ok())
      return Fail("nuclide "+nuclide.symbol+" borrowed a proton/alpha rate");
  return Pass("global dt and per-species W=Ndot*dt/N use independent inputs, the AMPS A*u mass convention maps to protons/alphas, and unsupported composition fails closed");
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

  // Parse the maintained unified example with the production reader; the
  // registered factory then resolves the real magnetic asset and constructs
  // the shared provider.  This exercises the actual parser/configuration/
  // adapter path; it is not a second reference-only event parser.
  RM::RunConfiguration3DOptions options;
  SEP3D::Core::Status status=RM::LoadConfigurationFile(
      (root/"srcSEP3D/examples/application-input/amps.in").string(),&options);
  if(!status.ok())return Fail("maintained unified input does not parse: "+
      status.message);
  options.meshMemoryBudgetBytes=1000000000000000ULL;
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
          options.particleNumerics.sourceNormalizationRadiusM);
  if(!flux.ok()||flux.value.acceptedAreaM2<=0||
      flux.value.electronRatePerS<=0||
      adapter->SharedProvider()->Current()!=nullptr)
    return Fail("maintained model cannot derive its non-mutating accepted-shock source rate");
  return Pass("maintained unified input resolves its magnetic asset, initializes one reduced provider, and derives a positive accepted-shock rate without advancing its epoch");
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
  RM::RunConfiguration3DOptions options;
  auto status=RM::LoadConfigurationFile(
      (root/"srcSEP3D/examples/application-input/amps.in").string(),&options);
  if(!status.ok())return Fail("maintained source fixture does not parse: "+
      status.message);
  options.meshMemoryBudgetBytes=1000000000000000ULL;
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
          *provider,options.particleNumerics.sourceNormalizationRadiusM);
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

  // The production entry point receives provider and epoch separately so that
  // normalization preflight need not mutate Provider::Current(). Prove that
  // this flexibility cannot mix one event's accepted faces with another
  // event's EOS/composition. Mutating only the immutable provenance label is
  // an independent negative control; the destination sentinel also checks the
  // documented transactional output contract.
  auto foreignEpoch=*epoch.value;
  foreignEpoch.eventIdentity="foreign-reduced-front-event";
  RM::SurfaceParticleRateDistribution provenanceSentinel;
  provenanceSentinel.generation=UINT64_C(0x123456789abcdef0);
  provenanceSentinel.physicalRatePerS=17.0;
  status=RM::BuildSurfaceParticleRateDistribution(
      *provider,foreignEpoch,electron,&provenanceSentinel);
  if(status.code!=SEP3D::Core::StatusCode::ConfigurationConflict||
      provenanceSentinel.generation!=UINT64_C(0x123456789abcdef0)||
      provenanceSentinel.physicalRatePerS!=17.0)
    return Fail("surface source accepted a foreign epoch or mutated output on provenance failure");

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
  RM::RunConfiguration3DOptions smokeOptions;
  status=RM::LoadConfigurationFile((root/
      "srcSEP3D/examples/application-input/amps-injection-smoke.in").string(),
      &smokeOptions);
  if(!status.ok())return Fail("positive native source fixture does not parse: "+
      status.message);
  smokeOptions.meshMemoryBudgetBytes=1000000000000000ULL;
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
  return Pass("accepted triangular rates reproduce the reference flux; foreign event epochs fail transactionally; keyed Poisson time/face/barycentric/momentum samples are deterministic and statistically convergent; reserved weighting fails closed");
}

Result RunCFG3D17() {
  std::ifstream input("examples/sep3d_analytic_parker.in");
  std::ostringstream buffer;
  buffer << input.rdbuf();
  if (!input || buffer.str().empty())
    return Fail("could not read the schema-6 fixture");

  auto replaceOnce = [](std::string* text, const std::string& from,
                        const std::string& to) {
    const std::size_t at = text == nullptr ? std::string::npos : text->find(from);
    if (at == std::string::npos) return false;
    text->replace(at, from.size(), to);
    return true;
  };
  std::string schema5 = buffer.str();
  if (!replaceOnce(&schema5,
          "spatial_diffusion_model = mean-free-path",
          "spatial_diffusion_model = parallel-diffusion-library"))
    return Fail("schema-6 fixture lost its spatial selector");
  const std::string shockMarker = "\n#begin physics: shock\n";
  const std::size_t shockAt = schema5.find(shockMarker);
  if (shockAt == std::string::npos)
    return Fail("schema-6 fixture lost its shock block");
  schema5.insert(shockAt, R"SEP3D(
#begin physics: parallel_diffusion
model = constant_kappa
kappa_parallel_m2_per_s = 1.25e18
#end physics
)SEP3D");

  RM::RunConfiguration3DOptions options;
  SEP3D::Core::Status status = RM::ParseConfigurationText(schema5, &options);
  if (!status.ok() || options.inputSchemaVersion != 6 ||
      options.particleMover != RM::ParticleMoverFamily::Parker ||
      !options.particleMoverExplicitlySelected ||
      options.spatialDiffusionModel !=
          RM::SpatialDiffusionModel::ParallelDiffusionLibrary ||
      options.parallelDiffusionModelId != "constant_kappa" ||
      options.parallelDiffusionParameters.size() != 1 ||
      options.parallelDiffusionParameters.front().name !=
          "kappa_parallel_m2_per_s" ||
      options.parallelDiffusionConfigurationFingerprint.empty())
    return Fail("valid parallel-diffusion block did not resolve: " +
                status.message);

  std::shared_ptr<const RM::RunConfiguration3D> configuration;
  status = RM::RunConfiguration3D::Create(options, &configuration);
  if (!status.ok() || !configuration ||
      configuration->resolved_manifest().find(
          "parallel_diffusion_model=constant_kappa") == std::string::npos ||
      configuration->restart_compatibility_manifest().find(
          options.parallelDiffusionConfigurationFingerprint) ==
              std::string::npos)
    return Fail("parallel model identity is absent from immutable/restart provenance");

  std::string changedParameter = schema5;
  if (!replaceOnce(&changedParameter,
          "kappa_parallel_m2_per_s = 1.25e18",
          "kappa_parallel_m2_per_s = 1.50e18"))
    return Fail("schema-6 fixture lost its coefficient parameter");
  RM::RunConfiguration3DOptions changedOptions;
  std::shared_ptr<const RM::RunConfiguration3D> changedConfiguration;
  if (!RM::ParseConfigurationText(changedParameter, &changedOptions).ok() ||
      !RM::RunConfiguration3D::Create(
          changedOptions, &changedConfiguration).ok() ||
      !changedConfiguration ||
      changedOptions.parallelDiffusionConfigurationFingerprint ==
          options.parallelDiffusionConfigurationFingerprint ||
      changedConfiguration->physics_fingerprint() ==
          configuration->physics_fingerprint() ||
      changedConfiguration->restart_compatibility_manifest() ==
          configuration->restart_compatibility_manifest())
    return Fail("parallel parameter change did not alter restart identity");

  // Runtime::Configure is the serial publication point for the library's
  // active configuration and model-specific function pointer. It must install
  // the same fingerprint that the parser/factory froze.
  RM::Runtime runtime;
  status = runtime.Configure(configuration);
  const auto active = SEP::ParallelDiffusion::GetActiveConfiguration();
  if (!status.ok() || SEP::ParallelDiffusion::ActiveModelFunction == nullptr ||
      SEP::ParallelDiffusion::ActiveParallelDiffusion !=
          SEP::ParallelDiffusion::BoundFunctionForModel(
              SEP::ParallelDiffusion::ModelId::ConstantKappa) ||
      SEP::ParallelDiffusion::ConfigurationFingerprint(active) !=
          options.parallelDiffusionConfigurationFingerprint)
    return Fail("Runtime did not install the validated shared model transactionally");

  std::string missingSection = schema5;
  const std::size_t sectionAt =
      missingSection.find("\n#begin physics: parallel_diffusion\n");
  const std::size_t nextAt =
      missingSection.find("\n#begin physics: shock\n", sectionAt + 1);
  missingSection.erase(sectionAt, nextAt - sectionAt);
  if (RM::ParseConfigurationText(missingSection, &options).ok())
    return Fail("library selector was accepted without physics: parallel_diffusion");

  std::string inactive = schema5;
  if (!replaceOnce(&inactive,
          "spatial_diffusion_model = parallel-diffusion-library",
          "spatial_diffusion_model = mean-free-path") ||
      RM::ParseConfigurationText(inactive, &options).ok())
    return Fail("inactive physics: parallel_diffusion block was accepted");

  std::string focused = schema5;
  if (!replaceOnce(&focused, "transport = parker", "transport = focused-diffusion"))
    return Fail("fixture lost its Parker transport spelling");
  if (RM::ParseConfigurationText(focused, &options).ok())
    return Fail("focused mover accepted the Parker spatial library");

  std::string unknown = schema5;
  unknown.insert(unknown.find("#end physics\n\n#begin physics: shock"),
                 "unknown_parameter = 7\n");
  if (RM::ParseConfigurationText(unknown, &options).ok())
    return Fail("shared model reader accepted an unknown parameter");

  std::string missingTurbulence = schema5;
  if (!replaceOnce(&missingTurbulence,
          "model = constant_kappa\nkappa_parallel_m2_per_s = 1.25e18",
          "model = qlt_slab_inertial") ||
      RM::ParseConfigurationText(missingTurbulence, &options).ok())
    return Fail("srcSEP3D accepted a spectral model without slab state");

  std::string rigidity = schema5;
  if (!replaceOnce(&rigidity,
          "model = constant_kappa\nkappa_parallel_m2_per_s = 1.25e18",
          "model = power_law_lambda\n"
          "lambda0_m = 1e10\n"
          "independent_variable = rigidity\n"
          "rigidity0_V = 1e9\n"
          "independent_exponent = 0.3333333333333333") ||
      !RM::ParseConfigurationText(rigidity, &options).ok())
    return Fail("case-sensitive SI rigidity key was not forwarded intact");
  std::string wrongCase = rigidity;
  if (!replaceOnce(&wrongCase, "rigidity0_V", "rigidity0_v") ||
      RM::ParseConfigurationText(wrongCase, &options).ok())
    return Fail("srcSEP3D lower-cased or accepted a misspelled library SI key");

  std::string nucleon = schema5;
  if (!replaceOnce(&nucleon,
          "model = constant_kappa\nkappa_parallel_m2_per_s = 1.25e18",
          "model = power_law_lambda\n"
          "lambda0_m = 1e10\n"
          "independent_variable = energy_per_nucleon\n"
          "energy_per_nucleon0_J = 1e-12\n"
          "independent_exponent = 0.3333333333333333") ||
      RM::ParseConfigurationText(nucleon, &options).ok())
    return Fail("srcSEP3D inferred an unavailable nucleon count");

  // The block body belongs to the library parser; '!' opens a comment (as
  // everywhere in the generic format), and its diagnostics carry the deck's
  // own line numbers.  The expected line is found by searching the fixture
  // text, not taken from the parser.
  std::string bangComment = schema5;
  if (!replaceOnce(&bangComment, "kappa_parallel_m2_per_s = 1.25e18",
                   "kappa_parallel_m2_per_s = 1.25e18 ! Fortran-style note") ||
      !RM::ParseConfigurationText(bangComment, &options).ok() ||
      options.parallelDiffusionParameters.size() != 1 ||
      options.parallelDiffusionParameters.front().value != "1.25e18")
    return Fail("library section grammar did not strip a '!' comment");
  std::string badValue = schema5;
  if (!replaceOnce(&badValue, "kappa_parallel_m2_per_s = 1.25e18",
                   "kappa_parallel_m2_per_s = -1"))
    return Fail("schema-6 fixture lost its coefficient parameter");
  const std::size_t badAt = badValue.find("kappa_parallel_m2_per_s = -1");
  const std::size_t badLine = 1 + static_cast<std::size_t>(
      std::count(badValue.begin(), badValue.begin() + badAt, '\n'));
  status = RM::ParseConfigurationText(badValue, &options);
  if (status.ok() ||
      status.message.find("line " + std::to_string(badLine) + ":") ==
          std::string::npos)
    return Fail("library diagnostic lost the deck line number: " +
                status.message);

  return Pass("the parallel-diffusion block selects the Parker mover family, calls the strict shared parser, freezes restart identity, installs one Parker-only active model, and rejects unavailable host state");
}

// CFG3D18 -- the shared srcSEP/srcSEP3D vocabulary (content schema 6,
// src/models/sep_common/sep_input_layout.h).
Result RunCFG3D18() {
  namespace Layout = SEP::InputLayout;
  // (1) The 3-D-only SWCME keys that srcSEP rejects are exactly the keys the
  // canonical SWCME3D resolver adds to SWCME1D.
  const std::set<std::string> table(Layout::SwcmeThreeDimensionalKeys().begin(),
                                    Layout::SwcmeThreeDimensionalKeys().end());
  if (table != swcme::input3d::detail::geometry_keys())
    return Fail("sep_input_layout.h 3-D SWCME key list differs from "
                "swcme::input3d::detail::geometry_keys()");

  // (2) D3: a srcSEP deck is read by srcSEP3D's parser, which stops at the
  // first srcSEP-only setting with the typed "not supported" diagnosis
  // instead of ignoring it or calling it unknown.
  const std::string srcSepDeck = ReadText("../srcSEP/examples/sep_parker_mesh.in");
  RM::RunConfiguration3DOptions options;
  const SEP3D::Core::Status srcSepStatus =
      RM::ParseConfigurationText(srcSepDeck, &options);
  if (srcSepDeck.empty() || srcSepStatus.ok() ||
      srcSepStatus.message.find("not supported by srcSEP3D") == std::string::npos)
    return Fail("a srcSEP schema-6 deck did not stop with a typed srcSEP3D "
                "not-supported diagnosis: " + srcSepStatus.message);

  // (3) Every key the table assigns to srcSEP3D is understood by the parser:
  // a key it does not know would be reported as "unknown configuration key".
  // Insert each into the complete fixture with a deliberately invalid value;
  // the diagnosis must then be a value or rule error, never "unknown".
  const std::string complete = CompleteInput();
  for (const Layout::KeyRule& rule : Layout::kKeys) {
    if ((rule.support & Layout::kSrcSep3d) == 0) continue;
    const std::string block = rule.block;
    const std::string qualified = block + "." + rule.key;
    std::string anchor;
    if (block == "observer.*") anchor = "#begin observer: earth\n";
    else if (block == "run") anchor = "#begin run\n";
    else if (block == "output") anchor = "#begin output\n";
    else continue;  // nested blocks are covered by the complete fixture
    std::string text = complete;
    if (!ReplaceFirst(&text, anchor, anchor + "  " + rule.key + " = @invalid@\n"))
      return Fail("fixture lost anchor " + anchor);
    const SEP3D::Core::Status status = RM::ParseConfigurationText(text, &options);
    if (status.ok() || status.message.find("unknown configuration key") !=
                           std::string::npos)
      return Fail("table key '" + qualified + "' is not understood by the "
                  "srcSEP3D parser: " + status.message);
  }

  // (4) Every maintained srcSEP3D example deck parses (the reduced-front
  // decks resolve their #included events and assets).
  const char* examples[] = {
      "examples/sep3d_analytic_parker.in",
      "examples/sep3d_analytic_parker_active_tube.in",
      "examples/sep3d_analytic_parker_corner_sphere.in",
      "examples/sep3d_swcme_20rs_1au.in",
      "examples/sep3d_swcme_sphere_mesh_background_20rs_1au.in",
      "examples/sep3d_swcme_sse_mesh_background_20rs_1au.in",
      "examples/application-input/amps.in",
      "examples/application-input/amps-injection-smoke.in",
      "examples/shock-front/positive_1au.in",
      "examples/shock-front/handoff_smoke.in",
      "examples/shock-front/corona_to_1au.in"};
  for (const char* example : examples) {
    const SEP3D::Core::Status status = RM::LoadConfigurationFile(example, &options);
    if (!status.ok())
      return Fail(std::string("example ") + example + " does not parse: " +
                  status.message);
  }
  return Pass("schema-6 table matches the SWCME3D key set, srcSEP decks stop "
              "with typed srcSEP3D errors, every srcSEP3D table key is "
              "understood, and every example deck parses");
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
      make("CFG3D01", "Input and CLI", "C01 schema 6, CLI, layout-independent identity, D1-D4 errors, and dry-run.", RunCFG3D01),
      make("CFG3D02", "Typed contracts", "C02 complete groups and fingerprint classes.", RunCFG3D02),
      make("CFG3D03", "Domain contract", "C03 presets, containment, and boundary direction.", RunCFG3D03),
      make("CFG3D04", "Parker geometry", "C04 shared polarity-independent geometry.", RunCFG3D04),
      make("CFG3D05", "Mesh preflight", "C05 composite refinement and memory planning.", RunCFG3D05),
      make("CFG3D06", "Finite field line", "Required finite Parker-line fields and fail-closed consistency checks.", RunCFG3D06),
      make("CFG3D07", "Compiled-species binding", "Complete generated AMPS table and fingerprint contract.", RunCFG3D07),
      make("CFG3D08", "Complete initialization", "Canonical SWCME, D4 Parker derivation, and exact per-step source contract.", RunCFG3D08),
      make("CFG3D09", "Background/turbulence selection", "Named prescribed slopes and reserved Python source.", RunCFG3D09),
      make("CFG3D10", "CME/Parker start linkage", "Canonical launch-apex linkage and fail-closed geometry checks.", RunCFG3D10),
      make("CFG3D11", "Transport/control schema", "Active corridor, population limits, mover coefficients, and dry-run output.", RunCFG3D11),
      make("CFG3D12", "Corner/sphere geometry", "Endpoint extent, input controls, identity and corner validation.", RunCFG3D12),
      make("CFG3D13", "Unified runtime input", "Generic block/include grammar, derived particle numerics keys, CLI spellings, immutable commit, and file:line diagnostics.", RunCFG3D13),
      make("CFG3D14", "Derived particle numerics", "Mesh CFL step and incident-flux per-species weight equations.", RunCFG3D14),
      make("CFG3D15", "Parsed reduced model", "Maintained unified input initializes the real reduced provider and source normalization.", RunCFG3D15),
      make("CFG3D16", "Reduced-front particle source", "Accepted-face rate sum, Poisson timing, triangular position, local spectrum, determinism, and reserved weighting.", RunCFG3D16),
      make("CFG3D17", "Parallel-diffusion binding", "Library block parser, provenance, Runtime installation, Parker-only ownership, and unavailable-state gates.", RunCFG3D17),
      make("CFG3D18", "Shared SEP vocabulary", "Schema-6 table vs SWCME3D keys, srcSEP decks give typed srcSEP3D errors, every example parses.", RunCFG3D18),
  };
}

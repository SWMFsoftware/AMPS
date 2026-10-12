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
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
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

// A deck with every '#include "FILE"' line replaced by that file's text
// (recursively, relative to the including file), so tests can edit any key
// of a multi-file deck in one string.
std::string InlineIncludes(const std::string& path) {
  const std::string text = ReadText(path);
  const std::string directory =
      std::filesystem::path(path).parent_path().string();
  std::istringstream input(text);
  std::string line, out;
  while (std::getline(input, line)) {
    const std::size_t at = line.find("#include \"");
    if (at != std::string::npos && line.find_first_not_of(" \t") == at) {
      const std::size_t begin = at + 10;
      const std::size_t end = line.find('"', begin);
      out += InlineIncludes(directory + "/" + line.substr(begin, end - begin));
      continue;
    }
    out += line + "\n";
  }
  return out;
}

// The complete schema-6 fixture is the maintained reduced-front application
// deck examples/application-input/amps.in (with sep3d.in and
// reduced-shock-surface.in inlined): every srcSEP3D key is explicit there, so
// the parser's completeness rules are exercised against the deck users
// actually start from.  The SWCME decks are not fixtures any more: the SWCME
// model is not ready and every physics: swcme block is rejected (CFG3D20).
// The inlined text is parsed from the srcSEP3D directory, so the reduced
// model's relative asset path is re-rooted there.
std::string CompleteInput() {
  std::string text = InlineIncludes("examples/application-input/amps.in");
  const std::string asset = "../shock-front/assets/";
  const std::size_t at = text.find(asset);
  if (at != std::string::npos)
    text.replace(at, asset.size(), "examples/shock-front/assets/");
  return text;
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
  const std::string retiredSpeciesField = CompleteInput() +
      "#begin species\n  mass_kg = 1.67262192369e-27\n#end species\n";
  if (RM::ParseConfigurationText(retiredSpeciesField, &parsed).ok())
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
  // SWCME is not ready (user decision 2026-10-12): a physics: swcme block and
  // source.shock.model = swcme stop with the typed diagnosis, and the source
  // model must be written explicitly.
  std::string retiredSchema = complete, srcSepKey = complete,
      srcSepTurbulence = complete, swcmeModel = complete, noModel = complete,
      retiredMover = complete, oldBlock = complete, fieldLineObserver = complete;
  const std::string swcmeBlock = complete +
      "#begin physics: swcme\n  ambient.wind_speed = 400 km/s\n#end physics\n";
  const bool edited =
      ReplaceFirst(&retiredSchema, "schema_version = 6", "schema_version = 4") &&
      ReplaceFirst(&srcSepKey, "kind = fixed-cartesian",
                   "kind = fixed-cartesian\n  heliocentric_radius_m = 1e11") &&
      ReplaceFirst(&srcSepTurbulence, "authority = prescribed",
                   "authority = prescribed\n    cascade = true") &&
      ReplaceFirst(&swcmeModel, "model = reduced-shock-front", "model = swcme") &&
      ReplaceFirst(&noModel, "  model = reduced-shock-front\n", "") &&
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
      {"SWCME block (not ready)", swcmeBlock,
       "the SWCME model is not ready yet to be used"},
      {"SWCME source model (not ready)", swcmeModel,
       "the SWCME model is not ready yet to be used"},
      {"implicit source model", noModel,
       "missing required key 'source.shock.model'"},
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
      "examples/application-input/amps.in", &exampleOptions);
  if (!exampleStatus.ok()) {
    return Fail("the annotated production example is invalid: " +
                exampleStatus.message);
  }
  std::shared_ptr<const RM::RunConfiguration3D> example;
  if (!RM::RunConfiguration3D::Create(exampleOptions, &example).ok() ||
      !RM::BuildDryRunSummary(*example, &summary).ok()) {
    return Fail("the annotated production example failed resource preflight");
  }
  if (exampleOptions.observers.size() != 1 ||
      exampleOptions.observers[0].allCompiledSpecies ||
      exampleOptions.observers[0].species != std::vector<int>{0} ||
      exampleOptions.source.model != RM::ParticleSourceModel::ReducedShockFront) {
    return Fail("the production example lost its explicit observer species "
                "list or its reduced shock-front source model");
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
  return Pass("schema-6 input, layout-independent identity, CLI normalization, D1-D3 and SWCME-not-ready rejections, early errors, and allocation-free dry-run passed");
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
  // The finite Parker-line definition is a required part of schema 6.  The
  // fixture writes the documented "derive" zeros; explicit values are set
  // here so the exact parse and the inner-boundary rule are exercised.
  RM::RunConfiguration3DOptions options;
  std::string complete = CompleteInput();
  if (!ReplaceFirst(&complete, "initial_x_m = 0\n", "initial_x_m = 1.3914e10\n") ||
      !ReplaceFirst(&complete, "length_m = 0\n", "length_m = 2.0e11\n") ||
      !ReplaceFirst(&complete, "point_count = 0\n", "point_count = 4001\n"))
    return Fail("fixture lost its field_line anchors");
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

// CFG3D08 -- the SWCME model is not ready (user decision 2026-10-12).
//
// The six SWCME example decks are kept as rejected fixtures: each must stop
// at its '#begin physics: swcme' line with the typed "not ready" diagnosis,
// whatever SWCME would have been used for (analytic-Parker decks with the
// SWCME shock source, SWCME background/propagation decks).  The maintained
// reduced-front deck still parses and selects the reduced shock-front source
// explicitly; an enabled source without the reduced front is rejected.
Result RunCFG3D08() {
  const char* decks[] = {
      "examples/sep3d_analytic_parker.in",
      "examples/sep3d_analytic_parker_active_tube.in",
      "examples/sep3d_analytic_parker_corner_sphere.in",
      "examples/sep3d_swcme_20rs_1au.in",
      "examples/sep3d_swcme_sphere_mesh_background_20rs_1au.in",
      "examples/sep3d_swcme_sse_mesh_background_20rs_1au.in"};
  for (const char* deck : decks) {
    // The expected line is found independently: the first '#begin physics:
    // swcme' of the deck file.
    std::istringstream lines(ReadText(deck));
    std::string line;
    std::size_t number = 0, swcmeLine = 0;
    while (std::getline(lines, line)) {
      ++number;
      if (line.rfind("#begin physics: swcme", 0) == 0) { swcmeLine = number; break; }
    }
    RM::RunConfiguration3DOptions options;
    const SEP3D::Core::Status status = RM::LoadConfigurationFile(deck, &options);
    const std::string where = std::string(deck) + ":" + std::to_string(swcmeLine) + ":";
    if (swcmeLine == 0 || status.ok() ||
        status.message.find("the SWCME model is not ready yet to be used") ==
            std::string::npos ||
        status.message.find(where) == std::string::npos)
      return Fail(std::string("SWCME deck ") + deck + " did not stop at " +
                  where + " with the not-ready diagnosis: " + status.message);
  }

  RM::RunConfiguration3DOptions options;
  const std::string complete = CompleteInput();
  if (!RM::ParseConfigurationText(complete, &options).ok() ||
      options.source.model != RM::ParticleSourceModel::ReducedShockFront ||
      !options.source.enabled || options.injectionCadenceSteps != 1)
    return Fail("the reduced-front fixture lost its explicit source model");
  std::string noFront = complete;
  const std::size_t frontBegin = noFront.find("#begin physics: shock_front");
  const std::size_t frontEnd = noFront.find("#end physics", frontBegin);
  if (frontBegin == std::string::npos || frontEnd == std::string::npos)
    return Fail("fixture lost its physics: shock_front block");
  noFront.erase(frontBegin, frontEnd + std::string("#end physics").size() - frontBegin);
  if (RM::ParseConfigurationText(noFront, &options).ok())
    return Fail("an enabled source without the reduced shock front was accepted");
  std::string skippedStep = complete;
  if (!ReplaceFirst(&skippedStep, "injection_cadence_steps = 1",
                    "injection_cadence_steps = 2") ||
      RM::ParseConfigurationText(skippedStep, &options).ok())
    return Fail("a source that skips simulation steps was accepted");
  return Pass("the six SWCME decks stop at their physics: swcme line with the "
              "not-ready diagnosis; the reduced-front source is explicit and "
              "requires its front");
}

Result RunCFG3D09() {
  const std::string complete = CompleteInput();
  if (complete.empty())
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

  // The future Python interpolation background is a typed reserved feature.
  // The reduced-front fixture's model_id and physics: shock_front belong to
  // provider = runtime-model, so the reservation is exercised at the typed
  // factory, which every parsed deck also passes.
  RM::RunConfiguration3DOptions python;
  if (!RM::ParseConfigurationText(complete, &python).ok())
    return Fail("fixture did not parse");
  python.background = RM::BackgroundAuthority::PythonInterpolator;
  std::shared_ptr<const RM::RunConfiguration3D> reserved;
  const SEP3D::Core::Status pythonStatus =
      RM::RunConfiguration3D::Create(python, &reserved);
  if (pythonStatus.code != SEP3D::Core::StatusCode::ReservedFeature ||
      pythonStatus.message.find("Python") == std::string::npos) {
    return Fail("reserved Python background did not fail with its typed status");
  }

  return Pass(
      "input selects validated spectral/amplitude turbulence models and reserves the future Python background explicitly");
}

// CFG3D10 -- Parker start linkage without SWCME.
//
// start_mode = cme-launch-point freezes the canonical SWCME launch apex as
// the Parker start.  With the SWCME model not ready (CFG3D08) that apex does
// not exist, so the linked mode must be rejected rather than silently fall
// back to an explicit point; the explicit mode of the reduced-front fixture
// carries no CME linkage into the fingerprint.
Result RunCFG3D10() {
  const std::string complete = CompleteInput();
  RM::RunConfiguration3DOptions linked;
  if (!RM::ParseConfigurationText(complete, &linked).ok() ||
      linked.cmeLaunchPointResolved ||
      linked.parkerSpiralStartMode != RM::ParkerSpiralStartMode::Explicit)
    return Fail("the explicit Parker start carried an implicit CME linkage");
  std::string cmeLinked = complete;
  if (!ReplaceFirst(&cmeLinked, "start_mode = explicit",
                    "start_mode = cme-launch-point"))
    return Fail("fixture lost its start mode");
  const SEP3D::Core::Status status = RM::ParseConfigurationText(cmeLinked, &linked);
  if (status.ok())
    return Fail("a CME-linked Parker start was accepted without an SWCME launch apex");
  return Pass("explicit Parker start has no CME linkage; cme-launch-point is "
              "rejected without SWCME (" + status.message + ")");
}

Result RunCFG3D11() {
  const std::string complete = CompleteInput();
  if (complete.empty())
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

  // The shipped active-tube example used the SWCME shock and is now a
  // rejected fixture (CFG3D08); the corridor cases below place the observer
  // on the exact finite Parker line of the reduced-front fixture instead.
  SEP3D::Core::ParkerSpiralGeometry fixtureGeometry;
  fixtureGeometry.sourceRadiusM = options.innerRadiusM;
  fixtureGeometry.sourceLongitudeRad = options.tubeLongitudeRad;
  fixtureGeometry.sourceColatitudeRad = options.tubeColatitudeRad;
  fixtureGeometry.solarWindSpeedMPerS = options.parker.solarWindSpeedMPerS;
  fixtureGeometry.solarRotationRateRadPerS =
      options.parker.solarRotationRateRadPerS;
  fixtureGeometry.rotationAxis = options.parker.rotationAxis;
  if (options.observers.size() != 1)
    return Fail("fixture should have exactly one observer");
  const double observerRadiusM =
      (options.observers[0].positionM - options.coordinateOriginM).Norm();
  const SEP3D::Core::Vec3 onLine = options.coordinateOriginM +
      SEP3D::Core::ParkerCurvePoint(observerRadiusM, fixtureGeometry);
  std::ostringstream onLineText;
  onLineText << std::setprecision(17) << "position_x_m = " << onLine.x
             << "\n  position_y_m = " << onLine.y
             << "\n  position_z_m = " << onLine.z;
  const std::string fixtureObserver =
      "position_x_m = 3.7399467675e10\n  position_y_m = 0\n  position_z_m = 0";

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
  // The corridor rules compare against the refinement tube, which the
  // fixture leaves disabled; enable it for the corridor variants.
  if (!replaceOnce(&corridor, "enabled = false\n    source_longitude_rad",
                   "enabled = true\n    source_longitude_rad") ||
      !replaceOnce(&corridor, "mode = full-domain", "mode = parker-tube") ||
      !replaceOnce(&corridor, "radius_at_reference_m = 0\n    radius_mode = constant-angular-width\n    buffer_blocks = 0",
                   "radius_at_reference_m = 7.479893535e9\n    radius_mode = constant-angular-width\n    buffer_blocks = 1") ||
      !replaceOnce(&corridor, fixtureObserver, onLineText.str()))
    return Fail("fixture lost an active-corridor anchor");
  const SEP3D::Core::Status corridorStatus =
      RM::ParseConfigurationText(corridor, &options);
  if (!corridorStatus.ok() ||
      options.activeRegion != RM::ActiveRegionMode::ParkerTube ||
      options.activeTubeBufferBlocks != 1) {
    return Fail("valid Parker active corridor did not parse: " +
                corridorStatus.message);
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
  // 90 degrees from the line at the same radius: far outside the corridor
  // plus its one-block buffer.
  if (!replaceOnce(&disconnectedObserver, onLineText.str(),
                   "position_x_m = 0\n  position_y_m = 3.7399467675e10\n  position_z_m = 0") ||
      RM::ParseConfigurationText(disconnectedObserver, &options).ok()) {
    return Fail("observer outside the stationary active corridor was accepted");
  }
  std::string truncatedBeforeObserver = corridor;
  if (!replaceOnce(&truncatedBeforeObserver, "length_m = 0\n",
                   "length_m = 1.0e10\n") ||
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

// The corner/sphere geometry of the former SWCME example
// sep3d_analytic_parker_corner_sphere.in (now a rejected fixture, CFG3D08),
// written into the reduced-front fixture: photospheric inner boundary, x-y
// corner cube with the outer radius at the field-line end point, the
// photosphere-anchored solar refinement and an active solar sphere.  Only the
// geometry keys change; the values are those of the former example.
std::string CornerSphereInput(std::string* why = nullptr) {
  std::string text = CompleteInput();
  const bool ok =
      ReplaceFirst(&text, "inner_radius_m = 1.3914e10", "inner_radius_m = 6.957e8") &&
      ReplaceFirst(&text, "outer_radius_mode = preset\n    outer_radius_m = 1.495978707e11",
                   "outer_radius_mode = field-line-endpoint\n"
                   "    box_geometry = field-line-xy-corner-cube\n"
                   "    corner_direction_x = 0\n    corner_direction_y = 0\n"
                   "    corner_direction_z = 0\n    corner_margin_m = 1.495978707e9\n"
                   "    outer_radius_m = 1.6455765777e11") &&
      ReplaceFirst(&text, "initial_x_m = 0\n", "initial_x_m = 6.957e8\n") &&
      ReplaceFirst(&text, "initial_z_m = 0\n",
                   "initial_z_m = 0\n    end_radius_m = 1.6455765777e11\n") &&
      ReplaceFirst(&text, "point_count = 0\n", "point_count = 4001\n") &&
      ReplaceFirst(&text, "#begin refinement: solar\n",
                   "#begin refinement: solar\n    anchor = photosphere\n") &&
      ReplaceFirst(&text, "surface_cell_size_m = 1.495978707e9",
                   "surface_cell_size_m = 3.4785e8") &&
      ReplaceFirst(&text, "transition_outer_radius_m = 3.7399467675e10",
                   "transition_outer_radius_m = 1.3914e10") &&
      // The photospheric cells need the former example's finer mesh limits.
      ReplaceFirst(&text, "minimum_cell_size_m = 1.495978707e9",
                   "minimum_cell_size_m = 3.4785e8") &&
      ReplaceFirst(&text, "maximum_level = 7", "maximum_level = 9") &&
      ReplaceFirst(&text, "memory_budget_bytes = 4294967296",
                   "memory_budget_bytes = 68719476736");
  if (!ok) { if (why) *why = "a geometry anchor is missing"; return std::string(); }
  // An active solar sphere connects to a stationary Parker corridor, which
  // needs the refinement tube and an observer on the finite Parker line.  The
  // line depends only on its source radius (the photosphere here), the tube
  // direction and the Parker wind/rotation, which the plain fixture supplies;
  // the observer is placed exactly on it (same construction as CFG3D11).
  RM::RunConfiguration3DOptions geometryOnly;
  const SEP3D::Core::Status first =
      RM::ParseConfigurationText(CompleteInput(), &geometryOnly);
  if (!first.ok() || geometryOnly.observers.size() != 1) {
    if (why) *why = "geometry pass: " + first.message;
    return std::string();
  }
  SEP3D::Core::ParkerSpiralGeometry line;
  line.sourceRadiusM = 6.957e8;
  line.sourceLongitudeRad = geometryOnly.tubeLongitudeRad;
  line.sourceColatitudeRad = geometryOnly.tubeColatitudeRad;
  line.solarWindSpeedMPerS = geometryOnly.parker.solarWindSpeedMPerS;
  line.solarRotationRateRadPerS = geometryOnly.parker.solarRotationRateRadPerS;
  line.rotationAxis = geometryOnly.parker.rotationAxis;
  const SEP3D::Core::Vec3 onLine = geometryOnly.coordinateOriginM +
      SEP3D::Core::ParkerCurvePoint(
          (geometryOnly.observers[0].positionM -
           geometryOnly.coordinateOriginM).Norm(), line);
  std::ostringstream observer;
  observer << std::setprecision(17) << "position_x_m = " << onLine.x
           << "\n  position_y_m = " << onLine.y
           << "\n  position_z_m = " << onLine.z;
  const bool corridor =
      ReplaceFirst(&text, "enabled = false\n    source_longitude_rad",
                   "enabled = true\n    source_longitude_rad") &&
      ReplaceFirst(&text, "#begin active_region\n",
                   "#begin active_region\n    solar_sphere_radius_m = 2.0871e10\n") &&
      ReplaceFirst(&text, "mode = full-domain", "mode = parker-tube") &&
      ReplaceFirst(&text, "radius_at_reference_m = 0\n    radius_mode = constant-angular-width\n    buffer_blocks = 0",
                   "radius_at_reference_m = 7.479893535e9\n    radius_mode = constant-angular-width\n    buffer_blocks = 1") &&
      ReplaceFirst(&text,
                   "position_x_m = 3.7399467675e10\n  position_y_m = 0\n  position_z_m = 0",
                   observer.str());
  if (!corridor && why) *why = "a corridor anchor is missing";
  return corridor ? text : std::string();
}

Result RunCFG3D12() {
  // The corner/sphere geometry selects the x-y corner box; the three-axis
  // corner mode is the same deck with the original geometry name.
  std::string why;
  const std::string xyDeck = CornerSphereInput(&why);
  std::string input = xyDeck;
  if (xyDeck.empty() ||
      !ReplaceFirst(&input, "box_geometry = field-line-xy-corner-cube",
                    "box_geometry = field-line-corner-cube"))
    return Fail("corner/sphere fixture unavailable: " + why);
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
      "  model = reduced-shock-front\n"
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
  const std::string fixture = CompleteInput();
  if (fixture.empty())
    return Fail("could not read the schema-6 fixture");

  auto replaceOnce = [](std::string* text, const std::string& from,
                        const std::string& to) {
    const std::size_t at = text == nullptr ? std::string::npos : text->find(from);
    if (at == std::string::npos) return false;
    text->replace(at, from.size(), to);
    return true;
  };
  std::string schema5 = fixture;
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
    if (block == "observer.*") anchor = "#begin observer: default\n";
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

  // (4) Every maintained (non-SWCME) srcSEP3D example deck parses (the
  // reduced-front decks resolve their #included events and assets).
  // (The six SWCME decks are rejected fixtures, checked by CFG3D08.)
  const char* examples[] = {
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

// CFG3D19 -- the point-sampling block ('#begin sampling: points').
//
// The complete fixture plus a sampling block whose points come from a point
// block and a point file (relative to the input file).  Checks: the typed
// records in input order (block points first, then file rows; an unnamed row
// is named by its zero-based file index), MeV channels and scale, immutable
// construction, the fingerprint change, and the fail-closed rules -- missing
// sub-block keys, a point file that does not exist, repeated point names, a
// point outside the radial domain, an unknown product, and a repeated
// species block.
Result RunCFG3D19() {
  namespace fs = std::filesystem;
  const fs::path directory = fs::temp_directory_path() /
      ("srcsep3d-cfg3d19-" + std::to_string(
          std::chrono::high_resolution_clock::now().time_since_epoch().count()));
  fs::create_directories(directory);
  struct Cleanup {
    fs::path path;
    ~Cleanup() { std::error_code e; fs::remove_all(path, e); }
  } cleanup{directory};
  {
    std::ofstream points(directory / "points.txt");
    points << "! x y z [name], metres\n1.0e11 0 0 far\n9.0e10 1.0e9 0\n";
  }
  const std::string complete = CompleteInput();
  const std::string block =
      "#begin sampling: points\n"
      "  products = flux, density, spectrum, fluence\n"
      "  point_file = points.txt\n"
      "  #begin species: H_PLUS\n"
      "    minimum_energy_mev = 1\n"
      "    maximum_energy_mev = 1000\n"
      "    energy_intervals = 30\n"
      "    energy_scale = log\n"
      "  #end species\n"
      "  #begin point: near\n"
      "    position_x_m = 7.5e10\n"
      "    position_y_m = 0\n"
      "    position_z_m = 0\n"
      "  #end point\n"
      "#end sampling\n";
  auto parse = [&](const std::string& text, RM::RunConfiguration3DOptions* o) {
    return RM::ParseConfigurationText(text, o, directory.string());
  };
  RM::RunConfiguration3DOptions options;
  const SEP3D::Core::Status status = parse(complete + block, &options);
  if (!status.ok())
    return Fail("a complete sampling block did not parse: " + status.message);
  const RM::PointSamplingOptions& sampling = options.pointSampling;
  if (!sampling.enabled || sampling.products != "flux, density, spectrum, fluence" ||
      sampling.species.size() != 1 || sampling.species[0].symbol != "H_PLUS" ||
      sampling.species[0].minimumEnergyMeV != 1.0 ||
      sampling.species[0].maximumEnergyMeV != 1000.0 ||
      sampling.species[0].energyIntervals != 30 ||
      sampling.species[0].energyScale != RM::EnergyChannelSpacing::Logarithmic)
    return Fail("sampling block or species channels parsed incorrectly");
  if (sampling.points.size() != 3 || sampling.points[0].name != "near" ||
      sampling.points[0].positionM.x != 7.5e10 ||
      sampling.points[1].name != "far" || sampling.points[1].positionM.x != 1.0e11 ||
      sampling.points[2].name != "1" || sampling.points[2].positionM.y != 1.0e9 ||
      sampling.points[1].origin.find("points.txt:2") == std::string::npos)
    return Fail("point order, names, positions or origins are wrong");
  std::shared_ptr<const RM::RunConfiguration3D> withSampling, without;
  RM::RunConfiguration3DOptions plain;
  if (!RM::RunConfiguration3D::Create(options, &withSampling).ok() ||
      !parse(complete, &plain).ok() ||
      !RM::RunConfiguration3D::Create(plain, &without).ok())
    return Fail("immutable construction failed");
  if (withSampling->physics_fingerprint() == without->physics_fingerprint())
    return Fail("the sampling definition does not enter the fingerprint");

  struct Negative { const char* label; std::string from, to, expect; };
  const std::vector<Negative> negatives = {
      {"missing channel key", "    energy_intervals = 30\n", "",
       "missing required key 'sampling.points.species.H_PLUS.energy_intervals'"},
      {"missing point coordinate", "    position_z_m = 0\n", "",
       "missing required key 'sampling.points.point.near.position_z_m'"},
      {"missing point file", "point_file = points.txt", "point_file = absent.txt",
       "cannot read sampling.points.point_file"},
      {"repeated point name", "#begin point: near", "#begin point: far",
       "sampling point name 'far'"},
      {"point outside the domain", "position_x_m = 7.5e10",
       "position_x_m = 1.0e9", "outside the heliocentric domain"},
      {"unknown product", "fluence\n", "intensity\n", "unknown product"},
      {"bad scale", "energy_scale = log", "energy_scale = logarithmic",
       "invalid value for 'sampling.points.species.H_PLUS.energy_scale'"},
  };
  for (const Negative& negative : negatives) {
    std::string text = block;
    if (!ReplaceFirst(&text, negative.from, negative.to))
      return Fail(std::string("fixture lost anchor for ") + negative.label);
    RM::RunConfiguration3DOptions bad;
    SEP3D::Core::Status result = parse(complete + text, &bad);
    std::shared_ptr<const RM::RunConfiguration3D> created;
    if (result.ok()) result = RM::RunConfiguration3D::Create(bad, &created);
    if (result.ok() || result.message.find(negative.expect) == std::string::npos)
      return Fail(std::string(negative.label) + " was not rejected as expected: " +
                  result.message);
  }
  std::string repeated = block;
  ReplaceFirst(&repeated, "  #end species\n",
               "  #end species\n  #begin species: H_PLUS\n  #end species\n");
  RM::RunConfiguration3DOptions bad;
  const SEP3D::Core::Status repeatedStatus = parse(complete + repeated, &bad);
  if (repeatedStatus.ok())
    return Fail("a repeated species block was accepted");
  return Pass("sampling block, point file, MeV channels, fingerprint and seven "
              "fail-closed rules");
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
      make("CFG3D08", "SWCME not ready", "The six SWCME decks stop at physics: swcme with the not-ready diagnosis; the reduced-front source is explicit and needs its front.", RunCFG3D08),
      make("CFG3D09", "Background/turbulence selection", "Named prescribed slopes and reserved Python source.", RunCFG3D09),
      make("CFG3D10", "Parker start without SWCME", "Explicit start has no CME linkage; cme-launch-point is rejected without an SWCME launch apex.", RunCFG3D10),
      make("CFG3D11", "Transport/control schema", "Active corridor, population limits, mover coefficients, and dry-run output.", RunCFG3D11),
      make("CFG3D12", "Corner/sphere geometry", "Endpoint extent, input controls, identity and corner validation.", RunCFG3D12),
      make("CFG3D13", "Unified runtime input", "Generic block/include grammar, derived particle numerics keys, CLI spellings, immutable commit, and file:line diagnostics.", RunCFG3D13),
      make("CFG3D14", "Derived particle numerics", "Mesh CFL step and incident-flux per-species weight equations.", RunCFG3D14),
      make("CFG3D15", "Parsed reduced model", "Maintained unified input initializes the real reduced provider and source normalization.", RunCFG3D15),
      make("CFG3D16", "Reduced-front particle source", "Accepted-face rate sum, Poisson timing, triangular position, local spectrum, determinism, and reserved weighting.", RunCFG3D16),
      make("CFG3D17", "Parallel-diffusion binding", "Library block parser, provenance, Runtime installation, Parker-only ownership, and unavailable-state gates.", RunCFG3D17),
      make("CFG3D18", "Shared SEP vocabulary", "Schema-6 table vs SWCME3D keys, srcSEP decks give typed srcSEP3D errors, every example parses.", RunCFG3D18),
      make("CFG3D19", "Point-sampling input", "sampling: points block, point file, MeV channels, fingerprint and fail-closed rules.", RunCFG3D19),
  };
}

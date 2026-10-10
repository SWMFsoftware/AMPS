// PDB01-PDB05: srcSEP schema-4 input: binding of the shared parallel-diffusion
// library (PDB01-PDB04) and the optional [run] particle_mover key (PDB05).
//
// Built without AMPS, PIC, or MPI by test/run_parallel_diffusion_binding_tests.sh
// from the production sources: util/sep_initialization.cpp (schema-4 section
// location), adapters/parallel_diffusion_adapter.cpp (selection cross-check,
// host-capability declaration, installation, SI state packing), the canonical
// sep_common registry, and the library objects.  The PIC sampling wrapper in
// coefficient_providers.cpp needs AMPS and is exercised only by a native build.
//
// Expected coefficients come from independent closed forms evaluated here
// (kappa = v lambda / 3 with v from p and m; lambda0 (r/r0)^alpha), never from
// the production evaluator.  The library's active model is process-global and
// cannot be uninstalled, so tests that must observe "nothing installed" run
// before the installing test; each test otherwise owns its configuration.

#include "../../adapters/parallel_diffusion_adapter.h"
#include "../../util/sep_initialization.h"

#include "parallel_diffusion/parallel_diffusion.h"
#include "sep_coefficient_registry.h"
#include "sep_test_registry.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

namespace Coefficient = SEP::Transport::Coefficient;
namespace Binding = SEP::ParallelDiffusionBinding;
namespace PD = SEP::ParallelDiffusion;
using Coefficient::SpatialKind;

std::string gExamplePath;

SEP::Testing::Result Pass(const std::string& message) {
  SEP::Testing::Result result;
  result.status = SEP::Testing::Status::Pass;
  result.message = message;
  result.metrics.push_back({"assertion_failures", 0.0, 0.0, "<=", "count"});
  return result;
}

SEP::Testing::Result Fail(const std::string& message) {
  SEP::Testing::Result result;
  result.status = SEP::Testing::Status::Fail;
  result.message = message;
  result.metrics.push_back({"assertion_failures", 1.0, 0.0, "<=", "count"});
  return result;
}

bool ReadFile(const std::string& path, std::string* text) {
  std::ifstream input(path.c_str());
  std::ostringstream buffer;
  buffer << input.rdbuf();
  *text = buffer.str();
  return input.good() || input.eof();
}

bool ReplaceOnce(std::string* text, const std::string& from,
                 const std::string& to) {
  const std::size_t at = text->find(from);
  if (at == std::string::npos || text->find(from, at + 1) != std::string::npos)
    return false;
  text->replace(at, from.size(), to);
  return true;
}

// Replace the example's [parallel_diffusion] body (everything after the header,
// which is the last section of the shipped example) with `body`.
bool WithSection(const std::string& example, const std::string& body,
                 std::string* deck) {
  const std::string header = "\n[parallel_diffusion]\n";
  const std::size_t at = example.find(header);
  if (at == std::string::npos) return false;
  *deck = example.substr(0, at + header.size()) + body;
  return true;
}

bool Parse(const std::string& deck, SEP::Initialization::Configuration* c,
           std::string* message) {
  const SEP::Transport::Status status = SEP::Initialization::ParseText(deck, c);
  if (message) *message = status.message;
  return status.ok();
}

bool Relative(double actual, double expected, double tolerance) {
  return std::fabs(actual - expected) <= tolerance * std::fabs(expected);
}

// Independent relativistic speed: v = p c^2 / sqrt((pc)^2 + (mc^2)^2).
double SpeedFromMomentum(double momentum, double mass) {
  const long double c = 299792458.0L;
  const long double pc = static_cast<long double>(momentum) * c;
  const long double rest = static_cast<long double>(mass) * c * c;
  return static_cast<double>(pc * c / std::sqrt(pc * pc + rest * rest));
}

// PDB01 -- schema-4 section location and legality in the srcSEP INI reader.
SEP::Testing::Result RunPDB01() {
  std::string example;
  if (!ReadFile(gExamplePath, &example) || example.empty())
    return Fail("cannot read shipped schema-4 example " + gExamplePath);
  SEP::Initialization::Configuration c;
  std::string message;
  if (!Parse(example, &c, &message))
    return Fail("shipped schema-4 example did not parse: " + message);
  // The body is kept raw (comments included, original line numbers) and is
  // not interpreted here: model ID and fingerprint stay unresolved.
  bool modelLine = false, valueLine = false;
  for (const SEP::Initialization::ParallelDiffusionLine& line :
       c.parallelDiffusionLines) {
    if (line.text == "model = constant_lambda") modelLine = line.line > 150;
    if (line.text.find("lambda_parallel_m = 1.495978707e10   ! [m]") == 0)
      valueLine = true;
  }
  if (c.schemaVersion != 4 || !c.hasParallelDiffusionSection || !modelLine ||
      !valueLine || !c.parallelDiffusionModelId.empty() ||
      SEP::Initialization::Fingerprint(c).empty())
    return Fail("schema-4 section was not recorded verbatim with line numbers");

  // Library-owned grammar: this reader must not lower-case or reject a
  // record that the library would accept (here a mixed-case SI key).
  std::string mixedCase;
  if (!WithSection(example, "model = power_law_lambda\nrigidity0_V = 1e9\n",
                   &mixedCase) || !Parse(mixedCase, &c, &message))
    return Fail("srcSEP reader interpreted library section records: " + message);
  bool caseKept = false;
  for (const auto& line : c.parallelDiffusionLines)
    if (line.text == "rigidity0_V = 1e9") caseKept = true;
  if (!caseKept) return Fail("srcSEP reader altered a case-sensitive SI key");

  std::string schema3 = example;
  if (!ReplaceOnce(&schema3, "schema_version = 4", "schema_version = 3") ||
      Parse(schema3, &c, &message))
    return Fail("schema 3 accepted [parallel_diffusion]");
  // Schemas 1 and 2 must reject the section as well (the check may not live
  // only inside schema-3 validation).  The schema-2 variant of the example is
  // otherwise invalid, so only the section diagnostic is asserted.
  std::string schema2 = example;
  ReplaceOnce(&schema2, "schema_version = 4", "schema_version = 2");
  if (Parse(schema2, &c, &message) ||
      message.find("[parallel_diffusion] requires run.schema_version=4") ==
          std::string::npos)
    return Fail("schema 2 did not reject [parallel_diffusion]: " + message);
  std::string duplicate = example + "\n[parallel_diffusion]\nmodel = bohm\n";
  if (Parse(duplicate, &c, &message))
    return Fail("duplicate [parallel_diffusion] section was accepted");

  // Without the section, schema 4 is exactly schema 3; the startup
  // fingerprint differs only through the schema number.
  std::string withoutSection = example.substr(0, example.find("\n[parallel_diffusion]\n") + 1);
  SEP::Initialization::Configuration plain;
  if (!Parse(withoutSection, &plain, &message) ||
      plain.hasParallelDiffusionSection)
    return Fail("schema 4 without [parallel_diffusion] was rejected: " + message);
  return Pass("schema-4 [parallel_diffusion] body is recorded verbatim with "
              "line numbers; schema 3 and duplicate sections are rejected");
}

// PDB02 -- selection cross-check and host-capability gate (nothing installed).
SEP::Testing::Result RunPDB02() {
  std::string example;
  if (!ReadFile(gExamplePath, &example)) return Fail("cannot read example");
  SEP::Initialization::Configuration c;
  std::string message;
  if (!Parse(example, &c, &message)) return Fail("example did not parse");

  struct Case {
    const char* name;
    bool librarySelected;
    const char* mover;
    std::string body;   // empty: use the example section unchanged
    bool withInput;
    bool withSection;
  };
  const std::vector<Case> cases = {
      {"section without CLI selection", false, "parker",
       "", true, true},
      {"CLI selection without --input", true,
       "parker", "", false, false},
      {"CLI selection without section", true,
       "parker", "", true, false},
      {"focused-transport mover", true,
       "fte-dmumu", "", true, true},
      {"turbulence-decomposition model", true, "parker",
       "model = qlt_slab_inertial\n", true, true},
      {"energy-per-nucleon law", true,
       "parker",
       "model = power_law_lambda\nlambda0_m = 1e10\n"
       "independent_variable = energy_per_nucleon\n"
       "energy_per_nucleon0_J = 1e-12\nindependent_exponent = 0.3\n",
       true, true},
      {"effective-field Bohm", true, "parker",
       "model = bohm\neta_B = 1\nfield_definition = effective_field\n",
       true, true},
      {"malformed section", true, "parker",
       "model = constant_lambda\nlambda_parallel_m = -1\n", true, true}};
  for (const Case& item : cases) {
    SEP::Initialization::Configuration candidate;
    if (item.withInput) {
      std::string deck = example;
      if (!item.withSection)
        deck = example.substr(0, example.find("\n[parallel_diffusion]\n") + 1);
      else if (!item.body.empty() && !WithSection(example, item.body, &deck))
        return Fail("could not build case deck");
      if (!Parse(deck, &candidate, &message))
        return Fail(std::string(item.name) + ": deck did not parse: " + message);
    }
    const SEP::Transport::Status status = Binding::Configure(
        item.librarySelected, item.mover, item.withInput ? &candidate : NULL);
    if (status.ok())
      return Fail(std::string("accepted ") + item.name);
    if (!candidate.parallelDiffusionModelId.empty())
      return Fail(std::string(item.name) + " wrote a model identity on failure");
  }
  // The malformed-section diagnostic must name the deck's own line.
  std::string bad;
  WithSection(example, "model = constant_lambda\nlambda_parallel_m = -1\n", &bad);
  SEP::Initialization::Configuration badConfig;
  Parse(bad, &badConfig, &message);
  const std::size_t valueLine = badConfig.parallelDiffusionLines.back().line;
  const SEP::Transport::Status badStatus = Binding::Configure(
      true, "parker", &badConfig);
  if (badStatus.message.find("line " + std::to_string(valueLine) + ":") ==
      std::string::npos)
    return Fail("library diagnostic lost the --input line number: " +
                badStatus.message);

  // Other providers are unaffected when no section is present.
  SEP::Initialization::Configuration plain;
  Parse(example.substr(0, example.find("\n[parallel_diffusion]\n") + 1),
        &plain, &message);
  if (!Binding::Configure(false, "parker", &plain).ok() ||
      !Binding::Configure(false, "fte-mfp", NULL).ok())
    return Fail("non-library provider was rejected by the binding");

  if (Binding::IsInstalled() || PD::HasActiveConfiguration())
    return Fail("a rejected configuration installed a library model");
  return Pass("section/CLI/mover cross-check and host-capability gate reject "
              "every inconsistent or unsupported selection without installing");
}

// PDB03 -- registry exposure of the new spatial provider.
SEP::Testing::Result RunPDB03() {
  SpatialKind kind = SpatialKind::FromPitchAngle;
  if (!Coefficient::ParseSpatial("Parallel-Diffusion-Library", &kind) ||
      kind != SpatialKind::ParallelDiffusionLibrary ||
      std::string(Coefficient::SpatialName(kind)) !=
          "parallel-diffusion-library")
    return Fail("parallel-diffusion-library is not a registry spatial provider");
  bool listed = false;
  for (const auto& descriptor : Coefficient::SpatialRegistry())
    if (descriptor.canonicalName == "parallel-diffusion-library") listed = true;
  Coefficient::Configuration configuration;
  configuration.spatial = SpatialKind::ParallelDiffusionLibrary;
  const bool parker =
      Coefficient::ValidateMoverCompatibility(configuration, "parker").ok();
  const bool dmumu =
      Coefficient::ValidateMoverCompatibility(configuration, "fte-dmumu").ok();
  configuration.meanFreePath = Coefficient::MeanFreePathKind::FromSpatial;
  const bool mfp =
      Coefficient::ValidateMoverCompatibility(configuration, "fte-mfp").ok();
  Coefficient::Configuration legacy;
  legacy.spatial = SpatialKind::FromMeanFreePath;
  const std::string legacyFingerprint =
      Coefficient::ConfigurationFingerprint(legacy);
  if (!listed || !parker || dmumu || mfp ||
      legacyFingerprint == Coefficient::ConfigurationFingerprint(configuration))
    return Fail("registry listing, Parker-only compatibility, or fingerprint "
                "distinction failed");
  return Pass("registry parses/lists the provider and admits it only for parker");
}

// PDB04 -- installation and mover-facing evaluation (installs a model).
SEP::Testing::Result RunPDB04() {
  std::string example;
  if (!ReadFile(gExamplePath, &example)) return Fail("cannot read example");
  SEP::Initialization::Configuration c;
  std::string message;
  if (!Parse(example, &c, &message)) return Fail("example did not parse");
  const std::string unresolved = SEP::Initialization::Fingerprint(c);
  const SEP::Transport::Status installed = Binding::Configure(
      true, "parker", &c);
  if (!installed.ok())
    return Fail("valid example section was not installed: " + installed.message);
  if (!Binding::IsInstalled() || c.parallelDiffusionModelId != "constant_lambda" ||
      c.parallelDiffusionConfigurationFingerprint.size() != 64 ||
      PD::ActiveParallelDiffusion !=
          PD::BoundFunctionForModel(PD::ModelId::ConstantLambda) ||
      SEP::Initialization::Fingerprint(c) == unresolved)
    return Fail("installation did not publish model, pointer, and fingerprint");

  // 10 MeV proton: pc = sqrt(T (T + 2 m c^2)), evaluated independently.
  const double mass = 1.67262192369e-27, charge = 1.602176634e-19;
  const double c0 = 299792458.0, kinetic = 10.0e6 * 1.602176634e-19;
  const double momentum =
      std::sqrt(kinetic * (kinetic + 2.0 * mass * c0 * c0)) / c0;
  const double position[3] = {1.0e11, 2.0e10, 0.0};
  const double field[3] = {5.0e-9, -1.0e-9, 0.0};
  const double lambda = 1.495978707e10;
  const Binding::KappaSample sample = Binding::EvaluateKappa(
      mass, charge, momentum, position, field, 3600.0, 7);
  const double expectedKappa = SpeedFromMomentum(momentum, mass) * lambda / 3.0;
  if (!sample.status.ok() || sample.modelId != "constant_lambda" ||
      !Relative(sample.lambdaParallelM, lambda, 1.0e-15) ||
      !Relative(sample.kappaParallelM2PerS, expectedKappa, 1.0e-13) ||
      sample.configurationFingerprint !=
          c.parallelDiffusionConfigurationFingerprint)
    return Fail("constant_lambda kappa differs from v*lambda/3: " +
                sample.status.message);

  // A second, radius-dependent model through the same pointer:
  // lambda = lambda0 (R/R0)^a (r/r0)^b with R = pc/|q|.
  std::string power;
  if (!WithSection(example,
          "model = power_law_lambda\n"
          "lambda0_m = 2.0e10\n"
          "independent_variable = rigidity\n"
          "rigidity0_V = 1.0e9\n"
          "independent_exponent = 0.5    ! a\n"
          "radius0_m = 1.495978707e11\n"
          "radial_exponent = 1.0         # b\n", &power) ||
      !Parse(power, &c, &message))
    return Fail("power-law deck did not parse: " + message);
  if (!Binding::Configure(
      true, "parker", &c).ok() ||
      PD::ActiveParallelDiffusion !=
          PD::BoundFunctionForModel(PD::ModelId::PowerLawLambda))
    return Fail("power-law section did not replace the active model");
  const Binding::KappaSample powerSample = Binding::EvaluateKappa(
      mass, charge, momentum, position, field, 3600.0, 7);
  const double rigidity = momentum * c0 / charge;
  const double radius = std::sqrt(position[0] * position[0] +
                                  position[1] * position[1]);
  const double expectedLambda = 2.0e10 * std::sqrt(rigidity / 1.0e9) *
      (radius / 1.495978707e11);
  if (!powerSample.status.ok() ||
      !Relative(powerSample.lambdaParallelM, expectedLambda, 1.0e-12) ||
      !Relative(powerSample.kappaParallelM2PerS,
                SpeedFromMomentum(momentum, mass) * expectedLambda / 3.0,
                1.0e-12))
    return Fail("power_law_lambda kappa differs from the closed form: " +
                powerSample.status.message);

  // A model failure is reported, never replaced by a value.
  const Binding::KappaSample neutral = Binding::EvaluateKappa(
      mass, 0.0, momentum, position, field, 3600.0, 7);
  if (neutral.status.ok() || neutral.kappaParallelM2PerS != 0.0)
    return Fail("library failure for a neutral particle was not propagated");
  return Pass("installed models are called through ActiveParallelDiffusion and "
              "match independent closed forms; failures propagate");
}

// PDB05 -- schema-4 [run] particle_mover (production mover from the input file).
// The expected mover for each spelling is written here, not obtained from the
// registry parser under test.
SEP::Testing::Result RunPDB05() {
  std::string example;
  if (!ReadFile(gExamplePath, &example)) return Fail("cannot read example");
  SEP::Initialization::Configuration c;
  std::string message;
  if (!Parse(example, &c, &message) || !c.particleMoverSpecified ||
      c.particleMover != SEP::Mover::ProductionMover::Parker ||
      c.particleMoverLine == 0)
    return Fail("shipped example did not select parker through "
                "run.particle_mover: " + message);
  const std::string parkerFingerprint = SEP::Initialization::Fingerprint(c);

  struct Accepted {
    const char* spelling;
    SEP::Mover::ProductionMover mover;
  };
  const Accepted accepted[] = {
      {"fte-dmumu", SEP::Mover::ProductionMover::FocusedTransportDiffusion},
      {"FTE-MFP", SEP::Mover::ProductionMover::FocusedTransportMeanFreePath}};
  for (const Accepted& item : accepted) {
    std::string deck = example;
    if (!ReplaceOnce(&deck, "particle_mover = parker",
                     std::string("particle_mover = ") + item.spelling) ||
        !Parse(deck, &c, &message) || c.particleMover != item.mover)
      return Fail(std::string("canonical mover '") + item.spelling +
                  "' was not accepted: " + message);
    if (SEP::Initialization::Fingerprint(c) == parkerFingerprint)
      return Fail("startup fingerprint does not depend on run.particle_mover");
  }

  // Rejections: deprecated CLI aliases and retired/unknown names, a schema-3
  // file, and a duplicate key.  The bad-value diagnostic names the line.
  const char* rejected[] = {"focused-diffusion", "ParticleMover_Parker_Dxx",
                            "boris", ""};
  for (const char* spelling : rejected) {
    std::string deck = example;
    ReplaceOnce(&deck, "particle_mover = parker",
                std::string("particle_mover = ") + spelling);
    if (Parse(deck, &c, &message))
      return Fail(std::string("accepted non-canonical mover '") + spelling + "'");
  }
  std::string schema3 = example.substr(0, example.find("\n[parallel_diffusion]\n") + 1);
  ReplaceOnce(&schema3, "schema_version = 4", "schema_version = 3");
  if (Parse(schema3, &c, &message) ||
      message.find("run.particle_mover requires run.schema_version=4") ==
          std::string::npos)
    return Fail("schema 3 accepted run.particle_mover: " + message);
  std::string duplicate = example;
  ReplaceOnce(&duplicate, "particle_mover = parker",
              "particle_mover = parker\nparticle_mover = fte-mfp");
  if (Parse(duplicate, &c, &message))
    return Fail("duplicate run.particle_mover was accepted");

  // Absent key: schema 4 stays valid and leaves the choice to CLI/default.
  std::string absent = example;
  ReplaceOnce(&absent, "particle_mover = parker\n", "");
  if (!Parse(absent, &c, &message) || c.particleMoverSpecified)
    return Fail("schema 4 without run.particle_mover was rejected: " + message);
  return Pass("run.particle_mover accepts canonical movers only, is schema-4 "
              "only, enters the startup fingerprint, and is optional");
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "usage: test_parallel_diffusion_binding <schema-4 example>\n";
    return 2;
  }
  gExamplePath = argv[1];
  auto make = [](const char* id, const char* name, const char* description,
                 SEP::Testing::TestCallback callback) {
    SEP::Testing::Descriptor descriptor;
    descriptor.id = id;
    descriptor.name = name;
    descriptor.group = "parallel-diffusion";
    descriptor.description = description;
    descriptor.initialization = SEP::Testing::InitializationLevel::None;
    descriptor.supportedBuildModes = "standalone-no-AMPS";
    descriptor.runtime = SEP::Testing::RuntimeClass::Routine;
    descriptor.seedPolicy = "deterministic; no RNG";
    descriptor.stateIsolation =
        "process-global library model; PDB01-PDB03 install nothing, PDB04 last";
    descriptor.callback = callback;
    return descriptor;
  };
  SEP::Testing::Registry registry({
      make("PDB01", "Schema-4 section", "srcSEP INI locates [parallel_diffusion] verbatim.", RunPDB01),
      make("PDB02", "Selection cross-check", "CLI/section/mover/host gate rejections.", RunPDB02),
      make("PDB03", "Registry provider", "parallel-diffusion-library spatial provider.", RunPDB03),
      make("PDB04", "Installed evaluation", "ActiveParallelDiffusion kappa vs closed forms.", RunPDB04),
      make("PDB05", "Schema-4 mover key", "run.particle_mover parsing, schema gate, fingerprint.", RunPDB05)});
  const SEP::Testing::Summary summary = registry.Run(
      registry.Select(std::vector<std::string>(), std::vector<std::string>(),
                      true), std::cout);
  return summary.ExitCode();
}

#include "tests.h"
#include "util/sep_acceptance_cases.h"
#include "util/sep_mover_validation.h"
#include "util/sep_scientific_validation.h"
#include "util/sep_swcme_validation.h"
#include "util/sep_initialization_validation.h"
#include "util/sep_turbulence_validation.h"
#include "validation/cases/CV01/cv01_model.h"
#include "validation/cases/controlled_transport_models.h"
#include "validation/cases/advanced_validation_models.h"
#include "validation/cases/integrated_validation_models.h"
#include "validation/cases/cross_model_validation_models.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace {

SEP::Testing::Result BooleanResult(bool passed, const char* passMessage,
                                   const char* failMessage) {
  SEP::Testing::Result result;
  result.status = passed ? SEP::Testing::Status::Pass
                         : SEP::Testing::Status::Fail;
  result.message = passed ? passMessage : failMessage;
  result.metrics.push_back(
      {"assertion_failures", passed ? 0.0 : 1.0, 0.0, "<=", "count"});
  return result;
}

SEP::Testing::Result RunDxx() {
  SEP::Testing::Result result = BooleanResult(DxxTest(),
      "analytical and quadrature spatial-diffusion checks passed",
      "spatial-diffusion analytical or quadrature check exceeded tolerance");
  result.configuration.push_back("speed_rounds_m_per_s=1e5,1e6,1e7,1e8");
  result.configuration.push_back("quadrature_panels=1000000");
  result.configuration.push_back("relative_tolerance=1e-5");
  return result;
}

SEP::Testing::Result RunFteConvection() {
  // Registry adapters use fixed seeds because their process exits after the
  // selected test set.  Re-seeding before each stochastic test makes results
  // independent of earlier selector overlap/order without perturbing any
  // subsequent production simulation (there is none in test-only mode).
  const std::uint64_t seed = 1002;
  rnd_seed(static_cast<int>(seed));
  SEP::Testing::Result result = BooleanResult(FTE_Convectoin(),
      "focused-transport convection preserved the expected displacement and velocity",
      "focused-transport convection exceeded its existing numerical tolerance");
  result.hasSeed = true;
  result.seed = seed;
  result.configuration.push_back("trials=10000");
  result.configuration.push_back("plasma_velocity_m_per_s=0,0,0");
  result.configuration.push_back("dt_s=1");
  return result;
}

SEP::Testing::Result RunParkerConvection() {
  const std::uint64_t seed = 1001;
  rnd_seed(static_cast<int>(seed));
  SEP::Testing::Result result = BooleanResult(ParkerModelMoverTest_convection(),
      "Parker convection and adiabatic momentum assertions passed",
      "Parker convection or adiabatic momentum assertion failed");
  result.hasSeed = true;
  result.seed = seed;
  result.configuration.push_back("trials=10000");
  result.configuration.push_back("density_jump=1_to_4");
  result.configuration.push_back("dt_s=1");
  return result;
}

SEP::Testing::Result RunParkerDistribution() {
  const std::uint64_t seed = 1003;
  rnd_seed(static_cast<int>(seed));
  SEP::Testing::Result result = BooleanResult(
      ParkerModelMoverTest_const_plasma_field(),
      "stochastic Parker campaign produced and wrote an in-range histogram",
      "stochastic Parker campaign produced no in-range histogram or output failed");
  result.hasSeed = true;
  result.seed = seed;
  result.configuration.push_back("trials=4000000");
  result.configuration.push_back("dt_s=3");
  result.configuration.push_back("histogram_bins=100");
  result.artifacts.push_back("dxParker.dat");
  return result;
}

SEP::Testing::Result RunScatteringBeyondOneAu() {
  // The historical routine performs 400,000 stochastic particle histories and
  // writes multiple diagnostics, but it defines no numerical acceptance rule
  // and includes an E=0 case whose mean-free-path expression is singular.  A
  // truthful registry must not infer PASS merely because such a routine
  // returned.  It remains discoverable/individually selectable and returns SKIP
  // until a later physics-validation step supplies an approved reference and
  // finite-energy domain; --run-test-manager retains the historical execution.
  SEP::Testing::Result result;
  result.status = SEP::Testing::Status::Skip;
  result.message =
      "legacy diagnostic has no authoritative acceptance criterion; use "
      "--run-test-manager for historical output execution";
  result.hasSeed = false;
  result.configuration.push_back("legacy_energy_grid_MeV=0,50,100,150");
  result.configuration.push_back("registry_execution=disabled_pending_reference");
  result.artifacts.push_back("rmax-E=<energy>MeV.dat");
  result.artifacts.push_back("time-E=<energy>MeV.dat");
  return result;
}

SEP::Testing::Result RunTurbulenceEnergyClosure() {
  // This assertion uses the public production helper and an independently
  // evaluated magnetic-energy expression.  It is intentionally initialization
  // free, deterministic, and non-mutating, so the turbulence group has one real
  // bounded component assertion without pretending to replace the later
  // transport/cascade/reflection validation campaign.
  const double magneticFieldTesla = 5.0e-9;
  const double fractionalFluctuation = 0.2;
  const double expected =
      std::pow(magneticFieldTesla * fractionalFluctuation, 2) /
      (2.0 * SEP::AlfvenTurbulence_Kolmogorov::WaveEnergyConstants::MU0);
  const double actual =
      SEP::AlfvenTurbulence_Kolmogorov::CalculateTypicalWaveEnergyDensity1AU(
          magneticFieldTesla, fractionalFluctuation);
  const double relativeError = std::fabs(actual - expected) / expected;
  const double tolerance = 32.0 * std::numeric_limits<double>::epsilon();

  SEP::Testing::Result result = BooleanResult(
      std::isfinite(actual) && relativeError <= tolerance,
      "wave-energy density matches the independent magnetic-pressure closure",
      "wave-energy density violates the magnetic-pressure closure");
  result.metrics.push_back(
      {"relative_error", relativeError, tolerance, "<=", "dimensionless"});
  result.configuration.push_back("B0_T=5e-9");
  result.configuration.push_back("deltaB_over_B0=0.2");
  return result;
}

SEP::Testing::Result RunCV01LinkedModel() {
  SEP::Testing::Result result;
  const SEP::Testing::ExecutionContext& context =
      SEP::Testing::GetExecutionContext();
  if (context.inputPath.empty() || context.artifactDirectory.empty()) {
    // Generic --all-tests/--test-group invocations have no case-specific
    // configuration by design. Report an explicit prerequisite SKIP so CV01
    // remains discoverable without breaking unrelated registry regression;
    // the validation runner always supplies both paths and therefore cannot
    // turn a missing linked-model run into PASS.
    result.status = SEP::Testing::Status::Skip;
    result.message =
        "CV01 requires --test-input and --test-output-dir; run it through "
        "test/run_tests.py --validation-case CV01 --amps /path/to/amps";
    return result;
  }

  // The Python case layer writes one UTF-8 argument per line. The first line
  // freezes the protocol version; all remaining lines are passed to the linked
  // model in --name/value pairs. This intentionally narrow format avoids a
  // second JSON implementation in C++ while retaining all reviewed JSON and
  // resolved-input evidence in the orchestration layer.
  std::ifstream input(context.inputPath.c_str());
  std::string line;
  if (!input.good() || !std::getline(input, line) ||
      line != "srcsep-cv01-native-args-v1") {
    result.status = SEP::Testing::Status::Error;
    result.message = "CV01 cannot read a supported native argument manifest";
    result.metrics.push_back({"native_model_execution_errors", 1.0, 0.0,
                              "<=", "count"});
    return result;
  }
  std::vector<std::string> arguments;
  while (std::getline(input, line)) {
    if (!line.empty() && line[line.size() - 1] == '\r')
      line.erase(line.size() - 1);
    if (line.empty()) {
      result.status = SEP::Testing::Status::Error;
      result.message = "CV01 native argument manifest contains an empty token";
      result.metrics.push_back({"native_model_execution_errors", 1.0, 0.0,
                                "<=", "count"});
      return result;
    }
    arguments.push_back(line);
  }
  if (!input.eof() || arguments.empty() || arguments.size() % 2 != 0) {
    result.status = SEP::Testing::Status::Error;
    result.message = "CV01 native argument manifest is truncated or unpaired";
    result.metrics.push_back({"native_model_execution_errors", 1.0, 0.0,
                              "<=", "count"});
    return result;
  }

  // Preserve the campaign seed reported by the reviewed input instead of
  // hard-coding the default. Although ballistic D_mumu=0 consumes no random
  // increment, recording the configured stream identity keeps provenance
  // stable when users run a reviewed input variant.
  std::uint64_t campaignSeed = 0;
  bool foundCampaignSeed = false;
  for (std::size_t i = 0; i + 1 < arguments.size(); i += 2) {
    if (arguments[i] != "--campaign-seed") continue;
    char* end = NULL;
    errno = 0;
    const unsigned long long parsed =
        std::strtoull(arguments[i + 1].c_str(), &end, 10);
    if (arguments[i + 1].empty() || arguments[i + 1][0] == '-' ||
        errno == ERANGE || !end || *end != '\0' ||
        parsed > std::numeric_limits<std::uint64_t>::max()) {
      result.status = SEP::Testing::Status::Error;
      result.message = "CV01 native manifest contains an invalid campaign seed";
      result.metrics.push_back({"native_model_execution_errors", 1.0, 0.0,
                                "<=", "count"});
      return result;
    }
    campaignSeed = static_cast<std::uint64_t>(parsed);
    foundCampaignSeed = true;
    break;
  }
  if (!foundCampaignSeed) {
    result.status = SEP::Testing::Status::Error;
    result.message = "CV01 native manifest does not declare a campaign seed";
    result.metrics.push_back({"native_model_execution_errors", 1.0, 0.0,
                              "<=", "count"});
    return result;
  }

  const std::string outputPath =
      context.artifactDirectory + "/CV01_model.csv";
  std::string error;
  const bool completed = SEP::Validation::CV01::RunModel(
      arguments, outputPath, &error);
  result.status = completed ? SEP::Testing::Status::Pass
                            : SEP::Testing::Status::Error;
  result.message = completed
      ? "linked srcSEP/AMPS application completed the controlled CV01 model stage"
      : "linked CV01 model stage failed: " + error;
  result.hasSeed = true;
  result.seed = campaignSeed;
  result.configuration.push_back("execution=linked-srcsep-amps");
  result.configuration.push_back("mover=fte-dmumu");
  result.configuration.push_back("background=analytic-uniform");
  result.configuration.push_back("Dmumu_s^-1=0");
  result.configuration.push_back("native_input=" + context.inputPath);
  result.configuration.push_back("native_artifact_directory=" +
                                 context.artifactDirectory);
  result.metrics.push_back({"native_model_execution_errors",
                            completed ? 0.0 : 1.0, 0.0, "<=", "count"});
  if (completed) result.artifacts.push_back(outputPath);
  return result;
}

SEP::Testing::Result RunLinkedControlledModel(const char* caseId) {
  SEP::Testing::Result result;
  const SEP::Testing::ExecutionContext& context =
      SEP::Testing::GetExecutionContext();
  if (context.inputPath.empty() || context.artifactDirectory.empty()) {
    // CV02-CV12, IV01-IV06, XM01-XM03, OV01-OV05, and EV01-EV02
    // require reviewed case-specific SI inputs. A generic registry
    // run cannot safely invent them, so discovery/all-tests gets an explicit
    // prerequisite SKIP. The external runner always supplies both paths and
    // requires the linked model CSV, so this cannot become a scientific PASS.
    result.status = SEP::Testing::Status::Skip;
    result.message = std::string(caseId) +
        " requires --test-input and --test-output-dir; run it through "
        "test/run_tests.py --validation-case " + caseId +
        " --amps /path/to/amps";
    return result;
  }

  // The narrow line protocol avoids a second permissive JSON parser in the
  // linked application. Its case-qualified header prevents a CV02 manifest,
  // for example, from being executed accidentally as CV03.
  std::ifstream input(context.inputPath.c_str());
  std::string line;
  const std::string expectedHeader =
      std::string("srcsep-controlled-native-args-v1:") + caseId;
  if (!input.good() || !std::getline(input, line) || line != expectedHeader) {
    result.status = SEP::Testing::Status::Error;
    result.message = std::string(caseId) +
        " cannot read its controlled native argument manifest";
    result.metrics.push_back({"native_model_execution_errors", 1.0, 0.0,
                              "<=", "count"});
    return result;
  }
  std::vector<std::string> arguments;
  while (std::getline(input, line)) {
    if (!line.empty() && line[line.size() - 1] == '\r')
      line.erase(line.size() - 1);
    if (line.empty()) {
      result.status = SEP::Testing::Status::Error;
      result.message = std::string(caseId) +
          " native argument manifest contains an empty token";
      result.metrics.push_back({"native_model_execution_errors", 1.0, 0.0,
                                "<=", "count"});
      return result;
    }
    arguments.push_back(line);
  }
  if (!input.eof() || arguments.empty() || arguments.size() % 2 != 0) {
    result.status = SEP::Testing::Status::Error;
    result.message = std::string(caseId) +
        " native argument manifest is truncated or unpaired";
    result.metrics.push_back({"native_model_execution_errors", 1.0, 0.0,
                              "<=", "count"});
    return result;
  }

  // A seed is mandatory even for deterministic characteristics. Preserving
  // the reviewed stream identity makes every case reproducible and prevents a
  // later stochastic extension from silently changing its evidence contract.
  std::uint64_t campaignSeed = 0;
  bool foundSeed = false;
  for (std::size_t i = 0; i + 1 < arguments.size(); i += 2) {
    if (arguments[i] != "--campaign-seed") continue;
    char* end = NULL;
    errno = 0;
    const unsigned long long parsed =
        std::strtoull(arguments[i + 1].c_str(), &end, 10);
    if (arguments[i + 1].empty() || arguments[i + 1][0] == '-' ||
        errno == ERANGE || !end || *end != '\0' ||
        parsed > std::numeric_limits<std::uint64_t>::max()) {
      result.status = SEP::Testing::Status::Error;
      result.message = std::string(caseId) +
          " native manifest contains an invalid campaign seed";
      result.metrics.push_back({"native_model_execution_errors", 1.0, 0.0,
                                "<=", "count"});
      return result;
    }
    campaignSeed = static_cast<std::uint64_t>(parsed);
    foundSeed = true;
    break;
  }
  if (!foundSeed) {
    result.status = SEP::Testing::Status::Error;
    result.message = std::string(caseId) +
        " native manifest does not declare a campaign seed";
    result.metrics.push_back({"native_model_execution_errors", 1.0, 0.0,
                              "<=", "count"});
    return result;
  }

  const std::string outputPath =
      context.artifactDirectory + "/" + caseId + "_model.csv";
  std::string error;
  // CV02-CV05 use the original controlled-transport collection; CV06-CV12
  // use the advanced collection. IV and XM/OV/EV cases have dedicated
  // collections.
  // Every collection enters through this native callback so
  // the Python campaign validates the linked srcSEP/AMPS application rather
  // than compiling and executing a replacement model.
  const std::string identifier(caseId);
  bool completed = false;
  if (identifier.size() >= 2 &&
      (identifier.substr(0, 2) == "XM" ||
       identifier.substr(0, 2) == "OV" ||
       identifier.substr(0, 2) == "EV"))
    completed = SEP::Validation::RunCrossModelValidationModel(
        identifier, arguments, outputPath, &error);
  else if (identifier.size() >= 2 && identifier.substr(0, 2) == "IV")
    completed = SEP::Validation::RunIntegratedValidationModel(
        identifier, arguments, outputPath, &error);
  else if (identifier >= "CV06")
    completed = SEP::Validation::RunAdvancedValidationModel(
        identifier, arguments, outputPath, &error);
  else
    completed = SEP::Validation::RunControlledTransportModel(
        identifier, arguments, outputPath, &error);
  result.status = completed ? SEP::Testing::Status::Pass
                            : SEP::Testing::Status::Error;
  result.message = completed
      ? std::string("linked srcSEP/AMPS application completed ") + caseId +
            " controlled model stage"
      : std::string("linked ") + caseId + " model stage failed: " + error;
  result.hasSeed = true;
  result.seed = campaignSeed;
  result.configuration.push_back("execution=linked-srcsep-amps");
  result.configuration.push_back(std::string("validation_case=") + caseId);
  result.configuration.push_back("units=SI");
  result.configuration.push_back("native_input=" + context.inputPath);
  result.configuration.push_back("native_artifact_directory=" +
                                 context.artifactDirectory);
  result.metrics.push_back({"native_model_execution_errors",
                            completed ? 0.0 : 1.0, 0.0, "<=", "count"});
  if (completed) result.artifacts.push_back(outputPath);
  return result;
}

SEP::Testing::Result RunCV02LinkedModel() {
  return RunLinkedControlledModel("CV02");
}

SEP::Testing::Result RunCV03LinkedModel() {
  return RunLinkedControlledModel("CV03");
}

SEP::Testing::Result RunCV04LinkedModel() {
  return RunLinkedControlledModel("CV04");
}

SEP::Testing::Result RunCV05LinkedModel() {
  return RunLinkedControlledModel("CV05");
}

SEP::Testing::Result RunCV06LinkedModel() { return RunLinkedControlledModel("CV06"); }
SEP::Testing::Result RunCV07LinkedModel() { return RunLinkedControlledModel("CV07"); }
SEP::Testing::Result RunCV08LinkedModel() { return RunLinkedControlledModel("CV08"); }
SEP::Testing::Result RunCV09LinkedModel() { return RunLinkedControlledModel("CV09"); }
SEP::Testing::Result RunCV10LinkedModel() { return RunLinkedControlledModel("CV10"); }
SEP::Testing::Result RunCV11LinkedModel() { return RunLinkedControlledModel("CV11"); }
SEP::Testing::Result RunCV12LinkedModel() { return RunLinkedControlledModel("CV12"); }
SEP::Testing::Result RunIV01LinkedModel() { return RunLinkedControlledModel("IV01"); }
SEP::Testing::Result RunIV02LinkedModel() { return RunLinkedControlledModel("IV02"); }
SEP::Testing::Result RunIV03LinkedModel() { return RunLinkedControlledModel("IV03"); }
SEP::Testing::Result RunIV04LinkedModel() { return RunLinkedControlledModel("IV04"); }
SEP::Testing::Result RunIV05LinkedModel() { return RunLinkedControlledModel("IV05"); }
SEP::Testing::Result RunIV06LinkedModel() { return RunLinkedControlledModel("IV06"); }
SEP::Testing::Result RunXM01LinkedModel() { return RunLinkedControlledModel("XM01"); }
SEP::Testing::Result RunXM02LinkedModel() { return RunLinkedControlledModel("XM02"); }
SEP::Testing::Result RunXM03LinkedModel() { return RunLinkedControlledModel("XM03"); }
SEP::Testing::Result RunOV01LinkedModel() { return RunLinkedControlledModel("OV01"); }
SEP::Testing::Result RunOV02LinkedModel() { return RunLinkedControlledModel("OV02"); }
SEP::Testing::Result RunOV03LinkedModel() { return RunLinkedControlledModel("OV03"); }
SEP::Testing::Result RunOV04LinkedModel() { return RunLinkedControlledModel("OV04"); }
SEP::Testing::Result RunOV05LinkedModel() { return RunLinkedControlledModel("OV05"); }
SEP::Testing::Result RunEV01LinkedModel() { return RunLinkedControlledModel("EV01"); }
SEP::Testing::Result RunEV02LinkedModel() { return RunLinkedControlledModel("EV02"); }

SEP::Testing::Descriptor MakeDescriptor(
    const char* id, const char* name, const char* group,
    const char* description, SEP::Testing::InitializationLevel initialization,
    SEP::Testing::RuntimeClass runtime, const char* seedPolicy,
    const char* stateIsolation, SEP::Testing::TestCallback callback) {
  SEP::Testing::Descriptor descriptor;
  descriptor.id = id;
  descriptor.name = name;
  descriptor.group = group;
  descriptor.description = description;
  descriptor.initialization = initialization;
  descriptor.supportedBuildModes = "serial and MPI standalone executable";
  descriptor.runtime = runtime;
  descriptor.seedPolicy = seedPolicy;
  descriptor.stateIsolation = stateIsolation;
  descriptor.callback = callback;
  return descriptor;
}

SEP::Testing::Descriptor MakeCV01Descriptor() {
  SEP::Testing::Descriptor descriptor = MakeDescriptor(
      "CV01", "Ballistic streaming on a uniform field line",
      "controlled-analytical",
      "Run the numerical stage inside the linked srcSEP/AMPS application; the Python case runner supplies the independent characteristic and final acceptance.",
      SEP::Testing::InitializationLevel::None,
      SEP::Testing::RuntimeClass::Routine, "fixed case seed from reviewed input",
      "isolated process and artifact directory; immutable synthetic input",
      RunCV01LinkedModel);
  // The callback writes one synthetic, non-decomposed particle CSV before MPI
  // initialization. The campaign therefore refuses --mpi-np and advertises
  // serial support until a future case defines rank-partitioned output.
  descriptor.supportedBuildModes = "serial linked srcSEP/AMPS executable";
  return descriptor;
}

SEP::Testing::Descriptor MakeControlledValidationDescriptor(
    const char* id, const char* name, const char* description,
    SEP::Testing::TestCallback callback) {
  SEP::Testing::Descriptor descriptor = MakeDescriptor(
      id, name, "controlled-analytical", description,
      SEP::Testing::InitializationLevel::None,
      SEP::Testing::RuntimeClass::Routine, "fixed case seed from reviewed input",
      "isolated process and artifact directory; immutable synthetic input",
      callback);
  // Model evidence is one non-decomposed CSV written before normal AMPS model
  // initialization. Serial is the only honest advertised mode until a case
  // defines rank partitioning plus deterministic reduction of its samples.
  descriptor.supportedBuildModes = "serial linked srcSEP/AMPS executable";
  return descriptor;
}

SEP::Testing::Descriptor MakeIntegratedValidationDescriptor(
    const char* id, const char* name, const char* description,
    SEP::Testing::TestCallback callback) {
  // Integrated manufactured cases exercise several production operators in
  // one controlled scenario.  Give them their own registry group so users can
  // select the IV campaign independently from the single-operator CV cases
  // with `--test-group integrated-manufactured`.
  SEP::Testing::Descriptor descriptor = MakeDescriptor(
      id, name, "integrated-manufactured", description,
      SEP::Testing::InitializationLevel::None,
      SEP::Testing::RuntimeClass::Routine, "fixed case seed from reviewed input",
      "isolated process and artifact directory; immutable synthetic input",
      callback);
  // Each IV callback writes a single deterministic evidence table before the
  // regular AMPS model initialization.  MPI execution is intentionally not
  // advertised until rank-local sampling and deterministic reduction are part
  // of the validation contract.
  descriptor.supportedBuildModes = "serial linked srcSEP/AMPS executable";
  return descriptor;
}

SEP::Testing::Descriptor MakeCrossModelValidationDescriptor(
    const char* id, const char* name, const char* description,
    SEP::Testing::TestCallback callback) {
  // Cross-model cases remain separate from analytical/manufactured evidence:
  // agreement with another numerical model is useful validation evidence but
  // does not establish correctness against an exact solution.
  SEP::Testing::Descriptor descriptor = MakeDescriptor(
      id, name, "cross-model", description,
      SEP::Testing::InitializationLevel::None,
      SEP::Testing::RuntimeClass::Extended, "fixed case seed from reviewed input",
      "isolated process; immutable reference and transactional model CSV",
      callback);
  descriptor.supportedBuildModes = "serial linked srcSEP/AMPS executable";
  return descriptor;
}

SEP::Testing::Descriptor MakeObservationalValidationDescriptor(
    const char* id, const char* name, const char* description,
    SEP::Testing::TestCallback callback) {
  // Observational/campaign cases are separate from cross-model tests because
  // their immutable baselines are spacecraft measurements.  The external
  // case descriptor further records whether the result is a release gate, a
  // diagnostic-only comparison, or a pilot-incomplete population study.
  SEP::Testing::Descriptor descriptor = MakeDescriptor(
      id, name, "observational-validation", description,
      SEP::Testing::InitializationLevel::None,
      SEP::Testing::RuntimeClass::Extended,
      "fixed case seed from reviewed observational input",
      "isolated linked process; immutable observations; transactional CSV",
      callback);
  descriptor.supportedBuildModes = "serial linked srcSEP/AMPS executable";
  return descriptor;
}

int StatusSeverity(SEP::Testing::Status status) {
  switch (status) {
    case SEP::Testing::Status::Pass: return 0;
    case SEP::Testing::Status::Skip: return 1;
    case SEP::Testing::Status::Fail: return 2;
    case SEP::Testing::Status::Error: return 3;
  }
  return 3;
}

std::string SerializeRankEvidence(const SEP::Testing::Result& result) {
  // The payload is intentionally plain text because it is embedded verbatim in
  // both JSON configuration evidence and JUnit system-out.  Structured report
  // writers own escaping, so rank-local diagnostics cannot corrupt either
  // output format.
  std::ostringstream evidence;
  evidence << "status=" << SEP::Testing::StatusName(result.status)
           << ";message=" << result.message;
  if (result.hasSeed) evidence << ";seed=" << result.seed;
  for (std::size_t i = 0; i < result.configuration.size(); ++i)
    evidence << ";configuration=" << result.configuration[i];
  for (std::size_t i = 0; i < result.metrics.size(); ++i)
    evidence << ";metric=" << result.metrics[i].name << ':'
             << result.metrics[i].value << ':'
             << result.metrics[i].comparison << ':'
             << result.metrics[i].tolerance << ':'
             << result.metrics[i].units;
  for (std::size_t i = 0; i < result.artifacts.size(); ++i)
    evidence << ";artifact=" << result.artifacts[i];
  return evidence.str();
}

SEP::Testing::Status StatusFromSeverity(int severity) {
  if (severity == 0) return SEP::Testing::Status::Pass;
  if (severity == 1) return SEP::Testing::Status::Skip;
  if (severity == 2) return SEP::Testing::Status::Fail;
  return SEP::Testing::Status::Error;
}

}  // namespace

const SEP::Testing::Registry& ComponentTestRegistry() {
  // Function-local static construction avoids cross-translation-unit
  // initialization ordering.  Registry construction also validates every
  // descriptor and deterministically sorts this deliberately unsorted source
  // list before it becomes observable through --list-tests.
  static const SEP::Testing::Registry registry([]() {
    std::vector<SEP::Testing::Descriptor> descriptors =
        SEP::Testing::AcceptanceCaseDescriptors();
    const std::vector<SEP::Testing::Descriptor> validation_descriptors =
        SEP::Testing::ScientificValidationDescriptors();
    descriptors.insert(descriptors.end(), validation_descriptors.begin(),
                       validation_descriptors.end());
    // The linked CLI and dependency-light runners intentionally share these
    // exact callbacks.  Keeping the controlled kernels in the production
    // catalog prevents a Make-only test from drifting away from --test ID.
    const std::vector<SEP::Testing::Descriptor> mover_descriptors =
        SEP::Testing::ControlledMoverDescriptors();
    descriptors.insert(descriptors.end(), mover_descriptors.begin(),
                       mover_descriptors.end());
    const std::vector<SEP::Testing::Descriptor> turbulence_descriptors =
        SEP::Testing::ControlledTurbulenceDescriptors();
    descriptors.insert(descriptors.end(), turbulence_descriptors.begin(),
                       turbulence_descriptors.end());
    // D01, D02, and the bounded D03 preflight are part of the same native C++
    // catalog as every other public ID.  test/run_tests.py --all discovers
    // this catalog through --list-tests and launches each descriptor in an
    // isolated linked process, so these checks can no longer be omitted merely
    // because their dependency-light Make targets were not selected.
    const std::vector<SEP::Testing::Descriptor> swcme_descriptors =
        SEP::Testing::SwcmeImprovementDescriptors();
    descriptors.insert(descriptors.end(), swcme_descriptors.begin(),
                       swcme_descriptors.end());
    // These AMPS-independent callbacks are linked into the production catalog
    // as routine tests, so `test/run_tests.py --all` discovers and executes
    // the same INIT01/INIT02 gates as the focused source-only build.
    const std::vector<SEP::Testing::Descriptor> initialization_descriptors =
        SEP::Testing::InitializationDescriptors();
    descriptors.insert(descriptors.end(), initialization_descriptors.begin(),
                       initialization_descriptors.end());
    const SEP::Testing::Descriptor legacy_descriptors[] = {
      MakeCV01Descriptor(),
      MakeControlledValidationDescriptor(
          "CV02", "Constant-coefficient spatial diffusion Green function",
          "Run production Parker diffusion in the linked application and compare packet moments and profiles with exact Gaussian bin integrals externally.",
          RunCV02LinkedModel),
      MakeControlledValidationDescriptor(
          "CV03", "Nonuniform diffusion and stochastic-calculus drift",
          "Run sinusoidal spatial diffusion in the linked application and compare equilibrium/transients with an independent conservative finite-volume reference.",
          RunCV03LinkedModel),
      MakeControlledValidationDescriptor(
          "CV04", "Adiabatic momentum change in expanding solar wind",
          "Run production Parker momentum updates in the linked application and compare constant-divergence and spherical-flow characteristics externally.",
          RunCV04LinkedModel),
      MakeControlledValidationDescriptor(
          "CV05", "Magnetic focusing in a prescribed field gradient",
          "Run production focused transport with zero scattering in the linked application and compare pitch-angle characteristics and invariants externally.",
          RunCV05LinkedModel),
      MakeControlledValidationDescriptor(
          "CV06", "Legendre-mode pitch-angle diffusion",
          "Evolve Legendre modes l=1..6 with the production Dmumu mover and compare decay rates, mode leakage, normalization, and boundary handling externally.",
          RunCV06LinkedModel),
      MakeControlledValidationDescriptor(
          "CV07", "Telegraph transport and diffusion limit",
          "Run persistent random flights with production exponential waiting times and compare fronts, support, moments, profiles, and the late diffusion limit externally.",
          RunCV07LinkedModel),
      MakeControlledValidationDescriptor(
          "CV08", "Absorbing-boundary first-passage distribution",
          "Run production Parker drift-diffusion to an absorbing boundary and compare censored arrival distributions with the exact inverse-Gaussian law externally.",
          RunCV08LinkedModel),
      MakeControlledValidationDescriptor(
          "CV09", "Planar diffusive-shock-acceleration spectrum",
          "Run a controlled shock-cycle benchmark for compression ratios 2, 3, and 4 and compare spectral slopes and acceleration times with planar DSA theory externally.",
          RunCV09LinkedModel),
      MakeControlledValidationDescriptor(
          "CV10", "Nonuniform turbulence advection",
          "Advect both spectral wave branches on fixed, variable-area, and remapped grids with the production conservative turbulence core and score invariants and convergence externally.",
          RunCV10LinkedModel),
      MakeControlledValidationDescriptor(
          "CV11", "Time-dependent resonant wave growth and damping",
          "Apply one-hot, time-dependent wave-energy source histories through the production turbulence ledger and compare exponential growth/damping histories externally.",
          RunCV11LinkedModel),
      MakeControlledValidationDescriptor(
          "CV12", "Controlled particle-wave total-energy exchange",
          "Exercise production wave-frame scattering and the turbulence exchange ledger in closed coupled and uncoupled controls, then score total-energy closure externally.",
          RunCV12LinkedModel),
      MakeIntegratedValidationDescriptor("IV01", "Streaming plus Parker-spiral focusing",
          "Couple analytic Parker geometry, production focusing, time of flight, weak scattering, and reversed field-line storage order.", RunIV01LinkedModel),
      MakeIntegratedValidationDescriptor("IV02", "Scattering-focusing equilibrium",
          "Evolve coupled focusing and Dmumu from isotropic and beam states toward the independent zero-flux equilibrium.", RunIV02LinkedModel),
      MakeIntegratedValidationDescriptor("IV03", "Manufactured spatial-pitch-momentum transport",
          "Evaluate a smooth full transport residual and refinement across advection, pitch, momentum, and diffusion terms.", RunIV03LinkedModel),
      MakeIntegratedValidationDescriptor("IV04", "Moving field-line grid conservation and GCL",
          "Remap a uniform state through stretching, compression, and sinusoidal segment motion and verify free-stream preservation.", RunIV04LinkedModel),
      MakeIntegratedValidationDescriptor("IV05", "Moving shock crossing and acceleration",
          "Compare stationary/moving shock frames, interpolated crossings, node coincidence, and the CV09 DSA reference.", RunIV05LinkedModel),
      MakeIntegratedValidationDescriptor("IV06", "Coupled self-generated turbulence feedback",
          "Compare frozen, one-way, and two-way scattering/growth controls with resonant and total-energy evidence.", RunIV06LinkedModel),
      MakeCrossModelValidationDescriptor("XM01", "Independent focused-transport PDE solver comparison",
          "Compare linked production characteristics with an independent conservative finite-volume solver for streaming, scattering, focusing, adiabatic momentum change, and their combination.", RunXM01LinkedModel),
      MakeCrossModelValidationDescriptor("XM02", "Published M-FLAMPA Parker-spiral comparison",
          "Run a controlled linked first-passage reconstruction for the three published mean free paths and compare with provenance-tracked M-FLAMPA sensitivity curves.", RunXM02LinkedModel),
      MakeCrossModelValidationDescriptor("XM03", "2013 April 11 Earth observation comparison",
          "Run event-informed Parker transport in the linked application and compare its Earth spectra with ACE, GOES-13, and SOHO observations digitized from Liu et al. Figure 12.", RunXM03LinkedModel),
      MakeObservationalValidationDescriptor("OV01", "2013 April 11 near-Earth observational benchmark",
          "Apply the release-gating Earth spectrum contract to linked Parker transport and ACE, GOES-13, and SOHO observations from Liu et al. Figure 12.", RunOV01LinkedModel),
      MakeObservationalValidationDescriptor("OV02", "2020 May 29 radial PSP/STEREO-A benchmark",
          "Run one unchanged transport setup at 0.33 and 0.96 AU and compare with PSP/EPI-Hi and STEREO-A/LET profiles from Cheng et al. Figures 3 and 6.", RunOV02LinkedModel),
      MakeObservationalValidationDescriptor("OV03", "2013 May 22 interacting-CME diagnostic",
          "Compare single- and twin-source linked transport hypotheses with GOES and STEREO-A profiles from Ding et al. Figure 1.", RunOV03LinkedModel),
      MakeObservationalValidationDescriptor("OV04", "2014 January 6 connectivity diagnostic",
          "Propagate three documented connection-delay realizations and compare their event spectra with the PAMELA spectrum in Bruno et al. Figure 4.", RunOV04LinkedModel),
      MakeObservationalValidationDescriptor("OV05", "September 2017 compound-event diagnostic",
          "Propagate publication-timed September 4, 6, and 10 injections and compare with STEREO-A profiles in Bruno et al. Figure 2.", RunOV05LinkedModel),
      MakeObservationalValidationDescriptor("EV01", "Calibration ensemble and parameter identifiability",
          "Generate linked event predictions for the sealed CCMC training/validation ensemble; calibration and observation scoring are performed by the campaign runner.", RunEV01LinkedModel),
      MakeObservationalValidationDescriptor("EV02", "Locked held-out validation campaign",
          "Generate linked event predictions without using held-out targets; the campaign runner applies the frozen split and scores untouched observations.", RunEV02LinkedModel),
      MakeDescriptor("TURB01", "Alfven wave-energy closure", "turbulence",
          "Compare the production 1-AU wave-energy helper with an independent magnetic-pressure expression.",
          SEP::Testing::InitializationLevel::None,
          SEP::Testing::RuntimeClass::Routine, "deterministic; no RNG",
          "read-only pure calculation", RunTurbulenceEnergyClosure),
      MakeDescriptor("PARKER02", "Parker stochastic displacement histogram", "parker",
          "Run the legacy constant-plasma Parker campaign and validate histogram/output completion.",
          SEP::Testing::InitializationLevel::FieldLineModel,
          SEP::Testing::RuntimeClass::Extended, "fixed registry seed 1003",
          "adapter reseeds; particle and segment lists are cleaned; process exits after tests",
          RunParkerDistribution),
      MakeDescriptor("DXX01", "Spatial diffusion coefficient", "diffusion",
          "Compare GetDxx with its constant-coefficient solution and independent numerical quadrature.",
          SEP::Testing::InitializationLevel::FieldLineModel,
          SEP::Testing::RuntimeClass::Routine, "deterministic; no RNG",
          "pitch-angle coefficient function pointer is restored before return", RunDxx),
      MakeDescriptor("SCAT01", "Scattering beyond 1 AU", "scattering",
          "Inventory the legacy return-probability diagnostic and its reference-acceptance limitation.",
          SEP::Testing::InitializationLevel::FieldLineModel,
          SEP::Testing::RuntimeClass::Extended, "uses the configured AMPS random stream; seed not exposed",
          "registry returns SKIP without mutation; legacy TestManager still executes it",
          RunScatteringBeyondOneAu),
      MakeDescriptor("FTE01", "Focused-transport convection", "transport",
          "Check field-line streaming displacement and velocity invariance in a static plasma fixture.",
          SEP::Testing::InitializationLevel::FieldLineModel,
          SEP::Testing::RuntimeClass::Routine, "fixed registry seed 1002",
          "adapter reseeds; plasma/coefficient/particle/list state is restored; process exits after tests", RunFteConvection),
      MakeDescriptor("PARKER01", "Parker convection and adiabatic momentum", "parker",
          "Check stationary line coordinates and the density-driven analytical momentum update.",
          SEP::Testing::InitializationLevel::FieldLineModel,
          SEP::Testing::RuntimeClass::Routine, "fixed registry seed 1001",
          "adapter reseeds; plasma/coefficient/particle/list state is restored; process exits after tests", RunParkerConvection),
    };
    descriptors.insert(descriptors.end(),
        legacy_descriptors,
        legacy_descriptors + sizeof(legacy_descriptors) /
            sizeof(legacy_descriptors[0]));
    return descriptors;
  }());
  return registry;
}

SEP::Testing::InitializationLevel RequiredInitializationLevel(
    const std::vector<const SEP::Testing::Descriptor*>& selected) {
  SEP::Testing::InitializationLevel required =
      SEP::Testing::InitializationLevel::None;
  for (const SEP::Testing::Descriptor* descriptor : selected) {
    if (descriptor &&
        descriptor->initialization ==
            SEP::Testing::InitializationLevel::FieldLineModel) {
      required = SEP::Testing::InitializationLevel::FieldLineModel;
    }
  }
  return required;
}

int RunSelectedComponentTests(
    const std::vector<const SEP::Testing::Descriptor*>& selected,
    std::ostream& out,
    const std::string& jsonReportPath,
    const std::string& junitReportPath) {
  SEP::Testing::Summary summary;
  const SEP::Testing::Registry& registry = ComponentTestRegistry();

  for (const SEP::Testing::Descriptor* descriptor : selected) {
    SEP::Testing::Result result = registry.RunOne(*descriptor);

    // Every required rank executes the callback.  Reduce the categorical result
    // by severity and duration by maximum so rank-local failure cannot be hidden
    // by a root PASS and the reported time reflects the slowest participant.
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (initialized) {
      int rank = 0;
      int ranks = 1;
      MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR, &rank);
      MPI_Comm_size(MPI_GLOBAL_COMMUNICATOR, &ranks);

      // Gather the full pre-reduction outcome from every rank.  A severity-only
      // reduction gives a correct exit status but would discard the message,
      // seed, metrics, and artifact paths from a failing non-root rank.  These
      // payloads make JSON/JUnit reports sufficient to diagnose MPI failures.
      const std::string localEvidence = SerializeRankEvidence(result);
      const int localLength = static_cast<int>(localEvidence.size());
      std::vector<int> lengths(rank == 0 ? ranks : 0, 0);
      MPI_Gather(&localLength, 1, MPI_INT,
                 rank == 0 ? lengths.data() : NULL, 1, MPI_INT, 0,
                 MPI_GLOBAL_COMMUNICATOR);
      std::vector<int> displacements(rank == 0 ? ranks : 0, 0);
      int totalLength = 0;
      if (rank == 0) {
        for (int i = 0; i < ranks; ++i) {
          displacements[i] = totalLength;
          totalLength += lengths[i];
        }
      }
      std::vector<char> gathered(rank == 0 ? totalLength : 0);
      MPI_Gatherv(localEvidence.data(), localLength, MPI_CHAR,
                  rank == 0 && totalLength != 0 ? gathered.data() : NULL,
                  rank == 0 ? lengths.data() : NULL,
                  rank == 0 ? displacements.data() : NULL,
                  MPI_CHAR, 0, MPI_GLOBAL_COMMUNICATOR);

      const int localSeverity = StatusSeverity(result.status);
      int globalSeverity = localSeverity;
      double globalDuration = result.elapsedSeconds;
      MPI_Allreduce(&localSeverity, &globalSeverity, 1, MPI_INT, MPI_MAX,
                    MPI_GLOBAL_COMMUNICATOR);
      MPI_Allreduce(&result.elapsedSeconds, &globalDuration, 1, MPI_DOUBLE,
                    MPI_MAX, MPI_GLOBAL_COMMUNICATOR);
      if (globalSeverity != localSeverity) {
        result.message = "one or more MPI ranks reported a more severe outcome";
      }
      result.status = StatusFromSeverity(globalSeverity);
      result.elapsedSeconds = globalDuration;
      if (rank == 0) {
        for (int i = 0; i < ranks; ++i) {
          result.configuration.push_back(
              "mpi_rank_" + std::to_string(i) + "=" +
              std::string(gathered.data() + displacements[i], lengths[i]));
        }
      }
    }

    if (PIC::ThisThread == 0) SEP::Testing::PrintResult(result, out);
    summary.Add(result);
  }

  if (PIC::ThisThread == 0) {
    out << "Component-test summary: PASS=" << summary.passed
        << " FAIL=" << summary.failed << " SKIP=" << summary.skipped
        << " ERROR=" << summary.errors << '\n';
  }

  // Only MPI root writes shared report paths.  A requested report is acceptance
  // evidence, so open/write failure is broadcast as ERROR to every rank and
  // changes the process status rather than being reduced to a console warning.
  int reportError = 0;
  if (PIC::ThisThread == 0) {
    std::string error;
    if (!jsonReportPath.empty() &&
        !SEP::Testing::WriteJsonSummary(summary, jsonReportPath, &error)) {
      out << "ERROR: " << error << '\n';
      reportError = 1;
    }
    error.clear();
    if (!junitReportPath.empty() &&
        !SEP::Testing::WriteJUnitSummary(summary, junitReportPath, &error)) {
      out << "ERROR: " << error << '\n';
      reportError = 1;
    }
  }
  int initialized = 0;
  MPI_Initialized(&initialized);
  if (initialized)
    MPI_Bcast(&reportError, 1, MPI_INT, 0, MPI_GLOBAL_COMMUNICATOR);
  return reportError ? 2 : summary.ExitCode();
}

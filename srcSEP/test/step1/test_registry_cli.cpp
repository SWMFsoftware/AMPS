#include "../../util/sep_cli.h"
#include "sep_test_registry.h"
#include "sep_transport_common.h"

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <cerrno>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <unistd.h>
#include <sys/stat.h>

namespace {

int failures = 0;

void Check(bool condition, const char* id, const std::string& message) {
  if (!condition) {
    std::cerr << '[' << id << "] FAIL: " << message << '\n';
    ++failures;
  }
}

bool Parse(const std::vector<std::string>& arguments,
           SEP::Util::CLI::Options& options, std::string& error) {
  std::vector<std::string> storage(arguments);
  std::vector<char*> argv;
  for (std::string& item : storage) argv.push_back(&item[0]);
  std::ostringstream out;
  std::ostringstream err;
  const bool ok = SEP::Util::CLI::ParseCommandLine(
      static_cast<int>(argv.size()), argv.data(), options, out, err);
  error = err.str();
  return ok;
}

SEP::Testing::Descriptor Descriptor(
    const char* id, const char* group, SEP::Testing::RuntimeClass runtime,
    SEP::Testing::Status status, int* callCounter) {
  SEP::Testing::Descriptor descriptor;
  descriptor.id = id;
  descriptor.name = std::string("fixture ") + id;
  descriptor.group = group;
  descriptor.description = "test-only registry fixture";
  descriptor.initialization = SEP::Testing::InitializationLevel::None;
  descriptor.supportedBuildModes = "test fixture";
  descriptor.runtime = runtime;
  descriptor.seedPolicy = "deterministic; no RNG";
  descriptor.stateIsolation = "test-local counter only";
  descriptor.callback = [status, callCounter]() {
    ++(*callCounter);
    SEP::Testing::Result result;
    result.status = status;
    result.message = "fixture outcome";
    return result;
  };
  return descriptor;
}

void TestCli01(const SEP::Testing::Registry& registry) {
  SEP::Util::CLI::Options help;
  SEP::Util::CLI::Options list;
  std::string error;
  Check(Parse({"sep", "--help"}, help, error) && help.printHelp,
        "CLI01", "--help must select the pre-initialization exit path");
  Check(Parse({"sep", "--list-tests"}, list, error) && list.listTests &&
            !SEP::Util::CLI::IsComponentTestExecutionRequested(list),
        "CLI01", "--list-tests must not become an execution selector");
  SEP::Util::CLI::Options input;
  Check(Parse({"sep", "--input=examples/sep_parker_mesh.in"}, input, error) &&
            input.inputPath == "examples/sep_parker_mesh.in",
        "CLI01", "--input must preserve its production initialization path");
  // -input is the AMPS-wide spelling of the same option (src/INPUT_FILE.md).
  SEP::Util::CLI::Options oneDash;
  Check(Parse({"sep", "-input", "examples/sep_parker_mesh.in"}, oneDash, error) &&
            oneDash.inputPath == "examples/sep_parker_mesh.in",
        "CLI01", "-input must be an exact synonym of --input");
  SEP::Util::CLI::Options reduced;
  Check(Parse({"sep", "--input", "examples/sep_parker_mesh.in",
               "--reduced-shock-event", "examples/reduced-shock/event.in"},
              reduced, error) &&
            reduced.reducedShockEventPath ==
                "examples/reduced-shock/event.in",
        "CLI01", "reduced event must preserve its production asset path");
  SEP::Util::CLI::Options initializationOnly;
  Check(Parse({"sep", "--input", "examples/sep_parker_mesh.in",
               "--initialization-only", "--initialization-output-dir",
               "preview"}, initializationOnly, error) &&
            initializationOnly.initializationOnly &&
            initializationOnly.initializationOutputDirectory == "preview",
        "CLI01", "initialization-only options must normalize exactly");

  // Default runtime input: without -input a production run reads ./amps.in
  // when it exists; without the file the legacy no-argument path (empty
  // inputPath) is kept, and help/list/test modes never pick the file up.
  // The checks run in a private scratch directory so no repository file is
  // touched; the original working directory is restored afterwards.
  {
    char originalDirectory[4096];
    char scratchTemplate[] = "/tmp/srcsep-cli01-XXXXXX";
    const bool haveOriginal =
        getcwd(originalDirectory, sizeof(originalDirectory)) != nullptr;
    const char* scratch = mkdtemp(scratchTemplate);
    if (haveOriginal && scratch != nullptr && chdir(scratch) == 0) {
      SEP::Util::CLI::Options absent;
      Check(Parse({"sep"}, absent, error) && absent.inputPath.empty() &&
                !absent.inputPathFromDefault,
            "CLI01", "without ./amps.in the legacy no-argument path must be kept");
      { std::ofstream create("amps.in"); create << "! fixture\n"; }
      SEP::Util::CLI::Options present;
      Check(Parse({"sep"}, present, error) && present.inputPath == "amps.in" &&
                present.inputPathFromDefault,
            "CLI01", "without -input an existing ./amps.in must become the input");
      SEP::Util::CLI::Options explicitInput;
      Check(Parse({"sep", "-input", "other.in"}, explicitInput, error) &&
                explicitInput.inputPath == "other.in" &&
                !explicitInput.inputPathFromDefault,
            "CLI01", "an explicit -input must take precedence over ./amps.in");
      SEP::Util::CLI::Options testMode;
      Check(Parse({"sep", "--test", "A01"}, testMode, error) &&
                testMode.inputPath.empty(),
            "CLI01", "component-test modes must not inherit ./amps.in");
      SEP::Util::CLI::Options listMode;
      Check(Parse({"sep", "--list-tests"}, listMode, error) &&
                listMode.inputPath.empty(),
            "CLI01", "list modes must not inherit ./amps.in");
      std::remove("amps.in");
      Check(chdir(originalDirectory) == 0, "CLI01",
            "could not restore the working directory");
      rmdir(scratch);
    } else {
      Check(false, "CLI01", "could not create the ./amps.in scratch directory");
    }
  }

  std::ostringstream listing;
  registry.PrintList(listing);
  const std::string listingText = listing.str();
  Check(listingText.find("A01") < listingText.find("B01"),
        "CLI01", "registry listing must be sorted by stable ID");
  Check(listingText.find("initialization") != std::string::npos &&
            listingText.find("test fixture") != std::string::npos,
        "CLI01", "listing must include prerequisites and descriptions");
}

void TestCli02And03(const SEP::Testing::Registry& registry,
                    int& passCalls, int& skipCalls) {
  SEP::Util::CLI::Options options;
  std::string error;
  Check(Parse({"sep", "--test", "a01", "--test=A01", "--test-group", "beta"},
              options, error),
        "CLI02", "both --test forms and repeated selectors must parse");
  const std::vector<const SEP::Testing::Descriptor*> selected =
      registry.Select(options.testIds, options.testGroups, options.runAllTests);
  Check(selected.size() == 2 && selected[0]->id == "A01" &&
            selected[1]->id == "B01",
        "CLI03", "overlapping selectors must de-duplicate and sort by ID");

  std::ostringstream output;
  const SEP::Testing::Summary summary = registry.Run(selected, output);
  Check(passCalls == 1 && skipCalls == 1,
        "CLI02", "each selected stable test must execute exactly once");
  Check(summary.passed == 1 && summary.skipped == 1 && summary.ExitCode() == 0,
        "CLI05", "PASS/SKIP summary and successful exit policy are incorrect");

  const std::vector<const SEP::Testing::Descriptor*> routine =
      registry.Select({}, {}, true);
  Check(routine.size() == 1 && routine[0]->id == "A01",
        "CLI03", "--all-tests policy must include routine but not extended tests");
}

void TestCli04(const SEP::Testing::Registry& registry) {
  SEP::Util::CLI::Options options;
  std::string error;
  Check(!Parse({"sep", "--test"}, options, error),
        "CLI04", "missing --test value must fail parsing");
  options = SEP::Util::CLI::Options();
  Check(!Parse({"sep", "--list-tests", "--test", "A01"}, options, error),
        "CLI04", "listing and execution must be rejected as conflicting modes");
  options = SEP::Util::CLI::Options();
  Check(!Parse({"sep", "--all-tests", "--test-group=alpha"}, options, error),
        "CLI04", "--all-tests and explicit selectors must be rejected");
  options = SEP::Util::CLI::Options();
  Check(!Parse({"sep", "--input", "run.in", "--all-tests"}, options, error),
        "CLI04", "production input must not leak into component-test mode");
  options = SEP::Util::CLI::Options();
  Check(!Parse({"sep", "--initialization-only"}, options, error),
        "CLI04", "initialization-only mode must require an explicit input");
  options = SEP::Util::CLI::Options();
  Check(!Parse({"sep", "--reduced-shock-event", "event.in"}, options, error),
        "CLI04", "reduced coupling must require an explicit field-line input");
  options = SEP::Util::CLI::Options();
  Check(!Parse({"sep", "--input", "run.in",
                "--initialization-output-dir", "preview"}, options, error),
        "CLI04", "the initialization output directory must require its mode");

  bool rejectedUnknown = false;
  try {
    (void)registry.Select({"unknown"}, {}, false);
  }
  catch (const std::invalid_argument&) {
    rejectedUnknown = true;
  }
  Check(rejectedUnknown, "CLI04", "unknown stable IDs must fail before initialization");

  options = SEP::Util::CLI::Options();
  Check(Parse({"sep", "--all-tests", "--test-json", "results.json",
               "--test-junit=results.xml"}, options, error) &&
            options.testJsonPath == "results.json" &&
            options.testJunitPath == "results.xml",
        "CLI04", "structured-report paths must support both value syntaxes");
  options = SEP::Util::CLI::Options();
  Check(!Parse({"sep", "--test-json", "results.json"}, options, error),
        "CLI04", "a report path without a test selector must be rejected");

  // Application-level validation paths are parsed by the same production CLI
  // and are deliberately restricted to one explicit test. This prevents a
  // case input from leaking into an unrelated callback selected by a group.
  SEP::Util::CLI::Options validation;
  Check(Parse({"sep", "--test", "CV01", "--test-input", "native.args",
               "--test-output-dir", "evidence/CV01"}, validation, error) &&
            validation.testInputPath == "native.args" &&
            validation.testArtifactDirectory == "evidence/CV01",
        "CLI04", "end-to-end paths must reach one explicit linked test");
  SEP::Util::CLI::Options missingArtifacts;
  Check(!Parse({"sep", "--test", "CV01", "--test-input", "native.args"},
               missingArtifacts, error),
        "CLI04", "end-to-end input and artifact paths must be atomic");
  SEP::Util::CLI::Options ambiguousValidation;
  Check(!Parse({"sep", "--test", "CV01", "--test", "PARK01",
                "--test-input", "native.args", "--test-output-dir", "out"},
               ambiguousValidation, error),
        "CLI04", "one case-specific input cannot configure multiple tests");
}

void TestCli05LegacyAndErrors() {
  std::string error;
  SEP::Util::CLI::Options legacy;
  Check(Parse({"sep", "--run-test-manager", "--no-test-manager",
               "--testmanager=YES"}, legacy, error) && legacy.runTestManager,
        "CLI05", "legacy TestManager aliases/Boolean values/last-wins changed");
  Check(!SEP::Util::CLI::IsComponentTestExecutionRequested(legacy),
        "CLI05", "legacy TestManager must not imply new test-only mode");
  for(const auto& selection:std::vector<std::vector<std::string>>{{"--test","A01"},
      {"--test-group","alpha"},{"--tests","A01"},{"--all-tests"}}) {
    auto arguments=std::vector<std::string>{"sep","--run-test-manager"};
    arguments.insert(arguments.end(),selection.begin(),selection.end());
    SEP::Util::CLI::Options mixed;
    Check(!Parse(arguments,mixed,error)&&error.find("TestManager")!=std::string::npos,
          "CLI05","enabled legacy diagnostics were silently ignored by component-only execution");
  }
  SEP::Util::CLI::Options disabledLegacy;
  Check(Parse({"sep","--run-test-manager","--no-test-manager","--test","A01"},disabledLegacy,error)&&
        !disabledLegacy.runTestManager&&SEP::Util::CLI::IsComponentTestExecutionRequested(disabledLegacy),
        "CLI05","explicit last-occurrence legacy disable must remain compatible with component execution");

  SEP::Util::CLI::Options wp30;
  Check(Parse({"sep", "--total-iterations=42", "--shock-model", "analytical",
               "--cme-scenario=slow", "--field-line-seed-area", "12.5",
               "--swcme-failure-policy=diagnostic-fallback",
               "--swcme-fallback-density", "7e6",
               "--swcme-fallback-speed=350000",
               "--swcme-fallback-divergence", "-2e-6",
               "--swcme-override", "ambient.wind_speed=450 km/s",
               "--swcme-override=cme.launch_speed=1400 km/s",
               "--shock-turbulence-efficiency", "0.04",
               "--shock-turbulence-plus-fraction=0.7",
               "--merge-minimum", "10", "--merge-maximum=20"},
              wp30,error) && wp30.totalIterations==42 &&
            wp30.totalIterationsProvided && wp30.analyticalShock &&
            wp30.slowCmeScenario && wp30.fieldLineSeedAreaM2==12.5 &&
            wp30.shockTurbulenceEfficiency==0.04 &&
            wp30.shockTurbulencePlusFraction==0.7 &&
            wp30.swcmeFailurePolicy==
                SEP::Util::CLI::Options::SwcmeFailurePolicy::DiagnosticFallback &&
            wp30.swcmeFailurePolicyProvided &&
            wp30.swcmeFallbackDensityM3==7.0e6 &&
            wp30.swcmeFallbackSpeedMPerS==350000.0 &&
            wp30.swcmeFallbackDivergencePerS==-2.0e-6 &&
            wp30.swcmeOverrides.size()==2 &&
            wp30.mergeMinimum==10 && wp30.mergeMaximum==20,
        "CLI05", "WP30 run controls or source-layer markers did not parse");
  SEP::Util::CLI::Options invalidWp30;
  Check(!Parse({"sep", "--shock-turbulence-efficiency", "1.1"},
               invalidWp30,error),
        "CLI05", "WP30 out-of-range source efficiency must fail preflight");
  SEP::Util::CLI::Options invalidPolicy;
  Check(!Parse({"sep", "--swcme-failure-policy", "silently-zero"},
               invalidPolicy,error),
        "CLI05", "unknown SWCME recovery policy must fail preflight");
  SEP::Util::CLI::Options invalidOverride;
  Check(!Parse({"sep", "--swcme-override", "ambient.wind_speed"},
               invalidOverride,error),
        "CLI05", "SWCME override without key=value must fail preflight");

  int failCalls = 0;
  int errorCalls = 0;
  SEP::Testing::Registry outcomes({
      Descriptor("FAIL01", "fixture", SEP::Testing::RuntimeClass::Routine,
                 SEP::Testing::Status::Fail, &failCalls),
      Descriptor("ERROR01", "fixture", SEP::Testing::RuntimeClass::Routine,
                 SEP::Testing::Status::Error, &errorCalls),
  });
  std::ostringstream output;
  SEP::Testing::Summary summary = outcomes.Run(
      outcomes.Select({"FAIL01"}, {}, false), output);
  Check(summary.ExitCode() == 1, "CLI05", "FAIL must return process status 1");
  output.str("");
  summary = outcomes.Run(outcomes.Select({"ERROR01"}, {}, false), output);
  Check(summary.ExitCode() == 2, "CLI05", "ERROR must return process status 2");
}

void TestRegistryCompleteness() {
  int calls = 0;
  bool duplicateRejected = false;
  try {
    SEP::Testing::Registry duplicate({
        Descriptor("A01", "alpha", SEP::Testing::RuntimeClass::Routine,
                   SEP::Testing::Status::Pass, &calls),
        Descriptor("a01", "alpha", SEP::Testing::RuntimeClass::Routine,
                   SEP::Testing::Status::Pass, &calls),
    });
  }
  catch (const std::invalid_argument&) {
    duplicateRejected = true;
  }
  Check(duplicateRejected, "REGISTRY", "case-folded duplicate IDs must be rejected");

  bool incompleteRejected = false;
  try {
    SEP::Testing::Descriptor incomplete;
    incomplete.id = "EMPTY01";
    SEP::Testing::Registry registry({incomplete});
  }
  catch (const std::invalid_argument&) {
    incompleteRejected = true;
  }
  Check(incompleteRejected, "REGISTRY", "incomplete metadata must be rejected");
}

void TestRefinementOrder() {
  // A fourfold error decrease for a twofold resolution decrease is exactly
  // second order.  This independent arithmetic fixture also checks that the
  // helper rejects zero-error and reversed-resolution studies rather than
  // manufacturing a finite order for unusable evidence.
  const SEP::Testing::RefinementOrderEstimate secondOrder =
      SEP::Testing::EstimateRefinementOrder(0.04, 0.01, 0.2, 0.1);
  Check(secondOrder.valid &&
            std::fabs(secondOrder.observedOrder - 2.0) <= 1.0e-14 &&
            std::fabs(secondOrder.coarseToFineErrorRatio - 4.0) <= 1.0e-14 &&
            std::fabs(secondOrder.coarseToFineResolutionRatio - 2.0) <= 1.0e-14,
        "REFINE", "automatic observed-order calculation is incorrect");
  Check(!SEP::Testing::EstimateRefinementOrder(0.0, 0.01, 0.2, 0.1).valid &&
            !SEP::Testing::EstimateRefinementOrder(0.04, 0.01, 0.1, 0.2).valid,
        "REFINE", "degenerate refinement evidence must be rejected");
}

std::string ReadFile(const char* path) {
  std::ifstream input(path);
  return std::string((std::istreambuf_iterator<char>(input)),
                     std::istreambuf_iterator<char>());
}

void TestStructuredReportsAndHiddenFailures() {
  int calls = 0;
  SEP::Testing::Descriptor hidden = Descriptor(
      "HIDDEN01", "fixture", SEP::Testing::RuntimeClass::Routine,
      SEP::Testing::Status::Pass, &calls);
  hidden.callback = [&calls]() {
    ++calls;
    SEP::Testing::Result result;
    result.status = SEP::Testing::Status::Pass;
    result.message = "legacy callback retained one failed assertion";
    result.hasSeed = true;
    result.seed = 1300;
    result.configuration.push_back("fixture=hidden-failure");
    result.metrics.push_back(
        {"assertion_failures", 1.0, 0.0, "<=", "count"});
    result.artifacts.push_back("fixture-evidence.txt");
    return result;
  };
  SEP::Testing::Registry registry({hidden});
  const SEP::Testing::Summary summary = registry.Run(
      registry.Select({}, {}, true), std::cout);
  Check(summary.failed == 1 && summary.ExitCode() == 1,
        "HIDDEN", "a positive assertion_failures metric must override PASS");

  std::string error;
  Check(SEP::Testing::WriteJsonSummary(summary, "step13-results.json", &error),
        "REPORT", "JSON report writer failed: " + error);
  error.clear();
  Check(SEP::Testing::WriteJUnitSummary(summary, "step13-results.xml", &error),
        "REPORT", "JUnit report writer failed: " + error);
  const std::string json = ReadFile("step13-results.json");
  const std::string junit = ReadFile("step13-results.xml");
  Check(json.find("srcsep-component-tests-v1") != std::string::npos &&
            json.find("\"id\": \"HIDDEN01\"") != std::string::npos &&
            json.find("fixture-evidence.txt") != std::string::npos,
        "REPORT", "JSON report omitted schema, result ID, or evidence");
  Check(junit.find("<failure") != std::string::npos &&
            junit.find("name=\"HIDDEN01\"") != std::string::npos &&
            junit.find("assertion_failures=1") != std::string::npos,
        "REPORT", "JUnit report omitted failure or metric evidence");
  std::remove("step13-results.json");
  std::remove("step13-results.xml");
}

// CLI06: mover precedence --particle-mover > application.sep.transport (input
// file) > default fte-dmumu.  The Parker-only parallel-diffusion provider
// is the probe: it is valid for parker and invalid for the fte-dmumu default.
// The file itself is not read by ParseCommandLine; main.cpp applies it and
// then calls ValidateMoverDependentOptions(), which is exercised directly.
void TestCli06MoverPrecedence() {
  namespace CLI = SEP::Util::CLI;
  using SEP::Mover::ProductionMover;
  std::string error;

  CLI::Options noInput;
  Check(!Parse({"sep", "--spatial-diffusion-provider",
                "parallel-diffusion-library"}, noInput, error) &&
            error.find("unsupported mover/coefficient combination") !=
                std::string::npos,
        "CLI06", "without --input the default fte-dmumu must be validated "
                 "immediately and reject a Parker-only provider");

  CLI::Options explicitParker;
  Check(Parse({"sep", "--particle-mover", "parker",
               "--spatial-diffusion-provider", "parallel-diffusion-library"},
              explicitParker, error) &&
            !CLI::IsMoverPendingInputFile(explicitParker),
        "CLI06", "an explicit parker mover must accept the Parker-only provider");

  CLI::Options deferred;
  const bool deferredParsed =
      Parse({"sep", "--input", "examples/sep_parker_mesh_parallel_diffusion.in",
             "--spatial-diffusion-provider", "parallel-diffusion-library"},
            deferred, error);
  Check(deferredParsed && CLI::IsMoverPendingInputFile(deferred) &&
            !deferred.particleMoverProvided &&
            deferred.particleMover == ProductionMover::FocusedTransportDiffusion,
        "CLI06", "--input without --particle-mover must defer mover checks: " +
                 error);
  std::ostringstream rejected;
  Check(!CLI::ValidateMoverDependentOptions(deferred, rejected),
        "CLI06", "deferred validation must still reject the fte-dmumu default");
  deferred.particleMover = ProductionMover::Parker;  // file value applied
  std::ostringstream accepted;
  Check(CLI::ValidateMoverDependentOptions(deferred, accepted),
        "CLI06", "deferred validation must accept the file's parker mover: " +
                 accepted.str());

  CLI::Options cliOverridesFile;
  Check(!Parse({"sep", "--input",
                "examples/sep_parker_mesh_parallel_diffusion.in",
                "--particle-mover", "fte-dmumu",
                "--spatial-diffusion-provider", "parallel-diffusion-library"},
               cliOverridesFile, error),
        "CLI06", "an explicit --particle-mover must be validated at parse "
                 "time even with --input");
}

void TestCampaignOptionsAndSeeds() {
  // Active native legacy keys must never map to the core multiplicative
  // generator's even fixed state 2^30. Historical defaults remain untouched.
  Check(SEP::Util::CLI::LegacyTestSeedFromKey(0)==1&&
        SEP::Util::CLI::LegacyTestSeedFromKey(UINT64_C(1073741823))==UINT32_C(2147483647)&&
        SEP::Util::CLI::LegacyTestSeedFromKey(UINT64_C(1073741824))==1,
        "CLI02","legacy odd31 seed bridge changed declared domain endpoints");
  for(std::uint64_t key:{UINT64_C(0),UINT64_C(1),UINT64_C(1073741823),
      UINT64_C(1073741824),UINT64_C(2147483647),std::numeric_limits<std::uint64_t>::max()}) {
    const auto native=SEP::Util::CLI::LegacyTestSeedFromKey(key);
    Check(native>0&&(native&1)==1&&native<=UINT32_C(2147483647),"CLI02",
          "active legacy test seed enters an even fixed/short-period state");
  }
  namespace T=SEP::Testing;
  SEP::Util::CLI::Options options;std::string error;
  Check(Parse({"sep","--tests","id:A01,tag:geometry","--test-seed=18446744073709551615",
      "--test-repeat","3","--test-output","CaseSensitive/Results","--test-format=JuNiT",
      "--fail-fast","--keep-test-artifacts=on"},options,error)&&
      options.testExpressions.size()==1&&options.testSeedProvided&&
      options.testSeed==std::numeric_limits<std::uint64_t>::max()&&options.testRepeat==3&&
      options.testOutputDirectory=="CaseSensitive/Results"&&
      options.testFormat==SEP::Util::CLI::Options::TestFormat::JUnit&&
      options.testFailFast&&options.keepTestArtifacts,
      "CLI02","additive campaign controls did not reach the existing Options object");
  for(const std::vector<std::string>& args:{
      std::vector<std::string>{"sep","--test","A01","--test-seed","18446744073709551616"},
      {"sep","--test","A01","--test-seed","-1"},
      {"sep","--test","A01","--test-seed","1.5"},
      {"sep","--test","A01","--test-repeat","0"},
      {"sep","--test","A01","--test-repeat","4294967296"},
      {"sep","--test","A01","--test-format","yaml"},
      {"sep","--test-seed","42"},
      {"sep","--list-tests","--fail-fast"}}) {
    SEP::Util::CLI::Options invalid;
    Check(!Parse(args,invalid,error),"CLI03","malformed/unscoped campaign request was accepted");
  }
  SEP::Util::CLI::Options scalar;
  Check(Parse({"sep","--test","A01","--test-seed","1","--test-seed=2",
      "--test-repeat=2","--test-repeat","4","--fail-fast","--fail-fast=off",
      "--keep-test-artifacts","no"},scalar,error)&&scalar.testSeed==2&&
      scalar.testRepeat==4&&!scalar.testFailFast&&!scalar.keepTestArtifacts,
      "CLI02","campaign scalar last-occurrence-wins/Boolean compatibility failed");

  int draws=0,fixedCalls=0,fakeCalls=0;
  auto seeded=Descriptor("A01","alpha",T::RuntimeClass::Routine,T::Status::Pass,&draws);
  seeded.tags={"geometry","mover:parker"};
  seeded.seedCapability=T::SeedCapability::Context;seeded.repeatSupported=true;
  seeded.callback=[&draws]() {
    ++draws;const auto seed=T::TestSeed(41);
    SEP::Transport::KeyedRandomStream random(seed,17,3,0);
    const double observed=random.UniformOpen01();
    const auto directory=T::GetExecutionContext().artifactDirectory;
    for(std::size_t length=directory.find('/');;length=directory.find('/',length+1)) {
      const std::string prefix=length==std::string::npos?directory:directory.substr(0,length);
      if(!prefix.empty()&&mkdir(prefix.c_str(),0700)!=0&&errno!=EEXIST)
        throw std::runtime_error("cannot create actual draw fixture directory");
      if(length==std::string::npos)break;
    }
    const auto artifact=directory+"/draw.txt";
    std::ofstream file(artifact.c_str());file<<std::setprecision(17)<<observed<<'\n';file.close();
    if(file.fail())throw std::runtime_error("cannot write actual draw fixture artifact");
    T::Result result;result.status=T::Status::Pass;result.message="actual keyed draw";
    result.hasSeed=true;result.seed=seed;
    result.metrics.push_back({"observed_uniform",observed,0,"observed","dimensionless"});
    result.artifacts.push_back(artifact);
    return result;
  };
  auto fixed=Descriptor("B01","beta",T::RuntimeClass::Extended,T::Status::Pass,&fixedCalls);
  auto fake=Descriptor("C01","gamma",T::RuntimeClass::Routine,T::Status::Pass,&fakeCalls);
  fake.seedCapability=T::SeedCapability::Context;fake.repeatSupported=true;
  const T::Registry registry({seeded,fixed,fake});
  Check(registry.Select({}, {},false,{"id:A01,tag:geometry"}).size()==1,
        "CLI02","expression union did not deduplicate an explicit ID and tag");
  Check(registry.Select({}, {},false,{"mover:parker"}).front()->id=="A01"&&
        registry.Select({}, {},false,{"all"}).size()==3,
        "CLI02","tag or explicit complete-suite selector did not resolve");
  bool rejected=false;
  try{registry.Select({}, {},false,{"tag:missing"});}catch(const std::invalid_argument&){rejected=true;}
  Check(rejected,"CLI03","unknown tag was silently ignored");
  T::CampaignOptions campaign;campaign.seedProvided=true;campaign.campaignSeed=123;
  campaign.repeatCount=3;campaign.artifactRoot="campaign-artifacts";
  std::ostringstream log;
  const auto first=registry.RunCampaign(registry.Select({"A01"},{},false),campaign,log);
  const auto second=registry.RunCampaign(registry.Select({}, {},false,{"tag:geometry"}),campaign,log);
  Check(first.results.size()==1&&first.results[0].invocations.size()==3&&draws==6,
        "CLI02","repetition lost a stable aggregate ID or an invocation");
  if(first.results.size()==1&&first.results[0].invocations.size()==3&&second.results.size()==1) {
    for(unsigned i=0;i<3;++i) {
      const auto& a=*first.results[0].invocations[i].result;
      const auto& b=*second.results[0].invocations[i].result;
      Check(a.seed==b.seed&&a.metrics[0].value==b.metrics[0].value,
            "CLI02","same campaign/test/replicate changed actual draws with selector form");
      Check(a.artifacts[0].find("replicate-"+std::to_string(i))!=std::string::npos,
            "REPORT","replicates did not receive distinct artifact namespaces");
      std::ifstream evidence(a.artifacts[0].c_str());double observed=0;evidence>>observed;
      Check(evidence.good()&&observed==a.metrics[0].value,"REPORT",
            "reported draw artifact did not retain the actual generated variate");
    }
    Check(first.results[0].invocations[0].result->seed!=first.results[0].invocations[1].result->seed&&
          first.results[0].invocations[0].result->metrics[0].value!=first.results[0].invocations[1].result->metrics[0].value,
          "CLI02","distinct replicates changed metadata but replayed the same actual draw");
  }
  campaign.repeatCount=1;
  const auto unsupported=registry.RunCampaign(registry.Select({"B01"},{},false),campaign,log);
  Check(unsupported.errors==1&&fixedCalls==0,"CLI03","fixed-contract seed override executed or falsely passed");
  const auto metadataOnly=registry.RunCampaign(registry.Select({"C01"},{},false),campaign,log);
  Check(metadataOnly.errors==1&&fakeCalls==1,"CLI04","metadata-only seeded callback was accepted without applying a stream seed");

  int a=0,b=0,d=0;
  auto pass=Descriptor("A","g",T::RuntimeClass::Routine,T::Status::Pass,&a);
  auto fail=Descriptor("B","g",T::RuntimeClass::Routine,T::Status::Fail,&b);
  auto untouched=Descriptor("D","g",T::RuntimeClass::Routine,T::Status::Pass,&d);
  pass.seedCapability=fail.seedCapability=untouched.seedCapability=T::SeedCapability::Deterministic;
  pass.repeatSupported=fail.repeatSupported=untouched.repeatSupported=true;
  const T::Registry stopRegistry({untouched,fail,pass});
  campaign.repeatCount=3;campaign.failFast=true;
  const auto stopped=stopRegistry.RunCampaign(stopRegistry.Select({}, {},true),campaign,log);
  Check(a==3&&b==1&&d==0&&stopped.failed==1&&stopped.notRun.size()==2&&stopped.notRunInvocationCount==5,
        "CLI04","fail-fast did not preserve completed replicates and explicit not-run evidence");
  Check(T::WriteJsonSummary(first,"campaign-results.json",&error)&&
        T::WriteJUnitSummary(first,"campaign-results.xml",&error),
        "REPORT","nested replicate reports could not be written transactionally");
  const auto json=ReadFile("campaign-results.json"),xml=ReadFile("campaign-results.xml");
  Check(json.find("\"replicate_index\": 0")!=std::string::npos&&
        json.find("\"replicate_index\": 2")!=std::string::npos&&
        json.find(std::to_string(first.results.front().invocations.back().result->seed))!=std::string::npos&&
        xml.find("tests=\"3\"")!=std::string::npos&&
        xml.find("A01[replicate=2]")!=std::string::npos,
        "REPORT","JSON/JUnit dropped replicate seed, index or unique testcase evidence");
  Check(!T::WriteJsonSummary(first,"campaign-results.json/missing",&error)&&
        ReadFile("campaign-results.json")==json,
        "REPORT","failed report creation damaged previously committed evidence");
  campaign.repeatCount=std::numeric_limits<std::uint32_t>::max();
  const auto bounded=stopRegistry.RunCampaign(stopRegistry.Select({"B"},{},false),campaign,log);
  Check(bounded.results.size()==1&&bounded.notRun.size()==1&&
        bounded.notRunInvocationCount==std::numeric_limits<std::uint32_t>::max()-UINT64_C(1),
        "CLI04","fail-fast materialized billions of skipped repetitions or lost their count");
  campaign.failFast=false;
  const auto unsupportedHuge=registry.RunCampaign(registry.Select({"B01"},{},false),campaign,log);
  Check(unsupportedHuge.errors==1&&unsupportedHuge.results.size()==1&&
        unsupportedHuge.results[0].invocations.size()==1&&fixedCalls==0&&
        unsupportedHuge.notRun==std::vector<std::string>{"B01#1..4294967294"}&&
        unsupportedHuge.notRunInvocationCount==UINT64_C(4294967294),
        "CLI04","unsupported UINT32_MAX repeat allocated duplicate errors or lost compact non-execution evidence");
  campaign.repeatCount=2;
  const auto fixedThenDeterministic=stopRegistry.RunCampaign({&fixed,&pass},campaign,log);
  Check(fixedThenDeterministic.errors==1&&fixedThenDeterministic.passed==1&&a==5&&
        fixedThenDeterministic.notRunInvocationCount==1,
        "CLI04","unsupported-contract rejection incorrectly stopped later supported descriptors without fail-fast");
}

}  // namespace

int main() {
  int passCalls = 0;
  int skipCalls = 0;
  SEP::Testing::Registry registry({
      Descriptor("B01", "beta", SEP::Testing::RuntimeClass::Extended,
                 SEP::Testing::Status::Skip, &skipCalls),
      Descriptor("A01", "alpha", SEP::Testing::RuntimeClass::Routine,
                 SEP::Testing::Status::Pass, &passCalls),
  });

  TestCli01(registry);
  TestCli02And03(registry, passCalls, skipCalls);
  TestCli04(registry);
  TestCli05LegacyAndErrors();
  TestRegistryCompleteness();
  TestRefinementOrder();
  TestStructuredReportsAndHiddenFailures();
  TestCli06MoverPrecedence();
  TestCampaignOptionsAndSeeds();

  if (failures != 0) {
    std::cerr << "Step 1 focused summary: PASS=0 FAIL=" << failures
              << " SKIP=0 ERROR=0 (unexpected outer assertion failures)\n";
    return EXIT_FAILURE;
  }
  std::cout << "CLI01 PASS\nCLI02 PASS\nCLI03 PASS\nCLI04 PASS\nCLI05 PASS\n"
            << "CLI06 PASS\n"
            << "REGISTRY PASS\nREFINE PASS\nHIDDEN PASS\nREPORT PASS\n"
            << "Negative control HIDDEN01: EXPECTED INTERNAL FAIL observed; "
               "the framework correctly rejected a callback that falsely "
               "reported PASS.\n"
            << "Step 1 focused summary: PASS=10 FAIL=0 SKIP=0 ERROR=0 "
               "(HIDDEN01 is expected behavior, not an outer failure)\n";
  return EXIT_SUCCESS;
}

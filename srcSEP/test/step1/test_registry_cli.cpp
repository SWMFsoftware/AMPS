#include "../../util/sep_cli.h"
#include "../../util/sep_test_registry.h"

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

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
}

void TestCli05LegacyAndErrors() {
  std::string error;
  SEP::Util::CLI::Options legacy;
  Check(Parse({"sep", "--run-test-manager", "--no-test-manager",
               "--testmanager=YES"}, legacy, error) && legacy.runTestManager,
        "CLI05", "legacy TestManager aliases/Boolean values/last-wins changed");
  Check(!SEP::Util::CLI::IsComponentTestExecutionRequested(legacy),
        "CLI05", "legacy TestManager must not imply new test-only mode");

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

  if (failures != 0) {
    std::cerr << "Step 1 focused tests: " << failures << " failure(s)\n";
    return EXIT_FAILURE;
  }
  std::cout << "CLI01 PASS\nCLI02 PASS\nCLI03 PASS\nCLI04 PASS\nCLI05 PASS\n"
            << "REGISTRY PASS\nREFINE PASS\nHIDDEN PASS\nREPORT PASS\n";
  return EXIT_SUCCESS;
}

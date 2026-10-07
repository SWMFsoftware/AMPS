// Native executable selectors and production CLI parsing. This translation
// unit has no AMPS/MPI/SWCME dependency: selection is validated before any
// model setup or mesh allocation. Suite membership belongs to the native
// descriptor registry, rather than a copied list of test IDs in the CLI.
#include "configuration_io.h"
#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <limits>
#include <map>

namespace SEP3D { namespace RuntimeModel {
namespace {
Core::Status Invalid(const std::string& message) {
  return Core::Status(Core::StatusCode::InvalidInput, message);
}
std::string Lower(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
      [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return value;
}
bool ParseUnsigned64(const std::string& text, std::uint64_t* value) {
  if (value == nullptr || text.empty() || text[0] == '-') return false;
  errno = 0;
  char* end = nullptr;
  const unsigned long long parsed = std::strtoull(text.c_str(), &end, 10);
  if (errno == ERANGE || end == text.c_str() || *end != '\0') return false;
  *value = static_cast<std::uint64_t>(parsed);
  return true;
}
template <typename Enum>
bool ParseEnum(const std::string& text,
               const std::map<std::string, Enum>& values, Enum* result) {
  const auto found = values.find(Lower(text));
  if (found == values.end()) return false;
  *result = found->second;
  return true;
}
} // namespace

const char* Name(LogVerbosity value) {
  switch (value) {
    case LogVerbosity::Quiet: return "quiet";
    case LogVerbosity::Normal: return "normal";
    case LogVerbosity::Verbose: return "verbose";
  }
  return "unknown";
}

Core::Status ParseStandaloneCommandLine(
    int argc, char* const argv[], StandaloneCommandLine* result) {
  if (result == nullptr) return Invalid("command-line output is null");
  StandaloneCommandLine candidate;
  for (int i = 1; i < argc; ++i) {
    const std::string argument(argv[i]);
    auto requireValue = [&](const char* option, std::string* value) {
      if (i + 1 >= argc) return false;
      *value = argv[++i];
      return !value->empty() && value->front() != '-' && *value != option;
    };
    if (argument == "-input" || argument == "--input") {
      if (!candidate.inputPath.empty())
        return Invalid("-input/--input may be specified only once");
      if (!requireValue("--input", &candidate.inputPath))
        return Invalid(argument + " requires a path");
      candidate.sectionInput = argument == "-input";
    } else if (argument == "--initialization-only") {
      candidate.initializationOnly = true;
    } else if (argument == "--initialization-output-dir") {
      if (!requireValue("--initialization-output-dir",
                        &candidate.initializationOutputDirectory)) {
        return Invalid("--initialization-output-dir requires a path");
      }
    } else if (argument == "--output-dir") {
      if (!requireValue("--output-dir", &candidate.outputDirectoryOverride))
        return Invalid("--output-dir requires a path");
    } else if (argument == "--restart") {
      if (!requireValue("--restart", &candidate.restartPath))
        return Invalid("--restart requires a path");
    } else if (argument == "--dry-run") {
      candidate.dryRun = true;
    } else if (argument == "--list-tests") {
      candidate.listTests = true;
    } else if (argument == "--all-tests") {
      candidate.allTests = true;
    } else if (argument == "--test-suite") {
      if (!candidate.testSuite.empty())
        return Invalid("--test-suite may be specified only once");
      if (!requireValue("--test-suite", &candidate.testSuite))
        return Invalid("--test-suite requires a suite name: sep-corona");
      candidate.testSuite = Lower(candidate.testSuite);
      if (candidate.testSuite != "sep-corona")
        return Invalid("unknown native test suite '" + candidate.testSuite + "'");
    } else if (argument == "--test") {
      std::string id;
      if (!requireValue("--test", &id)) return Invalid("--test requires an ID");
      candidate.tests.push_back(id);
    } else if (argument == "--test-input") {
      if (!requireValue("--test-input", &candidate.testInputPath))
        return Invalid("--test-input requires a path");
    } else if (argument == "--test-json") {
      if (!requireValue("--test-json", &candidate.testJsonPath))
        return Invalid("--test-json requires a path");
    } else if (argument == "--artifact-directory") {
      if (!requireValue("--artifact-directory",
                        &candidate.testArtifactDirectory))
        return Invalid("--artifact-directory requires a path");
    } else if (argument == "--test-steps") {
      std::string value;
      if (!requireValue("--test-steps", &value) ||
          !ParseUnsigned64(value, &candidate.testSteps))
        return Invalid("--test-steps requires a non-negative integer");
    } else if (argument == "--expect-mpi-ranks") {
      std::string value;
      std::uint64_t ranks = 0;
      if (!requireValue("--expect-mpi-ranks", &value) ||
          !ParseUnsigned64(value, &ranks) || ranks == 0 ||
          ranks > static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
        return Invalid("--expect-mpi-ranks requires a positive integer");
      candidate.expectedMpiRanks = static_cast<int>(ranks);
    } else if (argument == "--log-level") {
      std::string value;
      if (!requireValue("--log-level", &value) ||
          !ParseEnum(value, {{"quiet", LogVerbosity::Quiet},
                             {"normal", LogVerbosity::Normal},
                             {"verbose", LogVerbosity::Verbose}},
                     &candidate.verbosity)) {
        return Invalid("--log-level accepts quiet, normal, or verbose");
      }
    } else {
      return Invalid("unknown srcSEP3D option '" + argument + "'");
    }
  }
  const int selectionModes = static_cast<int>(candidate.listTests) +
      static_cast<int>(candidate.allTests) +
      static_cast<int>(!candidate.testSuite.empty()) +
      static_cast<int>(!candidate.tests.empty());
  if (selectionModes > 1)
    return Invalid("--list-tests, --all-tests, --test-suite, and --test are mutually exclusive");
  if (candidate.initializationOnly && candidate.dryRun) {
    return Invalid("--initialization-only and --dry-run are mutually exclusive: "
                   "the former builds the AMPS mesh, while the latter forbids "
                   "AMPS allocation");
  }
  if (!candidate.initializationOutputDirectory.empty() &&
      !candidate.initializationOnly) {
    return Invalid("--initialization-output-dir requires --initialization-only");
  }
  if (candidate.initializationOnly && selectionModes != 0) {
    return Invalid("--initialization-only cannot be combined with test selection");
  }
  const bool nativeTestRun = candidate.allTests ||
      !candidate.testSuite.empty() || !candidate.tests.empty();
  const bool hasNativeTestOptions = !candidate.testInputPath.empty() ||
      candidate.testJsonPath != "test_output/native/native.json" ||
      candidate.testArtifactDirectory != "test_output/native/artifacts" ||
      candidate.testSteps != 1 || candidate.expectedMpiRanks != 0;
  if (candidate.listTests && hasNativeTestOptions)
    return Invalid("--list-tests does not initialize AMPS and cannot be "
                   "combined with native-test execution options");
  if (!nativeTestRun && !candidate.listTests && hasNativeTestOptions)
    return Invalid("--test-input, --test-json, --artifact-directory, "
                   "--test-steps, and --expect-mpi-ranks require --test or "
                   "--all-tests or --test-suite");
  if (nativeTestRun) {
    if (candidate.sectionInput)
      return Invalid("linked native tests require the maintained --test-input/--input schema deck, not -input section mode");
    if (!candidate.inputPath.empty() && !candidate.testInputPath.empty() &&
        candidate.inputPath != candidate.testInputPath)
      return Invalid("--input and --test-input name different decks");
    if (candidate.testInputPath.empty())
      candidate.testInputPath = candidate.inputPath;
    if (candidate.testInputPath.empty())
      return Invalid("a linked native test requires --test-input PATH");
    if (candidate.dryRun)
      return Invalid("native tests require AMPS initialization and cannot be "
                     "combined with --dry-run");
  }
  // The shared runtime-file convention is intentionally useful without any
  // command-line arguments: batch systems can stage ``amps.in`` in each run
  // directory.  Native-test modes retain their explicit reviewed schema deck
  // because they must never inherit an unrelated working-directory file.
  if (candidate.inputPath.empty() && selectionModes == 0 &&
      !candidate.listTests) {
    candidate.inputPath = "amps.in";
    candidate.sectionInput = true;
  }
  *result = candidate;
  return Core::Status::OK();
}

} } // namespace SEP3D::RuntimeModel

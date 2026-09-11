#include "test_framework.hpp"

#include <cstring>
#include <exception>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

void test_cfg01(swcme_test::Context& context);
void test_cfg02(swcme_test::Context& context);
void test_den01(swcme_test::Context& context);
void test_par01(swcme_test::Context& context);
void test_par02(swcme_test::Context& context);
void test_par03(swcme_test::Context& context);
void test_1d_ambient_at_one_au(swcme_test::Context& context);

namespace {

// CLI modes are mutually exclusive. Keeping parsing separate from execution
// makes invalid combinations fail before any validation changes process state.
enum class Mode { All, Single, List, Help };

struct CommandLine {
  Mode mode = Mode::All;
  const char* test_id = nullptr;
};

// The registry is the single source used by --list, --all, and --test. Its
// declaration order is therefore also the deterministic --all execution order.
const swcme_test::TestCase kTests[] = {
    {"CFG01", "COMMON", "Configuration rejection and physical-range validation",
     test_cfg01},
    {"CFG02", "COMMON", "Unit-conversion and dimensional-consistency test",
     test_cfg02},
    {"DEN01", "COMMON",
     "Leblanc density normalization at the reference distance", test_den01},
    {"PAR01", "3D", "Parker-spiral equatorial vector orientation", test_par01},
    {"PAR02", "3D", "Parker field at arbitrary latitude and rotation axis", test_par02},
    {"PAR03", "3D", "Parker-field polar-limit regularity", test_par03},
    {"1D_AMBIENT_01", "1D", "Ambient values at 1 AU",
     test_1d_ambient_at_one_au},
};

constexpr std::size_t kTestCount = sizeof(kTests) / sizeof(kTests[0]);

void print_help(const char* executable) {
  std::cout
      << "usage: " << executable << " [--all | --list | --help | --test TEST_ID]\n\n"
      << "Options:\n"
      << "  --all           Run every registered test in registry order.\n"
      << "  --list          List registered tests without running them.\n"
      << "  --test TEST_ID  Run exactly one registered test.\n"
      << "  --help          Show this help text.\n\n"
      << "No arguments are equivalent to --all.\n\n"
      << "Examples:\n"
      << "  " << executable << " --list\n"
      << "  " << executable << " --all\n"
      << "  " << executable << " --test CFG02\n";
}

void print_registry() {
  std::cout << "Registered SWCME validation tests\n\n"
            << std::left << std::setw(16) << "ID"
            << std::setw(14) << "Class" << "Name\n"
            << "--------------------------------------------------------------------------\n";
  for (const auto& test : kTests) {
    std::cout << std::left << std::setw(16) << test.id
              << std::setw(14) << test.classification << test.name << '\n';
  }
  std::cout << "\nTotal: " << kTestCount << " tests\n";
}

bool parse_command_line(int argc, char** argv, CommandLine& command) {
  if (argc == 1) {
    command.mode = Mode::All;
    return true;
  }

  if (argc == 2) {
    if (std::strcmp(argv[1], "--all") == 0) {
      command.mode = Mode::All;
      return true;
    }
    if (std::strcmp(argv[1], "--list") == 0) {
      command.mode = Mode::List;
      return true;
    }
    if (std::strcmp(argv[1], "--help") == 0) {
      command.mode = Mode::Help;
      return true;
    }
    if (std::strcmp(argv[1], "--test") == 0) {
      std::cerr << "ERROR: --test requires a TEST_ID.\n"
                << "Use --list to display registered tests.\n";
      return false;
    }

    std::cerr << "ERROR: unknown option '" << argv[1] << "'.\n"
              << "Use --help for supported options.\n";
    return false;
  }

  if (argc == 3 && std::strcmp(argv[1], "--test") == 0) {
    command.mode = Mode::Single;
    command.test_id = argv[2];
    return true;
  }

  // More than one mode is deliberately rejected instead of applying an
  // order-dependent interpretation to combinations such as --all --test.
  std::cerr << "ERROR: command-line modes cannot be combined.\n"
            << "Use --help for supported syntax.\n";
  return false;
}

const swcme_test::TestCase* find_test(const char* id) {
  for (const auto& test : kTests) {
    if (std::strcmp(id, test.id) == 0) {
      return &test;
    }
  }
  return nullptr;
}

enum class Status { Pass, Fail, Skip };

struct ExecutionResult {
  const swcme_test::TestCase* test;
  Status status;
};

ExecutionResult run_test(const swcme_test::TestCase& test) {
  std::cout << "[ RUN      ] " << test.id << " [" << test.classification
            << "] " << test.name << '\n';
  swcme_test::Context context;

  try {
    test.function(context);
  } catch (const std::exception& error) {
    context.expect_true(false, std::string("unexpected exception: ") + error.what());
  } catch (...) {
    context.expect_true(false, "unexpected non-standard exception");
  }

  Status status = Status::Pass;
  if (context.failures() != 0) {
    status = Status::Fail;
    std::cout << "[  FAILED  ] " << test.id << '\n';
  } else if (context.passes() == 0 && context.skips() != 0) {
    status = Status::Skip;
    std::cout << "[  SKIPPED ] " << test.id << '\n';
  } else {
    std::cout << "[       OK ] " << test.id << '\n';
  }

  return {&test, status};
}

const char* status_name(Status status) {
  if (status == Status::Pass) return "PASS";
  if (status == Status::Fail) return "FAIL";
  return "SKIP";
}

int print_summary(const std::vector<ExecutionResult>& results) {
  int passed = 0;
  int failed = 0;
  int skipped = 0;
  for (const auto& result : results) {
    if (result.status == Status::Pass) ++passed;
    if (result.status == Status::Fail) ++failed;
    if (result.status == Status::Skip) ++skipped;
  }

  std::cout << "============================================================\n"
            << "SWCME validation summary\n"
            << "============================================================\n";
  for (const auto& result : results) {
    std::cout << std::left << std::setw(16) << result.test->id
              << std::setw(11) << result.test->classification
              << std::setw(58) << result.test->name
              << status_name(result.status) << '\n';
  }
  std::cout << "------------------------------------------------------------\n"
            << "Tests run:     " << results.size() << '\n'
            << "Passed:        " << passed << '\n'
            << "Failed:        " << failed << '\n'
            << "Skipped:       " << skipped << '\n'
            << "RESULT: " << (failed == 0 ? "PASS" : "FAIL") << '\n'
            << "============================================================\n";
  return failed == 0 ? 0 : 1;
}

}  // namespace

int main(int argc, char** argv) {
  CommandLine command;
  if (!parse_command_line(argc, argv, command)) {
    return 2;
  }

  if (command.mode == Mode::Help) {
    print_help(argv[0]);
    return 0;
  }
  if (command.mode == Mode::List) {
    print_registry();
    return 0;
  }

  std::vector<ExecutionResult> results;
  if (command.mode == Mode::Single) {
    const swcme_test::TestCase* test = find_test(command.test_id);
    if (test == nullptr) {
      std::cerr << "ERROR: unknown test ID '" << command.test_id << "'.\n"
                << "Use --list to display registered tests.\n";
      return 2;
    }
    results.push_back(run_test(*test));
  } else {
    for (const auto& test : kTests) {
      results.push_back(run_test(test));
    }
  }

  return print_summary(results);
}

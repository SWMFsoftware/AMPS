#include "test_framework.hpp"

#include <cstring>
#include <exception>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

void test_cfg01(swcme_test::Context& context);
void test_cfg02(swcme_test::Context& context);
void test_pst01(swcme_test::Context& context);
void test_pst02(swcme_test::Context& context);
void test_pst03(swcme_test::Context& context);
void test_pst06(swcme_test::Context& context);
void test_pst04(swcme_test::Context& context);
void test_out02(swcme_test::Context& context);
void test_out03(swcme_test::Context& context);
void test_out05(swcme_test::Context& context);
void test_out04(swcme_test::Context& context);
void test_out06(swcme_test::Context& context);
void test_out01(swcme_test::Context& context);
void test_out07(swcme_test::Context& context);
void test_pst07(swcme_test::Context& context);
void test_pst05(swcme_test::Context& context);
void test_pst08(swcme_test::Context& context);
void test_out08(swcme_test::Context& context);
void test_cfg03(swcme_test::Context& context);
void test_cfg04(swcme_test::Context& context);
void test_den01(swcme_test::Context& context);
void test_par01(swcme_test::Context& context);
void test_par02(swcme_test::Context& context);
void test_par03(swcme_test::Context& context);
void test_geo01(swcme_test::Context& context);
void test_geo02(swcme_test::Context& context);
void test_geo03(swcme_test::Context& context);
void test_geo04(swcme_test::Context& context);
void test_geo05(swcme_test::Context& context);
void test_geo06(swcme_test::Context& context);
void test_geo07(swcme_test::Context& context);
void test_geo08(swcme_test::Context& context);
void test_msh01(swcme_test::Context& context);
void test_msh02(swcme_test::Context& context);
void test_msh03(swcme_test::Context& context);
void test_msh04(swcme_test::Context& context);
void test_msh05(swcme_test::Context& context);
void test_kin01(swcme_test::Context& context);
void test_kin02(swcme_test::Context& context);
void test_kin03(swcme_test::Context& context);
void test_kin04(swcme_test::Context& context);
void test_kin05(swcme_test::Context& context);
void test_kin06(swcme_test::Context& context);
void test_kin07(swcme_test::Context& context);
void test_kin08(swcme_test::Context& context);
void test_kin09(swcme_test::Context& context);
void test_shk01(swcme_test::Context& context);
void test_shk02(swcme_test::Context& context);
void test_shk03(swcme_test::Context& context);
void test_shk04(swcme_test::Context& context);
void test_shk05(swcme_test::Context& context);
void test_shk06(swcme_test::Context& context);
void test_shk07(swcme_test::Context& context);
void test_shk08(swcme_test::Context& context);
void test_shk09(swcme_test::Context& context);
void test_shk10(swcme_test::Context& context);
void test_shk11(swcme_test::Context& context);
void test_shk12(swcme_test::Context& context);
void test_shk13(swcme_test::Context& context);
void test_shk14(swcme_test::Context& context);
void test_con01(swcme_test::Context& context);
void test_con02(swcme_test::Context& context);
void test_con03(swcme_test::Context& context);
void test_con04(swcme_test::Context& context);
void test_con05(swcme_test::Context& context);
void test_con06(swcme_test::Context& context);
void test_con07(swcme_test::Context& context);
void test_con08(swcme_test::Context& context);
void test_con09(swcme_test::Context& context);
void test_con10(swcme_test::Context& context);
void test_reg01(swcme_test::Context& context);
void test_reg02(swcme_test::Context& context);
void test_reg03(swcme_test::Context& context);
void test_reg04(swcme_test::Context& context);
void test_reg05(swcme_test::Context& context);
void test_acc01(swcme_test::Context& context);
void test_acc02(swcme_test::Context& context);
void test_acc03(swcme_test::Context& context);
void test_acc04(swcme_test::Context& context);
void test_acc05(swcme_test::Context& context);
void test_err01(swcme_test::Context& context);
void test_err02(swcme_test::Context& context);
void test_err03(swcme_test::Context& context);
void test_err04(swcme_test::Context& context);
void test_err05(swcme_test::Context& context);
void test_div01(swcme_test::Context& context);
void test_div02(swcme_test::Context& context);
void test_div03(swcme_test::Context& context);
void test_def01(swcme_test::Context& context);
void test_def02(swcme_test::Context& context);
void test_def03(swcme_test::Context& context);
void test_def04(swcme_test::Context& context);
void test_sep01(swcme_test::Context& context);
void test_sep02(swcme_test::Context& context);
void test_sep03(swcme_test::Context& context);
void test_sep04(swcme_test::Context& context);
void test_sep05(swcme_test::Context& context);
void test_sep06(swcme_test::Context& context);
void test_1d_ambient_at_one_au(swcme_test::Context& context);
void test_1d3d01(swcme_test::Context& context);
void test_1d3d02(swcme_test::Context& context);
void test_1d3d03(swcme_test::Context& context);

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
    {"PST01", "COMMON", "Prepared-state immutability", test_pst01},
    {"PST02", "COMMON", "Cross-model prepared-state rejection", test_pst02},
    {"PST03", "COMMON", "Cross-configuration prepared-state rejection", test_pst03},
    {"PST06", "COMMON", "Prepared-state record integrity", test_pst06},
    {"PST04", "COMMON", "Concurrent prepared-state evaluation", test_pst04},
    {"OUT02", "OUTPUT", "Write failure detection and propagation", test_out02},
    {"OUT03", "OUTPUT", "Transactional output commit", test_out03},
    {"OUT05", "OUTPUT", "Model-domain output preflight", test_out05},
    {"OUT04", "OUTPUT", "Box specification validation", test_out04},
    {"OUT06", "OUTPUT", "Mesh output validation", test_out06},
    {"OUT01", "OUTPUT", "Independent output parsing", test_out01},
    {"OUT07", "OUTPUT", "Demonstration program execution", test_out07},
    {"PST07", "INTEGRATION", "AMPS adapter equivalence", test_pst07},
    {"PST05", "COMMON", "Prepared-state lifetime contract", test_pst05},
    {"PST08", "PERFORMANCE", "State ownership performance", test_pst08},
    {"OUT08", "OUTPUT", "Strict warning writer build", test_out08},
    {"CFG03", "COMMON", "Smoothing width policy", test_cfg03},
    {"CFG04", "COMMON", "Configured radius domain", test_cfg04},
    {"KIN09", "COMMON", "Kinematic extrapolation domain", test_kin09},
    {"CON09", "3D", "Connectivity resolution-limit contract", test_con09},
    {"CON10", "3D", "Observer domain classification", test_con10},
    {"CFG01", "COMMON", "Configuration rejection and physical-range validation",
     test_cfg01},
    {"CFG02", "COMMON", "Unit-conversion and dimensional-consistency test",
     test_cfg02},
    {"DEF01", "1D<->3D", "Canonical shared default-configuration equivalence", test_def01},
    {"DEF02", "COMMON", "Science-scope, geometry, and Parker conventions", test_def02},
    {"DEF03", "COMMON", "Observer-local pre-shock model-scope validity", test_def03},
    {"DEF04", "COMMON", "Deterministic complete resolved-configuration manifest", test_def04},
    {"SEP01", "INTEGRATION", "AMPS background adapter matches direct production query", test_sep01},
    {"SEP02", "COMMON", "SEP source spectrum units and DSA slope conversion", test_sep02},
    {"SEP03", "1D<->3D", "AMPS-facing 1-D/3-D source-record identity", test_sep03},
    {"SEP04", "3D", "Shock-surface source patch area normalization", test_sep04},
    {"SEP05", "3D", "Observer cobpoint-to-SEP-source consistency", test_sep05},
    {"SEP06", "COMMON", "Resolved-compression source-disable contract", test_sep06},
    {"DEN01", "COMMON",
     "Leblanc density normalization at the reference distance", test_den01},
    {"PAR01", "3D", "Parker-spiral equatorial vector orientation", test_par01},
    {"PAR02", "3D", "Parker field at arbitrary latitude and rotation axis", test_par02},
    {"PAR03", "3D", "Parker-field polar-limit regularity", test_par03},
    {"GEO01", "3D", "Spherical shock geometry reference", test_geo01},
    {"GEO02", "3D", "True SSE front at the CME apex", test_geo02},
    {"GEO03", "3D", "SSE tangent-flank boundary behavior", test_geo03},
    {"GEO04", "3D", "Strict enforcement of finite SSE angular width", test_geo04},
    {"GEO05", "3D", "SSE surface level-set residual", test_geo05},
    {"GEO06", "3D", "Analytical shock normal versus level-set gradient", test_geo06},
    {"GEO07", "3D", "Shock normal speed versus finite-difference surface motion", test_geo07},
    {"GEO08", "3D", "Rotational covariance of geometry and local scalar physics", test_geo08},
    {"MSH01", "3D", "Shock-mesh nondegeneracy", test_msh01},
    {"MSH02", "3D", "Mesh-cell orientation relative to analytical shock normal", test_msh02},
    {"MSH03", "3D", "Shock-surface area convergence under angular refinement", test_msh03},
    {"MSH04", "3D", "Unique apex and periodic-seam topology", test_msh04},
    {"MSH05", "3D", "Area-weighted stochastic source-patch sampling", test_msh05},
    {"KIN01", "COMMON", "Fast-CME drag-based-model closed-form solution", test_kin01},
    {"KIN02", "COMMON", "Slow-CME sign-aware DBM branch", test_kin02},
    {"KIN03", "COMMON", "Zero-drag ballistic limit", test_kin03},
    {"KIN04", "COMMON", "Small-Gamma continuity across DBM/ballistic limit", test_kin04},
    {"KIN05", "COMMON", "Long-time DBM asymptotic stability", test_kin05},
    {"KIN06", "COMMON", "Data-driven PCHIP knot exactness", test_kin06},
    {"KIN07", "COMMON", "Data-driven monotonicity and no overshoot", test_kin07},
    {"KIN08", "COMMON", "Explicit data-driven time-domain policy", test_kin08},
    {"SHK01", "COMMON", "Fast-shock existence and no-shock threshold", test_shk01},
    {"SHK02", "COMMON", "Shock obliquity angle and polarity invariance", test_shk02},
    {"SHK03", "COMMON", "Parallel fast-shock limiting solution", test_shk03},
    {"SHK04", "COMMON", "Perpendicular MHD-shock benchmark", test_shk04},
    {"SHK05", "COMMON", "Oblique-MHD benchmark grid and branch continuity", test_shk05},
    {"SHK06", "COMMON", "Rankine-Hugoniot mass-flux conservation", test_shk06},
    {"SHK07", "COMMON", "Rankine-Hugoniot normal magnetic-field continuity", test_shk07},
    {"SHK08", "COMMON", "Rankine-Hugoniot tangential electric-field conservation", test_shk08},
    {"SHK09", "COMMON", "Rankine-Hugoniot momentum-flux conservation", test_shk09},
    {"SHK10", "COMMON", "Rankine-Hugoniot total-energy-flux conservation", test_shk10},
    {"SHK11", "COMMON", "Physical admissibility and entropy increase", test_shk11},
    {"SHK12", "COMMON", "Near-Mach-one weak-shock conditioning", test_shk12},
    {"SHK13", "3D", "Shock-state independence from arbitrary query radius", test_shk13},
    {"SHK14", "3D", "Canonical shock-state consistency across diagnostics and mesh", test_shk14},
    {"CON01", "3D", "Zero-solar-rotation radial connectivity limit", test_con01},
    {"CON02", "3D", "Parker-spiral field line intersecting a spherical shock", test_con02},
    {"CON03", "3D", "No connection to a finite-width shock", test_con03},
    {"CON04", "3D", "Tangent and near-tangent magnetic connection", test_con04},
    {"CON05", "3D", "Multiple intersections and deterministic root selection", test_con05},
    {"CON06", "3D", "Time-continuous cobpoint history", test_con06},
    {"CON07", "3D", "Cobpoint-to-ShockState consistency", test_con07},
    {"CON08", "3D", "Connectivity Parker path-length accuracy", test_con08},
    {"REG01", "COMMON", "SHOCK_ONLY upstream-field identity", test_reg01},
    {"REG02", "COMMON", "FULL_ICME resolved-shock inner Rankine-Hugoniot boundary", test_reg02},
    {"REG03", "1D", "Magnetic-ejecta density and velocity factors", test_reg03},
    {"REG04", "3D", "Self-similar nested shock/leading/trailing surfaces", test_reg04},
    {"REG05", "COMMON", "Continuity and smoothness across region transitions", test_reg05},
    {"ACC01", "COMMON", "SOURCE mode explicit source and SHOCK_ONLY transport flow", test_acc01},
    {"ACC02", "COMMON", "Acceleration-mode mutual-exclusion validation", test_acc02},
    {"ACC03", "COMMON", "Resolved-compression C1 shock smoothing and RH endpoints", test_acc03},
    {"ACC04", "1D<->3D", "Resolved shock-profile identity in spherical radial limit", test_acc04},
    {"ACC05", "COMMON", "Resolved mode disables prescribed DSA source", test_acc05},
    {"ERR01", "COMMON", "1-D outside-domain status and no radius clipping", test_err01},
    {"ERR02", "3D", "Cartesian non-finite input propagation", test_err02},
    {"ERR03", "COMMON", "Explicit Rankine-Hugoniot solver outcome status", test_err03},
    {"ERR04", "3D", "Degenerate direction rejection without +X fallback", test_err04},
    {"ERR05", "3D", "Tecplot writer rejects non-finite physics data", test_err05},
    {"DIV01", "COMMON", "Analytical divergence of constant radial solar wind", test_div01},
    {"DIV02", "COMMON", "Manufactured radial-flow divergence", test_div02},
    {"DIV03", "3D", "General Cartesian divergence convergence", test_div03},
    {"1D_AMBIENT_01", "1D", "Ambient values at 1 AU",
     test_1d_ambient_at_one_au},
    {"1D3D01", "1D<->3D", "Common-core upstream solar-wind state identity", test_1d3d01},
    {"1D3D02", "1D<->3D", "Common-core 1-D versus 3-D shock-state identity", test_1d3d02},
    {"1D3D03", "1D<->3D", "SOURCE acceleration-record identity", test_1d3d03},
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

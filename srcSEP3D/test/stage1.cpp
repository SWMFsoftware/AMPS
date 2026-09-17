// ============================================================================
// srcSEP3D/test/stage1.cpp
//
// AMPS-independent standalone test runner.
//
// LAYER: test binary only — not part of any library.
//   Links: CORE_OBJ + BG_OBJ + RUNTIME_OBJ + test objects.
//   Must NOT link any AMPS or MPI object.  The linker enforcing this is the
//   whole point of having a separate binary.
//
// WHAT "STANDALONE" MEANS:
//   These tests exercise core/, mesh/, background/, turbulence/, and runtime/
//   without any AMPS or MPI dependency.  The complete Stage-1 suite runs in
//   under two minutes on a laptop.  This is what makes the suite useful
//   during development rather than only in CI — a developer can run it
//   between every commit, without needing an AMPS installation.
//
//   R0/R1 source/ABI gates are orchestrated by test/run_tests.py because
//   they inspect the production tree or an external AMPS checkout.  Later
//   linked and validation tests must run through the production executable.
//
// HOW THE BINARY ENFORCES THE LAYERING BOUNDARY:
//   The linker sees AMPS-independent model/runtime and test objects — the Stage-1 link line in
//   the makefile does not include any AMPS archive.  If a developer adds
//   #include "pic.h" to a core/ file, the linker will produce an undefined-
//   symbol error here, immediately.  That error is BLD01 (test group BLD)
//   checked in positive form — see test/individual-test/test_build.cpp.
//
// FRAMEWORK:
//   Uses SEP::Testing::Registry from sep_common/sep_test_registry.{h,cpp},
//   accessed through the srcSEP3D alias header core/sep3d_test_registry.h.
//   The alias keeps test code independent of the installed sibling layout.
//
// CLI INTERFACE:
//   ./test/stage1 --help
//   ./test/stage1 --list-tests
//   ./test/stage1 --list-groups
//   ./test/stage1 --all-tests
//   ./test/stage1 --test LAY01
//   ./test/stage1 --test-group HARN
//   ./test/stage1 --all-tests --test-json report.json --test-junit report.xml
//
// EXIT CODES (match srcSEP and the SEP::Testing::Summary::ExitCode() contract):
//   0  all selected tests passed or were skipped (Skip is not a failure)
//   1  at least one test FAILED
//   2  at least one test produced an ERROR, or a usage error occurred
//
// HOW TO ADD A NEW TEST GROUP (summary; full procedure in test/README.md):
//   1. Create test/individual-test/test_<group>.cpp
//   2. Implement the Register<Group>Tests() function declared in
//      core/sep3d_test_registry.h
//   3. Uncomment the forward declaration and the append() call below
//   4. Uncomment the make target in the makefile
//   5. Verify: make test-stage1 lists the new group
// ============================================================================

#include "sep3d_test_registry.h"   // core/sep3d_test_registry.h via -Icore

#include <iostream>
#include <map>
#include <string>
#include <vector>

// ============================================================================
// Registration order
//
// Groups are appended in the order below.  The Registry sorts descriptors by
// ID before any run, so the run order is always alphabetical by test ID
// regardless of the append order.  The append order controls only which
// group's descriptors appear first in the raw vector before sorting.
//
// Convention: put the harness self-tests (HARN) before physics tests so that
// a broken harness is reported first, with a clear diagnosis, rather than
// being buried under physics failures.
// ============================================================================

// ---- Step 2: active groups -------------------------------------------------
// All three are declared in core/sep3d_test_registry.h.
// Implementations are in test/individual-test/test_harness.cpp,
// test/individual-test/test_layering.cpp, and test/individual-test/test_build.cpp.

// ---- Later inactive groups (uncomment as steps are completed) --------------
// std::vector<SEP3D::Testing::Descriptor> RegisterTimestepTests();    // Step 20
// std::vector<SEP3D::Testing::Descriptor> RegisterParkerMoverTests(); // Step 21
// std::vector<SEP3D::Testing::Descriptor> RegisterFTEMoverTests();    // Step 22


// ============================================================================
// Print help text
// ============================================================================
static void PrintHelp(const char* argv0) {
  std::cout
    << "srcSEP3D AMPS-independent standalone test runner\n"
    << "Uses SEP::Testing::Registry (same harness as srcSEP)\n"
    << "\n"
    << "Usage:\n"
    << "  " << argv0 << " [options]\n"
    << "\n"
    << "Options:\n"
    << "  --help                    Print this message and exit\n"
    << "  --list-tests              List all registered tests\n"
    << "  --list-groups             List test groups with counts\n"
    << "  --all-tests               Run all routine tests\n"
    << "  --test ID                 Run the test with this ID (repeatable;\n"
    << "                              case-insensitive; unknown IDs are errors)\n"
    << "  --test-group GROUP        Run all tests in GROUP (repeatable;\n"
    << "                              case-insensitive; unknown groups are errors)\n"
    << "  --test-json PATH          Write JSON summary to PATH\n"
    << "  --test-junit PATH         Write JUnit XML summary to PATH\n"
    << "  --test-input PATH         Set input file path for end-to-end tests\n"
    << "  --artifact-directory DIR  Set artifact directory for end-to-end tests\n"
    << "\n"
    << "Exit codes:\n"
    << "  0  all selected tests passed or were skipped\n"
    << "  1  at least one FAILED\n"
    << "  2  at least one ERROR, or a usage error (unknown option, bad ID)\n"
    << "\n"
    << "Groups through Phase O:   ADP3D  BGP3D  BLD  COEF3D  FTE3D  HARN  LAY  LIFE3D  MSH3D  NAT3D  PRK3D  RNG3D  RST3D  SHK3D  SNAP3D  TUR3D  UTIL\n"
    << "Frozen records:    test/frozen/\n"
    << "Test artifacts:    test/individual-test/\n"
    << "Full procedure:    test/README.md\n";
}


// ============================================================================
// main
// ============================================================================
int main(int argc, char** argv) {

  // argv[0] is passed to tests that need to invoke the runner as a subprocess
  // (HARN04, BLD01).  If argc == 0, a defensive empty string is used.
  const std::string selfPath = (argc > 0) ? argv[0] : "";

  // ---- Collect descriptors from all active groups -------------------------
  //
  // Each Register*Tests() call returns a fully populated vector of Descriptors.
  // We concatenate them all, then pass the combined vector to the Registry
  // constructor which validates (no empty fields, no duplicate IDs), de-
  // duplicates, and sorts by ID.
  //
  // The Registry constructor throws std::invalid_argument on any violation.
  // That is a programming error, not a test failure; we let it propagate so
  // the developer sees the full message rather than a silent error return.
  std::vector<SEP3D::Testing::Descriptor> all;

  auto append = [&](std::vector<SEP3D::Testing::Descriptor> v) {
    for (auto& d : v) all.push_back(std::move(d));
  };

  // HARN must come first so a broken harness is visible before physics tests.
  append(RegisterHarnessTests());
  append(RegisterLayeringTests(selfPath));
  append(RegisterBuildTests(selfPath));
  append(RegisterKernelTests());   // UTIL — shared-kernel frozen record (Step 3)
  append(RegisterRuntimeTests());  // LIFE3D — Phase R2 lifecycle
  append(RegisterMeshTests());     // MSH3D — Phase M mesh/storage
  append(RegisterBackgroundTests());  // BGP3D/SNAP3D — Phase B
  append(RegisterTurbulenceTests());  // TUR3D/COEF3D — Phase T
  append(RegisterTransportTests());   // PRK3D/FTE3D/RNG3D — Phase P
  append(RegisterAdapterTests());     // ADP3D/NAT3D/SHK3D — Phase A
  append(RegisterOutputTests());      // NAT3D/RST3D — Phase O

  // Linked multi-rank validation groups are added only when they execute the
  // configured AMPS binary; they must not be represented by standalone mocks.

  SEP3D::Testing::Registry registry(std::move(all));


  // ---- Parse arguments ----------------------------------------------------
  //
  // Parsing happens AFTER the Registry is built so that --list-tests can show
  // the full list even if no run is requested.  This mirrors srcSEP's parser
  // contract: discovery options exit before initialisation, but the descriptor
  // vector is populated first so the list is complete.
  //
  // Malformed options exit with code 2, not 1, to distinguish a usage error
  // from a test failure.  The SEP::Testing::Summary::ExitCode() contract
  // uses the same three-value scheme.

  if (argc == 1) { PrintHelp(argv[0]); return 0; }

  bool runAll = false;
  std::vector<std::string> runIds, runGroups;
  std::string jsonPath, junitPath, inputPath, artifactDir;

  for (int i = 1; i < argc; ++i) {
    std::string a(argv[i]);

    if (a == "--help") {
      PrintHelp(argv[0]); return 0;
    }
    else if (a == "--all-tests") {
      runAll = true;
    }
    else if (a == "--list-tests") {
      // Print the full descriptor table and exit before any initialisation.
      registry.PrintList(std::cout);
      return 0;
    }
    else if (a == "--list-groups") {
      // Summarise by group; print sorted for reproducibility.
      std::map<std::string,int> counts;
      for (const auto& d : registry.Descriptors()) counts[d.group]++;
      for (const auto& kv : counts)
        std::cout << kv.first << "  (" << kv.second << " tests)\n";
      return 0;
    }
    else if (a == "--test" && i+1 < argc) {
      runIds.push_back(argv[++i]);
    }
    else if (a == "--test-group" && i+1 < argc) {
      runGroups.push_back(argv[++i]);
    }
    else if (a == "--test-json" && i+1 < argc) {
      jsonPath = argv[++i];
    }
    else if (a == "--test-junit" && i+1 < argc) {
      junitPath = argv[++i];
    }
    else if (a == "--test-input" && i+1 < argc) {
      inputPath = argv[++i];
    }
    else if (a == "--artifact-directory" && i+1 < argc) {
      artifactDir = argv[++i];
    }
    else {
      // Unknown option: exit 2 (usage error) with a clear message.
      // srcSEP uses the same policy: "fail before initialisation".
      std::cerr << "[stage1] Unknown option '" << a
                << "'.  Run with --help for usage.\n";
      return 2;
    }
  }


  // ---- Install execution context ------------------------------------------
  //
  // The execution context carries the input-file path and artifact directory
  // to end-to-end tests that need them.  It is set once, before any test
  // runs, and is read-only after that (preventing one test from redirecting
  // a later test's I/O).
  SEP3D::Testing::SetExecutionContext({inputPath, artifactDir});


  // ---- Select tests -------------------------------------------------------
  //
  // Registry::Select() validates every requested ID and group name before
  // selecting anything.  An unknown ID or group is a usage error (exit 2),
  // not a test failure (exit 1).  This prevents a misspelled test name from
  // silently running nothing and reporting success.
  std::vector<const SEP3D::Testing::Descriptor*> selected;
  try {
    selected = registry.Select(runIds, runGroups, runAll);
  }
  catch (const std::invalid_argument& ex) {
    std::cerr << "[stage1] Selection error: " << ex.what() << "\n";
    return 2;
  }

  if (selected.empty()) {
    // This can only happen if no options that trigger a run were given.
    std::cerr << "[stage1] No tests selected.  "
                 "Use --all-tests, --test ID, or --test-group GROUP.\n";
    return 2;
  }


  // ---- Run ----------------------------------------------------------------
  SEP3D::Testing::Summary summary = registry.Run(selected, std::cout);


  // ---- Write structured reports -------------------------------------------
  //
  // A write failure is reported to stderr and reflected in the exit code via
  // a synthetic Error result added to the summary.  The process must never
  // exit 0 after a report write fails, because the caller (CI, campaign
  // runner) depends on the report file being valid.
  if (!jsonPath.empty()) {
    std::string err;
    if (!SEP3D::Testing::WriteJsonSummary(summary, jsonPath, &err)) {
      std::cerr << "[stage1] JSON report error: " << err << "\n";
      // Force a non-zero exit even if all tests passed.
      return 2;
    }
  }
  if (!junitPath.empty()) {
    std::string err;
    if (!SEP3D::Testing::WriteJUnitSummary(summary, junitPath, &err)) {
      std::cerr << "[stage1] JUnit report error: " << err << "\n";
      return 2;
    }
  }

  return summary.ExitCode();
}

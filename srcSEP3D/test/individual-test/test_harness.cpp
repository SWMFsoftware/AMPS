// ============================================================================
// srcSEP3D/test/individual-test/test_harness.cpp
//
// Test group HARN — test-harness self-verification.
//
// WHY THIS GROUP EXISTS:
//   The four HARN tests verify properties of the runner itself before any
//   physics test runs.  A broken harness produces misleading results: a
//   runner that ignores failures would report PASS for every test; a runner
//   that treats Skip as Fail would break every future group that uses Skip
//   for "not yet implemented" guards.  Catching those defects first, with a
//   clear diagnosis, saves time.
//
//   These tests are intentionally simple.  They do not test physics; they
//   test the infrastructure everything else depends on.
//
// TESTS (match the plan document Section 12.3 Step 2 exactly):
//
//   HARN01 — Empty registry: lists nothing, exits 0.
//     The runner must treat an empty test set as a valid (if trivial) outcome,
//     not as an error.  This matters because 'make test-stage1' will be run
//     throughout development before all groups are populated; it must always
//     exit 0 for a clean tree.
//     HOW: invoke the stage1 binary with a non-existent group (no tests in
//     it) via --test-group EMPTY, which the registry will reject.  So instead
//     we build a fresh Registry with zero descriptors and assert it lists
//     nothing and its Select(…,…,false) on empty IDs returns an empty vector
//     without throwing.
//
//   HARN02 — A Fail result produces ExitCode()==1; an Error result produces
//     ExitCode()==2; a Pass result produces ExitCode()==0.
//     HOW: build throwaway one-descriptor registries whose callbacks return
//     Fail / Error / Pass, run each, and assert the exit code and counters.
//     Returns Pass on success (does NOT itself fail).
//
//   HARN03 — A Skip result produces ExitCode()==0, skipped==1, passed==0;
//     a mixed Pass+Skip run produces exit 0 with passed==1 and skipped==1.
//     HOW: build throwaway registries and assert the counters and exit code.
//     Returns Pass on success (does NOT itself skip).
//
//   HARN04 — JSON and JUnit writers produce parseable files; a write to an
//     unwritable path is reported as ERROR with a non-zero exit, not a warning.
//     HOW (happy path): run a Pass stub, write both report files, check for
//     schema markers.
//     HOW (failure path): attempt to write to /dev/null/impossible (a path
//     that cannot be created on any POSIX system) and verify WriteJsonSummary
//     returns false with a non-empty error string.
//
// FROZEN RECORD: none — no numerical output.
//
// DESIGN — negative controls that do not leave the suite red:
//   HARN02 and HARN03 are negative controls: they deliberately drive Fail,
//   Error, and Skip results through throwaway in-process registries and
//   assert the runner reacts correctly.  Crucially, they verify this
//   INTERNALLY and then return Pass, so a healthy --all-tests run exits 0.
//
//   The OUTER-process exit codes (does the whole binary exit 1 on a Fail-only
//   run, and 0 on a Skip-only run?) are verified independently by the shell
//   scripts run_harn02.sh and run_harn03.sh, which the Python runner reports
//   as HARN02-EXITCODE and HARN03-EXITCODE.  Between the in-process assertions
//   here and those two shell checks, the exit-code contract is fully covered
//   without any test being permanently red.
// ============================================================================

#include "sep3d_test_registry.h"   // core/sep3d_test_registry.h via -Icore

#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

// ============================================================================
// Convenience helpers that build a SEP3D::Testing::Result with minimal
// boilerplate.  Mirrors the idiom used in srcSEP test files.
// ============================================================================
using R = SEP3D::Testing::Result;
using S = SEP3D::Testing::Status;

R Pass(std::string msg = "ok") {
  R r; r.status = S::Pass; r.message = std::move(msg); return r;
}
R Fail(std::string msg) {
  R r; r.status = S::Fail; r.message = std::move(msg); return r;
}
R Skip(std::string msg) {
  R r; r.status = S::Skip; r.message = std::move(msg); return r;
}
R Err(std::string msg) {
  R r; r.status = S::Error; r.message = std::move(msg); return r;
}

// ============================================================================
// HARN01 — empty registry lists nothing and exits 0
//
// WHAT IS BEING TESTED:
//   Registry({}) with zero descriptors must:
//     (a) construct without throwing
//     (b) PrintList() produce no test lines
//     (c) Select({}, {}, false) return an empty vector, not throw
//     (d) Summary{}.ExitCode() return 0
//
// WHY THESE PROPERTIES MATTER:
//   During development, test/stage1.cpp will compile with future Register*()
//   calls still commented out.  'make test-stage1' runs --all-tests, which
//   on an empty registry should exit 0 cleanly.  If it exited non-zero
//   because the selection was empty, every developer would need to uncomment
//   a placeholder test just to get a green build, which defeats the purpose
//   of incremental development.
// ============================================================================
SEP3D::Testing::Result run_HARN01() {
  // Build a registry with zero descriptors.
  std::vector<SEP3D::Testing::Descriptor> empty;
  SEP3D::Testing::Registry r(std::move(empty));

  // (a) Construction succeeded — we reached this point.

  // (b) PrintList produces no test descriptor lines.
  //
  // Capture PrintList output first.
  std::ostringstream buf;
  r.PrintList(buf);
  const std::string listing = buf.str();
  //
  // The real PrintList always emits a two-line header:
  //   "Registered srcSEP standalone component tests\n"
  //   "ID | group | class | initialization | build modes | description\n"
  // Then for each descriptor it emits: "<ID> | <group> | ..."
  // followed by two indented lines for seed policy and state isolation.
  //
  // With an empty registry there are no descriptors, so nothing after the
  // two header lines.  We check that no line looks like a descriptor line,
  // i.e. no line starts with a non-whitespace character AND is not one of
  // the two known header lines.
  {
    std::istringstream lines(listing);
    std::string line;
    int lineNum = 0;
    const std::string header1 = "Registered srcSEP standalone component tests";
    const std::string header2 = "ID | group | class | initialization";
    while (std::getline(lines, line)) {
      ++lineNum;
      if (line.empty()) continue;
      // Header lines are expected and benign
      if (line.find(header1) != std::string::npos) continue;
      if (line.find(header2) != std::string::npos) continue;
      // Any other non-empty line is unexpected in an empty registry listing
      return Fail(
          "PrintList() on an empty registry produced an unexpected line "
          "(line " + std::to_string(lineNum) + "): '" + line + "'");
    }
  }

  // (c) Select on an empty registry with no selectors returns an empty vector.
  const std::vector<const SEP3D::Testing::Descriptor*> sel =
      r.Select({}, {}, false);
  if (!sel.empty()) {
    return Fail("Select({},{},false) on an empty registry returned "
                + std::to_string(sel.size()) + " descriptors, expected 0.");
  }

  // (d) An empty Summary exits 0.
  SEP3D::Testing::Summary s;
  if (s.ExitCode() != 0) {
    return Fail("Empty Summary::ExitCode() returned "
                + std::to_string(s.ExitCode()) + ", expected 0.");
  }

  return Pass("empty registry: PrintList has no test lines, "
              "Select returns empty vector, empty Summary exits 0.");
}


// ============================================================================
// HARN02 — a Fail result must produce ExitCode()==1; an Error result must
//          produce ExitCode()==2.  Verified entirely through a throwaway
//          in-process registry; this test returns Pass on success.
//
// WHAT IS BEING TESTED:
//   The runner's exit-code contract for the two non-success outcomes:
//     - a Status::Fail  result  =>  Summary::ExitCode() == 1, failed == 1
//     - a Status::Error result  =>  Summary::ExitCode() == 2, errors == 1
//   and that a Pass result on its own gives ExitCode() == 0.
//
// WHY THIS MATTERS:
//   CI scripts, campaign runners, and Make all rely on the exit code to know
//   whether tests passed.  A runner that exits 0 on failure would hide broken
//   tests from every automated system — the single most dangerous defect a
//   test harness can have.
//
// DESIGN NOTE — negative control without an always-red result:
//   This is a negative control: it deliberately drives a Fail (and an Error)
//   through a throwaway registry and asserts the runner reacts correctly.
//   The assertion is what makes it a control — if the runner ever stopped
//   distinguishing Fail from Pass, this test would itself fail.
//
//   Earlier this test also returned Fail itself, so that the OUTER runner
//   would demonstrate a non-zero exit live.  That made --all-tests always
//   exit 1, which is noisy.  The outer exit-code behaviour is instead
//   verified independently by test/individual-test/run_harn02.sh (which the
//   Python runner reports as HARN02-EXITCODE).  This test therefore returns
//   Pass on success, and the suite exits 0 when healthy.
// ============================================================================
SEP3D::Testing::Result run_HARN02() {

  // Helper: build a one-descriptor registry whose single test returns the
  // given status, run it, and return the resulting Summary.
  auto run_one = [](SEP3D::Testing::Status status,
                    const char* id) -> SEP3D::Testing::Summary {
    SEP3D::Testing::Descriptor d;
    d.id                  = id;
    d.name                = "stub for HARN02 exit-code check";
    d.group               = "HARN_INTERNAL";
    d.description         = "internal throwaway stub for HARN02";
    d.supportedBuildModes = "all";
    d.seedPolicy          = "deterministic-no-rng";
    d.stateIsolation      = "none";
    d.callback = [status]() -> SEP3D::Testing::Result {
      SEP3D::Testing::Result r;
      r.status  = status;
      r.message = "deliberate stub result for HARN02";
      return r;
    };
    SEP3D::Testing::Registry reg({d});
    auto selected = reg.Select({id}, {}, false);
    std::ostringstream devnull;
    return reg.Run(selected, devnull);
  };

  // ---- (1) A Fail result must give ExitCode()==1 and failed==1 ------------
  {
    auto s = run_one(SEP3D::Testing::Status::Fail, "HARN02_FAIL_STUB");
    if (s.ExitCode() != 1)
      return Fail("A Fail result gave ExitCode()=="
                  + std::to_string(s.ExitCode()) + "; expected 1.  "
                  "The runner is not reporting failures as a non-zero exit — "
                  "this is the most dangerous possible harness defect.");
    if (s.failed != 1)
      return Fail("A Fail result gave failed=="
                  + std::to_string(s.failed) + "; expected 1.");
    if (s.passed != 0 || s.skipped != 0 || s.errors != 0)
      return Fail("A Fail result contaminated the other counters "
                  "(passed/skipped/errors should all be 0).");
  }

  // ---- (2) An Error result must give ExitCode()==2 and errors==1 ----------
  {
    auto s = run_one(SEP3D::Testing::Status::Error, "HARN02_ERROR_STUB");
    if (s.ExitCode() != 2)
      return Fail("An Error result gave ExitCode()=="
                  + std::to_string(s.ExitCode()) + "; expected 2.  "
                  "Error must be distinguishable from Fail by exit code.");
    if (s.errors != 1)
      return Fail("An Error result gave errors=="
                  + std::to_string(s.errors) + "; expected 1.");
  }

  // ---- (3) A Pass result on its own must give ExitCode()==0 ---------------
  {
    auto s = run_one(SEP3D::Testing::Status::Pass, "HARN02_PASS_STUB");
    if (s.ExitCode() != 0)
      return Fail("A lone Pass result gave ExitCode()=="
                  + std::to_string(s.ExitCode()) + "; expected 0.");
    if (s.passed != 1)
      return Fail("A Pass result gave passed=="
                  + std::to_string(s.passed) + "; expected 1.");
  }

  return Pass(
      "exit-code contract verified: Fail->1 (failed==1), Error->2 (errors==1), "
      "Pass->0 (passed==1).  Outer-process exit code is checked separately by "
      "run_harn02.sh (reported as HARN02-EXITCODE).");
}


// ============================================================================
// HARN03 — a Skip result must give ExitCode()==0, skipped==1, passed==0.
//          Verified through a throwaway in-process registry; returns Pass.
//
// WHAT IS BEING TESTED:
//   When a test returns Status::Skip:
//     (a) ExitCode() must return 0 (Skip is not a failure)
//     (b) Summary::skipped must be 1
//     (c) Summary::passed must be 0 (Skip must NOT be relabelled Pass)
//   And, in a mixed run of one Pass + one Skip:
//     (d) ExitCode() == 0, passed == 1, skipped == 1, failed == 0
//         (a Skip alongside a Pass must neither fail the run nor be counted
//          as a second pass)
//
// WHY (c) AND (d) MATTER:
//   Later groups use Skip heavily for "prerequisite not available" — BLD01
//   skips when nm is absent; the external-field tests (BGX01-07, Step 28)
//   are registered now but skip until the launcher exists.  If Skip were
//   silently promoted to Pass, those would report a false green, which is
//   worse than not running them.  If Skip counted as Fail, they would turn
//   the build red for no reason and people would learn to ignore failures.
//
// DESIGN NOTE:
//   Like HARN02, this is a negative-control-style check that no longer
//   returns Skip itself.  The outer-process exit code for a Skip-only run is
//   verified independently by run_harn03.sh (reported as HARN03-EXITCODE).
//   This test returns Pass on success so the suite exits 0 when healthy.
// ============================================================================
SEP3D::Testing::Result run_HARN03() {

  // Helper: run a registry built from the given (id, status) pairs and
  // return the Summary.
  auto run_set =
      [](const std::vector<std::pair<const char*, SEP3D::Testing::Status>>& specs)
      -> SEP3D::Testing::Summary {
    std::vector<SEP3D::Testing::Descriptor> ds;
    std::vector<std::string> ids;
    for (const auto& sp : specs) {
      SEP3D::Testing::Descriptor d;
      d.id                  = sp.first;
      d.name                = "stub for HARN03 skip check";
      d.group               = "HARN_INTERNAL";
      d.description         = "internal throwaway stub for HARN03";
      d.supportedBuildModes = "all";
      d.seedPolicy          = "deterministic-no-rng";
      d.stateIsolation      = "none";
      const SEP3D::Testing::Status status = sp.second;
      d.callback = [status]() -> SEP3D::Testing::Result {
        SEP3D::Testing::Result r;
        r.status  = status;
        r.message = "deliberate stub result for HARN03";
        return r;
      };
      ds.push_back(std::move(d));
      ids.push_back(sp.first);
    }
    SEP3D::Testing::Registry reg(std::move(ds));
    auto selected = reg.Select(ids, {}, false);
    std::ostringstream devnull;
    return reg.Run(selected, devnull);
  };

  // ---- (a)-(c) A lone Skip result -----------------------------------------
  {
    auto s = run_set({{"HARN03_SKIP_STUB", SEP3D::Testing::Status::Skip}});
    if (s.ExitCode() != 0)
      return Fail("A Skip result gave ExitCode()=="
                  + std::to_string(s.ExitCode())
                  + "; expected 0.  Skip must not be treated as a failure.");
    if (s.skipped != 1)
      return Fail("A Skip result gave skipped=="
                  + std::to_string(s.skipped) + "; expected 1.");
    if (s.passed != 0)
      return Fail("A Skip result gave passed=="
                  + std::to_string(s.passed) + "; expected 0.  "
                  "Skip must NOT be relabelled as Pass.");
    if (s.failed != 0 || s.errors != 0)
      return Fail("A Skip result contaminated the failed/errors counters.");
  }

  // ---- (d) A mixed Pass + Skip run ----------------------------------------
  {
    auto s = run_set({
      {"HARN03_PASS_STUB", SEP3D::Testing::Status::Pass},
      {"HARN03_SKIP_STUB", SEP3D::Testing::Status::Skip},
    });
    if (s.ExitCode() != 0)
      return Fail("A Pass+Skip run gave ExitCode()=="
                  + std::to_string(s.ExitCode()) + "; expected 0.");
    if (s.passed != 1)
      return Fail("A Pass+Skip run gave passed=="
                  + std::to_string(s.passed) + "; expected 1.");
    if (s.skipped != 1)
      return Fail("A Pass+Skip run gave skipped=="
                  + std::to_string(s.skipped) + "; expected 1.");
    if (s.failed != 0 || s.errors != 0)
      return Fail("A Pass+Skip run produced spurious failed/errors counts.");
  }

  return Pass(
      "skip contract verified: lone Skip -> exit 0, skipped==1, passed==0; "
      "mixed Pass+Skip -> exit 0, passed==1, skipped==1.  Outer-process exit "
      "code is checked separately by run_harn03.sh (HARN03-EXITCODE).");
}


// ============================================================================
// HARN04 — JSON and JUnit writers produce parseable files; unwritable path
//          is reported as ERROR with non-zero exit, not as a warning.
//
// WHAT IS BEING TESTED (two sub-cases):
//
//   Happy path:
//     Build a one-descriptor registry, run it, call WriteJsonSummary and
//     WriteJUnitSummary.  The output files must exist and contain:
//       JSON  : the schema marker "srcsep-component-tests-v1", the test ID,
//               and the status string "PASS".
//       JUnit : the <testsuite> and <testcase> elements, and the test ID.
//
//   Failure path:
//     Call WriteJsonSummary with a path that cannot be created on any POSIX
//     system ("/dev/null/impossible").  It must return false with a non-empty
//     error string.  This verifies that the writer fails hard rather than
//     silently producing partial output and returning true.
//
// WHY THE FAILURE PATH MATTERS:
//   CI systems read the JSON report to determine test outcomes.  A writer
//   that returns true after failing would cause the CI system to read a
//   non-existent or incomplete report and produce a confusing result.
// ============================================================================
SEP3D::Testing::Result run_HARN04() {
  const std::string jsonPath  = "/tmp/sep3d_harn04_test.json";
  const std::string junitPath = "/tmp/sep3d_harn04_test.xml";
  const std::string badPath   = "/dev/null/impossible/file.json";

  // ---- Build a single-descriptor registry with a always-passing callback --
  SEP3D::Testing::Descriptor d;
  d.id                  = "HARN04_STUB";
  d.name                = "stub that returns Pass for HARN04 writer test";
  d.group               = "HARN_INTERNAL";
  d.description         = "internal stub for HARN04";
  d.supportedBuildModes = "all";
  d.seedPolicy          = "deterministic-no-rng";
  d.stateIsolation      = "none";
  d.callback = []() -> SEP3D::Testing::Result {
    SEP3D::Testing::Result r;
    r.status  = SEP3D::Testing::Status::Pass;
    r.message = "pass from HARN04 stub";
    return r;
  };

  SEP3D::Testing::Registry reg({d});
  auto selected = reg.Select({"HARN04_STUB"}, {}, false);
  std::ostringstream devnull;
  SEP3D::Testing::Summary summary = reg.Run(selected, devnull);

  // ---- Happy path: JSON ---------------------------------------------------
  {
    std::string err;
    if (!SEP3D::Testing::WriteJsonSummary(summary, jsonPath, &err))
      return Err("WriteJsonSummary happy path failed: " + err);

    std::ifstream f(jsonPath);
    if (!f) return Err("JSON file not created at " + jsonPath);

    std::string content((std::istreambuf_iterator<char>(f)), {});

    if (content.find("srcsep-component-tests-v1") == std::string::npos)
      return Fail("JSON missing schema field 'srcsep-component-tests-v1'");
    if (content.find("HARN04_STUB") == std::string::npos)
      return Fail("JSON does not contain the test ID 'HARN04_STUB'");
    if (content.find("\"PASS\"") == std::string::npos)
      return Fail("JSON does not contain status '\"PASS\"'");

    std::remove(jsonPath.c_str());
  }

  // ---- Happy path: JUnit --------------------------------------------------
  {
    std::string err;
    if (!SEP3D::Testing::WriteJUnitSummary(summary, junitPath, &err))
      return Err("WriteJUnitSummary happy path failed: " + err);

    std::ifstream f(junitPath);
    if (!f) return Err("JUnit file not created at " + junitPath);

    std::string content((std::istreambuf_iterator<char>(f)), {});

    if (content.find("<testsuite") == std::string::npos)
      return Fail("JUnit XML missing '<testsuite' element");
    if (content.find("<testcase") == std::string::npos)
      return Fail("JUnit XML missing '<testcase' element");
    if (content.find("HARN04_STUB") == std::string::npos)
      return Fail("JUnit XML does not contain the test ID 'HARN04_STUB'");

    std::remove(junitPath.c_str());
  }

  // ---- Failure path: unwritable path --------------------------------------
  //
  // /dev/null/impossible is not writable on any POSIX system because /dev/null
  // is a character device, not a directory, so creating a file inside it is
  // always ENOTDIR.
  {
    std::string err;
    const bool wrote = SEP3D::Testing::WriteJsonSummary(summary, badPath, &err);

    if (wrote) {
      // If it somehow wrote, something is very wrong.
      return Fail("WriteJsonSummary returned true for path '" + badPath
                  + "' which should be unwritable.  "
                  "The writer must fail hard on a bad path, not silently "
                  "succeed.");
    }
    if (err.empty()) {
      return Fail("WriteJsonSummary returned false but left the error string "
                  "empty for path '" + badPath + "'.  "
                  "The error string must be non-empty so the caller can "
                  "report the problem.");
    }
  }

  return Pass("JSON and JUnit happy-path files contain schema, ID, and status; "
              "unwritable path returns false with non-empty error string.");
}

} // anonymous namespace


// ============================================================================
// RegisterHarnessTests
//
// Returns the four HARN descriptors.  Called from test/stage1.cpp::main().
//
// METADATA NOTES:
//   initialization = None: no model state is needed; these tests run against
//   the harness infrastructure only.
//
//   runtime = Routine: these tests are fast (< 10 ms each) and must always
//   run with --all-tests.  They must not be Extended.
//
//   seedPolicy = "deterministic-no-rng": none of these tests use random
//   numbers, so there is no seed to report.
//
//   stateIsolation = "independent": each test builds its own Registry from
//   scratch; no shared state can leak between them.
// ============================================================================
// ============================================================================
// Exit-code beacons for the shell scripts.
//
// The two shell scripts run_harn02.sh and run_harn03.sh verify the runner's
// OUTER-process exit code — something no in-process test can do about itself.
// To do that they need a test that genuinely returns Fail (so the process
// exits 1) and one that genuinely returns Skip (so the process exits 0).
//
// These beacons provide exactly that.  They are registered as EXTENDED tests,
// which means:
//   - they are NOT included in --all-tests (which runs only Routine tests),
//     so they never make a normal run red; and
//   - they ARE individually selectable by ID, so the shell scripts can invoke
//     them with --test HARN_FAIL_BEACON / --test HARN_SKIP_BEACON.
//
// This is the mechanism that lets the exit-code contract be fully verified
// (in-process by HARN02/HARN03, out-of-process by the shell scripts) without
// any test in the routine suite being permanently red.
// ============================================================================
SEP3D::Testing::Result run_HARN_FAIL_BEACON() {
  return Fail(
      "exit-code beacon: this EXTENDED test always returns Fail so that "
      "run_harn02.sh can confirm the outer process exits 1.  It is not part "
      "of --all-tests.");
}

SEP3D::Testing::Result run_HARN_SKIP_BEACON() {
  return Skip(
      "exit-code beacon: this EXTENDED test always returns Skip so that "
      "run_harn03.sh can confirm the outer process exits 0.  It is not part "
      "of --all-tests.");
}


// ============================================================================
// RegisterHarnessTests
// ============================================================================
std::vector<SEP3D::Testing::Descriptor> RegisterHarnessTests() {
  using D  = SEP3D::Testing::Descriptor;
  using IL = SEP3D::Testing::InitializationLevel;
  using RC = SEP3D::Testing::RuntimeClass;

  // Shared metadata for all HARN tests
  const char* builds   = "all";
  const char* seeds    = "deterministic-no-rng";
  const char* isolation= "independent: each builds its own Registry from scratch";

  auto make = [&](const char* id, const char* name, const char* desc,
                  SEP3D::Testing::TestCallback cb) {
    D d;
    d.id                  = id;
    d.name                = name;
    d.group               = "HARN";
    d.description         = desc;
    d.initialization      = IL::None;
    d.supportedBuildModes = builds;
    d.runtime             = RC::Routine;
    d.seedPolicy          = seeds;
    d.stateIsolation      = isolation;
    d.callback            = std::move(cb);
    return d;
  };

  // Like make(), but marks the test EXTENDED and puts it in its own BEACON
  // group.  EXTENDED excludes it from --all-tests; the separate group keeps
  // it out of --group HARN too, so the only way to run a beacon is by its
  // explicit ID (which is exactly how the shell scripts invoke it).
  auto make_beacon = [&](const char* id, const char* name, const char* desc,
                         SEP3D::Testing::TestCallback cb) {
    D d = make(id, name, desc, std::move(cb));
    d.runtime = RC::Extended;
    d.group   = "BEACON";
    return d;
  };

  return {
    make("HARN01",
         "Empty registry: lists nothing and exits 0",
         "Constructs Registry({}), asserts PrintList has no test lines, "
         "Select returns empty vector, empty Summary exits 0.",
         run_HARN01),

    make("HARN02",
         "Exit-code contract: Fail->1, Error->2, Pass->0",
         "Drives Fail, Error, and Pass results through throwaway registries "
         "and asserts ExitCode() and counters; returns Pass on success.",
         run_HARN02),

    make("HARN03",
         "Skip contract: Skip->exit 0, skipped==1, passed==0",
         "Asserts a lone Skip and a mixed Pass+Skip run give the correct exit "
         "code and counters; returns Pass on success.",
         run_HARN03),

    make("HARN04",
         "JSON and JUnit writers produce parseable files; unwritable path fails hard",
         "Happy path: files contain schema marker, test ID, status.  "
         "Failure path: /dev/null/impossible returns false with non-empty error.",
         run_HARN04),

    // Exit-code beacons — EXTENDED, in their own BEACON group, excluded from
    // both --all-tests and --group HARN.  Used only by the shell scripts
    // run_harn02.sh / run_harn03.sh, which select them by explicit ID.
    make_beacon("HARN_FAIL_BEACON",
         "Exit-code beacon: always Fail (for run_harn02.sh)",
         "EXTENDED. Always returns Fail so the shell script can confirm the "
         "outer process exits 1.  Not part of --all-tests or --group HARN.",
         run_HARN_FAIL_BEACON),

    make_beacon("HARN_SKIP_BEACON",
         "Exit-code beacon: always Skip (for run_harn03.sh)",
         "EXTENDED. Always returns Skip so the shell script can confirm the "
         "outer process exits 0.  Not part of --all-tests or --group HARN.",
         run_HARN_SKIP_BEACON),
  };
}

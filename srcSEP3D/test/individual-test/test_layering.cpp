// ============================================================================
// srcSEP3D/test/individual-test/test_layering.cpp
//
// Test group LAY — layering-boundary tests.
//
// WHY THIS GROUP EXISTS:
//   The Stage-1 test strategy depends entirely on the guarantee that every
//   model directory outside amps/ contains no AMPS or MPI symbols. If that boundary erodes,
//   the Stage-1 binary silently starts requiring AMPS headers and the "test
//   without building AMPS" property is lost — without any obvious error
//   message.  These tests enforce the boundary actively, at every build,
//   rather than relying on code review.
//
// TESTS (match plan document Section 12.3 Step 2 exactly):
//
//   LAY01 — Positive case: all AMPS-independent directories are free of pic.h, mpi.h,
//            and PIC:: on the current working tree.
//
//            HOW: invokes grep -rn -E over all three directories for each
//            forbidden pattern.  If any match is found the test fails and
//            prints the offending lines.
//
//            WHAT IT CATCHES: a developer inadvertently adds #include "pic.h"
//            to a core file.  The Stage-1 linker would also catch it, but
//            the grep catches it earlier (without even running a build) and
//            produces a more targeted message.
//
//   LAY02 — Negative control: a scratch file with '#include "pic.h"' planted
//            in core/ IS detected by the same grep, then cleaned up
//            unconditionally.
//
//            HOW: writes a temp file to core/, runs the grep, verifies a
//            match is found, removes the file regardless of the outcome.
//
//            WHY THIS TEST MUST EXIST:
//            A test that cannot be made to fail is not evidence.  LAY02
//            confirms that LAY01's grep command is actually pointing at the
//            right directory, using the right patterns, and would catch a
//            real violation.  Without LAY02, a misconfigured grep (wrong
//            directory, wrong pattern) would let LAY01 pass even when
//            violations exist.
//
//            NEGATIVE CONTROL CONTRACT:
//            LAY02 is the canonical example of a negative control.  If LAY02
//            ever starts PASSING when it should FAIL (i.e. the planted
//            violation is NOT detected), that is a regression — the same as
//            LAY01 incorrectly failing.  The runner treats both equally.
//
// FROZEN RECORD: none — no numerical output.
//
// REFERENCE SOLUTIONS: none — pass/fail is determined by grep match presence.
//
// STATE ISOLATION:
//   LAY01 reads the filesystem; it does not write anything.
//   LAY02 writes and removes a single temp file in core/.  The file is
//   removed in the cleanup block that runs regardless of outcome; if the
//   remove fails, the test returns Error (not Fail) and names the dirty file.
// ============================================================================

#include "sep3d_test_registry.h"   // core/sep3d_test_registry.h via -Icore

#include <cstdio>
#include <string>
#include <vector>

namespace {

// ============================================================================
// Convenience helpers (same as test_harness.cpp — will move to a shared
// utility header at Step 3 when sep_common.a is introduced).
// ============================================================================
using R = SEP3D::Testing::Result;
using S = SEP3D::Testing::Status;

R Pass(std::string msg) {
  R r; r.status = S::Pass; r.message = std::move(msg); return r;
}
R Fail(std::string msg) {
  R r; r.status = S::Fail; r.message = std::move(msg); return r;
}
R Err(std::string msg) {
  R r; r.status = S::Error; r.message = std::move(msg); return r;
}

// ============================================================================
// GrepFound — run a shell command, return true iff exit code is 0.
//
// grep exits 0 if at least one match was found, 1 if no match, 2 on error.
// We capture and discard stdout so test output is not cluttered; the caller
// knows the patterns and can reconstruct the meaning from the return value.
// ============================================================================
bool GrepFound(const std::string& cmd) {
  FILE* p = popen(cmd.c_str(), "r");
  if (!p) return false;
  char buf[4096];
  while (fgets(buf, sizeof(buf), p)) { /* discard */ }
  return (pclose(p) == 0);  // 0 = match found
}

// Patterns that must not appear in L0/L1 source files.
// These are POSIX extended-regexp strings safe for grep -E.
// The test strips false matches that come from comment text that describes
// what is forbidden — e.g. this very file contains the string "pic.h" in
// a comment.  The -E version with word-boundary anchors is more precise
// than a plain substring search.
struct Pattern {
  const char* re;
  const char* description;
};

const Pattern FORBIDDEN[] = {
  { "#include.*pic\\.h",   "#include of pic.h" },
  { "#include.*mpi\\.h",   "#include of mpi.h" },
  { "\\bPIC::",            "use of PIC:: namespace" },
};
const int N_FORBIDDEN = 3;

const char* DIRS[] = {
    "core", "background", "runtime", "mesh", "turbulence", "transport",
    "adapters", "output", "validation"};
const int N_DIRS = 9;

// ============================================================================
// LAY01 — positive case: no forbidden symbols in core/ or background/
// ============================================================================
SEP3D::Testing::Result run_LAY01() {
  std::string violations;

  for (int di = 0; di < N_DIRS; ++di) {
    const char* dir = DIRS[di];
    for (int pi = 0; pi < N_FORBIDDEN; ++pi) {
      const Pattern& pat = FORBIDDEN[pi];

      // Build grep command:
      //   -r  recursive
      //   -n  print line numbers
      //   -E  extended regexp
      //   --include="*.h" --include="*.cpp"  source files only
      std::string cmd =
          std::string("grep -rn -E --include=\"*.h\" --include=\"*.cpp\" '")
          + pat.re + "' " + dir + "/ 2>/dev/null";

      if (GrepFound(cmd)) {
        violations += std::string("  [") + dir + "/] " + pat.description + "\n";
      }
    }
  }

  if (!violations.empty()) {
    return Fail(
        "Layering violations found in an AMPS-independent source directory:\n" + violations
        + "These files must not include pic.h, mpi.h, or use PIC::.\n"
        + "See README.md (Layering rules) and MIGRATION_MANIFEST.md.");
  }

  return Pass(
      "core/, background/, runtime/, mesh/, turbulence/, transport/, adapters/, output/, and validation/ contain no #include pic.h, #include mpi.h, "
      "or PIC:: references.");
}


// ============================================================================
// LAY02 — negative control: a deliberate violation is detected
// ============================================================================
SEP3D::Testing::Result run_LAY02() {
  // The temp file must have a name that is:
  //   (a) impossible to confuse with a real source file
  //   (b) visibly abnormal if it survives a failed cleanup
  const std::string tmpPath = "core/__LAY02_violation_temp_DELETE_IF_SEEN.h";

  // ---- Write the deliberately violating file --------------------------------
  {
    FILE* f = std::fopen(tmpPath.c_str(), "w");
    if (!f) {
      return Err("Cannot create temp file " + tmpPath + " — check that "
                 "core/ is writable.  LAY02 requires write access to plant "
                 "a controlled violation.");
    }
    std::fprintf(f,
        "// TEMPORARY FILE — created by LAY02 — must be deleted automatically\n"
        "// If you see this file, the LAY02 cleanup step failed.  Delete it.\n"
        "#include \"pic.h\"   // deliberate L0 violation for LAY02\n");
    std::fclose(f);
  }

  // ---- Run the same grep as LAY01, on the same directory -------------------
  const std::string cmd =
      "grep -rn -E --include=\"*.h\" '#include.*pic\\.h' core/ 2>/dev/null";
  const bool detected = GrepFound(cmd);

  // ---- Clean up unconditionally -------------------------------------------
  //
  // The remove must happen before any return so the tree is never left dirty.
  // If the remove fails we return Error (not Fail) because a dirty tree
  // interferes with LAY01 on the next run.
  const bool cleaned = (std::remove(tmpPath.c_str()) == 0);

  if (!cleaned) {
    return Err("Could not remove temp file " + tmpPath
               + " — tree is dirty.  Delete this file manually before "
               "running LAY01 again.");
  }

  // ---- Evaluate result ----------------------------------------------------
  if (!detected) {
    return Fail(
        "The grep used in LAY01 did NOT detect a deliberate "
        "#include \"pic.h\" planted in " + tmpPath + ".  "
        "This means LAY01's grep is not examining core/ correctly — "
        "check the working directory when make test-stage1 is run "
        "(it must be the srcSEP3D/ root, not a subdirectory).");
  }

  return Pass(
      "Deliberate #include \"pic.h\" in core/ was detected by the LAY01 "
      "grep, and the temp file was removed cleanly.");
}

} // anonymous namespace


// ============================================================================
// RegisterLayeringTests
//
// Returns the LAY descriptors.  The selfPath argument is passed by
// test/stage1.cpp but is not used by this group (it is needed only by BLD01
// which invokes nm on the binary).  The parameter is kept for API symmetry
// with RegisterBuildTests.
// ============================================================================
std::vector<SEP3D::Testing::Descriptor> RegisterLayeringTests(
    const std::string& /*selfPath*/) {

  using D  = SEP3D::Testing::Descriptor;
  using IL = SEP3D::Testing::InitializationLevel;
  using RC = SEP3D::Testing::RuntimeClass;

  auto make = [](const char* id, const char* name, const char* desc,
                 SEP3D::Testing::TestCallback cb) {
    D d;
    d.id                  = id;
    d.name                = name;
    d.group               = "LAY";
    d.description         = desc;
    d.initialization      = IL::None;
    d.supportedBuildModes = "all";
    d.runtime             = RC::Routine;
    d.seedPolicy          = "deterministic-no-rng";
    d.stateIsolation      =
        "LAY01: read-only filesystem scan.  "
        "LAY02: writes then removes one temp file in core/; "
        "Error if remove fails.";
    d.callback            = std::move(cb);
    return d;
  };

  return {
    make("LAY01",
         "AMPS-independent directories contain no AMPS symbols",
         "grep -rn -E for #include pic.h, #include mpi.h, PIC:: in "
         "core/, background/, runtime/, mesh/, turbulence/, transport/, adapters/, output/, and validation/; fails if any match is found.",
         run_LAY01),

    make("LAY02",
         "Negative control: planted pic.h include is detected (and removed)",
         "Creates core/__LAY02_violation_temp_DELETE_IF_SEEN.h with "
         "#include pic.h, verifies detection, removes it unconditionally.",
         run_LAY02),
  };
}

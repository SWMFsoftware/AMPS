// ============================================================================
// srcSEP3D/test/individual-test/test_build.cpp
//
// Test group BLD — build-structure tests.
//
// WHY THIS GROUP EXISTS:
//   LAY01 checks the source files.  BLD01 checks the binary.  These are
//   different things: a file could pass the LAY01 grep (no forbidden text
//   in its own source) while still pulling in a forbidden symbol through a
//   chain of transitive includes that grep does not follow.
//
//   BLD01 closes that gap by inspecting the symbol table of the compiled,
//   linked Stage-1 binary with nm.  If any AMPS or MPI symbol appears there,
//   the AMPS-independent model/runtime boundary has been breached at the binary level, regardless of
//   what the source says.
//
// TESTS:
//
//   BLD01 — The Stage-1 binary's symbol table contains no AMPS or MPI symbols.
//
//            HOW: runs  nm <binary>  (or  nm -u <binary>  to list undefined
//            references only) and checks for the absence of patterns that
//            identify AMPS and MPI symbols.
//
//            WHAT IT CATCHES:
//              - A core/ file that includes a header which (after preprocessing)
//                declares a symbol from AMPS — the linker would then require
//                the AMPS archive and would fail, but nm reports it without
//                needing to re-link.
//              - A transitive include chain:
//                  core/foo.h -> some_util.h -> pic.h
//                that grep misses because only core/foo.h is in core/.
//
//            SKIP CONDITION:
//              nm is not available on every system (e.g. some cross-
//              compilation environments).  If which nm returns non-zero
//              the test returns Skip rather than Error, because the absence
//              of nm is an environment limitation, not a code defect.
//
//            BINARY PATH:
//              Provided by test/stage1.cpp via argv[0].  If the path is empty
//              or points to a file that does not exist, the test returns Skip
//              with an explanation.
//
// FROZEN RECORD: none.
// REFERENCE SOLUTIONS: none — pass/fail from nm output.
// ============================================================================

#include "sep3d_test_registry.h"   // core/sep3d_test_registry.h via -Icore

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sys/stat.h>
#include <vector>

namespace {

using R = SEP3D::Testing::Result;
using S = SEP3D::Testing::Status;

R Pass(std::string msg) {
  R r; r.status = S::Pass; r.message = std::move(msg); return r;
}
R Fail(std::string msg) {
  R r; r.status = S::Fail; r.message = std::move(msg); return r;
}
R Skip(std::string msg) {
  R r; r.status = S::Skip; r.message = std::move(msg); return r;
}

// ============================================================================
// Helper: check whether a command is available.
// Returns true if 'which <cmd>' exits 0.
// ============================================================================
bool CommandAvailable(const char* cmd) {
  std::string check = std::string("which ") + cmd + " >/dev/null 2>&1";
  return (std::system(check.c_str()) == 0);
}

// ============================================================================
// Helper: run a command and return all output as a string.
// Returns empty string if the command fails or produces no output.
// ============================================================================
std::string RunCapture(const std::string& cmd) {
  FILE* p = popen(cmd.c_str(), "r");
  if (!p) return {};
  std::string out;
  char buf[4096];
  while (fgets(buf, sizeof(buf), p)) out += buf;
  pclose(p);
  return out;
}

// Patterns in nm output that indicate an AMPS or MPI symbol.
// These are substring matches (not regexps) for portability — nm output
// format varies between platforms but AMPS and MPI function names are
// distinctive enough to identify without a regexp engine.
struct BinaryPattern {
  const char* substring;
  const char* description;
};

const BinaryPattern FORBIDDEN_BINARY[] = {
  // MPI runtime functions
  { "MPI_Init",        "MPI initialisation symbol MPI_Init" },
  { "MPI_Comm_rank",   "MPI symbol MPI_Comm_rank" },
  { "MPI_Barrier",     "MPI symbol MPI_Barrier" },
  { "MPI_Allreduce",   "MPI symbol MPI_Allreduce" },
  // AMPS mesh / PIC symbols
  { "cTreeNodeAMR",    "AMPS AMR tree symbol cTreeNodeAMR" },
  { "PIC_Mesh",        "AMPS PIC mesh symbol PIC_Mesh" },
  { "PIC_ParticleBuffer", "AMPS particle buffer symbol" },
  { "amps_init",       "AMPS application entry point amps_init" },
};
const int N_FORBIDDEN_BINARY = 8;

// ============================================================================
// BLD01 — Stage-1 binary symbol table contains no AMPS or MPI symbols
// ============================================================================
SEP3D::Testing::Result run_BLD01(const std::string& binaryPath) {

  // ---- Skip if no binary path was provided --------------------------------
  if (binaryPath.empty()) {
    return Skip(
        "No binary path provided.  Run via 'make test-stage1' or pass argv[0] "
        "to RegisterBuildTests().  BLD01 requires the path to the stage1 "
        "binary to inspect its symbol table with nm.");
  }

  // ---- Skip if the binary does not exist ----------------------------------
  {
    struct ::stat st;
    if (::stat(binaryPath.c_str(), &st) != 0) {
      return Skip(
          "Binary not found at '" + binaryPath + "'.  "
          "Build first with 'make test-stage1'.");
    }
  }

  // ---- Skip if nm is not available ----------------------------------------
  if (!CommandAvailable("nm")) {
    return Skip(
        "nm is not available in PATH.  BLD01 requires nm to inspect the "
        "binary symbol table.  On Debian/Ubuntu: sudo apt install binutils.");
  }

  // ---- Run nm on the binary -----------------------------------------------
  //
  // nm -u lists only undefined (external) symbols — exactly what we care
  // about: symbols the binary needs from the linker but did not define
  // itself.  AMPS and MPI symbols would appear here if any core/ or
  // background/ file transitively included their declarations.
  //
  // On macOS nm -u works the same way; on ELF systems nm -u is also
  // standard.  We do not use --demangle because C linkage symbols (MPI_*,
  // amps_*) are not mangled anyway, and demangled AMPS C++ names vary by
  // ABI.
  const std::string output = RunCapture("nm -u '" + binaryPath + "' 2>&1");

  if (output.empty()) {
    // nm ran but produced nothing; this can happen on a fully static binary
    // or when the binary has been stripped.  Skip rather than falsely pass.
    return Skip(
        "nm -u produced no output for '" + binaryPath + "'.  "
        "The binary may be stripped or fully static.  "
        "BLD01 cannot verify the symbol table in this case.");
  }

  // ---- Check for forbidden symbols ----------------------------------------
  std::string violations;
  for (int i = 0; i < N_FORBIDDEN_BINARY; ++i) {
    const BinaryPattern& pat = FORBIDDEN_BINARY[i];
    if (output.find(pat.substring) != std::string::npos) {
      violations += std::string("  ") + pat.description + "\n";
    }
  }

  if (!violations.empty()) {
    return Fail(
        "Stage-1 binary '" + binaryPath + "' contains AMPS or MPI symbols:\n"
        + violations
        + "This means core/ or background/ has a transitive dependency on "
        "AMPS or MPI headers.  Run  make check-layering  to find the source "
        "file; also inspect the transitive include chain from that file.\n"
        "nm output excerpt (first 1000 chars):\n"
        + output.substr(0, 1000));
  }

  return Pass(
      "nm -u '" + binaryPath + "' contains no AMPS or MPI symbols.  "
      "The core/mesh/background/turbulence/runtime boundary is intact at the binary level.");
}

} // anonymous namespace


// ============================================================================
// RegisterBuildTests
// ============================================================================
std::vector<SEP3D::Testing::Descriptor> RegisterBuildTests(
    const std::string& binaryPath) {

  using D  = SEP3D::Testing::Descriptor;
  using IL = SEP3D::Testing::InitializationLevel;
  using RC = SEP3D::Testing::RuntimeClass;

  D d;
  d.id                  = "BLD01";
  d.name                = "Stage-1 binary has no AMPS or MPI symbols";
  d.group               = "BLD";
  d.description         =
      "Runs 'nm -u <binary>' and verifies no AMPS or MPI symbol names appear.  "
      "Skips if nm is unavailable or the binary path is unknown.";
  d.initialization      = IL::None;
  d.supportedBuildModes = "all";
  d.runtime             = RC::Routine;
  d.seedPolicy          = "deterministic-no-rng";
  d.stateIsolation      = "read-only: runs nm on the binary, no filesystem writes";
  d.callback            = [binaryPath]() -> SEP3D::Testing::Result {
    return run_BLD01(binaryPath);
  };

  return { d };
}

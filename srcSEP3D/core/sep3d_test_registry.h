// ============================================================================
// srcSEP3D/core/sep3d_test_registry.h
//
// srcSEP3D wrapper over SEP::Testing::Registry.
//
// LAYER: L0 (numerical core).
//   This header must NEVER include pic.h, mpi.h, or any AMPS symbol.
//   The make check-layering target enforces this at every build.
//
// PURPOSE:
//   Re-exports the SEP::Testing types under the SEP3D::Testing alias and
//   declares the per-group registration functions whose implementations live
//   in test/individual-test/.
//
//   Having a single srcSEP3D-owned header for all of this means:
//     - test/stage1.cpp includes exactly ONE header to get everything it needs
//     - Each test/individual-test/test_*.cpp includes ONE header and knows its
//       own types
//     - sep_common owns the implementation; every test file is insulated from
//       the surrounding AMPS checkout layout
//
// RELATIONSHIP TO srcSEP:
//   SEP::Testing::Registry is the test harness both applications consume from
//   sep_common/sep_test_registry.{h,cpp}. srcSEP3D aliases that public API
//   rather than maintaining a private registry.
//
//   Using the same Registry ensures that:
//     (a) Descriptor metadata requirements are identical in both models
//     (b) The JSON/JUnit report schema is identical in both models
//     (c) srcSEP's test tooling (campaign runner, CI scripts) can consume
//         srcSEP3D Stage-1 reports without modification
//
// PHASE R0 STATE:
//   HARN, LAY, BLD, and UTIL are active.  The header is resolved from the
//   shared sep_common installation/source directory supplied by the build; this
//   application does not compile a private registry implementation.
// ============================================================================

#ifndef SEP3D_CORE_SEP3D_TEST_REGISTRY_H
#define SEP3D_CORE_SEP3D_TEST_REGISTRY_H

// ---- Pull in the real SEP::Testing types -----------------------------------
//
// The runner/makefile adds the selected shared-header directory to the include
// path, so the source remains independent of the surrounding checkout layout.
#include "sep_test_registry.h"   // found via the canonical sep_common include path

#include <string>
#include <vector>

// ============================================================================
// SEP3D::Testing namespace
//
// A thin alias of SEP::Testing so that srcSEP3D code can write
// SEP3D::Testing::Descriptor etc. and the Step 3 change is completely
// transparent.
// ============================================================================
namespace SEP3D {
namespace Testing {

  // Re-export every public name from SEP::Testing.
  // The common library owns the definitions; these aliases keep the 3-D model
  // namespace readable without adding wrapper implementations.
  using Descriptor          = SEP::Testing::Descriptor;
  using Result              = SEP::Testing::Result;
  using Status              = SEP::Testing::Status;
  using Metric              = SEP::Testing::Metric;
  using Summary             = SEP::Testing::Summary;
  using Registry            = SEP::Testing::Registry;
  using ExecutionContext     = SEP::Testing::ExecutionContext;
  using InitializationLevel = SEP::Testing::InitializationLevel;
  using RuntimeClass        = SEP::Testing::RuntimeClass;
  using TestCallback        = SEP::Testing::TestCallback;
  using RefinementOrderEstimate = SEP::Testing::RefinementOrderEstimate;

  // Free functions — re-exported by bringing them into this namespace so
  // callers can write SEP3D::Testing::WriteJsonSummary(…) throughout.
  using SEP::Testing::SetExecutionContext;
  using SEP::Testing::GetExecutionContext;
  using SEP::Testing::StatusName;
  using SEP::Testing::PrintResult;
  using SEP::Testing::WriteJsonSummary;
  using SEP::Testing::WriteJUnitSummary;
  using SEP::Testing::EstimateRefinementOrder;

} // namespace Testing
} // namespace SEP3D


// ============================================================================
// Per-group registration functions
//
// Each function returns a std::vector<SEP3D::Testing::Descriptor> containing
// the fully populated descriptors for that group.  test/stage1.cpp calls
// each in order and feeds the result to the Registry constructor.
//
// STATUS at Phase R0:
//   Active: HARN, LAY, BLD, UTIL.  Later physics groups remain declarations in
//   the development plan until their production implementations exist.
//
// TO ADD A NEW GROUP:
//   1. Create test/individual-test/test_<group>.cpp
//   2. Implement the function declared below (remove the // comment)
//   3. Uncomment the call in test/stage1.cpp
//   4. Uncomment the make target in the makefile
//
// NAMING CONVENTION:
//   Group IDs are two-to-five characters, all-caps, no underscores.
//   Test IDs are the group ID followed by a two-digit number: LAY01, MSH3D04.
//   The runner sorts by ID, so group prefixes also control display order.
// ============================================================================

// ---- Phase R0 standalone groups (active) -----------------------------------

// HARN — harness self-verification tests
// Verifies the runner's own exit-code, JSON/JUnit writer, and empty-registry
// behaviour.  These run before any physics test so a broken harness is
// caught immediately.
std::vector<SEP3D::Testing::Descriptor> RegisterHarnessTests();

// LAY — layering-boundary tests
// Verifies core/ and background/ contain no AMPS symbols (LAY01) and that
// the grep used by LAY01 actually detects a violation (LAY02 negative control).
std::vector<SEP3D::Testing::Descriptor> RegisterLayeringTests(
    const std::string& selfPath = {});

// BLD — build-structure tests
// BLD01: the Stage-1 binary's symbol table contains no AMPS or MPI symbols.
// Uses nm to inspect the binary at the path provided.
std::vector<SEP3D::Testing::Descriptor> RegisterBuildTests(
    const std::string& binaryPath = {});

// UTIL - shared-kernel frozen record (Step 3)
// UTIL02: reproduce a fixed set of sep_common.a kernel calls and compare
// byte-for-byte with test/frozen/S03_kernels.txt.
std::vector<SEP3D::Testing::Descriptor> RegisterKernelTests();


// ---- Future groups (uncomment at the step indicated) -----------------------

// Step 7  — mesh resolution law
// std::vector<SEP3D::Testing::Descriptor> RegisterMeshTests();

// Step 12 — background field providers
// std::vector<SEP3D::Testing::Descriptor> RegisterBackgroundTests();

// Step 18 — diffusion coefficient registry and tensor assembly
// std::vector<SEP3D::Testing::Descriptor> RegisterCoefficientTests();

// Step 20 — substep selector and timestep limits
// std::vector<SEP3D::Testing::Descriptor> RegisterTimestepTests();

// Step 21 — Parker 3-D mover core
// std::vector<SEP3D::Testing::Descriptor> RegisterParkerMoverTests();

// Step 22 — focused-transport 3-D mover core
// std::vector<SEP3D::Testing::Descriptor> RegisterFTEMoverTests();

// Step 25 — sampling products
// std::vector<SEP3D::Testing::Descriptor> RegisterSamplingTests();

// Step 27 — reproducibility across MPI configurations (Stage 2 only)
// std::vector<SEP3D::Testing::Descriptor> RegisterReproducibilityTests();

#endif // SEP3D_CORE_SEP3D_TEST_REGISTRY_H

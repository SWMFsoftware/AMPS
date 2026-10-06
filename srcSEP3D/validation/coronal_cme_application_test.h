// ============================================================================
// Native AMPS test boundary for the shared SEP coronal-CME model.
//
// This header deliberately contains only neutral value records.  main_lib.cpp
// captures authoritative AMPS state after the normal production lifecycle has
// initialized it; coronal_cme_application_test.cpp then evaluates named tests
// through public sep_coronal_cme APIs and writes portable JSON evidence.
// Production execution never enters this boundary unless --test or
// --all-tests or --test-suite was explicitly selected.
// ============================================================================

#ifndef SEP3D_VALIDATION_CORONAL_CME_APPLICATION_TEST_H
#define SEP3D_VALIDATION_CORONAL_CME_APPLICATION_TEST_H

#include "../core/sep3d_types.h"

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace SEP3D { namespace Validation {

enum class NativeTestStatus { Pass, Fail, Skip, Error };

struct NativeTestDescriptor {
  std::string id;
  std::string name;
  std::string description;
  // Explicit membership makes future suite selection independent of ID naming.
  std::string suite = "amps";
};

struct NativeSpeciesState {
  int compiledSlot = -1;
  std::string chemicalSymbol;
  double massKg = 0.0;
  double chargeC = 0.0;
  double timeStepS = 0.0;
  double particleWeight = 0.0;
};

// Every field is captured from the live AMPS application after amps_init().
// Values that are already collective are prefixed global; booleans are reduced
// with logical AND unless their name explicitly describes an optional mode.
struct NativeApplicationState {
  int mpiRankCount = 0;
  int expectedMpiRanks = 0;
  std::uint64_t globalAllocatedBlocks = 0;
  std::uint64_t globalPhysicalCells = 0;
  std::uint64_t plannedActiveLeaves = 0;
  std::uint64_t plannedInactiveLeaves = 0;
  std::uint64_t plannedSolarInteriorLeaves = 0;
  std::uint64_t backgroundGeneration = 0;
  std::uint64_t shockBackgroundGeneration = 0;
  std::uint64_t completedSteps = 0;
  std::uint32_t initializationMask = 0;
  bool solarBoundaryRegistered = false;
  // Installation includes a verified identity plan. Actual deactivation and
  // successful owner-block allocation are recorded independently below.
  bool activeMaskInstalled = false;
  bool activeRegionPruningApplied = false;
  bool activeRegionAllocationVerified = false;
  bool backgroundReady = false;
  bool turbulenceReady = false;
  bool shockRequired = false;
  bool shockReady = false;
  bool sourceEnabled = false;
  bool finiteBackgroundAndTurbulence = false;
  bool finiteBackgroundDerivatives = false;
  bool initializationProductsExist = false;
  bool initializationProductsFinite = false;
  bool mpiFingerprintConsistent = false;
  bool restartConfigured = false;
  // Read from the live provider, actual center-node bytes and received remote
  // blocks. These are collective checks, never synthetic mesh evidence.
  bool runtimeMeshOwnedFieldsMatch = false; // All owner cells, compared fields.
  bool runtimeMeshGhostFieldsMatch = false; // One center per received block.
  bool runtimeMeshProviderMatch = false; // Prepared epoch/tag and live readiness.
  std::uint64_t runtimeMeshGhostCellsChecked = 0; // Global representative count.
  std::uint64_t runtimeMeshPublishedUpdates = 0; // Completed refresh commits.
  std::uint64_t runtimeMeshExpectedUpdates = 0; // tick/backgroundCadenceSteps.
  bool runtimeCollectiveRollbackVerified = false;
  // Reduced-provider evidence is captured from the same adapter/epoch that
  // filled owner and received-ghost storage.  The front is replicated shared
  // analytical state; plasma arrays remain the undisturbed ambient reference.
  bool reducedProviderSelected = false;
  bool reducedGeometricEndpointReached = false;
  bool reducedApexShockAccepted = false;
  bool reducedEndpointObserverGeometricHit = false;
  bool reducedEndpointObserverShockAccepted = false;
  std::uint64_t reducedFrontGeneration = 0;
  std::uint64_t reducedAmbientGeneration = 0;
  double reducedEpochS = 0.0;
  double reducedApexRadiusM = 0.0;
  double reducedApexSpeedMPerS = 0.0;
  double reducedEndpointTimeS = 0.0;
  double reducedAcceptedAreaM2 = 0.0;
  double reducedNumericalFailureAreaM2 = 0.0;
  std::string reducedEventIdentity;
  std::string reducedPhase;
  std::string reducedEndpointObserverStatus;
  std::uint64_t globalParticleCount = 0;
  std::uint64_t globalInjectedParticleCount = 0;
  bool zeroParticleAllocationRequested = false;
  std::string configurationFingerprint;
  // Captured from the configured production authorities. Generic readiness
  // checks cannot prove that a particular coronal model was selected: the
  // provider names must accompany the evidence instead of being inferred
  // from the SCCM test prefix or an input filename.
  int inputSchemaVersion = 0;
  std::string backgroundAuthority = "unreported";
  std::string shockAuthority = "unreported";
  std::string activeRegionMode;
  std::vector<NativeSpeciesState> species;
};

struct NativeTestResult {
  std::string id;
  std::string name;
  NativeTestStatus status = NativeTestStatus::Error;
  std::string message;
  std::vector<std::pair<std::string, double>> metrics;
  std::vector<std::string> artifacts;
};

const std::vector<NativeTestDescriptor>& CoronalCmeNativeTests();
// Implemented by main_lib.cpp because it alone may read the private AMPS
// application state.  The capture performs collective reductions and must be
// called by every rank at the same joined lifecycle boundary.
Core::Status CaptureNativeApplicationState(
    int expectedMpiRanks, NativeApplicationState* state);
Core::Status SelectCoronalCmeNativeTests(
    bool allTests, const std::vector<std::string>& requested,
    std::vector<NativeTestDescriptor>* selected,
    const std::string& suite = "");
std::vector<NativeTestResult> EvaluateCoronalCmeNativeTests(
    const NativeApplicationState& state,
    const std::vector<NativeTestDescriptor>& selected);
Core::Status WriteNativeTestJson(
    const std::string& path, const NativeApplicationState& state,
    const std::vector<NativeTestResult>& results);
int NativeTestExitCode(const std::vector<NativeTestResult>& results);
const char* Name(NativeTestStatus status) noexcept;

} }  // namespace SEP3D::Validation

#endif  // SEP3D_VALIDATION_CORONAL_CME_APPLICATION_TEST_H

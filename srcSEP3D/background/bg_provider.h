// ============================================================================
// srcSEP3D/background/bg_provider.h
//
// Abstract background-field provider contract.
//
// LAYER: L1 (background models).
//   May include L0 (core/) headers.
//   Must NOT include pic.h, mpi.h, or any AMPS symbol.
//
// PHASE B STATUS:
//   The interface is active. Analytic Parker and imported SWMF/AWSoM
//   implementations use it to build immutable BackgroundSnapshot objects.
//
// PURPOSE:
//   Every source of background magnetic field and plasma state (analytic
//   Parker spiral, SWCME, external Python script) implements this interface.
//   The movers, the cell-filling code, and the tests all speak only this
//   language; swapping one provider for another is a one-line change in the
//   registry.
//
// KEY DESIGN CHOICES:
//
//   1. prepare() / evaluate() split.
//      prepare(t) is called once per background-update cadence; evaluate()
//      is called once per particle per substep.  This mirrors SWCME's
//      prepared-state lifecycle and prevents prepare() from accidentally
//      being called in a hot loop.
//
//   2. BackgroundSample carries a 'valid' flag.
//      A point may have no field value (e.g. inside the inner boundary,
//      or a cell not found by an external script).  The flag is carried
//      by all providers so no downstream code has a special case.
//
//   3. Batch evaluation.
//      EvaluateBatch() evaluates N points in one call.  This is the path
//      used for cell-centre filling, where N is typically 4^3 or 8^3 per
//      AMR block.  The single-point Evaluate() is for particle substeps.
//
// FULL SPECIFICATION: see BACKGROUND_FIELD.md.
// ============================================================================

#ifndef SEP3D_BG_PROVIDER_H
#define SEP3D_BG_PROVIDER_H

#include "../core/sep3d_types.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace SEP3D {
namespace Background {

// ----------------------------------------------------------------------------
// BackgroundSample — the values the provider delivers at one point.
//
// All quantities are in SI.  Fields marked (optional) are filled only when
// the provider advertises the corresponding capability flag; callers must
// check the capability before reading them.
// ----------------------------------------------------------------------------
struct BackgroundSample {
  Core::Status status;

  // true  = field is defined at this point; all scalar/vector fields below
  //         are valid.
  // false = no field at this point; the caller must apply the configured
  //         missing-value policy (mark-invalid, parker-fallback, fail, zero).
  bool  valid = false;

  // Magnetic field
  Core::Vec3   B;          // [T]
  double       absB = 0.0; // |B| [T]
  Core::Vec3   bHat;       // B/|B|, zero vector when absB == 0

  // Magnetic-field gradient tensor ∂B_i/∂x_j [T/m]   (optional)
  Core::Tensor3 gradB;

  // Precomputed scalars derived from the field (optional; the provider
  // may supply analytic values more accurate than the central-difference
  // stencil would produce)
  double divBhat      = 0.0;  // ∇·b̂  = -1/L_B  [1/m]
  double focusingLenM = 0.0;  // L_B = -(dln|B|/ds)^-1 [m]; +inf for uniform B
  Core::Vec3 curvature;       // (b̂·∇)b̂  [1/m]

  // Plasma bulk velocity [m/s]
  Core::Vec3   U;
  Core::Tensor3 gradU;              // dU_i/dx_j [1/s]
  double       divU = 0.0;         // ∇·U [1/s]
  double       fieldAlignedStrain = 0.0;  // b̂b̂:∇U [1/s]

  // Plasma state
  double numberDensityM3 = 0.0;  // n [m^-3]
  double temperatureK    = 0.0;  // T [K]
  double pressurePa      = 0.0;  // p [Pa]
  double alfvenSpeedMpS  = 0.0;  // v_A [m/s]

  // Snapshot identity — immutable once the provider is prepared.
  // The mover checks this to detect a stale sample.
  uint64_t generation = 0;

  // Configuration digest — hash of the provider's resolved configuration.
  // Stored in the run manifest; compared on restart.
  uint64_t configurationDigest = 0;
};


// ----------------------------------------------------------------------------
// ProviderCapabilities — bit flags advertising what the provider fills.
// Check these before reading optional fields in BackgroundSample.
// ----------------------------------------------------------------------------
struct ProviderCapabilities {
  bool hasAnalyticGradB   = false; // ∂B_i/∂x_j computed analytically
  bool hasAnalyticDivBhat = false; // ∇·b̂ and L_B computed analytically
  bool hasAnalyticCurvature = false; // (b̂·∇)b̂ computed analytically
  bool hasAnalyticDivU    = false; // ∇·U computed analytically
  bool hasFieldAlignedStrain = false; // b̂b̂:∇U computed analytically
  bool hasPlasmaState     = false; // n, T, p, v_A filled
  bool supportsBatchEval  = false; // EvaluateBatch() is implemented
};

// PythonInterpolator is a reserved provenance value for the future batched
// external-model bridge.  No current provider may publish it; recognizing it
// in the type system prevents a future implementation from masquerading as an
// analytic or coupled snapshot while its protocol is being introduced.
enum class ProviderKind { AnalyticParker, PythonInterpolator, SwmfAwsom };
enum class StorageOwnership { ModelOwned, ImportedReadOnly };

// Metadata are frozen by Prepare() and copied into every snapshot.  Keeping
// frame, ownership, epoch, and configuration identity next to the numerical
// values prevents an imported field from being silently reused in the wrong
// coordinate system or after the coupled epoch changes.
struct SnapshotMetadata {
  ProviderKind provider = ProviderKind::AnalyticParker;
  StorageOwnership ownership = StorageOwnership::ModelOwned;
  double epochS = 0.0;
  double validFromS = 0.0;
  double validUntilS = 0.0;
  std::uint64_t generation = 0;
  std::string coordinateFrame;
  std::string providerIdentity;
  std::string configurationFingerprint;
};


// ----------------------------------------------------------------------------
// BackgroundProvider — the abstract interface every provider implements.
//
// Concrete implementations (Phase B):
//   bg_parker.cpp   — analytic Parker spiral
//   bg_swmf.cpp     — read-only SWMF/AWSoM import adapter
//   external-script provider — RESERVED until its own release gate
// ----------------------------------------------------------------------------
class BackgroundProvider {
public:
  virtual ~BackgroundProvider() {}

  // Human-readable name used in log messages and the run manifest.
  virtual const char* CanonicalName() const = 0;

  // Validate the configuration before any preparation.
  // Returns an error Status if the configuration is inconsistent.
  virtual Core::Status Validate() const = 0;

  // Advance the provider's internal state to simulation time timeS.
  // Called once per background-update cadence, on every MPI rank.
  // May be expensive; must never be called from a particle loop.
  virtual Core::Status Prepare(double timeS) = 0;

  // Available only after successful Prepare().  Implementations must leave
  // their previous metadata and values unchanged when a new candidate fails.
  virtual const SnapshotMetadata* PreparedMetadata() const = 0;

  // Evaluate the background at one point x [m] in the heliocentric frame.
  // Returns immediately; the caller must have called Prepare() first.
  virtual BackgroundSample Evaluate(const Core::Vec3& x) const = 0;

  // Evaluate N points.  x_m, y_m, z_m are parallel arrays of length n.
  // out must be pre-allocated to length n.
  // Default implementation calls Evaluate() in a loop; providers override
  // this to use SIMD or batched analytic evaluation.
  virtual Core::Status EvaluateBatch(
      const double* x_m, const double* y_m, const double* z_m,
      std::size_t n, BackgroundSample* out) const {
    std::vector<Core::Status> statuses(n);
    return EvaluateBatchDetailed(x_m, y_m, z_m, n, out, statuses.data());
  }

  // Detailed batch evaluation has transactional per-sample output.  A valid
  // candidate replaces out[i]; a failed candidate writes only status[i] and
  // leaves out[i] byte-for-byte unchanged.  The aggregate return is the first
  // failed status, allowing a builder to reject the whole snapshot without
  // losing the identity of each bad sample.
  virtual Core::Status EvaluateBatchDetailed(
      const double* x_m, const double* y_m, const double* z_m,
      std::size_t n, BackgroundSample* out, Core::Status* status) const;

  // Returns the resolved configuration as a JSON-like manifest string.
  // Stored in the run manifest; compared on restart.
  virtual std::string ResolvedManifest() const = 0;

  // Returns the set of capabilities this provider fills.
  virtual ProviderCapabilities Capabilities() const = 0;
};

} // namespace Background
} // namespace SEP3D

#endif // SEP3D_BG_PROVIDER_H

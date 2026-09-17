// ============================================================================
// srcSEP3D/background/bg_provider.h
//
// Abstract background-field provider contract.
//
// LAYER: L1 (background models).
//   May include L0 (core/) headers.
//   Must NOT include pic.h, mpi.h, or any AMPS symbol.
//
// STEP 1 STATUS:
//   This is the forward declaration / stub that Step 1 puts in place so
//   that SEP3D.h can include it and the layering guard can verify it has
//   no AMPS dependency.  The full implementation is added in Steps 12–14.
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
// FULL SPECIFICATION: see BACKGROUND_FIELD.md (written at Step 17).
// ============================================================================

#ifndef SEP3D_BG_PROVIDER_H
#define SEP3D_BG_PROVIDER_H

#include "../core/sep3d_types.h"

#include <cstddef>
#include <string>

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


// ----------------------------------------------------------------------------
// BackgroundProvider — the abstract interface every provider implements.
//
// Concrete implementations (Steps 12–14):
//   bg_parker.cpp   — analytic Parker spiral
//   bg_swcme.cpp    — SWCME prepared-state adapter
//   bg_external_script.cpp — external Python script (RESERVED: Step WP28)
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
    for (std::size_t i = 0; i < n; ++i)
      out[i] = Evaluate({x_m[i], y_m[i], z_m[i]});
    return Core::Status::OK();
  }

  // Returns the resolved configuration as a JSON-like manifest string.
  // Stored in the run manifest; compared on restart.
  virtual std::string ResolvedManifest() const = 0;

  // Returns the set of capabilities this provider fills.
  virtual ProviderCapabilities Capabilities() const = 0;
};

} // namespace Background
} // namespace SEP3D

#endif // SEP3D_BG_PROVIDER_H

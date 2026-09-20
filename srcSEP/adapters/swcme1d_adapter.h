// ============================================================================
// srcSEP -> canonical SWCME adapter.
//
// This public boundary is deliberately provider-neutral: it contains no
// SWCME include and exposes no canonical SWCME type. AMPS translation units
// may therefore include it without inheriting src/models/swcme as a public
// include dependency. Provider-specific types remain confined to the .cpp.
//
// D01 fail-closed contract
// ------------------------
// A background query returns one immutable sample *or* an explicit failure.
// It never writes through caller-owned output references and never substitutes
// zero for a malformed physical field. Radius clamping and diagnostic fallback
// are named, opt-in policies; every recovered query is counted.
// ============================================================================

#ifndef SEP_ADAPTERS_SWCME1D_ADAPTER_H
#define SEP_ADAPTERS_SWCME1D_ADAPTER_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace SEP {
namespace SW1DAdapter {

enum class Scenario { Fast, Slow };

// Strict is the production default. ClampRadius may only repair a radius below
// the canonical inner domain. DiagnosticFallback may replace any failed query
// with a user-supplied, validated sample and is intended for controlled
// sensitivity experiments, not production forecasting.
enum class FailurePolicy { Strict, ClampRadius, DiagnosticFallback };

enum class StatusCode {
  Success,
  RecoveredWithRadiusClamp,
  RecoveredWithDiagnosticFallback,
  ModelNotPrepared,
  NonFiniteInput,
  RadiusOutsideDomain,
  NonFiniteOutput,
  NonPositiveDensity,
  InvalidSpeed,
  InvalidDivergence,
  CanonicalFailure,
  ConfigurationFailure,
  PreparationFailure,
  InvalidDiagnosticFallback
};

const char* StatusCodeName(StatusCode code);
const char* FailurePolicyName(FailurePolicy policy);

// Provider-neutral status used for state preparation and policy setup.
struct Status {
  StatusCode code = StatusCode::Success;
  std::string detail;
  bool ok() const { return code == StatusCode::Success; }
};

// Immutable physical value object. All fields use SI units. There are no
// mutators, so a successful query can be copied into a snapshot or source
// calculation without observing a partially updated record.
class BackgroundSample {
 public:
  BackgroundSample() = default;
  BackgroundSample(double number_density_m3, double speed_m_s,
                   double divergence_s_inv)
      : number_density_m3_(number_density_m3),
        speed_m_s_(speed_m_s),
        divergence_s_inv_(divergence_s_inv) {}

  double number_density_m3() const { return number_density_m3_; }
  double speed_m_s() const { return speed_m_s_; }
  double divergence_s_inv() const { return divergence_s_inv_; }

 private:
  double number_density_m3_ = 0.0;
  double speed_m_s_ = 0.0;
  double divergence_s_inv_ = 0.0;
};

struct QueryResult {
  StatusCode code = StatusCode::ModelNotPrepared;
  StatusCode original_code = StatusCode::ModelNotPrepared;
  BackgroundSample sample;
  double requested_radius_m = 0.0;
  double evaluated_radius_m = 0.0;
  double epoch_seconds = 0.0;
  std::uint64_t source_state_id = 0;
  std::string failed_field;
  std::string detail;

  bool ok() const {
    return code == StatusCode::Success ||
           code == StatusCode::RecoveredWithRadiusClamp ||
           code == StatusCode::RecoveredWithDiagnosticFallback;
  }
  bool recovered() const { return ok() && code != StatusCode::Success; }
};

struct Diagnostics {
  std::uint64_t successful_queries = 0;
  std::uint64_t failed_queries = 0;
  std::uint64_t radius_clamps = 0;
  std::uint64_t diagnostic_fallbacks = 0;
  std::uint64_t prepared_state_id = 0;
  double prepared_epoch_seconds = 0.0;
};

// Provider-neutral assignment transport used by both the PARAM parser and
// parser-free coupled hosts. Keys and units are interpreted only by the
// canonical resolver in src/models/swcme/swcme1d_input.hpp.
struct ParameterAssignment {
  std::string key;
  std::string value;
  std::string origin;
  std::size_t line = 0;
};

struct ConfigurationRequest {
  Scenario preset = Scenario::Fast;
  std::vector<ParameterAssignment> input_assignments;
  std::vector<ParameterAssignment> command_line_assignments;
};

// Stable, provider-neutral view of the frozen canonical configuration. It is
// sufficient for run fingerprints and for the legacy srcSEP injection bridge;
// the complete canonical Params/SpectrumConfig types remain private.
struct ConfigurationSummary {
  std::string preset;
  std::string fingerprint;
  std::string normalized_manifest;
  double launch_epoch_s = 0.0;
  double valid_from_s = 0.0;
  double valid_until_s = 0.0;
  double source_particle_mass_kg = 0.0;
  int source_charge_number = 0;
  double source_energy_min_MeV = 0.0;
  double source_energy_max_MeV = 0.0;
  double source_reference_energy_MeV = 0.0;
  double source_injection_efficiency = 0.0;
  double relative_source_weight_per_area = 0.0;
  // Startup consistency fields: the finite Parker mesh and SWCME background
  // must use exactly one wind/rotation/latitude law.  The signed radial field
  // is derived from SWCME's total 1-AU magnitude using its documented
  // reference latitude and polarity; application code must not reinterpret
  // the total magnitude as Br.
  double ambient_wind_speed_m_per_s = 0.0;
  double solar_rotation_rate_rad_per_s = 0.0;
  double parker_source_radius_m = 0.0;
  double parker_reference_sin_theta = 0.0;
  int parker_radial_polarity = 0;
  double parker_radial_field_at_one_au_t = 0.0;
};

// Configure the process-wide 1-D SW+CME model. Call PrepareState() after
// configuration and once for each physical background epoch.
void Configure(Scenario scenario);

// General D02 path. It applies the named preset, then the staged/input-file
// layer, then command-line/programmatic assignments. Failure leaves the prior
// model/configuration untouched.
Status Configure(const ConfigurationRequest& request);
ConfigurationSummary GetConfigurationSummary();

// The legacy PIC post-compile parser stages raw assignments here. Coupled hosts
// should skip staging and call Configure(request) directly, keeping their
// library path parser-free.
void ClearStagedInputAssignments();
Status StageInputAssignment(const ParameterAssignment& assignment);

// Preparation is transactional. If the canonical model rejects the new epoch,
// the last valid state and its ID remain unchanged and no background snapshot
// should be published by the caller.
Status PrepareState(double epoch_seconds);

// Configure the explicit D01 recovery policy. DiagnosticFallback additionally
// requires a valid strictly-positive density and speed sample.
Status SetFailurePolicy(FailurePolicy policy);
Status SetDiagnosticFallback(const BackgroundSample& sample);
FailurePolicy GetFailurePolicy();

// Optional monotonic stabilization applied only inside the modeled sheath.
// This is part of the physical profile and is distinct from radius clamping.
void EnableSheathClamp(bool enabled = true);

// Provider-neutral scalar views used by the rest of srcSEP. These accessors
// throw when no state is prepared; returning 0/1 would silently manufacture a
// shock state and violate the same fail-closed contract as QueryAtRadius().
double ShockRadiusM();
double ShockSpeedMPerS();
double CompressionRatio();
double DlnB_Dr_at_r(double radius_m);
double RelativeSourceWeightPerArea();

// Query density [m^-3], radial speed [m/s], and divergence [s^-1]. The result
// owns all three fields, so failure cannot partially mutate caller state.
QueryResult QueryAtRadius(double radius_m,
                          bool apply_sheath_clamp = true);

// Public validation seam used by focused tests and by the adapter after every
// canonical evaluation. It distinguishes the physical field responsible for
// rejection without requiring a deliberately corrupted canonical model.
Status ValidateSample(const BackgroundSample& sample);

Diagnostics GetDiagnostics();
void ResetDiagnostics();

// Include rank, epoch, location, field, state ID, and canonical detail in one
// deterministic diagnostic suitable for an AMPS fatal error or test log.
std::string FormatFailure(const QueryResult& result, int mpi_rank);

// Diagnostics remain owned by the adapter so their SWCME StepState argument
// never becomes part of the srcSEP application interface.
void WriteRadialProfileFromR(const double* radii_m, int count,
                             const char* file_name, double epoch_seconds);
void WriteShockVsTime(double duration_seconds, int sample_count,
                      const char* file_name);

}  // namespace SW1DAdapter
}  // namespace SEP

#endif  // SEP_ADAPTERS_SWCME1D_ADAPTER_H

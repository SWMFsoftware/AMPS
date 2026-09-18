#include "sep_swcme_validation.h"

#include "../adapters/swcme1d_adapter.h"
#include "sep_production_mover.h"

// The canonical resolver owns SWCME keys, units, precedence, and physics
// validation.  The srcSEP makefile provides src/models/swcme as an include
// root; no copy or compatibility header is kept in the application tree.
#include "swcme1d_input.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

// Accumulate all assertions so a single registry record reports the complete
// failure set.  This is preferable to throwing at the first discrepancy: D01
// and D02 are contract audits, and seeing every violated field/policy in one
// isolated process substantially shortens diagnosis.  assertion_failures is
// also understood by the shared registry and can never be hidden behind an
// accidentally returned PASS status.
class CheckSet {
 public:
  void Expect(bool condition, const std::string& description) {
    ++total_;
    if (condition) return;
    failures_.push_back(description);
  }

  SEP::Testing::Result Finish(const std::string& pass_message,
                              const std::string& failure_prefix) const {
    SEP::Testing::Result result;
    result.status = failures_.empty() ? SEP::Testing::Status::Pass
                                      : SEP::Testing::Status::Fail;
    if (failures_.empty()) {
      result.message = pass_message;
    } else {
      std::ostringstream message;
      message << failure_prefix << ": ";
      for (std::size_t i = 0; i < failures_.size(); ++i) {
        if (i != 0) message << "; ";
        message << failures_[i];
      }
      result.message = message.str();
    }
    result.metrics.push_back({"assertions_total",
                              static_cast<double>(total_), 0.0, ">=", "count"});
    result.metrics.push_back({"assertion_failures",
                              static_cast<double>(failures_.size()), 0.0,
                              "<=", "count"});
    return result;
  }

 private:
  std::size_t total_ = 0;
  std::vector<std::string> failures_;
};

bool Close(double left, double right, double relative = 1.0e-13) {
  return std::fabs(left - right) <=
      relative * std::max(1.0, std::max(std::fabs(left), std::fabs(right)));
}

swcme::input1d::Assignment CanonicalAssignment(
    const std::string& key, const std::string& value, std::size_t line = 1) {
  swcme::input1d::Assignment assignment;
  assignment.key = key;
  assignment.value = value;
  assignment.origin = "D02 native-registry fixture";
  assignment.line = line;
  return assignment;
}

SEP::SW1DAdapter::ParameterAssignment AdapterAssignment(
    const std::string& key, const std::string& value, std::size_t line = 1) {
  SEP::SW1DAdapter::ParameterAssignment assignment;
  assignment.key = key;
  assignment.value = value;
  assignment.origin = "D02 parser-free native-registry fixture";
  assignment.line = line;
  return assignment;
}

// Every callback changes process-owned SWCME state.  Python `--all` executes
// one registry ID per process, but direct users may select several IDs at once.
// Restore a documented fast/strict state so later callbacks do not inherit a
// fallback policy, staged PARAM assignment, or validity interval from a test.
// Cleanup is best-effort only in the destructor; each callback explicitly
// checks restoration before returning so a cleanup defect becomes a FAIL.
bool RestoreDefaultSwcmeState(std::string* detail) {
  try {
    SEP::SW1DAdapter::ClearStagedInputAssignments();
    SEP::SW1DAdapter::Configure(SEP::SW1DAdapter::Scenario::Fast);
    const SEP::SW1DAdapter::Status prepared =
        SEP::SW1DAdapter::PrepareState(0.0);
    if (!prepared.ok()) {
      if (detail) *detail = prepared.detail;
      return false;
    }
    const SEP::SW1DAdapter::Status strict =
        SEP::SW1DAdapter::SetFailurePolicy(
            SEP::SW1DAdapter::FailurePolicy::Strict);
    if (!strict.ok()) {
      if (detail) *detail = strict.detail;
      return false;
    }
    SEP::SW1DAdapter::ResetDiagnostics();
    return true;
  } catch (const std::exception& exception) {
    if (detail) *detail = exception.what();
    return false;
  } catch (...) {
    if (detail) *detail = "non-standard exception during SWCME restoration";
    return false;
  }
}

SEP::Testing::Result RunD01FailClosedBackground() {
  using namespace SEP::SW1DAdapter;
  CheckSet checks;

  try {
    // A D01 process begins before standalone model initialization.  A query at
    // that point must retain its typed failure; manufacturing a zero-valued
    // sample would let shock injection silently proceed with invalid physics.
    QueryResult query = QueryAtRadius(1.0e9);
    checks.Expect(!query.ok() && query.code == StatusCode::ModelNotPrepared,
                  "unprepared model was not rejected explicitly");

    // Exercise corruption classes through the public validation seam.  Valid
    // canonical SWCME output cannot normally generate these values, so direct
    // samples are the deterministic way to cover each fail-closed branch.
    checks.Expect(ValidateSample(BackgroundSample(-1.0, 4.0e5, 0.0)).code ==
                      StatusCode::NonPositiveDensity,
                  "negative density was accepted");
    checks.Expect(ValidateSample(BackgroundSample(5.0e6, 0.0, 0.0)).code ==
                      StatusCode::InvalidSpeed,
                  "zero solar-wind speed was accepted");
    checks.Expect(ValidateSample(BackgroundSample(
                      std::numeric_limits<double>::quiet_NaN(), 4.0e5, 0.0)).code ==
                      StatusCode::NonFiniteOutput,
                  "NaN density was accepted");
    checks.Expect(ValidateSample(BackgroundSample(
                      5.0e6, 4.0e5,
                      std::numeric_limits<double>::infinity())).code ==
                      StatusCode::InvalidDivergence,
                  "infinite velocity divergence was accepted");

    Configure(Scenario::Fast);
    Status status = PrepareState(0.0);
    checks.Expect(status.ok(), "fast preset did not prepare at launch epoch");
    const double fast_radius = status.ok() ? ShockRadiusM() : -1.0;
    const Diagnostics fast_state = GetDiagnostics();
    checks.Expect(status.ok() && Close(ShockSpeedMPerS(), 1.9e6),
                  "fast preset did not retain 1900 km/s launch speed");

    // The exact canonical radius is valid.  Its immediately adjacent floating
    // point below the boundary must fail under Strict rather than being
    // silently projected back onto the domain.
    const double inner_boundary = 1.05 * 6.957e8;
    query = QueryAtRadius(inner_boundary);
    checks.Expect(query.ok() && query.code == StatusCode::Success &&
                      query.sample.number_density_m3() > 0.0 &&
                      query.sample.speed_m_s() > 0.0,
                  "exact 1.05-R_sun boundary was not a valid sample");
    query = QueryAtRadius(std::nextafter(inner_boundary, 0.0));
    checks.Expect(!query.ok() && query.code == StatusCode::RadiusOutsideDomain,
                  "strict policy did not reject the adjacent sub-boundary radius");

    // Failed state preparation is transactional.  The last prepared physical
    // state and its monotonically increasing identity must remain unchanged.
    status = PrepareState(std::numeric_limits<double>::quiet_NaN());
    const Diagnostics after_bad_prepare = GetDiagnostics();
    checks.Expect(!status.ok() && status.code == StatusCode::PreparationFailure &&
                      Close(ShockRadiusM(), fast_radius) &&
                      after_bad_prepare.prepared_state_id ==
                          fast_state.prepared_state_id,
                  "failed epoch preparation replaced the valid state or ID");

    status = SetFailurePolicy(FailurePolicy::ClampRadius);
    query = QueryAtRadius(0.5 * inner_boundary);
    checks.Expect(status.ok() && query.ok() &&
                      query.code == StatusCode::RecoveredWithRadiusClamp &&
                      Close(query.evaluated_radius_m, inner_boundary) &&
                      GetDiagnostics().radius_clamps == 1,
                  "explicit radius clamp was not typed and counted");

    const BackgroundSample fallback(7.0e6, 3.5e5, -2.0e-6);
    status = SetDiagnosticFallback(fallback);
    checks.Expect(status.ok(), "valid diagnostic fallback was rejected");
    status = SetFailurePolicy(FailurePolicy::DiagnosticFallback);
    query = QueryAtRadius(std::numeric_limits<double>::quiet_NaN());
    checks.Expect(status.ok() && query.ok() &&
                      query.code == StatusCode::RecoveredWithDiagnosticFallback &&
                      query.original_code == StatusCode::NonFiniteInput &&
                      Close(query.sample.number_density_m3(), 7.0e6) &&
                      GetDiagnostics().diagnostic_fallbacks == 1,
                  "diagnostic fallback was not explicit, immutable, and counted");
    const std::string formatted = FormatFailure(query, 3);
    checks.Expect(formatted.find("rank=3") != std::string::npos &&
                      formatted.find("epoch_s=") != std::string::npos &&
                      formatted.find("source_state_id=") != std::string::npos &&
                      formatted.find("field=radius_m") != std::string::npos,
                  "failure diagnostic omitted rank, epoch, field, or state ID");

    // Reconfiguration creates a new canonical model instead of mutating a
    // prepared one.  Probe several physical epochs because the slow scenario
    // can be at a model transition exactly at one selected time.
    Configure(Scenario::Slow);
    const double slow_epochs[] = {0.0, 3600.0, 6.0 * 3600.0, 24.0 * 3600.0};
    for (double epoch : slow_epochs) {
      status = PrepareState(epoch);
      if (status.ok()) break;
    }
    const double slow_speed = status.ok() ? ShockSpeedMPerS() : -1.0;
    checks.Expect(status.ok() && slow_speed > 3.8e5 && slow_speed <= 9.5e5,
                  "slow preset did not retain its decelerating 950 km/s trajectory");
    checks.Expect(GetFailurePolicy() == FailurePolicy::Strict &&
                      GetDiagnostics().radius_clamps == 0 &&
                      GetDiagnostics().diagnostic_fallbacks == 0,
                  "new event retained a recovery policy or recovery counters");
  } catch (const std::exception& exception) {
    checks.Expect(false, std::string("unexpected D01 exception: ") +
                             exception.what());
  } catch (...) {
    checks.Expect(false, "unexpected non-standard D01 exception");
  }

  std::string restore_detail;
  checks.Expect(RestoreDefaultSwcmeState(&restore_detail),
                "default SWCME restoration failed: " + restore_detail);
  SEP::Testing::Result result = checks.Finish(
      "D01 fail-closed query, recovery, diagnostics, and transaction checks passed",
      "D01 fail-closed background contract failed");
  result.configuration.push_back("failure_policy_default=strict");
  result.configuration.push_back("background_units=SI");
  result.configuration.push_back("adapter_state_restored=fast@0s;strict");
  return result;
}

SEP::Testing::Result RunD02CanonicalConfiguration() {
  namespace I = swcme::input1d;
  CheckSet checks;

  try {
    const I::ResolveResult fast = I::Resolve(I::Preset::Fast, {});
    checks.Expect(fast.ok() && Close(fast.configuration.model.V_sw_kms, 400.0) &&
                      Close(fast.configuration.model.V0_sh_kms, 1900.0) &&
                      Close(fast.configuration.model.Gamma_kmInv, 8.0e-8) &&
                      Close(fast.configuration.source.spectrum.
                                kinetic_energy_min_MeV, 0.1) &&
                      Close(fast.configuration.source.spectrum.
                                kinetic_energy_max_MeV, 500.0),
                  "fast preset no longer reproduces its background/source values");
    const I::ResolveResult slow = I::Resolve(I::Preset::Slow, {});
    checks.Expect(slow.ok() && Close(slow.configuration.model.V_sw_kms, 380.0) &&
                      Close(slow.configuration.model.V0_sh_kms, 950.0) &&
                      Close(slow.configuration.model.Gamma_kmInv, 3.0e-8),
                  "slow preset no longer reproduces its background values");

    I::Layer input;
    input.name = "input-file";
    input.assignments.push_back(
        CanonicalAssignment("ambient.wind_speed", "450000 m/s"));
    input.assignments.push_back(
        CanonicalAssignment("cme.launch_radius", "0.01 AU", 2));
    input.assignments.push_back(
        CanonicalAssignment("source.energy_max", "0.8 GeV", 3));
    I::Layer cli;
    cli.name = "command-line";
    cli.assignments.push_back(
        CanonicalAssignment("ambient.wind_speed", "500 km/s"));
    I::ResolveResult layered = I::Resolve(I::Preset::Fast, {input, cli});
    checks.Expect(!layered.ok() && layered.status.code == I::Code::InvalidUnit &&
                      layered.status.key == "source.energy_max" &&
                      layered.status.layer == "input-file" &&
                      layered.status.origin == "D02 native-registry fixture" &&
                      layered.status.line == 3,
                  "invalid unit diagnostic lost key/layer/origin/line provenance");
    input.assignments.back() =
        CanonicalAssignment("source.energy_max", "800 MeV", 3);
    layered = I::Resolve(I::Preset::Fast, {input, cli});
    checks.Expect(layered.ok() &&
                      Close(layered.configuration.model.V_sw_kms, 500.0) &&
                      Close(layered.configuration.model.r0_Rs,
                            0.01 * swcme::constants::AU_M /
                                swcme::constants::SOLAR_RADIUS_M) &&
                      Close(layered.configuration.source.spectrum.
                                kinetic_energy_max_MeV, 800.0),
                  "preset/input/CLI precedence or one-time unit conversion failed");

    I::Layer duplicate;
    duplicate.name = "input-file";
    duplicate.assignments = {
        CanonicalAssignment("cme.launch_speed", "1200 km/s"),
        CanonicalAssignment("CME.LAUNCH_SPEED", "1300 km/s", 2)};
    const I::ResolveResult duplicate_result =
        I::Resolve(I::Preset::Fast, {duplicate});
    checks.Expect(!duplicate_result.ok() &&
                      duplicate_result.status.code == I::Code::DuplicateKey,
                  "canonicalized duplicate keys were accepted in one layer");

    I::Layer integer_overflow;
    integer_overflow.name = "input-file";
    integer_overflow.assignments.push_back(CanonicalAssignment(
        "source.charge_number", "999999999999999999999"));
    const I::ResolveResult integer_overflow_result =
        I::Resolve(I::Preset::Fast, {integer_overflow});
    checks.Expect(!integer_overflow_result.ok() &&
                      integer_overflow_result.status.code == I::Code::InvalidValue,
                  "integer overflow reached narrowed canonical storage");

    I::Layer disabled_source;
    disabled_source.name = "input-file";
    disabled_source.assignments.push_back(
        CanonicalAssignment("source.injection_efficiency", "0"));
    const I::ResolveResult disabled_source_result =
        I::Resolve(I::Preset::Fast, {disabled_source});
    checks.Expect(disabled_source_result.ok() &&
                      Close(disabled_source_result.configuration.source.
                                injection_efficiency, 0.0),
                  "zero source efficiency did not explicitly disable injection");

    I::Layer unknown;
    unknown.name = "input-file";
    unknown.assignments.push_back(
        CanonicalAssignment("shock.magic_compression", "4"));
    const I::ResolveResult unknown_result =
        I::Resolve(I::Preset::Fast, {unknown});
    checks.Expect(!unknown_result.ok() &&
                      unknown_result.status.code == I::Code::UnknownKey,
                  "unknown canonical configuration key was accepted");

    I::Layer deprecated;
    deprecated.name = "input-file";
    deprecated.assignments.push_back(
        CanonicalAssignment("sheath.compression_floor", "2"));
    const I::ResolveResult deprecated_result =
        I::Resolve(I::Preset::Fast, {deprecated});
    checks.Expect(!deprecated_result.ok() &&
                      deprecated_result.status.code == I::Code::UnsupportedField,
                  "deprecated no-effect field was not explicitly unsupported");

    I::Layer incompatible;
    incompatible.name = "input-file";
    incompatible.assignments.push_back(
        CanonicalAssignment("shock.region_mode", "full_icme"));
    incompatible.assignments.push_back(
        CanonicalAssignment("shock.acceleration_mode", "source", 2));
    const I::ResolveResult incompatible_result =
        I::Resolve(I::Preset::Fast, {incompatible});
    checks.Expect(!incompatible_result.ok() &&
                      incompatible_result.status.code ==
                          I::Code::CanonicalValidationFailure,
                  "incompatible shock representation passed canonical validation");

    I::Layer data;
    data.name = "input-file";
    data.assignments.push_back(
        CanonicalAssignment("cme.kinematics", "data_driven"));
    data.assignments.push_back(
        CanonicalAssignment("cme.data_times", "0 s, 1 h, 2 h", 2));
    data.assignments.push_back(
        CanonicalAssignment("cme.data_radii", "1.05 Rs, 5 Rs, 10 Rs", 3));
    const I::ResolveResult data_result = I::Resolve(I::Preset::Fast, {data});
    checks.Expect(data_result.ok() &&
                      data_result.configuration.model.data_time_s.size() == 3 &&
                      Close(data_result.configuration.model.data_time_s[2], 7200.0),
                  "data-driven height-time table unit parsing failed");

    // Equivalent assignments must serialize in schema order, not input order;
    // otherwise restarts would reject physically identical configurations.
    I::Layer order_a;
    order_a.name = "input";
    order_a.assignments = {
        CanonicalAssignment("ambient.wind_speed", "450 km/s"),
        CanonicalAssignment("ambient.density_1au", "7 cm^-3", 2)};
    I::Layer order_b;
    order_b.name = "input";
    order_b.assignments = {
        CanonicalAssignment("ambient.density_1au", "7 cm^-3"),
        CanonicalAssignment("ambient.wind_speed", "450 km/s", 2)};
    const I::ResolveResult first = I::Resolve(I::Preset::Fast, {order_a});
    const I::ResolveResult second = I::Resolve(I::Preset::Fast, {order_b});
    checks.Expect(first.ok() && second.ok() &&
                      first.configuration.normalized_manifest ==
                          second.configuration.normalized_manifest &&
                      first.configuration.fingerprint ==
                          second.configuration.fingerprint,
                  "manifest/fingerprint depends on assignment order");

    // Exercise the public parser-free bridge used by coupled hosts.  The
    // command/programmatic layer intentionally overrides the staged PARAM
    // preset, while every public summary remains in documented SI/MeV units.
    SEP::SW1DAdapter::ClearStagedInputAssignments();
    SEP::SW1DAdapter::Status adapter_status =
        SEP::SW1DAdapter::StageInputAssignment(
            AdapterAssignment("preset", "slow"));
    SEP::SW1DAdapter::ConfigurationRequest request;
    request.command_line_assignments.push_back(
        AdapterAssignment("preset", "fast"));
    request.command_line_assignments.push_back(
        AdapterAssignment("event.launch_epoch", "1 h", 2));
    request.command_line_assignments.push_back(
        AdapterAssignment("event.valid_until", "2 h", 3));
    request.command_line_assignments.push_back(
        AdapterAssignment("source.energy_min", "250 keV", 4));
    request.command_line_assignments.push_back(
        AdapterAssignment("source.injection_efficiency", "0.001", 5));
    adapter_status = SEP::SW1DAdapter::Configure(request);
    const SEP::SW1DAdapter::ConfigurationSummary summary =
        SEP::SW1DAdapter::GetConfigurationSummary();
    checks.Expect(adapter_status.ok() && summary.preset == "FAST" &&
                      Close(summary.launch_epoch_s, 3600.0) &&
                      Close(summary.valid_from_s, 3600.0) &&
                      Close(summary.valid_until_s, 7200.0) &&
                      Close(summary.source_energy_min_MeV, 0.25) &&
                      Close(summary.source_injection_efficiency, 0.001) &&
                      summary.fingerprint.size() == 16 &&
                      summary.normalized_manifest.find("preset=FAST") !=
                          std::string::npos,
                  "parser-free adapter did not expose frozen values/provenance");

    SEP::SW1DAdapter::ConfigurationRequest invalid_adapter_request;
    invalid_adapter_request.command_line_assignments.push_back(
        AdapterAssignment("ambient.wind_speed", "400 furlong/day", 12));
    const SEP::SW1DAdapter::Status invalid_adapter_status =
        SEP::SW1DAdapter::Configure(invalid_adapter_request);
    checks.Expect(!invalid_adapter_status.ok() &&
                      invalid_adapter_status.detail.find(
                          "origin=D02 parser-free native-registry fixture") !=
                          std::string::npos &&
                      invalid_adapter_status.detail.find("line=12") !=
                          std::string::npos &&
                      SEP::SW1DAdapter::GetConfigurationSummary().fingerprint ==
                          summary.fingerprint,
                  "failed adapter update lost provenance or mutated configuration");
    adapter_status = SEP::SW1DAdapter::PrepareState(3600.0);
    const std::uint64_t state_id =
        SEP::SW1DAdapter::GetDiagnostics().prepared_state_id;
    const SEP::SW1DAdapter::Status outside =
        SEP::SW1DAdapter::PrepareState(3599.0);
    checks.Expect(adapter_status.ok() && !outside.ok() &&
                      SEP::SW1DAdapter::GetDiagnostics().prepared_state_id ==
                          state_id,
                  "pre-launch epoch replaced the valid prepared state");

    SEP::SW1DAdapter::ConfigurationRequest unsupported_source;
    unsupported_source.command_line_assignments.push_back(AdapterAssignment(
        "source.normalization", "reference_differential_intensity"));
    unsupported_source.command_line_assignments.push_back(AdapterAssignment(
        "source.reference_intensity_si", "1e10 SI", 2));
    const SEP::SW1DAdapter::Status unsupported_status =
        SEP::SW1DAdapter::Configure(unsupported_source);
    checks.Expect(!unsupported_status.ok() &&
                      unsupported_status.detail.find("UNSUPPORTED_FIELD") !=
                          std::string::npos,
                  "unimplemented absolute source normalization was ignored");
  } catch (const std::exception& exception) {
    checks.Expect(false, std::string("unexpected D02 exception: ") +
                             exception.what());
  } catch (...) {
    checks.Expect(false, "unexpected non-standard D02 exception");
  }

  std::string restore_detail;
  checks.Expect(RestoreDefaultSwcmeState(&restore_detail),
                "default SWCME restoration failed: " + restore_detail);
  SEP::Testing::Result result = checks.Finish(
      "D02 canonical schema, units, precedence, provenance, and fingerprint checks passed",
      "D02 SWCME configuration contract failed");
  result.configuration.push_back("authority=preset<input<command-line/programmatic");
  result.configuration.push_back("canonical_owner=src/models/swcme");
  result.configuration.push_back("adapter_state_restored=fast@0s;strict");
  return result;
}

SEP::Testing::Result RunD03NativePreflight() {
  CheckSet checks;

  try {
    // D03's outer campaign appends one of these canonical names to every site
    // launch.  Verify the linked executable exposes exactly the intended 1-D
    // field-line algorithms and that each declares the coefficient contract
    // the campaign is supposed to exercise.
    const std::vector<SEP::Mover::Descriptor>& movers = SEP::Mover::Registry();
    const char* expected_names[] = {"parker", "fte-dmumu", "fte-mfp"};
    const SEP::Mover::CoefficientContract expected_coefficients[] = {
        SEP::Mover::CoefficientContract::SpatialDiffusion,
        SEP::Mover::CoefficientContract::PitchAngleDiffusion,
        SEP::Mover::CoefficientContract::MeanFreePath};
    checks.Expect(movers.size() == 3,
                  "linked production mover registry does not contain three entries");
    for (std::size_t i = 0; i < movers.size() && i < 3; ++i) {
      checks.Expect(std::string(movers[i].canonicalName) == expected_names[i],
                    "production mover ordering/name differs from D03 manifest contract");
      checks.Expect(movers[i].capabilities.requiresFieldLineAttachment,
                    "D03 mover does not require field-line attachment");
      checks.Expect(movers[i].capabilities.coefficientContract ==
                        expected_coefficients[i],
                    "D03 mover coefficient contract is inconsistent");
    }

    // Prepare two physical epochs through the same adapter used by the driver.
    // The campaign later parses final_source_state_id>=2 on every rank; this
    // preflight makes a build lacking monotonic refresh identity fail before a
    // costly site launch.  It deliberately does not claim MPI or restart
    // equivalence, which only the outer D03 campaign can establish.
    SEP::SW1DAdapter::ClearStagedInputAssignments();
    SEP::SW1DAdapter::Configure(SEP::SW1DAdapter::Scenario::Fast);
    const SEP::SW1DAdapter::Status first =
        SEP::SW1DAdapter::PrepareState(0.0);
    const std::uint64_t first_id =
        SEP::SW1DAdapter::GetDiagnostics().prepared_state_id;
    const SEP::SW1DAdapter::Status second =
        SEP::SW1DAdapter::PrepareState(1.0);
    const SEP::SW1DAdapter::Diagnostics diagnostics =
        SEP::SW1DAdapter::GetDiagnostics();
    checks.Expect(first.ok() && second.ok() && first_id == 1 &&
                      diagnostics.prepared_state_id == 2 &&
                      Close(diagnostics.prepared_epoch_seconds, 1.0),
                  "SWCME refresh did not expose monotonic state ID and epoch");
    const SEP::SW1DAdapter::ConfigurationSummary configuration =
        SEP::SW1DAdapter::GetConfigurationSummary();
    checks.Expect(!configuration.fingerprint.empty(),
                  "SWCME refresh has no canonical configuration fingerprint");
  } catch (const std::exception& exception) {
    checks.Expect(false, std::string("unexpected D03 preflight exception: ") +
                             exception.what());
  } catch (...) {
    checks.Expect(false, "unexpected non-standard D03 preflight exception");
  }

  std::string restore_detail;
  checks.Expect(RestoreDefaultSwcmeState(&restore_detail),
                "default SWCME restoration failed: " + restore_detail);
  SEP::Testing::Result result = checks.Finish(
      "D03 linked mover/refresh preflight passed; external MPI/restart campaign remains required",
      "D03 native integration preflight failed");
  result.configuration.push_back("claim=linked-native-preflight-only");
  result.configuration.push_back("required_external_gate=test-d03-native-integration");
  result.configuration.push_back("adapter_state_restored=fast@0s;strict");
  return result;
}

SEP::Testing::Descriptor MakeDescriptor(
    const char* id, const char* name, const char* group,
    const char* description, SEP::Testing::TestCallback callback) {
  SEP::Testing::Descriptor descriptor;
  descriptor.id = id;
  descriptor.name = name;
  descriptor.group = group;
  descriptor.description = description;
  descriptor.initialization = SEP::Testing::InitializationLevel::None;
  // These tests mutate the process-owned adapter and are therefore excluded
  // from native --all-tests (the routine subset).  Python --all discovers all
  // public IDs and executes each in a fresh process, which is the authoritative
  // complete-registry path requested for D01--D03.
  descriptor.runtime = SEP::Testing::RuntimeClass::Extended;
  descriptor.supportedBuildModes =
      "serial linked srcSEP/AMPS executable; isolated process";
  descriptor.seedPolicy = "deterministic; no RNG";
  descriptor.stateIsolation =
      "callback restores fast SWCME state; complete --all uses process isolation";
  descriptor.callback = callback;
  return descriptor;
}

}  // namespace

std::vector<SEP::Testing::Descriptor>
SEP::Testing::SwcmeImprovementDescriptors() {
  std::vector<Descriptor> descriptors;
  descriptors.push_back(MakeDescriptor(
      "D01", "Fail-closed SWCME background adapter", "swcme-background",
      "Verify typed query/preparation failures, explicit recovery policies, immutable SI samples, diagnostics, and transactional state.",
      RunD01FailClosedBackground));
  descriptors.push_back(MakeDescriptor(
      "D02", "Canonical SWCME configuration", "swcme-configuration",
      "Verify canonical presets, units, authority layers, provenance, validity, source controls, and deterministic fingerprints.",
      RunD02CanonicalConfiguration));
  descriptors.push_back(MakeDescriptor(
      "D03PRE", "Native integration evidence preflight", "native-integration",
      "Verify the linked three-mover catalog and monotonic SWCME refresh/fingerprint hooks required by the external MPI/restart D03 campaign.",
      RunD03NativePreflight));
  return descriptors;
}

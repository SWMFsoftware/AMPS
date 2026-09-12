#pragma once

// ============================================================================
// swcme_status.hpp
// ----------------------------------------------------------------------------
// Shared numerical/runtime status used by SWCME physics evaluators.
//
// Design rule
// -----------
// A physics routine must never turn an invalid/non-finite intermediate value
// into a plausible physical number (0 density, ambient flow, +X direction,
// etc.) merely to keep execution moving.  Instead, the first failing operation
// returns a ModelStatus that records WHAT failed and, for batch evaluators,
// WHICH sample failed.  Existing void APIs remain source-compatible wrappers;
// they convert a non-OK status into an exception rather than silently repairing
// the output.  New AMPS-facing code can use the *_checked APIs directly and
// propagate ModelStatus without exceptions.
//
// NoSurface, NoConnection, and SourceInactive are expected physical/geometric
// outcomes, not numerical failures.  They remain explicit status codes so an
// AMPS/SEP caller can distinguish "nothing to inject here" from an actual
// geometry, configuration, or shock-solver failure without relying on sentinel
// numbers or exceptions.
// ============================================================================

#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace swcme {

using ModelIdentity = std::uint64_t;
using ConfigurationDigest = std::uint64_t;

// Allocation-free FNV-1a builder used to fingerprint every field that can
// influence a prepared state.  Values are serialized explicitly in little-
// endian order instead of hashing object memory, so padding, host endianness,
// and enum storage width cannot make diagnostics vary between builds.
class ConfigurationDigestBuilder {
 public:
  // Configuration fingerprints use numerical equivalence for signed zero and
  // malformed NaN payloads.  Prepared-record integrity passes false so even a
  // bit-level cached-value change is detectable.
  explicit ConfigurationDigestBuilder(
      bool canonicalize_doubles=true) noexcept
      : canonicalize_doubles_(canonicalize_doubles) {}

  void add_byte(std::uint8_t value) noexcept {
    digest_ ^= value;
    digest_ *= 1099511628211ULL;
  }

  void add_uint64(std::uint64_t value) noexcept {
    for (unsigned shift=0; shift<64; shift+=8)
      add_byte(static_cast<std::uint8_t>((value >> shift) & 0xffU));
  }

  void add_bool(bool value) noexcept { add_byte(value ? 1U : 0U); }

  void add_double(double value) noexcept {
    // Equal numerical configurations must hash equally: normalize signed
    // zero and collapse every NaN payload to one diagnostic representation.
    // Validation still rejects NaNs; canonicalization only keeps an error
    // digest reproducible when malformed input reaches the ownership guard.
    if (canonicalize_doubles_ && value==0.0) value=0.0;
    std::uint64_t bits=0;
    if (canonicalize_doubles_ && std::isnan(value)) {
      bits=0x7ff8000000000000ULL;
    } else {
      static_assert(sizeof(bits)==sizeof(value),
                    "configuration digest requires 64-bit double");
      std::memcpy(&bits,&value,sizeof(bits));
    }
    add_uint64(bits);
  }

  void add_string(const char* value) noexcept {
    // Prefixing the length prevents ambiguous concatenations such as
    // ("ab","c") and ("a","bc") from producing the same byte stream.
    const std::size_t length=value ? std::strlen(value) : 0;
    add_uint64(static_cast<std::uint64_t>(length));
    for (std::size_t i=0; i<length; ++i)
      add_byte(static_cast<std::uint8_t>(value[i]));
  }

  ConfigurationDigest value() const noexcept { return digest_; }

 private:
  // Standard 64-bit FNV-1a offset basis.  This is a stable diagnostic digest,
  // not a cryptographic authenticator; model ownership remains independently
  // enforced by ModelIdentity.
  ConfigurationDigest digest_=14695981039346656037ULL;
  bool canonicalize_doubles_=true;
};

// Return a process-unique, nonzero identity for one logical Model instance.
// The identity is deliberately independent of parameter values: PST02 must
// reject a state prepared by a different instance even when both instances
// were constructed from byte-for-byte equivalent configurations.  A relaxed
// atomic is sufficient because uniqueness, rather than synchronization of any
// model data, is the only cross-thread property required here.
inline ModelIdentity next_model_identity() noexcept {
  static std::atomic<ModelIdentity> next{1};
  ModelIdentity identity=next.fetch_add(1,std::memory_order_relaxed);

  // Identity zero is reserved for an unprepared/default-constructed state.
  // Wrapping a 64-bit counter is practically unreachable, but skipping zero
  // makes the invariant explicit and keeps diagnostics unambiguous.
  if (identity==0) identity=next.fetch_add(1,std::memory_order_relaxed);
  return identity;
}

enum class StatusCode {
  Ok = 0,
  NoSurface,
  NoConnection,
  SourceInactive,
  InvalidConfiguration,
  NullPointer,
  NonFiniteInput,
  OutsideModelDomain,
  DegenerateVector,
  InvalidNumericalStep,
  GeometryFailure,
  ShockSolverFailure,
  NonFiniteResult,
  InvalidMesh,
  FileOpenFailure,
  FileWriteFailure,
  // Appended rather than inserted among earlier failures so the numeric values
  // of the pre-existing StatusCode entries remain source/binary-log compatible.
  StateModelMismatch,
  // Appended for the same compatibility reason as StateModelMismatch.  This
  // code identifies a state prepared by the same model before its mutable
  // configuration changed.
  StateConfigurationMismatch,
  // PST06 integrity failure.  Appending preserves numeric compatibility for
  // every status introduced before prepared-state record sealing.
  StalePreparedState,
  // OUT03 distinguishes a fully written temporary product that could not be
  // atomically installed from failures that occurred while writing its bytes.
  // Appending again preserves every pre-existing numeric status value.
  FileCommitFailure,
  // CON09 propagates a connectivity accuracy request that exceeds the
  // production scan budget through integration adapters.  This is distinct
  // from NoConnection: no physical connected/disconnected classification was
  // attempted because the requested resolution could not be achieved.
  ResolutionLimit
};

inline const char* status_code_name(StatusCode code) {
  switch (code) {
    case StatusCode::Ok: return "OK";
    case StatusCode::NoSurface: return "NO_SURFACE";
    case StatusCode::NoConnection: return "NO_CONNECTION";
    case StatusCode::SourceInactive: return "SOURCE_INACTIVE";
    case StatusCode::InvalidConfiguration: return "INVALID_CONFIGURATION";
    case StatusCode::StateModelMismatch: return "STATE_MODEL_MISMATCH";
    case StatusCode::StateConfigurationMismatch:
      return "STATE_CONFIGURATION_MISMATCH";
    case StatusCode::StalePreparedState: return "STALE_PREPARED_STATE";
    case StatusCode::NullPointer: return "NULL_POINTER";
    case StatusCode::NonFiniteInput: return "NONFINITE_INPUT";
    case StatusCode::OutsideModelDomain: return "OUTSIDE_MODEL_DOMAIN";
    case StatusCode::DegenerateVector: return "DEGENERATE_VECTOR";
    case StatusCode::InvalidNumericalStep: return "INVALID_NUMERICAL_STEP";
    case StatusCode::GeometryFailure: return "GEOMETRY_FAILURE";
    case StatusCode::ShockSolverFailure: return "SHOCK_SOLVER_FAILURE";
    case StatusCode::NonFiniteResult: return "NONFINITE_RESULT";
    case StatusCode::InvalidMesh: return "INVALID_MESH";
    case StatusCode::FileOpenFailure: return "FILE_OPEN_FAILURE";
    case StatusCode::FileWriteFailure: return "FILE_WRITE_FAILURE";
    case StatusCode::FileCommitFailure: return "FILE_COMMIT_FAILURE";
    case StatusCode::ResolutionLimit: return "RESOLUTION_LIMIT";
  }
  return "UNKNOWN_STATUS";
}

struct ModelStatus {
  static constexpr std::size_t npos = std::numeric_limits<std::size_t>::max();

  StatusCode code = StatusCode::Ok;
  std::size_t sample_index = npos;
  const char* context = "";
  double offending_value = 0.0;
  bool has_offending_value = false;
  ModelIdentity expected_model_identity = 0;
  ModelIdentity supplied_model_identity = 0;
  bool has_model_identities = false;
  ConfigurationDigest expected_configuration_digest = 0;
  ConfigurationDigest supplied_configuration_digest = 0;
  bool has_configuration_digests = false;
  ConfigurationDigest expected_state_integrity = 0;
  ConfigurationDigest computed_state_integrity = 0;
  bool has_state_integrity = false;
  // OUT02 reports where an output stream first stopped accepting data.  This
  // offset counts bytes successfully accepted by the writer backend before
  // the failed formatted write, flush, stream check, or close operation.
  std::size_t io_byte_offset = 0;
  bool has_io_byte_offset = false;

  constexpr bool ok() const noexcept { return code == StatusCode::Ok; }
  constexpr bool no_surface() const noexcept {
    return code == StatusCode::NoSurface;
  }
  constexpr bool no_connection() const noexcept {
    return code == StatusCode::NoConnection;
  }
  constexpr bool source_inactive() const noexcept {
    return code == StatusCode::SourceInactive;
  }
  constexpr bool failure() const noexcept {
    return code != StatusCode::Ok && code != StatusCode::NoSurface &&
           code != StatusCode::NoConnection && code != StatusCode::SourceInactive;
  }

  static constexpr ModelStatus success() noexcept { return {}; }

  static constexpr ModelStatus make(StatusCode c, const char* where,
                                    std::size_t index = npos) noexcept {
    ModelStatus s;
    s.code = c;
    s.context = where;
    s.sample_index = index;
    return s;
  }

  static constexpr ModelStatus make_value(StatusCode c, const char* where,
                                          double value,
                                          std::size_t index = npos) noexcept {
    ModelStatus s = make(c, where, index);
    s.offending_value = value;
    s.has_offending_value = true;
    return s;
  }

  // Build the dedicated ownership failure used by every state-consuming
  // model API.  Recording both identities makes a mixed-state failure
  // diagnosable without inspecting addresses or reproducing the calculation.
  static constexpr ModelStatus state_model_mismatch(
      const char* where, ModelIdentity expected,
      ModelIdentity supplied, ConfigurationDigest expected_configuration=0,
      ConfigurationDigest supplied_configuration=0,
      bool include_configuration=false) noexcept {
    ModelStatus s=make(StatusCode::StateModelMismatch,where);
    s.expected_model_identity=expected;
    s.supplied_model_identity=supplied;
    s.has_model_identities=true;
    s.expected_configuration_digest=expected_configuration;
    s.supplied_configuration_digest=supplied_configuration;
    s.has_configuration_digests=include_configuration;
    return s;
  }

  // Build the PST03 failure returned when the owner identity is correct but
  // the receiving model's current Params no longer match the configuration
  // from which the state was prepared.
  static constexpr ModelStatus state_configuration_mismatch(
      const char* where, ConfigurationDigest expected,
      ConfigurationDigest supplied) noexcept {
    ModelStatus s=make(StatusCode::StateConfigurationMismatch,where);
    s.expected_configuration_digest=expected;
    s.supplied_configuration_digest=supplied;
    s.has_configuration_digests=true;
    return s;
  }

  // Construct the explicit PST06 rejection.  "Expected" is the private seal
  // written by prepare_step(); "computed" is the digest of the record supplied
  // to the consumer.  Keeping both values makes corruption diagnosable without
  // exposing any API that can rewrite the private seal.
  static constexpr ModelStatus stale_prepared_state(
      const char* where, ConfigurationDigest expected,
      ConfigurationDigest computed) noexcept {
    ModelStatus s=make(StatusCode::StalePreparedState,where);
    s.expected_state_integrity=expected;
    s.computed_state_integrity=computed;
    s.has_state_integrity=true;
    return s;
  }

  // Construct the explicit OUT02 status without overloading offending_value,
  // which represents invalid physics input.  sample_index remains available
  // for a row/cell identifier while io_byte_offset locates the stream failure.
  static constexpr ModelStatus file_write_failure(
      const char* where, std::size_t byte_offset,
      std::size_t item_index=npos) noexcept {
    ModelStatus s=make(StatusCode::FileWriteFailure,where,item_index);
    s.io_byte_offset=byte_offset;
    s.has_io_byte_offset=true;
    return s;
  }

  // Construct the OUT03 status returned after the temporary output has been
  // written and closed but atomic replacement of the destination fails.  The
  // byte count records the complete staged product size and remains separate
  // from offending_value, which is reserved for invalid numerical inputs.
  static constexpr ModelStatus file_commit_failure(
      const char* where, std::size_t byte_offset) noexcept {
    ModelStatus s=make(StatusCode::FileCommitFailure,where);
    s.io_byte_offset=byte_offset;
    s.has_io_byte_offset=true;
    return s;
  }

  std::string summary() const {
    std::ostringstream out;
    out << status_code_name(code);
    if (context && context[0] != '\0') out << " in " << context;
    if (sample_index != npos) out << " at sample " << sample_index;
    if (has_offending_value) out << " (value=" << offending_value << ')';
    if (has_model_identities) {
      out << " (expected_model_identity=" << expected_model_identity
          << ", supplied_model_identity=" << supplied_model_identity << ')';
    }
    if (has_configuration_digests) {
      out << " (expected_configuration_digest=0x" << std::hex
          << std::setw(16) << std::setfill('0')
          << expected_configuration_digest
          << ", supplied_configuration_digest=0x" << std::setw(16)
          << supplied_configuration_digest << std::dec << ')';
    }
    if (has_state_integrity) {
      out << " (expected_state_integrity=0x" << std::hex
          << std::setw(16) << std::setfill('0') << expected_state_integrity
          << ", computed_state_integrity=0x" << std::setw(16)
          << computed_state_integrity << std::dec << ')';
    }
    if (has_io_byte_offset)
      out << " (io_byte_offset=" << io_byte_offset << ')';
    return out.str();
  }
};

// Source-compatible public wrappers use this helper so an error cannot be
// ignored merely because an older call site uses a void evaluator.  Expected
// NO_SURFACE results are handled by geometry APIs before reaching this helper.
inline void throw_if_error(const ModelStatus& status) {
  if (status.failure()) throw std::runtime_error(status.summary());
}

}  // namespace swcme

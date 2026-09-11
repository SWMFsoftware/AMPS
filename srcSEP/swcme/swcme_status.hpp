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
#include <cstddef>
#include <cstdint>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace swcme {

using ModelIdentity = std::uint64_t;

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
  StateModelMismatch
};

inline const char* status_code_name(StatusCode code) {
  switch (code) {
    case StatusCode::Ok: return "OK";
    case StatusCode::NoSurface: return "NO_SURFACE";
    case StatusCode::NoConnection: return "NO_CONNECTION";
    case StatusCode::SourceInactive: return "SOURCE_INACTIVE";
    case StatusCode::InvalidConfiguration: return "INVALID_CONFIGURATION";
    case StatusCode::StateModelMismatch: return "STATE_MODEL_MISMATCH";
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
      ModelIdentity supplied) noexcept {
    ModelStatus s=make(StatusCode::StateModelMismatch,where);
    s.expected_model_identity=expected;
    s.supplied_model_identity=supplied;
    s.has_model_identities=true;
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

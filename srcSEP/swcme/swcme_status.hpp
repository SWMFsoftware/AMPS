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
// NoSurface is an expected geometric outcome, not a numerical failure.  It is
// therefore represented explicitly and can be distinguished from GeometryError
// or ShockSolverFailure by callers that query finite shock surfaces.
// ============================================================================

#include <cstddef>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace swcme {

enum class StatusCode {
  Ok = 0,
  NoSurface,
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
  FileWriteFailure
};

inline const char* status_code_name(StatusCode code) {
  switch (code) {
    case StatusCode::Ok: return "OK";
    case StatusCode::NoSurface: return "NO_SURFACE";
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

  constexpr bool ok() const noexcept { return code == StatusCode::Ok; }
  constexpr bool no_surface() const noexcept {
    return code == StatusCode::NoSurface;
  }
  constexpr bool failure() const noexcept {
    return code != StatusCode::Ok && code != StatusCode::NoSurface;
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

  std::string summary() const {
    std::ostringstream out;
    out << status_code_name(code);
    if (context && context[0] != '\0') out << " in " << context;
    if (sample_index != npos) out << " at sample " << sample_index;
    if (has_offending_value) out << " (value=" << offending_value << ')';
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

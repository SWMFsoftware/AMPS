#ifndef SEP_COMMON_SEP_STATUS_H
#define SEP_COMMON_SEP_STATUS_H

#include <string>
#include <utility>

// This header is the lowest-level error contract shared by the SEP models.
// It intentionally depends only on the C++ standard library: physics kernels
// must be testable without AMPS, MPI, a mesh, or an application fatal-error
// handler.  Application adapters may translate a non-OK status at their
// boundary, but a shared kernel never aborts a process or invents a fallback.
namespace SEP {
namespace Core {

enum class StatusCode {
  Ok,
  InvalidConfiguration,
  UnsupportedCapability,
  InvalidState,
  OutOfDomain,
  DataIntegrityFailure,
  NumericalFailure,
  NotImplemented
};

struct Status {
  StatusCode code = StatusCode::Ok;
  std::string message;

  bool ok() const noexcept { return code == StatusCode::Ok; }

  static Status Success() { return {}; }

  static Status Failure(StatusCode failureCode, std::string diagnostic) {
    Status status;
    status.code = failureCode;
    status.message = std::move(diagnostic);
    return status;
  }
};

// Result<T> prevents a failed preparation step from publishing a partially
// initialized physical object.  A value is meaningful only when status.ok().
template <typename T>
struct Result {
  Status status;
  T value{};

  bool ok() const noexcept { return status.ok(); }

  static Result Success(T resolvedValue) {
    Result result;
    result.value = std::move(resolvedValue);
    return result;
  }

  static Result Failure(StatusCode code, std::string diagnostic) {
    Result result;
    result.status = Status::Failure(code, std::move(diagnostic));
    return result;
  }
};

}  // namespace Core
}  // namespace SEP

#endif  // SEP_COMMON_SEP_STATUS_H

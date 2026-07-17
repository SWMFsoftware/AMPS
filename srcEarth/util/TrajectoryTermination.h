#ifndef _SRC_EARTH_UTIL_TRAJECTORY_TERMINATION_H_
#define _SRC_EARTH_UTIL_TRAJECTORY_TERMINATION_H_

namespace Earth {
namespace GridlessMode {

// Explicit result of one backward trajectory.  The first three values are
// resolved physical classifications; all remaining values are numerical/unresolved
// outcomes that must be reported, retried, or excluded from physical denominators.
enum class TrajectoryTermination {
  OuterBoundaryAllowed = 0,
  InnerBoundaryForbidden,
  MagneticallyTrappedForbidden,
  TimeLimit,
  StepLimit,
  DistanceLimit,
  InvalidTimeStep,
  InvalidField,
  NumericalFailure,
  Count
};

inline const char* TrajectoryTerminationName(TrajectoryTermination t) {
  switch (t) {
    case TrajectoryTermination::OuterBoundaryAllowed: return "OUTER_BOUNDARY_ALLOWED";
    case TrajectoryTermination::InnerBoundaryForbidden: return "INNER_BOUNDARY_FORBIDDEN";
    case TrajectoryTermination::MagneticallyTrappedForbidden: return "MAGNETICALLY_TRAPPED_FORBIDDEN";
    case TrajectoryTermination::TimeLimit: return "TIME_LIMIT";
    case TrajectoryTermination::StepLimit: return "STEP_LIMIT";
    case TrajectoryTermination::DistanceLimit: return "DISTANCE_LIMIT";
    case TrajectoryTermination::InvalidTimeStep: return "INVALID_TIME_STEP";
    case TrajectoryTermination::InvalidField: return "INVALID_FIELD";
    case TrajectoryTermination::NumericalFailure: return "NUMERICAL_FAILURE";
    default: return "UNKNOWN";
  }
}

inline bool IsResolvedTermination(TrajectoryTermination t) {
  return t==TrajectoryTermination::OuterBoundaryAllowed ||
         t==TrajectoryTermination::InnerBoundaryForbidden ||
         t==TrajectoryTermination::MagneticallyTrappedForbidden;
}

inline bool IsAllowedTermination(TrajectoryTermination t) {
  return t==TrajectoryTermination::OuterBoundaryAllowed;
}

// Numerical safety limits are intentionally not part of IsResolvedTermination():
// density/transmission calculations such as F3 must continue to report and exclude
// them as unresolved samples.  Legacy Boolean cutoff searches, however, have no
// third state and historically interpreted "no escape before a configured limit" as
// FORBIDDEN.  Keep that caller-specific policy explicit rather than changing the
// structured trajectory result.
inline bool IsTraceLimitTermination(TrajectoryTermination t) {
  return t==TrajectoryTermination::TimeLimit ||
         t==TrajectoryTermination::StepLimit ||
         t==TrajectoryTermination::DistanceLimit;
}

inline bool IsPhysicalForbiddenTermination(TrajectoryTermination t) {
  return t==TrajectoryTermination::InnerBoundaryForbidden ||
         t==TrajectoryTermination::MagneticallyTrappedForbidden;
}

inline bool IsCutoffForbiddenTermination(TrajectoryTermination t) {
  return IsPhysicalForbiddenTermination(t) || IsTraceLimitTermination(t);
}

inline bool IsRetryableNumericalTermination(TrajectoryTermination t) {
  return t==TrajectoryTermination::InvalidTimeStep ||
         t==TrajectoryTermination::InvalidField ||
         t==TrajectoryTermination::NumericalFailure;
}

} // namespace GridlessMode
} // namespace Earth

#endif

#ifndef _SRC_EARTH_UTIL_TRAJECTORY_CONTRACT_H_
#define _SRC_EARTH_UTIL_TRAJECTORY_CONTRACT_H_

//======================================================================================
// TrajectoryContract.h
//======================================================================================
// Backend-neutral request/result contract for SEP-in-geospace characteristics.
//
// Why this contract exists
// ------------------------
// Cutoff, directional access, flux, and energy-spectrum products all depend on the
// same physical question: what happens to one phase-space seed when it is traced
// through one immutable field snapshot?  Before Roadmap Step 4, several pieces of
// that question were ambient process state (the mover), implicit convention (how
// "backward" time was represented), or backend-specific output (the boundary state).
// That made it possible for gridless and Mode3D/SWMF callers to ask subtly different
// questions while using nominally identical inputs.
//
// This header makes those choices explicit and shared.  It intentionally depends on
// neither AMPS nor MPI, Geopack, SPICE, or the mesh, so its validation and retry rules
// can be tested as ordinary C++11 code.  Production adapters may retain legacy wrapper
// functions, but new trajectory work enters and leaves through Request and Result.
//
// Released physics scope
// ----------------------
// Production Step-4 traces use a frozen magnetic snapshot.  StaticMagneticAntiparticle
// is the normal time-reversal construction; StaticMagneticSameCharge exists only for
// reproducibility of historical AMPS cutoff products.  PhysicalBackwardTime reserves
// the semantics required by a future E(x,t),B(x,t) integrator.  Merely selecting that
// enum does not release such an integrator: every current production backend rejects
// electric or explicitly time-dependent requests before taking a trajectory step.
//
// Unit contract
// -------------
//   position          m
//   momentum          kg m/s
//   rigidity          GV
//   charge            C
//   rest mass         kg
//   trace time        s
//   trace distance    m
//======================================================================================

#include "TrajectoryTermination.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>

namespace Earth {
namespace Trajectory {

// Mover identity belongs to a trajectory request, not only to a process-global parser
// setting.  GridlessParticleMovers.h aliases its historical MoverType name to this enum
// so old call sites remain source-compatible while request-based calls are re-entrant.
enum class Mover {
  BORIS,
  HC4,
  RK2,
  RK4,
  RK6,
  GC2,
  GC4,
  GC6,
  HYBRID
};

inline bool IsKnownMover(Mover mover) {
  switch (mover) {
    case Mover::BORIS:
    case Mover::HC4:
    case Mover::RK2:
    case Mover::RK4:
    case Mover::RK6:
    case Mover::GC2:
    case Mover::GC4:
    case Mover::GC6:
    case Mover::HYBRID:
      return true;
  }
  return false;
}

inline const char* MoverName(Mover mover) {
  switch (mover) {
    case Mover::BORIS:  return "BORIS";
    case Mover::HC4:    return "HC4";
    case Mover::RK2:    return "RK2";
    case Mover::RK4:    return "RK4";
    case Mover::RK6:    return "RK6";
    case Mover::GC2:    return "GC2";
    case Mover::GC4:    return "GC4";
    case Mover::GC6:    return "GC6";
    case Mover::HYBRID: return "HYBRID";
  }
  return "UNKNOWN";
}

inline bool IsFullOrbitMover(Mover mover) {
  return mover==Mover::BORIS || mover==Mover::HC4 || mover==Mover::RK2 ||
         mover==Mover::RK4 || mover==Mover::RK6;
}

inline bool IsReducedOrbitMover(Mover mover) {
  // HYBRID can enter a guiding-centre branch and therefore needs the same explicit
  // validity declaration as a pure GC mover.
  return mover==Mover::GC2 || mover==Mover::GC4 || mover==Mover::GC6 ||
         mover==Mover::HYBRID;
}

enum class BackwardTimeMode {
  // Historical compatibility: reverse the launch direction while retaining the
  // physical charge and integrate with positive numerical time.
  StaticMagneticSameCharge,

  // Static-field antiparticle construction: reverse the launch direction and charge,
  // then integrate with positive numerical time.  It is valid only for frozen B and
  // E=0, and is the released physical backtrace convention.
  StaticMagneticAntiparticle,

  // Same-charge integration backward in physical time.  This is the required public
  // convention for a future dynamic electromagnetic characteristic.  Current
  // production adapters reject requests that need that unreleased capability.
  PhysicalBackwardTime
};

inline bool IsKnownBackwardTimeMode(BackwardTimeMode mode) {
  switch (mode) {
    case BackwardTimeMode::StaticMagneticSameCharge:
    case BackwardTimeMode::StaticMagneticAntiparticle:
    case BackwardTimeMode::PhysicalBackwardTime:
      return true;
  }
  return false;
}

inline const char* BackwardTimeModeName(BackwardTimeMode mode) {
  switch (mode) {
    case BackwardTimeMode::StaticMagneticSameCharge:
      return "STATIC_MAGNETIC_SAME_CHARGE";
    case BackwardTimeMode::StaticMagneticAntiparticle:
      return "STATIC_MAGNETIC_ANTIPARTICLE";
    case BackwardTimeMode::PhysicalBackwardTime:
      return "PHYSICAL_BACKWARD_TIME";
  }
  return "UNKNOWN";
}

enum class RequestStatus {
  Valid,
  NonFiniteState,
  InvalidDirection,
  InvalidSpecies,
  InvalidRigidity,
  InvalidBudget,
  InvalidMover,
  InvalidBackwardTimeMode,
  MissingSnapshotIdentity,
  DynamicFieldNeedsPhysicalBackwardTime,
  ReducedOrbitNeedsStaticMagneticField,
  ReducedOrbitValidityNotDeclared,
  ElectromagneticMoverNotImplemented
};

inline const char* RequestStatusName(RequestStatus status) {
  switch (status) {
    case RequestStatus::Valid: return "VALID";
    case RequestStatus::NonFiniteState: return "NON_FINITE_STATE";
    case RequestStatus::InvalidDirection: return "INVALID_DIRECTION";
    case RequestStatus::InvalidSpecies: return "INVALID_SPECIES";
    case RequestStatus::InvalidRigidity: return "INVALID_RIGIDITY";
    case RequestStatus::InvalidBudget: return "INVALID_BUDGET";
    case RequestStatus::InvalidMover: return "INVALID_MOVER";
    case RequestStatus::InvalidBackwardTimeMode: return "INVALID_BACKWARD_TIME_MODE";
    case RequestStatus::MissingSnapshotIdentity: return "MISSING_SNAPSHOT_IDENTITY";
    case RequestStatus::DynamicFieldNeedsPhysicalBackwardTime:
      return "DYNAMIC_FIELD_NEEDS_PHYSICAL_BACKWARD_TIME";
    case RequestStatus::ReducedOrbitNeedsStaticMagneticField:
      return "REDUCED_ORBIT_NEEDS_STATIC_MAGNETIC_FIELD";
    case RequestStatus::ReducedOrbitValidityNotDeclared:
      return "REDUCED_ORBIT_VALIDITY_NOT_DECLARED";
    case RequestStatus::ElectromagneticMoverNotImplemented:
      return "ELECTROMAGNETIC_MOVER_NOT_IMPLEMENTED";
  }
  return "UNKNOWN";
}

// Stable, non-cryptographic digest used to carry FieldProvider snapshot identity
// through the compact trajectory result.  Empty means "no identity" and maps to zero;
// zero is otherwise reserved.  Backends must compare a nonzero request fingerprint
// with the active snapshot rather than simply copying caller input into the result.
inline std::uint64_t SnapshotFingerprint(const std::string& snapshotId) {
  if (snapshotId.empty()) return 0;
  std::uint64_t hash=14695981039346656037ULL;
  for (std::string::const_iterator it=snapshotId.begin();it!=snapshotId.end();++it) {
    hash^=static_cast<unsigned char>(*it);
    hash*=1099511628211ULL;
  }
  return hash==0 ? 1 : hash;
}

// A nonzero fingerprint in a request is always an assertion, even when identity is
// not marked mandatory.  When identity is mandatory, both sides must be nonzero.
// Keeping this logic in the common contract gives the direct and mesh adapters exactly
// the same stale/mismatched-snapshot behavior.
inline bool SnapshotIdentityMatches(std::uint64_t requestedFingerprint,
                                    bool identityRequired,
                                    std::uint64_t activeFingerprint) {
  if (identityRequired &&
      (requestedFingerprint==0 || activeFingerprint==0)) return false;
  if (requestedFingerprint!=0)
    return activeFingerprint!=0 && requestedFingerprint==activeFingerprint;
  return true;
}

struct Request {
  double x0_m[3]{0.0,0.0,0.0};
  double direction0_unit[3]{1.0,0.0,0.0};
  double rigidity_GV{0.0};
  double charge_C{0.0};
  double restMass_kg{0.0};

  Mover mover{Mover::BORIS};
  BackwardTimeMode backwardTimeMode{BackwardTimeMode::StaticMagneticAntiparticle};

  double maxTraceTime_s{0.0};
  double maxTraceDistance_m{0.0}; // <=0 deliberately disables only the path cap
  int maxSteps{0};
  bool captureExitState{false};

  // These flags describe fields actually used by the characteristic.  A coupled SWMF
  // snapshot may contain an E array while a released magnetic-only calculation leaves
  // electricFieldEnabled=false.
  bool electricFieldEnabled{false};
  bool fieldTimeDependent{false};
  bool electromagneticMoverImplemented{false};
  bool reducedOrbitValidityDeclared{false};

  std::uint64_t snapshotFingerprint{0};
  bool requireSnapshotIdentity{false};
};

inline RequestStatus ValidateRequest(const Request& request) {
  for (int d=0;d<3;++d) {
    if (!std::isfinite(request.x0_m[d]) ||
        !std::isfinite(request.direction0_unit[d]))
      return RequestStatus::NonFiniteState;
  }

  const double n2=request.direction0_unit[0]*request.direction0_unit[0]+
                  request.direction0_unit[1]*request.direction0_unit[1]+
                  request.direction0_unit[2]*request.direction0_unit[2];
  // The request is a contract, so silently normalizing a malformed direction would
  // hide a caller error and make cross-backend reproduction harder.
  if (!std::isfinite(n2) || std::fabs(n2-1.0)>1.0e-8)
    return RequestStatus::InvalidDirection;
  if (!(request.restMass_kg>0.0) || !std::isfinite(request.restMass_kg) ||
      request.charge_C==0.0 || !std::isfinite(request.charge_C))
    return RequestStatus::InvalidSpecies;
  if (!(request.rigidity_GV>0.0) || !std::isfinite(request.rigidity_GV))
    return RequestStatus::InvalidRigidity;
  if (!(request.maxTraceTime_s>0.0) || request.maxSteps<=0 ||
      !std::isfinite(request.maxTraceTime_s) ||
      !std::isfinite(request.maxTraceDistance_m))
    return RequestStatus::InvalidBudget;
  if (!IsKnownMover(request.mover)) return RequestStatus::InvalidMover;
  if (!IsKnownBackwardTimeMode(request.backwardTimeMode))
    return RequestStatus::InvalidBackwardTimeMode;
  if (request.requireSnapshotIdentity && request.snapshotFingerprint==0)
    return RequestStatus::MissingSnapshotIdentity;

  const bool electromagnetic=request.electricFieldEnabled || request.fieldTimeDependent;
  if (electromagnetic &&
      request.backwardTimeMode!=BackwardTimeMode::PhysicalBackwardTime)
    return RequestStatus::DynamicFieldNeedsPhysicalBackwardTime;
  if (electromagnetic && IsReducedOrbitMover(request.mover))
    return RequestStatus::ReducedOrbitNeedsStaticMagneticField;
  if (IsReducedOrbitMover(request.mover) && !request.reducedOrbitValidityDeclared)
    return RequestStatus::ReducedOrbitValidityNotDeclared;
  if (electromagnetic && !request.electromagneticMoverImplemented)
    return RequestStatus::ElectromagneticMoverNotImplemented;
  return RequestStatus::Valid;
}

// Capability gate used by every currently released production adapter.  It is kept
// separate from ValidateRequest because the common schema deliberately describes a
// future physical-backward electromagnetic request, while Step 4 implements only
// frozen magnetic characteristics.  A caller cannot enable unreleased physics merely
// by setting electromagneticMoverImplemented in its own request record.
inline bool IsReleasedFrozenMagneticRequest(const Request& request) {
  return !request.electricFieldEnabled && !request.fieldTimeDependent &&
         request.backwardTimeMode!=BackwardTimeMode::PhysicalBackwardTime;
}

// Complete phase-space record at the first outer-boundary event.  Position and
// momentum are evaluated at the same event fraction.  `valid` becomes true only after
// a finite boundary magnetic field has also supplied cosAlpha; every non-allowed
// result and every non-capturing request leaves it false.
struct ExitState {
  double x_exit_m[3]{0.0,0.0,0.0};
  double p_exit_SI[3]{0.0,0.0,0.0};
  double v_exit_unit[3]{0.0,0.0,0.0};
  double cosAlpha{0.0};
  double traceTimeAtExit_s{0.0};
  double rigidityAtExit_GV{0.0};
  bool valid{false};
};

// Populate the backend-independent kinematic part of ExitState.  The event position
// is supplied by the shared boundary locator; momentum is linearly interpolated over
// the same accepted numerical step.  This helper is used by both production backends
// and directly tested against an analytic interpolation reference.
inline bool PopulateExitKinematics(ExitState& state,
                                   const double eventPosition_m[3],
                                   const double momentumBefore_SI[3],
                                   const double momentumAfter_SI[3],
                                   double eventFraction,
                                   double traceTimeAtStepEnd_s,
                                   double acceptedStep_s,
                                   double absoluteCharge_C) {
  state=ExitState();
  if (eventPosition_m==nullptr || momentumBefore_SI==nullptr ||
      momentumAfter_SI==nullptr || !std::isfinite(eventFraction) ||
      eventFraction<0.0 || eventFraction>1.0 ||
      !std::isfinite(traceTimeAtStepEnd_s) || traceTimeAtStepEnd_s<0.0 ||
      !std::isfinite(acceptedStep_s) || acceptedStep_s<0.0 ||
      !std::isfinite(absoluteCharge_C) || !(absoluteCharge_C>0.0)) return false;
  const double timeTolerance=16.0*std::numeric_limits<double>::epsilon()*
      std::max(1.0,traceTimeAtStepEnd_s);
  if (acceptedStep_s>traceTimeAtStepEnd_s+timeTolerance) return false;

  double p2=0.0;
  for (int d=0;d<3;++d) {
    if (!std::isfinite(eventPosition_m[d]) ||
        !std::isfinite(momentumBefore_SI[d]) ||
        !std::isfinite(momentumAfter_SI[d])) return false;
    state.x_exit_m[d]=eventPosition_m[d];
    state.p_exit_SI[d]=momentumBefore_SI[d]+
        eventFraction*(momentumAfter_SI[d]-momentumBefore_SI[d]);
    p2+=state.p_exit_SI[d]*state.p_exit_SI[d];
  }
  if (!(p2>0.0) || !std::isfinite(p2)) return false;

  const double p=std::sqrt(p2);
  for (int d=0;d<3;++d) state.v_exit_unit[d]=state.p_exit_SI[d]/p;
  const double stepStart=traceTimeAtStepEnd_s-acceptedStep_s;
  state.traceTimeAtExit_s=std::max(0.0,stepStart+eventFraction*acceptedStep_s);
  state.rigidityAtExit_GV=(p*299792458.0/absoluteCharge_C)*1.0e-9;
  return std::isfinite(state.traceTimeAtExit_s) &&
         std::isfinite(state.rigidityAtExit_GV);
}

// Finish the exit record using B evaluated just inside the valid field domain.  The
// pitch cosine is clamped only to its exact mathematical range to remove roundoff-sized
// overshoot; a zero or non-finite field invalidates the entire requested exit record.
inline bool CompleteExitPitchAngle(ExitState& state,
                                   const double magneticField_T[3]) {
  state.valid=false;
  if (magneticField_T==nullptr) return false;
  double b2=0.0;
  double v2=0.0;
  for (int d=0;d<3;++d) {
    if (!std::isfinite(magneticField_T[d]) ||
        !std::isfinite(state.v_exit_unit[d])) return false;
    b2+=magneticField_T[d]*magneticField_T[d];
    v2+=state.v_exit_unit[d]*state.v_exit_unit[d];
  }
  if (!(b2>0.0) || !std::isfinite(b2) || !std::isfinite(v2) ||
      std::fabs(v2-1.0)>1.0e-10) return false;
  const double invB=1.0/std::sqrt(b2);
  double cosine=0.0;
  for (int d=0;d<3;++d)
    cosine+=state.v_exit_unit[d]*magneticField_T[d]*invB;
  if (!std::isfinite(cosine)) return false;
  state.cosAlpha=std::max(-1.0,std::min(1.0,cosine));
  state.valid=true;
  return true;
}

struct Result {
  GridlessMode::TrajectoryTermination termination{
      GridlessMode::TrajectoryTermination::NumericalFailure};
  ExitState exitState{};
  double traceTime_s{0.0};
  double traceDistance_m{0.0};
  int steps{0};
  int retryCount{0};

  // The primary fields preserve the normal-budget classification before any optional
  // time-convergence extension.  Extensions restart at the original phase-space seed;
  // they are distinct from numerical-failure retries and are counted separately.
  int primaryTerminationCode{-1};
  double primaryTraceTime_s{0.0};
  int traceExtensionCount{0};
  double initialTraceLimit_s{0.0};
  double finalTraceLimit_s{0.0};

  int mirrorPoints{0};
  int bounceCycles{0};
  int driftRevolutions{0};
  double driftAngle_rad{0.0};
  double driftMeanRadiusChange_m{0.0};
  int trapMechanism{0};
  double momentumRelativeSpread{0.0};

  Mover mover{Mover::BORIS};
  BackwardTimeMode backwardTimeMode{BackwardTimeMode::StaticMagneticAntiparticle};
  std::uint64_t snapshotFingerprint{0};

  bool allowed() const {
    return GridlessMode::IsAllowedTermination(termination);
  }
  bool resolved() const {
    return GridlessMode::IsResolvedTermination(termination);
  }
};

// One policy implementation is shared by gridless and mesh-backed traces.  A genuine
// numerical failure may receive a bounded smaller-step retry.  A valid TIME_LIMIT or
// STEP_LIMIT is never routed through that recovery path; optional convergence
// extensions are separate and preserve their own provenance in Result.
struct RetryPolicy {
  int maximumNumericalRetries{1};
  double numericalDtScale{0.5};
  double numericalTimeScale{2.0};
  double numericalStepScale{2.0};
  int unresolvedExtensionPasses{0};
  double unresolvedExtensionFactor{2.0};
};

inline bool IsValidRetryPolicy(const RetryPolicy& policy) {
  return policy.maximumNumericalRetries>=0 &&
         policy.unresolvedExtensionPasses>=0 &&
         policy.numericalDtScale>0.0 && policy.numericalDtScale<1.0 &&
         policy.numericalTimeScale>=1.0 && policy.numericalStepScale>=1.0 &&
         policy.unresolvedExtensionFactor>1.0 &&
         std::isfinite(policy.numericalDtScale) &&
         std::isfinite(policy.numericalTimeScale) &&
         std::isfinite(policy.numericalStepScale) &&
         std::isfinite(policy.unresolvedExtensionFactor);
}

inline bool ShouldRetryNumerical(const Result& result,int retriesUsed,
                                 const RetryPolicy& policy) {
  return IsValidRetryPolicy(policy) && retriesUsed>=0 &&
         retriesUsed<policy.maximumNumericalRetries &&
         GridlessMode::IsRetryableNumericalTermination(result.termination);
}

inline bool IsUnresolvedExtensionCandidate(
    GridlessMode::TrajectoryTermination termination) {
  // DISTANCE_LIMIT is a separately declared path budget and must never be relaxed by a
  // trace-time convergence study.
  return termination==GridlessMode::TrajectoryTermination::TimeLimit ||
         termination==GridlessMode::TrajectoryTermination::StepLimit;
}

inline bool ShouldExtendUnresolved(const Result& result,int extensionsUsed,
                                   const RetryPolicy& policy) {
  return IsValidRetryPolicy(policy) && extensionsUsed>=0 &&
         extensionsUsed<policy.unresolvedExtensionPasses &&
         IsUnresolvedExtensionCandidate(result.termination);
}

inline double ExtensionTimeBudget(double primaryTime_s,int extensionNumber,
                                  const RetryPolicy& policy) {
  if (!IsValidRetryPolicy(policy) || !(primaryTime_s>0.0) ||
      !std::isfinite(primaryTime_s) || extensionNumber<1)
    return std::numeric_limits<double>::quiet_NaN();
  const double value=primaryTime_s*std::pow(
      policy.unresolvedExtensionFactor,static_cast<double>(extensionNumber));
  return std::isfinite(value) ? value : std::numeric_limits<double>::quiet_NaN();
}

inline int ScaledStepBudget(int primarySteps,double timeFactor,
                            const Result& previous,double targetTime_s) {
  if (primarySteps<=0 || !(timeFactor>0.0) || !std::isfinite(timeFactor) ||
      !(targetTime_s>0.0) || !std::isfinite(targetTime_s)) return 0;
  long double desired=primarySteps;
  desired=std::max(desired,
      std::ceil(static_cast<long double>(primarySteps)*timeFactor*1.25L));
  if (previous.steps>0 && previous.traceTime_s>0.0 &&
      std::isfinite(previous.traceTime_s)) {
    const long double meanDt=static_cast<long double>(previous.traceTime_s)/
                             static_cast<long double>(previous.steps);
    if (meanDt>0.0L && std::isfinite(static_cast<double>(meanDt)))
      desired=std::max(desired,
          std::ceil(static_cast<long double>(targetTime_s)/meanDt*1.25L));
  }
  const long double maxInt=static_cast<long double>(
      std::numeric_limits<int>::max());
  return static_cast<int>(std::max(1.0L,std::min(maxInt,desired)));
}

// Exact dp/dt=qE update used to test forward/backward convention algebra.  Position
// integration remains mover-specific and no current production backend calls this
// helper to imply that electromagnetic tracing has been released.
inline void AdvanceUniformElectricMomentum(double p_SI[3],double charge_C,
                                            const double E_V_m[3],double dt_s) {
  for (int d=0;d<3;++d) p_SI[d]+=charge_C*E_V_m[d]*dt_s;
}

inline double MomentumMagnitude(const double p_SI[3]) {
  return std::sqrt(p_SI[0]*p_SI[0]+p_SI[1]*p_SI[1]+p_SI[2]*p_SI[2]);
}

// For differential intensity J=p^2 f, Liouville invariance of f gives the factor
// (p_local/p_boundary)^2.  Static-B traces normally return one, while the explicit
// expression is required by a future electric-field characteristic.
inline double PhaseSpaceIntensityFactor(const double pLocal_SI[3],
                                        const double pBoundary_SI[3]) {
  const double pl=MomentumMagnitude(pLocal_SI);
  const double pb=MomentumMagnitude(pBoundary_SI);
  if (!(pb>0.0) || !std::isfinite(pl) || !std::isfinite(pb))
    return std::numeric_limits<double>::quiet_NaN();
  const double ratio=pl/pb;
  return ratio*ratio;
}

} // namespace Trajectory
} // namespace Earth

// Compatibility aliases keep existing gridless and Mode3D signatures source-compatible
// while ownership of the records moves to the common backend-neutral contract.
namespace Earth {
namespace GridlessMode {
using TrajectoryRequest=Earth::Trajectory::Request;
using TrajectoryExitState=Earth::Trajectory::ExitState;
using TrajectoryResult=Earth::Trajectory::Result;
} // namespace GridlessMode
} // namespace Earth

#endif // _SRC_EARTH_UTIL_TRAJECTORY_CONTRACT_H_

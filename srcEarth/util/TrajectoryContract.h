#ifndef _SRC_EARTH_UTIL_TRAJECTORY_CONTRACT_H_
#define _SRC_EARTH_UTIL_TRAJECTORY_CONTRACT_H_

//======================================================================================
// TrajectoryContract.h
//======================================================================================
// Backend-neutral request/result and retry-policy contract for SEP-in-geospace
// backward characteristics.
//
// This header is deliberately independent of AMPS, MPI, Geopack, and the mesh.  The
// gridless, compact-Mode3D, and SWMF-snapshot adapters all receive the same request and
// return the same result.  The contract prevents three historically dangerous forms of
// implicit state:
//
//   * the selected mover cannot be inferred only from a process-global variable;
//   * a trace limit cannot be collapsed into a physical inner-boundary loss; and
//   * the static-B reversed-trajectory convention cannot be reused silently when an
//     electric or explicitly time-dependent field is enabled.
//
// Phase-1 production traces are frozen-snapshot, magnetic-only characteristics.  The
// PhysicalBackwardTime enum value is nevertheless part of the public record now so a
// later dynamic-field implementation has an unambiguous compatibility gate rather than
// inheriting a static-field shortcut by accident.
//======================================================================================

#include "TrajectoryTermination.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace Earth {
namespace Trajectory {

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

inline const char* MoverName(Mover mover) {
  switch (mover) {
    case Mover::BORIS: return "BORIS";
    case Mover::HC4: return "HC4";
    case Mover::RK2: return "RK2";
    case Mover::RK4: return "RK4";
    case Mover::RK6: return "RK6";
    case Mover::GC2: return "GC2";
    case Mover::GC4: return "GC4";
    case Mover::GC6: return "GC6";
    case Mover::HYBRID: return "HYBRID";
  }
  return "UNKNOWN";
}

inline bool IsFullOrbitMover(Mover mover) {
  return mover==Mover::BORIS || mover==Mover::HC4 || mover==Mover::RK2 ||
         mover==Mover::RK4 || mover==Mover::RK6;
}

inline bool IsReducedOrbitMover(Mover mover) {
  return mover==Mover::GC2 || mover==Mover::GC4 || mover==Mover::GC6 ||
         mover==Mover::HYBRID;
}

enum class BackwardTimeMode {
  // Historical AMPS convention: reverse the launch velocity while retaining the
  // physical charge.  It is preserved only for explicit regression reproduction.
  StaticMagneticSameCharge,

  // Static-B antiparticle construction: reverse launch velocity and charge while
  // integrating with a positive numerical step.  This is valid only when E=0 and the
  // frozen field has no explicit time dependence.
  StaticMagneticAntiparticle,

  // Integrate the same-charge characteristic backward in physical time.  This is the
  // required convention for a future released E(x,t),B(x,t) solver.
  PhysicalBackwardTime
};

enum class RequestStatus {
  Valid,
  NonFiniteState,
  InvalidDirection,
  InvalidSpecies,
  InvalidRigidity,
  InvalidBudget,
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

struct Request {
  double x0_m[3]{0.0,0.0,0.0};
  double direction0_unit[3]{1.0,0.0,0.0};
  double rigidity_GV{0.0};
  double charge_C{0.0};
  double restMass_kg{0.0};

  Mover mover{Mover::BORIS};
  BackwardTimeMode backwardTimeMode{BackwardTimeMode::StaticMagneticSameCharge};

  double maxTraceTime_s{0.0};
  double maxTraceDistance_m{0.0}; // <=0 disables only the cumulative-distance cap
  int maxSteps{0};
  bool captureExitState{false};

  // These flags describe fields actively used by the characteristic, not merely
  // arrays available in the snapshot.  A coupled snapshot may contain derived E while
  // the released calculation remains explicitly magnetic-only.
  bool electricFieldEnabled{false};
  bool fieldTimeDependent{false};
  bool electromagneticMoverImplemented{false};
  bool reducedOrbitValidityDeclared{false};

  // FNV/hash or other stable 64-bit digest of the immutable FieldSnapshot ID.  Zero is
  // reserved for legacy wrappers that have not yet attached provenance.
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

// Complete asymptotic state at the first outer-boundary event.  Position and momentum
// are interpolated to the same event fraction; direction and pitch angle are derived
// from that state.  `valid` is false for every non-allowed termination and whenever a
// requested boundary field could not be evaluated.
struct ExitState {
  double x_exit_m[3]{0.0,0.0,0.0};
  double p_exit_SI[3]{0.0,0.0,0.0};
  double v_exit_unit[3]{0.0,0.0,0.0};
  double cosAlpha{0.0};
  double traceTimeAtExit_s{0.0};
  double rigidityAtExit_GV{0.0};
  bool valid{false};
};

struct Result {
  GridlessMode::TrajectoryTermination termination{
      GridlessMode::TrajectoryTermination::NumericalFailure};
  ExitState exitState{};
  double traceTime_s{0.0};
  double traceDistance_m{0.0};
  int steps{0};
  int retryCount{0};

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
  BackwardTimeMode backwardTimeMode{BackwardTimeMode::StaticMagneticSameCharge};
  std::uint64_t snapshotFingerprint{0};

  bool allowed() const {
    return GridlessMode::IsAllowedTermination(termination);
  }
  bool resolved() const {
    return GridlessMode::IsResolvedTermination(termination);
  }
};

// One policy shared by gridless and mesh-backed traces.  Numerical failures receive
// at most one smaller-step retry.  Valid TIME/STEP limit results are never sent through
// that retry path; optional unresolved extensions are separately counted and restart
// from the original phase-space seed.
struct RetryPolicy {
  int maximumNumericalRetries{1};
  double numericalDtScale{0.5};
  double numericalTimeScale{2.0};
  double numericalStepScale{2.0};
  int unresolvedExtensionPasses{0};
  double unresolvedExtensionFactor{2.0};
};

inline bool ShouldRetryNumerical(const Result& result,int retriesUsed,
                                 const RetryPolicy& policy) {
  return retriesUsed<policy.maximumNumericalRetries &&
         GridlessMode::IsRetryableNumericalTermination(result.termination);
}

inline bool IsUnresolvedExtensionCandidate(
    GridlessMode::TrajectoryTermination termination) {
  // Distance is a separately declared physical path budget and is never increased by
  // the time-convergence study.
  return termination==GridlessMode::TrajectoryTermination::TimeLimit ||
         termination==GridlessMode::TrajectoryTermination::StepLimit;
}

inline bool ShouldExtendUnresolved(const Result& result,int extensionsUsed,
                                   const RetryPolicy& policy) {
  return extensionsUsed<policy.unresolvedExtensionPasses &&
         IsUnresolvedExtensionCandidate(result.termination);
}

inline double ExtensionTimeBudget(double primaryTime_s,int extensionNumber,
                                  const RetryPolicy& policy) {
  if (!(primaryTime_s>0.0) || extensionNumber<1 ||
      !(policy.unresolvedExtensionFactor>1.0))
    return std::numeric_limits<double>::quiet_NaN();
  return primaryTime_s*std::pow(policy.unresolvedExtensionFactor,
                                static_cast<double>(extensionNumber));
}

inline int ScaledStepBudget(int primarySteps,double timeFactor,
                            const Result& previous,double targetTime_s) {
  long double desired=std::max(1,primarySteps);
  desired=std::max(desired,
      std::ceil(static_cast<long double>(primarySteps)*timeFactor*1.25L));
  if (previous.steps>0 && previous.traceTime_s>0.0) {
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

// Exact uniform-electric-field momentum update used by the U-F11 convention test and
// by the future full electromagnetic characteristic.  The momentum equation is
// dp/dt=qE; applying +dt and then -dt must therefore recover the initial state to
// floating-point roundoff.  Position integration remains mover-specific.
inline void AdvanceUniformElectricMomentum(double p_SI[3],double charge_C,
                                            const double E_V_m[3],double dt_s) {
  for (int d=0;d<3;++d) p_SI[d]+=charge_C*E_V_m[d]*dt_s;
}

inline double MomentumMagnitude(const double p_SI[3]) {
  return std::sqrt(p_SI[0]*p_SI[0]+p_SI[1]*p_SI[1]+p_SI[2]*p_SI[2]);
}

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

// Compatibility aliases keep the existing public GridlessMode and Mode3D signatures
// source-compatible while moving ownership of the records to the common contract.
namespace Earth {
namespace GridlessMode {
using TrajectoryRequest=Earth::Trajectory::Request;
using TrajectoryExitState=Earth::Trajectory::ExitState;
using TrajectoryResult=Earth::Trajectory::Result;
} // namespace GridlessMode
} // namespace Earth

#endif

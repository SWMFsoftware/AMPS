#include "../../gridless/GridlessParticleMovers.h"
#include "../../util/FluxNumerics.h"
#include "../../util/TrajectoryContract.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

namespace {

int failures=0;

void Check(bool condition,const std::string& message) {
  if (!condition) {
    ++failures;
    std::cerr << "FAIL: " << message << "\n";
  }
}

bool Near(double actual,double expected,double relative=1.0e-12,
          double absolute=1.0e-12) {
  return std::fabs(actual-expected)<=absolute+
      relative*std::max(std::fabs(actual),std::fabs(expected));
}

// A uniform magnetic field has a closed-form relativistic helix.  This is an
// independent reference solution, not a comparison with a second implementation of
// the same numerical algorithm.  It catches both dispatcher mistakes and changes in
// the expected convergence properties of the production movers.
class UniformMagneticField : public IGridlessFieldEvaluator {
public:
  explicit UniformMagneticField(double bz) : bz_(bz) {}
  void GetB_T(const V3&,V3& b) const override { b=V3{0.0,0.0,bz_}; }
private:
  double bz_;
};

struct HelixError {
  double phaseSpace;
  double relativeMomentum;
};

HelixError IntegrateUniformHelix(MoverType mover,double dt) {
  const double charge=1.0;
  const double mass=1.0;
  const double bz=1.0;
  const double totalTime=1.0;
  const double p0=1.0;
  UniformMagneticField field(bz);
  V3 x{0.0,0.0,0.0};
  V3 p{p0,0.0,0.0};
  const double initialMomentum=norm(p);
  const int steps=static_cast<int>(std::lround(totalTime/dt));
  for (int i=0;i<steps;++i)
    StepParticle(mover,x,p,charge,mass,dt,field);

  const double gamma=std::sqrt(1.0+
      p0*p0/(mass*mass*SpeedOfLight*SpeedOfLight));
  const double omega=charge*bz/(gamma*mass);
  const double velocity=p0/(gamma*mass);
  const double angle=omega*totalTime;
  // dp/dt=q v x B gives clockwise momentum rotation for q>0 and Bz>0.
  const V3 pExact{p0*std::cos(angle),-p0*std::sin(angle),0.0};
  const V3 xExact{velocity*std::sin(angle)/omega,
                  velocity*(std::cos(angle)-1.0)/omega,0.0};
  return HelixError{
      std::sqrt(dot(sub(x,xExact),sub(x,xExact))+
                dot(sub(p,pExact),sub(p,pExact))),
      std::fabs(norm(p)-initialMomentum)/initialMomentum};
}

Earth::Trajectory::Request ValidRequest() {
  Earth::Trajectory::Request request;
  request.x0_m[0]=2.0;
  request.direction0_unit[0]=1.0;
  request.restMass_kg=1.0;
  request.charge_C=1.0;
  request.rigidity_GV=1.0;
  request.maxTraceTime_s=10.0;
  request.maxTraceDistance_m=100.0;
  request.maxSteps=1000;
  request.mover=Earth::Trajectory::Mover::BORIS;
  request.backwardTimeMode=
      Earth::Trajectory::BackwardTimeMode::StaticMagneticAntiparticle;
  request.snapshotFingerprint=42;
  request.requireSnapshotIdentity=true;
  return request;
}

} // namespace

int main() {
  using namespace Earth::Trajectory;
  using Earth::GridlessMode::TrajectoryTermination;

  // U-F10: strict analytic uniform-B reference.  Both released full-orbit movers
  // must converge monotonically at their expected order.  Boris must additionally
  // preserve |p| to roundoff in a static magnetic field.
  const HelixError b1=IntegrateUniformHelix(Mover::BORIS,0.1);
  const HelixError b2=IntegrateUniformHelix(Mover::BORIS,0.05);
  const HelixError b3=IntegrateUniformHelix(Mover::BORIS,0.025);
  Check(b2.phaseSpace<b1.phaseSpace && b3.phaseSpace<b2.phaseSpace,
        "BORIS uniform-helix error decreases under timestep refinement");
  Check(b1.phaseSpace/b2.phaseSpace>1.8 && b2.phaseSpace/b3.phaseSpace>1.8,
        "BORIS retains at least second-order-like uniform-helix convergence");
  Check(b3.phaseSpace<5.0e-4,
        "BORIS finest-step phase-space error stays below the absolute reference gate");
  Check(b1.relativeMomentum<5.0e-14 && b2.relativeMomentum<5.0e-14 &&
        b3.relativeMomentum<5.0e-14,
        "BORIS preserves momentum magnitude in a static magnetic field");

  const HelixError r1=IntegrateUniformHelix(Mover::RK4,0.2);
  const HelixError r2=IntegrateUniformHelix(Mover::RK4,0.1);
  const HelixError r3=IntegrateUniformHelix(Mover::RK4,0.05);
  Check(r2.phaseSpace<r1.phaseSpace && r3.phaseSpace<r2.phaseSpace,
        "RK4 uniform-helix error decreases under timestep refinement");
  Check(r1.phaseSpace/r2.phaseSpace>8.0 && r2.phaseSpace/r3.phaseSpace>8.0,
        "RK4 retains fourth-order-like uniform-helix convergence");
  Check(r3.phaseSpace<1.0e-6,
        "RK4 finest-step phase-space error stays below the absolute reference gate");

  // U-F11: request validation and explicit backward-time semantics.  Negative cases
  // are intentional hard gates: an adapter must not normalize bad input, infer a
  // missing snapshot, or accept an unreleased electromagnetic combination.
  Request request=ValidRequest();
  Check(ValidateRequest(request)==RequestStatus::Valid,
        "static magnetic full-orbit request is valid");
  for (int code=static_cast<int>(Mover::BORIS);
       code<=static_cast<int>(Mover::HYBRID);++code) {
    const Mover mover=static_cast<Mover>(code);
    Check(IsKnownMover(mover) && std::string(MoverName(mover))!="UNKNOWN" &&
          (IsFullOrbitMover(mover)!=IsReducedOrbitMover(mover)),
          "every released mover has one unambiguous orbit class and name");
  }

  request.direction0_unit[0]=0.9;
  Check(ValidateRequest(request)==RequestStatus::InvalidDirection,
        "non-unit launch direction is rejected rather than normalized");
  request=ValidRequest();
  request.x0_m[2]=std::numeric_limits<double>::quiet_NaN();
  Check(ValidateRequest(request)==RequestStatus::NonFiniteState,
        "non-finite launch state is rejected");
  request=ValidRequest();
  request.charge_C=0.0;
  Check(ValidateRequest(request)==RequestStatus::InvalidSpecies,
        "zero charge is rejected");
  request=ValidRequest();
  request.rigidity_GV=0.0;
  Check(ValidateRequest(request)==RequestStatus::InvalidRigidity,
        "non-positive rigidity is rejected");
  request=ValidRequest();
  request.maxSteps=0;
  Check(ValidateRequest(request)==RequestStatus::InvalidBudget,
        "zero step budget is rejected");
  request=ValidRequest();
  request.mover=static_cast<Mover>(999);
  Check(ValidateRequest(request)==RequestStatus::InvalidMover,
        "unknown mover enum is rejected");
  request=ValidRequest();
  request.backwardTimeMode=static_cast<BackwardTimeMode>(999);
  Check(ValidateRequest(request)==RequestStatus::InvalidBackwardTimeMode,
        "unknown backward-time enum is rejected");
  request=ValidRequest();
  request.snapshotFingerprint=0;
  Check(ValidateRequest(request)==RequestStatus::MissingSnapshotIdentity,
        "required snapshot identity cannot be omitted");

  request=ValidRequest();
  request.electricFieldEnabled=true;
  Check(ValidateRequest(request)==RequestStatus::DynamicFieldNeedsPhysicalBackwardTime,
        "electric field rejects a static-magnetic backtrace convention");
  request.backwardTimeMode=BackwardTimeMode::PhysicalBackwardTime;
  Check(ValidateRequest(request)==RequestStatus::ElectromagneticMoverNotImplemented,
        "unreleased electromagnetic mover fails fast");
  request.electromagneticMoverImplemented=true;
  Check(ValidateRequest(request)==RequestStatus::Valid,
        "declared full-orbit electromagnetic capability satisfies the abstract gate");
  Check(!IsReleasedFrozenMagneticRequest(request),
        "released production capability still rejects an abstract EM request");
  request.mover=Mover::GC4;
  request.reducedOrbitValidityDeclared=true;
  Check(ValidateRequest(request)==RequestStatus::ReducedOrbitNeedsStaticMagneticField,
        "guiding-centre mover is rejected for E/time-dependent fields");
  request=ValidRequest();
  request.mover=Mover::HYBRID;
  Check(ValidateRequest(request)==RequestStatus::ReducedOrbitValidityNotDeclared,
        "hybrid reduced-orbit branch requires an explicit validity declaration");
  request=ValidRequest();
  request.backwardTimeMode=BackwardTimeMode::PhysicalBackwardTime;
  Check(ValidateRequest(request)==RequestStatus::Valid &&
        !IsReleasedFrozenMagneticRequest(request),
        "physical-backward enum is representable but unreleased in production");

  double momentum[3]={1.0,-2.0,0.5};
  const double original[3]={momentum[0],momentum[1],momentum[2]};
  const double electric[3]={0.25,-0.5,0.75};
  AdvanceUniformElectricMomentum(momentum,2.0,electric,0.125);
  AdvanceUniformElectricMomentum(momentum,2.0,electric,-0.125);
  for (int d=0;d<3;++d)
    Check(Near(momentum[d],original[d],0.0,2.0e-15),
          "uniform-E forward/backward momentum update is reversible");
  const double localP[3]={2.0,0.0,0.0};
  const double boundaryP[3]={1.0,0.0,0.0};
  Check(Near(PhaseSpaceIntensityFactor(localP,boundaryP),4.0),
        "Liouville intensity factor uses the squared momentum ratio");
  const double zeroP[3]={0.0,0.0,0.0};
  Check(std::isnan(PhaseSpaceIntensityFactor(localP,zeroP)),
        "Liouville factor rejects zero boundary momentum");

  // U-F12: complete boundary state is reconstructed from one exact event fraction.
  // These expected values are analytic and exercise position, momentum, direction,
  // event time, rigidity, pitch angle, validity, and failure behavior.
  ExitState exitState;
  const double xEvent[3]={10.0,20.0,30.0};
  const double pBefore[3]={3.0,4.0,0.0};
  const double pAfter[3]={5.0,0.0,0.0};
  Check(PopulateExitKinematics(exitState,xEvent,pBefore,pAfter,0.25,
                               12.0,4.0,1.0),
        "valid boundary interpolation produces kinematic exit state");
  const double pExpected=std::sqrt(21.25);
  Check(Near(exitState.x_exit_m[1],20.0) &&
        Near(exitState.p_exit_SI[0],3.5) &&
        Near(exitState.p_exit_SI[1],3.0) &&
        Near(exitState.v_exit_unit[0],3.5/pExpected) &&
        Near(exitState.v_exit_unit[1],3.0/pExpected) &&
        Near(exitState.traceTimeAtExit_s,9.0) &&
        Near(exitState.rigidityAtExit_GV,pExpected*0.299792458),
        "exit state matches the analytic event-fraction reference");
  const double parallelB[3]={7.0,6.0,0.0};
  Check(CompleteExitPitchAngle(exitState,parallelB) && exitState.valid &&
        Near(exitState.cosAlpha,1.0,0.0,2.0e-15),
        "pitch angle is evaluated from the same interpolated momentum state");
  const double zeroB[3]={0.0,0.0,0.0};
  Check(!CompleteExitPitchAngle(exitState,zeroB) && !exitState.valid,
        "zero boundary field invalidates a requested complete exit state");
  Check(!PopulateExitKinematics(exitState,xEvent,pBefore,pAfter,1.01,
                                12.0,4.0,1.0),
        "out-of-segment event fraction is rejected rather than clamped");
  Check(!PopulateExitKinematics(exitState,xEvent,pBefore,pAfter,0.5,
                                3.0,4.0,1.0),
        "accepted step cannot predate the beginning of the trace");
  ExitState incompleteExit;
  const double unitB[3]={1.0,0.0,0.0};
  Check(!CompleteExitPitchAngle(incompleteExit,unitB),
        "pitch completion rejects a missing unit exit direction");

  // A fixed digest protects cross-run provenance from implementation-dependent
  // std::hash behavior.  The expected constant is an independent FNV-1a reference.
  Check(SnapshotFingerprint("field-v1-reference")==0xfe2b5d633cd201c7ULL,
        "snapshot fingerprint matches the fixed FNV-1a reference");
  Check(SnapshotFingerprint("")==0,
        "empty snapshot identity remains the reserved zero fingerprint");
  Check(SnapshotIdentityMatches(42,true,42) &&
        !SnapshotIdentityMatches(42,true,43) &&
        !SnapshotIdentityMatches(42,false,0) &&
        SnapshotIdentityMatches(0,false,43),
        "snapshot identity assertion accepts only the declared active generation");

  // U-F13: numerical retry and unresolved extension are bounded, distinct, and leave
  // DISTANCE_LIMIT untouched.  Count closure covers every terminal category used by
  // the cutoff and flux accumulators.
  RetryPolicy policy;
  policy.maximumNumericalRetries=1;
  policy.unresolvedExtensionPasses=2;
  Result result;
  result.termination=TrajectoryTermination::OuterBoundaryAllowed;
  Check(result.allowed() && result.resolved(),
        "outer-boundary termination is resolved allowed");
  result.termination=TrajectoryTermination::TimeLimit;
  Check(!result.allowed() && !result.resolved(),
        "time-limit termination remains unresolved rather than forbidden");
  result.termination=TrajectoryTermination::NumericalFailure;
  Check(ShouldRetryNumerical(result,0,policy) &&
        !ShouldRetryNumerical(result,1,policy),
        "numerical retry is allowed once and then stops");
  result.termination=TrajectoryTermination::TimeLimit;
  Check(!ShouldRetryNumerical(result,0,policy) &&
        ShouldExtendUnresolved(result,0,policy) &&
        ShouldExtendUnresolved(result,1,policy) &&
        !ShouldExtendUnresolved(result,2,policy),
        "trace-limit extension is bounded and separate from numerical retry");
  result.termination=TrajectoryTermination::DistanceLimit;
  Check(!ShouldExtendUnresolved(result,0,policy),
        "distance budget is never relaxed by a time-convergence extension");
  Check(Near(ExtensionTimeBudget(10.0,2,policy),40.0),
        "extension time follows the declared geometric policy");
  result.steps=50;
  result.traceTime_s=10.0;
  Check(ScaledStepBudget(100,4.0,result,40.0)==500,
        "step budget includes the fixed 25-percent convergence margin");
  RetryPolicy invalidPolicy=policy;
  invalidPolicy.unresolvedExtensionFactor=1.0;
  Check(!IsValidRetryPolicy(invalidPolicy) &&
        !ShouldExtendUnresolved(result,0,invalidPolicy) &&
        std::isnan(ExtensionTimeBudget(10.0,1,invalidPolicy)),
        "invalid retry policy fails closed");

  Earth::FluxNumerics::AccessAccumulator counts;
  for (int code=0;code<static_cast<int>(TrajectoryTermination::Count);++code)
    counts.Record(static_cast<TrajectoryTermination>(code),0.0);
  long long closed=0;
  for (std::size_t i=0;i<counts.terminationCounts.size();++i)
    closed+=counts.terminationCounts[i];
  Check(closed==counts.sampled &&
        counts.sampled==static_cast<long long>(TrajectoryTermination::Count),
        "termination counts close exactly over every terminal category");

  if (failures!=0) {
    std::cerr << "UTrajectoryCore: " << failures << " failure(s)\n";
    return EXIT_FAILURE;
  }
  std::cout << "UTrajectoryCore: PASS (U-F10..U-F13)\n";
  return EXIT_SUCCESS;
}

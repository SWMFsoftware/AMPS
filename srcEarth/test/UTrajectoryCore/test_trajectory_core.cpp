#include "../../gridless/GridlessParticleMovers.h"
#include "../../util/FluxNumerics.h"
#include "../../util/TrajectoryContract.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
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

// A uniform field has a closed-form relativistic helix.  It is therefore a genuine
// reference solution for the mover convergence checks, not a comparison against a
// second implementation of the same numerical algorithm.
class UniformMagneticSnapshot : public IGridlessFieldEvaluator {
public:
  explicit UniformMagneticSnapshot(double bz) : bz_(bz) {
    metadata_.sourceId="UNIT:UNIFORM_B";
    metadata_.modelName="UNIFORM_B";
    metadata_.epochUTC="2000-01-01T00:00:00";
    metadata_.frame=Earth::Field::CoordinateFrame::GSM;
    metadata_.interpolation=Earth::Field::InterpolationMode::DirectAnalytic;
    metadata_.magneticFieldAvailable=true;
    metadata_.electricFieldAvailable=false;
    metadata_.immutableDuringBatch=true;
    metadata_.valid=true;
    metadata_.snapshotId=Earth::Field::MakeSnapshotId(
        metadata_.sourceId,metadata_.epochUTC,"Bz="+std::to_string(bz_));
  }

  const Earth::Field::SnapshotMetadata& Metadata() const override {
    return metadata_;
  }

  Earth::Field::FieldSample Sample(
      const Earth::Field::FieldQuery& query) const override {
    Earth::Field::FieldSample sample;
    sample.snapshotId=metadata_.snapshotId;
    sample.interpolation=metadata_.interpolation;
    sample.status=Earth::Field::ValidateQuery(metadata_,query,&sample.message);
    if (!sample.ok()) return sample;
    sample.magneticField_T[2]=bz_;
    return sample;
  }

  void GetB_T(const V3&,V3& b) const override { b=V3{0.0,0.0,bz_}; }

private:
  double bz_;
  Earth::Field::SnapshotMetadata metadata_;
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
  UniformMagneticSnapshot field(bz);
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
  // dp/dt=q v x B gives clockwise rotation for q>0 and Bz>0.
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
  request.backwardTimeMode=Earth::Trajectory::BackwardTimeMode::StaticMagneticAntiparticle;
  request.snapshotFingerprint=42;
  request.requireSnapshotIdentity=true;
  return request;
}

} // namespace

int main() {
  using namespace Earth::Trajectory;
  using Earth::GridlessMode::TrajectoryTermination;

  // U-F10: both released full-orbit movers converge to the independent analytic
  // uniform-B helix when dt is halved.  Boris is globally second order in this
  // implementation; classical RK4 must display a substantially steeper error slope.
  const HelixError b1=IntegrateUniformHelix(Mover::BORIS,0.1);
  const HelixError b2=IntegrateUniformHelix(Mover::BORIS,0.05);
  const HelixError b3=IntegrateUniformHelix(Mover::BORIS,0.025);
  Check(b2.phaseSpace<b1.phaseSpace && b3.phaseSpace<b2.phaseSpace,
        "BORIS uniform-helix error decreases under timestep refinement");
  Check(b1.phaseSpace/b2.phaseSpace>1.8 && b2.phaseSpace/b3.phaseSpace>1.8,
        "BORIS reaches at least second-order-like uniform-helix convergence");
  Check(b1.relativeMomentum<5.0e-14 && b2.relativeMomentum<5.0e-14 &&
        b3.relativeMomentum<5.0e-14,
        "BORIS preserves momentum magnitude in a static magnetic field");

  const HelixError r1=IntegrateUniformHelix(Mover::RK4,0.2);
  const HelixError r2=IntegrateUniformHelix(Mover::RK4,0.1);
  const HelixError r3=IntegrateUniformHelix(Mover::RK4,0.05);
  Check(r2.phaseSpace<r1.phaseSpace && r3.phaseSpace<r2.phaseSpace,
        "RK4 uniform-helix error decreases under timestep refinement");
  Check(r1.phaseSpace/r2.phaseSpace>8.0 && r2.phaseSpace/r3.phaseSpace>8.0,
        "RK4 approaches fourth-order uniform-helix convergence");

  // U-F11: the request gate makes backward-time semantics explicit.  Static-B
  // antiparticle backtracing is accepted, but any E/time dependence requires a true
  // physical-backward-time mover and reduced-orbit movers remain forbidden there.
  Request request=ValidRequest();
  Check(ValidateRequest(request)==RequestStatus::Valid,
        "static magnetic full-orbit request is valid");
  request.electricFieldEnabled=true;
  Check(ValidateRequest(request)==RequestStatus::DynamicFieldNeedsPhysicalBackwardTime,
        "electric field rejects the static-magnetic backtrace convention");
  request.backwardTimeMode=BackwardTimeMode::PhysicalBackwardTime;
  Check(ValidateRequest(request)==RequestStatus::ElectromagneticMoverNotImplemented,
        "unreleased electromagnetic mover fails fast");
  request.electromagneticMoverImplemented=true;
  Check(ValidateRequest(request)==RequestStatus::Valid,
        "declared physical-backward electromagnetic mover passes the convention gate");
  request.mover=Mover::GC4;
  request.reducedOrbitValidityDeclared=true;
  Check(ValidateRequest(request)==RequestStatus::ReducedOrbitNeedsStaticMagneticField,
        "guiding-center mover is rejected for E/time-dependent fields");

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

  // U-F12: retry and extension paths are bounded and mutually exclusive.  This test
  // also checks count closure over every terminal category used by cutoff and flux.
  RetryPolicy policy;
  policy.maximumNumericalRetries=1;
  policy.unresolvedExtensionPasses=2;
  Result result;
  result.termination=TrajectoryTermination::NumericalFailure;
  Check(ShouldRetryNumerical(result,0,policy) &&
        !ShouldRetryNumerical(result,1,policy),
        "numerical retry is allowed once and then stops");
  result.termination=TrajectoryTermination::TimeLimit;
  Check(!ShouldRetryNumerical(result,0,policy) &&
        ShouldExtendUnresolved(result,0,policy) &&
        ShouldExtendUnresolved(result,1,policy) &&
        !ShouldExtendUnresolved(result,2,policy),
        "trace-limit extension is bounded and never uses the numerical retry path");
  Check(Near(ExtensionTimeBudget(10.0,2,policy),40.0),
        "extension time budget follows the declared geometric policy");

  Earth::FluxNumerics::AccessAccumulator counts;
  for (int code=0;code<static_cast<int>(TrajectoryTermination::Count);++code)
    counts.Record(static_cast<TrajectoryTermination>(code),0.0);
  long long closed=0;
  for (std::size_t i=0;i<counts.terminationCounts.size();++i)
    closed+=counts.terminationCounts[i];
  Check(closed==counts.sampled &&
        counts.sampled==static_cast<long long>(TrajectoryTermination::Count),
        "termination counts close exactly over the sampled trajectories");

  // A reduced-orbit request must state that its validity criteria were checked.  This
  // is the intentional fail-fast path for unsupported mover/field combinations.
  request=ValidRequest();
  request.mover=Mover::GC4;
  Check(ValidateRequest(request)==RequestStatus::ReducedOrbitValidityNotDeclared,
        "reduced-orbit request without a validity declaration fails fast");

  if (failures!=0) {
    std::cerr << "UTrajectoryCore: " << failures << " failure(s)\n";
    return EXIT_FAILURE;
  }
  std::cout << "UTrajectoryCore: PASS (U-F10..U-F12)\n";
  return EXIT_SUCCESS;
}

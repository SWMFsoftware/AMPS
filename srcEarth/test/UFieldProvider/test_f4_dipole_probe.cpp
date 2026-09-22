#include "../../gridless/DipoleInterface.h"
#include "../../gridless/GridlessParticleMovers.h"
#include "../../util/TrajectoryBoundary.h"
#include "../../util/TrajectoryTrapDetector.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

//====================================================================================
// F4 analytic-dipole viability regression for the Step-3 field snapshot
//====================================================================================
//
// F4 is an algebraic closure test: every saved energy node must have a finite nominal
// transmission before the Python runner can reconstruct density and flux.  The nominal
// transmission is deliberately undefined when *all* sampled directions end at a
// numerical time/step/distance limit.
//
// The former setup used RK4 with the trap detector's production-tight relative
// momentum tolerance (1e-4).  In the exact static dipole, |p| is an invariant, but
// explicit RK4 accumulated roughly 0.3--0.8 percent drift over the long low-energy
// F4 traces.  The detector correctly refused to turn those numerically degraded paths
// into a physical TRAPPED verdict.  At several latitude/energy pairs every one of the
// 16 paths was therefore unresolved and F4 emitted NaN.
//
// The closure test now uses the magnetic Boris mover, which preserves |p| to roundoff,
// and an explicit 1-Re bounce-envelope tolerance appropriate to this deliberately
// sparse angular sample.  This test mirrors the resulting F4 setup against the analytic
// dipole for all 5 x 32 energy/location bins and asserts that each bin contains at least
// one physical OUTER, INNER, or recurrent-TRAPPED classification.  Numerical limits
// are still counted as unresolved and are never relabeled.
//
// This source stays independent of the AMPS executable, MPI, Geopack, and SWMF.  It
// links the production mover, boundary-event logic, trap detector, and stateless
// Dipole::Params sampling path introduced by the Step-3 fix.

namespace {

constexpr double kPi=3.141592653589793238462643383279502884;
constexpr double kRe=Earth::GridlessMode::Dipole::Re_m;
constexpr double kQ=Earth::GridlessMode::Dipole::q_e;
constexpr double kMass=Earth::GridlessMode::Dipole::m0;

class FrozenDipoleField : public IGridlessFieldEvaluator {
 public:
  FrozenDipoleField()
      : params_(Earth::GridlessMode::Dipole::MakeParams(1.0,0.0)) {
    metadata_.snapshotId="f4-dipole-probe";
    metadata_.sourceId="TEST:STATELESS_DIPOLE";
    metadata_.modelName="DIPOLE";
    metadata_.epochUTC="2026-07-10T00:00:00";
    metadata_.frame=Earth::Field::CoordinateFrame::GSM;
    metadata_.interpolation=Earth::Field::InterpolationMode::DirectAnalytic;
    metadata_.magneticFieldAvailable=true;
    metadata_.immutableDuringBatch=true;
    metadata_.valid=true;
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
    Earth::GridlessMode::Dipole::GetB_Tesla(
        query.position_m,sample.magneticField_T,params_);
    return sample;
  }

  void GetB_T(const V3& x,V3& b) const override {
    const double position[3]={x.x,x.y,x.z};
    double field[3]={0.0,0.0,0.0};
    Earth::GridlessMode::Dipole::GetB_Tesla(position,field,params_);
    b=V3{field[0],field[1],field[2]};
  }

 private:
  const Earth::GridlessMode::Dipole::Params params_;
  Earth::Field::SnapshotMetadata metadata_;
};

double MomentumFromEnergy(double energyMeV) {
  const double kineticJ=energyMeV*1.0e6*kQ;
  const double gamma=1.0+kineticJ/(kMass*SpeedOfLight*SpeedOfLight);
  return kMass*SpeedOfLight*std::sqrt(gamma*gamma-1.0);
}

enum class Termination { Outer, Inner, Trapped, Time, Steps, Distance, Bad };

Termination Trace(const FrozenDipoleField& field,const V3& x0,
                  const V3& direction,double energyMeV) {
  V3 x=x0;
  V3 p=mul(MomentumFromEnergy(energyMeV),mul(-1.0,direction));
  const Earth::TrajectoryBoundary::Box box{
      {-20.0*kRe,-20.0*kRe,-20.0*kRe},
      { 20.0*kRe, 20.0*kRe, 20.0*kRe},
      1.01*kRe};

  Earth::TrajectoryTrap::Config cfg;
  cfg.enabled=true;
  cfg.minMirrorPoints=8;
  cfg.minBounceCycles=4;
  cfg.outerMargin_m=1.0*kRe;
  cfg.radialEnvelopeTolerance_m=1.0*kRe;
  cfg.energyRelativeTolerance=1.0e-4;
  cfg.parallelDeadband=1.0e-6;
  Earth::TrajectoryTrap::Detector detector(cfg,box);

  double time=0.0;
  double distance=0.0;
  int steps=0;

  V3 b;
  field.GetB_T(x,b);
  double xa[3]={x.x,x.y,x.z};
  double pa[3]={p.x,p.y,p.z};
  double ba[3]={b.x,b.y,b.z};
  detector.Update(xa,pa,ba);

  while (steps<100000 && time<120.0 && distance<300.0*kRe) {
    xa[0]=x.x; xa[1]=x.y; xa[2]=x.z;
    if (Earth::TrajectoryBoundary::InsideInnerSphere(xa,box))
      return Termination::Inner;
    if (!Earth::TrajectoryBoundary::InsideBox(xa,box))
      return Termination::Outer;

    const double gamma=std::sqrt(
        1.0+dot(p,p)/std::pow(kMass*SpeedOfLight,2));
    field.GetB_T(x,b);
    const double omega=kQ*norm(b)/(gamma*kMass);
    double dt=std::min(0.1,0.15/omega);
    dt=std::min(dt,120.0-time);
    if (!(dt>0.0) || !std::isfinite(dt)) return Termination::Bad;

    const V3 before=x;
    if (!StepParticleChecked(
            MoverType::BORIS,x,p,kQ,kMass,dt,field,box.innerRadius))
      return Termination::Inner;

    const double beforeArray[3]={before.x,before.y,before.z};
    const double afterArray[3]={x.x,x.y,x.z};
    const Earth::TrajectoryBoundary::Event event=
        Earth::TrajectoryBoundary::FindFirstEvent(
            beforeArray,afterArray,box,1.0);
    distance+=norm(sub(x,before));
    time+=dt;
    ++steps;
    if (event.type==Earth::TrajectoryBoundary::EventType::InnerSphere)
      return Termination::Inner;
    if (event.type==Earth::TrajectoryBoundary::EventType::OuterBox)
      return Termination::Outer;

    field.GetB_T(x,b);
    xa[0]=x.x; xa[1]=x.y; xa[2]=x.z;
    pa[0]=p.x; pa[1]=p.y; pa[2]=p.z;
    ba[0]=b.x; ba[1]=b.y; ba[2]=b.z;
    if (detector.Update(xa,pa,ba)) return Termination::Trapped;
  }

  if (steps>=100000) return Termination::Steps;
  if (distance>=300.0*kRe) return Termination::Distance;
  if (time>=120.0) return Termination::Time;
  return Termination::Bad;
}

std::vector<V3> F4Directions() {
  std::vector<V3> full;
  full.reserve(24*48);
  for (int i=0;i<24;++i) {
    const double mu=-1.0+2.0*(i+0.5)/24.0;
    const double sinTheta=std::sqrt(1.0-mu*mu);
    for (int j=0;j<48;++j) {
      const double phi=2.0*kPi*(j+0.5)/48.0;
      full.push_back(V3{sinTheta*std::cos(phi),sinTheta*std::sin(phi),mu});
    }
  }

  // Same deterministic 16-of-1152 selection used when F4 sets
  // DS_MAX_PARTICLES=512 over its 32-point energy grid.
  std::vector<V3> selected;
  selected.reserve(16);
  for (int k=0;k<16;++k) {
    const double fraction=static_cast<double>(k)/15.0;
    const int index=static_cast<int>(
        std::floor(fraction*(full.size()-1)+0.5));
    selected.push_back(full[static_cast<std::size_t>(index)]);
  }
  return selected;
}

bool Resolved(Termination termination) {
  return termination==Termination::Outer ||
         termination==Termination::Inner ||
         termination==Termination::Trapped;
}

} // namespace

int main() {
  const FrozenDipoleField field;
  const std::vector<V3> directions=F4Directions();
  const double latitudes[]={-60.0,-30.0,0.0,30.0,60.0};
  const double radius=(6371.2+9000.0)*1000.0;
  int checkedBins=0;
  int minimumResolved=static_cast<int>(directions.size());
  int failures=0;

  for (double latitude:latitudes) {
    const double latitudeRad=latitude*kPi/180.0;
    const V3 x0{radius*std::cos(latitudeRad),0.0,
                radius*std::sin(latitudeRad)};
    for (int energyIndex=0;energyIndex<32;++energyIndex) {
      const double energyMeV=std::pow(
          1000.0,static_cast<double>(energyIndex)/31.0);
      int resolved=0;
      for (const V3& direction:directions) {
        if (Resolved(Trace(field,x0,direction,energyMeV))) ++resolved;
      }
      ++checkedBins;
      minimumResolved=std::min(minimumResolved,resolved);
      if (resolved==0) {
        ++failures;
        std::cerr << "FAIL: F4 analytic-dipole bin has no resolved trajectories: "
                  << "latitude=" << latitude
                  << " deg, energy=" << energyMeV << " MeV\n";
      }
    }
  }

  if (failures!=0) {
    std::cerr << "UFieldProvider F4 probe: " << failures
              << " all-unresolved energy bin(s)\n";
    return EXIT_FAILURE;
  }
  std::cout << "UFieldProvider F4 probe: PASS (" << checkedBins
            << " bins, minimum resolved directions=" << minimumResolved << ")\n";
  return EXIT_SUCCESS;
}

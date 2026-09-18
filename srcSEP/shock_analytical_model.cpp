#include "sep.h"
#include "adapters/swcme1d_adapter.h"
#include "util/sep_shock_source_core.h"
#include "turbulence_production_adapter.h"
#include "util/sep_run_configuration.h"
//analytic model of a shock wave (Tenishev-2005-AIAA-4928

#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
#include "../srcInterface/amps2swmf.h"
#endif

double SEP::ParticleSource::ShockWave::Tenishev2005::rShock =1.0E-5*_SUN__RADIUS_;
bool SEP::ParticleSource::ShockWave::Tenishev2005::InitFlag=false;
double SEP::ParticleSource::ShockWave::Tenishev2005::MinFieldLineHeliocentricDistance=-1.0;

namespace {

SEP::Shock::TrajectoryState gAnalyticalShockState;

SEP::Shock::TrajectoryConfiguration AnalyticalShockConfiguration(
    double launchRadiusM) {
  SEP::Shock::TrajectoryConfiguration configuration;
  configuration.launchEpochS=0.0;
  configuration.launchRadiusM=launchRadiusM;
  const double radiiAu[]={0.1,0.15,0.3,0.5,1.3};
  const double speedsKmS[]={1800.0,1500.0,1500.0,1100.0,900.0};
  for (std::size_t i=0;i<sizeof(radiiAu)/sizeof(radiiAu[0]);++i) {
    SEP::Shock::SpeedKnot knot;
    knot.radiusM=radiiAu[i]*_AU_;
    knot.speedMPerS=speedsKmS[i]*1.0e3;
    configuration.knots.push_back(knot);
  }
  return configuration;
}

}  // namespace


void SEP::ParticleSource::ShockWave::Tenishev2005::Init() {
  InitFlag=true;

  //determine the  initial location of the shock that is the minimum helpocentric distance of the beginning of the simulated field lines
  //loop through all simulated field lines
  for (int i = 0; i < PIC::FieldLine::nFieldLine; i++) {
    double r, *x;

    x = PIC::FieldLine::FieldLinesAll[i].GetFirstSegment()->GetBegin()->GetX();
    r = Vector3D::Length(x);

    if ((MinFieldLineHeliocentricDistance < 0.0) || (MinFieldLineHeliocentricDistance > r))
      MinFieldLineHeliocentricDistance = r;
  }

  rShock=MinFieldLineHeliocentricDistance;

  // Establish a complete analytical state at the physical launch epoch.  The
  // state is advanced from this datum with exact piecewise integration, so a
  // restart at nonzero time and a differently partitioned timestep reproduce
  // the same radius instead of depending on a function-static last time.
  const SEP::Transport::Status trajectoryStatus=SEP::Shock::StateAtEpoch(
      AnalyticalShockConfiguration(rShock),0.0,&gAnalyticalShockState);
  if (!trajectoryStatus.ok())
    exit(__LINE__,__FILE__,trajectoryStatus.message.c_str());

}

double SEP::ParticleSource::ShockWave::Tenishev2005::GetCompressionRatio() {
  double r = rShock / _AU_;
  double res;

  if (r<0.04) {
    res=1.7;
  }
  else {
    res=2.0+(1.4-2.0)/(0.14-0.04)*(r-0.04);
    if (res<1.0) res=1.0;
  }

  return res;
}

double SEP::ParticleSource::ShockWave::Tenishev2005::GetShockSpeed() {
  if (SEP::ShockModelType==SEP::cShockModelType::SwCme1d) {
    return SEP::SW1DAdapter::ShockSpeedMPerS();
  }
  if (InitFlag==false) Init();
  const SEP::Transport::ScalarResult speed=SEP::Shock::SpeedAtRadius(
      gAnalyticalShockState.configuration,rShock);
  if (!speed.status.ok()) exit(__LINE__,__FILE__,speed.status.message.c_str());
  return speed.value;
}

void SEP::ParticleSource::ShockWave::Tenishev2005::UpdateShockLocation() {
  if (InitFlag==false) Init();

  // PIC owns the only simulation clock.  Read it through the same srcSEP clock
  // adapter used by background snapshots and output, rather than accessing the
  // mutable TimeCounter implementation detail or maintaining another elapsed
  // clock in this shock model.
  const double simulation_time = SEP::Background::SimulationTimeSeconds();

  if (SEP::ShockModelType==SEP::cShockModelType::SwCme1d) {
    rShock=SEP::SW1DAdapter::ShockRadiusM();
    return;
  }
  const SEP::Transport::Status status=
      SEP::Shock::AdvanceToEpoch(&gAnalyticalShockState,simulation_time);
  if (!status.ok()) exit(__LINE__,__FILE__,status.message.c_str());
  rShock=gAnalyticalShockState.radiusM;
}

double SEP::ParticleSource::ShockWave::Tenishev2005::GetSolarWindDensity() {
  double t=_AU_/rShock;
  return 5.0E6*t*t;
}

double SEP::ParticleSource::ShockWave::Tenishev2005::GetInjectionRate() {
  double r_sh,density;

  if (InitFlag==false) Init();

  switch (SEP::ShockModelType) {
  case SEP::cShockModelType::Analytic1D:
    density=GetSolarWindDensity();
    break;
  case SEP::cShockModelType::SwCme1d:
    r_sh=SEP::SW1DAdapter::ShockRadiusM();

    // D01: consume one transactional sample. A rejected provider query is a
    // fatal background-consistency error, not a zero-density interval. The
    // formatted message records the MPI rank, physical epoch, radius, failed
    // field, and immutable SWCME state ID needed to reproduce the failure.
    {
      const SEP::SW1DAdapter::QueryResult query =
          SEP::SW1DAdapter::QueryAtRadius(r_sh, /*applySheathClamp=*/true);
      if (!query.ok()) {
        const std::string diagnostic =
            SEP::SW1DAdapter::FormatFailure(query, PIC::ThisThread);
        exit(__LINE__,__FILE__,diagnostic.c_str());
      }
      density=query.sample.number_density_m3();
    }

    break;
  default:
    exit(__LINE__,__FILE__,"Error: the case is not known");
  }

  // Apply the single configured source efficiency used by every background
  // provider.  The former analytic-only (compression-1)/compression factor
  // made identical analytic, SWCME, and SWMF shock states inject different
  // physical particle counts.
  const double sourceWeight=SEP::ShockModelType==SEP::cShockModelType::SwCme1d
      ? SEP::SW1DAdapter::RelativeSourceWeightPerArea() : 1.0;
  return sourceWeight*
      SEP::FieldLine::FluxTubeGeometryCore::InjectedPhysicalParticleCount(
      SEP::Units::NumberDensityPerM3(density), SEP::Units::VolumeM3(1.0),
      SEP::FieldLine::InjectionParameters::InjectionEfficiency);
}

int SEP::ParticleSource::ShockWave::Tenishev2005::GetInjectionLocation(int iFieldLine,double &S,double *xInjection) {
  // Update the authoritative trajectory before deriving the surface radius.
  // The former order used a stale rShock^2 for analytical runs.
  UpdateShockLocation();
  const double radius=SEP::ShockModelType==SEP::cShockModelType::SwCme1d
      ? SEP::SW1DAdapter::ShockRadiusM() : rShock;
  std::vector<SEP::Shock::Vector3> vertices;
  PIC::FieldLine::cFieldLineSegment* segment=
      PIC::FieldLine::FieldLinesAll[iFieldLine].GetFirstSegment();
  if (segment==NULL) { S=-1.0; return -1; }
  for (;segment!=NULL;segment=segment->GetNext()) {
    if (vertices.empty()) {
      const double* x=segment->GetBegin()->GetX();
      SEP::Shock::Vector3 point; point.x=x[0];point.y=x[1];point.z=x[2];
      vertices.push_back(point);
    }
    const double* x=segment->GetEnd()->GetX();
    SEP::Shock::Vector3 point; point.x=x[0];point.y=x[1];point.z=x[2];
    vertices.push_back(point);
  }
  const SEP::Shock::IntersectionResult result=SEP::Shock::IntersectSphere(
      vertices,radius,std::max(1.0,1.0e-12*radius),
      SEP::Shock::IntersectionPolicy::FirstOutward);
  if (result.status!=SEP::Shock::IntersectionStatus::Ok) { S=-1.0; return -1; }
  const SEP::Shock::Intersection& hit=result.intersections.front();
  xInjection[0]=hit.positionM.x; xInjection[1]=hit.positionM.y;
  xInjection[2]=hit.positionM.z;
  S=static_cast<double>(hit.segment)+hit.fraction;
  return static_cast<int>(hit.segment);
}


namespace SEP {
namespace ParticleSource {
namespace ShockWave {

// Function to increment integrated wave energy due to shock passing
// r0, r1: initial and final heliocentric distances of the shock
// dt: time needed for shock to move from r0 to r1
void ShockTurbulenceEnergyInjection(double r0, double r1, double dt) {
    using namespace PIC::FieldLine;
    if (!(r1>r0) || !(r0>=0.0) || !(dt>0.0) || !std::isfinite(dt)) {
      SEP::Turbulence::PICAdapter::RecordShockSourceRejection(true);
      return;
    }
    const double shockSpeed=(r1-r0)/dt;
    const SEP::Run::Configuration& run=SEP::Run::Active().get();
    const double efficiency=run.shockTurbulenceEfficiency.value;
    const double plusFraction=run.shockTurbulencePlusFraction.value;

    for (int iFieldLine=0;iFieldLine<nFieldLine;++iFieldLine) {
      PIC::FieldLine::cFieldLine* line=&FieldLinesAll[iFieldLine];
      for (int iSegment=0;iSegment<line->GetTotalSegmentNumber();++iSegment) {
        cFieldLineSegment* segment=line->GetSegment(iSegment);
        if (segment==NULL) continue;
        const double* begin=segment->GetBegin()->GetX();
        const double* end=segment->GetEnd()->GetX();
        std::vector<SEP::Shock::Vector3> chord(2);
        chord[0].x=begin[0];chord[0].y=begin[1];chord[0].z=begin[2];
        chord[1].x=end[0];chord[1].y=end[1];chord[1].z=end[2];
        const double tolerance=std::max(1.0,1.0e-12*r1);
        SEP::Shock::IntersectionResult inner;
        if (r0>0.0) inner=SEP::Shock::IntersectSphere(
            chord,r0,tolerance,SEP::Shock::IntersectionPolicy::All);
        const SEP::Shock::IntersectionResult outer=SEP::Shock::IntersectSphere(
            chord,r1,tolerance,SEP::Shock::IntersectionPolicy::All);
        if (inner.status==SEP::Shock::IntersectionStatus::InvalidGeometry ||
            outer.status==SEP::Shock::IntersectionStatus::InvalidGeometry) {
          SEP::Turbulence::PICAdapter::RecordShockSourceRejection(true);
          continue;
        }
        // Sphere roots partition the straight PIC segment into intervals that
        // are wholly inside or outside the swept shell.  Midpoint classification
        // handles inward, outward, tangent, and two-crossing chords uniformly.
        std::vector<double> boundaries;
        boundaries.push_back(0.0);boundaries.push_back(1.0);
        for (std::size_t hit=0;hit<inner.intersections.size();++hit)
          boundaries.push_back(inner.intersections[hit].fraction);
        for (std::size_t hit=0;hit<outer.intersections.size();++hit)
          boundaries.push_back(outer.intersections[hit].fraction);
        std::sort(boundaries.begin(),boundaries.end());
        boundaries.erase(std::unique(boundaries.begin(),boundaries.end(),
            [tolerance,segment](double a,double b) {
              return std::fabs(a-b)*segment->GetLength()<=tolerance;
            }),boundaries.end());

        for (std::size_t interval=0;interval+1<boundaries.size();++interval) {
          const double start=boundaries[interval];
          const double finish=boundaries[interval+1];
          if (!(finish>start)) continue;
          double testPoint[3]={0.0,0.0,0.0};
          segment->GetCartesian(testPoint,0.5*(start+finish));
          const double testRadius=Vector3D::Length(testPoint);
          if (testRadius<r0 || testRadius>=r1) continue;

          double density0=0.0,density1=0.0;
          segment->GetPlasmaDensity(start,density0);
          segment->GetPlasmaDensity(finish,density1);
          double velocity0[3]={0.0,0.0,0.0};
          double velocity1[3]={0.0,0.0,0.0};
          segment->GetBegin()->GetPlasmaVelocity(velocity0);
          segment->GetEnd()->GetPlasmaVelocity(velocity1);
          const double radius=Vector3D::Length(testPoint);
          const double upstreamNormalSpeed=radius>0.0 ?
              0.5*((velocity0[0]+velocity1[0])*testPoint[0]+
                   (velocity0[1]+velocity1[1])*testPoint[1]+
                   (velocity0[2]+velocity1[2])*testPoint[2])/radius : 0.0;

          SEP::Shock::TurbulenceSourceInput input;
          input.sweptVolumeM3=SEP::FieldLine::FluxTubeGeometry::
              PartialSegmentVolumeM3(segment,iFieldLine,start,finish);
          input.upstreamMassDensityKgPerM3=
              0.5*(density0+density1)*PIC::CPLR::SWMF::MeanPlasmaAtomicMass;
          input.shockNormalSpeedMPerS=shockSpeed;
          input.upstreamNormalSpeedMPerS=upstreamNormalSpeed;
          input.efficiency=efficiency;
          input.plusBranchFraction=plusFraction;
          const SEP::Shock::TurbulenceSourceEnergy energy=
              SEP::Shock::ComputeTurbulenceSource(input);
          if (!energy.status.ok()) {
            SEP::Turbulence::PICAdapter::RecordShockSourceRejection(false);
            continue;
          }
          SEP::Turbulence::PICAdapter::ShockContribution contribution;
          contribution.fieldLine=iFieldLine;
          contribution.segment=iSegment;
          contribution.plusJ=energy.plusJ;
          contribution.minusJ=energy.minusJ;
          contribution.provenance="WP25 sphere-shell overlap; upstream normal-relative speed";
          SEP::Turbulence::PICAdapter::QueueShockContribution(contribution);
        }
      }
    }
}

} // namespace ShockWave
} // namespace ParticleSource
} // namespace SEP

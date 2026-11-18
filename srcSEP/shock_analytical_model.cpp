#include "sep.h"
#include "adapters/swcme1d_adapter.h"
#include "util/sep_shock_source_core.h"
#include "util/sep_analytic_shock_epoch.h"
#include "util/sep_analytic_radial_background.h"
#include "util/sep_background_runtime.h"
#include "turbulence_production_adapter.h"
#include "util/sep_run_configuration.h"
#include "physics_checkpoint_state.h"
#include "transport_common.h"
#include <utility>
#include <limits>
#include <cstdio>
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

double CompressionAtRadius(double radiusM) {
  // Preserve the existing declared radius-dependent analytic profile, including
  // its compression-one (no acceleration) branch. This is the model's existing
  // piecewise law, not a numerical repair of a bad input or a new thermal/RH
  // closure. Temporal consumers now pass the radius of their explicit view.
  const double radiusAu=radiusM/_AU_;
  if(radiusAu<0.04)return 1.7;
  const double value=2.0+(1.4-2.0)/(0.14-0.04)*(radiusAu-0.04);
  return value<1.0?1.0:value;
}
double DensityAtRadius(double radiusM) {
  const double inverseRadiusAu=_AU_/radiusM;
  return 5.0e6*inverseRadiusAu*inverseRadiusAu;
}

SEP::Transport::Status PrepareTrajectoryCandidate(double epochS,bool mayInitialize,
    SEP::PhysicsCheckpoint::AnalyticShockCache* candidate) {
  namespace T=SEP::Transport;
  if(SEP::ShockModelType!=SEP::cShockModelType::Analytic1D||
      SEP::Background::ConfiguredProvider()!=SEP::Background::Provider::Analytic)
    return T::Status::Error(T::StatusCode::UnsupportedConfiguration,
        "explicit analytic trajectory owner requires the analytic background/shock authority; SWCME uses its prepared provider state");
  if(SEP::Background::AnalyticRadial::DeclaredNativeWind().originM!=std::array<double,3>{{0,0,0}})
    return T::Status::Error(T::StatusCode::UnsupportedConfiguration,
        "legacy analytic shock sphere/density are global-zero-centered; a translated source/background count contract has not been declared");
  SEP::Background::SnapshotStore::Instance().AssertProviderMayWrite(SEP::Background::Provider::Analytic);
  if(!std::isfinite(epochS))return T::Status::Error(T::StatusCode::InvalidArgument,"analytic trajectory owner epoch must be finite [s]");
  SEP::Shock::TrajectoryConfiguration configuration;
  if(SEP::ParticleSource::ShockWave::Tenishev2005::InitFlag) {
    const auto captured=SEP::ParticleSource::ShockWave::Tenishev2005::CaptureAnalyticShockCache(candidate);
    if(!captured.ok())return captured;
    if(epochS<candidate->trajectory.epochS)
      return T::Status::Error(T::StatusCode::InvalidParticleState,"analytic trajectory owner cannot advance backwards");
    configuration=candidate->trajectory.configuration;
  } else {
    if(!mayInitialize)return T::Status::Error(T::StatusCode::InvalidParticleState,
        "analytic trajectory must be explicitly initialized before completing a host interval");
    double minimumRadius=-1;
    for(int line=0;line<PIC::FieldLine::nFieldLine;++line) {
      auto* first=PIC::FieldLine::FieldLinesAll[line].GetFirstVertex();
      if(!first)return T::Status::Error(T::StatusCode::InvalidParticleState,"analytic trajectory launch needs nonempty native field lines");
      const double* point=first->GetX();
      const double radius=std::hypot(std::hypot(point[0],point[1]),point[2]);
      if(!std::isfinite(radius)||radius<=0)return T::Status::Error(T::StatusCode::InvalidParticleState,
          "analytic trajectory launch radius must be finite and positive");
      if(minimumRadius<0||radius<minimumRadius)minimumRadius=radius;
    }
    if(minimumRadius<=0)return T::Status::Error(T::StatusCode::InvalidParticleState,"analytic trajectory launch requires a native field line");
    candidate->minimumLineRadiusM=minimumRadius;configuration=AnalyticalShockConfiguration(minimumRadius);
  }
  SEP::Shock::SourceEpochView view;
  const auto evaluated=SEP::Shock::EvaluateSourceAtEpoch(configuration,epochS,&view);
  if(!evaluated.ok())return evaluated;
  candidate->initialized=true;candidate->radiusM=view.trajectory.radiusM;
  candidate->trajectory=std::move(view.trajectory);
  return SEP::ParticleSource::ShockWave::Tenishev2005::ValidateAnalyticShockCache(*candidate);
}

SEP::Transport::Status PrepareAndCommitTrajectory(double epochS,bool mayInitialize) {
  namespace T=SEP::Transport;
  SEP::PhysicsCheckpoint::AnalyticShockCache candidate;T::Status local;
  std::vector<double> values,minimum,maximum;
  try {
    local=PrepareTrajectoryCandidate(epochS,mayInitialize,&candidate);
    if(local.ok()) {
      const auto& configuration=candidate.trajectory.configuration;
      if(configuration.knots.size()>(static_cast<std::size_t>(std::numeric_limits<int>::max())-5)/2)
        local=T::Status::Error(T::StatusCode::InvalidArgument,"analytic trajectory exceeds the MPI count representation");
      else {
        // Compare the complete governing trajectory, not merely one point at
        // which distinct knot tables could happen to intersect. All temporary
        // storage is allocated before collective readiness, so no rank can
        // publish a cache while another fails preparing its comparison data.
        values={configuration.launchEpochS,configuration.launchRadiusM,
                candidate.trajectory.epochS,candidate.radiusM,candidate.minimumLineRadiusM};
        for(const auto& knot:configuration.knots) {
          values.push_back(knot.radiusM);values.push_back(knot.speedMPerS);
        }
        minimum.resize(values.size());maximum.resize(values.size());
      }
    }
  }
  catch(const std::exception& error){local=T::Status::Error(T::StatusCode::InvalidParticleState,error.what());}
  const int bad=local.ok()?std::numeric_limits<int>::max():PIC::ThisThread;int first=0;
  MPI_Allreduce(&bad,&first,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  if(first!=std::numeric_limits<int>::max()) {
    int code=static_cast<int>(local.code);char message[2048]={};
    if(PIC::ThisThread==first)std::snprintf(message,sizeof(message),"analytic trajectory rank %d: %s",first,local.message.c_str());
    MPI_Bcast(&code,1,MPI_INT,first,MPI_GLOBAL_COMMUNICATOR);
    MPI_Bcast(message,static_cast<int>(sizeof(message)),MPI_CHAR,first,MPI_GLOBAL_COMMUNICATOR);
    return T::Status::Error(static_cast<T::StatusCode>(code),message);
  }
  const int count=static_cast<int>(values.size());int minCount=0,maxCount=0;
  MPI_Allreduce(&count,&minCount,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&count,&maxCount,1,MPI_INT,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  if(minCount!=maxCount)return T::Status::Error(T::StatusCode::InvalidParticleState,
      "MPI ranks disagree on analytic trajectory configuration size");
  MPI_Allreduce(values.data(),minimum.data(),count,MPI_DOUBLE,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(values.data(),maximum.data(),count,MPI_DOUBLE,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
  for(std::size_t i=0;i<values.size();++i)if(minimum[i]!=maximum[i])return T::Status::Error(T::StatusCode::InvalidParticleState,
      "MPI ranks disagree on the explicit analytic trajectory realization");
  // All validation/allocation precedes agreement. Commit only this prepared
  // native cache after every rank is ready and outside the particle read phase.
  SEP::ParticleSource::ShockWave::Tenishev2005::RestoreAnalyticShockCache(&candidate);
  return T::Status::Ok();
}

SEP::Transport::Status FindLocationAtRadius(double radiusM,int fieldLine,
    int* outputSegment,double* outputCoordinate,double outputPosition[3]) {
  namespace T=SEP::Transport;
  if(!outputSegment||!outputCoordinate||!outputPosition||fieldLine<0||fieldLine>=PIC::FieldLine::nFieldLine||
      !std::isfinite(radiusM)||radiusM<=0)
    return T::Status::Error(T::StatusCode::InvalidArgument,"analytic shock intersection requires a finite radius, native line and outputs");
  std::vector<SEP::Shock::Vector3> vertices;
  for(auto* segment=PIC::FieldLine::FieldLinesAll[fieldLine].GetFirstSegment();segment;segment=segment->GetNext()) {
    if(vertices.empty()) {
      const double* point=segment->GetBegin()->GetX();
      vertices.push_back({point[0],point[1],point[2]});
    }
    const double* point=segment->GetEnd()->GetX();vertices.push_back({point[0],point[1],point[2]});
  }
  if(vertices.empty())return T::Status::Error(T::StatusCode::InvalidParticleState,"analytic source line is empty");
  const auto intersection=SEP::Shock::IntersectSphere(vertices,radiusM,std::max(1.0,1e-12*radiusM),
      SEP::Shock::IntersectionPolicy::FirstOutward);
  if(intersection.status==SEP::Shock::IntersectionStatus::NoIntersection) {
    *outputSegment=-1;*outputCoordinate=-1;return T::Status::Ok();
  }
  if(intersection.status!=SEP::Shock::IntersectionStatus::Ok)
    return T::Status::Error(T::StatusCode::InvalidParticleState,"analytic source intersection rejected its native geometry: "+intersection.message);
  const auto& hit=intersection.intersections.front();
  *outputSegment=static_cast<int>(hit.segment);*outputCoordinate=static_cast<double>(hit.segment)+hit.fraction;
  outputPosition[0]=hit.positionM.x;outputPosition[1]=hit.positionM.y;outputPosition[2]=hit.positionM.z;
  return T::Status::Ok();
}

SEP::Transport::Status NativeSourceDeparture(double* departureEpochS) {
  namespace T=SEP::Transport;
  if(SEP::ShockModelType!=SEP::cShockModelType::Analytic1D||
      SEP::Background::ConfiguredProvider()!=SEP::Background::Provider::Analytic||
      SEP::Background::AnalyticRadial::DeclaredNativeWind().originM!=std::array<double,3>{{0,0,0}})
    return T::Status::Error(T::StatusCode::UnsupportedConfiguration,
        "native analytic source view requires the declared zero-centered analytic provider/shock owner");
  if(!SEP::ParticleSource::ShockWave::Tenishev2005::InitFlag)
    return T::Status::Error(T::StatusCode::InvalidParticleState,"native analytic source owner is unprepared");
  T::Birth::DriverWindow window;
  const auto frozen=T::PICAdapter::GetNativeDriverWindow(&window);
  if(!frozen.ok())return frozen;
  // The immutable host window supplies source time. A source-specific API must
  // not be usable as an arbitrary-time mathematical query or as a substitute
  // for opening the real driver window. No mutable PIC clock is consulted.
  *departureEpochS=window.startEpochS;return T::Status::Ok();
}
SEP::Transport::Status ValidateNativeSourceView(const SEP::Shock::SourceEpochView& view) {
  double departure=0;const auto frozen=NativeSourceDeparture(&departure);
  return frozen.ok()?SEP::Shock::ValidateSourceDepartureView(view,&gAnalyticalShockState,departure):frozen;
}

}  // namespace

SEP::Transport::Status SEP::ParticleSource::ShockWave::Tenishev2005::PrepareAnalyticTrajectoryAtEpoch(double epochS) {
  return PrepareAndCommitTrajectory(epochS,true);
}
SEP::Transport::Status SEP::ParticleSource::ShockWave::Tenishev2005::AdvanceAnalyticTrajectoryToEpoch(double epochS) {
  return PrepareAndCommitTrajectory(epochS,false);
}
SEP::Transport::Status SEP::ParticleSource::ShockWave::Tenishev2005::QueryAnalyticSourceAtEpoch(
    double epochS,SEP::Shock::SourceEpochView* output) {
  if(!output)return SEP::Transport::Status::Error(SEP::Transport::StatusCode::InvalidArgument,"native source view requires an output");
  double departure=0;const auto frozen=NativeSourceDeparture(&departure);if(!frozen.ok())return frozen;
  if(epochS!=departure)return SEP::Transport::Status::Error(SEP::Transport::StatusCode::InvalidParticleState,
      "native source query differs from the frozen driver departure");
  SEP::Shock::SourceEpochView candidate;
  const auto evaluated=SEP::Shock::EvaluateSourceAtEpoch(gAnalyticalShockState.configuration,epochS,&candidate);
  if(!evaluated.ok())return evaluated;
  const auto owned=SEP::Shock::ValidateSourceDepartureView(candidate,&gAnalyticalShockState,departure);
  if(!owned.ok())return owned;
  *output=std::move(candidate);return SEP::Transport::Status::Ok();
}
SEP::Transport::Status SEP::ParticleSource::ShockWave::Tenishev2005::InjectionLocationAt(
    const SEP::Shock::SourceEpochView& view,int line,int* segment,double* coordinate,double positionM[3]) {
  const auto valid=ValidateNativeSourceView(view);
  return valid.ok()?FindLocationAtRadius(view.trajectory.radiusM,line,segment,coordinate,positionM):valid;
}
SEP::Transport::ScalarResult SEP::ParticleSource::ShockWave::Tenishev2005::CompressionRatioAt(
    const SEP::Shock::SourceEpochView& view) {
  SEP::Transport::ScalarResult result;result.status=ValidateNativeSourceView(view);
  if(result.status.ok())result.value=CompressionAtRadius(view.trajectory.radiusM);
  return result;
}
SEP::Transport::ScalarResult SEP::ParticleSource::ShockWave::Tenishev2005::InjectionRateAt(
    const SEP::Shock::SourceEpochView& view) {
  SEP::Transport::ScalarResult result;result.status=ValidateNativeSourceView(view);
  if(!result.status.ok())return result;
  const double density=DensityAtRadius(view.trajectory.radiusM);
  if(!std::isfinite(density)||density<=0) {
    result.status=SEP::Transport::Status::Error(SEP::Transport::StatusCode::UnresolvedCoefficient,
        "analytic source density at the declared departure is not representable");return result;
  }
  result.value=SEP::FieldLine::FluxTubeGeometryCore::InjectedPhysicalParticleCount(
      SEP::Units::NumberDensityPerM3(density),SEP::Units::VolumeM3(1),
      SEP::FieldLine::InjectionParameters::InjectionEfficiency);
  if(!std::isfinite(result.value)||result.value<0)result.status=SEP::Transport::Status::Error(
      SEP::Transport::StatusCode::UnresolvedCoefficient,"analytic source departure rate is invalid");
  return result;
}

SEP::Transport::Status SEP::ParticleSource::ShockWave::Tenishev2005::CaptureAnalyticShockCache(
    SEP::PhysicsCheckpoint::AnalyticShockCache* cache) {
  if(!cache)return SEP::Transport::Status::Error(SEP::Transport::StatusCode::InvalidArgument,
      "null analytic shock checkpoint destination");
  cache->initialized=InitFlag;cache->radiusM=rShock;
  cache->minimumLineRadiusM=MinFieldLineHeliocentricDistance;
  if(InitFlag)cache->trajectory=gAnalyticalShockState;
  return ValidateAnalyticShockCache(*cache);
}
SEP::Transport::Status SEP::ParticleSource::ShockWave::Tenishev2005::ValidateAnalyticShockCache(
    const SEP::PhysicsCheckpoint::AnalyticShockCache& cache) {
  if(!std::isfinite(cache.radiusM)||cache.radiusM<0.0||
      !std::isfinite(cache.minimumLineRadiusM))
    return SEP::Transport::Status::Error(SEP::Transport::StatusCode::InvalidParticleState,
        "nonfinite analytic shock checkpoint cache");
  if(cache.initialized) {
    const auto valid=SEP::Shock::ValidateTrajectoryConfiguration(cache.trajectory.configuration);
    if(!valid.ok())return valid;
    if(!(cache.minimumLineRadiusM>0.0)||cache.radiusM!=cache.trajectory.radiusM||
        !std::isfinite(cache.trajectory.epochS))
      return SEP::Transport::Status::Error(SEP::Transport::StatusCode::InvalidParticleState,
          "analytic shock checkpoint radius/history mismatch");
  }
  return SEP::Transport::Status::Ok();
}
void SEP::ParticleSource::ShockWave::Tenishev2005::RestoreAnalyticShockCache(
    SEP::PhysicsCheckpoint::AnalyticShockCache* cache) {
  InitFlag=cache->initialized;rShock=cache->radiusM;
  MinFieldLineHeliocentricDistance=cache->minimumLineRadiusM;
  gAnalyticalShockState=std::move(cache->trajectory);
}


void SEP::ParticleSource::ShockWave::Tenishev2005::Init() {
  if(InitFlag)return;
  // Compatibility-only local initialization outside particle reads. The real
  // host uses the collective explicit-epoch owner before opening its window;
  // injection queries never initialize or mutate this cache.
  SEP::PhysicsCheckpoint::AnalyticShockCache candidate;
  const auto status=PrepareTrajectoryCandidate(0,true,&candidate);
  if(!status.ok())exit(__LINE__,__FILE__,status.message.c_str());
  RestoreAnalyticShockCache(&candidate);
}

double SEP::ParticleSource::ShockWave::Tenishev2005::GetCompressionRatio() {
  return CompressionAtRadius(rShock);
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
  if (SEP::ShockModelType==SEP::cShockModelType::SwCme1d) {
    // SWCME's prepared provider owns its radius/epoch. Do not initialize an
    // unrelated analytic trajectory or substitute a mutable PIC-clock query.
    rShock=SEP::SW1DAdapter::ShockRadiusM();
    return;
  }
  if(!InitFlag)Init();
  SEP::PhysicsCheckpoint::AnalyticShockCache candidate;
  const auto status=PrepareTrajectoryCandidate(SEP::Background::SimulationTimeSeconds(),false,&candidate);
  if(!status.ok())exit(__LINE__,__FILE__,status.message.c_str());
  RestoreAnalyticShockCache(&candidate);
}

double SEP::ParticleSource::ShockWave::Tenishev2005::GetSolarWindDensity() {
  return DensityAtRadius(rShock);
}

double SEP::ParticleSource::ShockWave::Tenishev2005::GetInjectionRate() {
  double r_sh,density;

  if(SEP::ShockModelType==SEP::cShockModelType::Analytic1D&&!InitFlag)Init();

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
  int segment=-1;const auto status=FindLocationAtRadius(radius,iFieldLine,&segment,&S,xInjection);
  if(!status.ok())exit(__LINE__,__FILE__,status.message.c_str());
  return segment;
}

int SEP::ParticleSource::ShockWave::Tenishev2005::GetCommittedInjectionLocation(int iFieldLine,double &S,double *xInjection) {
  // Read-only: the analytic trajectory is advanced (and rShock committed) only after the particle step
  // (AdvanceAnalyticTrajectoryToEpoch); advancing it here, inside the particle read phase, is a background write.
  // Before the trajectory is initialized there is nothing to read, and the updating query initializes it.
  if(SEP::ShockModelType!=SEP::cShockModelType::SwCme1d&&!InitFlag)return GetInjectionLocation(iFieldLine,S,xInjection);
  const double radius=SEP::ShockModelType==SEP::cShockModelType::SwCme1d
      ? SEP::SW1DAdapter::ShockRadiusM() : rShock;
  int segment=-1;const auto status=FindLocationAtRadius(radius,iFieldLine,&segment,&S,xInjection);
  if(!status.ok())exit(__LINE__,__FILE__,status.message.c_str());
  return segment;
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

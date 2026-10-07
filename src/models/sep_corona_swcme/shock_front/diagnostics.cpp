#include "diagnostics.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace SEP { namespace CoronaSwcme { namespace ShockFront { namespace {

bool Finite(CoronalCME::Vec3 value) {
  return std::isfinite(value.x)&&std::isfinite(value.y)&&std::isfinite(value.z);
}

struct SphereState {
  CoronalCME::Vec3 centerM;
  double radiusM = 0.0;
};

Core::Result<SphereState> Sphere(const Provider& provider,double timeS) {
  const auto trajectory=provider.Trajectory(timeS);
  if(!trajectory.ok())return Core::Result<SphereState>::Failure(
      trajectory.status.code,trajectory.status.message);
  const double sine=std::sin(provider.Event().halfWidthRad);
  const double center=trajectory.value.apexRadiusM/(1+sine);
  return Core::Result<SphereState>::Success(
      {center*provider.Event().direction,center*sine});
}

double SignedDistance(const Provider& provider,CoronalCME::Vec3 point,
    double timeS) {
  // The generating sphere supplies a smooth level set for temporal root
  // finding.  A zero of this full-sphere distance is only a *candidate*:
  // PassageAt re-queries the finite leading SSE ray support so the rear root
  // and a root beyond the tangent edge cannot be published as an arrival.
  const auto sphere=Sphere(provider,timeS);
  if(!sphere.ok())return std::numeric_limits<double>::quiet_NaN();
  return CoronalCME::Norm(point-sphere.value.centerM)-sphere.value.radiusM;
}

CoronalCME::Vec3 ObserverPosition(CoronalCME::Vec3 reference,
    CoronalCME::Vec3 velocity,double referenceTime,double time) {
  return reference+(time-referenceTime)*velocity;
}

bool SameRoot(const FrontIntersection& a,const FrontIntersection& b,
    double tolerance) {
  return CoronalCME::Norm(a.positionM-b.positionM)<=tolerance;
}

Core::Result<ObserverPassage> PassageAt(const Provider& provider,
    CoronalCME::Vec3 position,CoronalCME::Vec3 velocity,double time,
    double positionTolerance,bool tangent) {
  using Return=Core::Result<ObserverPassage>;
  const double radius=CoronalCME::Norm(position);
  if(radius==0)return Return::Failure(Core::StatusCode::OutOfDomain,
      "observer is at the HCI origin rather than on finite front support");
  const auto geometry=provider.EvaluateRay(position/radius,time,0);
  if(!geometry.ok())return Return::Failure(geometry.status.code,geometry.status.message);
  const double mismatch=CoronalCME::Norm(geometry.value.positionM-position);
  if(mismatch>positionTolerance)return Return::Failure(Core::StatusCode::OutOfDomain,
      "observer sphere root is outside the supported leading front");
  const auto record=provider.EvaluateFrontPoint(geometry.value.positionM,time,0);
  if(!record.ok())return Return::Failure(record.status.code,record.status.message);
  ObserverPassage out;out.timeS=time;out.positionM=position;out.front=record.value;
  // Relative normal speed, rather than apex speed or |V-vobs|, determines
  // whether the observer crosses transversely.  The front record still owns
  // the independent fast-shock decision: a geometric passage can therefore
  // be sub-fast/non-forward and is never promoted to a shock arrival here.
  out.relativeNormalSpeedMPerS=record.value.geometry.normalSpeedMPerS-
      CoronalCME::Dot(velocity,record.value.geometry.outwardNormal);
  out.kind=tangent||std::abs(out.relativeNormalSpeedMPerS)<=
      positionTolerance/std::max(1.0,std::abs(time))?
      ObserverPassageKind::Grazing:ObserverPassageKind::Transverse;
  return Return::Success(std::move(out));
}

} // namespace

const char* Name(ObserverPassageKind value) noexcept {
  return value==ObserverPassageKind::Transverse?"transverse":"grazing";
}
const char* Name(EnsembleOutcome value) noexcept {
  switch(value) {
    case EnsembleOutcome::KnownAcceptedShock:return "known-accepted-shock";
    case EnsembleOutcome::KnownNoShock:return "known-no-shock";
    case EnsembleOutcome::NumericalUnknown:return "numerical-unknown";
  }
  return "unknown";
}

Core::Result<std::vector<FrontIntersection>> IntersectPolylineWithFront(
    const Provider& provider,const std::vector<CoronalCME::Vec3>& points,
    double epoch,double tolerance) {
  using Return=Core::Result<std::vector<FrontIntersection>>;
  if(points.size()<2||!std::isfinite(epoch)||!std::isfinite(tolerance)||
      tolerance<=0||!std::all_of(points.begin(),points.end(),Finite))
    return Return::Failure(Core::StatusCode::InvalidConfiguration,
        "front intersection requires a finite polyline and positive distance tolerance");
  const auto sphere=Sphere(provider,epoch);
  if(!sphere.ok())return Return::Failure(sphere.status.code,sphere.status.message);
  std::vector<FrontIntersection> out;double prefix=0;
  for(std::size_t segment=0;segment+1<points.size();++segment) {
    const CoronalCME::Vec3 delta=points[segment+1]-points[segment];
    const double length=CoronalCME::Norm(delta);
    if(length==0)continue;
    const CoronalCME::Vec3 offset=points[segment]-sphere.value.centerM;
    // Solve |p+s*d-C|^2=a_sphere^2 in the segment parameter s.  Arc length is
    // retained separately because branch ordering along a traced field line is
    // physical, whereas the generating-sphere quadratic has no such ordering.
    const double a=CoronalCME::Dot(delta,delta);
    const double b=2*CoronalCME::Dot(offset,delta);
    const double c=CoronalCME::Dot(offset,offset)-
        sphere.value.radiusM*sphere.value.radiusM;
    double discriminant=b*b-4*a*c;
    const double discScale=std::max({b*b,std::abs(4*a*c),1.0});
    if(discriminant<-64*std::numeric_limits<double>::epsilon()*discScale) {
      prefix+=length;continue;
    }
    const bool tangent=std::abs(discriminant)<=
        64*std::numeric_limits<double>::epsilon()*discScale;
    // At a double root, carrying a roundoff-sized positive discriminant into
    // sqrt creates an O(sqrt(epsilon))*radius displacement (metres to tens of
    // metres at coronal scale).  Collapse the same symmetric roundoff band
    // used to identify tangency before taking the square root.
    discriminant=tangent?0.0:std::max(0.0,discriminant);
    const double root=std::sqrt(discriminant);
    const double fractions[2]={(-b-root)/(2*a),(-b+root)/(2*a)};
    const int count=tangent?1:2;
    for(int index=0;index<count;++index) {
      double fraction=fractions[index];
      const double fractionTolerance=tolerance/length;
      if(fraction<-fractionTolerance||fraction>1+fractionTolerance)continue;
      fraction=std::max(0.0,std::min(1.0,fraction));
      const CoronalCME::Vec3 point=points[segment]+fraction*delta;
      const double pointRadius=CoronalCME::Norm(point);if(pointRadius==0)continue;
      const auto geometry=provider.EvaluateRay(point/pointRadius,epoch,
          static_cast<std::uint64_t>(out.size()+1));
      if(!geometry.ok())continue; // full-sphere rear/support misses are expected
      if(CoronalCME::Norm(geometry.value.positionM-point)>tolerance)continue;
      const auto record=provider.EvaluateFrontPoint(geometry.value.positionM,epoch,
          static_cast<std::uint64_t>(out.size()+1));
      if(!record.ok())return Return::Failure(record.status.code,record.status.message);
      FrontIntersection hit;hit.segmentIndex=segment;
      hit.segmentFraction=fraction;hit.arcLengthM=prefix+fraction*length;
      hit.positionM=geometry.value.positionM;hit.tangent=tangent;hit.front=record.value;
      hit.branchId=static_cast<std::uint64_t>(segment+1)<<32|
          static_cast<std::uint64_t>(index+1);
      const double along=CoronalCME::Dot(record.value.upstream.magneticFieldT,
          delta/length);
      hit.signedPolarityAlongTrace=along>0?1:(along<0?-1:0);
      if(out.empty()||!SameRoot(out.back(),hit,tolerance))out.push_back(std::move(hit));
    }
    prefix+=length;
  }
  std::sort(out.begin(),out.end(),[](const auto& a,const auto& b){
    return a.arcLengthM<b.arcLengthM;});
  return Return::Success(std::move(out));
}

Core::Result<std::vector<ObserverPassage>> FindObserverPassages(
    const Provider& provider,CoronalCME::Vec3 reference,
    CoronalCME::Vec3 velocity,double referenceTime,double begin,double end,
    double scanStep,double timeTolerance,double distanceTolerance) {
  using Return=Core::Result<std::vector<ObserverPassage>>;
  if(!Finite(reference)||!Finite(velocity)||!std::isfinite(referenceTime)||
      !std::isfinite(begin)||!std::isfinite(end)||begin>=end||
      !std::isfinite(scanStep)||scanStep<=0||!std::isfinite(timeTolerance)||
      timeTolerance<=0||!std::isfinite(distanceTolerance)||distanceTolerance<=0)
    return Return::Failure(Core::StatusCode::InvalidConfiguration,
        "observer search requires finite ordered times and positive tolerances");
  auto value=[&](double time){return SignedDistance(provider,
      ObserverPosition(reference,velocity,referenceTime,time),time);};
  std::vector<double> times;std::vector<double> values;
  const int intervals=std::max(1,static_cast<int>(std::ceil((end-begin)/scanStep)));
  for(int i=0;i<=intervals;++i) {
    const double time=begin+(end-begin)*i/intervals;
    const double g=value(time);if(!std::isfinite(g))return Return::Failure(
        Core::StatusCode::NumericalFailure,"observer level-set evaluation failed");
    times.push_back(time);values.push_back(g);
  }
  std::vector<std::pair<double,bool>> roots;
  for(int i=0;i<intervals;++i) {
    if(values[i]==0)roots.push_back({times[i],false});
    if(values[i]*values[i+1]<0) {
      double left=times[i],right=times[i+1];double fleft=values[i];
      while(right-left>timeTolerance) {
        const double middle=0.5*(left+right);const double fm=value(middle);
        if(fleft*fm<=0)right=middle;else {left=middle;fleft=fm;}
      }
      roots.push_back({0.5*(left+right),false});
    }
  }
  // Minimize |g| around every sampled local minimum.  This catches a double
  // root without weakening the signed crossing logic or inventing a hit from
  // a merely close miss: the final dimensional distance tolerance still gates.
  for(int i=1;i<intervals;++i)if(std::abs(values[i])<=std::abs(values[i-1])&&
      std::abs(values[i])<=std::abs(values[i+1])) {
    double left=times[i-1],right=times[i+1];
    for(int iteration=0;iteration<100&&right-left>timeTolerance;++iteration) {
      const double third=(right-left)/3;
      const double x1=left+third,x2=right-third;
      if(std::abs(value(x1))<std::abs(value(x2)))right=x2;else left=x1;
    }
    const double time=0.5*(left+right);
    if(std::abs(value(time))<=distanceTolerance)roots.push_back({time,true});
  }
  std::sort(roots.begin(),roots.end());
  std::vector<ObserverPassage> out;
  for(const auto& root:roots) {
    if(!out.empty()&&std::abs(out.back().timeS-root.first)<=2*timeTolerance)continue;
    const auto point=ObserverPosition(reference,velocity,referenceTime,root.first);
    const auto passage=PassageAt(provider,point,velocity,root.first,
        2*distanceTolerance,root.second);
    if(passage.ok())out.push_back(passage.value);
    else if(passage.status.code!=Core::StatusCode::OutOfDomain)
      return Return::Failure(passage.status.code,passage.status.message);
  }
  return Return::Success(std::move(out));
}

Core::Result<IncidentParticleFlux> EvaluateIncidentParticleFluxAtApexRadius(
    const Provider& provider,double requestedRadius) {
  using Return=Core::Result<IncidentParticleFlux>;
  if(!std::isfinite(requestedRadius)||requestedRadius<=0)
    return Return::Failure(Core::StatusCode::InvalidConfiguration,
        "incident-flux normalization radius must be finite and positive");

  const auto& event=provider.Event();
  double left=event.ambient.support.startS;
  double right=event.ambient.support.endS;
  const auto first=provider.Trajectory(left);
  const auto last=provider.Trajectory(right);
  if(!first.ok())return Return::Failure(first.status.code,first.status.message);
  if(!last.ok())return Return::Failure(last.status.code,last.status.message);
  if(last.value.apexRadiusM<first.value.apexRadiusM)
    return Return::Failure(Core::StatusCode::InvalidState,
        "incident-flux normalization requires a monotonically outward apex history");
  const double radialScale=std::max({1.0,std::abs(first.value.apexRadiusM),
      std::abs(last.value.apexRadiusM),std::abs(requestedRadius)});
  const double radialTolerance=128*std::numeric_limits<double>::epsilon()*radialScale;
  if(requestedRadius<first.value.apexRadiusM-radialTolerance||
      requestedRadius>last.value.apexRadiusM+radialTolerance)
    return Return::Failure(Core::StatusCode::OutOfDomain,
        "incident-flux normalization radius is outside the front history");

  // Bisection uses only the analytical trajectory, while the final flux uses
  // the full curved surface.  Stop on a dimensional radius tolerance rather
  // than a fixed iteration count alone so the root error is meaningful from
  // the low corona through AU scales.
  for(int iteration=0;iteration<256&&right-left>
      16*std::numeric_limits<double>::epsilon()*
      std::max({1.0,std::abs(left),std::abs(right)});++iteration) {
    const double middle=0.5*(left+right);
    const auto state=provider.Trajectory(middle);
    if(!state.ok())return Return::Failure(state.status.code,state.status.message);
    if(state.value.apexRadiusM<requestedRadius)left=middle;else right=middle;
    if(std::abs(state.value.apexRadiusM-requestedRadius)<=radialTolerance) {
      left=right=middle;
      break;
    }
  }
  const double epoch=0.5*(left+right);
  const auto surface=provider.EvaluateEpoch(epoch,1);
  if(!surface.ok())return Return::Failure(surface.status.code,
      "cannot evaluate incident-flux surface: "+surface.status.message);

  long double proton=0,electron=0;
  for(const ShockRecord& record:surface.value->records) {
    if(record.status!=FrontStatus::SolvedFastShock)continue;
    if(!std::isfinite(record.inflowMPerS)||record.inflowMPerS<=0||
        !std::isfinite(record.geometry.areaM2)||record.geometry.areaM2<=0||
        !std::isfinite(record.upstream.plasma.protonNumberDensityM3)||
        !std::isfinite(record.upstream.plasma.electronNumberDensityM3)||
        record.upstream.plasma.protonNumberDensityM3<0||
        record.upstream.plasma.electronNumberDensityM3<0)
      return Return::Failure(Core::StatusCode::InvalidState,
          "accepted shock record has an invalid upstream incident flux");
    const long double measure=static_cast<long double>(record.inflowMPerS)*
        static_cast<long double>(record.geometry.areaM2);
    proton+=measure*record.upstream.plasma.protonNumberDensityM3;
    electron+=measure*record.upstream.plasma.electronNumberDensityM3;
  }
  IncidentParticleFlux out;
  out.epochS=epoch;
  out.apexRadiusM=surface.value->trajectory.apexRadiusM;
  out.acceptedAreaM2=surface.value->area.acceptedShockM2;
  out.excludedPhysicalAreaM2=surface.value->area.geometricSupportM2-
      surface.value->area.acceptedShockM2-
      surface.value->area.numericalFailureM2;
  out.numericalFailureAreaM2=surface.value->area.numericalFailureM2;
  out.protonRatePerS=static_cast<double>(proton);
  out.electronRatePerS=static_cast<double>(electron);
  out.alphaRatePerS=event.ambient.composition.alphaToProtonNumberRatio*
      out.protonRatePerS;
  if(out.numericalFailureAreaM2>0)
    return Return::Failure(Core::StatusCode::NumericalFailure,
        "incident-flux normalization surface contains numerically unresolved area");
  if(!std::isfinite(out.protonRatePerS)||!std::isfinite(out.electronRatePerS)||
      !std::isfinite(out.alphaRatePerS)||out.protonRatePerS<=0||
      out.electronRatePerS<=0||out.alphaRatePerS<0||out.acceptedAreaM2<=0)
    return Return::Failure(Core::StatusCode::InvalidState,
        "incident-flux normalization has no positive accepted-shock source");
  return Return::Success(out);
}

Core::Result<ReducedRestartState> MakeRestartState(const Provider& provider) {
  using Return=Core::Result<ReducedRestartState>;
  const auto current=provider.Current();if(!current)return Return::Failure(
      Core::StatusCode::InvalidState,"cannot checkpoint an uncommitted provider");
  ReducedRestartState out;out.eventIdentity=current->eventIdentity;
  out.epochS=current->trajectory.timeS;out.generation=current->generation;
  out.phase=current->trajectory.phase;out.handoffTimeS=current->handoff.timeS;
  out.handoffRadiusM=current->handoff.apexRadiusM;
  out.handoffSpeedMPerS=current->handoff.apexSpeedMPerS;
  return Return::Success(std::move(out));
}

Core::Result<std::shared_ptr<const Epoch>> RestoreRestartState(
    Provider* provider,const ReducedRestartState& state) {
  using Return=Core::Result<std::shared_ptr<const Epoch>>;
  if(!provider)return Return::Failure(Core::StatusCode::InvalidConfiguration,
      "restart target provider is null");
  // Restarts store no duplicate surface or ambient arrays: this analytical
  // model reconstructs them from the clock.  Reconstruction is safe only when
  // the complete event/asset identity and the phase-defining handoff receipt
  // match.  Validate all of that before Prepare can replace Current().
  if(state.eventIdentity!=provider->Event().physicsFingerprint)
    return Return::Failure(Core::StatusCode::DataIntegrityFailure,
        "restart event/asset fingerprint differs from the resolved provider");
  if(state.generation==0||!std::isfinite(state.epochS))return Return::Failure(
      Core::StatusCode::InvalidState,"restart clock/generation is invalid");
  const auto handoff=provider->Handoff();if(!handoff.ok())return Return::Failure(
      handoff.status.code,handoff.status.message);
  const double scale=std::max({1.0,std::abs(handoff.value.apexRadiusM),
      std::abs(handoff.value.apexSpeedMPerS)});
  if(std::abs(state.handoffTimeS-handoff.value.timeS)>1e-8||
      std::abs(state.handoffRadiusM-handoff.value.apexRadiusM)>1e-12*scale||
      std::abs(state.handoffSpeedMPerS-handoff.value.apexSpeedMPerS)>1e-12*scale)
    return Return::Failure(Core::StatusCode::DataIntegrityFailure,
        "restart handoff authority differs from the resolved provider");
  const auto trajectory=provider->Trajectory(state.epochS);
  if(!trajectory.ok()||trajectory.value.phase!=state.phase)return Return::Failure(
      Core::StatusCode::DataIntegrityFailure,"restart phase/clock is inconsistent");
  return provider->Prepare(state.epochS,state.generation);
}

std::string SerializeEpochJson(const Epoch& epoch) {
  std::size_t solved=0,subfast=0,nonforward=0,numerical=0;
  for(const auto& record:epoch.records) {
    solved+=record.status==FrontStatus::SolvedFastShock;
    subfast+=record.status==FrontStatus::SubfastFront;
    nonforward+=record.status==FrontStatus::NonForwardInflow;
    numerical+=record.status==FrontStatus::AmbientUnavailable||
      record.status==FrontStatus::NumericallyUnresolvedWeakShock||
      record.status==FrontStatus::WrongBranch||record.status==FrontStatus::InvalidJump;
  }
  // Areas are absolute SI measures over the complete geometric support.  Do
  // not normalize away sub-fast or numerical-unknown patches: consumers need
  // both the physical no-shock fraction and unresolved coverage to interpret
  // any connectivity or observer statistic.
  std::ostringstream out;out<<std::setprecision(17);
  out<<"{\"event_identity\":\""<<epoch.eventIdentity<<"\",\"time_s\":"
     <<epoch.trajectory.timeS<<",\"generation\":"<<epoch.generation
     <<",\"ambient_generation\":"<<epoch.ambientGeneration
     <<",\"phase\":\""<<Name(epoch.trajectory.phase)<<"\",\"area_m2\":{"
     <<"\"geometric\":"<<epoch.area.geometricSupportM2
     <<",\"superfast\":"<<epoch.area.superfastCandidateM2
     <<",\"accepted\":"<<epoch.area.acceptedShockM2
     <<",\"numerical_unknown\":"<<epoch.area.numericalFailureM2
     <<"},\"counts\":{\"solved\":"<<solved<<",\"subfast\":"<<subfast
     <<",\"nonforward\":"<<nonforward<<",\"numerical_unknown\":"<<numerical
     <<"},\"mesh\":{\"vertices\":"<<epoch.vertices.size()
     <<",\"triangles\":"<<epoch.triangles.size()
     <<"},\"geometric_endpoint_reached\":"
     <<(epoch.geometricEndpointReached?"true":"false")
     <<",\"apex_shock_accepted\":"<<(epoch.apexShockAccepted?"true":"false")
     <<"}";
  return out.str();
}

std::string SerializeSurfaceCsv(const Epoch& epoch) {
  std::ostringstream out;out<<std::setprecision(17);
  out<<"time_s,generation,stable_id,vertex0,vertex1,vertex2,x_m,y_m,z_m,"
      "nx,ny,nz,vn_m_s,curved_area_m2,planar_area_m2,status,"
      "rho1_kg_m3,p1_pa,u1x_m_s,u1y_m_s,u1z_m_s,b1x_t,b1y_t,b1z_t,mf,"
      "downstream_valid,compression,rho2_kg_m3,p2_pa,cb_valid,cb,ht_status\n";
  for(std::size_t face=0;face<epoch.records.size();++face) {
    const auto& record=epoch.records[face];
    const SurfaceTriangle* triangle=face<epoch.triangles.size()?
        &epoch.triangles[face]:nullptr;
    out<<epoch.trajectory.timeS<<','<<epoch.generation<<','
       <<record.geometry.stableId<<',';
    if(triangle)out<<triangle->vertex[0]<<','<<triangle->vertex[1]<<','
       <<triangle->vertex[2]<<',';
    else out<<",,,";
    out<<record.geometry.positionM.x<<','
       <<record.geometry.positionM.y<<','<<record.geometry.positionM.z<<','
       <<record.geometry.outwardNormal.x<<','<<record.geometry.outwardNormal.y<<','
       <<record.geometry.outwardNormal.z<<','<<record.geometry.normalSpeedMPerS<<','
       <<record.geometry.areaM2<<','
       <<(triangle?triangle->planarAreaM2:0)<<','<<Name(record.status)<<',';
    if(record.status==FrontStatus::BelowPhysicalInnerBoundary||
       record.status==FrontStatus::AmbientUnavailable)out<<",,,,,,,,";
    else out<<record.upstream.plasma.massDensityKgM3<<','
      <<record.upstream.plasma.pressurePa<<','<<record.upstream.velocityMPerS.x<<','
      <<record.upstream.velocityMPerS.y<<','<<record.upstream.velocityMPerS.z<<','
      <<record.upstream.magneticFieldT.x<<','<<record.upstream.magneticFieldT.y<<','
      <<record.upstream.magneticFieldT.z<<','<<record.fastMach;
    out<<','<<(record.downstreamValid?"true":"false")<<',';
    if(record.downstreamValid)out<<record.jump.compressionRatio<<','
      <<record.jump.downstream.massDensityKgM3<<','<<record.jump.downstream.pressurePa;
    else out<<",,";
    out<<','<<(record.diagnostics.magneticCompressionValid?"true":"false")<<',';
    if(record.diagnostics.magneticCompressionValid)
      out<<record.diagnostics.magneticCompression;
    out<<','<<Name(record.diagnostics.htStatus)<<'\n';
  }
  return out.str();
}

ConnectedBranchDerivative DifferentiateConnectedBranch(
    const ConnectedBranchSample& previous,const ConnectedBranchSample& current,
    const ConnectedBranchSample& next) {
  ConnectedBranchDerivative out;
  if(!previous.valid||!current.valid||!next.valid) {
    out.reason="connected derivative requires three valid branch samples";return out;
  }
  if(previous.branchId!=current.branchId||next.branchId!=current.branchId) {
    out.reason="branch identity changed; derivative is absent at merger/graze/selection";
    return out;
  }
  if(!(previous.timeS<current.timeS&&current.timeS<next.timeS)) {
    out.reason="connected derivative times are not strictly ordered";return out;
  }
  // Positions are already evaluated at the moving intersection.  Their
  // centered difference therefore contains both explicit surface evolution
  // and motion of the root along the supplied polyline (the xdot.grad term).
  // Refusing a changed branch is essential: differencing across a tangency or
  // merger would manufacture an arbitrarily large connection velocity.
  const double span=next.timeS-previous.timeS;
  out.intersectionVelocityMPerS=(next.positionM-previous.positionM)/span;
  out.obliquityRateRadPerS=(next.obliquityRad-previous.obliquityRad)/span;
  out.fastMachRatePerS=(next.fastMach-previous.fastMach)/span;
  out.valid=Finite(out.intersectionVelocityMPerS)&&
      std::isfinite(out.obliquityRateRadPerS)&&std::isfinite(out.fastMachRatePerS);
  if(!out.valid)out.reason="connected derivative produced a nonfinite result";
  return out;
}

Core::Result<EnsembleLedger> SummarizeEnsemble(
    const std::vector<EnsembleMember>& members,double tolerance) {
  using Return=Core::Result<EnsembleLedger>;
  if(members.empty()||!std::isfinite(tolerance)||tolerance<=0)
    return Return::Failure(Core::StatusCode::InvalidConfiguration,
        "ensemble requires members and a positive normalization tolerance");
  EnsembleLedger out;const std::string label=members.front().commonLabel;
  if(label.empty())return Return::Failure(Core::StatusCode::InvalidConfiguration,
      "ensemble common label is empty");
  for(const auto& member:members) {
    if(member.commonLabel!=label||!std::isfinite(member.weight)||member.weight<0)
      return Return::Failure(Core::StatusCode::InvalidConfiguration,
          "ensemble labels differ or a weight is invalid");
    switch(member.outcome) {
      case EnsembleOutcome::KnownAcceptedShock:out.acceptedWeight+=member.weight;break;
      case EnsembleOutcome::KnownNoShock:out.noShockWeight+=member.weight;break;
      case EnsembleOutcome::NumericalUnknown:out.unknownWeight+=member.weight;break;
    }
  }
  // Unknown numerical weight is neither discarded nor guessed.  It widens a
  // rigorous probability interval [accepted, accepted+unknown]; known
  // no-shock weight stays in the common declared measure rather than causing
  // survivor renormalization.
  const double total=out.acceptedWeight+out.noShockWeight+out.unknownWeight;
  if(std::abs(total-1)>tolerance)return Return::Failure(
      Core::StatusCode::InvalidConfiguration,
      "ensemble weights are not normalized under the declared measure");
  out.acceptedProbabilityLower=out.acceptedWeight;
  out.acceptedProbabilityUpper=out.acceptedWeight+out.unknownWeight;
  return Return::Success(out);
}

} } } // namespace SEP::CoronaSwcme::ShockFront

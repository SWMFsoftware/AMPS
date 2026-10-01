#include "sep_coronal_cme/research_extensions.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <tuple>

namespace SEP { namespace CoronalCME { namespace Research {
namespace {
bool finite(Vec3 v) { return std::isfinite(v.x)&&std::isfinite(v.y)&&std::isfinite(v.z); }
template<class T> Core::Result<T> bad(const std::string& s) {
  return Core::Result<T>::Failure(Core::StatusCode::InvalidConfiguration,s);
}
double energy(double mass, Vec3 p) {
  constexpr double c=299792458.0;
  return std::hypot(mass*c*c,Norm(p)*c);
}
double kinetic(double mass, Vec3 p) {
  constexpr double c=299792458.0;
  const double pc=Norm(p)*c;
  return pc*(pc/(energy(mass,p)+mass*c*c));
}
}
Core::Status ValidateCapability(const CapabilityIdentity& i) {
  if (i.schemaMajor!=ResearchSchemaMajor||i.capabilityId.empty()||i.algorithmVersion.empty()||
      i.coefficientFingerprint.empty()||!i.backgroundGeneration)
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,"incomplete research capability identity");
  return Core::Status::Success();
}
Core::Result<ForeshockEvaluation> EvaluateForeshockProxy(
    const ForeshockProxy& p,double d,std::uint64_t generation,double ambient) {
  const auto id=ValidateCapability(p.identity);
  if (!id.ok()||!std::isfinite(d)||!std::isfinite(ambient)||ambient<=0)
    return bad<ForeshockEvaluation>("invalid proxy identity/coefficient/distance");
  ForeshockEvaluation out; out.meanFreePathM=ambient; out.upstream=d>=0; out.shockGeneration=generation;
  if (!p.identity.enabled) return Core::Result<ForeshockEvaluation>::Success(out);
  if (!generation||generation!=p.identity.shockGeneration||p.referenceDistanceM<=0||
      p.supportWidthM<=0||!std::isfinite(p.referenceDistanceM)||!std::isfinite(p.supportWidthM)||
      !std::isfinite(p.minimumFactor)||p.minimumFactor<=0||p.minimumFactor>1)
    return bad<ForeshockEvaluation>("proxy generation/support/reduction range mismatch");
  if (d<p.referenceDistanceM)
    return bad<ForeshockEvaluation>("unresolved front-to-reference/downstream proxy evaluation");
  const double x=(d-p.referenceDistanceM)/p.supportWidthM;
  // A compact C2 quintic recovery has zero derivatives at either endpoint.
  // No wave-energy object is touched: this is a coefficient sensitivity only.
  const double smooth=x>=1?1:x*x*x*(10+x*(-15+6*x));
  out.factor=p.minimumFactor+(1-p.minimumFactor)*smooth;
  out.meanFreePathM=ambient*out.factor;
  return Core::Result<ForeshockEvaluation>::Success(out);
}
Core::Result<DynamicEllipsoid> DynamicEllipsoid::Build(const DynamicGeometryState& s) {
  if (!finite(s.centerM)||!finite(s.centerVelocityMPerS)||!finite(s.rotationAxis)||
      !std::isfinite(s.angleRad)||!std::isfinite(s.angularRateRadPerS)||
      !std::isfinite(s.epochS)||s.frameId.empty()||(!s.fullAttitude&&(!(Norm(s.rotationAxis)>0)||!std::isfinite(Norm(s.rotationAxis)))))
    return bad<DynamicEllipsoid>("invalid center/attitude/frame history state");
  for (int j=0;j<3;j++) if (!(s.axesM[j]>0)||!std::isfinite(s.axesM[j])||!std::isfinite(s.axisRatesMPerS[j]))
    return bad<DynamicEllipsoid>("invalid dynamic ellipsoid axes");
  if(s.fullAttitude) {
    if(!finite(s.inertialAngularVelocityPerS)) return bad<DynamicEllipsoid>("nonfinite angular velocity");
    const auto& r=s.bodyToInertial;
    for(double v:r) if(!std::isfinite(v)) return bad<DynamicEllipsoid>("nonfinite attitude matrix");
    Vec3 a{r[0],r[3],r[6]},b{r[1],r[4],r[7]},c{r[2],r[5],r[8]};
    if(std::abs(Dot(a,a)-1)>1e-10||std::abs(Dot(b,b)-1)>1e-10||std::abs(Dot(c,c)-1)>1e-10||
       std::abs(Dot(a,b))>1e-10||std::abs(Dot(a,c))>1e-10||std::abs(Dot(b,c))>1e-10||
       std::abs(Dot(a,Cross(b,c))-1)>1e-10) return bad<DynamicEllipsoid>("attitude is not a proper SO(3) rotation");
  }
  DynamicEllipsoid out; out.state_=s;
  if(!s.fullAttitude) out.state_.rotationAxis=Unit(s.rotationAxis);
  return Core::Result<DynamicEllipsoid>::Success(out);
}
Vec3 DynamicEllipsoid::Rotate(Vec3 b) const {
  if(state_.fullAttitude) {const auto& r=state_.bodyToInertial;return {r[0]*b.x+r[1]*b.y+r[2]*b.z,r[3]*b.x+r[4]*b.y+r[5]*b.z,r[6]*b.x+r[7]*b.y+r[8]*b.z};}
  const Vec3 n=state_.rotationAxis; const double c=std::cos(state_.angleRad),s=std::sin(state_.angleRad);
  return c*b+s*Cross(n,b)+(1-c)*Dot(n,b)*n;
}
Vec3 DynamicEllipsoid::InverseRotate(Vec3 b) const {
  if(state_.fullAttitude) {const auto& r=state_.bodyToInertial;return {r[0]*b.x+r[3]*b.y+r[6]*b.z,r[1]*b.x+r[4]*b.y+r[7]*b.z,r[2]*b.x+r[5]*b.y+r[8]*b.z};}
  const Vec3 n=state_.rotationAxis; const double c=std::cos(state_.angleRad),s=std::sin(state_.angleRad);
  return c*b-s*Cross(n,b)+(1-c)*Dot(n,b)*n;
}
Vec3 DynamicEllipsoid::Omega() const {
  return state_.fullAttitude?state_.inertialAngularVelocityPerS:state_.angularRateRadPerS*state_.rotationAxis;
}
Vec3 DynamicEllipsoid::Point(double t,double p) const {
  return state_.centerM+Rotate({state_.axesM[0]*std::cos(t),state_.axesM[1]*std::sin(t)*std::cos(p),
                              state_.axesM[2]*std::sin(t)*std::sin(p)});
}
Vec3 DynamicEllipsoid::SurfaceVelocity(double t,double p) const {
  const Vec3 relative=Point(t,p)-state_.centerM;
  return state_.centerVelocityMPerS+Cross(Omega(),relative)+
      Rotate({state_.axisRatesMPerS[0]*std::cos(t),state_.axisRatesMPerS[1]*std::sin(t)*std::cos(p),
              state_.axisRatesMPerS[2]*std::sin(t)*std::sin(p)});
}
Core::Result<SurfaceEvaluation> DynamicEllipsoid::Evaluate(Vec3 x) const {
  if (!finite(x)) return bad<SurfaceEvaluation>("nonfinite dynamic surface point");
  const Vec3 body=InverseRotate(x-state_.centerM);
  const Vec3 g={body.x/(state_.axesM[0]*state_.axesM[0]),body.y/(state_.axesM[1]*state_.axesM[1]),body.z/(state_.axesM[2]*state_.axesM[2])};
  if (Norm(g)==0) return bad<SurfaceEvaluation>("ellipsoid center has no surface normal");
  SurfaceEvaluation out; out.implicitValue=Dot(body,g)-1; out.outwardNormal=Unit(Rotate(g));
  const Vec3 expand=Rotate({body.x*state_.axisRatesMPerS[0]/state_.axesM[0],
      body.y*state_.axisRatesMPerS[1]/state_.axesM[1],body.z*state_.axisRatesMPerS[2]/state_.axesM[2]});
  // This is -F_t/|grad F|, including translation and every Q_dot term.
  out.normalSpeedMPerS=Dot(out.outwardNormal,state_.centerVelocityMPerS+
      Cross(Omega(),x-state_.centerM)+expand);
  const double a=1/(state_.axesM[0]*state_.axesM[0]),b=1/(state_.axesM[1]*state_.axesM[1]),c=1/(state_.axesM[2]*state_.axesM[2]);
  out.meanCurvaturePerM=(Dot(g,g)*(a+b+c)-(a*g.x*g.x+b*g.y*g.y+c*g.z*g.z))/(2*std::pow(Norm(g),3));
  out.gaussianCurvaturePerM2=a*b*c*Dot(body,g)/std::pow(Norm(g),4);
  return Core::Result<SurfaceEvaluation>::Success(out);
}
Core::Result<ReferenceFamilyMember> SolveReferenceFamilyMember(const ReferenceStratum& s,
    double target,double upper,double tol,const std::function<double(double)>& ratio,
    const std::function<Core::Result<double>(double)>& geometry) {
  if (s.speciesId.empty()||!(s.momentumLowerSi>=0)||!(s.momentumUpperSi>s.momentumLowerSi)||
      !(s.timeUpperS>s.timeLowerS)||!s.geometryGeneration||!(target>0)||!(upper>0)||!(tol>0)||
      !std::isfinite(s.momentumLowerSi)||!std::isfinite(s.momentumUpperSi)||!std::isfinite(s.timeLowerS)||!std::isfinite(s.timeUpperS)||
      !std::isfinite(target)||!std::isfinite(upper)||!std::isfinite(tol)||!ratio||!geometry)
    return bad<ReferenceFamilyMember>("invalid family stratum/target/support");
  auto integral=[&](double length)->Core::Result<double> {
    double previous=0;
    for (int n=32;n<=131072;n*=2) {
      double sum=0;
      for (int j=0;j<n;j++) {
        const double f=ratio(length*(j+0.5)/n);
        if (!(f>0)||!std::isfinite(f)) return bad<double>("nonpositive/nonfinite inflow/diffusivity");
        sum+=f*length/n;
      }
      if (n>32&&std::abs(sum-previous)<=tol*std::max(target,std::abs(sum))) return Core::Result<double>::Success(sum);
      previous=sum;
    }
    return bad<double>("Peclet quadrature did not converge");
  };
  auto full=integral(upper); if (!full.ok()) return bad<ReferenceFamilyMember>(full.status.message);
  if (full.value<target) return bad<ReferenceFamilyMember>("Peclet root outside geometry-certified support");
  double lo=0,hi=upper,depth=0;
  for (int i=0;i<80;i++) {
    const double mid=(lo+hi)/2; auto value=integral(mid);
    if (!value.ok()) return bad<ReferenceFamilyMember>(value.status.message);
    depth=value.value;
    if (depth<target) lo=mid; else hi=mid;
    if (std::abs(depth-target)<=tol*target) { lo=hi=mid; break; }
  }
  const double offset=(lo+hi)/2;
  auto measured=geometry(offset);
  if (!measured.ok()||!(measured.value>0)||!std::isfinite(measured.value))
    return bad<ReferenceFamilyMember>("family geometry/mask/clearance/joint measure rejected");
  auto exact=integral(offset); if (!exact.ok()) return bad<ReferenceFamilyMember>(exact.status.message);
  if(std::abs(exact.value-target)>tol*target) return bad<ReferenceFamilyMember>("Peclet root tolerance not achieved");
  return Core::Result<ReferenceFamilyMember>::Success({s,offset,exact.value,measured.value,false});
}
Core::Status ValidateReferenceFamily(const std::vector<ReferenceFamilyMember>& members) {
  if (members.empty()) return Core::Status::Failure(Core::StatusCode::InvalidState,"empty reference family");
  std::map<std::tuple<std::string,std::uint64_t,double,double>,std::vector<ReferenceFamilyMember>> groups;
  const double target=members.front().integratedPeclet;
  for (const auto& m:members) {
    if (m.sourceIsSeparable||m.offsetM<=0||m.jointMeasure<=0||!std::isfinite(m.offsetM)||!std::isfinite(m.jointMeasure)||!m.stratum.geometryGeneration||
        !(target>0)||!std::isfinite(m.integratedPeclet)||std::abs(m.integratedPeclet-target)>1e-8*target)
      return Core::Status::Failure(Core::StatusCode::InvalidState,"invalid/non-joint source family");
    const auto& s=m.stratum;
    groups[std::make_tuple(s.speciesId,s.patchId,s.timeLowerS,s.timeUpperS)].push_back(m);
  }
  for (auto& group:groups) {
    auto& values=group.second;
    std::sort(values.begin(),values.end(),[](const auto& a,const auto& b){return a.stratum.momentumLowerSi<b.stratum.momentumLowerSi;});
    for (std::size_t i=1;i<values.size();i++) if (values[i-1].stratum.momentumUpperSi!=values[i].stratum.momentumLowerSi||
        values[i-1].stratum.geometryGeneration!=values[i].stratum.geometryGeneration)
      return Core::Status::Failure(Core::StatusCode::InvalidState,"family gap/overlap/generation mismatch");
  }
  return Core::Status::Success();
}
Core::Result<RenewalBalance> EvaluateRenewal(double number,double mass,Vec3 incoming,
    const std::vector<RenewalBranch>& branches,bool authority,double tol) {
  if (!(number>0)||!(mass>0)||!std::isfinite(number)||!std::isfinite(mass)||!finite(incoming)||branches.empty()||!authority||!(tol>0)||!std::isfinite(tol))
    return bad<RenewalBalance>("renewal lacks validated conditional work authority");
  RenewalBalance out; out.incomingNumber=number; out.incomingEnergyJ=number*energy(mass,incoming); out.incomingMomentum=number*incoming;
  if(!std::isfinite(out.incomingEnergyJ)||!finite(out.incomingMomentum)) return bad<RenewalBalance>("nonfinite physical renewal inventory");
  double probability=0;
  for (const auto& b:branches) {
    if (!(b.probability>=0)||!std::isfinite(b.probability)||!finite(b.outgoingMomentumKgMPerS)||
        !finite(b.declaredShockImpulseKgMPerS)||!std::isfinite(b.declaredShockWorkJ)||!std::isfinite(b.residenceTimeS)||b.residenceTimeS<0||
        (b.outcome!=RenewalOutcome::Absorbed&&b.outcome!=RenewalOutcome::Downstream&&b.outcome!=RenewalOutcome::ReReleased)) return bad<RenewalBalance>("invalid renewal branch");
    if (b.outcome==RenewalOutcome::ReReleased&&(b.residenceTimeS<=0||
        Norm(b.outgoingMomentumKgMPerS-incoming)==0)) return bad<RenewalBalance>("unchanged/zero-delay re-emission double counts first passage");
    // Subtract kinetic energies in a cancellation-safe form. Rest energy is
    // unchanged; using it as the work-error scale would accept large relative
    // errors for nonrelativistic particles.
    const double work=kinetic(mass,b.outgoingMomentumKgMPerS)-kinetic(mass,incoming);
    if(!std::isfinite(work)||!std::isfinite(energy(mass,b.outgoingMomentumKgMPerS))) return bad<RenewalBalance>("nonfinite branch energy");
    const Vec3 impulse=b.outgoingMomentumKgMPerS-incoming;
    if (std::abs(work-b.declaredShockWorkJ)>tol*std::max(kinetic(mass,incoming),std::max(kinetic(mass,b.outgoingMomentumKgMPerS),std::abs(work)))||
        Norm(impulse-b.declaredShockImpulseKgMPerS)>tol*std::max(Norm(incoming),Norm(impulse)))
      return bad<RenewalBalance>("renewal four-momentum/shock-work authority mismatch");
    probability+=b.probability; const double weight=number*b.probability;
    if (b.outcome==RenewalOutcome::Absorbed) out.absorbedNumber+=weight;
    if (b.outcome==RenewalOutcome::Downstream) out.downstreamNumber+=weight;
    if (b.outcome==RenewalOutcome::ReReleased) out.reReleasedNumber+=weight;
    out.outgoingEnergyJ+=weight*energy(mass,b.outgoingMomentumKgMPerS);
    out.outgoingMomentum=out.outgoingMomentum+weight*b.outgoingMomentumKgMPerS;
    out.shockWorkJ+=weight*work; out.shockImpulse=out.shockImpulse+weight*impulse;
  }
  if (std::abs(probability-1)>tol) return bad<RenewalBalance>("conditional renewal probability not normalized");
  return Core::Result<RenewalBalance>::Success(out);
}
Core::Status ImpulsiveLedger::Birth(const ImpulsiveBirth& b) {
  if (!b.immutableBirthId||births_.count(b.immutableBirthId)||!(b.representedNumber>0)||
      !(b.kineticEnergyJ>=0)||!std::isfinite(b.representedNumber)||!std::isfinite(b.kineticEnergyJ)||
      b.sourceOrigin!="ImpulsiveCoronalRelease"||b.frontRelation==ImpulsiveFrontRelation::Downstream||
      (b.frontRelation==ImpulsiveFrontRelation::Upstream&&!b.shockGeneration))
    return Core::Status::Failure(Core::StatusCode::InvalidState,"invalid impulsive birth/mapping/front relation");
  births_[b.immutableBirthId]=b; bornNumber_+=b.representedNumber; bornEnergy_+=b.representedNumber*b.kineticEnergyJ; return Core::Status::Success();
}
Core::Status ImpulsiveLedger::AbsorbFrontEncounter(std::uint64_t id,std::uint64_t generation) {
  if (!births_.count(id)||removed_.count(id)||!generation||generation<births_[id].shockGeneration)
    return Core::Status::Failure(Core::StatusCode::InvalidState,"duplicate/stale impulsive front encounter");
  removed_[id]=true; frontNumber_+=births_[id].representedNumber;
  frontEnergy_+=births_[id].representedNumber*births_[id].kineticEnergyJ;
  return Core::Status::Success();
}
} } }

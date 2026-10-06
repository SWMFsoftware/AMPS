#include "provider.h"
#include "diagnostics.h"

#include "sep_coronal_cme/constants.h"
#include "sep_coronal_cme/ellipsoid_geometry.h"
#include "sep_coronal_cme/research_extensions.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace SF=SEP::CoronaSwcme::ShockFront;
namespace CME=SEP::CoronalCME;

namespace {

struct Counts { int pass=0,fail=0,skip=0,error=0; } counts;

void Check(bool condition,const std::string& id,const std::string& message) {
  if(condition) { ++counts.pass;std::cout<<"PASS "<<id<<" "<<message<<'\n'; }
  else { ++counts.fail;std::cerr<<"FAIL "<<id<<" "<<message<<'\n'; }
}

bool Near(double value,double expected,double relative,double absolute=0) {
  return std::abs(value-expected)<=absolute+relative*std::max(std::abs(expected),1.0);
}

CME::Vec3 RotateZ(CME::Vec3 value,double angle) {
  const double c=std::cos(angle),s=std::sin(angle);
  return {c*value.x-s*value.y,s*value.x+c*value.y,value.z};
}

double RelativeVectorError(CME::Vec3 value,CME::Vec3 expected) {
  return CME::Norm(value-expected)/std::max(1e-300,CME::Norm(expected));
}

double ReferenceEnergyFlux(const CME::MhdPrimitiveState& state,
    CME::Vec3 normal,double shockSpeed,double gamma) {
  // Independent conservative-flux evaluation used by RSH17/RSH18.  This is
  // intentionally not obtained from the solver's published residual, so a
  // shared normalization or reconstruction defect cannot self-validate.
  const CME::Vec3 u=state.velocityMPerS-shockSpeed*normal;
  const double un=CME::Dot(u,normal);
  const double bn=CME::Dot(state.magneticFieldT,normal);
  const double b2=CME::Dot(state.magneticFieldT,state.magneticFieldT);
  return un*(0.5*state.massDensityKgM3*CME::Dot(u,u)+
      gamma*state.pressurePa/(gamma-1)+
      b2/CME::Constants::kVacuumPermeabilityHPerM)-
      bn*CME::Dot(u,state.magneticFieldT)/
      CME::Constants::kVacuumPermeabilityHPerM;
}

double NormalizedDivergence(const SEP::CoronaSwcme::AmbientState& state,
    double lengthM) {
  const double divergence=state.gradientB[0]+state.gradientB[4]+state.gradientB[8];
  return std::abs(divergence)*lengthM/
      std::max(1e-300,CME::Norm(state.primitive.magneticFieldT));
}

std::string Bytes(const std::filesystem::path& path) {
  std::ifstream input(path,std::ios::binary);
  if(!input)throw std::runtime_error("cannot read "+path.string());
  std::ostringstream output;output<<input.rdbuf();return output.str();
}

std::filesystem::path Root() {
  for(const auto& candidate:{std::filesystem::path("."),
      std::filesystem::path("../../..")})
    if(std::filesystem::exists(candidate/"CODEX_REDUCED_SHOCK_TASK.txt"))
      return std::filesystem::canonical(candidate);
  throw std::runtime_error("cannot locate AMPS root");
}

std::shared_ptr<const SF::Configuration> Load(const std::filesystem::path& event) {
  const auto directory=event.parent_path();
  const auto resolved=SF::ResolveConfiguration(Bytes(event),[&](const std::string& name) {
    try { return SEP::Core::Result<std::string>::Success(Bytes(directory/name)); }
    catch(const std::exception& e) { return SEP::Core::Result<std::string>::Failure(
        SEP::Core::StatusCode::OutOfDomain,e.what()); }
  });
  if(!resolved.ok())throw std::runtime_error(resolved.status.message);
  return resolved.value;
}

} // namespace

int main() {
  try {
    const auto root=Root();
    const auto examples=root/"srcSEP3D/examples/shock-front";
    const auto smoke=Load(examples/"handoff_smoke.event");
    const auto longEvent=Load(examples/"corona_to_1au.event");
    const auto positiveEvent=Load(examples/"positive_1au.event");
    Check(smoke->physicsFingerprint.size()==64&&longEvent->physicsFingerprint.size()==64&&
        smoke->physicsFingerprint!=longEvent->physicsFingerprint&&
        smoke->initialApexRadiusM==13844430000.0&&
        smoke->initialApexSpeedMPerS==1000000.0&&
        smoke->coordinateFrame=="HCI"&&smoke->particleMode=="disabled","RSH01",
        "strict SI/HCI resolved events bind exact values, transitive asset and profile bytes");
    std::string bad=Bytes(examples/"handoff_smoke.event")+"unknown.key=1\n";
    auto rejected=SF::ResolveConfiguration(bad,[&](const std::string& name){
      return SEP::Core::Result<std::string>::Success(Bytes(examples/name));});
    Check(!rejected.ok(),"RSH01","unknown event key is rejected");
    bad=Bytes(examples/"handoff_smoke.event")+
        "history.initial_apex_speed_m_s=1000000\n";
    rejected=SF::ResolveConfiguration(bad,[&](const std::string& name){
      return SEP::Core::Result<std::string>::Success(Bytes(examples/name));});
    Check(!rejected.ok(),"RSH01","duplicate event key is rejected before construction");
    bad=Bytes(examples/"handoff_smoke.event");
    const std::string validEpoch="run.reference_epoch=2026-10-04T00:00:00Z";
    auto epochAt=bad.find(validEpoch);
    bad.replace(epochAt,validEpoch.size(),
        "run.reference_epoch=2026-02-30T00:00:00Z");
    rejected=SF::ResolveConfiguration(bad,[&](const std::string& name){
      return SEP::Core::Result<std::string>::Success(Bytes(examples/name));});
    Check(!rejected.ok(),"RSH02","invalid UTC calendar label is rejected");
    bad=Bytes(examples/"handoff_smoke.event");
    const auto checksum=bad.find("assets.harmonics_sha256=");
    bad[checksum+sizeof("assets.harmonics_sha256=")-1]='0';
    rejected=SF::ResolveConfiguration(bad,[&](const std::string& name){
      return SEP::Core::Result<std::string>::Success(Bytes(examples/name));});
    Check(!rejected.ok()&&rejected.status.code==SEP::Core::StatusCode::DataIntegrityFailure,
        "RSH01","changed harmonic bytes/identity are rejected");

    const auto made=SF::Provider::Create(smoke);
    Check(made.ok(),"RSH00","shared provider constructs without contact/sheath/particle state");
    if(!made.ok())throw std::runtime_error(made.status.message);
    const auto provider=made.value;
    Check(provider->RequireCapability(SF::VolumeCapability::AmbientReference).ok()&&
        !provider->RequireCapability(SF::VolumeCapability::PhysicalDownstreamVolume).ok(),
        "RSH01","ambient is supported and downstream volume is typed unsupported");

    const double lambda=smoke->halfWidthRad;
    const auto apex=provider->EvaluateRay(smoke->direction,0,1);
    const CME::Vec3 edge={std::cos(lambda),std::sin(lambda),0};
    const auto flank=provider->EvaluateRay(edge,0,2);
    Check(apex.ok()&&Near(CME::Norm(apex.value.positionM),smoke->initialApexRadiusM,2e-15)&&
        Near(apex.value.normalSpeedMPerS,1e6,2e-15),"RSH07",
        "finite SSE apex radius and normal speed are exact");
    Check(flank.ok()&&std::abs(flank.value.normalSpeedMPerS)<1e-8&&
        flank.value.supportEdge,"RSH07",
        "fixed-width tangent flank is geometric support with zero normal speed");
    const CME::Vec3 miss={std::cos(lambda+0.01),std::sin(lambda+0.01),0};
    Check(!provider->EvaluateRay(miss,0).ok(),"RSH07",
        "ray outside finite SSE support is absent");

    // RSH03: independent translating-plane level-set reference.  Its normal
    // speed is n.V, not the magnitude of a chosen apex/radial velocity.
    const CME::Vec3 planeNormal=CME::Unit({1,2,-1});
    const CME::Vec3 planeVelocity={8,-3,5};
    const CME::Vec3 planePoint0={4,-2,1};const double planeTime=3.0;
    const CME::Vec3 planePoint=planePoint0+planeTime*planeVelocity;
    const double planeG=CME::Dot(planeNormal,planePoint-
        (planePoint0+planeTime*planeVelocity));
    const double planeSpeed=CME::Dot(planeNormal,planeVelocity);
    Check(std::abs(planeG)<1e-14&&
        std::abs(planeSpeed-CME::Norm(planeVelocity))>1&&
        std::abs(planeSpeed-planeVelocity.x)>1,"RSH03",
        "translating-plane normal speed is analytic and rejects apex/radial substitutions");

    // RSH04: the maintained dynamic-ellipsoid kernel reduces exactly to a
    // Sun-centred expanding sphere when all axes/rates agree.  A conventional
    // theta-midpoint area integral is independently second-order convergent.
    CME::Research::DynamicGeometryState sphereState;
    sphereState.centerM={0,0,0};sphereState.centerVelocityMPerS={0,0,0};
    sphereState.rotationAxis={0,0,1};sphereState.axesM={3,3,3};
    sphereState.axisRatesMPerS={2,2,2};sphereState.frameId="HCI";
    const auto sphere=CME::Research::DynamicEllipsoid::Build(sphereState);
    const CME::Vec3 spherePoint=3*CME::Unit({1,2,3});
    const auto sphereEvaluation=sphere.ok()?sphere.value.Evaluate(spherePoint):
        SEP::Core::Result<CME::SurfaceEvaluation>::Failure(
          SEP::Core::StatusCode::NumericalFailure,"sphere unavailable");
    double sphereErrors[3]={};int sphereN[3]={8,16,32};
    for(int level=0;level<3;++level) {
      double area=0;const int n=sphereN[level];
      for(int i=0;i<n;++i)for(int j=0;j<2*n;++j) {
        const double th=(i+0.5)*CME::Constants::kPi/n;
        area+=9*std::sin(th)*(CME::Constants::kPi/n)*
          (2*CME::Constants::kPi/(2*n));
      }
      sphereErrors[level]=std::abs(area-4*CME::Constants::kPi*9);
    }
    Check(sphere.ok()&&sphereEvaluation.ok()&&
        RelativeVectorError(sphereEvaluation.value.outwardNormal,CME::Unit(spherePoint))<1e-14&&
        Near(sphereEvaluation.value.normalSpeedMPerS,2,2e-14)&&
        sphereErrors[0]/sphereErrors[1]>3.9&&sphereErrors[1]/sphereErrors[2]>3.9,
        "RSH04","Sun-centred sphere normal/speed are exact and area quadrature converges at second order");

    // RSH05: translating/expanding triaxial level-set derivatives are checked
    // by rebuilding independent +/- histories and differencing g at a fixed
    // Cartesian point.  |g|/|grad g| is also tested as a local distance, not
    // mistaken for a dimensionless position error.
    const auto ellipsoidBasis=CME::BuildRadialPrincipalBasis(0.2,0.4,0.3);
    auto ellipsoidAt=[&](double time) {
      CME::EllipsoidKinematics k;
      k.centerDistanceM={10+2*time,2,0};
      k.radialSemiAxisM={3+0.4*time,0.4,0};
      k.firstLateralSemiAxisM={2+0.2*time,0.2,0};
      k.secondLateralSemiAxisM={1+0.1*time,0.1,0};
      return CME::FixedOrientationEllipsoid::FromCenter(ellipsoidBasis.value,k);
    };
    const double ellipsoidTime=0.5,ellipsoidH=1e-5;
    const auto ellipsoid=ellipsoidAt(ellipsoidTime);
    const auto ellipsoidPlus=ellipsoidAt(ellipsoidTime+ellipsoidH);
    const auto ellipsoidMinus=ellipsoidAt(ellipsoidTime-ellipsoidH);
    const CME::Vec3 ellipsoidPoint=ellipsoid.value.Point(1.1,0.7);
    const auto ellipsoidEval=ellipsoid.value.Evaluate(ellipsoidPoint);
    const double gt=(ellipsoidPlus.value.Evaluate(ellipsoidPoint).value.implicitValue-
        ellipsoidMinus.value.Evaluate(ellipsoidPoint).value.implicitValue)/(2*ellipsoidH);
    const double spaceH=1e-5;
    const double gn=(ellipsoid.value.Evaluate(ellipsoidPoint+
        spaceH*ellipsoidEval.value.outwardNormal).value.implicitValue-
        ellipsoid.value.Evaluate(ellipsoidPoint-
        spaceH*ellipsoidEval.value.outwardNormal).value.implicitValue)/(2*spaceH);
    const double offset=1e-5;
    const auto displaced=ellipsoid.value.Evaluate(ellipsoidPoint+
        offset*ellipsoidEval.value.outwardNormal);
    const double displacedGradient=(ellipsoid.value.Evaluate(ellipsoidPoint+
        (offset+spaceH)*ellipsoidEval.value.outwardNormal).value.implicitValue-
        ellipsoid.value.Evaluate(ellipsoidPoint+
        (offset-spaceH)*ellipsoidEval.value.outwardNormal).value.implicitValue)/(2*spaceH);
    Check(ellipsoidBasis.ok()&&ellipsoid.ok()&&ellipsoidEval.ok()&&displaced.ok()&&
        Near(ellipsoidEval.value.normalSpeedMPerS,-gt/gn,2e-9)&&
        Near(std::abs(displaced.value.implicitValue/displacedGradient),offset,2e-5),
        "RSH05","maintained triaxial level set matches independent time differentiation and distance residual");

    // RSH06: general rotating/expanding triaxial motion.  Independent point
    // differentiation recovers the full surface velocity; deleting omega x r
    // is a controlled negative that changes its normal component.
    auto rotatingState=[](double time) {
      CME::Research::DynamicGeometryState s;s.centerM={10+2*time,3*time,0};
      s.centerVelocityMPerS={2,3,0};s.rotationAxis={0,0,1};
      s.angleRad=0.4*time;s.angularRateRadPerS=0.4;
      s.axesM={3+time,2+0.2*time,1+0.1*time};
      s.axisRatesMPerS={1,0.2,0.1};s.epochS=time;s.frameId="HCI";return s;
    };
    const double rotateTime=0.5,rotateH=1e-5;
    const auto rotating=CME::Research::DynamicEllipsoid::Build(rotatingState(rotateTime));
    const auto rotatingPlus=CME::Research::DynamicEllipsoid::Build(rotatingState(rotateTime+rotateH));
    const auto rotatingMinus=CME::Research::DynamicEllipsoid::Build(rotatingState(rotateTime-rotateH));
    const CME::Vec3 rotatingPoint=rotating.value.Point(1.1,0.7);
    const auto rotatingEval=rotating.value.Evaluate(rotatingPoint);
    const CME::Vec3 rotatingFd=(rotatingPlus.value.Point(1.1,0.7)-
        rotatingMinus.value.Point(1.1,0.7))/(2*rotateH);
    auto omitted=rotatingState(rotateTime);omitted.angularRateRadPerS=0;
    const auto withoutRotation=CME::Research::DynamicEllipsoid::Build(omitted);
    const auto omittedEval=withoutRotation.value.Evaluate(rotatingPoint);
    Check(rotating.ok()&&rotatingEval.ok()&&CME::Norm(rotatingFd-
        rotating.value.SurfaceVelocity(1.1,0.7))<2e-8&&
        std::abs(rotatingEval.value.normalSpeedMPerS-
          omittedEval.value.normalSpeedMPerS)>1e-2,"RSH06",
        "rotating triaxial normal speed includes the independently detected attitude derivative");

    // Complete SSE geometry checks beyond apex/tangent: leading versus rear
    // root, radial/normal speed distinction, regular width limits, and the
    // uniform Ma=2, lambda=45-degree accepted half-angle reference.
    SF::SseKinematicState sse{smoke->initialApexRadiusM,
        smoke->initialApexSpeedMPerS,smoke->direction,{},lambda,0};
    const double psi=0.5*lambda;
    const CME::Vec3 midRay={std::cos(psi),std::sin(psi),0};
    const auto midSse=SF::EvaluateSseRay(sse,midRay,7);
    const double ss=std::sin(lambda),cc=smoke->initialApexRadiusM/(1+ss);
    const double aa=cc*ss;
    const double leading=cc*std::cos(psi)+std::sqrt(aa*aa-
        cc*cc*std::sin(psi)*std::sin(psi));
    const double rear=cc*std::cos(psi)-std::sqrt(aa*aa-
        cc*cc*std::sin(psi)*std::sin(psi));
    auto narrow=sse;narrow.halfWidthRad=1e-5;
    auto wide=sse;wide.halfWidthRad=0.5*CME::Constants::kPi-1e-5;
    const auto narrowApex=SF::EvaluateSseRay(narrow,smoke->direction);
    const auto wideApex=SF::EvaluateSseRay(wide,smoke->direction);
    double lowAngle=0,highAngle=lambda;
    for(int iteration=0;iteration<80;++iteration) {
      const double test=0.5*(lowAngle+highAngle);
      const auto q=SF::EvaluateSseRay(sse,{std::cos(test),std::sin(test),0});
      if(q.value.normalSpeedMPerS>0.5*smoke->initialApexSpeedMPerS)
        lowAngle=test;else highAngle=test;
    }
    const double halfAngleDeg=0.5*(lowAngle+highAngle)*180/CME::Constants::kPi;
    const double radialSpeed=midSse.ok()?smoke->initialApexSpeedMPerS*
        CME::Norm(midSse.value.positionM)/smoke->initialApexRadiusM:0;
    Check(midSse.ok()&&Near(CME::Norm(midSse.value.positionM),leading,2e-14)&&
        leading>rear&&std::abs(radialSpeed-midSse.value.normalSpeedMPerS)>1e3&&
        narrowApex.ok()&&wideApex.ok()&&
        Near(CME::Norm(narrowApex.value.positionM),smoke->initialApexRadiusM,2e-14)&&
        Near(CME::Norm(wideApex.value.positionM),smoke->initialApexRadiusM,2e-14)&&
        std::abs(halfAngleDeg-32.4)<0.1,"RSH07",
        "production SSE selects the leading root, distinguishes radial/normal speed, retains regular width limits and matches the uniform shock half-angle");

    const double crossingExpected=69.57;
    const auto handoff=provider->HandoffTimeS();
    const auto before=provider->Trajectory(handoff.value-1e-4);
    const auto at=provider->Trajectory(handoff.value);
    const auto after=provider->Trajectory(handoff.value+1e-4);
    Check(handoff.ok()&&Near(handoff.value,crossingExpected,0,2e-9)&&before.ok()&&
        at.ok()&&after.ok()&&Near(at.value.apexRadiusM,smoke->handoffApexRadiusM,0,1e-6)&&
        Near(before.value.apexSpeedMPerS,after.value.apexSpeedMPerS,2e-9),"RSH20",
        "19.9-to-20 Rs crossing is bracketed at 69.57 s without radius/speed reset");
    Check(Near(after.value.apexAccelerationMPerS2,-7.2,2e-5),"RSH20",
        "outer one-sided acceleration records the declared C1, non-C2 handoff");
    const auto smokeHandoff=provider->Handoff();
    Check(smokeHandoff.ok()&&Near(smokeHandoff.value.timeS,69.57,0,2e-9)&&
        smokeHandoff.value.c1Matched&&!smokeHandoff.value.c2Matched&&
        !smokeHandoff.value.authorityReset&&
        Near(smokeHandoff.value.outerAccelerationMPerS2,-7.2,2e-14)&&
        smokeHandoff.value.maximumSurfacePositionMismatchM<=1e-6&&
        smokeHandoff.value.maximumNormalSpeedMismatchMPerS<=1e-9,"RSH20",
        "immutable handoff receipt records exact root, one authority, C1 match and the declared acceleration jump");

    const auto longMade=SF::Provider::Create(longEvent);
    Check(longMade.ok(),"RSH08","quintic-pulse long provider constructs");
    if(!longMade.ok())throw std::runtime_error(longMade.status.message);
    const auto longProvider=longMade.value;
    const auto longHandoff=longProvider->Handoff();
    const double handoffDelta=1e-3;
    const auto handoffBefore=longProvider->Prepare(
        longHandoff.value.timeS-handoffDelta,20);
    const auto handoffCrossing=longProvider->Prepare(longHandoff.value.timeS,21);
    const auto handoffAfter=longProvider->Prepare(
        longHandoff.value.timeS+handoffDelta,22);
    bool jumpContinuous=handoffBefore.ok()&&handoffCrossing.ok()&&handoffAfter.ok();
    double maximumCompressionChange=0;int matchedJumps=0;
    if(jumpContinuous)for(std::size_t i=0;i<handoffBefore.value->records.size();++i) {
      const auto& left=handoffBefore.value->records[i];
      const auto& right=handoffAfter.value->records[i];
      if(left.status!=SF::FrontStatus::SolvedFastShock||
          right.status!=SF::FrontStatus::SolvedFastShock)continue;
      ++matchedJumps;maximumCompressionChange=std::max(maximumCompressionChange,
          std::abs(left.jump.compressionRatio-right.jump.compressionRatio)/
          std::max(1.0,left.jump.compressionRatio));
    }
    Check(longHandoff.ok()&&handoffBefore.ok()&&handoffCrossing.ok()&&
        handoffAfter.ok()&&longHandoff.value.c1Matched&&
        longHandoff.value.canonicalWindValid&&
        handoffCrossing.value->handoff.timeS==longHandoff.value.timeS&&
        handoffBefore.value->eventIdentity==handoffCrossing.value->eventIdentity&&
        handoffCrossing.value->eventIdentity==handoffAfter.value->eventIdentity&&
        handoffBefore.value->trajectory.phase==SF::Phase::CoronalHistory&&
        handoffAfter.value->trajectory.phase==SF::Phase::SwcmeOuter&&
        matchedJumps>0&&maximumCompressionChange<2e-5,"RSH20",
        "pre/crossing/post full-surface epochs preserve identity and continuous local jumps across the exact matched state");

    // RSH02: the map is anchored in the rotating solar frame while every
    // returned vector is HCI.  Co-rotating both the query and epoch must rotate
    // B and U rigidly; it must not change scalar plasma state or time order.
    const double rotationTime=1234.0;
    const double theta=0.37,phi=0.28,radius=10*longEvent->ambient.support.solarRadiusM;
    const CME::Vec3 x0={radius*std::sin(theta)*std::cos(phi),
        radius*std::sin(theta)*std::sin(phi),radius*std::cos(theta)};
    const double angle=longEvent->ambient.ambient.rotationRateRadPerS*rotationTime;
    const auto frame0=longProvider->QueryAmbient(x0,0,10);
    const auto frame1=longProvider->QueryAmbient(RotateZ(x0,angle),rotationTime,11);
    Check(frame0.ok()&&frame1.ok()&&
        RelativeVectorError(frame1.value.primitive.magneticFieldT,
          RotateZ(frame0.value.primitive.magneticFieldT,angle))<2e-10&&
        RelativeVectorError(frame1.value.primitive.velocityMPerS,
          RotateZ(frame0.value.primitive.velocityMPerS,angle))<2e-10&&
        Near(frame1.value.primitive.plasma.massDensityKgM3,
          frame0.value.primitive.plasma.massDensityKgM3,2e-12)&&
        frame1.value.epochS>frame0.value.epochS,"RSH02",
        "HCI rotation, elapsed-SI time ordering and scalar invariance agree with an independent rigid rotation");

    // RSH10: independent l=1,m=0 PFSS source-surface limit.  For the real
    // orthonormal dipole, Br(Rss)=g10*Y10*3/[1+2(Rss/Rsun)^3].
    const auto& ambientDefinition=longEvent->ambient;
    const double rs=ambientDefinition.support.solarRadiusM;
    const double rss=ambientDefinition.ambient.sourceSurfaceRadiusM;
    const double g10=0.0002046653415892977;
    const double y10=std::sqrt(3.0/(4.0*CME::Constants::kPi));
    const double expectedSourcePole=g10*y10*3.0/
        (1.0+2.0*std::pow(rss/rs,3));
    const auto north=longProvider->QueryAmbient({0,0,rss},0,12);
    const auto south=longProvider->QueryAmbient({0,0,-rss},0,13);
    Check(north.ok()&&south.ok()&&Near(north.value.primitive.magneticFieldT.z,
        expectedSourcePole,3e-12,1e-18)&&Near(south.value.primitive.magneticFieldT.z,
        expectedSourcePole,3e-12,1e-18)&&
        std::hypot(north.value.primitive.magneticFieldT.x,
          north.value.primitive.magneticFieldT.y)<1e-18,"RSH10",
        "production PFSS matches the independent orthonormal dipole source-surface and pole limits");

    // RSH11: midpoint solid-angle quadrature of the signed source-surface
    // flux.  An even theta count avoids sampling the true dipole null itself;
    // the zero signed result is normalized by nonzero unsigned flux.
    double signedFlux=0,unsignedFlux=0;bool fluxQueries=true;
    const int nt=32,np=64;
    for(int i=0;i<nt;++i)for(int j=0;j<np;++j) {
      const double th=(i+0.5)*CME::Constants::kPi/nt;
      const double ph=(j+0.5)*2*CME::Constants::kPi/np;
      const CME::Vec3 n={std::sin(th)*std::cos(ph),std::sin(th)*std::sin(ph),std::cos(th)};
      const auto sample=longProvider->QueryAmbient(rss*n,0,14);
      if(!sample.ok()){fluxQueries=false;continue;}
      const double element=rss*rss*std::sin(th)*(CME::Constants::kPi/nt)*
          (2*CME::Constants::kPi/np);
      const double flux=CME::Dot(sample.value.primitive.magneticFieldT,n)*element;
      signedFlux+=flux;unsignedFlux+=std::abs(flux);
    }
    Check(fluxQueries&&unsignedFlux>0&&std::abs(signedFlux)/unsignedFlux<2e-15&&
        longEvent->harmonicAssetSha256==
          "b5f864f33b85034125a8c8ad162c3880d82985c1a97aafa8c7745a9737c751f9",
        "RSH11","signed source-surface flux closes while raw asset identity and nonzero unsigned flux remain explicit");

    // RSH12: exercise the production interpolated wind and Cartesian
    // derivative path on both the axisymmetric fixture and a nonaxisymmetric
    // harmonic definition.  A manufactured position-dependent scalar applied
    // after evaluation is a negative control: it introduces div B and is
    // correctly distinguishable from a globally calibrated field.
    const auto exterior=longProvider->QueryAmbient(x0,0,15);
    auto nonaxisDefinition=longEvent->ambient;
    nonaxisDefinition.ambient.harmonics="2,1,3e-5,-2e-5";
    nonaxisDefinition.physicsFingerprint="rsh12-nonaxisymmetric-reference";
    const auto nonaxisModel=SEP::CoronaSwcme::AmbientModel::Create(nonaxisDefinition);
    const auto nonaxis=nonaxisModel.ok()?nonaxisModel.value->EvaluateWithDerivatives(
        x0,0,1):SEP::Core::Result<SEP::CoronaSwcme::AmbientState>::Failure(
        SEP::Core::StatusCode::NumericalFailure,"nonaxis model unavailable");
    const double axisDiv=exterior.ok()?NormalizedDivergence(exterior.value,radius):1;
    const double nonaxisDiv=nonaxis.ok()?NormalizedDivergence(nonaxis.value,radius):1;
    // div[(1+alpha*x/r0)B] = scale*div(B)+alpha*Bx/r0.  This deliberately
    // illegal local rescaling must be visibly worse than the production field.
    const double alpha=0.2;
    const double badDiv=nonaxis.ok()?std::abs(alpha*nonaxis.value.primitive.magneticFieldT.x/radius)*
        radius/std::max(1e-300,CME::Norm(nonaxis.value.primitive.magneticFieldT)):0;
    Check(exterior.ok()&&nonaxis.ok()&&axisDiv<2e-7&&nonaxisDiv<3e-6&&
        badDiv>100*std::max(nonaxisDiv,1e-12),"RSH12",
        "interpolated axisymmetric/nonaxisymmetric exterior fields are divergence-consistent and a local-amplitude negative control is detected");

    // RSH13: the exterior is seeded from the one-sided PFSS radial state.
    // B and normal flux are continuous; their derivatives are allowed to jump
    // because PFSS and Parker are different analytical regions.
    const CME::Vec3 matchDirection={std::sin(0.31),0,std::cos(0.31)};
    const double matchEpsilon=1e-7;
    const auto matchIn=longProvider->QueryAmbient(rss*(1-matchEpsilon)*matchDirection,0,16);
    const auto matchOut=longProvider->QueryAmbient(rss*(1+matchEpsilon)*matchDirection,0,17);
    const double matchScale=matchIn.ok()&&matchOut.ok()?std::max(
        CME::Norm(matchIn.value.primitive.magneticFieldT),
        CME::Norm(matchOut.value.primitive.magneticFieldT)):1;
    Check(matchIn.ok()&&matchOut.ok()&&
        CME::Norm(matchIn.value.primitive.magneticFieldT-
          matchOut.value.primitive.magneticFieldT)/matchScale<8e-7&&
        matchIn.value.primitive.magneticSector==matchOut.value.primitive.magneticSector&&
        std::abs(CME::Dot(matchIn.value.primitive.magneticFieldT,matchDirection)-
          CME::Dot(matchOut.value.primitive.magneticFieldT,matchDirection))/matchScale<5e-7,
        "RSH13","PFSS/Parker matching preserves vector field, normal flux and polarity without requiring derivative continuity");

    // RSH14: pure hydrogen with Tp=Te has rho=ne*mp and total thermal
    // pressure ne*kB*(Te+Tp).  Test the production ambient values rather than
    // invoking the EOS twice, then independently reconstruct all wave limits.
    const auto eos=longProvider->QueryAmbient({0,0,2*rs},0,18);
    bool eosGood=eos.ok();
    if(eos.ok()) {
      const auto& p=eos.value.primitive;
      const double ne=p.plasma.electronNumberDensityM3;
      const double rho=ne*CME::Constants::kProtonMassKg;
      const double pressure=ne*CME::Constants::kBoltzmannJPerK*
          (p.electronTemperatureK+p.protonTemperatureK);
      const double cs=std::sqrt(longEvent->ambient.composition.gammaAdiabatic*pressure/rho);
      const double va=CME::Norm(p.magneticFieldT)/std::sqrt(
          CME::Constants::kVacuumPermeabilityHPerM*rho);
      CME::MhdPrimitiveState primitive{rho,pressure,p.velocityMPerS,p.magneticFieldT};
      const auto parallel=CME::EvaluateMhdCharacteristics(primitive,{0,0,1},
          longEvent->ambient.composition.gammaAdiabatic);
      const auto perpendicular=CME::EvaluateMhdCharacteristics(primitive,{1,0,0},
          longEvent->ambient.composition.gammaAdiabatic);
      eosGood=Near(p.plasma.massDensityKgM3,rho,3e-14)&&
          Near(p.plasma.pressurePa,pressure,3e-14)&&Near(p.plasma.soundSpeedMPerS,cs,3e-14)&&
          Near(p.plasma.alfvenSpeedMPerS,va,3e-14)&&parallel.ok()&&perpendicular.ok()&&
          Near(parallel.value.fastSpeedMPerS,std::max(cs,va),2e-14)&&
          Near(perpendicular.value.fastSpeedMPerS,std::hypot(cs,va),2e-14);
    }
    Check(eosGood,"RSH14","production ambient retains electron+ion pressure, hydrogen mass normalization and independent parallel/perpendicular wave limits");

    // RSH31/RSH37: vary one physical authority at a time.  Density changes
    // Alfven speed as rho^-1/2, global harmonic scaling changes B and vA
    // linearly without changing topology, and temperature changes cs.  The
    // last comparison constructs the forbidden independent exterior factor as
    // a discontinuity negative control rather than offering it as an option.
    auto densityDefinition=longEvent->ambient;
    densityDefinition.ambient.electronDensityAtReferenceM3*=2;
    densityDefinition.physicsFingerprint="rsh31-density";
    auto fieldDefinition=longEvent->ambient;
    fieldDefinition.ambient.harmonics="1,0,0.0004093306831785954,0";
    fieldDefinition.physicsFingerprint="rsh37-global-field";
    auto temperatureDefinition=longEvent->ambient;
    temperatureDefinition.composition.electronTemperatureK*=0.5;
    temperatureDefinition.composition.protonTemperatureK*=0.5;
    temperatureDefinition.physicsFingerprint="rsh31-temperature";
    const auto densityModel=SEP::CoronaSwcme::AmbientModel::Create(densityDefinition);
    const auto fieldModel=SEP::CoronaSwcme::AmbientModel::Create(fieldDefinition);
    const auto temperatureModel=SEP::CoronaSwcme::AmbientModel::Create(temperatureDefinition);
    const CME::Vec3 polarPoint={0,0,longEvent->ambient.ambient.referenceRadiusM};
    const auto nominal=longProvider->QueryAmbient(polarPoint,0,19);
    const auto dense=densityModel.ok()?densityModel.value->Evaluate(polarPoint,0):
        SEP::Core::Result<SEP::CoronaSwcme::AmbientPrimitive>::Failure(
          SEP::Core::StatusCode::NumericalFailure,"density model unavailable");
    const auto strong=fieldModel.ok()?fieldModel.value->Evaluate(polarPoint,0):
        SEP::Core::Result<SEP::CoronaSwcme::AmbientPrimitive>::Failure(
          SEP::Core::StatusCode::NumericalFailure,"field model unavailable");
    const auto cool=temperatureModel.ok()?temperatureModel.value->Evaluate(polarPoint,0):
        SEP::Core::Result<SEP::CoronaSwcme::AmbientPrimitive>::Failure(
          SEP::Core::StatusCode::NumericalFailure,"temperature model unavailable");
    const bool sensitivity=nominal.ok()&&dense.ok()&&strong.ok()&&cool.ok()&&
        Near(dense.value.plasma.massDensityKgM3/
          nominal.value.primitive.plasma.massDensityKgM3,2,2e-12)&&
        Near(dense.value.plasma.alfvenSpeedMPerS/
          nominal.value.primitive.plasma.alfvenSpeedMPerS,1/std::sqrt(2.0),2e-12)&&
        Near(CME::Norm(strong.value.magneticFieldT)/
          CME::Norm(nominal.value.primitive.magneticFieldT),2,3e-12)&&
        Near(strong.value.plasma.alfvenSpeedMPerS/
          nominal.value.primitive.plasma.alfvenSpeedMPerS,2,3e-12)&&
        Near(cool.value.plasma.soundSpeedMPerS/
          nominal.value.primitive.plasma.soundSpeedMPerS,1/std::sqrt(2.0),2e-12);
    Check(sensitivity,"RSH31","density, field and temperature sensitivity preserves EOS scalings without speed or Mach floors");
    const double forbiddenInterfaceMismatch=strong.ok()&&nominal.ok()?
        CME::Norm(1.1*strong.value.magneticFieldT-strong.value.magneticFieldT)/
          CME::Norm(strong.value.magneticFieldT):0;
    Check(sensitivity&&forbiddenInterfaceMismatch>0.099&&
        nominal.value.primitive.magneticSector==strong.value.magneticSector,
        "RSH37","one global harmonic scale preserves topology while an independent exterior radial scale fails flux continuity");
    const auto pulse=longProvider->Trajectory(600);
    Check(pulse.ok()&&Near(pulse.value.apexRadiusM,1.130055e9,0,1e-3)&&
        Near(pulse.value.apexSpeedMPerS,1e6,0,1e-8)&&
        std::abs(pulse.value.apexAccelerationMPerS2)<1e-12,"RSH08",
        "quintic history reaches the analytical pulse endpoint");
    const double h=0.01;const auto l=longProvider->Trajectory(300-h);
    const auto c=longProvider->Trajectory(300);const auto r=longProvider->Trajectory(300+h);
    Check(l.ok()&&c.ok()&&r.ok()&&Near((r.value.apexRadiusM-l.value.apexRadiusM)/(2*h),
        c.value.apexSpeedMPerS,2e-9),"RSH08",
        "independent centered differentiation recovers pulse speed");
    bool historyDerivatives=true;
    for(double time:{100.0,300.0,500.0}) {
      const double dh=1e-2;const auto hm=longProvider->Trajectory(time-dh);
      const auto hc=longProvider->Trajectory(time);const auto hp=longProvider->Trajectory(time+dh);
      historyDerivatives=historyDerivatives&&hm.ok()&&hc.ok()&&hp.ok()&&
          hc.value.apexRadiusM>0&&hc.value.apexSpeedMPerS>0&&
          Near((hp.value.apexRadiusM-hm.value.apexRadiusM)/(2*dh),
            hc.value.apexSpeedMPerS,3e-9)&&
          Near((hp.value.apexSpeedMPerS-hm.value.apexSpeedMPerS)/(2*dh),
            hc.value.apexAccelerationMPerS2,3e-8);
    }
    Check(historyDerivatives,"RSH08",
        "quintic launch remains positive between knots and independent derivatives recover speed and acceleration");

    // RSH09: all signs and smooth limits of the exact drag production
    // function, plus fourth-order convergence of an independent RK4 solve.
    const auto dragAbove=SF::EvaluateQuadraticDrag(2000,0,10,1000,400,2e-4);
    const auto dragBelow=SF::EvaluateQuadraticDrag(2000,0,10,200,400,2e-4);
    const auto dragEqual=SF::EvaluateQuadraticDrag(2000,0,10,400,400,2e-4);
    const auto dragZero=SF::EvaluateQuadraticDrag(2000,0,10,1000,400,0);
    auto integrate=[&](double step) {
      double position=10,speed=1000;
      auto acceleration=[](double u){return -2e-4*(u-400)*std::abs(u-400);};
      for(double time=0;time<2000-0.5*step;time+=step) {
        const double k1r=speed,k1v=acceleration(speed);
        const double k2r=speed+0.5*step*k1v,k2v=acceleration(speed+0.5*step*k1v);
        const double k3r=speed+0.5*step*k2v,k3v=acceleration(speed+0.5*step*k2v);
        const double k4r=speed+step*k3v,k4v=acceleration(speed+step*k3v);
        position+=step*(k1r+2*k2r+2*k3r+k4r)/6;
        speed+=step*(k1v+2*k2v+2*k3v+k4v)/6;
      }
      return std::pair<double,double>{position,speed};
    };
    const auto rk0=integrate(2),rk1=integrate(1),rk2=integrate(0.5);
    const double e0=std::abs(rk0.first-dragAbove.value.apexRadiusM);
    const double e1=std::abs(rk1.first-dragAbove.value.apexRadiusM);
    const double e2=std::abs(rk2.first-dragAbove.value.apexRadiusM);
    Check(dragAbove.ok()&&dragBelow.ok()&&dragEqual.ok()&&dragZero.ok()&&
        dragAbove.value.apexSpeedMPerS>400&&dragAbove.value.apexSpeedMPerS<1000&&
        dragBelow.value.apexSpeedMPerS<400&&dragBelow.value.apexSpeedMPerS>200&&
        dragEqual.value.apexSpeedMPerS==400&&dragZero.value.apexSpeedMPerS==1000&&
        e0/e1>12&&e1/e2>12,"RSH09",
        "quadratic drag handles above/below/equal/zero limits and matches an independently converged RK4 trajectory");

    // RSH32: the production SSE function carries width and direction-rate
    // terms.  Differentiate the same physical ray intersection while R, lambda
    // and d evolve, then remove each term separately as causal controls.
    const double dynamicOmega=2e-4,dynamicLambdaRate=1e-4;
    const double dynamicLambda=0.6,dynamicPsi=0.3;
    const CME::Vec3 dynamicRay={std::cos(dynamicPsi),std::sin(dynamicPsi),0};
    auto dynamicState=[&](double time) {
      SF::SseKinematicState state;
      state.apexRadiusM=1e10+8e5*time;state.apexSpeedMPerS=8e5;
      state.direction={std::cos(dynamicOmega*time),std::sin(dynamicOmega*time),0};
      state.directionRatePerS=dynamicOmega*CME::Vec3{-std::sin(dynamicOmega*time),
          std::cos(dynamicOmega*time),0};
      state.halfWidthRad=dynamicLambda+dynamicLambdaRate*time;
      state.halfWidthRateRadPerS=dynamicLambdaRate;return state;
    };
    const auto dynamic=SF::EvaluateSseRay(dynamicState(0),dynamicRay);
    const double dynamicH=1e-3;
    const auto dynamicPlus=SF::EvaluateSseRay(dynamicState(dynamicH),dynamicRay);
    const auto dynamicMinus=SF::EvaluateSseRay(dynamicState(-dynamicH),dynamicRay);
    const CME::Vec3 intersectionVelocity=(dynamicPlus.value.positionM-
        dynamicMinus.value.positionM)/(2*dynamicH);
    auto omitWidth=dynamicState(0);omitWidth.halfWidthRateRadPerS=0;
    auto omitDirection=dynamicState(0);omitDirection.directionRatePerS={};
    const auto noWidth=SF::EvaluateSseRay(omitWidth,dynamicRay);
    const auto noDirection=SF::EvaluateSseRay(omitDirection,dynamicRay);
    const CME::Vec3 dynamicEdgeRay={std::cos(dynamicLambda),std::sin(dynamicLambda),0};
    auto widthOnly=dynamicState(0);widthOnly.directionRatePerS={};
    const auto dynamicEdge=SF::EvaluateSseRay(widthOnly,dynamicEdgeRay);
    const double edgeReference=widthOnly.apexRadiusM*std::cos(dynamicLambda)*
        dynamicLambdaRate/(1+std::sin(dynamicLambda));
    auto invalidWidth=dynamicState(0);invalidWidth.halfWidthRad=0;
    auto invalidDirection=dynamicState(0);invalidDirection.directionRatePerS={1,0,0};
    std::ostringstream dynamicEvidence;dynamicEvidence<<std::setprecision(17)<<" fd="<<
        CME::Dot(intersectionVelocity,dynamic.value.outwardNormal)<<" exact="<<
        dynamic.value.normalSpeedMPerS<<" fd_delta="<<std::abs(CME::Dot(
        intersectionVelocity,dynamic.value.outwardNormal)-dynamic.value.normalSpeedMPerS)<<
        " omit_width_delta="<<
        std::abs(noWidth.value.normalSpeedMPerS-dynamic.value.normalSpeedMPerS)<<
        " omit_direction_delta="<<std::abs(noDirection.value.normalSpeedMPerS-
        dynamic.value.normalSpeedMPerS)<<" edge="<<dynamicEdge.value.normalSpeedMPerS<<
        " edge_ref="<<edgeReference;
    Check(dynamic.ok()&&dynamicPlus.ok()&&dynamicMinus.ok()&&noWidth.ok()&&
        noDirection.ok()&&dynamicEdge.ok()&&Near(CME::Dot(intersectionVelocity,
          dynamic.value.outwardNormal),dynamic.value.normalSpeedMPerS,2e-8)&&
        std::abs(noWidth.value.normalSpeedMPerS-dynamic.value.normalSpeedMPerS)>1e4&&
        std::abs(noDirection.value.normalSpeedMPerS-dynamic.value.normalSpeedMPerS)>1e4&&
        Near(dynamicEdge.value.normalSpeedMPerS,edgeReference,2e-12)&&
        !SF::EvaluateSseRay(invalidWidth,dynamicRay).ok()&&
        !SF::EvaluateSseRay(invalidDirection,dynamicRay).ok(),"RSH32",
        "complete SSE width/direction derivatives match independent motion; omitted terms and interval/tangency controls remain live "+dynamicEvidence.str());

    const auto endpoint=longProvider->EndpointTimeS();
    const auto endpointState=longProvider->Trajectory(endpoint.value);
    Check(endpoint.ok()&&Near(endpoint.value,203884.494,0,0.2)&&endpointState.ok()&&
        Near(endpointState.value.apexRadiusM,longEvent->endpointRadiusM,0,2.0),"RSH28",
        "analytical drag trajectory is evaluated at the actual 1-AU root");

    const double twoSolar=2*longEvent->ambient.support.solarRadiusM;
    const auto polar=longProvider->QueryAmbient({0,0,twoSolar},0,1);
    std::ostringstream polarEvidence;
    if(polar.ok())polarEvidence<<" measured_U="
        <<CME::Norm(polar.value.primitive.velocityMPerS)<<" measured_vA="
        <<polar.value.primitive.plasma.alfvenSpeedMPerS;
    Check(polar.ok()&&Near(CME::Norm(polar.value.primitive.velocityMPerS),49334,4e-3)&&
        Near(polar.value.primitive.plasma.alfvenSpeedMPerS,334155,5e-3),"RSH31",
        "canonical polar ambient reproduces Parker/PFSS wind and Alfven scales"+
        polarEvidence.str());
    if(polar.ok()) {
      CME::MhdPrimitiveState state{polar.value.primitive.plasma.massDensityKgM3,
          polar.value.primitive.plasma.pressurePa,polar.value.primitive.velocityMPerS,
          polar.value.primitive.magneticFieldT};
      const auto waves=CME::EvaluateMhdCharacteristics(state,{0,0,1},
          longEvent->ambient.composition.gammaAdiabatic);
      Check(waves.ok()&&Near(waves.value.fastSpeedMPerS,
          std::max(waves.value.soundSpeedMPerS,waves.value.alfvenSpeedMPerS),1e-12),
          "RSH14","parallel polar fast speed uses max(cs,vA), not perpendicular sum");
    }

    // RSH15: a Galilean boost changes U and the scalar normal surface speed
    // together, leaving shock-frame inflow and all jump invariants unchanged.
    // A rigid HCI rotation likewise only rotates vector results.  Merely
    // reversing the designated outward normal without swapping the upstream
    // side is not an invariance: it must be rejected as non-forward.
    const double gamma=5.0/3.0;
    const double mu0=CME::Constants::kVacuumPermeabilityHPerM;
    const CME::Vec3 frameNormal=CME::Unit({1,2,3});
    CME::MhdPrimitiveState frameUpstream{2,6,{4,-3,2},
        std::sqrt(mu0)*CME::Vec3{0.3,0.5,-0.2}};
    const double frameSpeed=20;
    const auto frameJump=CME::SolveObliqueFastShock(frameUpstream,frameNormal,
        frameSpeed,gamma,1e-9);
    const CME::Vec3 boost={8,-5,3};
    auto boostedUpstream=frameUpstream;
    boostedUpstream.velocityMPerS=boostedUpstream.velocityMPerS+boost;
    const auto boostedJump=CME::SolveObliqueFastShock(boostedUpstream,frameNormal,
        frameSpeed+CME::Dot(boost,frameNormal),gamma,1e-9);
    const double frameAngle=0.63;
    auto rotatedUpstream=frameUpstream;
    rotatedUpstream.velocityMPerS=RotateZ(frameUpstream.velocityMPerS,frameAngle);
    rotatedUpstream.magneticFieldT=RotateZ(frameUpstream.magneticFieldT,frameAngle);
    const auto rotatedJump=CME::SolveObliqueFastShock(rotatedUpstream,
        RotateZ(frameNormal,frameAngle),frameSpeed,gamma,1e-9);
    const auto reversed=SF::EvaluateLocalJump(frameUpstream,-1*frameNormal,
        -frameSpeed,gamma,1e-6,1e-9,0);
    Check(frameJump.ok()&&boostedJump.ok()&&rotatedJump.ok()&&reversed.ok()&&
        Near(boostedJump.value.upstreamFastMach,frameJump.value.upstreamFastMach,2e-13)&&
        Near(boostedJump.value.compressionRatio,frameJump.value.compressionRatio,2e-12)&&
        RelativeVectorError(boostedJump.value.downstream.velocityMPerS,
          frameJump.value.downstream.velocityMPerS+boost)<2e-12&&
        RelativeVectorError(rotatedJump.value.downstream.velocityMPerS,
          RotateZ(frameJump.value.downstream.velocityMPerS,frameAngle))<2e-12&&
        RelativeVectorError(rotatedJump.value.downstream.magneticFieldT,
          RotateZ(frameJump.value.downstream.magneticFieldT,frameAngle))<2e-12&&
        reversed.value.status==SF::FrontStatus::NonForwardInflow&&
        !reversed.value.downstreamValid,"RSH15",
        "shock invariants survive Galilean/rigid-frame changes while an unswapped normal reversal is typed non-forward");

    CME::MhdPrimitiveState gas{1,1/gamma,{0,0,0},{0,0,0}};
    const auto gasJump=CME::SolveObliqueFastShock(gas,{1,0,0},3,gamma,1e-9);
    Check(gasJump.ok()&&Near(gasJump.value.compressionRatio,3,1e-8)&&
        Near(gasJump.value.downstream.pressurePa/gas.pressurePa,11,1e-8)&&
        Near(gasJump.value.downstream.velocityMPerS.x,2,1e-8)&&
        Near((gasJump.value.downstream.pressurePa/gasJump.value.downstream.massDensityKgM3)/
          (gas.pressurePa/gas.massDensityKgM3),11.0/3.0,1e-8),"RSH16",
        "production RH solver recovers X=3, p2/p1=11, U2n=2 and T2/T1=11/3 for Ms=3");
    CME::MhdPrimitiveState perpendicular{1,0.1,{0,0,0},{0,std::sqrt(mu0),0}};
    const auto mhdJump=CME::SolveObliqueFastShock(perpendicular,{1,0,0},2,gamma,1e-9);
    const double energy1=ReferenceEnergyFlux(perpendicular,{1,0,0},2,gamma);
    const double energy2=mhdJump.ok()?ReferenceEnergyFlux(mhdJump.value.downstream,
        {1,0,0},2,gamma):0;
    const double expectedEntropy=std::log(6.0)-gamma*std::log(2.0);
    Check(mhdJump.ok()&&Near(mhdJump.value.compressionRatio,2,2e-8)&&
        Near(mhdJump.value.downstream.pressurePa/perpendicular.pressurePa,6,2e-8)&&
        Near(CME::Norm(mhdJump.value.downstream.magneticFieldT)/
        CME::Norm(perpendicular.magneticFieldT),2,2e-8)&&
        Near(mhdJump.value.entropyLogIncrement,expectedEntropy,2e-8)&&
        std::abs(energy2-energy1)/std::max(std::abs(energy1),1.0)<1e-9,
        "RSH17","finite-beta perpendicular reference closes independent total energy and entropy as well as X/B/p ratios");

    // RSH18: cover exact parallel, oblique and perpendicular limiting
    // families.  The negative control deliberately multiplies the *normal*
    // magnetic component by density compression; it must violate normal-flux
    // conservation even though that shortcut is valid for perpendicular Bt.
    CME::MhdPrimitiveState parallel{1,1/gamma,{0,0.2,-0.1},
        {0.5*std::sqrt(mu0),0,0}};
    CME::MhdPrimitiveState oblique{1,1/gamma,{0,0.2,-0.1},
        {0.4*std::sqrt(mu0),0.7*std::sqrt(mu0),0}};
    const auto parallelJump=CME::SolveObliqueFastShock(parallel,{1,0,0},3,gamma,1e-9);
    const auto obliqueJump=CME::SolveObliqueFastShock(oblique,{1,0,0},3.5,gamma,1e-9);
    double badNormalResidual=0,obliqueMassResidual=1,obliqueElectricResidual=1;
    if(obliqueJump.ok()) {
      const auto& down=obliqueJump.value.downstream;
      const CME::Vec3 n={1,0,0};
      const double m1=oblique.massDensityKgM3*CME::Dot(
          oblique.velocityMPerS-3.5*n,n);
      const double m2=down.massDensityKgM3*CME::Dot(
          down.velocityMPerS-3.5*n,n);
      obliqueMassResidual=std::abs(m2-m1)/std::max(std::abs(m1),std::abs(m2));
      const CME::Vec3 e1=CME::Cross(oblique.velocityMPerS-3.5*n,
          oblique.magneticFieldT);
      const CME::Vec3 e2=CME::Cross(down.velocityMPerS-3.5*n,
          down.magneticFieldT);
      // Faraday's jump condition constrains the tangential electric field.
      // The normal component -(u_t x B_t).n is not one of the 1-D RH fluxes
      // and need not agree across a discontinuity.
      const CME::Vec3 de=e2-e1;
      const CME::Vec3 deTangential=de-CME::Dot(de,n)*n;
      const CME::Vec3 e1Tangential=e1-CME::Dot(e1,n)*n;
      const CME::Vec3 e2Tangential=e2-CME::Dot(e2,n)*n;
      obliqueElectricResidual=CME::Norm(deTangential)/
          std::max(CME::Norm(e1Tangential),CME::Norm(e2Tangential));
      const CME::Vec3 badField=obliqueJump.value.compressionRatio*
          oblique.magneticFieldT;
      badNormalResidual=std::abs(CME::Dot(badField-oblique.magneticFieldT,n)) /
          std::max(CME::Norm(badField),CME::Norm(oblique.magneticFieldT));
    }
    std::ostringstream jumpEvidence;jumpEvidence<<std::setprecision(17)<<
        " mass="<<obliqueMassResidual<<" Et="<<obliqueElectricResidual<<
        " bad_Bn="<<badNormalResidual<<" parallel_ok="<<parallelJump.ok()<<
        " oblique_ok="<<obliqueJump.ok();
    Check(parallelJump.ok()&&obliqueJump.ok()&&mhdJump.ok()&&
        parallelJump.value.branch=="regular-parallel-fast"&&
        Near(CME::Dot(parallelJump.value.downstream.magneticFieldT,{1,0,0}),
          CME::Dot(parallel.magneticFieldT,{1,0,0}),2e-13)&&
        obliqueJump.value.branch=="regular-oblique-fast"&&
        obliqueMassResidual<1e-12&&obliqueElectricResidual<1e-9&&
        obliqueJump.value.residuals.maximum<1e-9&&badNormalResidual>0.1,
        "RSH18","parallel/oblique/perpendicular production branches conserve independent fluxes and reject full-vector B multiplication by X"+jumpEvidence.str());

    // RSH19 uses the same production classifier called by Prepare.  A small
    // positive Mach excess whose physical root lies below the maintained
    // compression bracket is numerical-unknown, not sub-fast and not a
    // fabricated finite jump.
    const auto gasWaves=CME::EvaluateMhdCharacteristics(gas,{1,0,0},gamma);
    const auto subfast=SF::EvaluateLocalJump(gas,{1,0,0},
        0.999*gasWaves.value.fastSpeedMPerS,gamma,1e-6,1e-9,0);
    const auto unresolved=SF::EvaluateLocalJump(gas,{1,0,0},
        (1+1e-10)*gasWaves.value.fastSpeedMPerS,gamma,1e-6,1e-9,0);
    const auto resolvedWeak=SF::EvaluateLocalJump(gas,{1,0,0},
        1.01*gasWaves.value.fastSpeedMPerS,gamma,1e-6,1e-8,0);
    const auto nonForward=SF::EvaluateLocalJump(gas,{1,0,0},-1,
        gamma,1e-6,1e-9,0);
    Check(subfast.ok()&&unresolved.ok()&&resolvedWeak.ok()&&nonForward.ok()&&
        subfast.value.status==SF::FrontStatus::SubfastFront&&
        !subfast.value.downstreamValid&&
        unresolved.value.fastMach>1&&
        unresolved.value.status==SF::FrontStatus::NumericallyUnresolvedWeakShock&&
        !unresolved.value.downstreamValid&&
        resolvedWeak.value.status==SF::FrontStatus::SolvedFastShock&&
        resolvedWeak.value.jump.compressionRatio>1&&
        nonForward.value.status==SF::FrontStatus::NonForwardInflow&&
        !nonForward.value.downstreamValid,"RSH19",
        "sub-fast/non-forward absence, resolvable weak jumps and positive-Mach unresolved failures remain distinct without floors");

    // RSH36: evaluate diagnostics through the production local-jump path.
    // Exact perpendicularity has no finite conventional HT boost, an exact
    // null has no magnetic diagnostics, and a near-perpendicular state keeps
    // its large condition number visible rather than imposing an angle floor.
    const auto perpendicularDiagnostic=SF::EvaluateLocalJump(perpendicular,
        {1,0,0},2,gamma,1e-6,1e-9,0);
    const auto parallelDiagnostic=SF::EvaluateLocalJump(parallel,{1,0,0},
        3,gamma,1e-6,1e-9,0);
    const auto obliqueDiagnostic=SF::EvaluateLocalJump(oblique,{1,0,0},
        3.5,gamma,1e-6,1e-9,0);
    const auto nullDiagnostic=SF::EvaluateLocalJump(gas,{1,0,0},3,
        gamma,1e-6,1e-9,0);
    CME::MhdPrimitiveState nearPerpendicular=perpendicular;
    nearPerpendicular.magneticFieldT.x=1e-8*std::sqrt(mu0);
    const auto nearDiagnostic=SF::EvaluateLocalJump(nearPerpendicular,
        {1,0,0},2,gamma,1e-6,1e-8,0);
    const CME::Vec3 diagnosticBoost={4,-7,2};
    auto shiftedOblique=oblique;
    shiftedOblique.velocityMPerS=shiftedOblique.velocityMPerS+diagnosticBoost;
    const auto shiftedDiagnostic=SF::EvaluateLocalJump(shiftedOblique,{1,0,0},
        3.5+diagnosticBoost.x,gamma,1e-6,1e-9,0);
    Check(perpendicularDiagnostic.ok()&&parallelDiagnostic.ok()&&
        obliqueDiagnostic.ok()&&nullDiagnostic.ok()&&nearDiagnostic.ok()&&
        shiftedDiagnostic.ok()&&
        perpendicularDiagnostic.value.diagnostics.magneticCompressionValid&&
        Near(perpendicularDiagnostic.value.diagnostics.magneticCompression,2,2e-8)&&
        perpendicularDiagnostic.value.diagnostics.htStatus==
          SF::HtDiagnosticStatus::PerpendicularNoFiniteBoost&&
        parallelDiagnostic.value.diagnostics.htStatus==SF::HtDiagnosticStatus::Valid&&
        Near(parallelDiagnostic.value.diagnostics.magneticCompression,1,2e-13)&&
        obliqueDiagnostic.value.diagnostics.htStatus==SF::HtDiagnosticStatus::Valid&&
        obliqueDiagnostic.value.diagnostics.upstreamElectricCancellation<1e-12&&
        obliqueDiagnostic.value.diagnostics.downstreamElectricCancellation<1e-9&&
        Near(obliqueDiagnostic.value.diagnostics.incidentHtSpeedMPerS,
          obliqueDiagnostic.value.inflowMPerS/
          obliqueDiagnostic.value.diagnostics.absoluteNormalFieldCosine,2e-12)&&
        nullDiagnostic.value.diagnostics.htStatus==
          SF::HtDiagnosticStatus::UpstreamMagneticNull&&
        !nullDiagnostic.value.diagnostics.magneticCompressionValid&&
        nearDiagnostic.value.diagnostics.htStatus==SF::HtDiagnosticStatus::Valid&&
        nearDiagnostic.value.diagnostics.absoluteNormalFieldCosine<2e-8&&
        shiftedDiagnostic.value.diagnostics.htStatus==SF::HtDiagnosticStatus::Valid&&
        RelativeVectorError(shiftedDiagnostic.value.diagnostics.htTangentialBoostMPerS,
          obliqueDiagnostic.value.diagnostics.htTangentialBoostMPerS+
          CME::Vec3{0,diagnosticBoost.y,diagnosticBoost.z})<1e-10&&
        RelativeVectorError(shiftedDiagnostic.value.diagnostics.upstreamVelocityHtMPerS,
          obliqueDiagnostic.value.diagnostics.upstreamVelocityHtMPerS)<1e-10,
        "RSH36","magnetic compression and HT vectors preserve frames/electric cancellation while perpendicular/null/ill-conditioned limits stay explicit");

    const auto first=provider->Prepare(0,1);
    const auto prior=provider->Current();
    const auto failed=provider->Prepare(-1,2);
    Check(first.ok()&&!failed.ok()&&provider->Current()==prior&&
        provider->Current()->generation==1,"RSH24",
        "failed candidate retains committed epoch and generation first_status="+
        first.status.message);
    // The production surface is a triangulated topological disk, not a set of
    // quadrature centroids connected only for plotting.  Reconstruct edge
    // incidence independently from connectivity: the finite SSE support rim
    // is the only boundary loop, the shared apex is interior, every chord is
    // outward, and physical records are in one-to-one correspondence with
    // exact-curved-area facets.
    bool triangleTopology=first.ok();
    std::map<std::pair<std::uint32_t,std::uint32_t>,int> edgeIncidence;
    std::vector<std::vector<std::uint32_t>> adjacency;
    std::set<std::uint32_t> boundaryVertices;
    std::size_t apexVertices=0;
    long double triangleArea=0;
    if(first.ok()) {
      const auto& mesh=*first.value;
      adjacency.resize(mesh.vertices.size());
      triangleTopology=mesh.vertices.size()==
          static_cast<std::size_t>(smoke->polarCells*smoke->azimuthCells+1)&&
          mesh.triangles.size()==static_cast<std::size_t>(
              (2*smoke->polarCells-1)*smoke->azimuthCells)&&
          mesh.records.size()==mesh.triangles.size();
      const double meshSine=std::sin(smoke->halfWidthRad);
      const double meshCenter=smoke->initialApexRadiusM/(1+meshSine);
      const double meshRadius=meshCenter*meshSine;
      const CME::Vec3 meshOrigin=meshCenter*smoke->direction;
      for(const auto& vertex:mesh.vertices) {
        apexVertices+=vertex.apex;
        triangleTopology=triangleTopology&&
            Near(CME::Norm(vertex.positionM-meshOrigin),meshRadius,4e-14);
      }
      for(std::size_t face=0;face<mesh.triangles.size();++face) {
        const auto& triangle=mesh.triangles[face];
        const auto& record=mesh.records[face];
        triangleTopology=triangleTopology&&triangle.stableId==face+1&&
            record.geometry.stableId==triangle.stableId&&
            record.geometry.areaM2==triangle.curvedAreaM2&&
            triangle.vertex[0]<mesh.vertices.size()&&
            triangle.vertex[1]<mesh.vertices.size()&&
            triangle.vertex[2]<mesh.vertices.size()&&
            triangle.vertex[0]!=triangle.vertex[1]&&
            triangle.vertex[1]!=triangle.vertex[2]&&
            triangle.vertex[2]!=triangle.vertex[0];
        if(!triangleTopology)break;
        const auto& x0=mesh.vertices[triangle.vertex[0]].positionM;
        const auto& x1=mesh.vertices[triangle.vertex[1]].positionM;
        const auto& x2=mesh.vertices[triangle.vertex[2]].positionM;
        const auto cross=CME::Cross(x1-x0,x2-x0);
        triangleTopology=triangleTopology&&
            CME::Dot(cross,record.geometry.outwardNormal)>0&&
            Near(0.5*CME::Norm(cross),triangle.planarAreaM2,3e-15);
        triangleArea+=triangle.curvedAreaM2;
        for(int edge=0;edge<3;++edge) {
          const std::uint32_t a=triangle.vertex[edge];
          const std::uint32_t b=triangle.vertex[(edge+1)%3];
          edgeIncidence[std::minmax(a,b)]++;
          adjacency[a].push_back(b);adjacency[b].push_back(a);
        }
      }
      for(const auto& edge:edgeIncidence) {
        triangleTopology=triangleTopology&&
            (edge.second==1||edge.second==2);
        if(edge.second==1) {
          boundaryVertices.insert(edge.first.first);
          boundaryVertices.insert(edge.first.second);
        }
      }
      std::vector<bool> visited(mesh.vertices.size(),false);
      std::vector<std::uint32_t> pending{0};visited[0]=true;
      for(std::size_t at=0;at<pending.size();++at)
        for(std::uint32_t next:adjacency[pending[at]])if(!visited[next]) {
          visited[next]=true;pending.push_back(next);
        }
      const std::size_t boundaryEdges=std::count_if(edgeIncidence.begin(),
          edgeIncidence.end(),[](const auto& item){return item.second==1;});
      const auto apexIndex=static_cast<std::uint32_t>(mesh.vertices.size()-1);
      triangleTopology=triangleTopology&&apexVertices==1&&
          boundaryEdges==static_cast<std::size_t>(smoke->azimuthCells)&&
          boundaryVertices.size()==static_cast<std::size_t>(smoke->azimuthCells)&&
          !boundaryVertices.count(apexIndex)&&
          std::all_of(visited.begin(),visited.end(),[](bool value){return value;})&&
          static_cast<long long>(mesh.vertices.size())-
              static_cast<long long>(edgeIncidence.size())+
              static_cast<long long>(mesh.triangles.size())==1&&
          Near(static_cast<double>(triangleArea),
              mesh.area.geometricSupportM2,3e-15);
    }
    Check(triangleTopology,"RSH23",
        "triangular SSE disk has one apex, one finite-support boundary, outward nondegenerate faces, exact patch/record ownership and Euler characteristic one");

    // At this production epoch one complete azimuthal face ring is only
    // M_f-1 ~= 1.06e-3 above the fast characteristic.  Triangular face
    // centroids sample that ring whereas the retired quadrilateral centers
    // did not.  Exercise the real positive event across the transition so a
    // flat-residual RH root cannot masquerade as topology-dependent coverage
    // loss.  The last epoch must change those faces to a physical sub-fast
    // classification, not retain or fabricate a downstream state.
    const auto positiveProvider=SF::Provider::Create(positiveEvent);
    const auto weakBefore=positiveProvider.ok()?
        positiveProvider.value->Prepare(56400,95):
        SEP::Core::Result<std::shared_ptr<const SF::Epoch>>::Failure(
          SEP::Core::StatusCode::NumericalFailure,"positive provider unavailable");
    const auto weakCrossing=positiveProvider.ok()?
        positiveProvider.value->Prepare(57000,96):
        SEP::Core::Result<std::shared_ptr<const SF::Epoch>>::Failure(
          SEP::Core::StatusCode::NumericalFailure,"positive provider unavailable");
    const auto weakAfter=positiveProvider.ok()?
        positiveProvider.value->Prepare(57600,97):
        SEP::Core::Result<std::shared_ptr<const SF::Epoch>>::Failure(
          SEP::Core::StatusCode::NumericalFailure,"positive provider unavailable");
    const auto acceptedCount=[](const std::shared_ptr<const SF::Epoch>& epoch) {
      return std::count_if(epoch->records.begin(),epoch->records.end(),
          [](const SF::ShockRecord& record) {
            return record.status==SF::FrontStatus::SolvedFastShock;
          });
    };
    Check(positiveProvider.ok()&&weakBefore.ok()&&weakCrossing.ok()&&
        weakAfter.ok()&&weakBefore.value->area.numericalFailureM2==0&&
        weakCrossing.value->area.numericalFailureM2==0&&
        weakAfter.value->area.numericalFailureM2==0&&
        acceptedCount(weakBefore.value)==1344&&
        acceptedCount(weakCrossing.value)==1344&&
        acceptedCount(weakAfter.value)==1296,"RSH23",
        "actual triangular production faces resolve the weak-shock ring and then classify its physical sub-fast transition without unknown area");
    auto strictRhConfiguration=std::make_shared<SF::Configuration>(*longEvent);
    strictRhConfiguration->rhResidualTolerance=1e-30;
    strictRhConfiguration->validityPolicy=
        SF::TrajectoryValidityPolicy::RequireDeclaredShockCoverage;
    strictRhConfiguration->requireFastShockAtObserver=false;
    const auto strictRhProvider=SF::Provider::Create(strictRhConfiguration);
    const auto strictRhInitial=strictRhProvider.ok()?strictRhProvider.value->Prepare(0,1):
        SEP::Core::Result<std::shared_ptr<const SF::Epoch>>::Failure(
          SEP::Core::StatusCode::NumericalFailure,"strict RH provider unavailable");
    const auto strictRhCommitted=strictRhProvider.ok()?strictRhProvider.value->Current():nullptr;
    const auto strictRhFailed=strictRhProvider.ok()?strictRhProvider.value->Prepare(600,2):
        SEP::Core::Result<std::shared_ptr<const SF::Epoch>>::Failure(
          SEP::Core::StatusCode::NumericalFailure,"strict RH provider unavailable");
    std::string brokenAsset=Bytes(examples/"handoff_smoke.event");
    const auto brokenChecksum=brokenAsset.find("assets.harmonics_sha256=");
    brokenAsset[brokenChecksum+sizeof("assets.harmonics_sha256=")-1]='0';
    const auto assetFailure=SF::ResolveConfiguration(brokenAsset,[&](const std::string& name){
      return SEP::Core::Result<std::string>::Success(Bytes(examples/name));});
    Check(strictRhProvider.ok()&&strictRhInitial.ok()&&!strictRhFailed.ok()&&
        strictRhProvider.value->Current()==strictRhCommitted&&
        strictRhProvider.value->Current()->generation==1&&!assetFailure.ok()&&
        provider->Current()==prior,"RSH24",
        "injected RH and asset-resolution failures cannot change either provider's committed epoch/clock");

    // RSH21: intersect independent straight-line field references with the
    // actual finite production front.  Quadratic segment roots retain all
    // supported crossings, catch a double tangent, deduplicate a shared
    // vertex, reject rear/support misses, and attach canonical polarity/jump
    // state rather than a reference-only placeholder.
    const double intersectionTime=1000;
    const CME::Vec3 intersectionRay=CME::Unit({std::sin(0.3),0,std::cos(0.3)});
    const auto intersectionGeometry=longProvider->EvaluateRay(
        intersectionRay,intersectionTime,701);
    const auto intersectionTrajectory=longProvider->Trajectory(intersectionTime);
    const double intersectionSphereRadius=intersectionTrajectory.ok()?
        intersectionTrajectory.value.apexRadiusM*std::sin(longEvent->halfWidthRad)/
        (1+std::sin(longEvent->halfWidthRad)):0;
    const CME::Vec3 intersectionNormal=intersectionGeometry.ok()?
        intersectionGeometry.value.outwardNormal:CME::Vec3{1,0,0};
    const CME::Vec3 tangentDirection=CME::Unit(CME::Cross(intersectionNormal,
        std::abs(intersectionNormal.z)<0.9?CME::Vec3{0,0,1}:CME::Vec3{0,1,0}));
    const CME::Vec3 intersectionPoint=intersectionGeometry.ok()?
        intersectionGeometry.value.positionM:CME::Vec3{};
    const auto transverseRoots=SF::IntersectPolylineWithFront(*longProvider,
        {intersectionPoint+2*intersectionSphereRadius*intersectionNormal,
         intersectionPoint,
         intersectionPoint-2*intersectionSphereRadius*intersectionNormal},
        intersectionTime,1e-3);
    const auto tangentRoots=SF::IntersectPolylineWithFront(*longProvider,
        {intersectionPoint-intersectionSphereRadius*tangentDirection,
         intersectionPoint+intersectionSphereRadius*tangentDirection},
        intersectionTime,1e-3);
    bool rootPayload=true;
    if(transverseRoots.ok())for(const auto& rootHit:transverseRoots.value)
      rootPayload=rootPayload&&rootHit.arcLengthM>=0&&
          rootHit.front.geometry.stableId>0&&rootHit.signedPolarityAlongTrace!=0&&
          rootHit.front.status!=SF::FrontStatus::OutsideFrontSupport;
    std::ostringstream intersectionEvidence;intersectionEvidence<<
        " transverse_ok="<<transverseRoots.ok()<<" transverse_n="<<
        (transverseRoots.ok()?transverseRoots.value.size():0)<<" tangent_ok="<<
        tangentRoots.ok()<<" tangent_n="<<(tangentRoots.ok()?tangentRoots.value.size():0)<<
        " payload="<<rootPayload<<" polarity="<<(transverseRoots.ok()&&
        !transverseRoots.value.empty()?transverseRoots.value.front().signedPolarityAlongTrace:99)<<
        " tangent_flag="<<(tangentRoots.ok()&&!tangentRoots.value.empty()?
        tangentRoots.value.front().tangent:false)<<" tangent_dx="<<
        (tangentRoots.ok()&&!tangentRoots.value.empty()?CME::Norm(
          tangentRoots.value.front().positionM-intersectionPoint):-1);
    Check(intersectionGeometry.ok()&&intersectionTrajectory.ok()&&
        transverseRoots.ok()&&!transverseRoots.value.empty()&&rootPayload&&
        tangentRoots.ok()&&tangentRoots.value.size()==1&&
        tangentRoots.value.front().tangent&&
        CME::Norm(tangentRoots.value.front().positionM-intersectionPoint)<1e-3,
        "RSH21","production finite-front intersections retain supported roots/polarity, deduplicate a shared vertex and recover an analytic tangent"+intersectionEvidence.str());

    // RSH22: a fixed observer on a regular ray gives a transverse crossing;
    // the fixed-width support edge gives a true geometric graze because Vn=0.
    // A moving observer reports Vn-vobs.n, and a point outside finite support
    // remains a miss even if the generating sphere has unrelated roots.
    const double observerHitTime=1200;
    const auto observerGeometry=longProvider->EvaluateRay(intersectionRay,
        observerHitTime,801);
    const CME::Vec3 observerVelocity=2e4*intersectionNormal;
    const auto fixedObserverCoarse=observerGeometry.ok()?
        SF::FindObserverPassages(*longProvider,observerGeometry.value.positionM,{},0,
          observerHitTime-100,observerHitTime+100,20,1e-6,2):
        SEP::Core::Result<std::vector<SF::ObserverPassage>>::Failure(
          SEP::Core::StatusCode::NumericalFailure,"observer geometry unavailable");
    const auto fixedObserverFine=observerGeometry.ok()?
        SF::FindObserverPassages(*longProvider,observerGeometry.value.positionM,{},0,
          observerHitTime-100,observerHitTime+100,10,2.5e-7,0.5):
        SEP::Core::Result<std::vector<SF::ObserverPassage>>::Failure(
          SEP::Core::StatusCode::NumericalFailure,"observer geometry unavailable");
    const auto movingObserver=observerGeometry.ok()?
        SF::FindObserverPassages(*longProvider,
          observerGeometry.value.positionM-observerHitTime*observerVelocity,
          observerVelocity,0,observerHitTime-100,observerHitTime+100,10,1e-6,2):
        SEP::Core::Result<std::vector<SF::ObserverPassage>>::Failure(
          SEP::Core::StatusCode::NumericalFailure,"observer geometry unavailable");
    const CME::Vec3 edgeRay=CME::Unit({std::sin(longEvent->halfWidthRad),0,
        std::cos(longEvent->halfWidthRad)});
    const auto edgeGeometry=longProvider->EvaluateRay(edgeRay,observerHitTime,802);
    const auto edgeObserver=edgeGeometry.ok()?SF::FindObserverPassages(*longProvider,
        edgeGeometry.value.positionM,{},0,observerHitTime-100,observerHitTime+100,
        10,1e-6,2):SEP::Core::Result<std::vector<SF::ObserverPassage>>::Failure(
          SEP::Core::StatusCode::NumericalFailure,"edge geometry unavailable");
    const CME::Vec3 missPosition=CME::Norm(observerGeometry.value.positionM)*
        CME::Vec3{0,1,0};
    const auto missedObserver=SF::FindObserverPassages(*longProvider,missPosition,{},0,
        observerHitTime-100,observerHitTime+100,10,1e-6,2);
    const bool observerGood=fixedObserverCoarse.ok()&&fixedObserverFine.ok()&&
        !fixedObserverCoarse.value.empty()&&!fixedObserverFine.value.empty()&&
        std::abs(fixedObserverCoarse.value.front().timeS-observerHitTime)<2e-5&&
        std::abs(fixedObserverFine.value.front().timeS-observerHitTime)<5e-6&&
        std::abs(fixedObserverFine.value.front().timeS-
          fixedObserverCoarse.value.front().timeS)<3e-5&&
        fixedObserverFine.value.front().kind==SF::ObserverPassageKind::Transverse&&
        movingObserver.ok()&&!movingObserver.value.empty()&&
        Near(movingObserver.value.front().relativeNormalSpeedMPerS,
          movingObserver.value.front().front.geometry.normalSpeedMPerS-
          CME::Dot(observerVelocity,
            movingObserver.value.front().front.geometry.outwardNormal),2e-13)&&
        edgeObserver.ok()&&!edgeObserver.value.empty()&&
        edgeObserver.value.front().kind==SF::ObserverPassageKind::Grazing&&
        missedObserver.ok()&&missedObserver.value.empty();
    std::ostringstream observerEvidence;observerEvidence<<" coarse_ok="<<
        fixedObserverCoarse.ok()<<" coarse_n="<<(fixedObserverCoarse.ok()?
        fixedObserverCoarse.value.size():0)<<" fine_ok="<<fixedObserverFine.ok()<<
        " fine_n="<<(fixedObserverFine.ok()?fixedObserverFine.value.size():0)<<
        " moving_ok="<<movingObserver.ok()<<" moving_n="<<(movingObserver.ok()?
        movingObserver.value.size():0)<<" edge_ok="<<edgeObserver.ok()<<
        " edge_n="<<(edgeObserver.ok()?edgeObserver.value.size():0)<<" miss_ok="<<
        missedObserver.ok()<<" miss_n="<<(missedObserver.ok()?missedObserver.value.size():0);
    Check(observerGood,"RSH22",
        "bracketed observer events converge, use relative normal speed, retain a support-edge graze and distinguish a finite-support miss"+observerEvidence.str());

    // RSH29: the reduced checkpoint contains only immutable identities and
    // analytical epoch state.  A restored provider must reproduce the later
    // front/jump/output bytes exactly; a changed physical identity is rejected
    // before it can alter the restored epoch.
    const auto restartA=SF::Provider::Create(longEvent);
    const auto checkpointEpoch=restartA.ok()?restartA.value->Prepare(600,11):
        SEP::Core::Result<std::shared_ptr<const SF::Epoch>>::Failure(
          SEP::Core::StatusCode::NumericalFailure,"restart provider unavailable");
    const auto checkpoint=restartA.ok()?SF::MakeRestartState(*restartA.value):
        SEP::Core::Result<SF::ReducedRestartState>::Failure(
          SEP::Core::StatusCode::NumericalFailure,"restart provider unavailable");
    const auto uninterrupted=restartA.ok()?restartA.value->Prepare(1800,31):
        SEP::Core::Result<std::shared_ptr<const SF::Epoch>>::Failure(
          SEP::Core::StatusCode::NumericalFailure,"restart provider unavailable");
    const auto restartB=SF::Provider::Create(longEvent);
    const auto restored=restartB.ok()&&checkpoint.ok()?
        SF::RestoreRestartState(restartB.value.get(),checkpoint.value):
        SEP::Core::Result<std::shared_ptr<const SF::Epoch>>::Failure(
          SEP::Core::StatusCode::NumericalFailure,"restart input unavailable");
    const auto resumed=restartB.ok()?restartB.value->Prepare(1800,31):
        SEP::Core::Result<std::shared_ptr<const SF::Epoch>>::Failure(
          SEP::Core::StatusCode::NumericalFailure,"restart provider unavailable");
    auto changedRestart=checkpoint.ok()?checkpoint.value:SF::ReducedRestartState{};
    changedRestart.eventIdentity="changed-physical-input";
    const auto restartPrior=restartB.ok()?restartB.value->Current():nullptr;
    const auto changedRejected=restartB.ok()?SF::RestoreRestartState(
        restartB.value.get(),changedRestart):
        SEP::Core::Result<std::shared_ptr<const SF::Epoch>>::Failure(
          SEP::Core::StatusCode::NumericalFailure,"restart provider unavailable");
    const std::string resumedJson=resumed.ok()?SF::SerializeEpochJson(*resumed.value):"";
    const std::string resumedCsv=resumed.ok()?SF::SerializeSurfaceCsv(*resumed.value):"";
    Check(restartA.ok()&&checkpointEpoch.ok()&&checkpoint.ok()&&uninterrupted.ok()&&
        restartB.ok()&&restored.ok()&&resumed.ok()&&!changedRejected.ok()&&
        restartB.value->Current()==restartPrior&&
        SF::SerializeEpochJson(*uninterrupted.value)==resumedJson&&
        SF::SerializeSurfaceCsv(*uninterrupted.value)==resumedCsv&&
        resumedJson.find("\"geometric_endpoint_reached\"")!=std::string::npos&&
        resumedCsv.find("downstream_valid,compression")!=std::string::npos,
        "RSH29","background-only restart reproduces later geometry/jumps/exports byte-for-byte and rejects a changed physics identity transactionally");

    // RSH39: samples are values at the moving connected root, so their total
    // derivative includes xdot.grad(F).  A fixed-label partial derivative is
    // a deliberate negative control.  Branch changes and grazing/invalid
    // samples return absence rather than a spurious finite rate.
    auto branchSample=[](double time,std::uint64_t branch,bool valid=true) {
      SF::ConnectedBranchSample sample;sample.branchId=branch;sample.timeS=time;
      sample.positionM={time*time,0,0};sample.obliquityRad=time+2*time*time;
      sample.fastMach=3*time-time*time;sample.valid=valid;return sample;
    };
    const auto branchDerivative=SF::DifferentiateConnectedBranch(
        branchSample(0.9,77),branchSample(1,77),branchSample(1.1,77));
    const auto branchChange=SF::DifferentiateConnectedBranch(
        branchSample(0.9,77),branchSample(1,77),branchSample(1.1,78));
    const auto branchGraze=SF::DifferentiateConnectedBranch(
        branchSample(0.9,77),branchSample(1,77),branchSample(1.1,77,false));
    const double forbiddenFixedLabelDerivative=1;
    Check(branchDerivative.valid&&Near(branchDerivative.intersectionVelocityMPerS.x,
        2,2e-14)&&Near(branchDerivative.obliquityRateRadPerS,5,2e-14)&&
        Near(branchDerivative.fastMachRatePerS,1,2e-14)&&
        std::abs(branchDerivative.obliquityRateRadPerS-
          forbiddenFixedLabelDerivative)>3.9&&!branchChange.valid&&!branchGraze.valid,
        "RSH39","connected-root differentiation includes intersection motion and becomes absent across branch changes/grazes");

    // RSH40: the declared measure includes scientific no-shock and numerical
    // unknown members.  Accepted probability is an interval; deleting unknown
    // members and renormalizing survivors is a detected negative control.
    const auto ensemble=SF::SummarizeEnsemble({
        {"observer-A",0.4,SF::EnsembleOutcome::KnownAcceptedShock},
        {"observer-A",0.35,SF::EnsembleOutcome::KnownNoShock},
        {"observer-A",0.25,SF::EnsembleOutcome::NumericalUnknown}},1e-12);
    const auto badLabel=SF::SummarizeEnsemble({
        {"observer-A",0.5,SF::EnsembleOutcome::KnownAcceptedShock},
        {"observer-B",0.5,SF::EnsembleOutcome::KnownNoShock}},1e-12);
    const auto badWeight=SF::SummarizeEnsemble({
        {"observer-A",0.4,SF::EnsembleOutcome::KnownAcceptedShock},
        {"observer-A",0.4,SF::EnsembleOutcome::KnownNoShock}},1e-12);
    const double forbiddenSurvivorRenormalization=0.4/(0.4+0.35);
    Check(ensemble.ok()&&Near(ensemble.value.acceptedProbabilityLower,0.4,1e-14)&&
        Near(ensemble.value.acceptedProbabilityUpper,0.65,1e-14)&&
        Near(ensemble.value.unknownWeight,0.25,1e-14)&&
        std::abs(forbiddenSurvivorRenormalization-
          ensemble.value.acceptedProbabilityLower)>0.1&&!badLabel.ok()&&!badWeight.ok(),
        "RSH40","ensemble measure retains no-shock/unknown weight and bounds probability without survivor renormalization");
    if(first.ok()) {
      const double s=std::sin(lambda);const double a=smoke->initialApexRadiusM*s/(1+s);
      const double exact=2*3.14159265358979323846*a*a*(1+s);
      Check(Near(first.value->area.geometricSupportM2,exact,2e-15),"RSH23",
          "curved SSE cap area uses exact generating-sphere metric");
    }
    std::vector<double> acceptedFractions,planarAreaErrors;
    bool refinedAreas=true;
    for(int n:{12,24,48}) {
      auto refinedConfiguration=std::make_shared<SF::Configuration>(*longEvent);
      refinedConfiguration->polarCells=n;refinedConfiguration->azimuthCells=2*n;
      const auto refinedProvider=SF::Provider::Create(refinedConfiguration);
      const auto refined=refinedProvider.ok()?refinedProvider.value->Prepare(14000,1):
          SEP::Core::Result<std::shared_ptr<const SF::Epoch>>::Failure(
            SEP::Core::StatusCode::NumericalFailure,"refined provider unavailable");
      if(!refined.ok()){refinedAreas=false;continue;}
      long double recordArea=0,planarArea=0;
      for(const auto& record:refined.value->records)
        recordArea+=record.geometry.areaM2;
      for(const auto& triangle:refined.value->triangles)
        planarArea+=triangle.planarAreaM2;
      refinedAreas=refinedAreas&&Near(static_cast<double>(recordArea),
          refined.value->area.geometricSupportM2,3e-15)&&
          refined.value->area.acceptedShockM2<=refined.value->area.superfastCandidateM2&&
          refined.value->area.superfastCandidateM2<=refined.value->area.geometricSupportM2&&
          refined.value->area.numericalFailureM2<=refined.value->area.geometricSupportM2;
      acceptedFractions.push_back(refined.value->area.acceptedShockM2/
          refined.value->area.geometricSupportM2);
      planarAreaErrors.push_back(std::abs(static_cast<double>(planarArea)-
          refined.value->area.geometricSupportM2)/
          refined.value->area.geometricSupportM2);
    }
    std::vector<SF::ShockRecord> disconnected(6);
    for(std::size_t i=0;i<disconnected.size();++i) {
      disconnected[i].geometry.stableId=i+1;disconnected[i].geometry.areaM2=1;
      disconnected[i].status=SF::FrontStatus::SubfastFront;
    }
    disconnected[0].status=disconnected[3].status=SF::FrontStatus::SolvedFastShock;
    disconnected[0].fastMach=disconnected[3].fastMach=2;
    disconnected[2].status=SF::FrontStatus::InvalidJump;disconnected[2].fastMach=1.2;
    disconnected[5].geometry.areaM2=0; // empty/degenerate chart cell
    const auto disconnectedArea=SF::SummarizeAreas(disconnected);
    const double coarseChange=acceptedFractions.size()==3?
        std::abs(acceptedFractions[1]-acceptedFractions[0]):1;
    const double fineChange=acceptedFractions.size()==3?
        std::abs(acceptedFractions[2]-acceptedFractions[1]):1;
    std::ostringstream refinementDetail;
    refinementDetail<<" refined="<<refinedAreas<<" fractions=";
    for(double value:acceptedFractions)refinementDetail<<value<<',';
    refinementDetail<<" changes="<<coarseChange<<','<<fineChange;
    const bool chordConverges=planarAreaErrors.size()==3&&
        planarAreaErrors[0]/planarAreaErrors[1]>3.5&&
        planarAreaErrors[1]/planarAreaErrors[2]>3.5;
    Check(refinedAreas&&acceptedFractions.size()==3&&acceptedFractions[2]>0&&
        acceptedFractions[2]<1&&fineChange<=coarseChange+1e-14&&chordConverges&&
        disconnectedArea.geometricSupportM2==5&&
        disconnectedArea.superfastCandidateM2==3&&
        disconnectedArea.acceptedShockM2==2&&
        disconnectedArea.numericalFailureM2==1,"RSH23",
        "actual curved area/accepted footprint refines at three angular levels; disconnected, unknown and zero-area accounting stays absolute"+
        refinementDetail.str());

    const auto finalEpoch=longProvider->Prepare(endpoint.value,100);
    bool anyNonForward=false,anyFalseDownstream=false;
    if(finalEpoch.ok())for(const auto& record:finalEpoch.value->records) {
      anyNonForward|=record.status==SF::FrontStatus::NonForwardInflow;
      anyFalseDownstream|=!record.downstreamValid&&
          record.status!=SF::FrontStatus::SolvedFastShock;
    }
    Check(finalEpoch.ok()&&finalEpoch.value->geometricEndpointReached&&anyNonForward,
        "RSH38","1-AU geometric arrival remains distinct from non-forward shock outcome");
    Check(finalEpoch.ok()&&anyFalseDownstream,"RSH19",
        "non-shock records carry absent rather than fabricated downstream state");
    auto requiredConfiguration=std::make_shared<SF::Configuration>(*longEvent);
    requiredConfiguration->validityPolicy=
        SF::TrajectoryValidityPolicy::RequireDeclaredShockCoverage;
    requiredConfiguration->requireFastShockAtObserver=true;
    const auto requiredProvider=SF::Provider::Create(requiredConfiguration);
    const auto requiredPrior=requiredProvider.ok()?requiredProvider.value->Prepare(600,1):
        SEP::Core::Result<std::shared_ptr<const SF::Epoch>>::Failure(
          SEP::Core::StatusCode::NumericalFailure,"required provider unavailable");
    const auto preserved=requiredProvider.ok()?requiredProvider.value->Current():nullptr;
    const auto requiredEndpoint=requiredProvider.ok()?requiredProvider.value->Prepare(
        endpoint.value,2):SEP::Core::Result<std::shared_ptr<const SF::Epoch>>::Failure(
          SEP::Core::StatusCode::NumericalFailure,"required provider unavailable");
    const auto unchangedTrajectory=requiredProvider.ok()?requiredProvider.value->Trajectory(
        endpoint.value):SEP::Core::Result<SF::TrajectoryState>::Failure(
          SEP::Core::StatusCode::NumericalFailure,"required provider unavailable");
    Check(longHandoff.ok()&&longHandoff.value.canonicalWindValid&&
        std::abs(longHandoff.value.effectiveMinusCanonicalWindMPerS)>1&&
        requiredProvider.ok()&&requiredPrior.ok()&&!requiredEndpoint.ok()&&
        requiredProvider.value->Current()==preserved&&
        requiredProvider.value->Current()->generation==1&&unchangedTrajectory.ok()&&
        Near(unchangedTrajectory.value.apexRadiusM,longEvent->endpointRadiusM,0,2),
        "RSH38","proxy/canonical wind discrepancy is retained and required shock coverage rejects transactionally without clipping the trajectory");
  } catch(const std::exception& error) {
    ++counts.error;std::cerr<<"ERROR RSH-HARNESS "<<error.what()<<'\n';
  }
  std::cout<<"SUMMARY PASS="<<counts.pass<<" FAIL="<<counts.fail
           <<" SKIP="<<counts.skip<<" ERROR="<<counts.error<<'\n';
  return counts.fail||counts.error?1:0;
}

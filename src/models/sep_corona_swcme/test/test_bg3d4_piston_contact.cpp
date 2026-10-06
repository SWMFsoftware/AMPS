#include "sep_corona_swcme/ambient_state.h"
#include "sep_corona_swcme/piston_contact.h"
#include "sep_coronal_cme/configuration_parser.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using SEP::CoronalCME::Dot;
using SEP::CoronalCME::Norm;
using SEP::CoronalCME::Unit;
using SEP::CoronalCME::Vec3;
using SEP::CoronaSwcme::EventConfiguration;
using SEP::CoronaSwcme::PistonContactModel;
using SEP::CoronaSwcme::PistonContactPhase;
using SEP::CoronaSwcme::PistonRayDisposition;
using SEP::Core::Result;

void Require(bool condition,const std::string& message) {
  if(!condition)throw std::runtime_error(message);
}

bool Close(double a,double b,double relative,double absolute=0) {
  return std::abs(a-b)<=absolute+relative*std::max(std::abs(a),std::abs(b));
}

std::string FileBytes(const std::string& path) {
  std::ifstream input(path,std::ios::binary);
  Require(static_cast<bool>(input),"cannot read fixture: "+path);
  std::ostringstream bytes;bytes<<input.rdbuf();
  Require(input.good()||input.eof(),"cannot finish reading fixture: "+path);
  return bytes.str();
}

Result<std::string> ReadFile(const std::string& path) {
  try {
    return Result<std::string>::Success(FileBytes(path));
  } catch(const std::exception& error) {
    return Result<std::string>::Failure(
        SEP::Core::StatusCode::DataIntegrityFailure,error.what());
  }
}

std::shared_ptr<const EventConfiguration> Event() {
  const auto event=SEP::CoronaSwcme::ResolveEventConfiguration(
      FileBytes("examples/bg3d4_piston/event.conf"),ReadFile);
  Require(event.ok(),"piston event rejected: "+event.status.message);
  return event.value;
}

std::shared_ptr<const PistonContactModel> Contact(
    const std::shared_ptr<const EventConfiguration>& event) {
  const auto model=PistonContactModel::Create(event);
  Require(model.ok(),"piston contact rejected: "+model.status.message);
  return model.value;
}

void TestIndependentIdentityAndAuthority() {
  const auto event=Event();
  Require(event->assets.size()==5&&event->pistonContact.enabled&&
      event->pistonRays.enabled&&event->pistonNumerics.enabled,
      "piston contact/rays/numerics are absent from immutable event identity");
  Require(std::any_of(event->assets.begin(),event->assets.end(),[](const auto& a) {
    return a.role=="piston_contact"&&a.sha256.size()==64;
  }),"piston contact role/checksum is missing");
  const auto model=Contact(event);
  const auto shock=event->At(event->support.startS);
  const auto piston=model->At(event->support.startS);
  Require(shock.ok()&&piston.ok(),"initial authorities cannot be queried");
  Require(std::abs(shock.value.apexRadiusM.value-
      piston.value.apexRadiusM.value)>1e8,
      "Level-A shock history was silently reused as the Level-B piston");

  std::string contact=FileBytes("examples/bg3d4_piston/contact.asset");
  const std::string oldRate="center_rate_m_s=700000";
  contact.replace(contact.find(oldRate),oldRate.size(),
      "center_rate_m_s=700001");
  const std::string digest=SEP::CoronalCME::ComputeContentChecksum(contact);
  std::string input=FileBytes("examples/bg3d4_piston/event.conf");
  const std::string oldDigest=event->assets.back().sha256;
  input.replace(input.find(oldDigest),oldDigest.size(),digest);
  const auto changed=SEP::CoronaSwcme::ResolveEventConfiguration(input,
      [&](const std::string& path) {
        if(path=="examples/bg3d4_piston/contact.asset")
          return Result<std::string>::Success(contact);
        return ReadFile(path);
      });
  Require(changed.ok(),"valid contact mutation was not resolvable");
  Require(changed.value->physicsFingerprint!=event->physicsFingerprint,
      "contact physics mutation did not change event fingerprint");

  std::string numericalInput=FileBytes("examples/bg3d4_piston/event.conf");
  const std::string oldCfl="bg3d4.cfl=0.2",newCfl="bg3d4.cfl=0.21";
  numericalInput.replace(numericalInput.find(oldCfl),oldCfl.size(),newCfl);
  const auto numericalChanged=SEP::CoronaSwcme::ResolveEventConfiguration(
      numericalInput,ReadFile);
  Require(numericalChanged.ok()&&
      numericalChanged.value->physicsFingerprint!=event->physicsFingerprint,
      "piston numerical mutation did not change event fingerprint");

  const auto stale=SEP::CoronaSwcme::ResolveEventConfiguration(
      FileBytes("examples/bg3d4_piston/event.conf"),
      [&](const std::string& path) {
        if(path=="examples/bg3d4_piston/contact.asset")
          return Result<std::string>::Success(contact);
        return ReadFile(path);
      });
  Require(!stale.ok()&&
      stale.status.code==SEP::Core::StatusCode::DataIntegrityFailure,
      "stale contact checksum was not rejected transactionally");
}

void TestRampHandoffAndCoverage() {
  const auto event=Event();
  const auto model=Contact(event);
  const auto& c=event->pistonContact;
  const double t0=c.startS,t1=t0+c.startupRampDurationS;
  const auto initial=model->At(t0),middle=model->At(t0+0.5*c.startupRampDurationS),
      end=model->At(t1),late=model->At(event->support.endS);
  Require(initial.ok()&&middle.ok()&&end.ok()&&late.ok(),
      "contact time law left support");
  Require(initial.value.phase==PistonContactPhase::StartupRamp&&
      initial.value.apexRadiusM.firstDerivative==0&&
      initial.value.apexRadiusM.secondDerivative==0,
      "startup does not begin from a C2-compatible rest state");
  const double prescribed=c.centerRateMPerS+c.radialRateMPerS;
  Require(Close(middle.value.apexRadiusM.secondDerivative,
      1.875*prescribed/c.startupRampDurationS,2e-14),
      "quintic startup does not attain its analytic acceleration maximum");
  Require(Close(end.value.apexRadiusM.firstDerivative,prescribed,1e-14)&&
      std::abs(end.value.apexRadiusM.secondDerivative)<1e-12,
      "startup ramp does not join the prescribed rate with C2 continuity");
  Require(late.value.apexRadiusM.value>149597870700.0&&
      late.value.apexRadiusM.value<event->support.coverageRadiusM,
      "generic contact lacks actual 1-AU propagation inside coverage");

  const double th=model->HandoffTimeS(),h=0.1;
  const auto before=model->At(th-h),at=model->At(th),after=model->At(th+h);
  Require(before.ok()&&at.ok()&&after.ok()&&
      at.value.phase==PistonContactPhase::CoronalAnalytic&&
      after.value.phase==PistonContactPhase::HandoffTransition,
      "contact handoff phase ownership is wrong");
  const double derivative=(after.value.apexRadiusM.value-
      before.value.apexRadiusM.value)/(2*h);
  Require(Close(derivative,at.value.apexRadiusM.firstDerivative,2e-9),
      "contact position and velocity are discontinuous at handoff");

  const double outerTime=th+c.handoffTransitionDurationS+1000;
  const auto outer=model->At(outerTime);
  Require(outer.ok()&&outer.value.phase==PistonContactPhase::SwcmeOuter,
      "contact did not enter outer self-similar continuation");
  const auto& reference=at.value.ellipsoid;
  const auto& scaled=outer.value.ellipsoid;
  const double ratios[]={
      scaled.centerDistanceM.value/reference.centerDistanceM.value,
      scaled.radialSemiAxisM.value/reference.radialSemiAxisM.value,
      scaled.firstLateralSemiAxisM.value/reference.firstLateralSemiAxisM.value,
      scaled.secondLateralSemiAxisM.value/reference.secondLateralSemiAxisM.value};
  for(double ratio:ratios)Require(Close(ratio,ratios[0],2e-14),
      "outer contact is not heliocentrically self-similar");
}

void TestRayGeometryAndAnalyticHistory() {
  const auto event=Event();
  const auto model=Contact(event);
  const Vec3 apex=model->Basis().radial;
  const double time=0.5*event->pistonContact.startupRampDurationS;
  const auto state=model->At(time);
  const auto ray=model->EvaluateRay(apex,time);
  Require(state.ok()&&ray.ok()&&ray.value.Supported(),
      "apex ray is not supported");
  Require(Close(ray.value.radiusM,state.value.apexRadiusM.value,2e-15)&&
      Close(ray.value.radialSpeedMPerS,
          state.value.apexRadiusM.firstDerivative,2e-14)&&
      Close(ray.value.radialAccelerationMPerS2,
          state.value.apexRadiusM.secondDerivative,2e-13)&&
      Close(ray.value.incidence,1,2e-14),
      "implicit ray solution disagrees with exact apex kinematics");

  const Vec3 oblique=Unit(apex+0.12*model->Basis().firstLateral+
      0.07*model->Basis().secondLateral);
  const auto center=model->EvaluateRay(oblique,time);
  Require(center.ok()&&center.value.Supported(),
      "chosen independent oblique ray is unsupported");
  const auto local=[&](Vec3 x) {
    return std::vector<double>{Dot(x,model->Basis().radial),
        Dot(x,model->Basis().firstLateral),
        Dot(x,model->Basis().secondLateral)};
  };
  const auto x=local(center.value.positionM);
  const auto& k=state.value.ellipsoid;
  const double implicit=std::pow((x[0]-k.centerDistanceM.value)/
      k.radialSemiAxisM.value,2)+
      std::pow(x[1]/k.firstLateralSemiAxisM.value,2)+
      std::pow(x[2]/k.secondLateralSemiAxisM.value,2)-1;
  Require(std::abs(implicit)<2e-14,
      "outermost ray intersection is not on the implicit contact");

  // Fourth-order differences of independently queried positions check the
  // analytic implicit derivative.  The acceleration comparison differentiates
  // position itself rather than subtracting a sheath velocity from a contact
  // velocity produced by the same routine.
  const double h=2;
  const auto m2=model->EvaluateRay(oblique,time-2*h),
      m1=model->EvaluateRay(oblique,time-h),
      p1=model->EvaluateRay(oblique,time+h),
      p2=model->EvaluateRay(oblique,time+2*h);
  Require(m2.ok()&&m1.ok()&&p1.ok()&&p2.ok(),
      "independent contact-history stencil left support");
  const double velocity=(m2.value.radiusM-8*m1.value.radiusM+
      8*p1.value.radiusM-p2.value.radiusM)/(12*h);
  const double acceleration=(-p2.value.radiusM+16*p1.value.radiusM-
      30*center.value.radiusM+16*m1.value.radiusM-m2.value.radiusM)/
      (12*h*h);
  Require(Close(velocity,center.value.radialSpeedMPerS,3e-9)&&
      Close(acceleration,center.value.radialAccelerationMPerS2,2e-5,2e-5),
      "analytic contact derivatives fail independent history differentiation");

  const auto reverse=model->EvaluateRay(-1*apex,time);
  Require(reverse.ok()&&!reverse.value.Supported(),
      "rearward nonintersection was admitted as a piston ray");
  bool sawFlank=false;
  for(int i=1;i<160&&!sawFlank;++i) {
    const double angle=0.5*3.14159265358979323846*i/160;
    const Vec3 q=std::cos(angle)*apex+
        std::sin(angle)*model->Basis().firstLateral;
    const auto candidate=model->EvaluateRay(q,time);
    Require(candidate.ok(),"flank scan failed");
    sawFlank=candidate.value.disposition==PistonRayDisposition::UnsupportedFlank;
  }
  Require(sawFlank,"non-grazing incidence gate did not type an unsupported flank");
}

void TestStartupAmbientCompatibility() {
  const auto event=Event();
  const auto model=Contact(event);
  const auto ambient=SEP::CoronaSwcme::AmbientModel::Create(event);
  Require(ambient.ok(),"ambient model for startup preflight failed");
  std::vector<Vec3> rays;
  for(const auto& ray:event->pistonRays.rays)rays.push_back(ray.direction);
  const auto diagnostics=model->CheckStartupCompatibility(
      *ambient.value,rays);
  Require(diagnostics.ok(),"startup preflight failed: "+
      diagnostics.status.message);
  Require(diagnostics.value.compatible&&diagnostics.value.supportedRays>0&&
      diagnostics.value.supportedRays+diagnostics.value.unsupportedRays==
          rays.size()&&
      diagnostics.value.incompatibleRays==0&&
      diagnostics.value.maximumAmbientRadialMach<=
          event->pistonContact.startupMachTolerance,
      "contact startup is incompatible with the qualified ambient");
  double supportedSolidAngle=0,unsupportedSolidAngle=0;
  for(const auto& ray:event->pistonRays.rays) {
    const auto state=model->EvaluateRay(ray.direction,
        event->pistonContact.startS);
    Require(state.ok(),state.status.message);
    (state.value.Supported()?supportedSolidAngle:unsupportedSolidAngle)+=
        ray.solidAngleSr;
  }
  Require(supportedSolidAngle>0&&unsupportedSolidAngle>0&&
      Close(supportedSolidAngle+unsupportedSolidAngle,
          4*3.14159265358979323846,2e-15),
      "ray support/unsupported solid-angle budget does not close 4*pi");
  std::cout<<"[CSWC0620-EVIDENCE] fingerprint="<<event->physicsFingerprint
           <<" contact_sha="<<event->assets.back().sha256
           <<" handoff_s="<<model->HandoffTimeS()
           <<" final_apex_m="<<model->At(event->support.endS).value.apexRadiusM.value
           <<" startup_max_mach="<<diagnostics.value.maximumAmbientRadialMach
           <<" supported_rays="<<diagnostics.value.supportedRays
           <<" unsupported_rays="<<diagnostics.value.unsupportedRays
           <<" supported_sr="<<supportedSolidAngle
           <<" unsupported_sr="<<unsupportedSolidAngle<<'\n';
}

} // namespace

int main() {
  try {
    TestIndependentIdentityAndAuthority();
    TestRampHandoffAndCoverage();
    TestRayGeometryAndAnalyticHistory();
    TestStartupAmbientCompatibility();
    std::cout<<"[CSWC0620] PASS independent Level-B piston contact authority\n";
    return 0;
  } catch(const std::exception& error) {
    std::cerr<<"BG3D-4 piston contact FAIL: "<<error.what()<<'\n';
    return 1;
  }
}

#include "sep_corona_swcme/cross_ray_contact.h"
#include "sep_corona_swcme/ambient_state.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

#include <type_traits>

static_assert(!std::is_default_constructible<SEP::CoronaSwcme::CrossRay::CorotatingPistonContact>::value,
    "invalid empty authority must not bypass its validated factory");

namespace {
namespace XR=SEP::CoronaSwcme::CrossRay;
using SEP::CoronalCME::Dot;
using SEP::CoronalCME::Norm;
using SEP::CoronalCME::Vec3;

void Require(bool ok,const std::string& message) {
  if(!ok)throw std::runtime_error(message);
}

std::shared_ptr<const SEP::CoronaSwcme::EventConfiguration> Event() {
  const auto read=[](const std::string& name) {
    std::ifstream f(name,std::ios::binary);std::ostringstream bytes;bytes<<f.rdbuf();
    if(!f)return SEP::Core::Result<std::string>::Failure(
        SEP::Core::StatusCode::DataIntegrityFailure,"cannot read contact fixture: "+name);
    return SEP::Core::Result<std::string>::Success(bytes.str());
  };
  const auto text=read("examples/bg3d4_piston/event.conf");Require(text.ok(),text.status.message);
  const auto event=SEP::CoronaSwcme::ResolveEventConfiguration(text.value,read);
  Require(event.ok(),event.status.message);return event.value;
}

std::shared_ptr<const XR::CorotatingPistonContact> Contact(
    std::shared_ptr<const SEP::CoronaSwcme::EventConfiguration> event,double duration=1200) {
  const auto legacy=SEP::CoronaSwcme::PistonContactModel::Create(event);
  Require(legacy.ok(),legacy.status.message);
  const auto contact=XR::CorotatingPistonContact::Create(legacy.value,duration);
  Require(contact.ok(),contact.status.message);return contact.value;
}

void Startup(const std::shared_ptr<const SEP::CoronaSwcme::EventConfiguration>& event) {
  const auto model=Contact(event);const auto ambient=SEP::CoronaSwcme::AmbientModel::Create(event);
  Require(ambient.ok(),ambient.status.message);
  const double omega=event->ambient.rotationRateRadPerS;
  const auto start=model->RotationAt(0),middle=model->RotationAt(600),end=model->RotationAt(1200),
      late=model->RotationAt(1300);
  Require(start.ok()&&middle.ok()&&end.ok()&&late.ok(),"rotation coverage failed");
  // Independent polynomial endpoint/midpoint values. Integral_0^.5 S5=5/64,
  // so theta(.5T)=27*Omega*T/64, omega(.5T)=Omega/2,
  // alpha(.5T)=-15*Omega/(8T); theta(T)=Omega*T/2.
  Require(start.value.value==0&&start.value.firstDerivative==omega&&
      start.value.secondDerivative==0,"startup is not exact ambient corotation");
  Require(std::abs(middle.value.value-27*omega*1200/64)<1e-16&&
      std::abs(middle.value.firstDerivative-omega/2)<1e-18&&
      std::abs(middle.value.secondDerivative+15*omega/(8*1200))<1e-20,
      "independent quintic rate/angle/acceleration reference failed");
  Require(std::abs(end.value.value-600*omega)<1e-16&&end.value.firstDerivative==0&&
      end.value.secondDerivative==0&&late.value.value==end.value.value&&
      late.value.firstDerivative==0,"corotation endpoint reset or residual spin");
  std::size_t supported=0;double largestFlux=0,largestBeta=0;
  for(const auto& ray:event->pistonRays.rays) {
    const auto old=model->Reference()->EvaluateRay(ray.direction,0);
    const auto now=model->EvaluateRay(ray.direction,0);
    Require(old.ok()&&now.ok(),"startup ray evaluation failed");
    if(!now.value.Supported())continue;
    ++supported;
    Require(Norm(now.value.positionM-old.value.positionM)<1e-6,
        "expanded startup reset the initial contact geometry");
    const auto sample=ambient.value->Evaluate(now.value.positionM,0);
    Require(sample.ok(),sample.status.message);
    const double slip=Dot(sample.value.velocityMPerS,now.value.outwardNormal)-now.value.normalSpeedMPerS;
    largestFlux=std::max(largestFlux,std::abs(sample.value.plasma.massDensityKgM3*slip));
    largestBeta=std::max(largestBeta,std::abs(Dot(sample.value.magneticFieldT,
        now.value.outwardNormal))/Norm(sample.value.magneticFieldT));
  }
  Require(supported==4,"contact startup silently omitted original supported rays");
  Require(largestFlux<1e-20,"corotating contact fails unchanged physical startup mass-flux bound");
  Require(largestBeta>0.3,"actual threaded field was projected into a closed/draped contact");
  std::cout<<"[CSWC0628-XRAY-CONTACT-EVIDENCE] supported="<<supported<<
      " maximum_startup_flux_kg_m2_s="<<largestFlux<<" max_abs_bn_over_b="<<largestBeta<<'\n';
}

void Derivatives(const std::shared_ptr<const SEP::CoronaSwcme::EventConfiguration>& event) {
  const auto model=Contact(event);
  // Independently differentiate actual positions on fixed material labels
  // and actual radii on FIXED HCI rays. These are different derivatives;
  // the latter includes rotation of the ray relative to the reference body.
  for(double t:{100.0,600.0,1190.0,1300.0,
      model->Reference()->HandoffTimeS()+300,event->support.endS-2}) {
    const double h=0.5;
    const auto point=model->SurfacePoint(0.6,0.7,t),v=model->MaterialVelocity(0.6,0.7,t);
    Require(point.ok()&&v.ok(),"material point/velocity evaluation failed");
    const auto p1=model->SurfacePoint(0.6,0.7,t-h),p2=model->SurfacePoint(0.6,0.7,t+h),
        p3=model->SurfacePoint(0.6,0.7,t-2*h),p4=model->SurfacePoint(0.6,0.7,t+2*h);
    Require(p1.ok()&&p2.ok()&&p3.ok()&&p4.ok(),"material difference support missing");
    const Vec3 independent=(8*(p2.value-p1.value)-(p4.value-p3.value))/(12*h);
    Require(Norm(independent-v.value)<0.002,"rotating material velocity is not derivative of position");
    const auto surface=model->Evaluate(point.value,t);Require(surface.ok(),surface.status.message);
    Require(std::abs(surface.value.implicitValue)<2e-13&&
        std::abs(Dot(v.value,surface.value.outwardNormal)-surface.value.normalSpeedMPerS)<1e-8,
        "rotated implicit surface and material normal velocity disagree");
    for(const auto& ray:event->pistonRays.rays) {
      const auto center=model->EvaluateRay(ray.direction,t);
      if(!center.ok()||!center.value.Supported())continue;
      const auto m=model->EvaluateRay(ray.direction,t-h),p=model->EvaluateRay(ray.direction,t+h),
          mm=model->EvaluateRay(ray.direction,t-2*h),pp=model->EvaluateRay(ray.direction,t+2*h);
      Require(m.ok()&&p.ok()&&mm.ok()&&pp.ok()&&m.value.Supported()&&p.value.Supported()&&
          mm.value.Supported()&&pp.value.Supported(),"ray stencil lost supported contact");
      const double speed=(8*(p.value.radiusM-m.value.radiusM)-
          (pp.value.radiusM-mm.value.radiusM))/(12*h);
      const double acceleration=(8*(p.value.radialSpeedMPerS-m.value.radialSpeedMPerS)-
          (pp.value.radialSpeedMPerS-mm.value.radialSpeedMPerS))/(12*h);
      Require(std::abs(speed-center.value.radialSpeedMPerS)<0.002&&
          std::abs(acceleration-center.value.radialAccelerationMPerS2)<2e-5,
          "fixed-HCI ray speed/acceleration omits orientation derivatives");
    }
  }
}

void Interface() {
  const XR::Thermodynamics gas{5.0/3.0,2};
  const XR::FaceGeometry face{{1,0,0},2};
  const XR::Primitive left{3,{2,4,-1},7,{0.5,0.3,-0.2}},
      right{9,{2,4,-1},7,{0.5,0.3,-0.2}};
  const auto good=XR::EvaluateIdealContactResiduals(left,right,face,gas);
  Require(good.ok(),good.status.message);const auto& g=good.value;
  Require(g.leftRelativeMassFluxKgM2S==0&&g.rightRelativeMassFluxKgM2S==0&&
      Norm(g.velocityJumpMPerS)==0&&Norm(g.magneticJumpT)==0&&g.pressureJumpPa==0&&
      Norm(g.tractionJumpPa)==0&&Norm(g.tangentialElectricJumpVPerM)==0&&
      std::abs(g.relativeTotalEnergyFluxJumpWPerM2)<1e-13,
      "density/entropy jump was mistaken for a non-material contact");
  auto changed=right;changed.velocityMPerS.x+=0.1;
  const auto leak=XR::EvaluateIdealContactResiduals(left,changed,face,gas);
  Require(leak.ok()&&std::abs(leak.value.rightRelativeMassFluxKgM2S-0.9)<1e-13,
      "independent signed mass leakage was hidden");
  changed=right;changed.velocityMPerS.y+=0.2;
  const auto shear=XR::EvaluateIdealContactResiduals(left,changed,face,gas);
  Require(shear.ok()&&std::abs(shear.value.velocityJumpMPerS.y-0.2)<1e-14&&
      std::abs(shear.value.tangentialElectricJumpVPerM.z-0.1)<1e-14,
      "threaded contact ignored tangential velocity/electric mismatch");
  // A proper oblique Cartesian rotation must rotate vector residuals and
  // leave scalar balances unchanged. This catches normal/tangent bookkeeping
  // that happens to work only for an axis-aligned contact.
  const auto rotate=[](Vec3 x) {return Vec3{0.6*x.x-0.8*x.y,0.8*x.x+0.6*x.y,x.z};};
  auto rotatedLeft=left,rotatedRight=changed;
  rotatedLeft.velocityMPerS=rotate(left.velocityMPerS);
  rotatedLeft.magneticFieldT=rotate(left.magneticFieldT);
  rotatedRight.velocityMPerS=rotate(changed.velocityMPerS);
  rotatedRight.magneticFieldT=rotate(changed.magneticFieldT);
  const auto oblique=XR::EvaluateIdealContactResiduals(rotatedLeft,rotatedRight,
      {rotate(face.normalLeftToRight),face.normalSpeedMPerS},gas);
  Require(oblique.ok()&&Norm(oblique.value.tangentialElectricJumpVPerM-
      rotate(shear.value.tangentialElectricJumpVPerM))<2e-14&&
      Norm(oblique.value.velocityJumpMPerS-rotate(shear.value.velocityJumpMPerS))<2e-14,
      "oblique contact residuals do not rotate as physical HCI vectors");
  changed=right;changed.magneticFieldT.y+=0.2;
  const auto magnetic=XR::EvaluateIdealContactResiduals(left,changed,face,gas);
  // Bn=.5, mu=2: transverse Maxwell traction change is -.5*.2/2=-.05 Pa.
  Require(magnetic.ok()&&std::abs(magnetic.value.tractionJumpPa.y+0.05)<1e-14,
      "physical Maxwell traction was reduced to scalar pressure");
  changed=right;changed.pressurePa+=0.4;
  const auto pressure=XR::EvaluateIdealContactResiduals(left,changed,face,gas);
  Require(pressure.ok()&&std::abs(pressure.value.tractionJumpPa.x-0.4)<1e-14&&
      std::abs(pressure.value.relativeTotalEnergyFluxJumpWPerM2-0.8)<1e-13,
      "moving-interface pressure work has wrong HCI/ALE sign");
  auto zeroLeft=left,zeroRight=right;zeroLeft.magneticFieldT.x=0;zeroRight.magneticFieldT.x=0;
  const auto zero=XR::EvaluateIdealContactResiduals(zeroLeft,zeroRight,face,gas);
  Require(zero.ok()&&Norm(zero.value.tractionJumpPa)==0,
      "continuous-field branch failed at zero normal B");
}

void Rejection(const std::shared_ptr<const SEP::CoronaSwcme::EventConfiguration>& event) {
  const auto model=Contact(event);const auto before=model->Fingerprint();
  Require(!XR::CorotatingPistonContact::Create(nullptr,1200).ok()&&
      !XR::CorotatingPistonContact::Create(model->Reference(),0).ok()&&
      !XR::CorotatingPistonContact::Create(model->Reference(),-1).ok()&&
      !XR::CorotatingPistonContact::Create(model->Reference(),event->support.endS).ok(),
      "corotation accepted missing authority or invalid transition support");
  Require(!model->RotationAt(-1).ok()&&!model->EvaluateRay({2,0,0},0).ok()&&
      !model->Evaluate({std::numeric_limits<double>::quiet_NaN(),0,0},0).ok(),
      "corotation silently extrapolated/normalized invalid input");
  Require(Contact(event,600)->Fingerprint()!=before,"transition duration missing from contact identity");
  std::size_t excludedWithGeometry=0;
  for(const auto& ray:event->pistonRays.rays) {
    const auto value=model->EvaluateRay(ray.direction,100);
    Require(value.ok(),value.status.message);
    if(!value.value.Supported()&&value.value.radiusM>0) {
      ++excludedWithGeometry;
      Require(Norm(value.value.positionM-value.value.radiusM*ray.direction)<2e-6&&
          std::abs(Dot(ray.direction,value.value.outwardNormal)-value.value.incidence)<2e-13,
          "excluded contact geometry escaped in reference-body rather than HCI coordinates");
    }
  }
  Require(excludedWithGeometry>0,"unsupported-area frame control did not execute");
  const XR::Primitive p{1,{1,0,0},1,{0.5,0.2,0}};
  Require(!XR::EvaluateIdealContactResiduals(p,p,{{2,0,0},1},{5.0/3.0,2}).ok(),
      "contact residuals silently normalized invalid physical normal");
  Require(model->Fingerprint()==before,"failed contact queries changed immutable identity");
}
} // namespace

int main() {
  try {
    const auto event=Event();
    Startup(event);std::cout<<"[CSWC0628-XRAY-CONTACT-STARTUP] PASS full-vector ambient corotation compatibility\n";
    Derivatives(event);std::cout<<"[CSWC0628-XRAY-CONTACT-DERIVATIVES] PASS independent moving-point and fixed-HCI-ray derivatives\n";
    Interface();std::cout<<"[CSWC0628-XRAY-CONTACT-INTERFACE] PASS manufactured threaded-contact dimensional residuals\n";
    Rejection(event);std::cout<<"[CSWC0628-XRAY-CONTACT-REJECTION] PASS explicit law identity and typed rejection\n";
  }catch(const std::exception& e){std::cerr<<"cross-ray contact FAIL: "<<e.what()<<'\n';return 1;}
  return 0;
}

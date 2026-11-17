#include "sep_analytic_radial_validation.h"
#include "sep_analytic_radial_background.h"
#include "sep_background_derivatives.h"
#include "sep_initialization.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_background_snapshot.h)

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace SEP { namespace Testing { namespace {
namespace R=Background::AnalyticRadial;
struct Checks {
  unsigned failures=0;std::string message;
  void Require(bool condition,const std::string& reason) {
    if(condition)return;
    ++failures;if(!message.empty())message+="; ";message+=reason;
  }
  bool Status(const Transport::Status& status) {Require(status.ok(),status.message);return status.ok();}
};
Result Finish(const Checks& check,const std::string& message) {
  Result result;result.status=check.failures?Status::Fail:Status::Pass;
  result.message=check.failures?check.message:message;
  result.metrics.push_back({"assertion_failures",static_cast<double>(check.failures),0,"<=","count"});
  result.configuration.push_back("scope=explicit-fixed-radial-velocity-and-geometry-values;no-native-startup-or-thermal-closure-claim");
  return result;
}
R::Wind Declaration() {
  R::Wind wind;wind.originM={{2.0e10,-1.0e10,3.0e9}};wind.speedMPerS=725000;
  return wind;
}
Result TranslatedVelocity() {
  Checks check;const auto wind=Declaration();
  const double position[3]={5.0e10,3.0e10,3.0e9};R::Sample sample;
  if(check.Status(R::Evaluate(wind,position,&sample))) {
    // The independent 3-4-5 reference does not use the helper's radius/unit
    // direction. Translating both origin and point leaves those values fixed.
    check.Require(sample.radiusM==5.0e10,"translated physical radius changed");
    check.Require(std::fabs(sample.velocityMPerS[0]-435000)<=1e-9&&
                  std::fabs(sample.velocityMPerS[1]-580000)<=1e-9&&sample.velocityMPerS[2]==0,
                  "explicit nondefault wind lost the translated 3-4-5 direction");
    R::Wind translated=wind;translated.originM={{-4.0e10,7.0e10,-9.0e9}};
    const double translatedPoint[3]={-1.0e10,11.0e10,-9.0e9};R::Sample second;
    if(check.Status(R::Evaluate(translated,translatedPoint,&second)))
      check.Require(second.velocityMPerS==sample.velocityMPerS&&second.gradientPerS==sample.gradientPerS,
                    "joint origin/point translation changed the physical radial field");
  }
  return Finish(check,"translated-origin nondefault-speed radial U matches independent 3-4-5 values");
}
Result TensorReferences() {
  Checks check;const auto wind=Declaration();const double position[3]={5e10,3e10,3e9};R::Sample sample;
  if(check.Status(R::Evaluate(wind,position,&sample))) {
    const double scale=725000.0/5e10;
    const double factors[9]={0.64,-0.48,0,-0.48,0.36,0,0,0,1};
    for(unsigned index=0;index<9;++index)
      check.Require(std::fabs(sample.gradientPerS[index]/scale-factors[index])<=2e-15,
                    "radial Cartesian tensor differs from the closed-form projection matrix");
    const double radial[3]={0.6,0.8,0},transverse[3]={0,0,1};
    Transport::VelocityDerivatives radialResult,transverseResult;
    if(check.Status(Transport::ContractVelocityGradient(sample.gradientPerS,radial,&radialResult))&&
        check.Status(Transport::ContractVelocityGradient(sample.gradientPerS,transverse,&transverseResult))) {
      check.Require(std::fabs(radialResult.divergencePerS/scale-2)<=2e-15,
                    "constant radial speed lost its two transverse expansion dimensions");
      check.Require(std::fabs(radialResult.fieldAlignedStrainPerS/scale)<=2e-15&&
                    std::fabs(transverseResult.fieldAlignedStrainPerS/scale-1)<=2e-15,
                    "field-aligned tensor contraction is not the exact radial/transverse strain");
    }
  }
  return Finish(check,"full radial tensor, divergence and directional strains match independent projection references");
}
Result DifferenceReference() {
  Checks check;const auto wind=Declaration();const double position[3]={5e10,3e10,3e9};
  R::Sample sample;double maximum=0;
  if(check.Status(R::Evaluate(wind,position,&sample)))for(unsigned axis=0;axis<3;++axis) {
    // A central velocity difference is an independent derivative reference.
    // h/r=1e-4 bounds the leading O(h^2/r^2) error; the literal 2e-8 gate also
    // covers binary64 subtraction, without adding an SI max(1,...) floor.
    constexpr double stepM=5.0e6;
    double left[3]={position[0],position[1],position[2]},right[3]={position[0],position[1],position[2]};
    left[axis]-=stepM;right[axis]+=stepM;R::Sample minus,plus;
    if(!check.Status(R::Evaluate(wind,left,&minus))||!check.Status(R::Evaluate(wind,right,&plus)))continue;
    for(unsigned component=0;component<3;++component) {
      const double reference=(plus.velocityMPerS[component]-minus.velocityMPerS[component])/(2*stepM);
      maximum=std::max(maximum,std::fabs(reference-sample.gradientPerS[3*component+axis])/(wind.speedMPerS/sample.radiusM));
    }
  }
  check.Require(maximum<=2.0e-8,"native Cartesian radial gradient failed the central-difference reference");
  auto result=Finish(check,"radial velocity differences independently recover the Cartesian gradient");
  result.metrics.push_back({"central_difference_scaled_error",maximum,2e-8,"<=",""});return result;
}
Result AtomicRejection() {
  Checks check;auto wind=Declaration();const double position[3]={5e10,3e10,3e9};
  R::Sample output;output.radiusM=17;output.velocityMPerS={{19,23,29}};output.gradientPerS.fill(31);
  const auto original=output;
  const auto unchanged=[&](){return output.radiusM==original.radiusM&&output.velocityMPerS==original.velocityMPerS&&output.gradientPerS==original.gradientPerS;};
  const double origin[3]={wind.originM[0],wind.originM[1],wind.originM[2]};
  check.Require(!R::Evaluate(wind,origin,&output).ok()&&unchanged(),"undefined-origin rejection modified a caller sample");
  wind.speedMPerS=0;
  check.Require(!R::Evaluate(wind,position,&output).ok()&&unchanged(),"invalid-speed rejection modified a caller sample");
  wind=Declaration();wind.originM[1]=std::numeric_limits<double>::quiet_NaN();
  check.Require(!R::Evaluate(wind,position,&output).ok()&&unchanged(),"nonfinite-origin rejection modified a caller sample");
  wind=Declaration();const double wrong[3]={400000,0,0};
  check.Require(!R::ValidateVelocity(wind,position,wrong,2e-12).ok(),"wrong legacy zero-origin velocity passed configured endpoint verification");
  const double correct[3]={435000,580000,0};
  check.Status(R::ValidateVelocity(wind,position,correct,2e-12));
  return Finish(check,"invalid/undefined radial inputs reject atomically and mismatched endpoint U never becomes a provider fallback");
}
Result GeometryOnlyBuilder() {
  Checks check;Initialization::Configuration geometry;
  geometry.parkerOriginM={2e10,-1e10,3e9};geometry.parkerInitialPointM={5e10,3e10,3e9};
  geometry.parkerLengthM=1e9;geometry.parkerPointCount=7;geometry.solarWindSpeedMPerS=725000;
  geometry.solarRotationRateRadPerS=2.86533e-6;geometry.innerRadiusM=5e10;
  std::vector<Initialization::Vec3> points;
  check.Require(!Initialization::Validate(geometry).ok(),"full runtime configuration unexpectedly accepts absent numerical/source/output contracts");
  if(check.Status(Initialization::BuildParkerLine(geometry,&points))) {
    check.Require(points.size()==7,"geometry-only builder did not preserve requested vertex count");
    check.Require(points.front().x==5e10&&points.front().y==3e10&&points.front().z==3e9,"geometry-only builder changed the initial point");
  }
  for(unsigned defect=0;defect<4;++defect) {
    auto invalid=geometry;
    if(defect==0)invalid.parkerPointCount=1;
    if(defect==1)invalid.solarRotationRateRadPerS=std::numeric_limits<double>::quiet_NaN();
    if(defect==2)invalid.parkerInitialPointM=invalid.parkerOriginM;
    if(defect==3)invalid.innerRadiusM*=2;
    const auto previous=points;
    const auto status=Initialization::BuildParkerLine(invalid,&points);
    bool same=previous.size()==points.size();
    for(std::size_t i=0;same&&i<points.size();++i)
      same=previous[i].x==points[i].x&&previous[i].y==points[i].y&&previous[i].z==points[i].z;
    check.Require(!status.ok()&&same,"invalid geometry changed output or bypassed geometry-domain validation");
  }
  return Finish(check,"geometry-only API needs no manufactured run numerics while complete-configuration and atomic invalid-geometry gates remain strict");
}
Result DeclarationIdentity() {
  Checks check;const auto original=Declaration();
  const auto identity=[](const R::Wind& wind) {
    return Background::FingerprintConfiguration("schema=srcsep-background-v3;analytic-radial="+R::CanonicalDeclaration(wind));
  };
  auto changedOrigin=original;changedOrigin.originM[1]=std::nextafter(original.originM[1],0.0);
  auto changedSpeed=original;changedSpeed.speedMPerS=std::nextafter(original.speedMPerS,std::numeric_limits<double>::infinity());
  check.Require(identity(original)!=identity(changedOrigin),"one-binary64-step origin change lost background identity");
  check.Require(identity(original)!=identity(changedSpeed),"one-binary64-step speed change lost background identity");
  auto positiveZero=original,negativeZero=original;positiveZero.originM[2]=0.0;negativeZero.originM[2]=-0.0;
  check.Require(identity(positiveZero)==identity(negativeZero),"signed-zero presentation changed an identical physical origin identity");
  check.Require(R::CanonicalDeclaration(original).find("7.25000000000000000e+05")!=std::string::npos,
                "canonical physical speed did not retain its declared precision/units");
  return Finish(check,"existing canonical background hash distinguishes adjacent representable origin/speed declarations and preserves signed-zero equivalence");
}
} // namespace
std::vector<Descriptor> AnalyticRadialDescriptors() {
  const char* names[]={"Explicit translated radial velocity","Exact radial Cartesian tensor","Independent velocity-difference gradient","Atomic radial-domain rejection","Geometry-only Parker construction contract","Canonical physical radial declaration identity"};
  const TestCallback callbacks[]={TranslatedVelocity,TensorReferences,DifferenceReference,AtomicRejection,GeometryOnlyBuilder,DeclarationIdentity};
  std::vector<Descriptor> result;
  for(unsigned index=0;index<6;++index) {
    Descriptor descriptor;descriptor.id="ARB0"+std::to_string(index+1);descriptor.name=names[index];descriptor.description=names[index];
    descriptor.group="analytic-radial";descriptor.initialization=InitializationLevel::None;
    descriptor.runtime=RuntimeClass::Routine;descriptor.seedCapability=SeedCapability::Deterministic;descriptor.repeatSupported=true;
    descriptor.seedPolicy="independent deterministic SI tensor/geometry fixtures; no random draws";
    descriptor.supportedBuildModes="strict dependency-light and native callback factory";
    descriptor.stateIsolation="private explicit declarations/configurations/output values; no active initialization install";
    descriptor.callback=callbacks[index];result.push_back(descriptor);
  }
  return result;
}
} }

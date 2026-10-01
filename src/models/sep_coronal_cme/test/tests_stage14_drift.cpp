#include "test_framework.h"
#include "sep_coherent_transport.h"
#include "sep_coronal_cme/vector_math.h"
#include <cmath>
#include <cstring>
#include <limits>

namespace SCCMTest {
namespace {
using namespace SEP::CoronalCME;
using namespace SEP::Coherent;
std::array<double,3> array(Vec3 v){return {v.x,v.y,v.z};}
Vec3 vec(std::array<double,3> a){return {a[0],a[1],a[2]};}
Vec3 parker(Vec3 x,double polarity){
  // The second manufacture has dimensional proton/alpha masses and charges,
  // a 10^7 m length unit and 10^-7 T field unit. The small-domain manufacture
  // below remains an independent dimensionless numerical reference.
  const bool physical=Norm(x)>1e6;
  if(physical)x=x/1e7;
  const double r=Norm(x);const Vec3 radial=x/r;
  // Divergence-free static Parker field away from the Sun and polar axis.
  return polarity*(physical?1e-7:1.0)*(100/(r*r))*(radial+Vec3{0.01*x.y,-0.01*x.x,0});
}
SmoothSnapshot snapshot(Vec3 x,double polarity){
  SmoothSnapshot s;s.frameId="manufactured-inertial";s.generation=1;
  const double scale=Norm(x)>1e6?1e7:1.;
  s.magneticFieldT=array(parker(x,polarity));s.gradientScaleM=scale;s.curvatureScaleM=scale;s.stepS=0;s.maximumOrderingError=0.01;
  const double h=1e-4*scale;double db[3][3]{};
  for(int j=0;j<3;j++){
    Vec3 step{};if(j==0)step.x=h;if(j==1)step.y=h;if(j==2)step.z=h;
    const Vec3 a=parker(x+step,polarity),b=parker(x-step,polarity);
    s.gradientMagnitudeTPerM[j]=(Norm(a)-Norm(b))/(2*h);
    const Vec3 derivative=(Unit(a)-Unit(b))/(2*h);
    db[j][0]=derivative.x;db[j][1]=derivative.y;db[j][2]=derivative.z;
  }
  s.curlUnitFieldPerM={db[1][2]-db[2][1],db[2][0]-db[0][2],db[0][1]-db[1][0]};return s;
}
struct Orbit { Vec3 x,p; };
Orbit add(Orbit a,Orbit b,double dt){return {a.x+dt*b.x,a.p+dt*b.p};}
Orbit rate(Orbit a,double mass,double charge,double polarity){
  constexpr double c=299792458.0;
  const Vec3 v=a.p/(mass*std::hypot(1.0,Norm(a.p)/(mass*c)));
  return {v,charge*Cross(v,parker(a.x,polarity))};
}
Orbit orbitStep(Orbit a,double mass,double q,double polarity,double dt){
  const auto k1=rate(a,mass,q,polarity),k2=rate(add(a,k1,dt/2),mass,q,polarity),
      k3=rate(add(a,k2,dt/2),mass,q,polarity),k4=rate(add(a,k3,dt),mass,q,polarity);
  return {a.x+dt/6*(k1.x+2*k2.x+2*k3.x+k4.x),a.p+dt/6*(k1.p+2*k2.p+2*k3.p+k4.p)};
}
Phase phaseAdd(Phase a,const FocusedCharacteristic& b,double dt){
  a.positionM=array(vec(a.positionM)+dt*vec(b.positionRateMPerS));a.momentumSi+=dt*b.momentumRateSiPerS;
  a.pitchCosine+=dt*b.pitchRatePerS;a.epochS+=dt;return a;
}
Phase focusedStep(Phase a,double polarity,double dt){
  auto eval=[&](const Phase& p){auto result=EvaluateFocused(snapshot(vec(p.positionM),polarity),p);Require(result.ok(),result.status.message);return result.value;};
  auto k1=eval(a),k2=eval(phaseAdd(a,k1,dt/2)),k3=eval(phaseAdd(a,k2,dt/2)),k4=eval(phaseAdd(a,k3,dt));
  FocusedCharacteristic sum;sum.positionRateMPerS=array((vec(k1.positionRateMPerS)+2*vec(k2.positionRateMPerS)+2*vec(k3.positionRateMPerS)+vec(k4.positionRateMPerS))/6);
  sum.momentumRateSiPerS=(k1.momentumRateSiPerS+2*k2.momentumRateSiPerS+2*k3.momentumRateSiPerS+k4.momentumRateSiPerS)/6;
  sum.pitchRatePerS=(k1.pitchRatePerS+2*k2.pitchRatePerS+2*k3.pitchRatePerS+k4.pitchRatePerS)/6;
  return phaseAdd(a,sum,dt);
}
Phase initial(){Phase p;p.positionM={10,2,3};p.momentumSi=0.001;p.pitchCosine=0.3;p.massKg=1;p.chargeC=1;p.speciesId="manufactured-proton";p.equationFrameId="manufactured-inertial";return p;}
}
void RegisterStage14Drift(Registry* tests){
  (*tests)["FTE3D10"]=[]{
    for(double mass:{1.0,4.0}){
      auto phase=initial();phase.massKg=mass;phase.chargeC=mass==1?1:2;
      const Vec3 b=Unit(parker(vec(phase.positionM),1));const Vec3 tangent=Unit(Cross(b,{0,0,1}));
      const Vec3 momentum=phase.momentumSi*(phase.pitchCosine*b+std::sqrt(1-phase.pitchCosine*phase.pitchCosine)*tangent);
      const Vec3 x=vec(phase.positionM)-Cross(momentum,b)/(phase.chargeC*Norm(parker(vec(phase.positionM),1)));
      Orbit orbit{x,momentum};
      // Independent full-orbit RK4 integrates the Lorentz equation, not the
      // guiding-center equations under test. Charge and mass are distinct.
      for(int j=0;j<20000;j++){orbit=orbitStep(orbit,mass,phase.chargeC,1,0.03);phase=focusedStep(phase,1,0.03);}
      const Vec3 center=orbit.x+Cross(orbit.p,Unit(parker(orbit.x,1)))/(phase.chargeC*Norm(parker(orbit.x,1)));
      Require(Norm(center-vec(phase.positionM))<2e-4,"Parker full-orbit/GC trajectory in magnetization domain");
    }
    for(double massMultiplier:{1.,4.}) {
      auto phase=initial();phase.positionM={1e8,2e7,3e7};phase.massKg=massMultiplier*1.67262192369e-27;
      phase.chargeC=(massMultiplier==1.?1.:2.)*1.602176634e-19;
      phase.momentumSi=massMultiplier*1e-22;phase.speciesId=massMultiplier==1.?"proton":"alpha";
      const Vec3 field=parker(vec(phase.positionM),1),b=Unit(field),tangent=Unit(Cross(b,{0,0,1}));
      const Vec3 momentum=phase.momentumSi*(phase.pitchCosine*b+std::sqrt(1-phase.pitchCosine*phase.pitchCosine)*tangent);
      Orbit orbit{vec(phase.positionM)-Cross(momentum,b)/(phase.chargeC*Norm(field)),momentum};
      for(int j=0;j<20000;j++){orbit=orbitStep(orbit,phase.massKg,phase.chargeC,1,.003);phase=focusedStep(phase,1,.003);}
      const Vec3 center=orbit.x+Cross(orbit.p,Unit(parker(orbit.x,1)))/(phase.chargeC*Norm(parker(orbit.x,1)));
      Require(Norm(center-vec(phase.positionM))<200.,"dimensional proton/alpha full-orbit reference at rho/L < 0.01");
    }
    auto p=initial();p.equation=Equation::ParkerIsotropic;p.pitchApplicable=false;p.pitchCosine=0;
    auto plus=EvaluateParker(snapshot(vec(p.positionM),1),p);p.chargeC=-1;
    auto negative=EvaluateParker(snapshot(vec(p.positionM),1),p);p.chargeC=1;
    auto polarity=EvaluateParker(snapshot(vec(p.positionM),-1),p);
    Require(plus.ok()&&negative.ok()&&polarity.ok(),"signed Parker operator");
    Require(Norm(vec(plus.value.diagnosticDriftVelocityMPerS)+vec(negative.value.diagnosticDriftVelocityMPerS))<1e-14,"charge reversal");
    Require(Norm(vec(plus.value.diagnosticDriftVelocityMPerS)+vec(polarity.value.diagnosticDriftVelocityMPerS))<1e-14,"polarity reversal");
  };
  (*tests)["FTE3D11"]=[]{
    auto p=initial();auto s=snapshot(vec(p.positionM),1);auto characteristic=EvaluateFocused(s,p);
    Require(characteristic.ok()&&characteristic.value.ownership==Ownership::FullFocusedDeterministicCharacteristic,"focused complete ownership");
    Require(characteristic.value.momentumRateSiPerS==0,"static magnetic field cannot invent electric/adiabatic work");
    auto original=characteristic.value;original.momentumRateSiPerS=0.123456789;int calls=0;
    auto baseline=[&]{calls++;return SEP::Core::Result<FocusedCharacteristic>::Success(original);};
    auto disabled=SelectFocused(false,s,p,baseline);Require(disabled.ok()&&calls==1,"disabled original authority selected once");
    Require(std::memcmp(&original.momentumRateSiPerS,&disabled.value.momentumRateSiPerS,sizeof(double))==0&&
        original.positionRateMPerS==disabled.value.positionRateMPerS,"drift-disabled characteristic bitwise recovery");
    for(std::uint32_t region:{Hcs,Separatrix,Null,Transition,Shock}){
      auto invalid=s;invalid.invalidRegion=region;Require(!EvaluateFocused(invalid,p).ok(),"invalid-region interface dispatch");
    }
    auto invalid=s;invalid.magneticFieldT={0,0,0};Require(!EvaluateFocused(invalid,p).ok(),"weak/null field guard");
    invalid=s;invalid.maximumOrderingError=1e-6;Require(!EvaluateFocused(invalid,p).ok(),"gyroradius ordering guard");
    invalid=s;invalid.stationaryFields=false;Require(!EvaluateFocused(invalid,p).ok(),"unqualified time-dependent Hamiltonian rejected");
    auto wrong=p;wrong.equation=Equation::ParkerIsotropic;Require(!EvaluateFocused(s,wrong).ok(),"wrong equation ownership rejected");
    wrong=p;wrong.pitchCosine=std::numeric_limits<double>::quiet_NaN();Require(!EvaluateFocused(s,wrong).ok(),"nonfinite phase cannot pass ordering comparisons");
    auto uniform=s;uniform.magneticFieldT={0,0,1};uniform.electricFieldVPerM={0,0,0.002};
    uniform.gradientMagnitudeTPerM={0,0,0};uniform.curlUnitFieldPerM={0,0,0};
    auto electric=EvaluateFocused(uniform,p);Require(electric.ok(),"stationary electric work domain");
    const double gamma=std::hypot(1.0,p.momentumSi/(p.massKg*299792458.0));
    const double kineticRate=p.momentumSi/(p.massKg*gamma)*electric.value.momentumRateSiPerS;
    const double physicalWork=p.chargeC*0.002*electric.value.positionRateMPerS[2];
    Require(std::abs(kineticRate-physicalWork)<1e-18,"single Hamiltonian electric work equals q E dot velocity");
    Require(std::abs(electric.value.momentumRateSiPerS-p.chargeC*0.002*p.pitchCosine)<1e-14,"uniform-field independent Lorentz momentum reference");
    auto a=p,b=p;for(int j=0;j<100;j++)a=focusedStep(a,1,0.1);for(int j=0;j<200;j++)b=focusedStep(b,1,0.05);
    Require(Norm(vec(a.positionM)-vec(b.positionM))<1e-9&&std::abs(a.pitchCosine-b.pitchCosine)<1e-9,"time convergence");
    const auto moment=[](const Phase& v){return v.momentumSi*v.momentumSi*(1-v.pitchCosine*v.pitchCosine)/(2*v.massKg*Norm(parker(vec(v.positionM),1)));};
    Require(std::abs(moment(a)-moment(p))<1e-12,"unified magnetic-moment invariant");
  };
}
}

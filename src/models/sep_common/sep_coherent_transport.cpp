#include "sep_coherent_transport.h"
#include <algorithm>
#include <cmath>
namespace SEP { namespace Coherent {
namespace {
using V=std::array<double,3>;
V add(V a,V b){for(int i=0;i<3;i++)a[i]+=b[i];return a;}
V mul(double s,V a){for(auto& x:a)x*=s;return a;}
double dot(V a,V b){double sum=0;for(int i=0;i<3;i++)sum+=a[i]*b[i];return sum;}
double norm(V a){return std::sqrt(dot(a,a));}
V cross(V a,V b){return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
bool finite(V a){for(double x:a)if(!std::isfinite(x))return false;return true;}
Core::Status validate(const SmoothSnapshot& s,const Phase& p,Equation e){
  if(p.equation!=e||p.equationFrameId!=s.frameId||p.speciesId.empty()||s.frameId.empty()||!s.generation||
      !(p.massKg>0)||!(p.momentumSi>0)||p.chargeC==0||!std::isfinite(p.chargeC)||
      !std::isfinite(p.massKg)||!std::isfinite(p.momentumSi)||!std::isfinite(p.pitchCosine)||
      !finite(p.positionM)||!finite(s.magneticFieldT)||!finite(s.electricFieldVPerM)||
      !finite(s.gradientMagnitudeTPerM)||!finite(s.curlUnitFieldPerM)||!std::isfinite(p.epochS))
    return Core::Status::Failure(Core::StatusCode::InvalidState,"coherent equation/frame/species mismatch");
  if((e==Equation::FocusedGyrotropic&&(!p.pitchApplicable||std::abs(p.pitchCosine)>1))||
      (e==Equation::ParkerIsotropic&&(p.pitchApplicable||p.pitchCosine!=0)))
    return Core::Status::Failure(Core::StatusCode::InvalidState,"pitch applicability mismatch");
  if(!s.stationaryFields||!s.stationaryInertialEquationFrame)
    return Core::Status::Failure(Core::StatusCode::UnsupportedCapability,"moving-frame Hamiltonian provider not qualified");
  const double b=norm(s.magneticFieldT),scale=std::min(s.gradientScaleM,s.curvatureScaleM);
  if(s.invalidRegion||!(b>0)||!std::isfinite(b)||!(scale>0)||
      !std::isfinite(s.gradientScaleM)||!std::isfinite(s.curvatureScaleM)||
      !(s.maximumOrderingError>0)||!std::isfinite(s.maximumOrderingError)||
      !(s.stepS>=0)||!std::isfinite(s.stepS))
    return Core::Status::Failure(Core::StatusCode::OutOfDomain,"coherent interface/weak-field dispatch");
  const double rho=p.momentumSi/(std::abs(p.chargeC)*b);
  constexpr double c=299792458.0;
  const double speed=p.momentumSi/(p.massKg*std::hypot(1.0,p.momentumSi/(p.massKg*c)));
  if(rho/scale>s.maximumOrderingError||speed*s.stepS/scale>s.maximumOrderingError)
    return Core::Status::Failure(Core::StatusCode::OutOfDomain,"magnetization/step ordering exceeded");
  return Core::Status::Success();
}
}
Core::Result<ParkerOperator> EvaluateParker(const SmoothSnapshot& s,const Phase& p){
  auto valid=validate(s,p,Equation::ParkerIsotropic);
  if(!valid.ok())return Core::Result<ParkerOperator>::Failure(valid.code,valid.message);
  constexpr double c=299792458.0;const double bmag=norm(s.magneticFieldT);
  const V b=mul(1/bmag,s.magneticFieldT);
  const double speed=p.momentumSi/(p.massKg*std::hypot(1.0,p.momentumSi/(p.massKg*c)));
  const double k=p.momentumSi*speed/(3*p.chargeC*bmag);
  ParkerOperator out;
  // K^A_ij=epsilon_ijk k_A b_k; the derived curl is diagnostic only.
  // The application adds this tensor OR its equivalent discrete operator,
  // never a second explicit diagnostic velocity to the same Parker equation.
  out.antisymmetricTensorM2PerS={0,k*b[2],-k*b[1],-k*b[2],0,k*b[0],k*b[1],-k*b[0],0};
  out.diagnosticDriftVelocityMPerS=mul(k,add(s.curlUnitFieldPerM,mul(-1/bmag,cross(s.gradientMagnitudeTPerM,b))));
  out.estimatedOrderingError=p.momentumSi/(std::abs(p.chargeC)*bmag*std::min(s.gradientScaleM,s.curvatureScaleM));
  return Core::Result<ParkerOperator>::Success(out);
}
Core::Result<FocusedCharacteristic> EvaluateFocused(const SmoothSnapshot& s,const Phase& p){
  auto valid=validate(s,p,Equation::FocusedGyrotropic);
  if(!valid.ok())return Core::Result<FocusedCharacteristic>::Failure(valid.code,valid.message);
  constexpr double c=299792458.0;const double bmag=norm(s.magneticFieldT);
  const V b=mul(1/bmag,s.magneticFieldT);const double parallel=p.momentumSi*p.pitchCosine;
  const double gamma=std::hypot(1.0,p.momentumSi/(p.massKg*c));
  const double moment=p.momentumSi*p.momentumSi*(1-p.pitchCosine*p.pitchCosine)/(2*p.massKg*bmag);
  const V starB=add(s.magneticFieldT,mul(parallel/p.chargeC,s.curlUnitFieldPerM));
  const double denominator=dot(b,starB);
  if(!(denominator>0)||std::abs(denominator-bmag)/bmag>s.maximumOrderingError)
    return Core::Result<FocusedCharacteristic>::Failure(Core::StatusCode::OutOfDomain,"B-star parallel ordering/dispatch");
  const V starE=add(s.electricFieldVPerM,mul(-moment/(p.chargeC*gamma),s.gradientMagnitudeTPerM));
  FocusedCharacteristic out;
  // One relativistic Hamiltonian owns streaming, mirror force, grad-B and
  // curvature drift and electric work. H=sqrt(m²c⁴+c²p_parallel²+2m mu B c²).
  // There is no separately appended adiabatic or focusing momentum increment.
  out.positionRateMPerS=mul(1/denominator,add(mul(parallel/(p.massKg*gamma),starB),cross(starE,b)));
  const double parallelRate=p.chargeC*dot(starE,starB)/denominator;
  const double energyRate=p.chargeC*dot(s.electricFieldVPerM,out.positionRateMPerS);
  out.momentumRateSiPerS=p.massKg*gamma*energyRate/p.momentumSi;
  out.pitchRatePerS=(parallelRate-p.pitchCosine*out.momentumRateSiPerS)/p.momentumSi;
  out.estimatedOrderingError=p.momentumSi/(std::abs(p.chargeC)*bmag*std::min(s.gradientScaleM,s.curvatureScaleM));
  return Core::Result<FocusedCharacteristic>::Success(out);
}
Core::Result<FocusedCharacteristic> SelectFocused(bool enabled,const SmoothSnapshot& s,const Phase& p,
    const std::function<Core::Result<FocusedCharacteristic>()>& baseline){
  // The disabled route returns the original authority directly, with no
  // recomputation, frame conversion or additional floating-point operation.
  if(!enabled){
    if(!baseline)return Core::Result<FocusedCharacteristic>::Failure(Core::StatusCode::InvalidState,"missing baseline characteristic authority");
    return baseline();
  }
  return EvaluateFocused(s,p);
}
} }

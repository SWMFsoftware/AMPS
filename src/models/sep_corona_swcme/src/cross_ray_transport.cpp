#include "sep_corona_swcme/cross_ray_transport.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace SEP { namespace CoronaSwcme { namespace CrossRay { namespace {

using CoronalCME::Dot;
using CoronalCME::Norm;
using CoronalCME::Vec3;
using Vector=std::array<double,8>;

bool Finite(Vec3 v) {
  return std::isfinite(v.x)&&std::isfinite(v.y)&&std::isfinite(v.z);
}

bool ValidGas(const Thermodynamics& gas) {
  return std::isfinite(gas.gammaAdiabatic)&&gas.gammaAdiabatic>1&&
      std::isfinite(gas.magneticPermeabilityNPerA2)&&
      gas.magneticPermeabilityNPerA2>0;
}

bool UnitNormal(Vec3 normal) {
  return Finite(normal)&&std::abs(Dot(normal,normal)-1)<=
      64*std::numeric_limits<double>::epsilon();
}

Vector Pack(const ConservedDensity& u) {
  return {u.massDensityKgM3,u.momentumDensityKgM2S.x,
      u.momentumDensityKgM2S.y,u.momentumDensityKgM2S.z,
      u.totalEnergyDensityJPerM3,u.magneticFieldT.x,
      u.magneticFieldT.y,u.magneticFieldT.z};
}

Vector Pack(const Flux& f) {
  return {f.massKgM2S,f.momentumPa.x,f.momentumPa.y,f.momentumPa.z,
      f.energyWPerM2,f.inductionTMPerS.x,f.inductionTMPerS.y,f.inductionTMPerS.z};
}

ConservedDensity Density(const Vector& u) {
  return {u[0],{u[1],u[2],u[3]},u[4],{u[5],u[6],u[7]}};
}

Flux FluxValue(const Vector& f) {
  return {f[0],{f[1],f[2],f[3]},f[4],{f[5],f[6],f[7]}};
}

bool Finite(const Vector& values) {
  for(double value:values)if(!std::isfinite(value))return false;
  return true;
}

Vector AleFlux(const Primitive& p,const ConservedDensity& u,
    const FaceGeometry& face,const Thermodynamics& gas) {
  const Vec3 n=face.normalLeftToRight;
  const double vn=Dot(p.velocityMPerS,n),bn=Dot(p.magneticFieldT,n);
  const double pt=p.pressurePa+Dot(p.magneticFieldT,p.magneticFieldT)/
      (2*gas.magneticPermeabilityNPerA2);
  Flux f;
  f.massKgM2S=p.densityKgM3*vn;
  f.momentumPa=f.massKgM2S*p.velocityMPerS+pt*n-
      (bn/gas.magneticPermeabilityNPerA2)*p.magneticFieldT;
  f.energyWPerM2=(u.totalEnergyDensityJPerM3+pt)*vn-
      Dot(p.magneticFieldT,p.velocityMPerS)*bn/gas.magneticPermeabilityNPerA2;
  f.inductionTMPerS=vn*p.magneticFieldT-bn*p.velocityMPerS;
  auto values=Pack(f);const auto conserved=Pack(u);
  // Reynolds transport at a face moving in the SAME HCI frame. Retaining the
  // full conserved energy here is essential: replacing only vn by vn-wn in
  // the stationary formulas loses pressure/Maxwell work and is not ALE MHD.
  for(std::size_t i=0;i<values.size();++i)
    values[i]-=face.normalSpeedMPerS*conserved[i];
  return values;
}

} // namespace

Core::Result<ConservedDensity> ToConserved(
    const Primitive& state,const Thermodynamics& gas) {
  using Return=Core::Result<ConservedDensity>;
  if(!ValidGas(gas))return Return::Failure(Core::StatusCode::InvalidConfiguration,
      "cross-ray thermodynamics require finite gamma>1 and mu0>0");
  if(!(std::isfinite(state.densityKgM3)&&state.densityKgM3>0&&
      std::isfinite(state.pressurePa)&&state.pressurePa>0&&
      Finite(state.velocityMPerS)&&Finite(state.magneticFieldT)))
    return Return::Failure(Core::StatusCode::InvalidState,
        "cross-ray primitive state is nonfinite or nonpositive");
  ConservedDensity out;
  out.massDensityKgM3=state.densityKgM3;
  out.momentumDensityKgM2S=state.densityKgM3*state.velocityMPerS;
  out.totalEnergyDensityJPerM3=state.pressurePa/(gas.gammaAdiabatic-1)+
      0.5*state.densityKgM3*Dot(state.velocityMPerS,state.velocityMPerS)+
      Dot(state.magneticFieldT,state.magneticFieldT)/(2*gas.magneticPermeabilityNPerA2);
  out.magneticFieldT=state.magneticFieldT;
  if(!Finite(Pack(out))||!(out.totalEnergyDensityJPerM3>0))
    return Return::Failure(Core::StatusCode::NumericalFailure,
      "cross-ray conservative energy/momentum overflowed");
  return Return::Success(out);
}

Core::Result<Primitive> ToPrimitive(
    const ConservedDensity& state,const Thermodynamics& gas) {
  using Return=Core::Result<Primitive>;
  if(!ValidGas(gas))return Return::Failure(Core::StatusCode::InvalidConfiguration,
      "cross-ray thermodynamics require finite gamma>1 and mu0>0");
  if(!Finite(Pack(state))||!(state.massDensityKgM3>0))return Return::Failure(
      Core::StatusCode::InvalidState,"cross-ray conserved density is invalid");
  Primitive out;
  out.densityKgM3=state.massDensityKgM3;
  out.velocityMPerS=state.momentumDensityKgM2S/state.massDensityKgM3;
  out.magneticFieldT=state.magneticFieldT;
  const double thermal=state.totalEnergyDensityJPerM3-
      0.5*Dot(state.momentumDensityKgM2S,state.momentumDensityKgM2S)/
          state.massDensityKgM3-
      Dot(state.magneticFieldT,state.magneticFieldT)/(2*gas.magneticPermeabilityNPerA2);
  out.pressurePa=(gas.gammaAdiabatic-1)*thermal;
  // A hypersonic/subthermal state whose internal energy is lost to rounding
  // is explicitly rejected. No total-energy reset or pressure floor hides it.
  if(!(Finite(out.velocityMPerS)&&std::isfinite(out.pressurePa)&&out.pressurePa>0))
    return Return::Failure(Core::StatusCode::InvalidState,
        "cross-ray conserved state has no resolved positive thermal pressure");
  return Return::Success(out);
}

Core::Result<double> FastSpeedMPerS(const Primitive& state,Vec3 n,
    const Thermodynamics& gas) {
  using Return=Core::Result<double>;
  const auto conserved=ToConserved(state,gas);
  if(!conserved.ok())return Return::Failure(
      conserved.status.code,conserved.status.message);
  if(!UnitNormal(n))return Return::Failure(Core::StatusCode::InvalidConfiguration,
      "cross-ray fast-wave normal must be finite and unit length");
  const double bn=Dot(state.magneticFieldT,n);
  const Vec3 bt=state.magneticFieldT-bn*n;
  const double a2=gas.gammaAdiabatic*state.pressurePa/state.densityKgM3;
  const double denominator=gas.magneticPermeabilityNPerA2*state.densityKgM3;
  const double br2=bn*bn/denominator,bt2=Dot(bt,bt)/denominator;
  const double scale=std::max({a2,br2,bt2});
  if(!(std::isfinite(scale)&&scale>0&&std::isfinite(a2)&&
      std::isfinite(br2)&&std::isfinite(bt2)))return Return::Failure(
          Core::StatusCode::NumericalFailure,"cross-ray canonical wave scale is unresolved");
  // (a^2+v_A^2)^2-4*a^2*v_An^2 = (a^2-v_A^2)^2+4*a^2*v_At^2.
  // Scale before squaring and use explicit tangential B, avoiding negative
  // roundoff discriminants near the parallel acoustic/Alfven degeneracy.
  const double a=a2/scale,b=br2/scale+bt2/scale,t=bt2/scale;
  const double fast2=0.5*scale*(a+b+std::sqrt((a-b)*(a-b)+4*a*t));
  if(!(std::isfinite(fast2)&&fast2>0))return Return::Failure(
      Core::StatusCode::NumericalFailure,"cross-ray fast speed is nonfinite or unresolved");
  return Return::Success(std::sqrt(fast2));
}

Core::Result<FaceFlux> HllFlux(const Primitive& left,const Primitive& right,
    const FaceGeometry& geometry,const Thermodynamics& gas) {
  using Return=Core::Result<FaceFlux>;
  if(!UnitNormal(geometry.normalLeftToRight)||!std::isfinite(geometry.normalSpeedMPerS))
    return Return::Failure(Core::StatusCode::InvalidConfiguration,
        "cross-ray face normal/speed is invalid");
  const auto ul=ToConserved(left,gas),ur=ToConserved(right,gas);
  if(!ul.ok())return Return::Failure(ul.status.code,ul.status.message);
  if(!ur.ok())return Return::Failure(ur.status.code,ur.status.message);
  const Vec3 n=geometry.normalLeftToRight;
  const double bnLeft=Dot(left.magneticFieldT,n),bnRight=Dot(right.magneticFieldT,n);
  const double magneticScale=std::max(Norm(left.magneticFieldT),Norm(right.magneticFieldT));
  if(std::abs(bnLeft-bnRight)>64*std::numeric_limits<double>::epsilon()*magneticScale)
    return Return::Failure(Core::StatusCode::DataIntegrityFailure,
        "cross-ray face states disagree on normal magnetic field; shared face flux required");
  const auto cl=FastSpeedMPerS(left,n,gas),cr=FastSpeedMPerS(right,n,gas);
  if(!cl.ok())return Return::Failure(cl.status.code,cl.status.message);
  if(!cr.ok())return Return::Failure(cr.status.code,cr.status.message);
  const double vl=Dot(left.velocityMPerS,n)-geometry.normalSpeedMPerS;
  const double vr=Dot(right.velocityMPerS,n)-geometry.normalSpeedMPerS;
  FaceFlux out;out.geometry=geometry;
  out.leftWaveMPerS=std::min({0.0,vl-cl.value,vr-cr.value});
  out.rightWaveMPerS=std::max({0.0,vl+cl.value,vr+cr.value});
  out.maximumRelativeSignalMPerS=std::max(-out.leftWaveMPerS,out.rightWaveMPerS);
  const auto fl=AleFlux(left,ul.value,geometry,gas),fr=AleFlux(right,ur.value,geometry,gas);
  Vector f{};
  if(out.leftWaveMPerS>=0)f=fl;
  else if(out.rightWaveMPerS<=0)f=fr;
  else {
    const auto l=Pack(ul.value),r=Pack(ur.value);
    const double sl=out.leftWaveMPerS,sr=out.rightWaveMPerS;
    for(std::size_t i=0;i<f.size();++i)
      f[i]=(sr*fl[i]-sl*fr[i]+sl*sr*(r[i]-l[i]))/(sr-sl);
  }
  if(!Finite(f)||!std::isfinite(out.maximumRelativeSignalMPerS)||
      !(out.maximumRelativeSignalMPerS>0))return Return::Failure(
          Core::StatusCode::NumericalFailure,"cross-ray HLL flux/waves overflowed");
  out.value=FluxValue(f);
  return Return::Success(out);
}

Core::Result<ExchangeCandidate> TransferStationaryFace(const Cell& left,
    const Cell& right,const FaceFlux& flux,double area,double step,
    const Thermodynamics& gas) {
  using Return=Core::Result<ExchangeCandidate>;
  if(!UnitNormal(flux.geometry.normalLeftToRight)||
      !std::isfinite(flux.geometry.normalSpeedMPerS))return Return::Failure(
          Core::StatusCode::InvalidConfiguration,"cross-ray transfer face geometry is invalid");
  if(flux.geometry.normalSpeedMPerS!=0)return Return::Failure(
      Core::StatusCode::UnsupportedCapability,
      "moving-face update requires unimplemented geometric-conservation authority");
  if(!(std::isfinite(left.volumeM3)&&left.volumeM3>0&&
      std::isfinite(right.volumeM3)&&right.volumeM3>0&&
      std::isfinite(area)&&area>0&&std::isfinite(step)&&step>=0))
    return Return::Failure(Core::StatusCode::InvalidConfiguration,
        "cross-ray shared transfer requires positive finite volumes/area and nonnegative step");
  const auto oldLeft=ToPrimitive(left.density,gas),oldRight=ToPrimitive(right.density,gas);
  if(!oldLeft.ok())return Return::Failure(oldLeft.status.code,oldLeft.status.message);
  if(!oldRight.ok())return Return::Failure(oldRight.status.code,oldRight.status.message);
  const auto values=Pack(flux.value);auto l=Pack(left.density),r=Pack(right.density);
  Vector amount{};
  for(std::size_t i=0;i<amount.size();++i) {
    amount[i]=area*step*values[i];
    l[i]-=amount[i]/left.volumeM3;
    r[i]+=amount[i]/right.volumeM3;
  }
  if(!Finite(amount)||!Finite(l)||!Finite(r))return Return::Failure(
      Core::StatusCode::NumericalFailure,"cross-ray integrated transfer overflowed");
  ExchangeCandidate out;
  out.left={left.volumeM3,Density(l)};out.right={right.volumeM3,Density(r)};
  out.transfer={amount[0],{amount[1],amount[2],amount[3]},amount[4],
      {amount[5],amount[6],amount[7]}};
  const auto newLeft=ToPrimitive(out.left.density,gas),
      newRight=ToPrimitive(out.right.density,gas);
  if(!newLeft.ok()||!newRight.ok())return Return::Failure(
      Core::StatusCode::NumericalFailure,
      "cross-ray private transfer loses positive density/thermal pressure; reduce qualified step");
  return Return::Success(std::move(out));
}

} } } // namespace SEP::CoronaSwcme::CrossRay

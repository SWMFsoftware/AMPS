#include "sep_corona_swcme/material_transport.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace SEP { namespace CoronaSwcme { namespace Material { namespace {

bool Finite(CoronalCME::Vec3 x) {
  return std::isfinite(x.x)&&std::isfinite(x.y)&&std::isfinite(x.z);
}

} // namespace

Core::Result<TransportedState> PushForwardAdiabatic(
    const CrossRay::Primitive& reference,const MapJet& map,
    const CrossRay::Thermodynamics& gas) {
  using Return=Core::Result<TransportedState>;
  const auto original=CrossRay::ToConserved(reference,gas);
  if(!original.ok())return Return::Failure(original.status.code,original.status.message);
  if(!Finite(map.positionHciM)||!Finite(map.velocityHciMPerS)||!Finite(map.accelerationHciMPerS2))
    return Return::Failure(Core::StatusCode::InvalidState,"material-map position/time derivatives are nonfinite");
  long double scale=0;
  for(double value:map.deformation) {
    if(!std::isfinite(value))return Return::Failure(Core::StatusCode::InvalidState,
        "material deformation gradient is nonfinite");
    scale=std::max(scale,std::abs(static_cast<long double>(value)));
  }
  for(double value:map.deformationRatePerS)if(!std::isfinite(value))return Return::Failure(
      Core::StatusCode::InvalidState,"material deformation-rate gradient is nonfinite");
  if(!(scale>0))return Return::Failure(Core::StatusCode::NumericalFailure,
      "material map has zero numerical rank");
  long double f[3][3]{},cofactor[3][3]{};
  for(int i=0;i<3;++i)for(int j=0;j<3;++j)f[i][j]=map.deformation[3*i+j]/scale;
  for(int i=0;i<3;++i)for(int j=0;j<3;++j) {
    const int a=(i+1)%3,b=(i+2)%3,c=(j+1)%3,d=(j+2)%3;
    cofactor[i][j]=f[a][c]*f[b][d]-f[a][d]*f[b][c];
  }
  const long double determinant=f[0][0]*cofactor[0][0]+f[0][1]*cofactor[0][1]+f[0][2]*cofactor[0][2];
  long double magnitude=0;
  for(int j=0;j<3;++j) {
    const int c=(j+1)%3,d=(j+2)%3;
    magnitude+=std::abs(f[0][j])*(std::abs(f[1][c]*f[2][d])+std::abs(f[1][d]*f[2][c]));
  }
  const long double numericalBound=64*std::numeric_limits<double>::epsilon()*magnitude;
  if(!std::isfinite(determinant)||!(determinant>numericalBound))return Return::Failure(
      Core::StatusCode::NumericalFailure,
      "material map is folded or determinant unresolved at input precision; no regularization applied");
  const long double jacobian=determinant*scale*scale*scale;
  TransportedState out;out.jacobian=static_cast<double>(jacobian);
  if(!std::isfinite(out.jacobian)||!(out.jacobian>0))return Return::Failure(
      Core::StatusCode::NumericalFailure,"positive material Jacobian is not representable");
  long double gradient[3][3]{},dilatation=0;
  for(int i=0;i<3;++i)for(int j=0;j<3;++j) {
    for(int k=0;k<3;++k)gradient[i][j]+=map.deformationRatePerS[3*i+k]*
        cofactor[j][k]/(determinant*scale);
    out.velocityGradientPerS[3*i+j]=static_cast<double>(gradient[i][j]);
    if(!std::isfinite(out.velocityGradientPerS[3*i+j]))return Return::Failure(
        Core::StatusCode::NumericalFailure,"material velocity gradient overflows SI output");
    if(i==j)dilatation+=gradient[i][j];
  }
  out.primitive.densityKgM3=static_cast<double>(reference.densityKgM3/jacobian);
  // The logarithmic evaluation avoids an unnecessary intermediate overflow
  // in J^gamma; the physical pressure must still be positive/representable.
  out.primitive.pressurePa=static_cast<double>(std::exp(std::log(
      static_cast<long double>(reference.pressurePa))-gas.gammaAdiabatic*std::log(jacobian)));
  out.primitive.velocityMPerS=map.velocityHciMPerS;
  const long double b0[]={reference.magneticFieldT.x,reference.magneticFieldT.y,reference.magneticFieldT.z};
  double magnetic[3]{},rate[3]{};
  for(int i=0;i<3;++i) {
    long double b=0,bd=0;
    for(int j=0;j<3;++j) {
      b+=static_cast<long double>(map.deformation[3*i+j])*b0[j]/jacobian;
      bd+=static_cast<long double>(map.deformationRatePerS[3*i+j])*b0[j]/jacobian;
    }
    magnetic[i]=static_cast<double>(b);rate[i]=static_cast<double>(bd-b*dilatation);
  }
  out.primitive.magneticFieldT={magnetic[0],magnetic[1],magnetic[2]};
  out.magneticMaterialRateTPerS={rate[0],rate[1],rate[2]};
  out.densityMaterialRateKgM3S=static_cast<double>(-out.primitive.densityKgM3*dilatation);
  out.pressureMaterialRatePaPerS=static_cast<double>(-gas.gammaAdiabatic*out.primitive.pressurePa*dilatation);
  const auto conserved=CrossRay::ToConserved(out.primitive,gas);
  if(!conserved.ok())return Return::Failure(conserved.status.code,
      "transported material state rejected: "+conserved.status.message);
  if(!Finite(out.magneticMaterialRateTPerS)||!std::isfinite(out.densityMaterialRateKgM3S)||
      !std::isfinite(out.pressureMaterialRatePaPerS))return Return::Failure(
          Core::StatusCode::NumericalFailure,"material thermodynamic/induction rates overflow SI output");
  return Return::Success(out);
}

} } } // namespace SEP::CoronaSwcme::Material

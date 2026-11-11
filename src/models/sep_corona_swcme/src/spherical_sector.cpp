#include "sep_corona_swcme/spherical_sector.h"

#include "sep_corona_swcme/spherical_shell.h"
#include "sep_coronal_cme/constants.h"

#include <cmath>
#include <limits>

namespace SEP { namespace CoronaSwcme { namespace CrossRay { namespace {

using CoronalCME::Cross;
using CoronalCME::Dot;
using CoronalCME::Vec3;

bool Finite(Vec3 v) {
  return std::isfinite(v.x)&&std::isfinite(v.y)&&std::isfinite(v.z);
}

bool ValidBasis(const SectorBasis& b) {
  constexpr double tolerance=64*std::numeric_limits<double>::epsilon();
  return Finite(b.x)&&Finite(b.y)&&Finite(b.z)&&
      std::abs(Dot(b.x,b.x)-1)<=tolerance&&
      std::abs(Dot(b.y,b.y)-1)<=tolerance&&
      std::abs(Dot(b.z,b.z)-1)<=tolerance&&
      std::abs(Dot(b.x,b.y))<=tolerance&&
      std::abs(Dot(b.x,b.z))<=tolerance&&
      std::abs(Dot(b.y,b.z))<=tolerance&&
      Dot(Cross(b.x,b.y),b.z)>=1-tolerance;
}

Vec3 Hci(const SectorBasis& b,long double x,long double y,long double z) {
  return static_cast<double>(x)*b.x+static_cast<double>(y)*b.y+
      static_cast<double>(z)*b.z;
}

long double DifferenceMinusSine(long double x) {
  // x-sin(x) loses its leading cubic term near a pole's small latitude span.
  // The alternating expansion through x^9 has omitted relative error below
  // 2e-23 for |x|<0.01, beneath both double and 80-bit accumulation roundoff.
  // This is a numerical evaluation of the exact geometric primitive, not a
  // thin-angle physical approximation or a minimum polar area.
  if(std::abs(x)<0.01L) {
    const long double x2=x*x;
    return x*x2*(1.0L/6-x2*(1.0L/120-x2*(1.0L/5040-x2/362880)));
  }
  return x-std::sin(x);
}

} // namespace

Core::Result<SectorMeasure> MeasureSphericalSector(
    const SectorBounds& b,const SectorBasis& basis) {
  using Return=Core::Result<SectorMeasure>;
  const double inputs[]={b.innerRadiusM,b.outerRadiusM,b.muLower,b.muUpper,
      b.phiLowerRad,b.phiUpperRad};
  for(double value:inputs)if(!std::isfinite(value))return Return::Failure(
      Core::StatusCode::InvalidConfiguration,"nonfinite spherical-sector coordinate");
  const long double dphi=static_cast<long double>(b.phiUpperRad)-b.phiLowerRad;
  if(!(b.innerRadiusM>=0&&b.outerRadiusM>b.innerRadiusM&&
      b.muLower>=-1&&b.muUpper<=1&&b.muUpper>b.muLower&&dphi>0&&
      dphi<=2*static_cast<long double>(CoronalCME::Constants::kPi))||!ValidBasis(basis))
    return Return::Failure(Core::StatusCode::InvalidConfiguration,
        "invalid spherical-sector radii/angular bounds or HCI basis");
  const long double lo=b.muLower,hi=b.muUpper;
  const long double lambdaLo=std::asin(lo),lambdaHi=std::asin(hi);
  const long double dlambda=lambdaHi-lambdaLo;
  const long double phiMid=(static_cast<long double>(b.phiLowerRad)+b.phiUpperRad)/2;
  const long double sineIntegral=2*std::cos(phiMid)*std::sin(dphi/2);
  const long double cosineIntegral=2*std::sin(phiMid)*std::sin(dphi/2);
  // Integral sqrt(1-mu^2) dmu = [Delta lambda + cos(lambdaR+lambdaL)
  // sin(Delta lambda)]/2. The alternative below remains well conditioned at
  // polar caps, where the two terms of that expression nearly cancel.
  const long double halfSum=(lambdaLo+lambdaHi)/2;
  const long double transverseIntegral=0.5L*(DifferenceMinusSine(dlambda)+
      2*std::cos(halfSum)*std::cos(halfSum)*std::sin(dlambda));
  const long double radialZ=0.5L*(hi-lo)*(hi+lo)*dphi;
  const long double omega=(hi-lo)*dphi;
  const Vec3 angularVector=Hci(basis,transverseIntegral*sineIntegral,
      transverseIntegral*cosineIntegral,radialZ);
  const long double left=b.innerRadiusM,right=b.outerRadiusM;
  const long double shellAreaFactor=(right-left)*(right+left)/2;
  SectorMeasure out;
  out.solidAngleSr=static_cast<double>(omega);
  out.volumeM3=out.solidAngleSr*SphericalShellVolumePerSrM3(
      b.innerRadiusM,b.outerRadiusM);
  out.volumeCentroidHciM=(SphericalShellCentroidM(b.innerRadiusM,b.outerRadiusM)/
      out.solidAngleSr)*angularVector;
  out.faces[0]={-static_cast<double>(left*left)*angularVector,
      static_cast<double>(left*left*omega)};
  out.faces[1]={static_cast<double>(right*right)*angularVector,
      static_cast<double>(right*right*omega)};
  const auto muFace=[&](long double mu,double sign) {
    const long double transverse=std::sqrt((1-mu)*(1+mu));
    const long double factor=sign*shellAreaFactor;
    return SectorFaceMeasure{
        Hci(basis,-factor*mu*transverse*sineIntegral,
            -factor*mu*transverse*cosineIntegral,
            factor*(1-mu)*(1+mu)*dphi),
        static_cast<double>(shellAreaFactor*transverse*dphi)};
  };
  out.faces[2]=muFace(lo,-1);out.faces[3]=muFace(hi,1);
  const double phiArea=static_cast<double>(shellAreaFactor*dlambda);
  out.faces[4]={Hci(basis,phiArea*std::sin(static_cast<long double>(b.phiLowerRad)),
      -phiArea*std::cos(static_cast<long double>(b.phiLowerRad)),0),phiArea};
  out.faces[5]={Hci(basis,-phiArea*std::sin(static_cast<long double>(b.phiUpperRad)),
      phiArea*std::cos(static_cast<long double>(b.phiUpperRad)),0),phiArea};
  if(!(std::isfinite(out.volumeM3)&&out.volumeM3>0&&
      std::isfinite(out.solidAngleSr)&&out.solidAngleSr>0&&Finite(out.volumeCentroidHciM)))
    return Return::Failure(Core::StatusCode::NumericalFailure,
        "spherical-sector volume/centroid is not representable");
  for(const auto& face:out.faces) {
    if(!(Finite(face.outwardAreaVectorM2)&&std::isfinite(face.physicalAreaM2)&&
        face.physicalAreaM2>=0))return Return::Failure(
            Core::StatusCode::NumericalFailure,"spherical-sector face measure is nonfinite");
    out.areaClosureResidualM2=out.areaClosureResidualM2+face.outwardAreaVectorM2;
  }
  return Return::Success(std::move(out));
}

} } } // namespace SEP::CoronaSwcme::CrossRay

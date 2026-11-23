#include "sep_reduced_source.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace SEP { namespace ReducedParticles {
namespace {
using ReducedShock::Position;
double Dot(const Position& a,const Position& b) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
Position Minus(const Position& a,const Position& b) {return {{a[0]-b[0],a[1]-b[1],a[2]-b[2]}};}
bool Finite(const Position& x){return std::isfinite(x[0])&&std::isfinite(x[1])&&std::isfinite(x[2]);}
}
Core::Result<double> NetProtonRate(const ReducedShock::FrontSample& front,
    double flux,double fraction) {
  using R=Core::Result<double>;
  if(!std::isfinite(flux)||flux<=0||!std::isfinite(fraction)||fraction<0||fraction>1)
    return R::Failure(Core::StatusCode::InvalidConfiguration,"net source requires positive tube magnetic flux [T m^2] and fraction in [0,1]");
  if(!front.acceptedShock)return R::Success(0);
  const auto& up=front.upstream;
  if(!front.upstreamValid||!Finite(up.magneticFieldT)||!Finite(front.outwardNormal)||
      !std::isfinite(up.protonNumberDensityM3)||up.protonNumberDensityM3<=0||
      !std::isfinite(front.inflowMPerS)||front.inflowMPerS<=0)
    return R::Failure(Core::StatusCode::InvalidState,"accepted source front lacks physical upstream proton state");
  const double b=std::sqrt(Dot(up.magneticFieldT,up.magneticFieldT));
  const double bn=std::abs(Dot(up.magneticFieldT,front.outwardNormal));
  const double normalNorm=std::sqrt(Dot(front.outwardNormal,front.outwardNormal));
  if(!std::isfinite(b)||std::abs(normalNorm-1)>64*std::numeric_limits<double>::epsilon())
    return R::Failure(Core::StatusCode::InvalidState,"source magnetic magnitude/unit normal is not finite normalized SI geometry");
  if(front.epochS!=up.epochS||front.generation!=up.generation||front.eventIdentity!=up.eventIdentity)
    return R::Failure(Core::StatusCode::DataIntegrityFailure,"source front and proton state belong to different core realizations");
  if(!(b>0)||!(bn>0))return R::Failure(Core::StatusCode::UnsupportedCapability,
      "tangent/zero-field tube has no finite shock source patch; no obliquity floor is permitted");
  // The tube's conserved flux [Wb] divided by |B dot n| [T] gives the
  // intercepted area on the shock [m^2]. Use upstream proton number
  // density and the shock-frame normal inflow so the product has units
  // protons/s. The configured fraction is an already accelerated yield.
  const double rate=fraction*up.protonNumberDensityM3*front.inflowMPerS*flux/bn;
  if(!std::isfinite(rate))return R::Failure(Core::StatusCode::InvalidState,"net proton source rate is not representable");
  return R::Success(rate);
}
Core::Result<std::vector<ReleasePoint>> ReferenceIntersections(
    const ReducedShock::IntervalView& view,const std::vector<Position>& vertices,
    double epoch,double distance,double tolerance) {
  using R=Core::Result<std::vector<ReleasePoint>>;
  if(vertices.size()<2||!std::isfinite(distance)||distance<=0||!std::isfinite(tolerance)||tolerance<=0)
    return R::Failure(Core::StatusCode::InvalidConfiguration,"reference surface requires a material line, positive L_ref and SI tolerance");
  const auto sphere=view.SphereAt(epoch);if(!sphere.ok())return R::Failure(sphere.status.code,sphere.status.message);
  const double radius=sphere.value.radiusM+distance;
  if(!std::isfinite(radius))return R::Failure(Core::StatusCode::InvalidState,"offset radius is not representable");
  std::vector<ReleasePoint> out;
  for(std::size_t i=0;i+1<vertices.size();++i) {
    if(!Finite(vertices[i])||!Finite(vertices[i+1]))return R::Failure(Core::StatusCode::InvalidState,"nonfinite material chord in reference placement");
    const Position v=Minus(vertices[i+1],vertices[i]),r=Minus(vertices[i],sphere.value.centerM);
    const double a=Dot(v,v),b=Dot(r,v),c=Dot(r,r)-radius*radius;
    if(!(a>0))return R::Failure(Core::StatusCode::InvalidState,"zero material chord in reference placement");
    const double discriminant=b*b-a*c;
    if(discriminant<0)continue;
    // Rationalized quadratic roots avoid cancellation for a surface near a
    // vertex. An exact double root is retained once; an interior shared vertex
    // belongs to the next segment, consistently with floor(S) in native PIC.
    const double q=-b-std::copysign(std::sqrt(discriminant),b);
    double roots[2]={q/a,q==0?-b/a:c/q};
    if(roots[1]<roots[0])std::swap(roots[0],roots[1]);
    for(int j=0;j<2;++j) {
      const double f=roots[j];
      if(f<0||f>1||(j==1&&roots[0]==roots[1])||(f==1&&i+2<vertices.size()))continue;
      ReleasePoint p;p.segment=i;p.fraction=f;
      Position foot{};
      for(int d=0;d<3;++d){p.xM[d]=vertices[i][d]+f*v[d];p.normal[d]=(p.xM[d]-sphere.value.centerM[d])/radius;
        foot[d]=sphere.value.centerM[d]+sphere.value.radiusM*p.normal[d];}
      if(std::abs(std::sqrt(Dot(Minus(p.xM,sphere.value.centerM),Minus(p.xM,sphere.value.centerM)))-radius)>tolerance)
        return R::Failure(Core::StatusCode::InvalidState,"reference intersection exceeds its absolute surface residual tolerance");
      // Query the normal foot through the core support classifier. The rear
      // sphere or an unsupported SSE flank is not an upstream release site.
      const auto front=view.FrontAt(foot,epoch);
      if(!front.ok()) {
        if(front.status.code==Core::StatusCode::OutOfDomain)continue;
        return R::Failure(front.status.code,front.status.message);
      }
      p.normalFoot=front.value;
      if(front.value.supportEdge)return R::Failure(Core::StatusCode::UnsupportedCapability,
          "release touches finite-front support edge; no endpoint extrapolation is permitted");
      out.push_back(p);
    }
  }
  return R::Success(std::move(out));
}
} }

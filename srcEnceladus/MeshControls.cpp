#include "ForwardModel.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace Enceladus { namespace ForwardModel {
namespace {
double Component(const Vec3& p,int d) { return d==0?p.x:d==1?p.y:p.z; }
double Radius2(const Vec3& p) { return p.x*p.x+p.y*p.y+p.z*p.z; }
double MinRadius2(const Vec3& a,const Vec3& b,const Vec3& c={}) {
  double r=0;
  for(int d=0;d<3;++d) {
    const double q=Component(c,d),lo=Component(a,d),hi=Component(b,d);
    const double delta=q<lo?lo-q:q>hi?q-hi:0;
    r+=delta*delta;
  }
  return r;
}
double MaxRadius2(const Vec3& a,const Vec3& b,const Vec3& c={}) {
  double r=0;
  for(int d=0;d<3;++d) {
    const double delta=std::max(std::abs(Component(a,d)-Component(c,d)),
                                std::abs(Component(b,d)-Component(c,d)));
    r+=delta*delta;
  }
  return r;
}
bool Overlap(const Vec3& a,const Vec3& b,const Vec3& c,const Vec3& e) {
  for(int d=0;d<3;++d)
    if(Component(a,d)>Component(e,d)||Component(b,d)<Component(c,d))return false;
  return true;
}
bool Intersects(const RefinementRegion& r,const Vec3& a,const Vec3& b) {
  if(r.geometry==RefinementRegion::Geometry::Box)return Overlap(a,b,r.minimumM,r.maximumM);
  const double minimum=MinRadius2(a,b,r.centerM);
  const double maximum=MaxRadius2(a,b,r.centerM);
  return minimum<=r.outerRadiusM*r.outerRadiusM && maximum>=r.innerRadiusM*r.innerRadiusM;
}
}

double MeshCellSizeAt(const MeshControls& m,const Vec3& p,bool surface) {
  double size=m.backgroundCellSizeM;
  if(surface) {
    const double requested=p.z>0?m.northSurfaceCellSizeM:m.southSurfaceCellSizeM;
    if(requested>0)size=std::min(size,requested);
  }
  if(m.polarCellSizeM>0&&p.z<0&&Radius2(p)>0) {
    const double tangent=std::tan(m.polarHalfAngleDeg*std::acos(-1.0)/180.0);
    if(p.x*p.x+p.y*p.y<=p.z*p.z*tangent*tangent)size=std::min(size,m.polarCellSizeM);
  }
  for(const auto& r:m.regions)if(Intersects(r,p,p))size=std::min(size,r.cellSizeM);
  return size;
}

double MeshCellSizeForBox(const MeshControls& m,const Vec3& a,const Vec3& b,double radius) {
  double size=m.backgroundCellSizeM;
  // A sphere intersects an AABB iff its radius lies between the nearest and
  // farthest point radii. Hemisphere requests are both considered at z=0.
  if(MinRadius2(a,b)<=radius*radius&&MaxRadius2(a,b)>=radius*radius) {
    if(b.z>=0&&m.northSurfaceCellSizeM>0)size=std::min(size,m.northSurfaceCellSizeM);
    if(a.z<=0&&m.southSurfaceCellSizeM>0)size=std::min(size,m.southSurfaceCellSizeM);
  }
  if(a.z<0&&m.polarCellSizeM>0) {
    const double tangent=std::tan(m.polarHalfAngleDeg*std::acos(-1.0)/180.0);
    const Vec3 xyA{a.x,a.y,0},xyB{b.x,b.y,0};
    if(MinRadius2(xyA,xyB)<=a.z*a.z*tangent*tangent)size=std::min(size,m.polarCellSizeM);
  }
  for(const auto& r:m.regions)if(Intersects(r,a,b))size=std::min(size,r.cellSizeM);
  return size;
}

double SegmentBoxDistanceSquared(const Vec3& a,const Vec3& b,const Vec3& lo,const Vec3& hi) {
  // Squared distance from a+t(b-a) to an AABB is continuous convex piecewise
  // quadratic in t. Face-crossing times partition [0,1]; in each interval
  // evaluate its stationary point and endpoints. This is exact geometry, not
  // a sampled approximation that could miss a thin flyby-refinement tube.
  std::vector<double> times{0,1};
  for(int d=0;d<3;++d) {
    const double delta=Component(b,d)-Component(a,d);
    if(delta==0)continue;
    for(double face:{Component(lo,d),Component(hi,d)}) {
      double t=(face-Component(a,d))/delta;
      if(t>0&&t<1)times.push_back(t);
    }
  }
  std::sort(times.begin(),times.end());
  auto distance=[&](double t) {
    Vec3 p{a.x+t*(b.x-a.x),a.y+t*(b.y-a.y),a.z+t*(b.z-a.z)};
    return MinRadius2(lo,hi,p);
  };
  double result=std::numeric_limits<double>::infinity();
  for(double t:times)result=std::min(result,distance(t));
  for(std::size_t i=1;i<times.size();++i) {
    const double middle=(times[i-1]+times[i])/2;
    double quadratic=0,linear=0;
    for(int d=0;d<3;++d) {
      const double start=Component(a,d),delta=Component(b,d)-start,p=start+middle*delta;
      if(p>=Component(lo,d)&&p<=Component(hi,d))continue;
      const double face=p<Component(lo,d)?Component(lo,d):Component(hi,d);
      quadratic+=delta*delta;linear+=delta*(start-face);
    }
    if(quadratic>0)result=std::min(result,distance(std::max(times[i-1],std::min(times[i],-linear/quadratic))));
  }
  return result;
}

Status ValidateMeshCapabilities(const MeshControls& m,double radius,const Vec3& counts,int cap,double extra) {
  if(!m.explicitConfiguration)return Status::Ok();
  double rootMax=0;
  for(int d=0;d<3;++d) {
    const double a=Component(m.minimumM,d),b=Component(m.maximumM,d),n=Component(counts,d);
    if(!std::isfinite(a)||!std::isfinite(b)||!(a<-radius&&b>radius)||!std::isfinite(b-a)||!(n>0))
      return Status::Error("runtime mesh domain must finitely contain the complete Enceladus sphere on every axis");
    rootMax=std::max(rootMax,(b-a)/n);
  }
  double requested=m.backgroundCellSizeM;
  if(!(requested>0)||!std::isfinite(requested))return Status::Error("mesh background cell size must be finite and positive");
  for(double x:{m.northSurfaceCellSizeM,m.southSurfaceCellSizeM,m.polarCellSizeM,extra})
    if(x>0)requested=std::min(requested,x);
  for(const auto& r:m.regions) {
    if(!Intersects(r,m.minimumM,m.maximumM))return Status::Error("refinement '"+r.id+"' does not intersect the runtime domain");
    requested=std::min(requested,r.cellSizeM);
  }
  // Compiled storage/layout still caps AMR depth. Refuse an unattainable
  // explicit request even when the legacy core enforcement option is off.
  if(cap<0||cap>30||requested<std::ldexp(rootMax,-cap))
    return Status::Error("runtime cell-size request exceeds the compiled maximum AMR level; increase the build capability explicitly");
  return Status::Ok();
}
} }

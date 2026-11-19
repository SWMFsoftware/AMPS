#include "TrajectoryExposure.h"
#include "TrajectoryIntervals.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>

namespace Moon { namespace Species { namespace TrajectoryExposure {
namespace {
namespace I=TrajectoryIntervals;
using Polynomial=I::Polynomial;
using IVector=std::array<I::Interval,3>;
// The polynomial describes the ideal continuous arc. The spatial evaluator
// describes the SAME geometric constraint at the represented callback point.
// Keeping both prevents a binary64 position projection from silently crossing
// a silhouette that was excluded by the continuous-arc root proof.
struct Constraint { Polynomial polynomial;std::function<I::Interval(const IVector&)> spatial; };
I::Interval Dot(const IVector&a,const IVector&b){I::Interval result;for(int d=0;d<3;++d)result=I::Add(result,I::Multiply(a[d],b[d]));return result;}
IVector Subtract(const IVector&a,const IVector&b){IVector result;for(int d=0;d<3;++d)result[d]=I::Subtract(a[d],b[d]);return result;}
int StrictSign(const I::Interval&value){if(!I::Finite(value))return 0;return value.lo>0?1:(value.hi<0?-1:0);}
bool Reject(std::string* error,const char* text){if(error)*error=text;return false;}
bool Finite(const Vector& x){return std::all_of(x.begin(),x.end(),[](double a){return std::isfinite(a);});}
IVector Exact(const Vector&x){return {{I::Point(x[0]),I::Point(x[1]),I::Point(x[2])}};}
IVector Sub(const Vector&a,const Vector&b){IVector c;for(int d=0;d<3;++d)c[d]=I::Subtract(I::Point(a[d]),I::Point(b[d]));return c;}
IVector Cross(const IVector&a,const IVector&b){return {{I::Subtract(I::Multiply(a[1],b[2]),I::Multiply(a[2],b[1])),I::Subtract(I::Multiply(a[2],b[0]),I::Multiply(a[0],b[2])),I::Subtract(I::Multiply(a[0],b[1]),I::Multiply(a[1],b[0]))}};}
Vector Position(const Arc&a,long double t){Vector p;for(int d=0;d<3;++d)p[d]=static_cast<double>(a.positionM[d]+t*a.velocityMS[d]+t*t*a.accelerationMS2[d]/2);return p;}
bool SphereVisible(const Vector&x,const Vector&sun,const Sphere&s){
  long double direction[3],offset[3],q=0,projection=0;
  for(int d=0;d<3;++d){direction[d]=static_cast<long double>(sun[d])-x[d];offset[d]=static_cast<long double>(x[d])-s.centreM[d];q+=direction[d]*direction[d];projection+=offset[d]*direction[d];}
  const long double u=std::max(0.0L,std::min(1.0L,-projection/q));long double distance2=0;
  for(int d=0;d<3;++d){const long double z=offset[d]+u*direction[d];distance2+=z*z;}
  return distance2>=static_cast<long double>(s.radiusM)*s.radiusM;
}
Polynomial SolarCoordinate(const Arc&a,const Field&f,int d){return {I::Subtract(I::Point(f.sunPositionM[d]),I::Point(a.positionM[d])),I::Point(-a.velocityMS[d]),I::Multiply(I::Point(-.5L),I::Point(a.accelerationMS2[d]))};}
// Analytic f'' bound, with outward arithmetic from the actual binary inputs.
// Interval coordinate boxes may be conservative; unresolved support refuses
// instead of replacing a positive distance by an arbitrary numerical floor.
bool FactorBounds(const Arc&a,const Field&f,long double lo,long double hi,
    long double* minimum,long double* maximum,long double* curvature){
  I::Interval qmin,qmax,qp,qpp;const I::Interval time{lo,hi};
  for(int d=0;d<3;++d){const auto coordinate=I::Value(SolarCoordinate(a,f,d),time);
    if(!I::Finite(coordinate))return false;
    const long double mx=std::max(std::abs(coordinate.lo),std::abs(coordinate.hi));
    const long double mn=coordinate.lo<=0&&coordinate.hi>=0?0:std::min(std::abs(coordinate.lo),std::abs(coordinate.hi));
    const auto velocity=I::Add(I::Point(a.velocityMS[d]),I::Multiply(I::Point(a.accelerationMS2[d]),time));
    const auto speed=I::Point(std::max(std::abs(velocity.lo),std::abs(velocity.hi)));
    qmin=I::Add(qmin,I::Multiply(I::Point(mn),I::Point(mn)));qmax=I::Add(qmax,I::Multiply(I::Point(mx),I::Point(mx)));
    qp=I::Add(qp,I::Multiply(I::Point(2),I::Multiply(I::Point(mx),speed)));
    qpp=I::Add(qpp,I::Multiply(I::Point(2),I::Add(I::Multiply(speed,speed),I::Multiply(I::Point(mx),I::Point(std::abs(a.accelerationMS2[d]))))));
  }
  if(!(qmin.lo>0)||!I::Finite(qmax))return false;
  const auto ref=I::Multiply(I::Point(f.referenceDistanceM),I::Point(f.referenceDistanceM));
  *minimum=I::Divide(ref,qmax).lo;*maximum=I::Divide(ref,qmin).hi;
  const auto term1=I::Divide(I::Multiply(I::Point(2),I::Multiply(qp,qp)),I::Multiply(qmin,I::Multiply(qmin,qmin)));
  const auto term2=I::Divide(qpp,I::Multiply(qmin,qmin));
  *curvature=I::Multiply(ref,I::Add(term1,term2)).hi;
  return std::isfinite(*minimum)&&std::isfinite(*maximum)&&std::isfinite(*curvature);
}
I::Interval Factor(const Arc&a,const Field&f,long double t){I::Interval q;for(int d=0;d<3;++d){const auto r=I::Value(SolarCoordinate(a,f,d),I::Point(t));q=I::Add(q,I::Multiply(r,r));}return I::Divide(I::Multiply(I::Point(f.referenceDistanceM),I::Point(f.referenceDistanceM)),q);}
long double Midpoint(const I::Interval&i){return i.lo+(i.hi-i.lo)/2;}
long double ProjectionError(const I::Interval&i,long double value){return I::Up(std::max(value-i.lo,i.hi-value));}

}
bool BuildPhotonHistory(const Arc&a,double duration,const Field&f,
    const std::vector<PhotonChannel>& photons,const std::vector<ConstantChannel>& constants,
    const Controls&c,History* output,std::string* error){
  if(!output||!error)return false;
  if(!Finite(a.positionM)||!Finite(a.velocityMS)||!Finite(a.accelerationMS2)||!Finite(f.sunPositionM)||
      f.frame.empty()||!std::isfinite(f.frozenGeometryEpochS)||!(duration>0)||!std::isfinite(duration)||
      !(f.referenceDistanceM>0)||!std::isfinite(f.referenceDistanceM)||!(c.rootTimeToleranceS>0)||
      !std::isfinite(c.rootTimeToleranceS)||c.rootTimeToleranceS>=duration||
      !std::isfinite(c.absoluteDepthTolerance)||c.absoluteDepthTolerance<0||
      !std::isfinite(c.relativeDepthTolerance)||c.relativeDepthTolerance<0||
      !(c.absoluteDepthTolerance>0||c.relativeDepthTolerance>0)||c.maximumSegments==0||c.maximumGeometryIntervals==0)
    return Reject(error,"trajectory exposure requires finite SI arc/frame/epoch and explicit positive numerical/work controls");
  if(f.terrainFaces.empty()!=!static_cast<bool>(f.terrainVisibility))return Reject(error,"terrain silhouette geometry and production ray callback must be supplied together");
  std::set<std::size_t> ids;I::Interval photonSum,constantSum;
  for(const auto&p:photons){if(!std::isfinite(p.rateAtReferenceS)||p.rateAtReferenceS<0||!ids.insert(p.reaction).second)return Reject(error,"invalid/duplicate photon channel");photonSum=I::Add(photonSum,I::Point(p.rateAtReferenceS));}
  for(const auto&p:constants){if(!std::isfinite(p.rateS)||p.rateS<0||!ids.insert(p.reaction).second)return Reject(error,"invalid/duplicate constant channel");constantSum=I::Add(constantSum,I::Point(p.rateS));}
  const long double totalPhotons=photonSum.hi;
  if(ids.empty())return Reject(error,"trajectory exposure requires the complete declared channel family");
  const long double H=duration;const auto durationInterval=I::Point(H);
  std::vector<I::Interval> brackets;std::size_t visits=0;History result;
  std::vector<Constraint> constraints;
  const long double unitTolerance=I::Divide(I::Point(c.rootTimeToleranceS),durationInterval).lo;
  auto roots=[&](const Polynomial& polynomial){return I::Isolate(polynomial,unitTolerance,c.maximumGeometryIntervals,&visits,&brackets,error);};
  auto plane=[&](const IVector&normal,const Vector&point){Polynomial polynomial(3);
    for(int d=0;d<3;++d){polynomial[0]=I::Add(polynomial[0],I::Multiply(normal[d],I::Subtract(I::Point(a.positionM[d]),I::Point(point[d]))));
      polynomial[1]=I::Add(polynomial[1],I::Multiply(normal[d],I::Multiply(I::Point(a.velocityMS[d]),durationInterval)));
      polynomial[2]=I::Add(polynomial[2],I::Multiply(normal[d],I::Multiply(I::Point(.5L*a.accelerationMS2[d]),I::Multiply(durationInterval,durationInterval))));}
    if(!roots(polynomial))return false;
    constraints.push_back({std::move(polynomial),[normal,point](const IVector&x){return Dot(normal,Subtract(x,Exact(point)));}});
    return true;};
  for(const auto&face:f.terrainFaces){for(const auto&v:face.verticesM)if(!Finite(v))return Reject(error,"nonfinite production terrain vertex");
    // Stable edge form, with interval coefficients and nearby plane origin.
    // Root-free complement intervals must be proven even for a distant Sun;
    // a zero/cancelled normal is never silently turned into no shadow.
    for(int edge=0;edge<3;++edge)if(!plane(Cross(Sub(face.verticesM[edge],f.sunPositionM),Sub(face.verticesM[(edge+1)%3],face.verticesM[edge])),face.verticesM[edge]))return false;
    if(!plane(Cross(Sub(face.verticesM[1],face.verticesM[0]),Sub(face.verticesM[2],face.verticesM[0])),face.verticesM[0]))return false;}
  for(const auto&s:f.opaqueSpheres){
    if(!Finite(s.centreM)||!(s.radiusM>0)||!std::isfinite(s.radiusM))return Reject(error,"opaque spheres need finite explicit centres and positive radii");
    Polynomial ray[3],offset[3],norm(1),projection(1),offsetNorm(1);
    for(int d=0;d<3;++d){const auto velocity=I::Multiply(I::Point(a.velocityMS[d]),durationInterval),acceleration=I::Multiply(I::Point(.5L*a.accelerationMS2[d]),I::Multiply(durationInterval,durationInterval));
      ray[d]={I::Subtract(I::Point(f.sunPositionM[d]),I::Point(a.positionM[d])),I::Negate(velocity),I::Negate(acceleration)};
      offset[d]={I::Subtract(I::Point(a.positionM[d]),I::Point(s.centreM[d])),velocity,acceleration};
      I::Accumulate(norm,I::Product(ray[d],ray[d]));I::Accumulate(projection,I::Product(offset[d],ray[d]));I::Accumulate(offsetNorm,I::Product(offset[d],offset[d]));}
    const auto radius2=I::Multiply(I::Point(s.radiusM),I::Point(s.radiusM));offsetNorm[0]=I::Subtract(offsetNorm[0],radius2);
    const auto base=Cross(Sub(a.positionM,s.centreM),Sub(f.sunPositionM,s.centreM));
    const auto velocity=Cross(Exact(a.velocityMS),Sub(f.sunPositionM,s.centreM)),acceleration=Cross(Exact(a.accelerationMS2),Sub(f.sunPositionM,s.centreM));
    Polynomial tangent(1);
    for(int d=0;d<3;++d){Polynomial cross{base[d],I::Multiply(velocity[d],durationInterval),I::Multiply(I::Multiply(I::Point(.5L),acceleration[d]),I::Multiply(durationInterval,durationInterval))};I::Accumulate(tangent,I::Product(cross,cross));}
    I::Accumulate(tangent,norm,I::Negate(radius2));
    if(!roots(tangent)||!roots(projection))return false;
    const IVector sun=Exact(f.sunPositionM),centre=Exact(s.centreM);
    // Identical tangent/projection/endpoint predicates are evaluated outward
    // at the callback's actual represented point, without a position epsilon.
    auto rayAt=[sun](const IVector&x){return Subtract(sun,x);};
    auto offsetAt=[centre](const IVector&x){return Subtract(x,centre);};
    constraints.push_back({tangent,[rayAt,offsetAt,radius2](const IVector&x){const auto ray=rayAt(x),cross=Cross(offsetAt(x),ray);return I::Subtract(Dot(cross,cross),I::Multiply(radius2,Dot(ray,ray)));}});
    constraints.push_back({projection,[rayAt,offsetAt](const IVector&x){return Dot(offsetAt(x),rayAt(x));}});
    I::Accumulate(projection,norm);
    if(!roots(projection)||!roots(offsetNorm))return false;
    constraints.push_back({projection,[rayAt,offsetAt](const IVector&x){const auto ray=rayAt(x);return I::Add(Dot(offsetAt(x),ray),Dot(ray,ray));}});
    constraints.push_back({offsetNorm,[offsetAt,radius2](const IVector&x){const auto offset=offsetAt(x);return I::Subtract(Dot(offset,offset),radius2);}});
  }
  // Union actual unresolved enclosures, not a count-times-epsilon heuristic.
  // Every complementary interval has a certified strict sign for EVERY
  // silhouette polynomial. A short or near-tangent shadow may exist inside
  // a bracket; its entire possible exposure is carried as numerical error.
  std::sort(brackets.begin(),brackets.end(),[](const I::Interval&a,const I::Interval&b){return a.lo<b.lo;});
  std::vector<I::Interval> merged;for(const auto&b:brackets){if(!merged.empty()&&b.lo<=merged.back().hi)merged.back().hi=std::max(merged.back().hi,b.hi);else merged.push_back(b);}
  std::vector<long double> knots{0,H};I::Interval uncertainDuration;
  std::vector<I::Interval> timeBrackets;
  for(const auto&b:merged){const auto lo=I::Multiply(I::Point(b.lo),durationInterval),hi=I::Multiply(I::Point(b.hi),durationInterval);
    const long double begin=std::max(0.0L,lo.lo),end=std::min(H,hi.hi);knots.push_back(begin);knots.push_back(end);
    timeBrackets.push_back({begin,end});uncertainDuration=I::Add(uncertainDuration,I::Subtract(I::Point(end),I::Point(begin)));}
  result.silhouetteRoots=merged.size();result.geometryIntervalsVisited=visits;
  if(knots.size()>c.maximumSegments+1)return Reject(error,"shadow enclosure count exceeds maximumSegments; reduce arc/mesh complexity or increase the explicit bound");
  std::sort(knots.begin(),knots.end());knots.erase(std::unique(knots.begin(),knots.end()),knots.end());
  long double globalMin=0,globalMax=0,globalCurvature=0;
  if(!FactorBounds(a,f,0,H,&globalMin,&globalMax,&globalCurvature))return Reject(error,"Sun-distance support reaches an unresolved/singular interval");
  const long double constantIntegralLower=I::Multiply(constantSum,durationInterval).lo;
  I::Interval visiblePhotonLower;
  struct Work{long double lo,hi;bool visible;};std::vector<Work> pending;
  for(std::size_t i=knots.size()-1;i>0;--i){
    const long double lo=knots[i-1],hi=knots[i];if(!(hi>lo))continue;
    const long double middle=lo+(hi-lo)/2;const Vector x=Position(a,middle);
    if(!Finite(x))return Reject(error,"visibility midpoint position is not representable in callback coordinates");
    // The cells partition the expanded root-bracket union. A cell already in
    // that union has its whole width charged. On every other cell prove that
    // position projection preserves ALL strict geometric predicate signs.
    // If not, charge the WHOLE cell, not a guessed spatial rounding length.
    const bool rootBracket=std::any_of(timeBrackets.begin(),timeBrackets.end(),[&](const I::Interval&b){return lo>=b.lo&&hi<=b.hi;});
    if(!rootBracket){
      const auto unitMiddle=I::Divide(I::Point(middle),durationInterval);const auto represented=Exact(x);bool certified=true;
      for(const auto&constraint:constraints){const int idealSign=StrictSign(I::Value(constraint.polynomial,unitMiddle));const int representedSign=StrictSign(constraint.spatial(represented));
        if(idealSign==0||representedSign!=idealSign){certified=false;break;}}
      if(!certified)uncertainDuration=I::Add(uncertainDuration,I::Subtract(I::Point(hi),I::Point(lo)));
    }
    bool visible=true;for(const auto&s:f.opaqueSpheres)visible=visible&&SphereVisible(x,f.sunPositionM,s);
    if(visible&&f.terrainVisibility&&!f.terrainVisibility(x,f.sunPositionM,&visible,error))return false;
    pending.push_back({lo,hi,visible});
    if(visible)visiblePhotonLower=I::Add(visiblePhotonLower,I::Multiply(I::Subtract(I::Point(hi),I::Point(lo)),I::Multiply(I::Point(globalMin),I::Point(photonSum.lo))));
  }
  result.shadowBoundaryDepthErrorBound=I::Multiply(uncertainDuration,I::Multiply(I::Point(totalPhotons),I::Point(globalMax))).hi;
  const long double lowerIntegral=I::Add(I::Point(constantIntegralLower),I::Point(std::max(0.0L,I::Subtract(visiblePhotonLower,I::Point(result.shadowBoundaryDepthErrorBound)).lo))).lo;
  const long double budget=I::Add(I::Point(c.absoluteDepthTolerance),I::Multiply(I::Point(c.relativeDepthTolerance),I::Point(lowerIntegral))).lo;
  if(result.shadowBoundaryDepthErrorBound>budget)return Reject(error,"shadow-root depth-error envelope exceeds declared exposure budget; reduce rootTimeToleranceS or maximum arc duration");
  const long double smoothBudget=I::Subtract(I::Point(budget),I::Point(result.shadowBoundaryDepthErrorBound)).lo;
  while(!pending.empty()){
    const Work w=pending.back();pending.pop_back();const auto durationBound=I::Subtract(I::Point(w.hi),I::Point(w.lo));const long double h=w.hi-w.lo;long double mn=0,mx=0,curve=0;
    if(!FactorBounds(a,f,w.lo,w.hi,&mn,&mx,&curve))return Reject(error,"unresolved solar-distance subinterval");
    const auto begin=w.visible?Factor(a,f,w.lo):I::Point(0),end=w.visible?Factor(a,f,w.hi):I::Point(0);
    std::vector<PiecewiseHazard::Channel> rates;I::Interval projectionBound;
    for(const auto&p:photons){const auto first=I::Multiply(I::Point(p.rateAtReferenceS),begin),last=I::Multiply(I::Point(p.rateAtReferenceS),end);
      if(!I::Finite(first)||!I::Finite(last))return Reject(error,"nonfinite bounded photon rate");
      const long double v0=Midpoint(first),v1=Midpoint(last);rates.push_back({p.reaction,v0,v1});
      projectionBound=I::Add(projectionBound,I::Multiply(I::Multiply(durationBound,I::Point(.5L)),I::Add(I::Point(ProjectionError(first,v0)),I::Point(ProjectionError(last,v1)))));}
    const auto hInterval=durationBound;
    const long double analytic=w.visible?I::Divide(I::Multiply(I::Point(totalPhotons),I::Multiply(I::Point(curve),I::Multiply(hInterval,I::Multiply(hInterval,hInterval)))),I::Point(12)).hi:0;
    const long double bound=I::Add(I::Point(analytic),projectionBound).hi;
    const long double allowed=I::Divide(I::Multiply(I::Point(smoothBudget),I::Point(durationBound.lo)),durationInterval).lo;
    if(bound>allowed){const long double middle=w.lo+h/2;if(!(middle>w.lo&&middle<w.hi))return Reject(error,"exposure refinement is no longer representable");if(pending.size()+result.segments.size()+2>c.maximumSegments)return Reject(error,"exposure maximumSegments exhausted before declared error budget; reduce physical substep or allocate more segments");pending.push_back({middle,w.hi,w.visible});pending.push_back({w.lo,middle,w.visible});continue;}
    PiecewiseHazard::Segment segment;segment.beginS=w.lo;segment.endS=w.hi;
    segment.channels=std::move(rates);
    for(const auto&p:constants)segment.channels.push_back({p.reaction,p.rateS,p.rateS});
    result.interpolationDepthErrorBound=I::Add(I::Point(result.interpolationDepthErrorBound),I::Point(bound)).hi;result.segments.push_back(std::move(segment));
    if(result.segments.size()>c.maximumSegments)return Reject(error,"exposure maximumSegments exhausted");
  }
  if(!PiecewiseHazard::ValidateHistory(result.segments,error))return false;
  if(I::Add(I::Point(result.interpolationDepthErrorBound),I::Point(result.shadowBoundaryDepthErrorBound)).hi>budget)
    return Reject(error,"accumulated outward exposure-error bound exceeds the declared budget");
  *output=std::move(result);error->clear();return true;
}
}}}

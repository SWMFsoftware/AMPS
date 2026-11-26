#include "sep_reduced_geometry.h"
#include "pic_field_line_geometry.h"
#include "sep_geometric_spacing.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <iomanip>
#include <sstream>

namespace SEP { namespace ReducedShock {
namespace G=PIC::FieldLine::Geometry;
namespace {
Core::Result<Position> Failure(const std::string& why) {
  return Core::Result<Position>::Failure(Core::StatusCode::NumericalFailure,why);
}
Core::Result<Position> RK4(const VelocityQuery& velocity,const Position& x,
    double t,double end) {
  const double h=end-t,middle=t+0.5*h;
  Position p,k[4];
  auto query=velocity(x,t);if(!query.ok())return query;k[0]=query.value;
  for(int stage=1;stage<4;++stage) {
    const double fraction=stage==3?1:0.5;
    for(int d=0;d<3;++d)p[d]=x[d]+fraction*h*k[stage-1][d];
    // The caller owns the exact closed endpoint. Reconstructing it by adding
    // two rounded half-steps can query one ULP outside a valid provider view.
    // Carry the endpoint itself through both full/refined RK maps instead;
    // the strict provider guard remains unchanged and no time is clamped.
    query=velocity(p,stage==3?end:middle);if(!query.ok())return query;k[stage]=query.value;
  }
  for(int d=0;d<3;++d)p[d]=x[d]+h*(k[0][d]+2*k[1][d]+2*k[2][d]+k[3][d])/6;
  return G::Finite(p)?Core::Result<Position>::Success(p):Failure("trajectory produced nonfinite position");
}
}
Core::Result<Position> IntegrateTrajectory(const VelocityQuery& velocity,
    const Position& x,double departure,double arrival,const GeometryControls& c) {
  if(!velocity||!G::Finite(x)||!std::isfinite(departure)||!std::isfinite(arrival)||
      !std::isfinite(arrival-departure)||!std::isfinite(c.positionToleranceM)||c.positionToleranceM<=0||
      !std::isfinite(c.maximumTimeSubstepS)||c.maximumTimeSubstepS<=0||c.maximumAttempts==0)
    return Failure("invalid trajectory interval/accuracy controls");
  Position current=x;double time=departure;
  const double total=std::abs(arrival-departure);
  if(total==0)return Core::Result<Position>::Success(current);
  // Partition the COMPLETE requested interval by bisection, retaining each
  // exact endpoint. Taking repeated maximum-size steps first can leave an
  // arbitrarily short tail. At a heliocentric coordinate, its proportional
  // error budget can fall below coordinate roundoff: rejecting/halving that
  // tail reduces its budget again and can consume the entire work limit.
  // Bisection avoids creating that unrelated tiny remainder. It changes only
  // the internal RK mesh; neither physical AMPS dt nor the requested endpoint,
  // spatial tolerance, maximum numerical step or work bound is increased.
  //
  // Process the left (departure-side, also for a reverse query) child first.
  // The parent's endpoint stays pending until both children complete. Thus
  // accepted intervals are adjacent and cover exactly [departure,arrival];
  // no accumulated subtraction decides whether a final remainder exists.
  std::vector<double> endpoints{arrival};
  // Failure diagnostics only: keep the actual controller history so a work
  // limit can be reproduced from its SI inputs instead of guessed from a
  // downstream mover trap. These counters never select a step or change the
  // accepted position/error budget, provider query or physical timestep.
  std::size_t accepted=0,rejected=0,subdivided=0;double acceptedErrorM=0,minimumStepS=total;
  for(std::size_t attempt=0;!endpoints.empty()&&attempt<c.maximumAttempts;++attempt) {
    const double end=endpoints.back();
    const double h=end-time;
    if(h==0)return Failure("trajectory time lost floating-point progress");
    const double middle=time+0.5*h;
    const double step=std::abs(h);minimumStepS=std::min(minimumStepS,step);
    // Even maximum-step partitioning consumes the existing work budget. A
    // failed refinement must not silently accept an unresolved interval or
    // replace its endpoint with an epsilon-shifted/clamped time.
    const auto bisect=[&]() {
      if(middle==time||middle==end)return false;
      endpoints.push_back(middle);return true;
    };
    if(step>c.maximumTimeSubstepS) {
      ++subdivided;
      if(!bisect())return Failure("trajectory maximum-step partition lost floating-point progress");
      continue;
    }
    const auto full=RK4(velocity,current,time,end);if(!full.ok())return full;
    const auto half=RK4(velocity,current,time,middle);if(!half.ok())return half;
    const auto refined=RK4(velocity,half.value,middle,end);if(!refined.ok())return refined;
    const double estimate=G::Norm(G::Difference(full.value,refined.value))/15;
    // Keep the original proportional local-error criterion. Adjacent accepted
    // intervals partition the total duration, so their allocated estimates
    // sum to the same global spatial budget; no unused-budget redistribution
    // or roundoff acceptance floor is introduced.
    if(estimate>c.positionToleranceM*(step/total)) {
      ++rejected;
      if(!bisect())return Failure("trajectory accuracy refinement lost floating-point progress");
      continue;
    }
    ++accepted;acceptedErrorM+=estimate;
    current=refined.value;time=end;
    endpoints.pop_back();
  }
  if(!endpoints.empty()) {
    std::ostringstream why;why<<std::setprecision(17)
      <<"trajectory exceeded its declared work limit: initial_m="<<x[0]<<','<<x[1]<<','<<x[2]
      <<" departure_s="<<departure<<" arrival_s="<<arrival<<" last_time_s="<<time
      <<" remaining_s="<<std::abs(arrival-time)<<" minimum_attempt_s="<<minimumStepS
      <<" accepted="<<accepted<<" rejected="<<rejected<<" accepted_error_estimate_m="<<acceptedErrorM
      <<" maximum_step_subdivisions="<<subdivided
      <<" position_tolerance_m="<<c.positionToleranceM<<" maximum_time_substep_s="<<c.maximumTimeSubstepS
      <<" maximum_attempts="<<c.maximumAttempts;
    return Failure(why.str());
  }
  return Core::Result<Position>::Success(current);
}
Core::Result<Position> MoveVertex(const IntervalView& view,const Position& x,
    double departure,double arrival,const GeometryControls& c) {
  if(departure<view.DepartureS()||departure>view.ArrivalS()||
      arrival<view.DepartureS()||arrival>view.ArrivalS())
    return Core::Result<Position>::Failure(Core::StatusCode::OutOfDomain,"trajectory exceeds captured physical interval");
  // Validate even a zero-duration query, which otherwise executes no stages.
  const auto initial=view.PrimitiveAt(x,departure);
  if(!initial.ok())return Core::Result<Position>::Failure(initial.status.code,initial.status.message);
  return IntegrateTrajectory([&](const Position& p,double epoch) {
    const auto sample=view.PrimitiveAt(p,epoch);
    if(!sample.ok())return Core::Result<Position>::Failure(sample.status.code,sample.status.message);
    return Core::Result<Position>::Success(sample.value.velocityMPerS);
  },x,departure,arrival,c);
}
Core::Result<Position> TraceMagnetic(const IntervalView& view,const Position& x,
    const Position& outwardChord,double arc,double epoch,const GeometryControls& c) {
  Position target{}; // Published only after the tracer's complete success.
  const auto status=G::TraceAlongField(x,outwardChord,arc,epoch,
      [&](const Position& p,double t,Position* b) {
        const auto sample=view.PrimitiveAt(p,t);
        if(!sample.ok())return G::Status::Error(sample.status.message);
        *b=sample.value.magneticFieldT;return G::Status{};
      },{c.maximumArcSubstepM,c.positionToleranceM,c.maximumAttempts},&target);
  return status.ok()?Core::Result<Position>::Success(target):Failure(status.message);
}
Core::Result<std::vector<Position>> BuildInitialLine(const IntervalView& view,
    const Position& initial,double length,std::size_t count,double epoch,
    const GeometryControls& c,double firstArc) {
  using Return=Core::Result<std::vector<Position>>;
  if(!std::isfinite(length)||length<=0||count<2)
    return Return::Failure(Core::StatusCode::InvalidConfiguration,"invalid core line length or vertex count");
  std::vector<double> arcs(count-1,length/static_cast<double>(count-1));
  if(firstArc!=0) {
    std::string error;
    if(!GeometricArcs(length,count,firstArc,&arcs,&error))
      return Return::Failure(Core::StatusCode::InvalidConfiguration,error);
  }
  std::vector<Position> points;points.reserve(count);points.push_back(initial);
  for(std::size_t i=1;i<count;++i) {
    const double arc=arcs[i-1];
    const auto state=view.PrimitiveAt(points.back(),epoch);
    if(!state.ok())return Return::Failure(state.status.code,state.status.message);
    if(state.value.magneticSector==0)
      return Return::Failure(Core::StatusCode::UnsupportedCapability,"material open-line initialization cannot infer orientation on a closed core tube");
    Position reverse;
    for(int d=0;d<3;++d)reverse[d]=-state.value.magneticSector*state.value.magneticFieldT[d];
    // TraceMagnetic follows opposite its supplied oriented chord; reversing
    // the outward orientation here builds the initial line outward along B.
    const auto next=TraceMagnetic(view,points.back(),reverse,arc,epoch,c);
    if(!next.ok())return Return::Failure(next.status.code,next.status.message);
    points.push_back(next.value);
  }
  return Return::Success(std::move(points));
}
Core::Result<Position> TraceSunwardToRadius(const IntervalView& view,const Position& anchor,
    double radius,double epoch,const GeometryControls& c,double radiusTolerance,std::size_t maximumSteps) {
  using Return=Core::Result<Position>;
  double r=G::Norm(anchor);
  if(!G::Finite(anchor)||!std::isfinite(radius)||!(radius>0)||!(r>radius)||
      !std::isfinite(c.positionToleranceM)||!(c.positionToleranceM>0)||
      !std::isfinite(radiusTolerance)||!(radiusTolerance>0))
    return Return::Failure(Core::StatusCode::InvalidConfiguration,
        "sunward trace requires a finite anchor outside a positive target radius");
  Position x=anchor;
  for(std::size_t step=0;step<maximumSteps;++step) {
    const double remaining=r-radius;
    if(remaining<=radiusTolerance)return Return::Success(x);
    const auto state=view.PrimitiveAt(x,epoch);
    if(!state.ok())return Return::Failure(state.status.code,state.status.message);
    if(state.value.magneticSector==0)
      return Return::Failure(Core::StatusCode::UnsupportedCapability,
          "sunward trace reached a closed core tube; the anchor line is not open");
    Position outward;
    for(int d=0;d<3;++d)outward[d]=state.value.magneticSector*state.value.magneticFieldT[d];
    // TraceMagnetic follows opposite its supplied outward chord: sunward.
    const auto next=TraceMagnetic(view,x,outward,remaining,epoch,c);
    if(!next.ok())return Return::Failure(next.status.code,next.status.message);
    const double rn=G::Norm(next.value);
    if(!(rn<r))return Return::Failure(Core::StatusCode::NumericalFailure,
        "sunward magnetic trace does not decrease heliocentric radius");
    if(rn<radius-c.positionToleranceM-radiusTolerance)return Return::Failure(Core::StatusCode::NumericalFailure,
        "sunward magnetic trace crossed the target sphere beyond its tolerance");
    x=next.value;r=rn;
  }
  return Return::Failure(Core::StatusCode::NumericalFailure,
      "sunward magnetic trace did not reach the target radius within its step limit");
}
RadialBounds BuildRadialBounds(const std::vector<Position>& departure,const Position& center) {
  RadialBounds bounds;bounds.center=center;
  const std::size_t segments=departure.size()>1?departure.size()-1:0;
  bounds.prefixMaximumHigh.resize(segments);bounds.suffixMinimumLow.resize(segments);
  bounds.finite=segments>0;
  std::vector<double> low(segments);
  bounds.arc.assign(departure.size(),0.0);
  for(std::size_t i=0;i<segments;++i) {
    const double a=G::Norm(G::Difference(departure[i],center));
    const double b=G::Norm(G::Difference(departure[i+1],center));
    const double chord=G::Norm(G::Difference(departure[i+1],departure[i]));
    const double high=std::max(a,b);low[i]=std::min(a,b)-0.5*chord;
    bounds.arc[i+1]=bounds.arc[i]+chord;
    bounds.finite=bounds.finite&&std::isfinite(high)&&std::isfinite(low[i]);
    bounds.prefixMaximumHigh[i]=i?std::max(bounds.prefixMaximumHigh[i-1],high):high;
  }
  for(std::size_t k=segments;k-->0;)
    bounds.suffixMinimumLow[k]=k+1<segments?std::min(bounds.suffixMinimumLow[k+1],low[k]):low[k];
  return bounds;
}
VertexRange CandidateVertices(const RadialBounds& bounds,double inner,double outer,double reach) {
  VertexRange range;
  const std::size_t segments=bounds.prefixMaximumHigh.size();
  if(segments==0)return range;
  // Nonfinite geometry is never excluded: the full checks must see it.
  if(!bounds.finite||!std::isfinite(inner)||!std::isfinite(outer)||!std::isfinite(reach)) {
    range.count=segments+1;return range;
  }
  // Segments before 'first' satisfy max(r_i,r_i+1)+reach < inner; segments
  // after 'last' satisfy min(r_i,r_i+1)-l_i/2-2 reach > outer: neither can
  // meet the band (moved chord radii: see CandidateVertices above).
  const auto& high=bounds.prefixMaximumHigh;const auto& low=bounds.suffixMinimumLow;
  const std::size_t first=static_cast<std::size_t>(
      std::lower_bound(high.begin(),high.end(),inner-reach)-high.begin());
  const std::size_t end=static_cast<std::size_t>(
      std::upper_bound(low.begin(),low.end(),outer+2*reach)-low.begin());
  if(first>=end)return range;
  range.first=first;range.count=end-first+1;
  return range;
}
VertexRange CandidateVertices(const std::vector<Position>& departure,
    const Position& center,double inner,double outer,double reach) {
  return CandidateVertices(BuildRadialBounds(departure,center),inner,outer,reach);
}
Core::Result<std::vector<Connection>> ConnectLine(const IntervalView& view,
    const std::vector<Position>& points,double epoch,double tolerance,
    std::size_t firstIndex,std::size_t lineVertexCount) {
  const std::size_t total=lineVertexCount?lineVertexCount:points.size();
  using Return=Core::Result<std::vector<Connection>>;
  if(points.size()<2||!std::isfinite(tolerance)||tolerance<=0)
    return Return::Failure(Core::StatusCode::InvalidConfiguration,"invalid reduced line intersection inputs");
  const auto sphere=view.SphereAt(epoch);if(!sphere.ok())return Return::Failure(sphere.status.code,sphere.status.message);
  std::vector<Connection> result;double arc=0;
  const double eps=std::numeric_limits<double>::epsilon();
  for(std::size_t segment=0;segment+1<points.size();++segment) {
    const auto& x=points[segment];const auto& y=points[segment+1];
    const auto delta=G::Difference(y,x);const double length=G::Norm(delta);
    if(!G::Finite(x)||!G::Finite(y)||!std::isfinite(length)||length<=tolerance)
      return Return::Failure(Core::StatusCode::InvalidState,"unresolved or nonfinite reduced line segment");
    // Dimensionless quadratic, accumulated in long double. The uncertainty
    // band accounts for the DOUBLE precision supplied geometry, rather than
    // mistaking long-double accumulation for more accurate physical inputs.
    long double a=0,b=0,c=-1;
    for(int d=0;d<3;++d) {
      const long double q=(static_cast<long double>(x[d])-sphere.value.centerM[d])/sphere.value.radiusM;
      const long double v=static_cast<long double>(delta[d])/sphere.value.radiusM;
      a+=v*v;b+=2*q*v;c+=q*q;
    }
    long double discriminant=b*b-4*a*c;
    const long double uncertainty=128*eps*(b*b+4*std::abs(a*c));
    if(discriminant < -uncertainty){arc+=length;continue;}
    const bool unresolvedTangent=std::abs(discriminant)<=uncertainty;
    if(unresolvedTangent)discriminant=0;
    const long double root=std::sqrt(std::max(0.L,discriminant));
    const long double q=-0.5L*(b+std::copysign(root,b));
    std::vector<long double> fractions;
    if(root==0)fractions.push_back(-b/(2*a));
    else {fractions.push_back(q/a);fractions.push_back(q==0?-b/a:c/q);std::sort(fractions.begin(),fractions.end());}
    for(long double value:fractions) {
      // Only coordinate-roundoff excursions coalesce at an endpoint. An
      // actually exterior root is never pulled into the computational line.
      if(value < -32*eps||value>1+32*eps)continue;
      double fraction=static_cast<double>(std::max(0.L,std::min(1.L,value)));
      if(std::abs(fraction)<32*eps)fraction=0;
      if(std::abs(fraction-1)<32*eps)fraction=1;
      Position point;for(int d=0;d<3;++d)point[d]=x[d]+fraction*delta[d];
      const double radius=G::Norm(point);
      if(radius==0)continue;
      const auto normalRadius=G::Difference(point,sphere.value.centerM);
      if(G::Dot(point,sphere.value.direction)/radius<sphere.value.minimumRayCosine-64*eps||
          G::Dot(point,normalRadius)<-64*eps*radius*sphere.value.radiusM)continue;
      const auto canonical=view.LeadingPointOnRay(point,epoch);
      if(!canonical.ok())return Return::Failure(canonical.status.code,canonical.status.message);
      const double residual=G::Norm(G::Difference(point,canonical.value));
      if(residual>tolerance)return Return::Failure(Core::StatusCode::NumericalFailure,
          "supported reduced intersection does not meet its spatial tolerance");
      const auto front=view.FrontAt(canonical.value,epoch);
      if(!front.ok())return Return::Failure(front.status.code,front.status.message);
      Connection hit;hit.segment=firstIndex+segment;hit.fraction=fraction;hit.arcLengthM=arc+fraction*length;
      hit.surfaceResidualM=residual;hit.front=front.value;
      hit.tangent=unresolvedTangent||std::abs(G::Dot(delta,front.value.outwardNormal))<=64*eps*length;
      if(fraction==1&&firstIndex+segment+2<total){hit.segment=firstIndex+segment+1;hit.fraction=0;}
      if(!result.empty()&&result.back().segment==hit.segment&&result.back().fraction==hit.fraction)continue;
      result.push_back(std::move(hit));
    }
    arc+=length;
  }
  return Return::Success(std::move(result));
}
Core::Result<std::vector<Connection>> ConnectMaterialLine(const IntervalView& view,
    const std::vector<Position>& departure,double epoch,double tolerance,const GeometryControls& controls,
    const RadialBounds* bounds) {
  using Return=Core::Result<std::vector<Connection>>;
  if(departure.size()<2||!std::isfinite(tolerance)||tolerance<=0)
    return Return::Failure(Core::StatusCode::InvalidConfiguration,"invalid reduced line intersection inputs");
  const auto sphere=view.SphereAt(epoch);if(!sphere.ok())return Return::Failure(sphere.status.code,sphere.status.message);
  // Every vertex moves at most (speed bound)*|elapsed| plus the trajectory's
  // own position tolerance; the intersection tolerance widens the band too.
  const double reach=view.SpeedBoundMPerS()*std::abs(epoch-view.DepartureS())+
      controls.positionToleranceM+tolerance;
  const bool reuse=bounds&&bounds->center==sphere.value.centerM&&
      bounds->prefixMaximumHigh.size()+1==departure.size();
  const auto range=reuse?CandidateVertices(*bounds,sphere.value.radiusM,sphere.value.radiusM,reach):
      CandidateVertices(departure,sphere.value.centerM,sphere.value.radiusM,sphere.value.radiusM,reach);
  if(range.count==0)return Return::Success({});
  std::vector<Position> moved;moved.reserve(range.count);
  for(std::size_t i=range.first;i<range.first+range.count;++i) {
    const auto target=MoveVertex(view,departure[i],view.DepartureS(),epoch,controls);
    if(!target.ok())return Return::Failure(target.status.code,target.status.message);
    moved.push_back(target.value);
  }
  auto connections=ConnectLine(view,moved,epoch,tolerance,range.first,departure.size());
  if(!connections.ok())return connections;
  // Material coordinate of each root on the departure chords.
  if(!connections.value.empty()) {
    const RadialBounds local=reuse?RadialBounds{}:BuildRadialBounds(departure,sphere.value.centerM);
    const auto& arc=reuse?bounds->arc:local.arc;
    for(auto& hit:connections.value)
      hit.arcLengthM=arc[hit.segment]+hit.fraction*(arc[hit.segment+1]-arc[hit.segment]);
  }
  return connections;
}
} }

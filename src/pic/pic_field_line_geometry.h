#ifndef AMPS_PIC_FIELD_LINE_GEOMETRY_H
#define AMPS_PIC_FIELD_LINE_GEOMETRY_H

#include <array>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace PIC { namespace FieldLine { namespace Geometry {
using Point=std::array<double,3>;
struct Status {
  std::string message;
  bool ok() const {return message.empty();}
  static Status Error(const std::string& reason) {Status s;s.message=reason;return s;}
};
inline double Dot(const Point& a,const Point& b) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
inline Point Difference(const Point& a,const Point& b) {return {{a[0]-b[0],a[1]-b[1],a[2]-b[2]}};}
inline bool Finite(const Point& p) {return std::isfinite(p[0])&&std::isfinite(p[1])&&std::isfinite(p[2]);}
inline double Norm(const Point& p) {return std::hypot(std::hypot(p[0],p[1]),p[2]);}

struct Context {
  double departureS=0,arrivalS=0;
  Point centerM{{0,0,0}};
  double minimumRadiusM=0,maximumRadiusM=0;
  // Explicit numerical controls, not hidden repairs. Exceeding a limit fails
  // the entire candidate; it never silently clips geometry or changes dt.
  double toleranceM=0;
  std::size_t maximumInsertions=0;
  std::uint64_t expectedGeneration=0;
};
// The geometry owner integrates motion over the full interval. A callback
// returns a position [m], not a speed that this container integrates again.
using Motion=std::function<Status(std::size_t,const Point&,double,double,Point*)>;
// Trace toward the beginning along the owner's target-epoch B. The oriented
// first chord resolves B polarity; length is requested magnetic arc length.
using InnerTrace=std::function<Status(const Point&,const Point&,double,double,Point*)>;
using MagneticField=std::function<Status(const Point&,double,Point*)>;
struct TraceControl {
  double maximumStepM=0,toleranceM=0;
  std::size_t maximumAttempts=0;
};
// Integrate dx/dl=+-B/|B| at one target epoch. The sign follows the material
// ordering, independently of magnetic polarity. RK4 step doubling controls
// spatial error without changing any physical particle timestep. Only a
// fully successful requested arc length is published to output.
inline Status TraceAlongField(const Point& first,const Point& outwardChord,
    double distanceM,double epochS,const MagneticField& field,
    const TraceControl& controls,Point* output) {
  if(!output||!field||!Finite(first)||!Finite(outwardChord)||
      !std::isfinite(distanceM)||distanceM<=0||!std::isfinite(epochS)||
      !std::isfinite(controls.maximumStepM)||controls.maximumStepM<=0||
      !std::isfinite(controls.toleranceM)||controls.toleranceM<=0||controls.maximumAttempts==0)
    return Status::Error("invalid inner magnetic tracing controls");
  Point initial;auto status=field(first,epochS,&initial);if(!status.ok())return status;
  const double projection=Dot(initial,outwardChord);
  if(!Finite(initial)||Norm(initial)==0||!std::isfinite(projection)||projection==0)
    return Status::Error("inner magnetic trace cannot resolve nonzero B and material orientation");
  const double sign=projection>0?-1:1;
  auto direction=[&](const Point& x,Point* d)->Status {
    Point b;const auto evaluated=field(x,epochS,&b);if(!evaluated.ok())return evaluated;
    const double magnitude=Norm(b);
    if(!Finite(b)||!std::isfinite(magnitude)||magnitude<=0)
      return Status::Error("inner magnetic trace encountered invalid/null B");
    for(int i=0;i<3;++i)(*d)[i]=sign*b[i]/magnitude;
    return Status{};
  };
  auto rk4=[&](const Point& x,double h,Point* end)->Status {
    Point k1,k2,k3,k4,p;
    auto s=direction(x,&k1);if(!s.ok())return s;
    for(int d=0;d<3;++d)p[d]=x[d]+h*k1[d]/2;
    s=direction(p,&k2);if(!s.ok())return s;
    for(int d=0;d<3;++d)p[d]=x[d]+h*k2[d]/2;
    s=direction(p,&k3);if(!s.ok())return s;
    for(int d=0;d<3;++d)p[d]=x[d]+h*k3[d];
    s=direction(p,&k4);if(!s.ok())return s;
    for(int d=0;d<3;++d)(*end)[d]=x[d]+h*(k1[d]+2*k2[d]+2*k3[d]+k4[d])/6;
    return Finite(*end)?Status{}:Status::Error("nonfinite RK magnetic trace result");
  };
  Point x=first;double remaining=distanceM,h=std::min(distanceM,controls.maximumStepM);
  for(std::size_t attempt=0;remaining>0&&attempt<controls.maximumAttempts;++attempt) {
    h=std::min(h,remaining);Point full,half,refined;
    status=rk4(x,h,&full);if(!status.ok())return status;
    status=rk4(x,h/2,&half);if(!status.ok())return status;
    status=rk4(half,h/2,&refined);if(!status.ok())return status;
    const double error=Norm(Difference(full,refined))/15;
    const double allowed=controls.toleranceM*(h/distanceM);
    if(error>allowed) {
      h/=2;
      if(h==0||Norm(Difference(refined,x))==0)
        return Status::Error("inner trace cannot resolve the requested spatial tolerance");
      continue;
    }
    x=refined;
    if(h==remaining)remaining=0;
    else {
      const double next=remaining-h;
      if(next==remaining)return Status::Error("inner trace arc length lost floating-point progress");
      remaining=next;
    }
    h=std::min(controls.maximumStepM,h*2);
  }
  if(remaining>0)return Status::Error("inner magnetic trace exceeded its explicit work limit");
  *output=x;return Status{};
}
struct Plan {
  Context context;
  std::vector<Point> oldPositions,movedPositions,insertedPositions;
  // insertedPositions is in insertion order (nearest old first first).
  // Retained old segment identities survive; indices shift by prefix count.
  std::size_t removedBack=0;
};

// A chord can cross the sphere even when both endpoint radii are outside.
// Use its actual closest point. Degenerate/nonfinite input is rejected by
// preparation; this helper also rejects it for independent direct callers.
inline double MinimumRadius(const Point& a,const Point& b,const Point& center) {
  const Point d=Difference(b,a),r=Difference(center,a);
  const double length=Norm(d);
  if(!Finite(a)||!Finite(b)||!Finite(center)||!std::isfinite(length)||length<=0)
    throw std::invalid_argument("invalid material segment in radial domain test");
  const Point direction{{d[0]/length,d[1]/length,d[2]/length}};
  const double along=std::max(0.0,std::min(length,Dot(r,direction)));
  return Norm(Point{{a[0]+along*direction[0]-center[0],
      a[1]+along*direction[1]-center[1],a[2]+along*direction[2]-center[2]}});
}

// Preparation has no observable geometry mutation. Callback errors, domain
// holes and invalid later vertices leave both the live line and output plan
// unchanged. The native wrapper additionally stages physical data/particles.
inline Status Prepare(const std::vector<Point>& old,const Context& c,
                      const Motion& move,const InnerTrace& trace,Plan* output) {
  if(!output||!move||!trace||old.size()<2||!Finite(c.centerM)||
      !std::isfinite(c.departureS)||!std::isfinite(c.arrivalS)||c.arrivalS<=c.departureS||
      !std::isfinite(c.arrivalS-c.departureS)||
      !std::isfinite(c.minimumRadiusM)||c.minimumRadiusM<=0||
      !std::isfinite(c.maximumRadiusM)||c.maximumRadiusM<=c.minimumRadiusM||
      !std::isfinite(c.toleranceM)||c.toleranceM<=0||c.toleranceM>=c.minimumRadiusM||
      c.maximumInsertions==0)
    return Status::Error("invalid material geometry interval/domain/callback controls");
  Plan candidate;candidate.context=c;candidate.oldPositions=old;
  candidate.movedPositions.resize(old.size());
  try {
    for(std::size_t i=0;i<old.size();++i) {
      if(!Finite(old[i]))return Status::Error("nonfinite departure material vertex");
      if(i&&(!std::isfinite(Norm(Difference(old[i],old[i-1])))||
          Norm(Difference(old[i],old[i-1]))<=c.toleranceM))
        return Status::Error("degenerate/nonfinite departure material segment");
      const auto status=move(i,old[i],c.departureS,c.arrivalS,&candidate.movedPositions[i]);
      if(!status.ok())return status;
      if(!Finite(candidate.movedPositions[i]))return Status::Error("motion returned a nonfinite material vertex");
      if(i&&(!std::isfinite(Norm(Difference(candidate.movedPositions[i],candidate.movedPositions[i-1])))||
          Norm(Difference(candidate.movedPositions[i],candidate.movedPositions[i-1]))<=c.toleranceM))
        return Status::Error("motion collapsed a material segment");
    }
    Point first=candidate.movedPositions[0],second=candidate.movedPositions[1];
    while(Norm(Difference(first,c.centerM))>c.minimumRadiusM+c.toleranceM) {
      if(candidate.insertedPositions.size()>=c.maximumInsertions)
        return Status::Error("material inner extension exceeded its explicit insertion limit");
      const Point chord=Difference(second,first);const double length=Norm(chord);
      Point next;
      const auto status=trace(first,chord,length,c.arrivalS,&next);
      if(!status.ok())return status;
      if(!Finite(next)||Norm(Difference(next,first))<=c.toleranceM)
        return Status::Error("inner magnetic trace returned invalid/degenerate geometry");
      if(Norm(Difference(next,c.centerM))>=Norm(Difference(first,c.centerM))-c.toleranceM)
        return Status::Error("inner magnetic trace made no inward progress or overshot the radial domain");
      candidate.insertedPositions.push_back(next);second=first;first=next;
    }
    // Trim the COMPLETE candidate, including any new prefix records still
    // beyond the outer sphere after a large displacement. All old material
    // may leave while newly entering material still covers the physical line;
    // that is a valid loss event, not permission to retain out-of-domain edges.
    std::vector<Point> all(candidate.insertedPositions.rbegin(),candidate.insertedPositions.rend());
    all.insert(all.end(),candidate.movedPositions.begin(),candidate.movedPositions.end());
    for(std::size_t i=all.size()-1;i>0;--i) {
      if(MinimumRadius(all[i-1],all[i],c.centerM)
          <=c.maximumRadiusM+c.toleranceM)break;
      ++candidate.removedBack;
    }
    if(candidate.removedBack==all.size()-1)
      return Status::Error("the complete candidate line left the domain; explicit line retirement is required");
  }catch(const std::exception& e){return Status::Error(std::string("material geometry callback failed: ")+e.what());}
  *output=std::move(candidate);return Status{};
}

// Counts measure sampling opportunities, not particle counts. Empty segments
// accumulate exactly the same exposure as occupied segments of the same age.
// Collecting/completed buffers have independent lifetimes and swap together.
struct Exposure {
  std::uint64_t collecting=0,completed=0;
  double collectingSeconds=0,completedSeconds=0;
  double fixedStepSeconds=0;
  void Sample(double dt) {
    if(!std::isfinite(dt)||dt<=0||collecting==std::numeric_limits<std::uint64_t>::max()||
        !std::isfinite(collectingSeconds+dt))
      throw std::invalid_argument("invalid/overflowing material sampling exposure");
    // The lifetime denominator assumes one unchanged configured physical dt.
    // Recording exposure is not permission to introduce remainder timesteps.
    if(fixedStepSeconds!=0&&fixedStepSeconds!=dt)
      throw std::invalid_argument("material sampling physical timestep changed");
    fixedStepSeconds=dt;
    ++collecting;collectingSeconds+=dt;
  }
  void Complete() {
    completed=collecting;completedSeconds=collectingSeconds;
    collecting=0;collectingSeconds=0;
  }
  // Absence is explicit. A newly created segment has no completed measurement,
  // which must not be confused with a sampled vacuum value of zero.
  bool HasCompletedMeasurement() const {return completed!=0;}
  double Average(double numerator) const {
    if(!HasCompletedMeasurement()||!std::isfinite(numerator))
      throw std::invalid_argument("material sample average has no valid completed exposure");
    return numerator/static_cast<double>(completed);
  }
};
} } }
#endif

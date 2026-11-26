#include "reduced_shock_boundary.h"
#include "reduced_shock_geometry.h"
#include "util/sep_run_configuration.h"
#include <algorithm>
#include <cmath>
#include <map>

namespace SEP { namespace ReducedParticles { namespace Native {
namespace {
using TStatus=Transport::Status;
TStatus Error(const std::string& why){return TStatus::Error(Transport::StatusCode::UnsupportedConfiguration,why);}
double EpochAt(const Transport::PICAdapter::ParticleContext& context,double elapsed,double offset=0) {
  // Use the birth clock's end-relative duration arithmetic. Adding a rounded
  // birth epoch, accumulated elapsed time and the last remainder can exceed
  // the immutable endpoint by one ULP. Subtract the remaining duration FIRST:
  // an accepted final remainder then queries the exact declared endpoint.
  // No tolerance, timestamp clamp, extrapolation or physical-dt change occurs.
  if(context.usesDriverWindow||context.consumesPendingBirth)
    return context.firstMove.parent.endEpochS-((context.moveDurationS-elapsed)-offset);
  return context.particleStepEpochS+elapsed+offset;
}
struct Root {bool active=false;double arc=0;int sign=1;};
Core::Result<Root> PolarRoot(const ReducedShock::Native::MaterialView& mesh,double epoch,
    double shock) {
  using R=Core::Result<Root>;Root root;
  const auto controls=ReducedShock::Native::Controls(mesh.interval.ArrivalS()-mesh.interval.DepartureS(),mesh.arc[1]);
  const auto height=[&](std::size_t i) {
    return ReducedShock::MoveVertex(mesh.interval,mesh.vertices[i],mesh.interval.DepartureS(),epoch,controls);
  };
  const auto first=height(0),last=height(mesh.vertices.size()-1);
  if(!first.ok())return R::Failure(first.status.code,first.status.message);
  if(!last.ok())return R::Failure(last.status.code,last.status.message);
  if(shock<first.value[2]||shock>last.value[2])return R::Success(root);
  std::size_t left=0,right=mesh.vertices.size()-1;
  auto a=first.value,b=last.value;
  // Use the stored departure geometry only as a SEARCH HINT. The HCI shock
  // height often lies in the same or a neighbouring departure segment, so
  // integrating every midpoint of a full-range binary search is redundant.
  // Retain the actual first/last flow queries above, then validate the hinted
  // vertices with the SAME controlled core trajectories at the EXACT epoch.
  // Only a proved bracket can replace the current range; otherwise ordinary
  // binary search continues in the remaining proved range. No U*dt estimate,
  // approximate root, temporal interpolation or omitted validity gate occurs.
  const auto upper=std::upper_bound(mesh.vertices.begin(),mesh.vertices.end(),shock,
      [](double z,const ReducedShock::Position& point){return z<point[2];});
  const std::size_t hintedRight=static_cast<std::size_t>(upper-mesh.vertices.begin());
  if(hintedRight>0&&hintedRight<mesh.vertices.size()) {
    const std::size_t hintedLeft=hintedRight-1;
    const auto lowerPoint=hintedLeft==0?first:height(hintedLeft);
    const auto upperPoint=hintedRight==right?last:height(hintedRight);
    if(!lowerPoint.ok())return R::Failure(lowerPoint.status.code,lowerPoint.status.message);
    if(!upperPoint.ok())return R::Failure(upperPoint.status.code,upperPoint.status.message);
    if(lowerPoint.value[2]>shock) {
      right=hintedLeft;b=lowerPoint.value;
      // A material displacement of less than one departure segment commonly
      // moves the true bracket one segment inward from the first hint. Probe
      // that neighbour with the full trajectory before paying for a binary
      // search. Accept it ONLY if its actual height proves the lower bound;
      // otherwise retain the narrowed proved range and continue normally.
      if(right>0) {
        const auto neighbour=height(right-1);
        if(!neighbour.ok())return R::Failure(neighbour.status.code,neighbour.status.message);
        if(neighbour.value[2]<=shock){left=right-1;a=neighbour.value;}
        else {--right;b=neighbour.value;}
      }
    }
    else if(upperPoint.value[2]<=shock&&hintedRight<right){left=hintedRight;a=upperPoint.value;}
    else {left=hintedLeft;right=hintedRight;a=lowerPoint.value;b=upperPoint.value;}
    // Equality at an internal vertex belongs to the departure-side bound,
    // matching the retained binary search and native floor(S) convention.
    // The terminal vertex can remain the upper bound of the final segment.
  }
  while(right-left>1) {
    const std::size_t middle=left+(right-left)/2;const auto point=height(middle);
    if(!point.ok())return R::Failure(point.status.code,point.status.message);
    if(point.value[2]<=shock){left=middle;a=point.value;}else{right=middle;b=point.value;}
  }
  // On the positive solar-rotation axis, the qualified Parker exterior
  // has radial U and B. Its scalar flow map preserves vertex ordering, so
  // binary search brackets the core's leading ray intersection on the
  // actual advected chords. The fitted cap may have a tilted axis, and its
  // leading root on this line is then not the apex.
  if(a[0]!=0||a[1]!=0||b[0]!=0||b[1]!=0||!(b[2]>a[2]))
    return R::Failure(Core::StatusCode::InvalidState,"qualified polar material flow lost its axis/order invariant");
  const auto front=mesh.interval.FrontAt({{0,0,shock}},epoch);
  if(!front.ok())return R::Failure(front.status.code,front.status.message);
  if(front.value.supportEdge||!(front.value.outwardNormal[2]>0))
    return R::Failure(Core::StatusCode::UnsupportedCapability,
        "polar ray reached a grazing or unsupported leading barrier");
  if(!front.value.acceptedShock)return R::Success(root);
  root.active=true;root.arc=mesh.arc[left]+(shock-a[2])/(b[2]-a[2])*(mesh.arc[right]-mesh.arc[left]);
  return R::Success(root);
}
Core::Result<Root> EvaluateRootAt(const ReducedShock::Native::MaterialView& mesh,double epoch) {
  using R=Core::Result<Root>;
  const auto sphere=mesh.interval.SphereAt(epoch);
  if(!sphere.ok())return R::Failure(sphere.status.code,sphere.status.message);
  const auto& front=sphere.value;
  // Strict interior support proves there is one forward leading root
  // on the +Z ray. Tangent and outside-support cases retain the general
  // intersection path and its typed validity checks.
  bool polar=front.direction[2]>front.minimumRayCosine;
  double previous=0;
  for(const auto& point:mesh.vertices) {
    polar=polar&&point[0]==0&&point[1]==0&&point[2]>previous;previous=point[2];
  }
  if(polar) {
    const auto ambient=mesh.interval.PrimitiveAt(mesh.vertices.front(),epoch);
    if(!ambient.ok())return R::Failure(ambient.status.code,ambient.status.message);
    polar=ambient.value.parkerExterior&&ambient.value.velocityMPerS[0]==0&&ambient.value.velocityMPerS[1]==0;
    if(polar) {
      if(front.direction[2]==1&&front.direction[0]==0&&front.direction[1]==0&&
         front.centerM[0]==0&&front.centerM[1]==0)
        return PolarRoot(mesh,epoch,front.centerM[2]+front.radiusM);
      // Reuse the production finite-SSE solver at this actual epoch.
      // It applies the same support and picks the same leading branch
      // as the general chord-intersection path, without visiting all
      // of the unrelated advected segments on a proved radial line.
      const auto leading=mesh.interval.LeadingPointOnRay({{0,0,1}},epoch);
      if(!leading.ok())return R::Failure(leading.status.code,leading.status.message);
      return PolarRoot(mesh,epoch,leading.value[2]);
    }
  }
  // Curved/nonpolar lines retain the complete finite-front intersection path
  // and its multiple-root/support checks. No geometrical shortcut is inferred
  // from the current particle's position or from a function-pointer identity.
  const auto connections=ReducedShock::Native::ConnectionsAt(epoch);
  if(!connections.ok())return R::Failure(connections.status.code,connections.status.message);
  Root root;
  for(const auto& c:connections.value)if(c.front.acceptedShock) {
    if(root.active)return R::Failure(Core::StatusCode::UnsupportedCapability,
        "multiple simultaneous reduced shock barriers need a joint first-passage operator; this first consumer admits one transverse leading connection");
    if(c.tangent||c.front.supportEdge)return R::Failure(Core::StatusCode::UnsupportedCapability,
        "tangent/support-edge shock barrier has no qualified material-coordinate first-passage map");
    root.active=true;root.arc=mesh.arc[c.segment]+c.fraction*(mesh.arc[c.segment+1]-mesh.arc[c.segment]);
    const auto a=ReducedShock::MoveVertex(mesh.interval,mesh.vertices[c.segment],mesh.interval.DepartureS(),epoch,
        ReducedShock::Native::Controls(mesh.interval.ArrivalS()-mesh.interval.DepartureS(),mesh.arc[1]));
    const auto b=ReducedShock::MoveVertex(mesh.interval,mesh.vertices[c.segment+1],mesh.interval.DepartureS(),epoch,
        ReducedShock::Native::Controls(mesh.interval.ArrivalS()-mesh.interval.DepartureS(),mesh.arc[1]));
    if(!a.ok())return R::Failure(a.status.code,a.status.message);
    if(!b.ok())return R::Failure(b.status.code,b.status.message);
    double projection=0;for(int d=0;d<3;++d)projection+=(b.value[d]-a.value[d])*c.front.outwardNormal[d];
    if(projection==0)return R::Failure(Core::StatusCode::UnsupportedCapability,"material chord tangent to absorbing front");
    root.sign=projection>0?1:-1;
  }
  return R::Success(root);
}
// Verified piecewise-linear shock-root table over the retained step (opt-in,
// --shock-boundary-root-table). Nodes are exact roots. An interval [a,b] is
// accepted only when the exact roots at its 1/4, 1/2 and 3/4 epochs have the
// same activity and upstream sign as both nodes and, when active, lie within
// boundaryToleranceM of the linear interpolant of the node arcs; otherwise it
// is bisected. An interval that cannot be resolved above a 1e-9 relative
// width (e.g. a connectivity change inside it) is flagged for exact queries,
// so the table never interpolates across a topology change.
struct RootTable {
  std::shared_ptr<const ReducedShock::Native::MaterialView> owner;
  std::vector<double> epochs;   // node epochs, increasing, covering the view
  std::vector<Root> roots;      // exact roots at the nodes
  std::vector<char> exact;      // interval i -> answer queries exactly
};
Core::Status BuildRootTable(const std::shared_ptr<const ReducedShock::Native::MaterialView>& mesh,RootTable* table) {
  const double a=mesh->interval.DepartureS(),b=mesh->interval.ArrivalS();
  const double tolerance=Run::Active().get().reducedParticles.boundaryToleranceM;
  struct Node{double t;Root r;};
  const auto exactAt=[&](double t)->Core::Result<Root>{return EvaluateRootAt(*mesh,t);};
  std::vector<Node> done;std::vector<char> flags;
  const auto first=exactAt(a);if(!first.ok())return Core::Status::Failure(first.status.code,first.status.message);
  done.push_back({a,first.value});
  if(!(b>a)){table->owner=mesh;table->epochs={a};table->roots={first.value};table->exact={};return Core::Status::Success();}
  // Depth-first refinement from left to right keeps the node list ordered.
  std::vector<Node> pending;
  const int initial=8;
  for(int k=initial;k>=1;--k) {
    const double t=k==initial?b:a+(b-a)*k/initial;
    const auto r=exactAt(t);if(!r.ok())return Core::Status::Failure(r.status.code,r.status.message);
    pending.push_back({t,r.value});
  }
  while(!pending.empty()) {
    const Node left=done.back(),right=pending.back();
    const double width=right.t-left.t;
    bool accept=left.r.active==right.r.active&&left.r.sign==right.r.sign;
    if(accept)for(double f:{0.25,0.5,0.75}) {
      const auto m=exactAt(left.t+f*width);if(!m.ok())return Core::Status::Failure(m.status.code,m.status.message);
      if(m.value.active!=left.r.active||(m.value.active&&m.value.sign!=left.r.sign)){accept=false;break;}
      if(m.value.active&&std::abs(m.value.arc-(left.r.arc+f*(right.r.arc-left.r.arc)))>tolerance){accept=false;break;}
    }
    const double middle=left.t+0.5*width;
    if(accept||!(width>1e-9*(b-a))||middle<=left.t||middle>=right.t) {
      flags.push_back(accept?0:1);done.push_back(right);pending.pop_back();
      continue;
    }
    const auto m=exactAt(middle);if(!m.ok())return Core::Status::Failure(m.status.code,m.status.message);
    pending.push_back({middle,m.value});
  }
  table->owner=mesh;table->epochs.clear();table->roots.clear();
  for(const auto& node:done){table->epochs.push_back(node.t);table->roots.push_back(node.r);}
  table->exact=std::move(flags);
  return Core::Status::Success();
}
Core::Result<Root> TableRootAt(const std::shared_ptr<const ReducedShock::Native::MaterialView>& mesh,double epoch) {
  using R=Core::Result<Root>;
  static thread_local RootTable table;
  if(!std::isfinite(epoch)||epoch<mesh->interval.DepartureS()||epoch>mesh->interval.ArrivalS())
    return EvaluateRootAt(*mesh,epoch);
  if(table.owner!=mesh) {
    const auto built=BuildRootTable(mesh,&table);
    if(!built.ok()){table.owner.reset();return R::Failure(built.code,built.message);}
  }
  const auto& t=table.epochs;
  if(t.size()<2)return R::Success(table.roots.front());
  std::size_t i=static_cast<std::size_t>(std::upper_bound(t.begin(),t.end(),epoch)-t.begin());
  i=i==0?0:std::min(i-1,t.size()-2);
  if(epoch==t[i])return R::Success(table.roots[i]);
  if(epoch==t[i+1])return R::Success(table.roots[i+1]);
  if(table.exact[i])return EvaluateRootAt(*mesh,epoch);
  Root root=table.roots[i];
  if(root.active) {
    const double f=(epoch-t[i])/(t[i+1]-t[i]);
    root.arc=table.roots[i].arc+f*(table.roots[i+1].arc-table.roots[i].arc);
  }
  return R::Success(root);
}
Core::Result<Root> RootAt(const std::shared_ptr<const ReducedShock::Native::MaterialView>& mesh,double epoch) {
  // The shock's material coordinate depends only on this immutable line view
  // and the EXACT epoch, never on the particle or its random stream. Focused
  // particles often share the same numerical times; the previous endpoint is
  // also the next interval's start. Reintegrating the same vertex trajectories
  // and resolving the same RH jump for every repetition is unnecessary work.
  // Memoize successful values without rounding times, interpolating roots,
  // omitting physical checks on a new query, or changing the accuracy budget.
  //
  // Each worker owns its cache. Retain the shared immutable owner itself, not
  // just its address/generation: otherwise allocator address reuse after a
  // geometry commit could make an old material coordinate appear current.
  // Changing the view clears every entry, including a shared physical endpoint
  // whose arc label changes after inner insertion. The finite entry limit is
  // only a memory bound; eviction repeats the full existing calculation.
  struct Cache {
    std::shared_ptr<const ReducedShock::Native::MaterialView> owner;
    std::map<double,Root> roots;
  };
  static thread_local Cache cache;
  if(cache.owner!=mesh){cache.roots.clear();cache.owner=mesh;}
  if(Run::Active().get().reducedParticles.boundaryRootTable)return TableRootAt(mesh,epoch);
  // NaN cannot be used as an ordered-map key: it is not strictly ordered and
  // could compare equivalent to an unrelated valid timestamp. Let the actual
  // provider reject every nonfinite/out-of-interval query before any lookup.
  if(!std::isfinite(epoch)||epoch<mesh->interval.DepartureS()||epoch>mesh->interval.ArrivalS())
    return EvaluateRootAt(*mesh,epoch);
  const auto found=cache.roots.find(epoch);
  if(found!=cache.roots.end())return Core::Result<Root>::Success(found->second);
  const auto result=EvaluateRootAt(*mesh,epoch);
  if(result.ok()) {
    if(cache.roots.size()==4096)cache.roots.clear();
    cache.roots.emplace(epoch,result.value);
  }
  return result;
}
}
Transport::Status AbsorbInitialReturn(const Transport::PICAdapter::ParticleContext& context,
    double momentum,bool* removed) {
  if(!removed)return Error("initial shock-return query requires an output flag");
  *removed=false;if(!Enabled())return TStatus::Ok();
  const auto retained=ReducedShock::Native::CaptureMaterialView();
  if(!retained.ok())return Error(retained.status.message);
  const auto& mesh=*retained.value;
  const auto root=RootAt(retained.value,EpochAt(context,0));
  if(!root.ok())return Error(root.status.message);
  if(!root.value.active)return TStatus::Ok();
  const auto segment=static_cast<std::size_t>(context.state.coordinate);
  if(segment+1>=mesh.arc.size())return Error("initial return particle is outside its material line");
  const double arc=mesh.arc[segment]+(context.state.coordinate-segment)*(mesh.arc[segment+1]-mesh.arc[segment]);
  if(root.value.sign*(arc-root.value.arc)<=0) {
    // A boundary condition is applied before reading downstream transport
    // coefficients or drawing a stochastic path. This provider supplies no
    // downstream volume; an already returned record cannot advance in ambient
    // fields for even one numerical interval before being removed.
    RecordTerminal(context,Terminal::ShockReturn,momentum);
    PIC::ParticleBuffer::DeleteParticle(context.pointer);*removed=true;
  }
  return TStatus::Ok();
}
Transport::Status PrepareBoundary(const Transport::PICAdapter::ParticleContext& context,
    double elapsed,double* interval,BoundaryInterval* output,Transport::BoundaryIntervalControl* history) {
  if(!interval||!output||!std::isfinite(*interval)||*interval<0)return Error("invalid numerical boundary interval");
  *output=BoundaryInterval{};
  if(!Enabled()||*interval==0)return TStatus::Ok();
  const auto retained=ReducedShock::Native::CaptureMaterialView();
  if(!retained.ok())return Error(retained.status.message);
  const auto& mesh=*retained.value;
  if(context.state.fieldLineId!=0)return Error("reduced absorbing boundary requested on a foreign material line");
  const auto segment=static_cast<std::size_t>(context.state.coordinate);
  if(segment+1>=mesh.arc.size())return Error("boundary particle coordinate exceeds its immutable material mesh");
  const double begin=EpochAt(context,elapsed);
  const auto initial=RootAt(retained.value,begin);if(!initial.ok())return Error(initial.status.message);
  const double tolerance=Run::Active().get().reducedParticles.boundaryToleranceM;
  // Warm-start only the numerical proposal. Every new interval is checked at
  // its actual explicit times and geometry before a particle/path variate is
  // drawn; past error estimates never authorize an unchecked future interval.
  return Transport::RefineBoundaryInterval(begin,*interval,
      Transport::ActiveNumericalTolerances().minimumStepS,[&](double candidate,bool* resolved) {
    const auto end=RootAt(retained.value,EpochAt(context,elapsed,candidate));if(!end.ok())return Error(end.status.message);
    *resolved=true;
    if(initial.value.active!=end.value.active)return Error(
        "an occupied particle interval crosses reduced front connectivity/criticality creation or loss; event-local topology transport is not qualified by stages 5-6");
    for(double fraction:{0.25,0.5,0.75}) {
      const auto middle=RootAt(retained.value,EpochAt(context,elapsed,fraction*candidate));
      if(!middle.ok())return Error(middle.status.message);
      if(middle.value.active!=initial.value.active||middle.value.sign!=initial.value.sign)
        return Error("reduced shock branch changes inside occupied interval; a joint topology-event operator is required");
      if(initial.value.active&&std::abs(middle.value.arc-(initial.value.arc+fraction*(end.value.arc-initial.value.arc)))>tolerance)
        *resolved=false;
    }
    if(*resolved) {
      output->active=initial.value.active;output->durationS=candidate;
      output->particleArcM=mesh.arc[segment]+(context.state.coordinate-segment)*(mesh.arc[segment+1]-mesh.arc[segment]);
      output->frontBeginArcM=initial.value.arc;output->frontEndArcM=end.value.arc;output->upstreamSign=initial.value.sign;
    }
    return TStatus::Ok();
  },interval,history);
}
Transport::AbsorbingBridgeResult CheckBoundary(const BoundaryInterval& b,double displacement0,
    double displacement1,double elapsed0,double elapsed1,double kdt,Transport::RandomStream& random) {
  Transport::AbsorbingBridgeResult result;result.status=TStatus::Ok();
  if(!b.active)return result;
  if(!(b.durationS>0)||elapsed0<0||elapsed1<elapsed0||elapsed1>b.durationS)
    {result.status=Error("first-passage query exceeds its prepared interval");return result;}
  const double speed=(b.frontEndArcM-b.frontBeginArcM)/b.durationS;
  const double d0=b.upstreamSign*(b.particleArcM+displacement0-b.frontBeginArcM-speed*elapsed0);
  const double d1=b.upstreamSign*(b.particleArcM+displacement1-b.frontBeginArcM-speed*elapsed1);
  return Transport::AbsorbingLinearBridge(d0,d1,kdt,random);
}
} } }

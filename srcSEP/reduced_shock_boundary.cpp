#include "reduced_shock_boundary.h"
#include "reduced_shock_geometry.h"
#include "util/sep_run_configuration.h"
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
    const ReducedShock::SurfaceSphere& sphere) {
  using R=Core::Result<Root>;Root root;
  const auto controls=ReducedShock::Native::Controls(mesh.interval.ArrivalS()-mesh.interval.DepartureS(),mesh.arc[1]);
  const auto height=[&](std::size_t i) {
    return ReducedShock::MoveVertex(mesh.interval,mesh.vertices[i],mesh.interval.DepartureS(),epoch,controls);
  };
  const auto first=height(0),last=height(mesh.vertices.size()-1);
  if(!first.ok())return R::Failure(first.status.code,first.status.message);
  if(!last.ok())return R::Failure(last.status.code,last.status.message);
  const double shock=sphere.centerM[2]+sphere.radiusM;
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
  // This is an exact specialization of the supplied core geometry, not an
  // origin-centred substitute shock: the finite-SSE apex is C_z+a. On the
  // positive solar-rotation axis, the qualified Parker exterior has radial U
  // and B. Its scalar flow map preserves vertex ordering, so a binary search
  // locates the unique leading intersection on the advected material chords.
  // We still call the actual trajectory and RH/front evaluators. No U*dt
  // approximation, cached front time, or copied compression is introduced.
  if(a[0]!=0||a[1]!=0||b[0]!=0||b[1]!=0||!(b[2]>a[2]))
    return R::Failure(Core::StatusCode::InvalidState,"qualified polar material flow lost its axis/order invariant");
  const auto front=mesh.interval.FrontAt({{0,0,shock}},epoch);
  if(!front.ok())return R::Failure(front.status.code,front.status.message);
  if(front.value.supportEdge)return R::Failure(Core::StatusCode::UnsupportedCapability,"polar apex unexpectedly lies on a finite-front support edge");
  if(!front.value.acceptedShock)return R::Success(root);
  root.active=true;root.arc=mesh.arc[left]+(shock-a[2])/(b[2]-a[2])*(mesh.arc[right]-mesh.arc[left]);
  return R::Success(root);
}
Core::Result<Root> EvaluateRootAt(const ReducedShock::Native::MaterialView& mesh,double epoch) {
  using R=Core::Result<Root>;
  const auto sphere=mesh.interval.SphereAt(epoch);
  if(!sphere.ok())return R::Failure(sphere.status.code,sphere.status.message);
  bool polar=sphere.value.centerM[0]==0&&sphere.value.centerM[1]==0&&
      sphere.value.direction[0]==0&&sphere.value.direction[1]==0&&sphere.value.direction[2]==1;
  double previous=0;
  for(const auto& point:mesh.vertices) {
    polar=polar&&point[0]==0&&point[1]==0&&point[2]>previous;previous=point[2];
  }
  if(polar) {
    const auto ambient=mesh.interval.PrimitiveAt(mesh.vertices.front(),epoch);
    if(!ambient.ok())return R::Failure(ambient.status.code,ambient.status.message);
    polar=ambient.value.parkerExterior&&ambient.value.velocityMPerS[0]==0&&ambient.value.velocityMPerS[1]==0;
    if(polar)return PolarRoot(mesh,epoch,sphere.value);
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

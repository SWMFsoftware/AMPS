// Force-aware first/second-order LOCAL PHYSICAL TIME particle tracers.
// See srcMoon/test/M3_FORCE_CUT_CELL_CONTRACT.md. ForceFreeCutSurface remains
// unchanged. Geometry is solved on each accepted parabolic position arc,
// rather than on a chord whose endpoints can both miss an intervening body.
#include "pic.h"
#include "meshAMRcutcell_watertight.h"
#include "pic_mover_force_cut_surface_math.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <sstream>
#include <vector>

namespace PIC { namespace Mover { namespace ForceCutSurface {
cSettings Settings;
fAcceleration Acceleration=nullptr;
namespace {
using Node=cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>;
using Math::State;
struct Record {
  cTimeContext time;
  double currentEpoch;
  unsigned int epochValid,active;
  // Global triangulation index+1 is identical on all ranks and survives
  // opaque MPI cloning. Zero is the freshly allocated/reset "no face" value.
  // A raw face pointer would be invalid in the receiving rank's address space.
  unsigned int lastFaceIndexPlusOne;
};
long int Offset=-1;
Record Read(PIC::ParticleBuffer::byte* data) {
  Record r;std::memcpy(&r,data+Offset,sizeof(r));return r;
}
void Write(PIC::ParticleBuffer::byte* data,const Record& r) {
  std::memcpy(data+Offset,&r,sizeof(r));
}
[[noreturn]] void Error(const char* what,const Record* r=nullptr) {
  std::ostringstream message;
  message<<"force-aware cut-surface tracer: "<<what;
  if(r)message<<"; requested="<<r->time.RequestedTimeS<<" s, elapsed="
      <<r->time.ElapsedTimeS<<" s, remaining="<<r->time.RemainingTimeS<<" s";
  message<<". Check the selected boundary/force/clock callback and explicit "
      "numerical tolerances. Reduce the caller interval or maximum substep; "
      "increase the declared iteration/resource limit only after reviewing cost. "
      "No particle or unused physical interval was silently discarded.";
  exit(__LINE__,__FILE__,message.str().c_str());
  std::abort();
}
bool Finite(const State& s) {
  for(int d=0;d<3;++d)if(!std::isfinite(s.x[d])||!std::isfinite(s.v[d]))return false;
  return true;
}
void Sync(PIC::ParticleBuffer::byte* data,const State& s) {
  PIC::ParticleBuffer::SetX(const_cast<double*>(s.x),data);
  PIC::ParticleBuffer::SetV(const_cast<double*>(s.v),data);
}
bool Within(const double* x,Node* n,double tolerance) {
  if(!n)return false;
  for(int d=0;d<3;++d)if(x[d]<n->xmin[d]-tolerance||x[d]>n->xmax[d]+tolerance)return false;
  return true;
}
// Query just to the requested side of a face, leaving physical x unchanged.
Node* Locate(const double* x,const double* direction,Node* hint,bool across=false) {
  if(!across && Within(x,hint,0))return hint;
  double probe[3];
  for(int d=0;d<3;++d) {
    probe[d]=x[d];
    if(direction && direction[d]!=0) {
      // nextafter(0) is subnormal and is rounded away by mesh coordinate
      // normalization. Probe by a bounded domain-scale rounding distance;
      // the physical position and flight time remain exactly unchanged.
      const double scale=std::max(std::fabs(x[d]),PIC::Mesh::mesh->xGlobalMax[d]-PIC::Mesh::mesh->xGlobalMin[d]);
      probe[d]+=std::copysign(16*std::numeric_limits<double>::epsilon()*scale,direction[d]);
    }
  }
  return PIC::Mesh::mesh->findTreeNode(probe,hint);
}
void Insert(long int ptr,PIC::ParticleBuffer::byte* data,Node* node,const State& state) {
#if _MESH_DIMENSION_ == 3
  if(!node||!node->block)Error("destination block is unavailable; verify allocated MPI halo");
  double probe[3];
  for(int d=0;d<3;++d) {
    const double centre=0.5*(node->xmin[d]+node->xmax[d]);
    probe[d]=state.x[d];
    if(probe[d]<=node->xmin[d]||probe[d]>=node->xmax[d]) {
      const double scale=std::max(std::fabs(probe[d]),node->xmax[d]-node->xmin[d]);
      probe[d]+=std::copysign(16*std::numeric_limits<double>::epsilon()*scale,centre-probe[d]);
    }
  }
  int i,j,k;
  if(PIC::Mesh::mesh->FindCellIndex(probe,i,j,k,node,false)<0)Error("cannot locate destination cell");
  auto block=node->block;
#if _PIC_MOVER__MPI_MULTITHREAD_ == _PIC_MODE_ON_
  // Match the existing temporary-list protocol in the other native movers.
  // Each worker owns a distinct particle; traversal/commit must follow the
  // worker barrier/join, since exchange publishes the head before its next
  // link is filled. Set prev before publication so a later inserter's update
  // of this particle's predecessor cannot be overwritten by this worker.
  PIC::ParticleBuffer::SetPrev(-1,data);
  const long int previous=std::atomic_exchange(
      block->tempParticleMovingListTable+i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k),ptr);
  PIC::ParticleBuffer::SetNext(previous,data);
  if(previous!=-1)PIC::ParticleBuffer::SetPrev(ptr,previous);
#elif _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
  long int& head=block->tempParticleMovingListTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
  PIC::ParticleBuffer::SetPrev(-1,data);PIC::ParticleBuffer::SetNext(head,data);
  if(head!=-1)PIC::ParticleBuffer::SetPrev(ptr,head);head=ptr;
#else
  (void)ptr;(void)data;(void)block;
  Error("this physical-time tracer requires the qualified MPI-only list mode");
#endif
#else
  // DIM<3 build: this is a 3D algorithm.  It uses all three components of
  // node->xmin/xmax (or of the cell indices), which have _MESH_DIMENSION_
  // elements, so in a DIM1/DIM2 build it would read past those arrays.  It is
  // not defined for DIM<3, so stop here instead of reading other memory.
  exit(__LINE__,__FILE__,"Error: PIC::Mover::ForceCutSurface (anonymous namespace)::Insert requires a 3D mesh (_MESH_DIMENSION_==3)");
#endif
}
void ChangeWeightBase(long int ptr,int species,Node* from,Node* to) {
  const double before=from->block->GetLocalParticleWeight(species);
  const double after=to->block->GetLocalParticleWeight(species);
  if(!(before>0)||!(after>0)||!std::isfinite(before)||!std::isfinite(after))Error("invalid local physical particle weight");
  if(before==after)return;
#if _INDIVIDUAL_PARTICLE_WEIGHT_MODE_ == _INDIVIDUAL_PARTICLE_WEIGHT_ON_
  const double correction=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ptr)*before/after;
  if(!(correction>0)||!std::isfinite(correction))Error("local weight conversion overflow");
  PIC::ParticleBuffer::SetIndividualStatWeightCorrection(correction,ptr);
#else
  Error("different local particle weights require individual corrections to conserve N; no W/dt resampling is applied");
#endif
}
struct Force {
  long int ptr;int species;Node* node;Record* record;
  void operator()(const State& s,double elapsed,double* a) {
    if(!Finite(s))Error("nonfinite force evaluation state",record);
    cTimeContext time=record->time;time.ElapsedTimeS=elapsed;
    time.RemainingTimeS=std::max(0.0,time.RequestedTimeS-elapsed);
    Node* evaluationNode=Within(s.x,node,0)?node:Locate(s.x,nullptr,node);
    if(!evaluationNode||!evaluationNode->block)Error("force substage needs an allocated interpolation block; reduce maximum substep",record);
    if(Acceleration)Acceleration(a,species,ptr,s.x,s.v,evaluationNode,time);
    else {
      double x[3],v[3];for(int d=0;d<3;++d){x[d]=s.x[d];v[d]=s.v[d];}
      using PIC::Mover::TotalParticleAcceleration_default;
      _PIC_PARTICLE_MOVER__TOTAL_PARTICLE_ACCELERATION_(a,species,ptr,x,v,evaluationNode);
    }
    for(int d=0;d<3;++d)if(!std::isfinite(a[d]))Error("nonfinite production acceleration",record);
  }
};
enum Kind {Moving=0,Triangle=1,Sphere=2,Domain=3,Block=4,Event=5,End=6};
struct Hit {
  double time;Kind kind=End;CutCell::cTriangleFace* face=nullptr;
  cInternalBoundaryConditionsDescriptor* sphere=nullptr;
  double position[3]{},normal[3]{},surfaceVelocity[3]{};
  Hit(double h):time(h){}
};
bool Earlier(double t,Kind kind,const Hit& hit) {
  return t<hit.time || (t==hit.time && kind<hit.kind);
}
bool InTriangle(double* point,CutCell::cTriangleFace* face) {
  // The adaptive exact 2D orientation signs are shared with the committed
  // watertight line/triangle primitive; closed edges belong to both faces.
  int omit=0;for(int d=1;d<3;++d)if(std::fabs(face->ExternalNormal[d])>std::fabs(face->ExternalNormal[omit]))omit=d;
  const int p=(omit+1)%3,q=(omit+2)%3;
  const double* a=face->x0Face;const double* b=face->x1Face;const double* c=face->x2Face;
  const int s0=CutCell::Watertight::Orient2DSign(b[p]-a[p],b[q]-a[q],point[p]-a[p],point[q]-a[q]);
  const int s1=CutCell::Watertight::Orient2DSign(c[p]-b[p],c[q]-b[q],point[p]-b[p],point[q]-b[q]);
  const int s2=CutCell::Watertight::Orient2DSign(a[p]-c[p],a[q]-c[q],point[p]-c[p],point[q]-c[q]);
  return (s0>=0&&s1>=0&&s2>=0)||(s0<=0&&s1<=0&&s2<=0);
}
double Derivative(const State& s,const double* a,const double* n,double t) {
  double value=0;for(int d=0;d<3;++d)value+=n[d]*(s.v[d]+t*a[d]);return value;
}
double BoxExit(const State& s,const double* a,Node* node,double h) {
  double result=HUGE_VAL;
  for(int d=0;d<3;++d)for(int side=0;side<2;++side) {
    double n[3]={},point[3]={};n[d]=side?1:-1;point[d]=side?node->xmax[d]:node->xmin[d];
    auto roots=Math::PlaneRoots(s,a,n,point,h,Settings.TimeToleranceS);
    if(std::fabs(s.x[d]-point[d])<=Settings.GeometryToleranceM)roots.insert(roots.begin(),0);
    for(auto root:roots) {
      const double t=static_cast<double>(root)*h;
      const double derivative=n[d]*(s.v[d]+t*a[d]);
      const bool out=derivative>0 || (t==0&&derivative==0&&n[d]*a[d]>0);
      if(out&&t<result)result=t;
    }
  }
  return result;
}
Hit Geometry(const State& s,const double* a,Node* node,double h,
    CutCell::cTriangleFace* lastFace,long int ptr,const Record& record) {
  Hit hit(h);
  const double box=BoxExit(s,a,node,h);
  if(box<=h) {hit.time=box;hit.kind=Block;}
  const bool straight=a[0]==0&&a[1]==0&&a[2]==0;
  if(straight) {
    // Reuse the committed mover's exact watertight edge/vertex decisions,
    // gas-side search origin, concave-corner rule and conservative AMR face
    // walk in the zero-force limit. The search may examine neighbour lists
    // but only a hit before this physical leaf exit can be committed here.
    double displacement[3];for(int d=0;d<3;++d)displacement[d]=s.v[d]*h;
    double fraction=0;
    auto face=PIC::Mover::CutSurface::FindFirstCrossing(s.x,displacement,node,lastFace,fraction,s.v);
    const double time=fraction*h;
    if(face&&Earlier(time,Triangle,hit)) {
      hit.time=time;hit.kind=Triangle;hit.face=face;
      double point[3];Math::Position(s,a,time,point);
      CutCell::Watertight::ClosestPointOnTriangle(point,face->x0Face,face->x1Face,face->x2Face,hit.position);
      double correction=0;for(int d=0;d<3;++d)correction+=(hit.position[d]-point[d])*(hit.position[d]-point[d]);
      if(std::sqrt(correction)>Settings.GeometryToleranceM)Error("watertight hit correction exceeds declared geometry tolerance",&record);
    }
  }
  for(auto descriptor=straight?nullptr:node->FirstTriangleCutFace;descriptor;descriptor=descriptor->next) {
    auto face=descriptor->TriangleFace;
    auto roots=Math::PlaneRoots(s,a,face->ExternalNormal,face->x0Face,h,Settings.TimeToleranceS);
    double signedDistance=0;for(int d=0;d<3;++d)signedDistance+=(s.x[d]-face->x0Face[d])*face->ExternalNormal[d];
    if(std::fabs(signedDistance)<=Settings.GeometryToleranceM)roots.insert(roots.begin(),0);
    for(auto root:roots) {
      const double t=static_cast<double>(root)*h;
      if(t>hit.time)continue;
      const double dn=Derivative(s,a,face->ExternalNormal,t);
      double an=0;for(int d=0;d<3;++d)an+=a[d]*face->ExternalNormal[d];
      if(!(dn<0 || (t==0&&dn==0&&an<0)))continue;
      if(lastFace&&face!=lastFace&&t<=Settings.TimeToleranceS) {
        double front=-HUGE_VAL;const double* vertices[3]={face->x0Face,face->x1Face,face->x2Face};
        for(auto vertex:vertices) {double value=0;for(int d=0;d<3;++d)value+=(vertex[d]-lastFace->x0Face[d])*lastFace->ExternalNormal[d];front=std::max(front,value);}
        if(!(front>Settings.GeometryToleranceM))continue;
      }
      double point[3];Math::Position(s,a,t,point);
      double distance=0;for(int d=0;d<3;++d)distance+=(point[d]-face->x0Face[d])*face->ExternalNormal[d];
      for(int d=0;d<3;++d)point[d]-=distance*face->ExternalNormal[d];
      if(!InTriangle(point,face)) {
        double nearest[3];CutCell::Watertight::ClosestPointOnTriangle(point,face->x0Face,face->x1Face,face->x2Face,nearest);
        double norm=0;for(int d=0;d<3;++d)norm+=(nearest[d]-point[d])*(nearest[d]-point[d]);
        if(std::sqrt(norm)>Settings.GeometryToleranceM)continue;
        for(int d=0;d<3;++d)point[d]=nearest[d];
      }
      bool choose=Earlier(t,Triangle,hit);
      if(t==hit.time&&hit.kind==Triangle) {
        const double previous=Derivative(s,a,hit.face->ExternalNormal,t);
        choose=dn<previous || (dn==previous&&face<hit.face);
      }
      if(choose) {hit.time=t;hit.kind=Triangle;hit.face=face;for(int d=0;d<3;++d)hit.position[d]=point[d];}
    }
  }
  for(auto descriptor=node->InternalBoundaryDescriptorList;descriptor;descriptor=descriptor->nextInternalBCelement) {
    if(descriptor->BondaryType==_INTERNAL_BOUNDARY_TYPE_NASTRAN_SURFACE_)continue;
    if(descriptor->BondaryType!=_INTERNAL_BOUNDARY_TYPE_SPHERE_)
      Error("unsupported analytic internal-boundary type; use a qualified sphere or triangulated surface",&record);
    auto sphere=static_cast<cInternalSphericalData*>(descriptor->BoundaryElement);
    // A NULL logical sphere is a sampling/inventory surface, not a physical
    // boundary (LOLA uses its moving terrain and a separate safety sphere).
    if(!sphere->ParticleSphereInteraction)continue;
    double* centre;double radius;sphere->GetSphereGeometricalParameters(centre,radius);
    double startRadius2=0;for(int d=0;d<3;++d)startRadius2+=(s.x[d]-centre[d])*(s.x[d]-centre[d]);
    if(!(radius>0)||!std::isfinite(radius)||std::sqrt(startRadius2)<radius-Settings.GeometryToleranceM)
      Error("invalid analytic sphere or particle starts inside its solid; review source placement",&record);
    auto roots=Math::SphereRoots(s,a,centre,radius,h,Settings.TimeToleranceS);
    for(auto root:roots) {
      const double t=static_cast<double>(root)*h;
      if(!Earlier(t,Sphere,hit))continue;
      double point[3],dot=0;Math::Position(s,a,t,point);
      for(int d=0;d<3;++d)dot+=(point[d]-centre[d])*(s.v[d]+t*a[d]);
      double second=0;for(int d=0;d<3;++d)second+=s.v[d]*s.v[d]+(point[d]-centre[d])*a[d];
      if(!(dot<0 || (t==0&&dot==0&&second<0)))continue;
      hit.time=t;hit.kind=Sphere;hit.sphere=descriptor;
    }
  }
  if(PIC::Mover::FindMovingSurfaceIntersection && hit.time>0) {
    double time=-1,position[3],normal[3],velocity[3];
    if(PIC::Mover::FindMovingSurfaceIntersection(ptr,s.x,s.v,a,
        record.time.ElapsedTimeS,hit.time,time,position,normal,velocity)) {
      if(!std::isfinite(time)||time<0||time>hit.time+Settings.TimeToleranceS)Error("moving-surface callback returned invalid time",&record);
      if(Earlier(time,Moving,hit)) {
        hit.time=time;hit.kind=Moving;
        for(int d=0;d<3;++d){hit.position[d]=position[d];hit.normal[d]=normal[d];hit.surfaceVelocity[d]=velocity[d];}
      }
    }
  }
  return hit;
}
bool EnteringVelocity(const State& state,const Hit& hit) {
  double normal[3]={},surfaceVelocity[3]={};
  if(hit.kind==Triangle)for(int d=0;d<3;++d)normal[d]=hit.face->ExternalNormal[d];
  else if(hit.kind==Moving)for(int d=0;d<3;++d){normal[d]=hit.normal[d];surfaceVelocity[d]=hit.surfaceVelocity[d];}
  else if(hit.kind==Sphere) {
    auto sphere=static_cast<cInternalSphericalData*>(hit.sphere->BoundaryElement);
    double* centre=nullptr,radius=0;sphere->GetSphereGeometricalParameters(centre,radius);
    for(int d=0;d<3;++d)normal[d]=(state.x[d]-centre[d])/radius;
  }
  else return true;
  double relative[3];for(int d=0;d<3;++d)relative[d]=state.v[d]-surfaceVelocity[d];
  return Math::IncomingVelocity(relative,normal);
}
int Trace(long int ptr,double interval,Node* node,int order) {
  ValidateConfiguration();
  if(!IsStorageRegistered())Error("trace storage must be registered before particle-buffer initialization");
  if(!std::isfinite(interval)||interval<0)Error("caller physical interval must be finite and nonnegative");
  auto data=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  Record record=Read(data);
  if(!record.active) {
    record.time=cTimeContext();record.time.InitialEpochS=record.epochValid?record.currentEpoch:NAN;
    record.time.RequestedTimeS=record.time.RemainingTimeS=interval;
    record.time.IntegrationOrder=order;record.active=interval>0;
  }
  else if(record.time.IntegrationOrder!=order || std::fabs(interval-record.time.RemainingTimeS)>Settings.TimeToleranceS)
    Error("resumed trace must retain its original order and physical remaining interval",&record);
  State state;PIC::ParticleBuffer::GetX(state.x,data);PIC::ParticleBuffer::GetV(state.v,data);
  if(!Finite(state)||!node||!node->block)Error("invalid initial phase space or owner block",&record);
  if(node->Thread!=PIC::ThisThread)Error("a trace must begin on its MPI owner",&record);
  if(record.lastFaceIndexPlusOne>static_cast<unsigned int>(CutCell::nBoundaryTriangleFaces))Error("corrupt continued cut-face identity",&record);
  CutCell::cTriangleFace* lastFace=record.lastFaceIndexPlusOne?
      CutCell::BoundaryTriangleFaces+record.lastFaceIndexPlusOne-1:nullptr;
  long int segments=0,interactions=0;int samePoint=0;
  const int eventCount=(QueryTrajectoryEvent!=nullptr)+(AdvanceTrajectoryEvent!=nullptr)+(ProcessTrajectoryEvent!=nullptr);
  if(eventCount!=0&&eventCount!=3)Error("trajectory event requires its complete three-callback set");
  while(record.time.RemainingTimeS>0) {
    if(++segments>Settings.MaximumSubsteps)Error("maximum substep count exhausted",&record);
    const int species=PIC::ParticleBuffer::GetI(data);
    if(species<0||species>=PIC::nTotalSpecies)Error("invalid particle species",&record);
    Write(data,record);Sync(data,state);
    Force force{ptr,species,node,&record};double acceleration[3];force(state,record.time.ElapsedTimeS,acceleration);
    double h=std::min(record.time.RemainingTimeS,Settings.MaximumSubstepS);
    h=std::min(h,BoxExit(state,acceleration,node,h));
    auto Preview=[&](double maximum) {
      Hit candidate=maximum>0?Geometry(state,acceleration,node,maximum,lastFace,ptr,record):Hit(0);
      if(maximum==0)candidate.kind=Block;
      if(eventCount) {
        double eventTime=-1;
        const int code=QueryTrajectoryEvent(ptr,species,node,node->block->GetLocalTimeStep(species),candidate.time,eventTime);
        if(code!=TrajectoryEventNone&&code!=TrajectoryEventPending)Error("invalid trajectory-event query result",&record);
        if(code==TrajectoryEventPending) {
          if(!std::isfinite(eventTime)||eventTime<0||eventTime>candidate.time+Settings.TimeToleranceS)Error("invalid trajectory-event time",&record);
          if(eventTime<candidate.time || (eventTime==candidate.time&&candidate.kind==End)) {candidate.time=eventTime;candidate.kind=Event;}
        }
      }
      return candidate;
    };
    // Clip BEFORE evaluating any midpoint/substage acceleration. A stage on
    // the original longer arc could otherwise enter a solid, or use the old
    // species beyond a due reaction. Refinement repeats the complete preview.
    Hit hit=Preview(h);h=hit.time;
    if(h>0) {
      bool accepted=false;
      for(int refinement=0;refinement<=Settings.MaximumRefinements;++refinement) {
        if(refinement!=0){hit=Preview(h);h=hit.time;}
        if(h==0){accepted=true;break;}
        const State whole=Math::Step(state,record.time.ElapsedTimeS,h,order,force);
        if(!Finite(whole))Error("nonfinite trial state before physical event",&record);
        // The quadratic position extension and midpoint velocity agree only
        // to their declared numerical order. Near a grazing accelerated hit,
        // never reflect an integrated velocity that is ALREADY leaving the
        // surface back into its solid. No clock, state or ledger is committed
        // until this consistency check and the force error criterion pass.
        if(!EnteringVelocity(whole,hit)) {
          if(!Settings.Adaptive)Error("quadratic arc predicts entry but integrated event velocity is leaving; enable adaptive mode or reduce maximum substep",&record);
          h*=.5;
          if(h<Settings.MinimumSubstepS)Error("grazing impact consistency reached minimum substep",&record);
          continue;
        }
        if(!Settings.Adaptive){accepted=true;break;}
        const State middle=Math::Step(state,record.time.ElapsedTimeS,h/2,order,force);
        double middleAcceleration[3];force(middle,record.time.ElapsedTimeS+h/2,middleAcceleration);
        // The two-half-step estimator has its OWN second parabolic arc.
        // Reject/refine if that trial arc hits a surface/leaf before its
        // endpoint, rather than sampling a midpoint inside the obstacle.
        Record middleRecord=record;middleRecord.time.ElapsedTimeS+=h/2;
        const auto middleHit=Geometry(middle,middleAcceleration,node,h/2,nullptr,ptr,middleRecord);
        if(middleHit.time<h/2-Settings.TimeToleranceS) {
          h*=.5;
          if(h<Settings.MinimumSubstepS)Error("surface-aware force error control reached minimum substep",&record);
          continue;
        }
        const State fine=Math::Step(middle,record.time.ElapsedTimeS+h/2,h/2,order,force,middleAcceleration);
        double error=0;
        for(int d=0;d<3;++d) {
          const double xs=Settings.AbsolutePositionToleranceM+Settings.RelativeTolerance*std::max(std::fabs(state.x[d]),std::fabs(fine.x[d]));
          const double vs=Settings.AbsoluteVelocityToleranceMPerS+Settings.RelativeTolerance*std::max(std::fabs(state.v[d]),std::fabs(fine.v[d]));
          const double divisor=order==1?0.5:0.75;
          error=std::max(error,std::fabs(whole.x[d]-fine.x[d])/(divisor*xs));
          error=std::max(error,std::fabs(whole.v[d]-fine.v[d])/(divisor*vs));
        }
        if(error<=1){accepted=true;break;}
        h*=0.5;
        if(h<Settings.MinimumSubstepS)Error("force error control reached minimum substep",&record);
      }
      if(!accepted)Error("force error control exhausted refinements",&record);
    }
    const double dt=hit.time;
    if(dt>0) {
      state=Math::Step(state,record.time.ElapsedTimeS,dt,order,force);
      if(eventCount)AdvanceTrajectoryEvent(ptr,species,node,node->block->GetLocalTimeStep(species),dt);
      record.time.ElapsedTimeS+=dt;record.time.RemainingTimeS-=dt;
      if(record.time.RemainingTimeS<0 && record.time.RemainingTimeS>=-Settings.TimeToleranceS)record.time.RemainingTimeS=0;
      if(record.time.RemainingTimeS<0)Error("physical clock advanced beyond requested interval",&record);
      samePoint=0;
    }
    else if(++samePoint>Settings.MaximumSamePointInteractions) {
      std::ostringstream diagnostic;diagnostic.precision(17);
      diagnostic<<"unresolved zero-time boundary/event cycle; ptr="<<ptr<<" kind="<<hit.kind
          <<" last_face="<<(lastFace?lastFace-CutCell::BoundaryTriangleFaces:-1)<<" x=("
          <<state.x[0]<<","<<state.x[1]<<","<<state.x[2]<<") v=("<<state.v[0]<<","<<state.v[1]<<","<<state.v[2]
          <<") leaf=["<<node->xmin[0]<<","<<node->xmax[0]<<"]x["<<node->xmin[1]<<","<<node->xmax[1]<<"]x["<<node->xmin[2]<<","<<node->xmax[2]<<"]";
      Error(diagnostic.str().c_str(),&record);
    }
    // A zero-time AMR hand-off is NOT a new geometric trajectory. Preserve
    // convex/concave same-point history through local/MPI ownership changes;
    // clear it only after actual positive flight (or a fresh-record reset).
    if(dt>0){lastFace=nullptr;record.lastFaceIndexPlusOne=0;}
    if(hit.kind==Moving||hit.kind==Triangle)for(int d=0;d<3;++d)state.x[d]=hit.position[d];
    if(hit.kind==Block) {
      // Root arithmetic can leave x=+5e-20 at an exact x=0 crossing.
      // A directional ownership probe then selects the negative-x leaf,
      // while strict containment would immediately select positive-x again.
      // Project ONLY outward crossed coordinates onto their exact shared
      // box plane, within the declared geometry tolerance. Never displace
      // an arbitrary interior coordinate or change elapsed physical time.
      for(int d=0;d<3;++d) {
        const double direction=state.v[d]+Settings.TimeToleranceS*acceleration[d];
        if(direction<0&&std::fabs(state.x[d]-node->xmin[d])<=Settings.GeometryToleranceM)state.x[d]=node->xmin[d];
        if(direction>0&&std::fabs(state.x[d]-node->xmax[d])<=Settings.GeometryToleranceM)state.x[d]=node->xmax[d];
      }
    }
    Sync(data,state);Write(data,record);
    bool internalEndpointEvent=false;
    if(hit.kind==Block && record.time.RemainingTimeS==0 && eventCount==3) {
      double direction[3];for(int d=0;d<3;++d)direction[d]=state.v[d]+Settings.TimeToleranceS*acceleration[d];
      Node* neighbour=Locate(state.x,direction,node,true);
      // AMR ownership is administrative, so an exact endpoint cannot hide
      // an already due physical reaction. External boundaries and material
      // surfaces keep their declared priority; no epsilon advances an event.
      if(neighbour && neighbour!=node) {
        double due=-1;const int code=QueryTrajectoryEvent(ptr,species,node,node->block->GetLocalTimeStep(species),0,due);
        if(code!=TrajectoryEventNone&&code!=TrajectoryEventPending)Error("invalid internal-endpoint event query result",&record);
        if(code==TrajectoryEventPending) {
          if(due!=0)Error("zero-horizon endpoint query must return an exactly due event",&record);
          internalEndpointEvent=true;hit.kind=Event;
        }
      }
    }
    if(hit.kind==Event) {
      const double remaining=record.time.RemainingTimeS;
      double callbackRemaining=remaining;int changedSpecies=species;
      const int result=ProcessTrajectoryEvent(ptr,state.x,state.v,changedSpecies,node,
          node->block->GetLocalTimeStep(species),callbackRemaining);
      if(!PIC::ParticleBuffer::IsParticleAllocated(data))Error("event callback deleted mover-owned current particle",&record);
      if(result==TrajectoryEventDelete){PIC::ParticleBuffer::DeleteParticle(ptr);return _PARTICLE_LEFT_THE_DOMAIN_;}
      if(result!=TrajectoryEventContinue || !Finite(state) ||
          changedSpecies!=PIC::ParticleBuffer::GetI(data) ||
          !std::isfinite(callbackRemaining)||std::fabs(callbackRemaining-remaining)>Settings.TimeToleranceS)
        Error("event callback must preserve physical remaining seconds and a valid current particle",&record);
      Sync(data,state);
      double again=-1;
      if(QueryTrajectoryEvent(ptr,changedSpecies,node,node->block->GetLocalTimeStep(changedSpecies),0,again)!=TrajectoryEventNone)
        Error("event callback did not consume its due event",&record);
      if(internalEndpointEvent)hit.kind=Block;
    }
    if(hit.kind==Triangle||hit.kind==Moving||hit.kind==Sphere) {
      if(++interactions>Settings.MaximumInteractions)Error("surface interaction limit exhausted",&record);
      int result=_PARTICLE_DELETED_ON_THE_FACE_;
      if(hit.kind==Triangle) {
        if(ProcessTriangleCutFaceIntersection)result=ProcessTriangleCutFaceIntersection(ptr,state.x,state.v,hit.face,node);
        lastFace=hit.face;
        record.lastFaceIndexPlusOne=static_cast<unsigned int>(hit.face-CutCell::BoundaryTriangleFaces)+1;
      }
      else if(hit.kind==Moving) {
        if(!ProcessMovingSurfaceIntersection)Error("moving intersection has no interaction callback");
        result=ProcessMovingSurfaceIntersection(ptr,state.x,state.v,hit.normal,hit.surfaceVelocity,node);
      }
      else {
        auto sphere=static_cast<cInternalSphericalData*>(hit.sphere->BoundaryElement);
        double remaining=record.time.RemainingTimeS;
        result=sphere->ParticleSphereInteraction(species,ptr,state.x,state.v,remaining,node,sphere);
        if(std::fabs(remaining-record.time.RemainingTimeS)>Settings.TimeToleranceS)Error("sphere callback changed physical tracing interval",&record);
      }
      if(result==_PARTICLE_DELETED_ON_THE_FACE_){PIC::ParticleBuffer::DeleteParticle(ptr);return _PARTICLE_LEFT_THE_DOMAIN_;}
      if(!Finite(state))Error("surface callback returned nonfinite state",&record);
    }
    else if(hit.kind==Block) {
      double direction[3];for(int d=0;d<3;++d)direction[d]=state.v[d]+Settings.TimeToleranceS*acceleration[d];
      Node* next=Locate(state.x,direction,node,true);
      if(!next) {
        bool crossed[3]={};int face=-1;
        for(int d=0;d<3;++d) {
          if(direction[d]>0&&state.x[d]>=PIC::Mesh::mesh->xGlobalMax[d]-Settings.GeometryToleranceM){crossed[d]=true;if(face<0)face=2*d+1;}
          if(direction[d]<0&&state.x[d]<=PIC::Mesh::mesh->xGlobalMin[d]+Settings.GeometryToleranceM){crossed[d]=true;if(face<0)face=2*d;}
        }
        if(face<0)Error("AMR face traversal failed inside computational domain",&record);
#if _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE_ == _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__DELETE_
        PIC::ParticleEvents::RecordDomainExit(ptr,node,state.x,state.v);
        PIC::ParticleBuffer::DeleteParticle(ptr);return _PARTICLE_LEFT_THE_DOMAIN_;
#elif _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE_ == _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__SPECULAR_REFLECTION_
        for(int d=0;d<3;++d)if(crossed[d])state.v[d]=-state.v[d];
#elif _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE_ == _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__USER_FUNCTION_
        if(!ProcessOutsideDomainParticles)Error("domain user callback is unavailable");
        if(ProcessOutsideDomainParticles(ptr,state.x,state.v,face,node)==_PARTICLE_DELETED_ON_THE_FACE_){PIC::ParticleBuffer::DeleteParticle(ptr);return _PARTICLE_LEFT_THE_DOMAIN_;}
#else
        Error("unsupported domain boundary policy");
#endif
      }
      else {
        if(next==node && !internalEndpointEvent)Error("AMR face probe did not leave current block",&record);
        if(!next->block)Error("adjacent destination MPI ghost block is not allocated",&record);
        // A Continue event may change species or reverse velocity exactly
        // here. Use its actual new species/base weight and allow it to stay
        // in the old closed leaf when its new direction points back inside.
        ChangeWeightBase(ptr,PIC::ParticleBuffer::GetI(data),node,next);node=next;
        if(node->Thread!=PIC::ThisThread) {
          if(record.time.RemainingTimeS==0) {
            record.active=0;
            if(record.epochValid)record.currentEpoch=record.time.InitialEpochS+record.time.ElapsedTimeS;
          }
          Sync(data,state);Write(data,record);Insert(ptr,data,node,state);
          return _PARTICLE_IN_NOT_IN_USE_NODE_;
        }
      }
    }
    if(!Finite(state))Error("nonfinite advanced state",&record);
    // A surface policy may apply a documented outward rounding offset. Use
    // its actual new leaf and conserve represented N before another force
    // evaluation; do not keep interpolating from the pre-impact leaf.
    if(!Within(state.x,node,0)) {
      Node* next=Locate(state.x,nullptr,node);
      if(!next||!next->block)Error("surface/event callback left the allocated physical domain",&record);
      ChangeWeightBase(ptr,PIC::ParticleBuffer::GetI(data),node,next);node=next;
      if(node->Thread!=PIC::ThisThread) {
        if(record.time.RemainingTimeS==0){record.active=0;if(record.epochValid)record.currentEpoch=record.time.InitialEpochS+record.time.ElapsedTimeS;}
        Sync(data,state);Write(data,record);Insert(ptr,data,node,state);return _PARTICLE_IN_NOT_IN_USE_NODE_;
      }
    }
  }
  record.time.ElapsedTimeS=record.time.RequestedTimeS;record.time.RemainingTimeS=0;record.active=0;
  if(record.epochValid)record.currentEpoch=record.time.InitialEpochS+record.time.ElapsedTimeS;
  Write(data,record);Sync(data,state);Insert(ptr,data,node,state);
  return _PARTICLE_MOTION_FINISHED_;
}
} // implementation details

bool IsStorageRegistered(){return Offset>=0;}
void RegisterStorage() {
  if(IsStorageRegistered())return;
  if(PIC::ParticleBuffer::ParticleDataBuffer)Error("cannot register tracing fields after particle-buffer initialization");
  ValidateConfiguration();PIC::ParticleBuffer::RequestDataStorage(Offset,sizeof(Record),alignof(Record));
}
void ValidateConfiguration() {
#if DIM != 3
  Error("force-aware cut-surface tracing requires DIM=3");
#endif
#if _COMPILATION_MODE_ != _COMPILATION_MODE__MPI_
  Error("force-aware cut-surface tracing requires qualified MPI-only execution");
#endif
#if _PIC_LOCAL_TIME_STEP_FLUX_RESAMPLING_MODE_ == _PIC_MODE_ON_
  Error("physical-seconds tracing cannot silently use W/dt flux resampling");
#endif
#if _CUDA_MODE_ == _ON_
  Error("force-aware cut-surface tracing is not implemented for CUDA");
#endif
#if _PIC_FIELD_LINE_MODE_ == _PIC_MODE_ON_
  Error("force-aware cut-surface tracing requires spatial particle lists, not field-line mode");
#endif
#if _PIC_BC__PERIODIC_MODE_ == _PIC_BC__PERIODIC_MODE_ON_
  Error("periodic wrapping is not implemented by this tracer; use a qualified periodic mover or an explicitly nonperiodic domain");
#endif
#if _PIC_PARTICLE_MOVER__BACKWARD_TIME_INTEGRATION_MODE_ != _PIC_PARTICLE_MOVER__BACKWARD_TIME_INTEGRATION_MODE__DISABLED_
  Error("backward integration is not implemented by this forward physical-time/event clock; use a qualified backward mover");
#endif
#if _PIC_FIELD_SOLVER_MODE_ != _PIC_FIELD_SOLVER_MODE__OFF_
  Error("self-consistent PIC field solving/current-charge deposition is not implemented by this tracer; disable the field solver for prescribed Lorentz forcing or select a qualified depositing mover");
#endif
  const double positive[]={Settings.RelativeTolerance,Settings.AbsolutePositionToleranceM,
      Settings.AbsoluteVelocityToleranceMPerS,Settings.GeometryToleranceM,Settings.TimeToleranceS,Settings.MinimumSubstepS};
  for(double x:positive)if(!(x>0)||!std::isfinite(x))Error("numerical tolerances must be finite and positive");
  if(!(Settings.MaximumSubstepS>0)||std::isnan(Settings.MaximumSubstepS)||
      Settings.MaximumSubstepS<Settings.MinimumSubstepS||Settings.MaximumSubsteps<1||
      Settings.MaximumRefinements<1||Settings.MaximumInteractions<1||Settings.MaximumSamePointInteractions<1)
    Error("invalid segment/iteration controls");
}
void ResetParticleState(PIC::ParticleBuffer::byte* data) {
  if(!IsStorageRegistered())return;
  if(!data)Error("cannot reset a null particle record");
  std::memset(data+Offset,0,sizeof(Record));
}
void SetParticleEpochS(long int ptr,double epoch) {
  if(!IsStorageRegistered()||!std::isfinite(epoch))Error("set epoch requires registered storage and finite physical epoch");
  auto data=PIC::ParticleBuffer::GetParticleDataPointer(ptr);auto record=Read(data);
  if(record.active)Error("cannot reset an active physical continuation epoch",&record);
  record.currentEpoch=epoch;record.epochValid=1;record.time=cTimeContext();record.time.InitialEpochS=epoch;Write(data,record);
}
double GetParticleEpochS(long int ptr) {
  if(!IsStorageRegistered())return NAN;
  const auto record=Read(PIC::ParticleBuffer::GetParticleDataPointer(ptr));
  return !record.epochValid?NAN:record.active?record.time.InitialEpochS+record.time.ElapsedTimeS:record.currentEpoch;
}
cTimeContext GetTimeContext(long int ptr) {
  if(!IsStorageRegistered())return cTimeContext();
  auto record=Read(PIC::ParticleBuffer::GetParticleDataPointer(ptr));
  if(!record.epochValid)record.time.InitialEpochS=NAN;
  return record.time;
}
void QueueDetachedParticle(long int ptr,double dt,Node* node,int order) {
  if(!IsStorageRegistered()||!node||!node->block||node->Thread!=PIC::ThisThread||
      !std::isfinite(dt)||dt<0||(order!=1&&order!=2))Error("invalid detached sibling continuation request");
  auto data=PIC::ParticleBuffer::GetParticleDataPointer(ptr);Record record=Read(data);
  State state;PIC::ParticleBuffer::GetX(state.x,data);PIC::ParticleBuffer::GetV(state.v,data);
  if(!PIC::ParticleBuffer::IsParticleAllocated(data)||record.active||!record.epochValid||record.lastFaceIndexPlusOne!=0||
      PIC::ParticleBuffer::GetNext(data)!=-1||PIC::ParticleBuffer::GetPrev(data)!=-1||
      !Finite(state)||!Within(state.x,node,Settings.GeometryToleranceM))Error("sibling must be allocated, unlinked, reset and assigned its event epoch");
  record.time=cTimeContext();record.time.InitialEpochS=record.currentEpoch;
  record.time.RequestedTimeS=record.time.RemainingTimeS=dt;record.time.IntegrationOrder=order;record.active=dt>0;
  Write(data,record);Insert(ptr,data,node,state);
}
long int CountPending() {
  if(!IsStorageRegistered())return 0;
  long int count=0;const int cells=_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;
  for(int thread=0;thread<PIC::Mesh::mesh->nTotalThreads;++thread)
    for(auto node=thread==PIC::ThisThread?PIC::Mesh::mesh->ParallelNodesDistributionList[thread]:PIC::Mesh::mesh->DomainBoundaryLayerNodesList[thread];node;node=node->nextNodeThisThread) {
      if(!node->block)continue;
      for(int cell=0;cell<cells;++cell)for(long int ptr=node->block->FirstCellParticleTable[cell];ptr!=-1;ptr=PIC::ParticleBuffer::GetNext(ptr)) {
        const auto r=Read(PIC::ParticleBuffer::GetParticleDataPointer(ptr));
        if(r.active) {if(!(r.time.RemainingTimeS>0)||!std::isfinite(r.time.RemainingTimeS))Error("corrupt pending physical clock",&r);++count;}
      }
    }
  return count;
}
long int MoveContinuations() {
  if(!IsStorageRegistered())return 0;
  long int count=0;const int cells=_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;
  for(unsigned int index=0;index<PIC::DomainBlockDecomposition::nLocalBlocks;++index) {
    auto node=PIC::DomainBlockDecomposition::BlockTable[index];if(!node->block)continue;
    for(int cell=0;cell<cells;++cell) {
      long int& head=node->block->FirstCellParticleTable[cell];long int ptr=head;
      while(ptr!=-1) {
        const long int next=PIC::ParticleBuffer::GetNext(ptr);const auto r=Read(PIC::ParticleBuffer::GetParticleDataPointer(ptr));
        if(r.active) {PIC::ParticleBuffer::ExcludeParticleFromList(ptr,head);Trace(ptr,r.time.RemainingTimeS,node,r.time.IntegrationOrder);++count;}
        ptr=next;
      }
    }
  }
  PIC::Mover::CommitTempParticleMovingListsToFirstCellParticleTable();return count;
}
}}}

int PIC::Mover::ForceCutSurfaceFirstOrder(long int ptr,double dt,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  return ForceCutSurface::Trace(ptr,dt,node,1);
}
int PIC::Mover::ForceCutSurfaceSecondOrder(long int ptr,double dt,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  return ForceCutSurface::Trace(ptr,dt,node,2);
}

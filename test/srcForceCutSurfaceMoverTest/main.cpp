// Native, analytical verification of the OPTIONAL production force tracers.
// These deliberately dimensioned metre/second fields are numerical fixtures,
// not observations and not alternative Moon force implementations.  Every
// trajectory, collision, list insertion and MPI hand-off uses the core mover.
#include "pic.h"
#include "AnalyticPolyhedra.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

namespace {
using Node=cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>;
namespace Trace=PIC::Mover::ForceCutSurface;
struct EventData {double remaining;int active;};
long int EventOffset=-1;
int Field=0,Order=1,Failures=0,Checks=0,Absorb=0;
long int TriangleHits=0,SphereHits=0,MovingHits=0,EventHits=0;
double HitTime=0,EventPosition=0;
bool WithSibling=false;
bool WithMoving=false;
bool DeleteCurrent=false;
FILE* Report=nullptr;
std::vector<Node*> Leaves;

void Collect(Node* node) {
  if(node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_){Leaves.push_back(node);return;}
  for(int i=0;i<8;++i)if(node->downNode[i])Collect(node->downNode[i]);
}
double Resolution(double*){return 1.01*std::sqrt(3.0)*0.25;}
double SurfaceResolution(CutCell::cTriangleFace*){return Resolution(nullptr);}
double UnitLoad(Node*){return 1;}
void Check(bool condition,const char* name,double value=0,double expected=0) {
  int local=condition?1:0,all=0;
  MPI_Allreduce(&local,&all,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  ++Checks;if(!all)++Failures;
  if(PIC::ThisThread==0){fprintf(Report,"%s order=%d %s value=%.17g expected=%.17g\n",all?"PASS":"FAIL",Order,name,value,expected);fflush(Report);}
}
double Sum(double value) {
  double result=0;MPI_Allreduce(&value,&result,1,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);return result;
}
EventData GetEvent(long int ptr){EventData result;std::memcpy(&result,PIC::ParticleBuffer::GetParticleDataPointer(ptr)+EventOffset,sizeof(result));return result;}
void SetEvent(long int ptr,EventData event){std::memcpy(PIC::ParticleBuffer::GetParticleDataPointer(ptr)+EventOffset,&event,sizeof(event));}
void Force(double* a,int,long int,const double* x,const double* v,Node*,const Trace::cTimeContext& time) {
  // Negative control: this callback REFUSES to supply a force inside either
  // obstacle. A passing collision test therefore proves the actual numerical
  // midpoint/error-estimator never sampled beyond the first physical impact.
  const bool inFloor=x[0]>-1.5+1e-10&&x[0]<1.5-1e-10&&x[1]>-1.5+1e-10&&x[1]<1.5-1e-10&&x[2]>-.9+1e-10&&x[2]<-.8-1e-10;
  const double sphereDistance=std::sqrt((x[0]-1)*(x[0]-1)+x[1]*x[1]+(x[2]-.5)*(x[2]-.5));
  const bool underMoving=WithMoving&&x[2]<.3+.1*time.ElapsedTimeS-1e-10;
  if(inFloor||sphereDistance<.1-1e-10||underMoving)exit(__LINE__,__FILE__,"verification force refuses illegal substage inside obstacle");
  a[0]=a[1]=a[2]=0;
  if(Field==1){a[0]=0.3;a[1]=-0.2;}
  if(Field==2)a[2]=2;
  if(Field==3)a[2]=-1;
  // A unit cyclotron/Coriolis frequency is an analytical, velocity-dependent
  // control of the actual callback/midpoint path, including transverse drift.
  if(Field==4){a[0]=v[1];a[1]=-v[0];}
  if(Field==5)a[0]=time.ElapsedTimeS;
  if(Field==6)a[0]=1;
  if(Field==7)a[2]=50*time.ElapsedTimeS;
}
int Triangle(long int ptr,double*,double* v,CutCell::cTriangleFace* face,Node*) {
  ++TriangleHits;HitTime=Trace::GetTimeContext(ptr).ElapsedTimeS;
  if(Absorb)return _PARTICLE_DELETED_ON_THE_FACE_;
  double dot=0;for(int d=0;d<3;++d)dot+=v[d]*face->ExternalNormal[d];
  for(int d=0;d<3;++d)v[d]-=2*dot*face->ExternalNormal[d];
  return _PARTICLE_REJECTED_ON_THE_FACE_;
}
int Sphere(int,long int,double* x,double* v,double&,void*,void* descriptor) {
  ++SphereHits;auto sphere=static_cast<cInternalSphericalData*>(descriptor);
  double* centre=nullptr,radius=0;sphere->GetSphereGeometricalParameters(centre,radius);
  double normal[3],dot=0;for(int d=0;d<3;++d){normal[d]=(x[d]-centre[d])/radius;dot+=v[d]*normal[d];}
  for(int d=0;d<3;++d)v[d]-=2*dot*normal[d];
  return _PARTICLE_REJECTED_ON_THE_FACE_;
}
bool MovingIntersection(long int,const double* x,const double* v,const double* a,double elapsed,
    double maximum,double& due,double* position,double* normal,double* velocity) {
  // Explicit analytical moving plane z=.3 m + (.1 m/s)*elapsed. This is a
  // test geometry adapter through the REAL production moving-surface API.
  // Root evaluation is independent of the fixture's closed-form final oracle.
  if(!WithMoving||maximum<=0)return false;
  const double gap=x[2]-(.3+.1*elapsed),relative=v[2]-.1;
  double time=HUGE_VAL;
  if(a[2]==0){if(relative<0)time=-gap/relative;}
  else {
    const double discriminant=relative*relative-2*a[2]*gap;
    if(discriminant>=0)for(double candidate:{(-relative-std::sqrt(discriminant))/a[2],(-relative+std::sqrt(discriminant))/a[2]})
      if(candidate>=0&&relative+a[2]*candidate<0)time=std::min(time,candidate);
  }
  if(time<0||time>maximum)return false;
  due=time;for(int d=0;d<3;++d){position[d]=x[d]+time*v[d]+.5*time*time*a[d];normal[d]=d==2?1:0;velocity[d]=d==2?.1:0;}
  return true;
}
int MovingInteraction(long int,double*,double* v,const double* normal,const double* surfaceVelocity,Node*) {
  ++MovingHits;double incoming=0;for(int d=0;d<3;++d)incoming+=(v[d]-surfaceVelocity[d])*normal[d];
  for(int d=0;d<3;++d)v[d]-=2*incoming*normal[d];
  return _PARTICLE_REJECTED_ON_THE_FACE_;
}
int Query(long int ptr,int,Node*,double,double maximum,double& due) {
  const auto event=GetEvent(ptr);
  if(event.active&&event.remaining<=maximum){due=std::max(0.0,event.remaining);return PIC::Mover::TrajectoryEventPending;}
  return PIC::Mover::TrajectoryEventNone;
}
void Advance(long int ptr,int,Node*,double,double dt) {
  auto event=GetEvent(ptr);if(event.active)event.remaining-=dt;SetEvent(ptr,event);
}
int Process(long int ptr,double* x,double* v,int& species,Node* node,double,double&) {
  auto event=GetEvent(ptr);event.active=0;SetEvent(ptr,event);
  ++EventHits;EventPosition=x[0];
  const double represented=node->block->GetLocalParticleWeight(species)*PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ptr);
  species=_NA_PLUS_SPEC_;PIC::ParticleBuffer::SetI(species,ptr);
  PIC::ParticleBuffer::SetIndividualStatWeightCorrection(represented/node->block->GetLocalParticleWeight(species),ptr);
  v[0]=2;
  if(WithSibling) {
    // The callback stages one deliberately independent control tracer.  This
    // is a work-queue/list test, not an atom-conservation reaction fixture.
    const long int sibling=PIC::ParticleBuffer::GetNewParticle();
    auto data=PIC::ParticleBuffer::GetParticleDataPointer(sibling);
    PIC::ParticleBuffer::CloneParticle(data,PIC::ParticleBuffer::GetParticleDataPointer(ptr));
    PIC::ParticleBuffer::SetNext(-1,data);PIC::ParticleBuffer::SetPrev(-1,data);
    PIC::ParticleBuffer::SetX(x,data);double velocity[3]={-1,0,0};PIC::ParticleBuffer::SetV(velocity,data);
    Trace::ResetParticleState(data);Trace::SetParticleEpochS(sibling,Trace::GetParticleEpochS(ptr));
    SetEvent(sibling,{0,0});
    Trace::QueueDetachedParticle(sibling,Trace::GetTimeContext(ptr).RemainingTimeS,node,Order);
  }
  return DeleteCurrent?PIC::Mover::TrajectoryEventDelete:PIC::Mover::TrajectoryEventContinue;
}
void Setup() {
  PIC::Init_BeforeParser();
  PIC::ParticleBuffer::RequestDataStorage(EventOffset,sizeof(EventData),alignof(EventData));
  char report[128];std::snprintf(report,sizeof(report),"force-cut-surface.np%d.report",PIC::nTotalThreads);
  if(PIC::ThisThread==0)Report=fopen(report,"w");
  AnalyticPolyhedra::cScene scene;auto floor=AnalyticPolyhedra::Box(3,3,.1,1);
  double identity[3][3]={{1,0,0},{0,1,0},{0,0,1}},shift[3]={-1.5,-1.5,-.9};floor.Transform(identity,shift);scene.Add(floor);
  if(PIC::ThisThread==0&&!scene.WriteNastran("force-cut-fixture.bdf","analytical force tracer fixture"))exit(__LINE__,__FILE__,"cannot write fixture surface");
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  CutCell::ReadNastranSurfaceMeshLongFormat("force-cut-fixture.bdf");
  PIC::BC::InternalBoundary::NastranSurface::RegisterInternalNastranSurface();
  PIC::Mesh::mesh->CutCellSurfaceLocalResolution=SurfaceResolution;
  SetGeneralSphereSurfaceMeshParameters(8,16);PIC::BC::InternalBoundary::Sphere::Init();
  auto descriptor=PIC::BC::InternalBoundary::Sphere::RegisterInternalSphere();
  auto sphere=static_cast<cInternalSphericalData*>(descriptor.BoundaryElement);
  double centre[3]={1,0,.5};sphere->SetSphereGeometricalParameters(centre,.1);sphere->ParticleSphereInteraction=Sphere;
  sphere->InjectionRate=nullptr;sphere->InjectionBoundaryCondition=nullptr;
  double xmin[3]={-2,-2,-2},xmax[3]={2,2,2};
  PIC::Mesh::mesh->AllowBlockAllocation=false;PIC::Mesh::mesh->init(xmin,xmax,Resolution);
  if(PIC::ThisThread==0){PIC::Mesh::mesh->buildMesh();PIC::Mesh::mesh->saveMeshFile("force-cut-fixture.mesh");}
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  if(PIC::ThisThread!=0)PIC::Mesh::mesh->readMeshFile("force-cut-fixture.mesh");
  PIC::Mesh::mesh->SetParallelLoadMeasure(UnitLoad);PIC::Mesh::mesh->CreateNewParallelDistributionLists();
  PIC::Mesh::initCellSamplingDataBuffer();PIC::Mesh::mesh->AllowBlockAllocation=true;PIC::Mesh::mesh->AllocateTreeBlocks();
  PIC::Mesh::IrregularSurface::InitExternalNormalVector();PIC::Init_AfterParser();PIC::Mover::Init();
  PIC::DomainBlockDecomposition::UpdateBlockTable();Collect(PIC::Mesh::mesh->rootTree);
  for(auto node:Leaves)if(node->block) {
    for(int species=0;species<PIC::nTotalSpecies;++species) {
      // Deliberately unequal neighbouring W and local dt expose accidentally
      // substituting destination dt for the caller's remaining physical time.
      node->block->SetLocalParticleWeight(1+species+.05*(node->xmin[0]+2),species);
      node->block->SetLocalTimeStep(.07+.03*(node->xmin[0]+2)+.02*species,species);
    }
    const double volume=(node->xmax[0]-node->xmin[0])*(node->xmax[1]-node->xmin[1])*(node->xmax[2]-node->xmin[2])/
      (_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_);
    for(int k=0;k<_BLOCK_CELLS_Z_;++k)for(int j=0;j<_BLOCK_CELLS_Y_;++j)for(int i=0;i<_BLOCK_CELLS_X_;++i){auto cell=node->block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k));if(cell)cell->Measure=volume;}
  }
  Trace::Acceleration=Force;PIC::Mover::ProcessTriangleCutFaceIntersection=Triangle;
  Trace::Settings.Adaptive=false;Trace::Settings.MaximumSubstepS=.01;
  Trace::Settings.GeometryToleranceM=1e-11;Trace::Settings.TimeToleranceS=1e-13;
}
void Start(const double* x,const double* v,double dt,bool events=false,double eventDelay=.4) {
  PIC::ParticleBuffer::DeleteAllParticles();TriangleHits=SphereHits=MovingHits=EventHits=0;HitTime=EventPosition=0;
  PIC::Mover::FindMovingSurfaceIntersection=WithMoving?MovingIntersection:nullptr;
  PIC::Mover::ProcessMovingSurfaceIntersection=WithMoving?MovingInteraction:nullptr;
  PIC::Mover::QueryTrajectoryEvent=events?Query:nullptr;PIC::Mover::AdvanceTrajectoryEvent=events?Advance:nullptr;
  PIC::Mover::ProcessTrajectoryEvent=events?Process:nullptr;
  PIC::Mover::SetUserDefinedParticleMover(Order==1?PIC::Mover::ForceCutSurfaceFirstOrder:PIC::Mover::ForceCutSurfaceSecondOrder);
  Node* node=PIC::Mesh::mesh->findTreeNode(const_cast<double*>(x));
  if(node->Thread==PIC::ThisThread) {
    int i,j,k;PIC::Mesh::mesh->FindCellIndex(const_cast<double*>(x),i,j,k,node,false);
    long int& head=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
    const long int ptr=PIC::ParticleBuffer::GetNewParticle(head);
    auto data=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
    PIC::ParticleBuffer::SetI(_NA_SPEC_,data);PIC::ParticleBuffer::SetX(const_cast<double*>(x),data);PIC::ParticleBuffer::SetV(const_cast<double*>(v),data);
    PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1/node->block->GetLocalParticleWeight(_NA_SPEC_),data);
    Trace::ResetParticleState(data);Trace::SetParticleEpochS(ptr,1000);SetEvent(ptr,{eventDelay,events?1:0});
    node->block->SetLocalTimeStep(dt,_NA_SPEC_);
  }
  double moving=0,exchange=0;PIC::TimeStepInternal::ExecutionTrackDefault(moving,exchange);
}
struct Result {double count=0,x[3]{},v[3]{},weight=0,epoch=0,clockError=0;};
Result CollectResult() {
  Result local,result;
  for(unsigned int index=0;index<PIC::DomainBlockDecomposition::nLocalBlocks;++index) {
    auto node=PIC::DomainBlockDecomposition::BlockTable[index];
    for(int cell=0;cell<_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;++cell)
      for(long int ptr=node->block->FirstCellParticleTable[cell];ptr!=-1;ptr=PIC::ParticleBuffer::GetNext(ptr)) {
        double x[3],v[3];PIC::ParticleBuffer::GetX(x,ptr);PIC::ParticleBuffer::GetV(v,ptr);
        ++local.count;for(int d=0;d<3;++d){local.x[d]+=x[d];local.v[d]+=v[d];}
        local.weight+=node->block->GetLocalParticleWeight(PIC::ParticleBuffer::GetI(ptr))*PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ptr);
        local.epoch+=Trace::GetParticleEpochS(ptr);
        const auto time=Trace::GetTimeContext(ptr);local.clockError+=std::fabs(time.RemainingTimeS)+std::fabs(time.RequestedTimeS-time.ElapsedTimeS);
      }
  }
  result.count=Sum(local.count);result.weight=Sum(local.weight);result.epoch=Sum(local.epoch);result.clockError=Sum(local.clockError);
  for(int d=0;d<3;++d){result.x[d]=Sum(local.x[d]);result.v[d]=Sum(local.v[d]);}
  return result;
}
void StateCheck(const char* name,const Result& result,const double* x,const double* v,double dt,double tolerance=1e-9,double count=1) {
  double error=std::fabs(result.count-count);
  for(int d=0;d<3;++d)error=std::max(error,std::max(std::fabs(result.x[d]-x[d]),std::fabs(result.v[d]-v[d])));
  Check(error<tolerance,name,error,0);
  Check(std::fabs(result.weight-count)<1e-12,"represented N conserved across local W",result.weight,count);
  Check(std::fabs(result.epoch-count*(1000+dt))<1e-9&&result.clockError<1e-12,"physical epoch and remaining clock",result.epoch,count*(1000+dt));
  Check(Trace::CountPending()==0,"no pending continuation after collective exchange");
  PIC::ParticleBuffer::CheckParticleList();
}
void Cases() {
  Field=0;Absorb=0;WithSibling=false;
  {const double previous=Trace::Settings.MaximumSubstepS;Trace::Settings.MaximumSubstepS=.2;
   double x[3]={.0001081,1.7,1.6},v[3]={-.0068821990473957543,0,0};Start(x,v,.2);
   double end[3]={x[0]+.2*v[0],x[1],x[2]};StateCheck("rounded x=0 face never reverses owner",CollectResult(),end,v,.2);
   Trace::Settings.MaximumSubstepS=previous;}
  {double x[3]={-1.5,1,1.5},v[3]={1,0,0};Start(x,v,3);double end[3]={1.5,1,1.5};StateCheck("arbitrary local interval and ownership",CollectResult(),end,v,3);}
  Field=1;
  {double x[3]={-1,.7,.5},v[3]={.2,0,0};Start(x,v,1);double end[3]={-.65,.6,.5},velocity[3]={.5,-.2,0};StateCheck("constant acceleration",CollectResult(),end,velocity,1);}
  Field=2;Absorb=1;
  {double x[3]={0,0,-.4},v[3]={0,0,-2};Start(x,v,2);auto result=CollectResult();const double expected=(2-std::sqrt(2.4))/2;
   Check(result.count==0&&Sum(TriangleHits)==1,"curved arc enters floor with clear chord endpoints",result.count,0);
   Check(std::fabs(Sum(HitTime)-expected)<1e-10,"first physical curved impact time",Sum(HitTime),expected);}
  Field=3;Absorb=0;
  {double x[3]={0,0,-.3},v[3]={0,0,0};Start(x,v,2);StateCheck("gravity and reflected remaining interval",CollectResult(),x,v,2);Check(Sum(TriangleHits)==1,"one floor reflection",Sum(TriangleHits),1);}
  Field=0;
  {double x[3]={.5,0,.5},v[3]={1,0,0};Start(x,v,1);double end[3]={.3,0,.5},velocity[3]={-1,0,0};StateCheck("analytic sphere reflection",CollectResult(),end,velocity,1);Check(Sum(SphereHits)==1,"one sphere reflection",Sum(SphereHits),1);}
  Field=6;
  {double x[3]={.5,0,.5},v[3]={.5,0,0};Start(x,v,1);
   const double impact=-.5+std::sqrt(1.05),incoming=.5+impact,remaining=1-impact;
   double end[3]={.9-incoming*remaining+.5*remaining*remaining,0,.5},velocity[3]={-incoming+remaining,0,0};
   StateCheck("accelerated quartic sphere reflection",CollectResult(),end,velocity,1);Check(Sum(SphereHits)==1,"one accelerated sphere reflection",Sum(SphereHits),1);}
  Field=0;WithMoving=true;
  {double x[3]={-.5,.5,.5},v[3]={0,0,-.2};Start(x,v,1);double end[3]={-.5,.5,.5},velocity[3]={0,0,.4};
   StateCheck("moving plane reflection in surface rest frame",CollectResult(),end,velocity,1);Check(Sum(MovingHits)==1,"one actual moving callback",Sum(MovingHits),1);}
  WithMoving=false;
  {double x[3]={1.5,1.8,1.5},v[3]={1,0,0};Start(x,v,1);double velocity[3]={-1,0,0};StateCheck("specular external domain",CollectResult(),x,velocity,1);}
  Field=4;
  {double x[3]={0,.5,.5},v[3]={1,0,0};Start(x,v,1);double end[3]={std::sin(1.),.5+std::cos(1.)-1,.5},velocity[3]={std::cos(1.),-std::sin(1.),0};StateCheck("velocity dependent acceleration",CollectResult(),end,velocity,1,Order==1?.02:2e-4);}
  Field=0;
  {double x[3]={-1,0,.5},v[3]={1,0,0};Start(x,v,1,true);double end[3]={.6,0,.5},velocity[3]={2,0,0};StateCheck("trajectory event species and residual interval",CollectResult(),end,velocity,1);
   Check(Sum(EventHits)==1&&std::fabs(Sum(EventPosition)+.6)<1e-10,"event at accumulated physical time",Sum(EventPosition),-.6);}
  WithSibling=true;
  {double x[3]={-1,0,.5},v[3]={1,0,0};Start(x,v,1,true);double summedEnd[3]={-.6,0,1},summedVelocity[3]={1,0,0};StateCheck("detached sibling queued exactly once",CollectResult(),summedEnd,summedVelocity,1,1e-9,2);Check(Sum(EventHits)==1,"sibling inherits no expired event",Sum(EventHits),1);}
  WithSibling=false;
  {const double previous=Trace::Settings.MaximumSubstepS;Trace::Settings.MaximumSubstepS=.125;
   double x[3]={-1,1.2,.5},v[3]={1,0,0};Start(x,v,1,true,1);
   double end[3]={0,1.2,.5},velocity[3]={2,0,0};StateCheck("internal AMR endpoint Continue event",CollectResult(),end,velocity,1);
   Check(Sum(EventHits)==1,"internal endpoint Continue processed once",Sum(EventHits),1);
   WithSibling=true;DeleteCurrent=true;Start(x,v,1,true,1);velocity[0]=-1;
   StateCheck("internal AMR endpoint Delete with zero-time child",CollectResult(),end,velocity,1);
   Check(Sum(EventHits)==1,"internal endpoint Delete processed once",Sum(EventHits),1);
   WithSibling=false;DeleteCurrent=false;Trace::Settings.MaximumSubstepS=previous;}
  if(Trace::Settings.Adaptive) {
    Field=7;const double previousMaximum=Trace::Settings.MaximumSubstepS;
    Trace::Settings.MaximumSubstepS=.11;
    double x[3]={0,0,-.79},v[3]={0,0,-.1};Start(x,v,.11);
    // z(t)=z0+v0*t+(50/6)*t^3 has a minimum ABOVE the
    // floor. Its coarse frozen-a quadratic spuriously predicts contact.
    double end[3]={0,0,-.79-.011+50*.001331/6},velocity[3]={0,0,.2025};
    StateCheck("accelerated grazing refinement with no physical impact",CollectResult(),end,velocity,.11,Order==1?1e-4:1e-5);
    Check(Sum(TriangleHits)==0,"grazing control never reflects a leaving velocity",Sum(TriangleHits),0);
    Trace::Settings.MaximumSubstepS=previousMaximum;
  }
}
}
int main(int argc,char** argv) {
  Setup();
  // The SAME cases and tolerances run with fixed/adaptive force quadrature;
  // changing this numerical execution mode does not change the oracle.
  for(int i=1;i<argc;++i) {
    if(std::strcmp(argv[i],"adaptive")==0)Trace::Settings.Adaptive=true;
    else if(std::strcmp(argv[i],"inconsistent-impact")==0) {
      // Expected-fatal fixed-step protocol: the caller must reduce its
      // numerical segment instead of receiving a nonphysical reflection.
      Order=2;Field=7;Trace::Settings.Adaptive=false;Trace::Settings.MaximumSubstepS=.11;
      double x[3]={0,0,-.79},v[3]={0,0,-.1};Start(x,v,.11);
      if(PIC::ThisThread==0)std::fprintf(stderr,"ERROR: inconsistent impact unexpectedly returned\n");MPI_Finalize();return 2;
    }
    else {if(PIC::ThisThread==0)std::fprintf(stderr,"usage: amps [adaptive]\n");MPI_Finalize();return 2;}
  }
  for(Order=1;Order<=2;++Order)Cases();
  if(PIC::ThisThread==0){fprintf(Report,"SUMMARY checks=%d passed=%d failed=%d\n",Checks,Checks-Failures,Failures);fclose(Report);printf("force-cut-surface: %d checks, %d failures\n",Checks,Failures);}
  MPI_Finalize();return Failures?1:0;
}

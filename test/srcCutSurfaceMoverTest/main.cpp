/*
 * Test of the force-free cut-surface mover PIC::Mover::ForceFreeCutSurface() (src/pic/pic_mover_force_free_cut_surface.cpp).
 *
 * Usage: amps [new|trajectory|simplified|direct] [nightly|stress] [noadversarial]
 *   new         PIC::Mover::ForceFreeCutSurface (the default)
 *   trajectory  PIC::Mover::TrajectoryTrackingMover                                        (baseline)
 *   simplified  PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory_BoundaryInjection_Simplified (baseline)
 *   direct      PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory_SecondOrder (baseline; MoverIntegratorMode=DIRECT)
 *   nightly     the nightly test (a few minutes on 8 processes)
 *   stress      more particles and steps, and (with one process) particles that cross the domain in one step
 *   noadversarial  without the adversarial class (a baseline mover that does not finish a step with it)
 *
 * The scene: six closed bodies in the box [-1,1]^3 written as one NASTRAN file and read through the usual surface
 * path (CutCell::ReadNastranSurfaceMeshLongFormat, RegisterInternalNastranSurface, InitExternalNormalVector): an
 * icosphere (level 3), a cone with a sharp tip, a cylinder, a box, an L-shaped body and a hollow cube with a closed
 * cavity (AnalyticPolyhedra.h), all rotated.  The mesh is refined near the surfaces (refinement jumps), the outer
 * boundary reflects specularly.  The inside test is the exact test of the triangulated polyhedra
 * (AnalyticPolyhedra::cScene::Depth()).
 *
 * The particles, one species per class (input: SpeciesList=H2O CO2 O N2 H):
 *   H2O uniform      uniform in the gas volume (also in the cavity), isotropic, one speed; specular reflection
 *   CO2 adversarial  aimed at vertices (many at the tip of the cone), at edges, grazing the faces at angles 0..1e-2,
 *                    starting on the surface (on faces, edges, vertices; inward, outward, random), slow to very fast;
 *                    specular reflection
 *   O   fast         isotropic, 0.88 of the smallest block per step (5 cells); stress mode with one process: 1.5 m per step
 *   N2  diffuse      Maxwellian (300 K), diffuse reflection (Surface::DiffuseReflection::Processor) at 300 K
 *   H   absorbing    removed at the first hit
 * The specular and absorbing classes do not use random numbers in the run, so their results do not depend on
 * the number of processes.
 *
 * Checks (one line "PASS|FAIL <check> ..." each, in cut-surface-mover.<mover>.<mode>.np<N>.report):
 *   geometry     every body is a closed oriented manifold with a positive volume; the faces read by AMPS are the
 *                faces of the file (vertices bit for bit, body number); every external normal of AMPS points
 *                out of its body
 *   face-lists   every face that intersects a block (exact triangle/box test) is in the cut-face list of the block
 *   initial      no particle starts inside a body
 *   leaks        after EVERY step, for every particle: no particle inside a body deeper than 1e-12 of the domain
 *                size (and none deeper than EPS)
 *   lost         no particle of the specular classes is lost; no mover error
 *   speed        |v| is conserved by specular reflection to 1e-12
 *   hits         every hit passed to the surface callback is on the plane of its face (1e-11 of the domain size),
 *                inside the closed face (barycentric 1e-9) and moves into the body
 *   reference    a part of the specular particles (all adversarial ones with a hash condition, 1/8 of the others)
 *                is integrated over every step by an independent brute-force tracer (all faces, Moeller-Trumbore in
 *                long double); unambiguous steps (no hit within 1e-9 of an edge, no grazing, no two hits within
 *                1e-9) must give the same position and velocity (1e-9)
 *   hit-rate     the hits of the uniform gas on each body agree with n<v>/4 per unit area (5 sigma)
 *   diffuse      the mean normal velocity of the diffusely reflected particles is sqrt(pi)/(2 beta) (5 sigma)
 *   absorbing    every absorbed particle is removed; their number equals the number of absorption events
 * cut-surface-mover.<mover>.<mode>.np<N>.digest: a digest of the final states (sorted) and of the hits per face of
 * the deterministic classes; it must not depend on the number of processes.
 */

#include "pic.h"
#include "Surface.h"
#include "AnalyticPolyhedra.h"

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <glob.h>
#include <signal.h>
#include <unistd.h>
#include <map>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

typedef cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> cNode;
using AnalyticPolyhedra::cBody;
using AnalyticPolyhedra::cScene;

enum {sUniform=0,sAdversarial=1,sFast=2,sDiffuse=3,sAbsorb=4,nClasses=5};
const char *ClassName[nClasses]={"uniform","adversarial","fast","diffuse","absorbing"};

//----------------------------------------- settings -----------------------------------------
std::string MoverName="new",Mode="nightly";
bool WithAdversarial=true;
int nSteps=120;
long int nUniform=40000,nFast=8000,nDiffuse=8000,nAbsorb=8000;
int AdversarialStride=3;   //every n-th vertex/edge/face gets adversarial particles

const double DomainMin=-1.0,DomainMax=1.0,DomainSize=2.0;
//The gas is enclosed by a closed shell (a body of the scene): [-0.94,0.94]^3 is the gas.  In a build with the specular
//reflection at the boundary of the domain (input/test/cut-surface-mover-specular-domain.input) there is no shell: the
//gas fills the domain and is reflected by its boundary (the mover's own treatment of the domain boundary is tested).
#if _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE_ == _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__SPECULAR_REFLECTION_
const bool SpecularDomain=true;
const double ContainerInner=1.0,ContainerOuter=1.0;
#else
const bool SpecularDomain=false;
const double ContainerInner=0.94,ContainerOuter=0.97;
#endif
const double FineCell=2.0/16/6,ShellCell=2.0/8/6,CoarseCell=2.0/4/6;  //cell sizes (levels 4, 3, 2) near the bodies, near the shell, elsewhere
double SpeedUniform=0.02,SpeedFast=0.11,SpeedAbsorb=0.03;
const double AdversarialSpeed[]={1.0E-4,2.0E-3,7.0E-3,0.02,0.05,0.11};
const int nAdversarialSpeeds=sizeof(AdversarialSpeed)/sizeof(double);
const double SurfaceTemperature=300.0;
double dtDiffuse=0.0;

const double TolStrict=1.0E-12*DomainSize;   //a particle deeper than this inside a body is a leak
const double TolPlane=1.0E-11*DomainSize;    //a reported hit point farther from the face plane is wrong

cScene Scene;
std::vector<int> AmpsFaceBody;         //body of each AMPS face
std::vector<double> FaceNormal;        //outward unit normal of each face (from the generator)

FILE *Report=NULL;
int nFailed=0;

void Check(bool pass,const char *name,const char *format,...) {
  char message[4000];
  va_list args;

  va_start(args,format);
  vsnprintf(message,sizeof(message),format,args);
  va_end(args);

  if (PIC::ThisThread==0) {
    fprintf(Report,"%s %s %s\n",pass ? "PASS" : "FAIL",name,message);
    printf("%s %s %s\n",pass ? "PASS" : "FAIL",name,message);
    fflush(Report);
    fflush(stdout);
  }

  if (!pass) nFailed++;
}

void Info(const char *format,...) {
  char message[4000];
  va_list args;

  va_start(args,format);
  vsnprintf(message,sizeof(message),format,args);
  va_end(args);

  if (PIC::ThisThread==0) {
    fprintf(Report,"INFO %s\n",message);
    printf("INFO %s\n",message);
    fflush(Report);
    fflush(stdout);
  }
}

long long SumLL(long long v) {long long g; MPI_Allreduce(&v,&g,1,MPI_LONG_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR); return g;}
double SumD(double v) {double g; MPI_Allreduce(&v,&g,1,MPI_DOUBLE,MPI_SUM,MPI_GLOBAL_COMMUNICATOR); return g;}
double MaxD(double v) {double g; MPI_Allreduce(&v,&g,1,MPI_DOUBLE,MPI_MAX,MPI_GLOBAL_COMMUNICATOR); return g;}

//------------------------------------------- scene -------------------------------------------
void BuildScene() {
  using namespace AnalyticPolyhedra;

  struct cPlacement {double center[3],euler[3];};
  const cPlacement placement[6]={
    {{-0.55,-0.45, 0.03},{0.31,0.77,1.13}},
    {{ 0.00,-0.45,-0.02},{0.13,0.41,0.27}},
    {{ 0.55,-0.45, 0.04},{0.91,1.21,0.37}},
    {{-0.55, 0.45,-0.03},{0.47,0.61,0.83}},
    {{ 0.00, 0.45, 0.02},{1.37,0.53,0.19}},
    {{ 0.55, 0.45,-0.01},{0.23,0.87,1.41}}
  };

  cBody body[6]={Icosphere(3,0.22),Cone(40,0.16,0.42),Cylinder(32,0.15,0.36,4),Box(0.34,0.26,0.30,3),LShape(0.2,3),HollowCube(0.11,2)};

  for (int b=0;b<6;b++) {
    double c[3],R[3][3],zero[3]={0.0,0.0,0.0},I[3][3]={{1,0,0},{0,1,0},{0,0,1}},shift[3];

    //center the bounding box at the origin, rotate, move to the place
    for (int i=0;i<3;i++) c[i]=-0.5*(body[b].xmin[i]+body[b].xmax[i]);
    body[b].Transform(I,c);
    RotationMatrix(placement[b].euler[0],placement[b].euler[1],placement[b].euler[2],R);
    body[b].Transform(R,zero);
    for (int i=0;i<3;i++) shift[i]=placement[b].center[i];
    body[b].Transform(I,shift);

    Scene.Add(body[b]);
  }

  //the container: a closed shell; the gas is inside it, and no particle reaches the boundary of the domain (which
  //deletes particles: a particle that passes through the shell is lost)
  if (!SpecularDomain) Scene.Add(HollowBox(ContainerInner,ContainerOuter,8));

  for (size_t b=0;b<Scene.Body.size();b++) for (int f=0;f<Scene.Body[b].nFaces();f++) {
    double n[3];
    Scene.Body[b].FaceNormal(f,n);
    FaceNormal.insert(FaceNormal.end(),n,n+3);
  }
}

//distance from the bounding spheres of the bodies (for the mesh resolution)
double DistanceToBodies(const double *x) {
  double dmin=HUGE_VAL;

  for (size_t b=0;b<Scene.Body.size();b++) {
    const cBody& body=Scene.Body[b];
    double c[3],r=0.0,d=0.0;

    if (body.Name=="shell") continue;

    for (int i=0;i<3;i++) c[i]=0.5*(body.xmin[i]+body.xmax[i]);
    for (int v=0;v<body.nVertices();v++) {
      double s=0.0;
      for (int i=0;i<3;i++) s+=pow(body.V(v)[i]-c[i],2);
      r=std::max(r,sqrt(s));
    }
    for (int i=0;i<3;i++) d+=pow(x[i]-c[i],2);
    dmin=std::min(dmin,fabs(sqrt(d)-r));
  }

  return dmin;
}

//AMPS refines a block until the diagonal of its cells is below the requested value
double Request(double cell) {return 1.01*sqrt(3.0)*cell;}

double Resolution(double *x) {
  if (DistanceToBodies(x)<0.08) return Request(FineCell);
  if (std::max(fabs(x[0]),std::max(fabs(x[1]),fabs(x[2])))>ContainerInner-0.06) return Request(ShellCell);
  return Request(CoarseCell);
}

//a face is exposed to the gas when the point just outside its center is inside the container
bool FaceExposed(const cBody& body,int f) {
  double c[3]={0,0,0},n[3];
  for (int k=0;k<3;k++) for (int i=0;i<3;i++) c[i]+=body.V(body.Face[3*f+k])[i]/3.0;
  body.FaceNormal(f,n);
  for (int i=0;i<3;i++) if (fabs(c[i]+1.0E-6*n[i])>=ContainerInner) return false;
  return true;
}
double SurfaceResolution(CutCell::cTriangleFace *f) {return Request((f->attribute==7) ? ShellCell : FineCell);}
double UnitLoad(cNode*) {return 1.0;}

void CollectLeaves(cNode* node,std::vector<cNode*>& leaves) {
  if (node->lastBranchFlag()!=_BOTTOM_BRANCH_TREE_) {
    for (int n=0;n<8;n++) if (node->downNode[n]!=NULL) CollectLeaves(node->downNode[n],leaves);
  }
  else leaves.push_back(node);
}

//triangle/box overlap (Akenine-Moeller): the box center c, half sizes h
bool TriangleBoxOverlap(const double *c,const double *h,const double *v0,const double *v1,const double *v2) {
  double a[3],b[3],d[3],e[3][3];
  for (int i=0;i<3;i++) a[i]=v0[i]-c[i],b[i]=v1[i]-c[i],d[i]=v2[i]-c[i];
  const double *v[3]={a,b,d};
  for (int i=0;i<3;i++) e[0][i]=b[i]-a[i],e[1][i]=d[i]-b[i],e[2][i]=a[i]-d[i];

  //the box faces
  for (int i=0;i<3;i++) {
    double mn=std::min(a[i],std::min(b[i],d[i])),mx=std::max(a[i],std::max(b[i],d[i]));
    if ((mn>h[i])||(mx<-h[i])) return false;
  }

  //the plane of the triangle
  double n[3];
  AnalyticPolyhedra::Cross(e[0],e[1],n);
  double r=h[0]*fabs(n[0])+h[1]*fabs(n[1])+h[2]*fabs(n[2]),s=AnalyticPolyhedra::Dot(n,a);
  if (fabs(s)>r) return false;

  //the 9 cross products of the edges with the axes
  for (int k=0;k<3;k++) for (int axis=0;axis<3;axis++) {
    double u[3]={0,0,0},ax[3];
    u[axis]=1.0;
    AnalyticPolyhedra::Cross(u,e[k],ax);
    double p0=AnalyticPolyhedra::Dot(ax,v[0]),p1=AnalyticPolyhedra::Dot(ax,v[1]),p2=AnalyticPolyhedra::Dot(ax,v[2]);
    double mn=std::min(p0,std::min(p1,p2)),mx=std::max(p0,std::max(p1,p2));
    double rr=h[0]*fabs(ax[0])+h[1]*fabs(ax[1])+h[2]*fabs(ax[2]);
    if ((mn>rr)||(mx<-rr)) return false;
  }

  return true;
}

//the present step; the state of every particle at the beginning of the step (one process: the buffer index identifies
//the particle; for the diagnostics of leaks and wrong hits)
int CurrentStep=-1;
double HitTolerance=0.0;   //the other movers: 2*EPS (set after the mesh)
struct cPreState {double x[3],v[3]; int spec;};
std::vector<cPreState> PreState;

//------------------------------------------- the surface interaction -------------------------------------------
std::vector<unsigned long long> HitCount[nClasses];  //per face
long long BadIncoming=0,OffPlane=0,OffFace=0;
double MaxPlaneDistance=0.0;
long long DiffuseEvents=0,AbsorptionEvents=0;
double DiffuseNormalSpeedSum=0.0,DiffuseNormalSpeedMin=HUGE_VAL;

int SurfaceInteraction(long int ptr,double* x,double* v,CutCell::cTriangleFace *face,cNode* node) {
  int spec=PIC::ParticleBuffer::GetI(ptr);
  int iface=face-CutCell::BoundaryTriangleFaces;
  const double *n=face->ExternalNormal;

  //the geometry of the reported hit
  double vn=v[0]*n[0]+v[1]*n[1]+v[2]*n[2];
  if (!(vn<0.0)) {
    BadIncoming++;
    if (BadIncoming<=5) {
      printf("BAD-HIT (thread %i) step %i class %i face %i v.n=%.3e x=(%.17e,%.17e,%.17e) v=(%.17e,%.17e,%.17e)\n",PIC::ThisThread,CurrentStep,spec,iface,vn,x[0],x[1],x[2],v[0],v[1],v[2]);
      if ((ptr<(long int)PreState.size())&&(PreState[ptr].spec==spec)) {
        const cPreState& p=PreState[ptr];
        printf("BAD-HIT-START (thread %i) step %i: replay with: %i %.17e %.17e %.17e %.17e %.17e %.17e\n",PIC::ThisThread,CurrentStep,spec,p.x[0],p.x[1],p.x[2],p.v[0],p.v[1],p.v[2]);
      }
    }
  }

  double dist=fabs((x[0]-face->x0Face[0])*n[0]+(x[1]-face->x0Face[1])*n[1]+(x[2]-face->x0Face[2])*n[2]);
  MaxPlaneDistance=std::max(MaxPlaneDistance,dist);
  if (dist>std::max(TolPlane,HitTolerance)) OffPlane++;

  //the other movers find a hit with the tolerance EPS (a face enlarged by EPS, a hit up to EPS behind the start)
  double xLocal[2];
  const double bary=std::max(1.0E-9,HitTolerance/std::min(face->e0Length,face->e1Length));
  face->GetProjectedLocalCoordinates(xLocal,x);
  if ((xLocal[0]<-bary)||(xLocal[1]<-bary)||(xLocal[0]+xLocal[1]>1.0+bary)) OffFace++;

  if ((spec>=0)&&(spec<nClasses)) HitCount[spec][iface]++;

  switch (spec) {
  case sDiffuse: {
    int res=Surface::DiffuseReflection::Processor(ptr,x,v,face,node);
    double vout=v[0]*n[0]+v[1]*n[1]+v[2]*n[2];
    DiffuseEvents++,DiffuseNormalSpeedSum+=vout;
    DiffuseNormalSpeedMin=std::min(DiffuseNormalSpeedMin,vout);
    return res;
  }
  case sAbsorb:
    AbsorptionEvents++;
    return _PARTICLE_DELETED_ON_THE_FACE_;
  default:
    return Surface::SpecularReflection::Processor(ptr,x,v,face,node);
  }
}

//------------------------------------------- the movers -------------------------------------------
int MoverSelection=0;
long int TracerIdOffset=-1;
unsigned long long NextTracerId=0;
void InitializeTracer(PIC::ParticleBuffer::byte* data) {
  PIC::Mover::ForceCutSurface::ResetParticleState(data);
  const unsigned long long id=(static_cast<unsigned long long>(PIC::ThisThread)<<48)|(++NextTracerId);
  std::memcpy(data+TracerIdOffset,&id,sizeof(id));
}

//debugging: AMPS_CUTMOVER_WATCHDOG=<s> stops the run when one mover call takes longer than <s> seconds and prints
//the particle (a mover that loops); the state is the one at the beginning of the call
int WatchdogSeconds=0;
double WatchX[3],WatchV[3];
int WatchSpec=-1;

void WatchdogHandler(int) {
  char buf[800];
  int n=snprintf(buf,sizeof(buf),"WATCHDOG (thread %i): a mover call took more than %i s: step %i class %i x=(%.17e,%.17e,%.17e) v=(%.17e,%.17e,%.17e)\n",
      PIC::ThisThread,WatchdogSeconds,CurrentStep,WatchSpec,WatchX[0],WatchX[1],WatchX[2],WatchV[0],WatchV[1],WatchV[2]);
  if (n>0) {ssize_t w=write(2,buf,n); (void)w;}
  _exit(3);
}

int nLostReports=0;

int Mover(long int ptr,double dt,cNode* node) {
  {
    PIC::ParticleBuffer::byte *data=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
    PIC::ParticleBuffer::GetX(WatchX,data);
    PIC::ParticleBuffer::GetV(WatchV,data);
    WatchSpec=PIC::ParticleBuffer::GetI(data);
  }

  if (WatchdogSeconds>0) alarm(WatchdogSeconds);

  int res;
  switch (MoverSelection) {
  case 0: res=PIC::Mover::ForceFreeCutSurface(ptr,dt,node); break;
  case 1: res=PIC::Mover::TrajectoryTrackingMover(ptr,dt,node); break;
  case 2: res=PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory_BoundaryInjection_Simplified(ptr,dt,node,false); break;
  case 4: res=PIC::Mover::ForceCutSurfaceFirstOrder(ptr,dt,node); break;
  case 5: res=PIC::Mover::ForceCutSurfaceSecondOrder(ptr,dt,node); break;
  default: res=PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory_SecondOrder(ptr,dt,node);
  }

  if (WatchdogSeconds>0) alarm(0);

  //a specular particle removed by the mover is lost (the 'lost' check): report the first ones with their start
  if ((res==_PARTICLE_LEFT_THE_DOMAIN_)&&((WatchSpec==sUniform)||(WatchSpec==sAdversarial)||(WatchSpec==sFast))&&(nLostReports<10)) {
    nLostReports++;
    printf("LOST (thread %i) step %i: the mover removed the particle; replay with: %i %.17e %.17e %.17e %.17e %.17e %.17e\n",PIC::ThisThread,CurrentStep,WatchSpec,
        WatchX[0],WatchX[1],WatchX[2],WatchV[0],WatchV[1],WatchV[2]);
  }

  return res;
}


void RecordPreState() {
  if (PIC::nTotalThreads!=1) return;
  PreState.assign(PIC::ParticleBuffer::GetMaxNPart(),cPreState{{0,0,0},{0,0,0},-1});

  for (cNode *node=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    PIC::Mesh::cDataBlockAMR *block=node->block;
    if (block==NULL) continue;
    for (int nd=0;nd<_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;nd++) {
      for (long int ptr=block->FirstCellParticleTable[nd];ptr!=-1;ptr=PIC::ParticleBuffer::GetNext(ptr)) {
        PIC::ParticleBuffer::byte *data=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
        cPreState& p=PreState[ptr];
        PIC::ParticleBuffer::GetX(p.x,data);
        PIC::ParticleBuffer::GetV(p.v,data);
        p.spec=PIC::ParticleBuffer::GetI(data);
      }
    }
  }
}

//replay: one particle, given by the command line, for a few steps
bool ReplayMode=false;
int ReplaySpec=0,ReplaySteps=1;
double ReplayX[3],ReplayV[3];

//the depth inside a body; the outside of the container (beyond the shell) counts as inside the shell: a particle
//there has passed through it
double LeakDepth(const double *x,int *body=NULL) {
  double depth=Scene.Depth(x,body);
  double escape=std::max(fabs(x[0]),std::max(fabs(x[1]),fabs(x[2])))-ContainerInner;
  if (escape>depth) {
    depth=escape;
    if (body!=NULL) *body=Scene.Body.size()-1;
  }
  return depth;
}

bool InsideScene(const double *x) {return LeakDepth(x)>TolStrict;}

//------------------------------------------- the reference tracer -------------------------------------------
struct cRefFace {double a[3],e1[3],e2[3],n[3],bmin[3],bmax[3];};
double Dot3(const double *a) {return a[0]*a[0]+a[1]*a[1]+a[2]*a[2];}
double AmbiguityLength=0.0;   //0: the relative scale 1e-9; the other movers: 2*EPS (set after the mesh)
std::vector<cRefFace> RefFace;

void InitReference() {
  for (size_t b=0;b<Scene.Body.size();b++) {
    const cBody& body=Scene.Body[b];
    for (int f=0;f<body.nFaces();f++) {
      cRefFace r;
      const double *v0=body.V(body.Face[3*f]),*v1=body.V(body.Face[3*f+1]),*v2=body.V(body.Face[3*f+2]);
      for (int i=0;i<3;i++) {
        r.a[i]=v0[i],r.e1[i]=v1[i]-v0[i],r.e2[i]=v2[i]-v0[i];
        r.bmin[i]=std::min(v0[i],std::min(v1[i],v2[i]))-1.0E-9;
        r.bmax[i]=std::max(v0[i],std::max(v1[i],v2[i]))+1.0E-9;
      }
      body.FaceNormal(f,r.n);
      RefFace.push_back(r);
    }
  }
}

//one step of a specular particle: brute force over all faces; 'ambiguous' when the result depends on
//rounding (a hit or a near miss within 1e-9 of an edge, grazing, two hits or a hit and a wall within 1e-9)
void ReferenceStep(double *x,double *v,double dt,bool& ambiguous) {
  double timeLeft=dt;
  ambiguous=false;

  for (int iter=0;(iter<10000)&&(timeLeft>0.0);iter++) {
    double D[3]={v[0]*timeLeft,v[1]*timeLeft,v[2]*timeLeft},lD=sqrt(D[0]*D[0]+D[1]*D[1]+D[2]*D[2]);
    if (lD==0.0) return;
    double smin[3],smax[3];
    for (int i=0;i<3;i++) smin[i]=std::min(x[i],x[i]+D[i])-AmbiguityLength,smax[i]=std::max(x[i],x[i]+D[i])+AmbiguityLength;

    double tBest=HUGE_VAL,tSecond=HUGE_VAL;
    int fBest=-1;

    for (size_t f=0;f<RefFace.size();f++) {
      const cRefFace& r=RefFace[f];
      if ((smax[0]<r.bmin[0])||(smin[0]>r.bmax[0])||(smax[1]<r.bmin[1])||(smin[1]>r.bmax[1])||(smax[2]<r.bmin[2])||(smin[2]>r.bmax[2])) continue;

      long double dn=(long double)D[0]*r.n[0]+(long double)D[1]*r.n[1]+(long double)D[2]*r.n[2];
      if (!(dn<0.0L)) continue;
      if (fabsl(dn)<1.0E-6L*lD) {
        //grazing: ambiguous when the line comes close to the face
        ambiguous=true;
      }

      //Moeller-Trumbore
      long double p[3],q[3],s[3];
      p[0]=(long double)D[1]*r.e2[2]-(long double)D[2]*r.e2[1];
      p[1]=(long double)D[2]*r.e2[0]-(long double)D[0]*r.e2[2];
      p[2]=(long double)D[0]*r.e2[1]-(long double)D[1]*r.e2[0];
      long double det=p[0]*r.e1[0]+p[1]*r.e1[1]+p[2]*r.e1[2];
      if (det==0.0L) continue;
      for (int i=0;i<3;i++) s[i]=(long double)x[i]-r.a[i];
      long double u=(s[0]*p[0]+s[1]*p[1]+s[2]*p[2])/det;
      q[0]=s[1]*r.e1[2]-s[2]*r.e1[1];
      q[1]=s[2]*r.e1[0]-s[0]*r.e1[2];
      q[2]=s[0]*r.e1[1]-s[1]*r.e1[0];
      long double w=((long double)D[0]*q[0]+(long double)D[1]*q[1]+(long double)D[2]*q[2])/det;
      long double t=(r.e2[0]*q[0]+r.e2[1]*q[1]+r.e2[2]*q[2])/det;

      //the scale of the ambiguity: 1e-9 for the new mover; the other movers decide with the tolerance EPS (a hit up to EPS
      //behind the start, a face enlarged by EPS): a step that comes within 2*EPS of an edge, or starts or ends within 2*EPS
      //of a face plane, is ambiguous for them
      const long double faceLength=sqrtl(std::max((long double)Dot3(r.e1),(long double)Dot3(r.e2)));
      const long double ambBary=std::max(1.0E-9L,(long double)AmbiguityLength/faceLength);
      const long double ambT=std::max(1.0E-9L,(long double)AmbiguityLength/lD);

      long double bmin=std::min(u,std::min(w,1.0L-u-w));
      if ((t<-ambT)||(t>1.0L+ambT)) continue;
      if (fabsl(bmin)<ambBary) ambiguous=true;   //on or near an edge
      if (bmin<0.0L) continue;
      if ((t<ambT)||(t>1.0L-ambT)) ambiguous=true;   //at the start (on the surface) or at the end
      if (t<0.0L) t=0.0L;

      if ((double)t<tBest) tSecond=tBest,tBest=(double)t,fBest=f;
      else if ((double)t<tSecond) tSecond=(double)t;
    }

    //the walls of the domain
    double tWall=HUGE_VAL,tAxis[3];
    for (int i=0;i<3;i++) {
      tAxis[i]=HUGE_VAL;
      if (D[i]>0.0) tAxis[i]=(DomainMax-x[i])/D[i];
      else if (D[i]<0.0) tAxis[i]=(DomainMin-x[i])/D[i];
      tWall=std::min(tWall,tAxis[i]);
    }

    if ((fBest>=0)&&(fabs(tSecond-tBest)*lD<std::max(1.0E-9,AmbiguityLength))) ambiguous=true;
    if ((fBest>=0)&&(fabs(tWall-tBest)*lD<std::max(1.0E-9,AmbiguityLength))) ambiguous=true;

    if ((fBest>=0)&&(tBest<=tWall)) {
      const double *n=RefFace[fBest].n;
      for (int i=0;i<3;i++) x[i]+=tBest*D[i];
      double c=v[0]*n[0]+v[1]*n[1]+v[2]*n[2];
      for (int i=0;i<3;i++) v[i]-=2.0*c*n[i];
      timeLeft*=(1.0-tBest);
    }
    else if (tWall<=1.0) {
      int nCrossed=0;
      for (int i=0;i<3;i++) {
        x[i]+=tWall*D[i];
        if (tAxis[i]<=tWall*(1.0+1.0E-12)) x[i]=(D[i]>0.0) ? DomainMax : DomainMin,v[i]=-v[i],nCrossed++;
      }
      if (nCrossed>1) ambiguous=true;
      timeLeft*=(1.0-tWall);
    }
    else {
      for (int i=0;i<3;i++) x[i]+=D[i];
      timeLeft=0.0;
    }
  }
}

//------------------------------------------- particles -------------------------------------------
std::mt19937_64 Rng(20261116);
double U() {return std::uniform_real_distribution<double>(0.0,1.0)(Rng);}

void RandomDirection(double *d) {
  double z=2.0*U()-1.0,phi=2.0*M_PI*U(),s=sqrt(std::max(0.0,1.0-z*z));
  d[0]=s*cos(phi),d[1]=s*sin(phi),d[2]=z;
}

//a random direction in the hemisphere around n
void RandomDirection(double *d,const double *n) {
  RandomDirection(d);
  if (AnalyticPolyhedra::Dot(d,n)<0.0) for (int i=0;i<3;i++) d[i]=-d[i];
}

long int nCreated[nClasses]={0,0,0,0,0},nGeneratedInside=0;

//create the particle on the process that owns its block (every process generates the same sequence)
void AddParticle(int spec,const double *xIn,const double *vIn) {
  double x[3]={xIn[0],xIn[1],xIn[2]},v[3]={vIn[0],vIn[1],vIn[2]};

  if (LeakDepth(x)>TolStrict) {nGeneratedInside++; return;}   //inside a body or outside the container

  cNode *node=PIC::Mesh::mesh->findTreeNode(x);
  if (node==NULL) return;

  nCreated[spec]++;
  if (node->Thread!=PIC::ThisThread) return;

  PIC::ParticleBuffer::InitiateParticle(x,v,NULL,&spec,NULL,_PIC_INIT_PARTICLE_MODE__ADD2LIST_,(void*)node,
      MoverSelection>=4?InitializeTracer:nullptr);
}

double SpeedOfClass(int spec) {
  switch (spec) {
  case sUniform: return SpeedUniform;
  case sFast: return SpeedFast;
  case sAbsorb: return SpeedAbsorb;
  default: return 0.0;
  }
}

void RandomGasPoint(double *x) {
  do {
    for (int i=0;i<3;i++) x[i]=ContainerInner*(2.0*U()-1.0);
  }
  while (Scene.Depth(x)>-1.0E-9);
}

void CreateParticles() {
  double x[3],v[3],d[3];

  //uniform, fast, absorbing: isotropic, one speed each
  const int iso[3]={sUniform,sFast,sAbsorb};
  const long int nIso[3]={nUniform,nFast,nAbsorb};

  for (int c=0;c<3;c++) for (long int n=0;n<nIso[c];n++) {
    RandomGasPoint(x);
    RandomDirection(d);
    for (int i=0;i<3;i++) v[i]=SpeedOfClass(iso[c])*d[i];
    AddParticle(iso[c],x,v);
  }

  //diffuse: Maxwellian
  {
    const double m=PIC::MolecularData::GetMass(sDiffuse),vth=sqrt(Kbol*SurfaceTemperature/m);
    std::normal_distribution<double> g(0.0,vth);
    for (long int n=0;n<nDiffuse;n++) {
      RandomGasPoint(x);
      for (int i=0;i<3;i++) v[i]=g(Rng);
      AddParticle(sDiffuse,x,v);
    }
  }

  //adversarial
  if (!WithAdversarial) return;
  int iSpeed=0;
  auto Speed=[&]() {iSpeed=(iSpeed+1)%nAdversarialSpeeds; return AdversarialSpeed[iSpeed];};

  if (SpecularDomain) {
    //aimed exactly at the corners and at points of the edges of the domain, and starting on its faces
    for (int c=0;c<8;c++) for (int k=0;k<20;k++) {
      double corner[3]={(c&1) ? DomainMax : DomainMin,(c&2) ? DomainMax : DomainMin,(c&4) ? DomainMax : DomainMin},n[3];
      for (int i=0;i<3;i++) n[i]=-corner[i]/sqrt(3.0);
      double sp=Speed(),s=FineCell*(0.2+2.3*U());
      RandomDirection(d,n);
      if (k==0) for (int i=0;i<3;i++) d[i]=n[i];
      for (int i=0;i<3;i++) x[i]=corner[i]+s*d[i],v[i]=-sp*d[i];
      AddParticle(sAdversarial,x,v);
    }

    for (int axis=0;axis<3;axis++) for (int c=0;c<4;c++) for (int k=0;k<10;k++) {
      double e[3],n[3]={0,0,0};
      int a1=(axis+1)%3,a2=(axis+2)%3;
      e[axis]=DomainMin+(DomainMax-DomainMin)*U();
      e[a1]=(c&1) ? DomainMax : DomainMin,e[a2]=(c&2) ? DomainMax : DomainMin;
      n[a1]=-e[a1],n[a2]=-e[a2];
      double l=sqrt(Dot3(n)),sp=Speed(),s=FineCell*(0.2+2.3*U());
      for (int i=0;i<3;i++) n[i]/=l;
      RandomDirection(d,n);
      for (int i=0;i<3;i++) x[i]=e[i]+s*d[i],v[i]=-sp*d[i];
      AddParticle(sAdversarial,x,v);
    }

    for (int face=0;face<6;face++) for (int k=0;k<20;k++) {
      int axis=face/2;
      double p[3];
      for (int i=0;i<3;i++) p[i]=DomainMin+(DomainMax-DomainMin)*U();
      p[axis]=(face%2==0) ? DomainMin : std::nextafter(DomainMax,0.0);
      double sp=Speed();
      RandomDirection(d);
      for (int i=0;i<3;i++) v[i]=sp*d[i];
      AddParticle(sAdversarial,p,v);
    }
  }

  for (size_t b=0;b<Scene.Body.size();b++) {
    const cBody& body=Scene.Body[b];

    //the shell: only its inner surface faces the gas
    auto GasSide=[&](int iv) {
      if (body.Name!="shell") return true;
      const double *V=body.V(iv);
      return std::max(fabs(V[0]),std::max(fabs(V[1]),fabs(V[2])))<=ContainerInner*(1.0+1.0E-12);
    };

    //the vertex normals and the edges
    std::vector<double> vn(3*body.nVertices(),0.0);
    std::map<std::pair<int,int>,std::vector<int> > edges;
    for (int f=0;f<body.nFaces();f++) {
      double n[3];
      body.FaceNormal(f,n);
      for (int k=0;k<3;k++) {
        for (int i=0;i<3;i++) vn[3*body.Face[3*f+k]+i]+=n[i];
        int a=body.Face[3*f+k],c=body.Face[3*f+(k+1)%3];
        edges[std::make_pair(std::min(a,c),std::max(a,c))].push_back(f);
      }
    }

    //aimed at the vertices; the vertex of the cone with the largest z in local coordinates is its tip: all
    //vertices of a cone with the largest number of faces get many particles
    std::vector<int> nVertexFaces(body.nVertices(),0);
    for (int f=0;f<body.nFaces();f++) for (int k=0;k<3;k++) nVertexFaces[body.Face[3*f+k]]++;

    for (int iv=0;iv<body.nVertices();iv++) {
      bool tip=(body.Name=="cone")&&(nVertexFaces[iv]>=20)&&(iv==0);
      if (!tip&&(iv%AdversarialStride!=0)) continue;
      if (!GasSide(iv)) continue;
      const double *V=body.V(iv);
      double n[3]={vn[3*iv],vn[3*iv+1],vn[3*iv+2]},l=sqrt(AnalyticPolyhedra::Dot(n,n));
      for (int i=0;i<3;i++) n[i]/=l;

      int nAim=tip ? 200 : 2;
      for (int k=0;k<nAim;k++) {
        for (int itry=0;itry<20;itry++) {
          double s=FineCell*(0.2+2.3*U()),sp=Speed();
          RandomDirection(d,n);
          if (tip&&(k==0)) for (int i=0;i<3;i++) d[i]=n[i];   //straight at the tip along its axis
          for (int i=0;i<3;i++) x[i]=V[i]+s*d[i],v[i]=-sp*d[i];
          if (LeakDepth(x)<=0.0) {AddParticle(sAdversarial,x,v); break;}
        }
      }

      //starting on the vertex: inward, outward, random
      double dirs[3][3]={{-n[0],-n[1],-n[2]},{n[0],n[1],n[2]},{0,0,0}};
      RandomDirection(dirs[2]);
      for (int k=0;k<3;k++) {
        double sp=Speed();
        for (int i=0;i<3;i++) v[i]=sp*dirs[k][i];
        AddParticle(sAdversarial,V,v);
      }
    }

    //aimed at the edges, and starting on them
    int iedge=0;
    for (std::map<std::pair<int,int>,std::vector<int> >::const_iterator it=edges.begin();it!=edges.end();it++,iedge++) {
      if (iedge%AdversarialStride!=0) continue;
      if (!GasSide(it->first.first)||!GasSide(it->first.second)) continue;
      const double *A=body.V(it->first.first),*B=body.V(it->first.second);
      double n[3]={0,0,0};
      for (size_t k=0;k<it->second.size();k++) {double nf[3]; body.FaceNormal(it->second[k],nf); for (int i=0;i<3;i++) n[i]+=nf[i];}
      double l=sqrt(AnalyticPolyhedra::Dot(n,n));
      if (l==0.0) continue;
      for (int i=0;i<3;i++) n[i]/=l;

      const double tau[3]={0.5,U(),1.0E-9};
      for (int k=0;k<3;k++) {
        double E[3];
        for (int i=0;i<3;i++) E[i]=A[i]+tau[k]*(B[i]-A[i]);
        for (int itry=0;itry<20;itry++) {
          double s=FineCell*(0.2+2.3*U()),sp=Speed();
          RandomDirection(d,n);
          for (int i=0;i<3;i++) x[i]=E[i]+s*d[i],v[i]=-sp*d[i];
          if (LeakDepth(x)<=0.0) {AddParticle(sAdversarial,x,v); break;}
        }
      }

      double E[3];
      for (int i=0;i<3;i++) E[i]=A[i]+tau[1]*(B[i]-A[i]);
      double sp=Speed();
      for (int i=0;i<3;i++) v[i]=-sp*n[i];
      AddParticle(sAdversarial,E,v);
      RandomDirection(d);
      for (int i=0;i<3;i++) v[i]=sp*d[i];
      AddParticle(sAdversarial,E,v);
    }

    //grazing, and starting on the faces
    for (int f=0;f<body.nFaces();f++) {
      if (f%AdversarialStride!=0) continue;
      if (!FaceExposed(body,f)) continue;
      const double *a=body.V(body.Face[3*f]),*bb=body.V(body.Face[3*f+1]),*c=body.V(body.Face[3*f+2]);
      double n[3],t1[3],t2[3],Q[3];
      body.FaceNormal(f,n);
      double r1=U(),r2=U();
      if (r1+r2>1.0) r1=1.0-r1,r2=1.0-r2;
      for (int i=0;i<3;i++) Q[i]=a[i]+r1*(bb[i]-a[i])+r2*(c[i]-a[i]);

      //a tangent direction
      RandomDirection(d);
      AnalyticPolyhedra::Cross(n,d,t1);
      double l=sqrt(AnalyticPolyhedra::Dot(t1,t1));
      for (int i=0;i<3;i++) t1[i]/=l;
      AnalyticPolyhedra::Cross(n,t1,t2);

      const double alpha[6]={0.0,1.0E-12,1.0E-8,1.0E-4,1.0E-2,-1.0E-8};   //>0: into the body
      const double height[4]={0.0,1.0E-14,1.0E-9,1.0E-6};
      for (int ka=0;ka<6;ka++) {
        int kh=(f/AdversarialStride+ka)%4;
        double sp=Speed();
        //start upstream so that the particle reaches Q within about one step
        for (int i=0;i<3;i++) {
          double dir=cos(alpha[ka])*t1[i]-sin(alpha[ka])*n[i];
          v[i]=sp*dir;
          x[i]=Q[i]+height[kh]*DomainSize*n[i]-0.5*sp*dir;
        }
        //keep the start above the plane (a large angle moved it up by 0.5*sp*sin(alpha))
        AddParticle(sAdversarial,x,v);
      }

      //starting on the face: inward, outward, random, tangent
      double dirs[4][3]={{-n[0],-n[1],-n[2]},{n[0],n[1],n[2]},{0,0,0},{t2[0],t2[1],t2[2]}};
      RandomDirection(dirs[2]);
      for (int k=0;k<4;k++) {
        double sp=Speed();
        for (int i=0;i<3;i++) v[i]=sp*dirs[k][i];
        AddParticle(sAdversarial,Q,v);
      }
    }
  }
}

//------------------------------------------- per-step checks -------------------------------------------
long long LeakStrict[nClasses]={0},LeakEPS[nClasses]={0};
int LeakReportLimit=10;
double MaxLeakDepth=0.0,MaxSpeedError=0.0;
int nLeakReports=0;

double SpeedError(int spec,const double *v) {
  double s=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]),e=HUGE_VAL;
  if (spec==sAdversarial) {
    for (int k=0;k<nAdversarialSpeeds;k++) e=std::min(e,fabs(s/AdversarialSpeed[k]-1.0));
  }
  else e=fabs(s/SpeedOfClass(spec)-1.0);
  return e;
}

//leaks (removed after counting) and the speeds; returns the number of particles of each class
void CheckParticles(int step,long long *count) {
  for (int s=0;s<nClasses;s++) count[s]=0;

  for (cNode *node=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    PIC::Mesh::cDataBlockAMR *block=node->block;
    if (block==NULL) continue;

    for (int nd=0;nd<_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;nd++) {
      long int ptr=block->FirstCellParticleTable[nd];

      while (ptr!=-1) {
        long int next=PIC::ParticleBuffer::GetNext(ptr);
        PIC::ParticleBuffer::byte *data=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
        double x[3],v[3];
        int spec=PIC::ParticleBuffer::GetI(data);

        PIC::ParticleBuffer::GetX(x,data);
        PIC::ParticleBuffer::GetV(v,data);

        int body;
        double depth=LeakDepth(x,&body);

        if (depth>TolStrict) {
          LeakStrict[spec]++;
          if (depth>PIC::Mesh::mesh->EPS) LeakEPS[spec]++;
          MaxLeakDepth=std::max(MaxLeakDepth,depth);

          if (nLeakReports<LeakReportLimit) {
            nLeakReports++;
            printf("LEAK (thread %i) step %i class %s body %s depth %.3e x=(%.17e,%.17e,%.17e) v=(%.17e,%.17e,%.17e)\n",PIC::ThisThread,step,ClassName[spec],
                Scene.Body[body].Name.c_str(),depth,x[0],x[1],x[2],v[0],v[1],v[2]);
            if ((ptr<(long int)PreState.size())&&(PreState[ptr].spec==spec)) {
              const cPreState& p=PreState[ptr];
              printf("LEAK-START (thread %i) step %i: replay with: %i %.17e %.17e %.17e %.17e %.17e %.17e\n",PIC::ThisThread,step,spec,p.x[0],p.x[1],p.x[2],p.v[0],p.v[1],p.v[2]);
            }
          }

          PIC::ParticleBuffer::DeleteParticle(ptr,block->FirstCellParticleTable[nd]);
        }
        else {
          count[spec]++;
          if ((spec==sUniform)||(spec==sAdversarial)||(spec==sFast)) MaxSpeedError=std::max(MaxSpeedError,SpeedError(spec,v));
        }

        ptr=next;
      }
    }
  }
}

//------------------------------------------- the reference comparison -------------------------------------------
struct cRefParticle {long int ptr; double x[3],v[3]; int spec;unsigned long long tracerId;};
std::vector<cRefParticle> RefParticles;
long long RefCompared=0,RefAmbiguous=0,RefMismatch=0;
//the position tolerance of the comparison: 1e-9 m for the new mover; the other movers move a hit point into its face
//(1e-5 of the face in local coordinates) and particles off the surface or into a block (by up to 10*EPS): 1e-4 m.
//The velocity must agree to 1e-9 for every mover: the same faces in the same order.
double RefPositionTolerance=1.0E-9,MaxRefPositionDifference=0.0;
int nRefReports=0;

unsigned long long HashState(const double *x) {
  unsigned long long h=1469598103934665603ULL;
  const unsigned char *p=(const unsigned char*)x;
  for (int i=0;i<24;i++) h=(h^p[i])*1099511628211ULL;
  return h;
}

void CollectReferenceParticles() {
  RefParticles.clear();

  for (cNode *node=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    PIC::Mesh::cDataBlockAMR *block=node->block;
    if (block==NULL) continue;

    for (int nd=0;nd<_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;nd++) {
      for (long int ptr=block->FirstCellParticleTable[nd];ptr!=-1;ptr=PIC::ParticleBuffer::GetNext(ptr)) {
        cRefParticle r;
        PIC::ParticleBuffer::byte *data=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
        r.spec=PIC::ParticleBuffer::GetI(data);
        if ((r.spec!=sUniform)&&(r.spec!=sAdversarial)&&(r.spec!=sFast)) continue;

        PIC::ParticleBuffer::GetX(r.x,data);
        PIC::ParticleBuffer::GetV(r.v,data);
        unsigned long long h=HashState(r.x);
        if ((r.spec==sAdversarial) ? (h%4!=0) : (h%8!=0)) continue;
        r.ptr=ptr;
        if(MoverSelection>=4)std::memcpy(&r.tracerId,data+TracerIdOffset,sizeof(r.tracerId));
        RefParticles.push_back(r);
      }
    }
  }
  if(MoverSelection>=4) {
    // Rank-local particle indices change under real MPI transfer. Preserve
    // an opaque fixture identity through production cloning/packing and
    // gather independent starting states, never compare reused old indices.
    const int bytes=static_cast<int>(RefParticles.size()*sizeof(cRefParticle));
    std::vector<int> counts(PIC::nTotalThreads),offsets(PIC::nTotalThreads);
    MPI_Allgather(&bytes,1,MPI_INT,counts.data(),1,MPI_INT,MPI_GLOBAL_COMMUNICATOR);
    int total=0;for(int i=0;i<PIC::nTotalThreads;++i){offsets[i]=total;total+=counts[i];}
    std::vector<cRefParticle> all(total/sizeof(cRefParticle));
    MPI_Allgatherv(RefParticles.data(),bytes,MPI_BYTE,all.data(),counts.data(),offsets.data(),MPI_BYTE,MPI_GLOBAL_COMMUNICATOR);
    RefParticles.swap(all);
  }
}

//after the move, before the particle exchange: the moved particles still have their buffer index
void CompareWithReference() {
  std::unordered_map<unsigned long long,long int> ownedTracers;
  long long found=0;
  if(MoverSelection>=4) {
    for(cNode* node=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];node;node=node->nextNodeThisThread)
      if(node->block)for(int cell=0;cell<_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;++cell)
        for(long int ptr=node->block->FirstCellParticleTable[cell];ptr!=-1;ptr=PIC::ParticleBuffer::GetNext(ptr)) {
          unsigned long long id=0;std::memcpy(&id,PIC::ParticleBuffer::GetParticleDataPointer(ptr)+TracerIdOffset,sizeof(id));ownedTracers[id]=ptr;
        }
  }
  for (size_t n=0;n<RefParticles.size();n++) {
    cRefParticle& r=RefParticles[n];
    long int ptr=r.ptr;
    if(MoverSelection>=4){auto current=ownedTracers.find(r.tracerId);if(current==ownedTracers.end())continue;ptr=current->second;++found;}
    if (!PIC::ParticleBuffer::IsParticleAllocated(ptr)) continue;

    double x[3],v[3];
    bool ambiguous;
    PIC::ParticleBuffer::byte *data=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
    PIC::ParticleBuffer::GetX(x,data);
    PIC::ParticleBuffer::GetV(v,data);

    double xr[3]={r.x[0],r.x[1],r.x[2]},vr[3]={r.v[0],r.v[1],r.v[2]};
    ReferenceStep(xr,vr,PIC::ParticleWeightTimeStep::GlobalTimeStep[r.spec],ambiguous);

    if (ambiguous) {RefAmbiguous++; continue;}
    RefCompared++;

    double s=sqrt(vr[0]*vr[0]+vr[1]*vr[1]+vr[2]*vr[2]),ex=0.0,ev=0.0;
    for (int i=0;i<3;i++) ex=std::max(ex,fabs(x[i]-xr[i])),ev=std::max(ev,fabs(v[i]-vr[i]));

    MaxRefPositionDifference=std::max(MaxRefPositionDifference,ex);
    if ((ex>RefPositionTolerance)||(ev>1.0E-9*s)) {
      RefMismatch++;
      if (nRefReports<10) {
        nRefReports++;
        printf("REFERENCE (thread %i) class %s start x=(%.17e,%.17e,%.17e) v=(%.17e,%.17e,%.17e): mover x=(%.17e,%.17e,%.17e) v=(%.17e,%.17e,%.17e), reference x=(%.17e,%.17e,%.17e) v=(%.17e,%.17e,%.17e)\n",
            PIC::ThisThread,ClassName[r.spec],r.x[0],r.x[1],r.x[2],r.v[0],r.v[1],r.v[2],x[0],x[1],x[2],v[0],v[1],v[2],xr[0],xr[1],xr[2],vr[0],vr[1],vr[2]);
      }
    }
  }
  if(MoverSelection>=4) {
    const long long globalFound=SumLL(found);
    // Every selected reference particle is specular, so disappearance is an
    // actual mismatch rather than a silently skipped MPI reference check.
    if(PIC::ThisThread==0&&globalFound!=static_cast<long long>(RefParticles.size()))
      RefMismatch+=std::llabs(globalFound-static_cast<long long>(RefParticles.size()));
  }
}

//------------------------------------------- the digest -------------------------------------------
unsigned long long Fnv(const void *p,size_t n,unsigned long long h=1469598103934665603ULL) {
  const unsigned char *c=(const unsigned char*)p;
  for (size_t i=0;i<n;i++) h=(h^c[i])*1099511628211ULL;
  return h;
}

void WriteDigest(const char *fname) {
  //the final states of the deterministic classes, sorted
  std::vector<double> local;

  for (cNode *node=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    PIC::Mesh::cDataBlockAMR *block=node->block;
    if (block==NULL) continue;

    for (int nd=0;nd<_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;nd++) {
      for (long int ptr=block->FirstCellParticleTable[nd];ptr!=-1;ptr=PIC::ParticleBuffer::GetNext(ptr)) {
        PIC::ParticleBuffer::byte *data=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
        int spec=PIC::ParticleBuffer::GetI(data);
        if (spec==sDiffuse) continue;
        double x[3],v[3];
        PIC::ParticleBuffer::GetX(x,data);
        PIC::ParticleBuffer::GetV(v,data);
        local.push_back(spec);
        local.insert(local.end(),x,x+3);
        local.insert(local.end(),v,v+3);
      }
    }
  }

  int n=local.size(),nAll=0;
  std::vector<int> counts(PIC::nTotalThreads),displ(PIC::nTotalThreads);
  MPI_Gather(&n,1,MPI_INT,counts.data(),1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);
  if (PIC::ThisThread==0) for (int t=0;t<PIC::nTotalThreads;t++) displ[t]=nAll,nAll+=counts[t];
  std::vector<double> all((PIC::ThisThread==0) ? std::max(nAll,1) : 1);
  MPI_Gatherv(local.data(),n,MPI_DOUBLE,all.data(),counts.data(),displ.data(),MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);

  //the hits per face of the deterministic classes
  std::vector<unsigned long long> hits[nClasses];
  for (int s=0;s<nClasses;s++) {
    hits[s].resize(HitCount[s].size());
    MPI_Reduce(HitCount[s].data(),hits[s].data(),HitCount[s].size(),MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
  }

  if (PIC::ThisThread==0) {
    std::vector<std::vector<double> > p(nAll/7);
    for (int i=0;i<nAll/7;i++) p[i].assign(all.begin()+7*i,all.begin()+7*i+7);
    std::sort(p.begin(),p.end());

    FILE *f=fopen(fname,"w");
    for (int s=0;s<nClasses;s++) {
      if (s==sDiffuse) continue;
      unsigned long long h=1469598103934665603ULL;
      long int np=0;
      for (size_t i=0;i<p.size();i++) if ((int)p[i][0]==s) h=Fnv(p[i].data(),7*sizeof(double),h),np++;
      unsigned long long hh=Fnv(hits[s].data(),hits[s].size()*sizeof(unsigned long long));
      unsigned long long total=0;
      for (size_t i=0;i<hits[s].size();i++) total+=hits[s][i];
      fprintf(f,"%s: %li particles, state digest %016llx; %llu hits, hits-per-face digest %016llx\n",ClassName[s],np,h,total,hh);
    }
    fclose(f);
  }
}

}  // namespace


int main(int argc,char **argv) {
  PIC::Init_BeforeParser();
  std::string ReuseMesh;

  for (int i=1;i<argc;i++) {
    std::string a=argv[i];
    if ((a=="new")||(a=="trajectory")||(a=="simplified")||(a=="direct")||(a=="force-first")||(a=="force-second")) MoverName=a;
    else if ((a=="nightly")||(a=="stress")||(a=="movers")) Mode=a;
    else if (a=="noadversarial") WithAdversarial=false;
    else if(a=="reuse-mesh"&&i+1<argc)ReuseMesh=argv[++i];
    else if ((a=="replay")&&(i+7<argc)) {
      //replay <class> x y z vx vy vz [steps]
      ReplayMode=true,ReplaySpec=atoi(argv[i+1]);
      for (int k=0;k<3;k++) ReplayX[k]=atof(argv[i+2+k]),ReplayV[k]=atof(argv[i+5+k]);
      i+=7;
      if ((i+1<argc)&&(isdigit(argv[i+1][0]))) ReplaySteps=atoi(argv[++i]);
    }
    else {
      if (PIC::ThisThread==0) printf("usage: %s [new|trajectory|simplified|direct|force-first|force-second] [nightly|movers|stress] [noadversarial] [replay class x y z vx vy vz [steps]]\n",argv[0]);
      MPI_Finalize();
      return 2;
    }
  }

  MoverSelection=(MoverName=="new") ? 0 : ((MoverName=="trajectory") ? 1 : ((MoverName=="simplified") ? 2 : 3));
  if(MoverName=="force-first")MoverSelection=4;
  if(MoverName=="force-second")MoverSelection=5;
  if (MoverSelection!=0&&MoverSelection<4) RefPositionTolerance=1.0E-4;
  if(MoverSelection>=4) {
    // Registration is runtime opt-in; default/old mover fixtures allocate
    // exactly their historical layout and keep their existing execution.
    PIC::Mover::ForceCutSurface::RegisterStorage();
    PIC::ParticleBuffer::RequestDataStorage(TracerIdOffset,sizeof(unsigned long long),alignof(unsigned long long));
    PIC::Mover::ForceCutSurface::Acceleration=[](double* a,int,long int,const double*,const double*,cNode*,const PIC::Mover::ForceCutSurface::cTimeContext&){a[0]=a[1]=a[2]=0;};
    PIC::Mover::ForceCutSurface::Settings.Adaptive=false;
    PIC::Mover::ForceCutSurface::Settings.GeometryToleranceM=2e-11;
  }

  if (Mode=="movers") {
    //the existing movers in the nightly test: fewer particles and steps
    nSteps=80,nUniform=20000,nFast=4000,nDiffuse=4000,nAbsorb=4000,AdversarialStride=6;
  }

  if (Mode.compare(0,6,"stress")==0) {
    nSteps=400,nUniform=150000,nFast=30000,nDiffuse=30000,nAbsorb=20000,AdversarialStride=1;
    if (PIC::nTotalThreads==1) SpeedFast=1.5;   //crosses the domain in one step: only one process has all blocks
  }

  char fname[400];
  if (!WithAdversarial) Mode+="-noadversarial";
  if (ReplayMode) Mode+="-replay",nSteps=ReplaySteps,LeakReportLimit=1000000;
  if (getenv("AMPS_CUTMOVER_WATCHDOG")!=NULL) {
    WatchdogSeconds=atoi(getenv("AMPS_CUTMOVER_WATCHDOG"));
    signal(SIGALRM,WatchdogHandler);
  }
  sprintf(fname,"cut-surface-mover.%s.%s.np%i.report",MoverName.c_str(),Mode.c_str(),PIC::nTotalThreads);
  if (PIC::ThisThread==0) Report=fopen(fname,"w");

  //1. the scene, written as a NASTRAN file and read as every surface
  BuildScene();
  InitReference();

  {
    bool ok=true;
    std::string msg,all;
    for (size_t b=0;b<Scene.Body.size();b++) {
      bool closed=Scene.Body[b].CheckClosedManifold(msg);
      double vol=Scene.Body[b].Volume();
      ok=ok&&closed&&(vol>0.0);
      char buf[200];
      sprintf(buf,"; volume %.6e area %.6e",vol,Scene.Body[b].Area());
      all+=msg+buf+((b+1<Scene.Body.size()) ? " | " : "");
    }
    Check(ok,"geometry-bodies","%s",all.c_str());
  }

  if (PIC::ThisThread==0) {
    Scene.WriteNastran("cut-surface-mover.bdf","AMPS cut-surface mover test scene (AnalyticPolyhedra.h)");

    //the external normals of an earlier run in this directory are not used
    glob_t g;
    if (glob("amr.sig=*.TriangulationExternalNormalVector.bin",0,NULL,&g)==0) {
      for (size_t i=0;i<g.gl_pathc;i++) remove(g.gl_pathv[i]);
      globfree(&g);
    }
  }
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

  CutCell::ReadNastranSurfaceMeshLongFormat("cut-surface-mover.bdf");
  PIC::BC::InternalBoundary::NastranSurface::RegisterInternalNastranSurface();
  PIC::Mesh::mesh->CutCellSurfaceLocalResolution=SurfaceResolution;

  //2. the mesh: generated by process 0 and read from the file by the others
  double xmin[3]={DomainMin,DomainMin,DomainMin},xmax[3]={DomainMax,DomainMax,DomainMax};

  PIC::Mesh::mesh->AllowBlockAllocation=false;
  PIC::Mesh::mesh->init(xmin,xmax,Resolution);

  if(!ReuseMesh.empty()) {
    // An explicit verification-runner optimization, never an implicit cache:
    // the runner preserves/checks the NASTRAN/config/mesh hashes from the
    // producing run. All geometry, per-face coverage and acceptance checks
    // below still execute. No particle state or reference result is loaded.
    PIC::Mesh::mesh->readMeshFile(ReuseMesh.c_str());
  }
  else if (PIC::ThisThread==0) {
    PIC::Mesh::mesh->buildMesh();
    PIC::Mesh::mesh->saveMeshFile("cut-surface-mover.mesh");
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  }
  else {
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    PIC::Mesh::mesh->readMeshFile("cut-surface-mover.mesh");
  }

  PIC::Mesh::mesh->SetParallelLoadMeasure(UnitLoad);
  PIC::Mesh::mesh->CreateNewParallelDistributionLists();
  PIC::Mesh::initCellSamplingDataBuffer();
  PIC::Mesh::mesh->AllowBlockAllocation=true;
  PIC::Mesh::mesh->AllocateTreeBlocks();

  PIC::Mesh::IrregularSurface::InitExternalNormalVector();

  PIC::Init_AfterParser();
  PIC::Mover::Init();

  //The cell measures: the existing movers check them in the debugger mode, which the nightly test mode turns on
  //(-amps-test=on, which an earlier test of the nightly run leaves set: the Table entry does not set it); the motion does
  //not use them.  Every cell gets its full volume, so that the test runs the same in both modes.
  {
    std::vector<cNode*> leaves;
    CollectLeaves(PIC::Mesh::mesh->rootTree,leaves);

    for (cNode *node : leaves) if (node->block!=NULL) {
      const double vCell=(node->xmax[0]-node->xmin[0])*(node->xmax[1]-node->xmin[1])*(node->xmax[2]-node->xmin[2])/
          (_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_);

      for (int k=0;k<_BLOCK_CELLS_Z_;k++) for (int j=0;j<_BLOCK_CELLS_Y_;j++) for (int i=0;i<_BLOCK_CELLS_X_;i++) {
        PIC::Mesh::cDataCenterNode *cell=node->block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k));
        if ((cell!=NULL)&&(cell->Measure<=0.0)) cell->Measure=vCell;
      }
    }
  }

  //3. the faces read by AMPS and their normals
  {
    long long badVertex=0,badBody=0,badNormal=0;
    int nAmps=CutCell::nBoundaryTriangleFaces,f=0;

    AmpsFaceBody.assign(nAmps,-1);
    for (size_t b=0;b<Scene.Body.size();b++) for (int k=0;k<Scene.Body[b].nFaces();k++,f++) {
      if (f>=nAmps) break;
      CutCell::cTriangleFace& face=CutCell::BoundaryTriangleFaces[f];
      const cBody& body=Scene.Body[b];
      const double *v[3]={face.x0Face,face.x1Face,face.x2Face};
      for (int j=0;j<3;j++) if (memcmp(v[j],body.V(body.Face[3*k+j]),3*sizeof(double))!=0) badVertex++;
      if (face.attribute!=(int)b+1) badBody++;
      if (AnalyticPolyhedra::Dot(face.ExternalNormal,&FaceNormal[3*f])<1.0-1.0E-12) badNormal++;
      AmpsFaceBody[f]=b;
    }

    Check((nAmps==Scene.nFaces())&&(badVertex==0)&&(badBody==0)&&(badNormal==0),"geometry-amps",
        "%i faces read by AMPS (%i written); %lli vertices differ from the file, %lli faces have another body, %lli external normals do not point out of the body",
        nAmps,Scene.nFaces(),badVertex,badBody,badNormal);
  }

  //4. the cut-face lists of the blocks (every process checks a part of the leaves)
  {
    std::vector<cNode*> leaves;
    CollectLeaves(PIC::Mesh::mesh->rootTree,leaves);
    long long missing=0,listed=0,checked=0;

    for (size_t n=PIC::ThisThread;n<leaves.size();n+=PIC::nTotalThreads) {
      cNode *node=leaves[n];
      std::vector<char> inList(CutCell::nBoundaryTriangleFaces,0);
      for (CutCell::cTriangleFaceDescriptor *t=node->FirstTriangleCutFace;t!=NULL;t=t->next) inList[t->TriangleFace-CutCell::BoundaryTriangleFaces]=1,listed++;

      double c[3],h[3];
      for (int i=0;i<3;i++) c[i]=0.5*(node->xmin[i]+node->xmax[i]),h[i]=0.5*(node->xmax[i]-node->xmin[i])*(1.0-1.0E-9);

      for (int f=0;f<CutCell::nBoundaryTriangleFaces;f++) {
        CutCell::cTriangleFace& face=CutCell::BoundaryTriangleFaces[f];
        checked++;
        if (TriangleBoxOverlap(c,h,face.x0Face,face.x1Face,face.x2Face)&&(inList[f]==0)) {
          missing++;
          if (missing<=5) printf("FACE-LIST (thread %i): face %i intersects the block [%e,%e]x[%e,%e]x[%e,%e] (level %i) and is not in its list\n",
              PIC::ThisThread,f,node->xmin[0],node->xmax[0],node->xmin[1],node->xmax[1],node->xmin[2],node->xmax[2],node->RefinmentLevel);
        }
      }
    }

    missing=SumLL(missing),listed=SumLL(listed);
    long long nLeaves=leaves.size();
    int maxLevel=0,minLevel=100;
    for (size_t n=0;n<leaves.size();n++) maxLevel=std::max(maxLevel,(int)leaves[n]->RefinmentLevel),minLevel=std::min(minLevel,(int)leaves[n]->RefinmentLevel);

    Check(missing==0,"face-lists","%lli faces missing from the cut-face lists of %lli blocks (levels %i..%i; %lli list entries)",missing,nLeaves,minLevel,maxLevel,listed);
  }

  //5. the movers, the surface interaction, the time steps
  if (MoverSelection!=0&&MoverSelection<4) AmbiguityLength=2.0*PIC::Mesh::mesh->EPS,HitTolerance=2.0*PIC::Mesh::mesh->EPS;
  //the gas is enclosed by the shell: a particle that reaches the boundary of the domain has passed through it
  if (MoverSelection==0) PIC::Mover::ForceFreeCutSurfaceDiagnostics::InsideBody=InsideScene,PIC::Mover::ForceFreeCutSurfaceDiagnostics::ReportDomainBoundary=true;
  PIC::Mover::SetUserDefinedParticleMover(Mover);
  PIC::Mover::ProcessTriangleCutFaceIntersection=SurfaceInteraction;
  Surface::Temperature::Isothremal::Temp=SurfaceTemperature;

  {
    const double m=PIC::MolecularData::GetMass(sDiffuse);
    dtDiffuse=FineCell/sqrt(2.0*Kbol*SurfaceTemperature/m);
  }

  for (int s=0;s<PIC::nTotalSpecies;s++) {
    PIC::ParticleWeightTimeStep::GlobalTimeStep[s]=(s==sDiffuse) ? dtDiffuse : 1.0;
    PIC::ParticleWeightTimeStep::GlobalParticleWeight[s]=1.0;
  }

  for (int s=0;s<nClasses;s++) HitCount[s].assign(CutCell::nBoundaryTriangleFaces,0);

  //the smallest block, the maximum displacement per step
  {
    std::vector<cNode*> leaves;
    CollectLeaves(PIC::Mesh::mesh->rootTree,leaves);
    double minBlock=HUGE_VAL;
    for (size_t n=0;n<leaves.size();n++) minBlock=std::min(minBlock,leaves[n]->xmax[0]-leaves[n]->xmin[0]);
    Info("mesh: %li blocks, smallest block %.4e m (%d cells of %.4e m), EPS=%.3e m; displacement per step: uniform %.3e, fast %.3e, absorbing %.3e, adversarial up to %.3e m",
        (long)leaves.size(),minBlock,_BLOCK_CELLS_X_,minBlock/_BLOCK_CELLS_X_,PIC::Mesh::mesh->EPS,SpeedUniform,SpeedFast,SpeedAbsorb,AdversarialSpeed[nAdversarialSpeeds-1]);
  }

  //6. the particles
  if (ReplayMode) AddParticle(ReplaySpec,ReplayX,ReplayV);
  else CreateParticles();

  long long count0[nClasses];
  CheckParticles(-1,count0);
  for (int s=0;s<nClasses;s++) count0[s]=SumLL(count0[s]);
  long long leakInitial=0;
  for (int s=0;s<nClasses;s++) leakInitial+=SumLL(LeakStrict[s]),LeakStrict[s]=0,LeakEPS[s]=0;
  MaxLeakDepth=0.0,nLeakReports=0;

  Check(leakInitial==0,"initial","%lli particles start inside a body (%li rejected before the creation); particles: uniform %lli, adversarial %lli, fast %lli, diffuse %lli, absorbing %lli",
      leakInitial,nGeneratedInside,count0[0],count0[1],count0[2],count0[3],count0[4]);

  //the gas volume: the inside of the container without the bodies (the cavity of the hollow cube is gas)
  double freeVolume=pow(2.0*ContainerInner,3);
  for (size_t b=0;b<Scene.Body.size();b++) if (Scene.Body[b].Name!="shell") freeVolume-=Scene.Body[b].Volume();

  //7. the steps
  PIC::Mover::ForceFreeCutSurfaceDiagnostics::Reset();
  double tStart=MPI_Wtime();
  long long count[nClasses];

  for (int step=0;step<nSteps;step++) {
    PIC::DomainBlockDecomposition::UpdateBlockTable();
    CollectReferenceParticles();
    RecordPreState();
    CurrentStep=step;

    if(MoverSelection>=4) {
      // Complete pending physical-time/MPI trajectories before evaluating
      // the existing independent brute-force reference. No extra full dt is
      // applied to particles received during continuation exchanges.
      double moving=0,exchange=0;PIC::TimeStepInternal::ExecutionTrackDefault(moving,exchange);
    }
    else PIC::Mover::MoveParticles();
    CompareWithReference();
    if(MoverSelection<4)PIC::Parallel::ExchangeParticleData();

    CheckParticles(step,count);

    if (ReplayMode&&(PIC::nTotalThreads==1)) {
      for (size_t ptr=0;ptr<PreState.size();ptr++) if ((PreState[ptr].spec>=0)&&PIC::ParticleBuffer::IsParticleAllocated(ptr)) {
        double x[3],v[3];
        PIC::ParticleBuffer::byte *data=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
        PIC::ParticleBuffer::GetX(x,data);
        PIC::ParticleBuffer::GetV(v,data);
        printf("REPLAY step %i: x=(%.17e,%.17e,%.17e) v=(%.17e,%.17e,%.17e) depth %.3e\n",step,x[0],x[1],x[2],v[0],v[1],v[2],LeakDepth(x));
      }
    }

    if ((PIC::ThisThread==0)&&((step%20==0)||(MoverSelection!=0))) {printf("step %i (%.1f s)\n",step,MPI_Wtime()-tStart); fflush(stdout);}
  }

  double tRun=MPI_Wtime()-tStart;
  for (int s=0;s<nClasses;s++) count[s]=SumLL(count[s]);

  //8. the results
  long long leak[nClasses],leakEPS[nClasses],leakTotal=0;
  for (int s=0;s<nClasses;s++) leak[s]=SumLL(LeakStrict[s]),leakEPS[s]=SumLL(LeakEPS[s]),leakTotal+=leak[s];
  double maxDepth=MaxD(MaxLeakDepth);

  PIC::Mover::ForceFreeCutSurfaceDiagnostics::cCounters c;
  PIC::Mover::ForceFreeCutSurfaceDiagnostics::Reduce(c);

  Check(leakTotal==0,"leaks","%lli particle-steps found inside a body deeper than %.1e m (%lli deeper than EPS=%.1e m; max depth %.3e m): uniform %lli, adversarial %lli, fast %lli, diffuse %lli, absorbing %lli",
      leakTotal,TolStrict,leakEPS[0]+leakEPS[1]+leakEPS[2]+leakEPS[3]+leakEPS[4],PIC::Mesh::mesh->EPS,maxDepth,leak[0],leak[1],leak[2],leak[3],leak[4]);

  Check((count[sUniform]==count0[sUniform])&&(count[sAdversarial]==count0[sAdversarial])&&(count[sFast]==count0[sFast])&&(c.Errors()==0),"lost",
      "specular particles: uniform %lli of %lli, adversarial %lli of %lli, fast %lli of %lli; mover errors %llu (not leaving %llu, too many interactions %llu, unallocated block %llu, outside %llu, walk %llu)",
      count[0],count0[0],count[1],count0[1],count[2],count0[2],c.Errors(),c.ErrorNotLeaving,c.ErrorTooManyInteractions,c.ErrorUnallocatedBlock,c.ErrorOutsideDomain,c.ErrorWalk);

  double speedError=MaxD(MaxSpeedError);
  Check(speedError<1.0E-12,"speed","max relative change of |v| of the specular particles %.3e",speedError);

  {
    long long bad=SumLL(BadIncoming),off=SumLL(OffPlane),offFace=SumLL(OffFace);
    double maxDist=MaxD(MaxPlaneDistance);
    Check((bad==0)&&(off==0)&&(offFace==0),"hits","%lli hits not moving into the body, %lli farther than %.1e m from the face plane (max %.3e m), %lli outside the closed face",
        bad,off,TolPlane,maxDist,offFace);
  }

  {
    long long compared=SumLL(RefCompared),ambiguous=SumLL(RefAmbiguous),mismatch=SumLL(RefMismatch);
    double maxDx=MaxD(MaxRefPositionDifference);
    Check((mismatch==0)&&(compared>0),"reference","%lli particle-steps compared with the brute-force tracer, %lli differ (position %.0e m, velocity 1e-9; max position difference %.3e m); %lli ambiguous steps not compared",
        compared,mismatch,RefPositionTolerance,maxDx,ambiguous);
  }

  //the hit rate of the uniform gas on each body
  {
    std::vector<unsigned long long> hits(HitCount[sUniform].size());
    MPI_Allreduce(HitCount[sUniform].data(),hits.data(),hits.size(),MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);

    const double n=count0[sUniform]/freeVolume,flux=n*SpeedUniform/4.0;   //per unit area and step
    bool ok=true;
    std::string s;
    double totalExpected=0.0,totalHits=0.0;
    int f=0;

    for (size_t b=0;b<Scene.Body.size();b++) {
      double h=0.0,area=0.0;
      for (int k=0;k<Scene.Body[b].nFaces();k++,f++) {
        h+=hits[f];
        if (FaceExposed(Scene.Body[b],k)) area+=Scene.Body[b].FaceArea(k);
      }
      double expected=flux*area*nSteps,dev=(h-expected)/sqrt(expected);
      if (fabs(dev)>5.0) ok=false;
      totalExpected+=expected,totalHits+=h;
      char buf[300];
      sprintf(buf,"%s%s %.0f/%.1f (%+.2f sigma)",(b==0) ? "" : ", ",Scene.Body[b].Name.c_str(),h,expected,dev);
      s+=buf;
    }

    double dev=(totalHits-totalExpected)/sqrt(totalExpected);
    if (fabs(dev)>5.0) ok=false;
    Check(ok,"hit-rate","uniform gas hits/expected n<v>/4*A*t: total %.0f/%.1f (%+.2f sigma, ratio %.4f); %s",totalHits,totalExpected,dev,totalHits/totalExpected,s.c_str());
  }

  //diffuse reflection
  {
    double events=SumD(DiffuseEvents),sum=SumD(DiffuseNormalSpeedSum),vmin=-MaxD(-DiffuseNormalSpeedMin);
    const double m=PIC::MolecularData::GetMass(sDiffuse),beta=sqrt(m/(2.0*Kbol*SurfaceTemperature));
    double expected=sqrt(M_PI)/(2.0*beta),mean=(events>0.0) ? sum/events : 0.0;
    double sigma=expected*sqrt((1.0-M_PI/4.0)/(M_PI/4.0))/sqrt(std::max(events,1.0));
    Check((events>100.0)&&(fabs(mean-expected)<5.0*sigma)&&(vmin>0.0),"diffuse","%.0f diffuse reflections: mean normal velocity %.4f m/s, expected sqrt(pi)/(2 beta)=%.4f m/s (%+.2f sigma); min %.3e m/s; particles %lli of %lli",
        events,mean,expected,(mean-expected)/sigma,vmin,count[sDiffuse],count0[sDiffuse]);
  }

  //absorption
  {
    long long events=SumLL(AbsorptionEvents);
    Check(count0[sAbsorb]-count[sAbsorb]==events+leak[sAbsorb],"absorbing","%lli absorption events, %lli of %lli particles left (%lli removed)",events,count[sAbsorb],count0[sAbsorb],count0[sAbsorb]-count[sAbsorb]);
  }

  double maxCorrection=MaxD(PIC::Mover::ForceFreeCutSurfaceDiagnostics::MaxHitPointCorrection);
  Info("mover %s: %.1f s for %i steps; moves %llu, segments %llu, blocks visited %llu, interactions %llu (%llu at the start of a segment), absorbed %llu, domain boundary %llu, segments starting on the surface %llu, hit points moved onto their face %llu (max %.3e m), inside before/after the move %llu/%llu",
      MoverName.c_str(),tRun,nSteps,c.Moves,c.Segments,c.BlocksVisited,c.Interactions,c.InteractionsAtSegmentStart,c.Absorbed,c.DomainBoundary,c.OnSurfaceStarts,c.HitPointsMovedOntoFace,maxCorrection,c.InsideBeforeMove,c.InsideAfterMove);

  sprintf(fname,"cut-surface-mover.%s.%s.np%i.digest",MoverName.c_str(),Mode.c_str(),PIC::nTotalThreads);
  WriteDigest(fname);

  if (PIC::ThisThread==0) {
    fprintf(Report,"%s %s %s: %i checks failed\n",(nFailed==0) ? "PASS" : "FAIL",MoverName.c_str(),Mode.c_str(),nFailed);
    fclose(Report);
  }

  MPI_Finalize();
  return (nFailed==0) ? 0 : 1;
}

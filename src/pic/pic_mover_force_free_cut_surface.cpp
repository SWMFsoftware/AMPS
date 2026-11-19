//$Id$
/*
 * PIC::Mover::ForceFreeCutSurface(): force-free (straight-line) particle mover whose purpose is the exact
 * interaction of particles with the triangulated internal surface (cut faces).
 *
 * In every time step the particle moves along the segment x+t*D, D=v*dt, 0<=t<=1.  The mover
 *  1. finds the EARLIEST crossing of the segment into the body: it walks through every block the segment
 *     crosses (any refinement level, also blocks that are neither the start nor the end block, also blocks of
 *     other processes: the tree and its cut-face lists are known on every process) and tests the cut faces
 *     listed in each of them with the watertight line/triangle test CutCell::Watertight::
 *     LineTriangleIntersection() (meshAMRcutcell_watertight.h).  Only crossings into the body count (v.n<0,
 *     n = the external normal of the face).  The walk stops as soon as the next block begins after the earliest
 *     crossing found so far, so a face listed in several blocks or hit in a block other than the one that
 *     lists it does not change the result;
 *  2. moves the particle to the crossing, applies the surface interaction through the usual callback
 *     PIC::Mover::ProcessTriangleCutFaceIntersection (specular/diffuse reflection, absorption, sampling of the
 *     hit are the business of the callback), and continues with the rest of the step: any number of
 *     interactions per step;
 *  3. handles the boundary of the computational domain as configured by
 *     _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE_ (delete, specular reflection, user function);
 *  4. puts the particle into the temporary list of the cell where the step ends (as the other movers: a block
 *     of another process is a ghost block of this one, and the particle exchange sends it).
 *
 * Robustness:
 *  - edges and vertices: the line/triangle decision is exact for the rounded coordinates and consistent
 *    between faces that share an edge or a vertex, so a crossing cannot fall between two faces; a crossing
 *    exactly on an edge or a vertex is reported by every face that contains it, the earliest (then the most
 *    head-on, then the lowest face index) is processed; at the point of an interaction another face counts
 *    only if part of it lies in front of the face just left (a concave edge or corner: several reflections at
 *    the same point), not at a convex edge or vertex (the particle leaves the body there);
 *  - the face just left is not hit again because only crossings INTO the body count: after a reflection the
 *    particle moves out of that face.  Nothing is excluded by identity and nothing is ignored for being
 *    close to the start: a particle that starts exactly on the surface (or up to BehindTolerance behind a face
 *    plane because of rounding) and moves into the body interacts at t=0;
 *  - a segment that starts on the surface (within BehindTolerance=1e-12 of the domain size of a face plane,
 *    over the face) is tested from a point BehindTolerance on the gas side: rounding can put such a start just
 *    behind the plane, and a particle sliding along the surface, or entering it at a grazing angle, would cross
 *    no face (at a concave edge it would pass beside the next face).  Only the tests start there;
 *  - grazing incidence: the decision whether the line passes through a face does not depend on the angle; a
 *    line in the plane of a face does not cross it.  The crossing point of a grazing line is ill conditioned
 *    along the face: a computed hit point beside the face is moved to the nearest point of the closed face
 *    (counted in HitPointsMovedOntoFace, the displacement in MaxHitPointCorrection: about 2e-12 m in the test);
 *  - very fast particles: the walk covers the whole segment, whatever the number of blocks; the step must end
 *    in a block allocated on this process (a local block or a ghost block), as for every AMPS mover;
 *  - nothing is fixed silently: a callback that does not send the particle out of the face (after
 *    MaxCallbackRetries attempts), too many interactions at one point or in one step, a step that ends in a
 *    block not allocated on this process are counted in ForceFreeCutSurfaceDiagnostics::Counters, reported
 *    (the first MaxReportedErrors of them) and the particle is removed.  With
 *    ForceFreeCutSurfaceDiagnostics::InsideBody set, the start and the end of every step are tested and a
 *    particle inside the body is counted and reported; it is never moved.
 *
 * The time step is the one of the start block for the whole step (uniform time step), forces are not applied.
 * Select it with "MoverIntegratorMode=ForceFreeCutSurface" in #main of the input file (ampsConfig.pl), or
 * define _PIC_PARTICLE_MOVER__MOVE_PARTICLE_TIME_STEP_(ptr,LocalTimeStep,node) as
 * PIC::Mover::ForceFreeCutSurface(ptr,LocalTimeStep,node).
 */

#include "pic.h"
#include "meshAMRcutcell_watertight.h"
#include "pic_mover_force_cut_surface_math.h"

#include <cfloat>
#include <cmath>
#include <cstdio>

namespace PIC {
namespace Mover {
namespace ForceFreeCutSurfaceDiagnostics {
  cCounters Counters;
  bool (*InsideBody)(const double *x)=NULL;
  int MaxReportedErrors=20;
  bool ReportDomainBoundary=false;
  double MaxHitPointCorrection=0.0;
  int MaxCallbackRetries=100;
  int MaxSamePointInteractions=32;
  long int MaxInteractionsPerStep=100000;

  void Reset() {Counters=cCounters(),MaxHitPointCorrection=0.0;}

  void Reduce(cCounters& global) {
    unsigned long long l[cCounters::nFields],g[cCounters::nFields];
    Counters.Pack(l);
    MPI_Allreduce(l,g,cCounters::nFields,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
    global.Unpack(g);
  }
}
}
}

//The mover and the crossing search work in three dimensions (block boxes, the cut-face lists, 3D positions); in the 1D
//and 2D builds, which compile this file too, they are stubs (see the end of the file).
#if DIM == 3

namespace {

typedef cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> cNode;
namespace Diag=PIC::Mover::ForceFreeCutSurfaceDiagnostics;

int nReportedErrors=0;

inline void Count(unsigned long long& c) {
#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  #pragma omp atomic
#endif
  c++;
}

void Report(const char *what,long int ptr,const double *x,const double *v,const double *x0=NULL,const double *v0=NULL) {
  bool print=false;

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  #pragma omp critical (ForceFreeCutSurfaceReport)
#endif
  {
    if (nReportedErrors<Diag::MaxReportedErrors) nReportedErrors++,print=true;
  }

  if (print) {
    if (x0==NULL) printf("$PREFIX: ForceFreeCutSurface (thread %i): %s: particle %ld x=(%.17e,%.17e,%.17e) v=(%.17e,%.17e,%.17e)\n",
        PIC::ThisThread,what,ptr,x[0],x[1],x[2],v[0],v[1],v[2]);
    else printf("$PREFIX: ForceFreeCutSurface (thread %i): %s: particle %ld x=(%.17e,%.17e,%.17e) v=(%.17e,%.17e,%.17e); the step started at x=(%.17e,%.17e,%.17e) v=(%.17e,%.17e,%.17e)\n",
        PIC::ThisThread,what,ptr,x[0],x[1],x[2],v[0],v[1],v[2],x0[0],x0[1],x0[2],v0[0],v0[1],v0[2]);
    fflush(stdout);
  }
}

//the distance behind a face plane up to which a start point counts as on the face (rounding of the
//coordinates of points on the surface); set from the size of the domain at the first call
double BehindTolerance=-1.0;

void InitTolerances() {
  if (BehindTolerance>0.0) return;

  double l=0.0;
  for (int idim=0;idim<3;idim++) l=std::max(l,std::max(std::fabs(PIC::Mesh::mesh->xGlobalMin[idim]),std::fabs(PIC::Mesh::mesh->xGlobalMax[idim])));
  for (int idim=0;idim<3;idim++) l=std::max(l,PIC::Mesh::mesh->xGlobalMax[idim]-PIC::Mesh::mesh->xGlobalMin[idim]);
  BehindTolerance=1.0E-12*l;
}

//a coordinate just below the upper boundary of the domain that findTreeNode places inside the domain: it works with
//floor((x-xGlobalMin)/dx), in which a point one rounding unit below xGlobalMax is rounded back onto the boundary (the
//particles reflected at the upper faces of the domain were removed as being outside of it: specular-domain test)
inline double UpperInside(int idim) {
  const double s=std::max(std::max(std::fabs(PIC::Mesh::mesh->xGlobalMin[idim]),std::fabs(PIC::Mesh::mesh->xGlobalMax[idim])),
      PIC::Mesh::mesh->xGlobalMax[idim]-PIC::Mesh::mesh->xGlobalMin[idim]);

  return PIC::Mesh::mesh->xGlobalMax[idim]-16.0*DBL_EPSILON*s;
}

//the leaf containing x; x on the upper boundary of the domain, or within a few rounding units of it, is looked up just
//inside (UpperInside)
cNode* LocateNode(const double *x,cNode *hint) {
  double xp[3]={x[0],x[1],x[2]};
  cNode *node=PIC::Mesh::mesh->findTreeNode(xp,hint);
  if (node!=NULL) return node;

  bool changed=false;
  for (int idim=0;idim<3;idim++) {
    if (xp[idim]>UpperInside(idim)) xp[idim]=UpperInside(idim),changed=true;
    if (xp[idim]<PIC::Mesh::mesh->xGlobalMin[idim]) xp[idim]=PIC::Mesh::mesh->xGlobalMin[idim],changed=true;
  }

  return (changed) ? PIC::Mesh::mesh->findTreeNode(xp,hint) : NULL;
}

struct cHit {
  double t,dn;
  CutCell::cTriangleFace *face;
};

inline bool Earlier(double t,double dn,CutCell::cTriangleFace *f,const cHit& hit) {
  if (hit.face==NULL) return true;
  if (t!=hit.t) return t<hit.t;
  if (dn!=hit.dn) return dn<hit.dn;
  return (f-CutCell::BoundaryTriangleFaces)<(hit.face-CutCell::BoundaryTriangleFaces);
}

enum eEvent {EventNone,EventFace,EventDomainBoundary,EventError};

//the earliest event of the segment x+t*D, 0<=t<=1: a crossing into the body (hit), or the exit from the
//domain (tDomain); 'node' contains x on input; on output it is the last block visited; 'lastFace' is the face
//of the interaction at x that ended the previous segment (NULL at the beginning of the step)
eEvent FindEarliestEvent(const double *xStart,const double *D,double lengthD,const CutCell::cTriangleFace *lastFace,cNode* &node,cHit& hit,double& tDomain,const double* physicalVelocity=nullptr) {
  const double tBehind=BehindTolerance/lengthD;
  cNode *cur=node;
  double tCur=0.0,x[3]={xStart[0],xStart[1],xStart[2]};

  hit.face=NULL,hit.t=HUGE_VAL,hit.dn=0.0;

  //A start point on the surface: within BehindTolerance of the plane of a face and over the face.  Rounding
  //can put such a point just behind the plane, and a particle that then moves along the surface (or into it at
  //a grazing angle) crosses no face: at a concave edge it would pass beside the next face into the body.  The
  //line tests of this segment therefore start BehindTolerance on the gas side of the surface (along the sum of
  //the normals of the faces the point is on).  Only the tests start there: the particle is not moved.
  {
    double N[3]={0.0,0.0,0.0};

    for (CutCell::cTriangleFaceDescriptor *t=cur->FirstTriangleCutFace;t!=NULL;t=t->next) {
      CutCell::cTriangleFace *f=t->TriangleFace;
      const double *n=f->ExternalNormal;
      const double sd=(x[0]-f->x0Face[0])*n[0]+(x[1]-f->x0Face[1])*n[1]+(x[2]-f->x0Face[2])*n[2];

      if (!(std::fabs(sd)<=BehindTolerance)) continue;

      double xLocal[2];
      const double beta=4.0*BehindTolerance/std::min(f->e0Length,f->e1Length);
      f->GetProjectedLocalCoordinates(xLocal,x);
      if ((xLocal[0]<-beta)||(xLocal[1]<-beta)||(xLocal[0]+xLocal[1]>1.0+beta)) continue;

      for (int idim=0;idim<3;idim++) N[idim]+=n[idim];
    }

    const double l=std::sqrt(N[0]*N[0]+N[1]*N[1]+N[2]*N[2]);
    if (l>1.0E-3) {
      for (int idim=0;idim<3;idim++) x[idim]+=BehindTolerance*N[idim]/l;
      Count(Diag::Counters.OnSurfaceStarts);
    }
  }

  for (long int nVisited=0;nVisited<10000000;nVisited++) {
    Count(Diag::Counters.BlocksVisited);

    for (CutCell::cTriangleFaceDescriptor *t=cur->FirstTriangleCutFace;t!=NULL;t=t->next) {
      CutCell::cTriangleFace *f=t->TriangleFace;
      const double* direction=physicalVelocity?physicalVelocity:D;
      const double dn=direction[0]*f->ExternalNormal[0]+direction[1]*f->ExternalNormal[1]+direction[2]*f->ExternalNormal[2];
      double tHit;

      if(physicalVelocity&&!PIC::Mover::ForceCutSurface::Math::IncomingVelocity(physicalVelocity,f->ExternalNormal))continue;
      if (!(dn<0.0)) continue;   //not moving into the body through this face
      if (!CutCell::Watertight::LineTriangleIntersection(x,D,f->x0Face,f->x1Face,f->x2Face,tHit)) continue;
      if ((tHit<-tBehind)||(tHit>1.0)) continue;
      if (tHit<0.0) tHit=0.0;

      if ((lastFace!=NULL)&&(f!=lastFace)&&(tHit*lengthD<=2.0*BehindTolerance)) {
        //A face met again at the point of the previous interaction (an edge or a vertex): the particle moves
        //into it only if the face lies (partly) in front of the plane of the face it has just left - a concave
        //edge or corner.  A face entirely behind that plane (a convex edge or vertex) cannot be entered: the
        //particle leaves the body through the face of the previous interaction.
        const double *nl=lastFace->ExternalNormal,*pl=lastFace->x0Face;
        const double *v[3]={f->x0Face,f->x1Face,f->x2Face};
        double front=-HUGE_VAL;
        for (int k=0;k<3;k++) front=std::max(front,(v[k][0]-pl[0])*nl[0]+(v[k][1]-pl[1])*nl[1]+(v[k][2]-pl[2])*nl[2]);
        if (!(front>BehindTolerance)) continue;
      }

      if (Earlier(tHit,dn,f,hit)) hit.t=tHit,hit.dn=dn,hit.face=f;
    }

    //the exit of the segment from this block
    double tExit=HUGE_VAL;
    for (int idim=0;idim<3;idim++) {
      double te;
      if (D[idim]>0.0) te=(cur->xmax[idim]-x[idim])/D[idim];
      else if (D[idim]<0.0) te=(cur->xmin[idim]-x[idim])/D[idim];
      else continue;
      if (te<tExit) tExit=te;
    }
    if (tExit<tCur) tExit=tCur;

    if ((hit.face!=NULL)&&(hit.t<=tExit)) {node=cur; return EventFace;}
    if (tExit>=1.0) {node=cur; return EventNone;}

    //the next block: probe just beyond the exit point (much less than the smallest block size)
    const double minBlock=PIC::Mesh::mesh->dxRootBlock[0]/(1<<_MAX_REFINMENT_LEVEL_);
    cNode *next=cur;
    double probeLength=1.0E-9*minBlock;

    for (int itry=0;(itry<12)&&(next==cur);itry++,probeLength*=10.0) {
      double tp=tExit+probeLength/lengthD,xp[3];
      for (int idim=0;idim<3;idim++) xp[idim]=x[idim]+tp*D[idim];
      next=PIC::Mesh::mesh->findTreeNode(xp,cur);
    }

    if (next==NULL) {
      //the segment leaves the domain: the exit parameter from the box of the domain
      tDomain=HUGE_VAL;
      for (int idim=0;idim<3;idim++) {
        double te;
        if (D[idim]>0.0) te=(PIC::Mesh::mesh->xGlobalMax[idim]-x[idim])/D[idim];
        else if (D[idim]<0.0) te=(PIC::Mesh::mesh->xGlobalMin[idim]-x[idim])/D[idim];
        else continue;
        if (te<tDomain) tDomain=te;
      }
      if (tDomain<0.0) tDomain=0.0;
      node=cur;
      return EventDomainBoundary;
    }

    if (next==cur) {node=cur; return EventError;}

    cur=next,tCur=tExit;
  }

  node=cur;
  return EventError;
}

//put the particle into the temporary list of its cell (as PIC::Mover::TrajectoryTrackingMover())
void AddToTempList(long int ptr,PIC::ParticleBuffer::byte *ParticleData,PIC::Mesh::cDataBlockAMR *block,int i,int j,int k) {
#if _PIC_MOVER__MPI_MULTITHREAD_ == _PIC_MODE_ON_
  PIC::ParticleBuffer::SetPrev(-1,ParticleData);

  long int tempFirstCellParticle=atomic_exchange(block->tempParticleMovingListTable+i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k),ptr);

  PIC::ParticleBuffer::SetNext(tempFirstCellParticle,ParticleData);
  if (tempFirstCellParticle!=-1) PIC::ParticleBuffer::SetPrev(ptr,tempFirstCellParticle);

#elif _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
  long int tempFirstCellParticle,*tempFirstCellParticlePtr;

  tempFirstCellParticlePtr=block->tempParticleMovingListTable+i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k);
  tempFirstCellParticle=(*tempFirstCellParticlePtr);

  PIC::ParticleBuffer::SetNext(tempFirstCellParticle,ParticleData);
  PIC::ParticleBuffer::SetPrev(-1,ParticleData);

  if (tempFirstCellParticle!=-1) PIC::ParticleBuffer::SetPrev(ptr,tempFirstCellParticle);
  *tempFirstCellParticlePtr=ptr;

#elif _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  PIC::Mesh::cDataBlockAMR::cTempParticleMovingListMultiThreadTable* ThreadTempParticleMovingData=block->GetTempParticleMovingListMultiThreadTable(omp_get_thread_num(),i,j,k);

  PIC::ParticleBuffer::SetNext(ThreadTempParticleMovingData->first,ParticleData);
  PIC::ParticleBuffer::SetPrev(-1,ParticleData);

  if (ThreadTempParticleMovingData->last==-1) ThreadTempParticleMovingData->last=ptr;
  if (ThreadTempParticleMovingData->first!=-1) PIC::ParticleBuffer::SetPrev(ptr,ThreadTempParticleMovingData->first);
  ThreadTempParticleMovingData->first=ptr;
#else
#error The option is unknown
#endif
}

//move x, if its projection onto the plane of the face is outside the closed face, to the nearest point of the
//face plus the distance of x from the plane; returns the displacement (0 when x is not moved)
double ClampToFace(double *x,CutCell::cTriangleFace *f) {
  double xLocal[2];
  f->GetProjectedLocalCoordinates(xLocal,x);
  if ((xLocal[0]>=0.0)&&(xLocal[1]>=0.0)&&(xLocal[0]+xLocal[1]<=1.0)) return 0.0;

  const double *a=f->x0Face,*b=f->x1Face,*c=f->x2Face,*n=f->ExternalNormal;
  double q[3];
  CutCell::Watertight::ClosestPointOnTriangle(x,a,b,c,q);

  const double sd=(x[0]-a[0])*n[0]+(x[1]-a[1])*n[1]+(x[2]-a[2])*n[2];
  double moved=0.0;
  for (int idim=0;idim<3;idim++) {
    double xn=q[idim]+sd*n[idim];
    moved+=(xn-x[idim])*(xn-x[idim]);
    x[idim]=xn;
  }
  return std::sqrt(moved);
}

int RemoveParticle(long int ptr) {
  PIC::ParticleBuffer::DeleteParticle(ptr);
  return _PARTICLE_LEFT_THE_DOMAIN_;
}

}  // namespace


//The earliest crossing into the body of the straight segment x0 -> x0+D (0<=t<=1) over the cut faces of all blocks the
//segment crosses, with the same rules as ForceFreeCutSurface() (see the description at the top); used by the other
//movers to verify the segments of their steps.
CutCell::cTriangleFace* PIC::Mover::CutSurface::FindFirstCrossing(const double *x0,const double *D,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node,
    const CutCell::cTriangleFace *lastFace,double& t) {
  return FindFirstCrossing(x0,D,node,lastFace,t,nullptr);
}

CutCell::cTriangleFace* PIC::Mover::CutSurface::FindFirstCrossing(const double *x0,const double *D,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node,
    const CutCell::cTriangleFace *lastFace,double& t,const double* physicalVelocity) {
  if ((CutCell::nBoundaryTriangleFaces==0)||(node==NULL)) return NULL;

  const double lengthD=std::sqrt(D[0]*D[0]+D[1]*D[1]+D[2]*D[2]);
  if (!(lengthD>0.0)||!std::isfinite(lengthD)) return NULL;

  InitTolerances();

  for (int idim=0;idim<3;idim++) if ((x0[idim]<node->xmin[idim])||(x0[idim]>node->xmax[idim])) {
    node=LocateNode(x0,node);
    break;
  }

  if (node==NULL) return NULL;

  cHit hit;
  double tDomain;

  if (FindEarliestEvent(x0,D,lengthD,lastFace,node,hit,tDomain,physicalVelocity)!=EventFace) return NULL;

  t=hit.t;
  return hit.face;
}


int PIC::Mover::ForceFreeCutSurface(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  PIC::ParticleBuffer::byte *ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  double x[3],v[3];
  int spec;

  PIC::ParticleBuffer::GetX(x,ParticleData);
  PIC::ParticleBuffer::GetV(v,ParticleData);
  spec=PIC::ParticleBuffer::GetI(ParticleData);

  const double xStart[3]={x[0],x[1],x[2]},vStart[3]={v[0],v[1],v[2]};

  InitTolerances();
  Count(Diag::Counters.Moves);

  if ((Diag::InsideBody!=NULL)&&Diag::InsideBody(x)) {
    Count(Diag::Counters.InsideBeforeMove);
    Report("the particle is inside the body at the beginning of the step",ptr,x,v);
  }

  //the block containing the particle (a point on a block face may be listed in the neighbor)
  cNode *node=startNode;
  for (int idim=0;idim<3;idim++) if ((x[idim]<node->xmin[idim])||(x[idim]>node->xmax[idim])) {node=LocateNode(x,startNode); break;}

  if (node==NULL) {
    Count(Diag::Counters.ErrorOutsideDomain);
    Report("the particle is outside of the domain",ptr,x,v);
    return RemoveParticle(ptr);
  }

  double timeLeft=dtTotal;
  long int nInteractions=0;
  int nSamePoint=0;
  CutCell::cTriangleFace *lastFace=NULL;

  while (timeLeft>0.0) {
    const double D[3]={v[0]*timeLeft,v[1]*timeLeft,v[2]*timeLeft};
    const double lengthD=std::sqrt(D[0]*D[0]+D[1]*D[1]+D[2]*D[2]);

    if (!(lengthD>0.0)) break;
    if (!std::isfinite(lengthD)) {
      Count(Diag::Counters.ErrorNonFinite);
      Report("non-finite velocity",ptr,x,v);
      return RemoveParticle(ptr);
    }

    Count(Diag::Counters.Segments);

    cHit hit;
    double tDomain=0.0;
    eEvent event=FindEarliestEvent(x,D,lengthD,lastFace,node,hit,tDomain);
    lastFace=NULL;

    if (event==EventNone) {
      for (int idim=0;idim<3;idim++) x[idim]+=D[idim];
      timeLeft=0.0;
      break;
    }

    if (event==EventError) {
      Count(Diag::Counters.ErrorWalk);
      Report("the walk through the blocks along the trajectory failed",ptr,x,v);
      return RemoveParticle(ptr);
    }

    if (event==EventFace) {
      CutCell::cTriangleFace *face=hit.face;

      for (int idim=0;idim<3;idim++) x[idim]+=hit.t*D[idim];
      timeLeft*=(1.0-hit.t);

      //The crossing of a grazing line with the plane is ill conditioned along the face: the decision that the
      //line passes through this face is exact, but the computed point can lie beside it.  Such a point is moved
      //to the nearest point of the closed face (keeping its distance from the plane): next to a concave edge a
      //point beside the face would be inside the body.
      {
        double moved=ClampToFace(x,face);
        if (moved>0.0) {
          Count(Diag::Counters.HitPointsMovedOntoFace);
#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
          #pragma omp critical (ForceFreeCutSurfaceMaxCorrection)
#endif
          if (moved>Diag::MaxHitPointCorrection) Diag::MaxHitPointCorrection=moved;
        }
      }

      Count(Diag::Counters.Interactions);
      if (hit.t==0.0) Count(Diag::Counters.InteractionsAtSegmentStart);

      if (++nInteractions>Diag::MaxInteractionsPerStep) {
        Count(Diag::Counters.ErrorTooManyInteractions);
        Report("too many surface interactions in one step",ptr,x,v);
        return RemoveParticle(ptr);
      }

      if (hit.t*lengthD<=BehindTolerance) {
        if (++nSamePoint>Diag::MaxSamePointInteractions) {
          Count(Diag::Counters.ErrorTooManyInteractions);
          Report("too many surface interactions at one point",ptr,x,v);
          return RemoveParticle(ptr);
        }
      }
      else nSamePoint=0;

      node=LocateNode(x,node);
      if ((node==NULL)||(node->block==NULL)) {
        Count(Diag::Counters.ErrorUnallocatedBlock);
        Report("the surface is hit in a block not allocated on this process (the time step is too long)",ptr,x,v);
        return RemoveParticle(ptr);
      }

      //the surface interaction: the callback returns the new velocity; it must leave the face
      int code,nTry=0;
      double vn;

      do {
        code=(PIC::Mover::ProcessTriangleCutFaceIntersection!=NULL) ? PIC::Mover::ProcessTriangleCutFaceIntersection(ptr,x,v,face,node) : _PARTICLE_DELETED_ON_THE_FACE_;

        if (code==_PARTICLE_DELETED_ON_THE_FACE_) {
          Count(Diag::Counters.Absorbed);
          return RemoveParticle(ptr);
        }

        vn=v[0]*face->ExternalNormal[0]+v[1]*face->ExternalNormal[1]+v[2]*face->ExternalNormal[2];
      }
      while ((!(vn>0.0))&&(++nTry<Diag::MaxCallbackRetries));

      if (!(vn>0.0)) {
        Count(Diag::Counters.ErrorNotLeaving);
        Report("the surface interaction does not send the particle out of the face",ptr,x,v);
        return RemoveParticle(ptr);
      }

      lastFace=face;
      continue;
    }

    //the boundary of the domain
    {
      if (tDomain>1.0) tDomain=1.0;

      bool crossed[3]={false,false,false};
      double tAxis[3];

      for (int idim=0;idim<3;idim++) {
        tAxis[idim]=HUGE_VAL;
        if (D[idim]>0.0) tAxis[idim]=(PIC::Mesh::mesh->xGlobalMax[idim]-x[idim])/D[idim];
        else if (D[idim]<0.0) tAxis[idim]=(PIC::Mesh::mesh->xGlobalMin[idim]-x[idim])/D[idim];
      }

      for (int idim=0;idim<3;idim++) {
        x[idim]+=tDomain*D[idim];
        if (tAxis[idim]<=tDomain*(1.0+1.0E-12)) {
          crossed[idim]=true;
          x[idim]=(D[idim]>0.0) ? PIC::Mesh::mesh->xGlobalMax[idim] : PIC::Mesh::mesh->xGlobalMin[idim];
        }
      }

      timeLeft*=(1.0-tDomain);
      Count(Diag::Counters.DomainBoundary);

#if _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE_ == _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__DELETE_
      if (Diag::ReportDomainBoundary) Report("the particle leaves the domain",ptr,x,v,xStart,vStart);
      return RemoveParticle(ptr);
#elif _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE_ == _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__SPECULAR_REFLECTION_
      for (int idim=0;idim<3;idim++) if (crossed[idim]) v[idim]=-v[idim];
#elif _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE_ == _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__USER_FUNCTION_
      {
        int nface=-1;
        for (int idim=0;idim<3;idim++) if (crossed[idim]) {nface=2*idim+((D[idim]>0.0) ? 1 : 0); break;}

        cNode *bnode=LocateNode(x,node);
        int code=(bnode!=NULL) ? ProcessOutsideDomainParticles(ptr,x,v,nface,bnode) : _PARTICLE_DELETED_ON_THE_FACE_;
        if (code==_PARTICLE_DELETED_ON_THE_FACE_) return RemoveParticle(ptr);
      }
#else
      exit(__LINE__,__FILE__,"Error: ForceFreeCutSurface supports the DELETE, SPECULAR_REFLECTION and USER_FUNCTION domain boundary modes");
#endif

      node=LocateNode(x,node);
      if (node==NULL) {
        Count(Diag::Counters.ErrorOutsideDomain);
        Report("the particle is outside of the domain after the interaction with its boundary",ptr,x,v);
        return RemoveParticle(ptr);
      }
    }
  }

  //the end of the step: the block and the cell
  cNode *newNode=LocateNode(x,node);

  if (newNode!=NULL) for (int idim=0;idim<3;idim++) {
    //a point on the upper boundary of the domain, or within a few rounding units of it, belongs to no block for
    //findTreeNode: move it just inside (UpperInside), so that findTreeNode finds it later
    if (x[idim]>UpperInside(idim)) x[idim]=UpperInside(idim);
  }

  if (newNode==NULL) {
    Count(Diag::Counters.ErrorOutsideDomain);
    Report("the particle is outside of the domain at the end of the step",ptr,x,v);
    return RemoveParticle(ptr);
  }

  if (newNode->block==NULL) {
    Count(Diag::Counters.ErrorUnallocatedBlock);
    Report("the step ends in a block not allocated on this process (the time step is too long)",ptr,x,v);
    return RemoveParticle(ptr);
  }

  if ((Diag::InsideBody!=NULL)&&Diag::InsideBody(x)) {
    Count(Diag::Counters.InsideAfterMove);
    Report("the particle is inside the body at the end of the step",ptr,x,v,xStart,vStart);
  }

  int i,j,k;
  if (PIC::Mesh::mesh->FindCellIndex(x,i,j,k,newNode,false)==-1) {
    Count(Diag::Counters.ErrorWalk);
    Report("cannot find the cell of the particle",ptr,x,v);
    return RemoveParticle(ptr);
  }

#if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
  PIC::ParticleTracker::RecordTrajectoryPoint(x,v,spec,ParticleData,(void*)newNode);

#if _PIC_PARTICLE_TRACKER__TRACKING_CONDITION_MODE__DYNAMICS_ == _PIC_MODE_ON_
  PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(x,v,spec,ParticleData,(void*)newNode);
#endif
#endif

  (void)spec;

  PIC::ParticleBuffer::SetV(v,ParticleData);
  PIC::ParticleBuffer::SetX(x,ParticleData);
  AddToTempList(ptr,ParticleData,newNode->block,i,j,k);

  return _PARTICLE_MOTION_FINISHED_;
}

#else //DIM == 3

CutCell::cTriangleFace* PIC::Mover::CutSurface::FindFirstCrossing(const double*,const double*,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*,const CutCell::cTriangleFace*,double&,const double*) {return NULL;}

CutCell::cTriangleFace* PIC::Mover::CutSurface::FindFirstCrossing(const double *x0,const double *D,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node,
    const CutCell::cTriangleFace *lastFace,double& t) {
  //there are no cut faces in the 1D and 2D meshes
  return NULL;
}

int PIC::Mover::ForceFreeCutSurface(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  exit(__LINE__,__FILE__,"Error: PIC::Mover::ForceFreeCutSurface is implemented for the 3D mesh only (DIM==3)");
  return _PARTICLE_LEFT_THE_DOMAIN_;
}

#endif //DIM == 3

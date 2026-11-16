//$Id$

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <limits>
#include <string.h>
#include <list>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <iostream>
#include <fstream>
#include <signal.h>

#include <sys/time.h>
#include <sys/resource.h>

using namespace std;

#include "constants.h"
#include "specfunc.h"
#include "ifileopr.h"
#include "meshAMRcutcell.h"

#include "pic.h"
#include "pic_geometry_kernels.h"

//init the ray tracking module
void PIC::RayTracing::Init() {

}

//calculate the exist point of the ray from a given block
bool PIC::RayTracing::GetBlockExitPoint(double *xBlockMin,double *xBlockMax,double *x0Ray,double *lRay,double *xBlockExit, double *xFaceExitLocal, int &nExitFace) {
  // The slab kernel selects the minimum forward face parameter.  It also
  // rejects zero directions and malformed boxes explicitly, avoiding the
  // historical collection of face-specific divisions.  The public boolean
  // contract is retained: only a geometrically valid forward exit is true.
  const PIC::GeometryKernels::RayBoxStatus status=
      PIC::GeometryKernels::GetRayBoxExit(
          xBlockMin,xBlockMax,x0Ray,lRay,PIC::Mesh::mesh->EPS,
          xBlockExit,xFaceExitLocal,nExitFace);
  return status==PIC::GeometryKernels::RayBoxStatus::Success;
}

/* Return true when no face of the surface triangulation intersects the ray
 * from xStart toward xTarget.  SetCutCellShadowAttribute() starts the ray at
 * the center of the face whose illumination it determines, and that face has
 * to be skipped (ExeptionFace): the origin lies on its plane, so whether
 * RayIntersection() reported the face itself as an obstacle depended only on
 * the sign of the rounding error of the origin-plane distance, which shadowed
 * about half of the sunlit faces at random; since RayIntersection() accepts
 * hits up to EPS behind the origin (as t=0) it shadowed all of them. */
bool PIC::RayTracing::TestDirectAccess(double *xStart,double *xTarget,void* ExeptionFace) {
	double c=0.0,l[3],x[3];
	int idim;

	//clculate the pointing vector
	for (idim=0;idim<3;idim++) {
		l[idim]=xTarget[idim]-xStart[idim];
		c+=l[idim]*l[idim];
    x[idim]=xStart[idim];
	}

	c=sqrt(c);
	for (idim=0;idim<3;idim++) l[idim]/=c;

	//determine the initial block
	cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=NULL;

	node=PIC::Mesh::mesh->findTreeNode(xStart);
	if (node==NULL) return false;

	//increment the call counter
  PIC::Mesh::IrregularSurface::CutFaceAccessCounter::IncrementCounter();

	//trace the ray
	PIC::Mesh::IrregularSurface::cTriangleFaceDescriptor *t;

	while (node!=NULL) {
	  PIC::Mesh::IrregularSurface::cTriangleFace *TriangleFace;
		double IntersectionTime;

		if ((t=node->FirstTriangleCutFace)!=NULL) {
      for (;t!=NULL;t=t->next) {
        TriangleFace=t->TriangleFace;

        if (PIC::Mesh::IrregularSurface::CutFaceAccessCounter::IsFirstAccecssWithAccessCounterUpdate(TriangleFace)==true) {
          if ((TriangleFace!=(PIC::Mesh::IrregularSurface::cTriangleFace*)ExeptionFace)&&(TriangleFace->RayIntersection(x,l,IntersectionTime,PIC::Mesh::mesh->EPS)==true)) return false;
        }
      }
		}

		//determine the new node
		double xNodeExit[3],xFaceExitLocal[2];
		int nExitFace,iFace,jFace;

		if (PIC::RayTracing::GetBlockExitPoint(node->xmin,node->xmax,x,l,xNodeExit,xFaceExitLocal,nExitFace)==true) {
      iFace=(xFaceExitLocal[0]<0.5) ? 0 : 1;
		  jFace=(xFaceExitLocal[1]<0.5) ? 0 : 1;

		  node=node->GetNeibFace(nExitFace,iFace,jFace,PIC::Mesh::mesh);
		}
		else {
		  PIC::RayTracing::GetBlockExitPoint(node->xmin,node->xmax,x,l,xNodeExit,xFaceExitLocal,nExitFace);
		  node=PIC::Mesh::mesh->findTreeNode(xNodeExit,node);
		}

		memcpy(x,xNodeExit,3*sizeof(double));
	}

	return true;
}

/* Count the crossings of the surface triangulation by the ray xStart+s*l,
 * s>=0, l=(xTarget-xStart)/|xTarget-xStart|; callers use the parity of the
 * count (CheckPointInsideDomain(), the orientation of the external normals in
 * InitExternalNormalVector()).
 *
 * cTriangleFace::RayIntersection() accepts a hit up to EPS outside the edges
 * of a triangle (and up to EPS behind the origin), so a ray that crosses the
 * surface within EPS of an edge shared by two triangles is reported by both of
 * them.  Counting each report flipped the parity for such rays: some external
 * normals were oriented into the body (and particles injected from those faces
 * started inside the nucleus), and points next to the surface were classified
 * on the wrong side.  As in the classification of the cut-cell volume
 * integration, the path length s of every hit is therefore recorded and hits
 * that coincide within the geometric tolerance count as one crossing. */
int PIC::RayTracing::CountFaceIntersectionNumber(double *xStart,double *xTarget,int MeshFileID,bool ParallelCheck,void* ExeptionFace) {
  double c=0.0,l[3],x[3];
  int cnt=-1,idim,IntersectionCounter=0;
  std::vector<double> CrossingPathLength; //path length from xStart of every reported hit
  double PathLengthToBlockEntry=0.0; //path length from xStart to 'x'


  //clculate the pointing vector
  for (idim=0;idim<3;idim++) {
    l[idim]=xTarget[idim]-xStart[idim];
    c+=l[idim]*l[idim];
    x[idim]=xStart[idim];
  }

  c=sqrt(c);
  for (idim=0;idim<3;idim++) l[idim]/=c;

  //determine the initial block
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=NULL;

  node=PIC::Mesh::mesh->findTreeNode(xStart);
  if (node==NULL) return false;

  //increment the call counter
  PIC::Mesh::IrregularSurface::CutFaceAccessCounter::IncrementCounter();

  //trace the ray
  PIC::Mesh::IrregularSurface::cTriangleFaceDescriptor *t;

  while (node!=NULL) {
    PIC::Mesh::IrregularSurface::cTriangleFace *TriangleFace;
    double IntersectionTime;
    bool code;

    if ((t=node->FirstTriangleCutFace)!=NULL) {
      for (;t!=NULL;t=t->next) {
        TriangleFace=t->TriangleFace;

        if (PIC::Mesh::IrregularSurface::CutFaceAccessCounter::IsFirstAccecssWithAccessCounterUpdate(TriangleFace)==true) {
          if ( ((MeshFileID<0)||(TriangleFace->MeshFileID==MeshFileID)) && ((ParallelCheck==false)||((++cnt)%PIC::nTotalThreads==PIC::ThisThread)) ) {
            code=TriangleFace->RayIntersection(x,l,IntersectionTime,PIC::Mesh::mesh->EPS);
            if ((TriangleFace!=(PIC::Mesh::IrregularSurface::cTriangleFace*)ExeptionFace) && (code==true)) CrossingPathLength.push_back(PathLengthToBlockEntry+IntersectionTime);
          }
        }
      }
    }

    //determine the new node
    double xNodeExit[3],xFaceExitLocal[2];
    int nExitFace,iFace,jFace;

    if (PIC::RayTracing::GetBlockExitPoint(node->xmin,node->xmax,x,l,xNodeExit,xFaceExitLocal,nExitFace)==true) {
      iFace=(xFaceExitLocal[0]<0.5) ? 0 : 1;
      jFace=(xFaceExitLocal[1]<0.5) ? 0 : 1;

      node=node->GetNeibFace(nExitFace,iFace,jFace,PIC::Mesh::mesh);
    }
    else {
      PIC::RayTracing::GetBlockExitPoint(node->xmin,node->xmax,x,l,xNodeExit,xFaceExitLocal,nExitFace);
      node=PIC::Mesh::mesh->findTreeNode(xNodeExit,node);
    }

    for (idim=0;idim<3;idim++) PathLengthToBlockEntry+=(xNodeExit[idim]-x[idim])*l[idim];
    memcpy(x,xNodeExit,3*sizeof(double));
  }

  //with ParallelCheck each process has tested a part of the faces: collect all hits
  if (ParallelCheck==true) {
    int nLocal=CrossingPathLength.size(),nTotal=0;
    std::vector<int> nHits(PIC::nTotalThreads),offset(PIC::nTotalThreads);

    MPI_Allgather(&nLocal,1,MPI_INT,nHits.data(),1,MPI_INT,MPI_GLOBAL_COMMUNICATOR);
    for (int thread=0;thread<PIC::nTotalThreads;thread++) offset[thread]=nTotal,nTotal+=nHits[thread];

    std::vector<double> all(nTotal>0 ? nTotal : 1);
    MPI_Allgatherv(CrossingPathLength.data(),nLocal,MPI_DOUBLE,all.data(),nHits.data(),offset.data(),MPI_DOUBLE,MPI_GLOBAL_COMMUNICATOR);
    all.resize(nTotal);
    CrossingPathLength.swap(all);
  }

  //hits closer than the tolerance of RayIntersection() along the ray are one crossing
  std::sort(CrossingPathLength.begin(),CrossingPathLength.end());

  const double tolerance=2.0*fabs(PIC::Mesh::mesh->EPS);
  double previous=0.0;

  for (std::vector<double>::const_iterator it=CrossingPathLength.begin();it!=CrossingPathLength.end();it++) {
    const double scale=std::max(1.0,fabs(*it));

    if ((IntersectionCounter==0)||(fabs(*it-previous)>std::max(tolerance,128.0*std::numeric_limits<double>::epsilon()*scale))) {
      previous=*it;
      IntersectionCounter++;
    }
  }

  return IntersectionCounter;
}

int PIC::RayTracing::FindFistIntersectedFace(double *x0Ray,double *lRay,double *xIntersection,bool ParallelCheck,void* ExeptionFace) {
  double x[3];
  int cnt=-1;

  //determine the initial block
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=NULL;

  //in case ParallelCheck==true -> broad casr the inital point and the search direction from the root
  if (ParallelCheck==true) {
    MPI_Bcast(x0Ray,3,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);
    MPI_Bcast(lRay,3,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);
  }

  //find the node
  node=PIC::Mesh::mesh->findTreeNode(x0Ray);
  if (node==NULL) return false;

  //increment the call counter
  PIC::Mesh::IrregularSurface::CutFaceAccessCounter::IncrementCounter();

  memcpy(x,x0Ray,3*sizeof(double));

  //trace the ray
  PIC::Mesh::IrregularSurface::cTriangleFaceDescriptor *t;
  PIC::Mesh::IrregularSurface::cTriangleFace *minIntersectionTimeTriangleFace=NULL;
  PIC::Mesh::IrregularSurface::cTriangleFace *TriangleFace;
  double IntersectionTime,minIntersectionTime=-1.0;
  double xIntersectionTemp[3],xMinTimeIntersection[3];

  while ((node!=NULL)&&(minIntersectionTimeTriangleFace==NULL)) {
    bool code;

    if ((t=node->FirstTriangleCutFace)!=NULL) {
      for (;t!=NULL;t=t->next) {
        TriangleFace=t->TriangleFace;

        if ((ParallelCheck==false)||((++cnt)%PIC::nTotalThreads==PIC::ThisThread)) {
          if (PIC::Mesh::IrregularSurface::CutFaceAccessCounter::IsFirstAccecssWithAccessCounterUpdate(TriangleFace)==true) {
            code=TriangleFace->RayIntersection(x0Ray,lRay,IntersectionTime,xIntersectionTemp,PIC::Mesh::mesh->EPS);

            if ((TriangleFace!=(PIC::Mesh::IrregularSurface::cTriangleFace*)ExeptionFace) && (code==true)) {
              if ((minIntersectionTimeTriangleFace==NULL) || (IntersectionTime>minIntersectionTime)) {
                minIntersectionTimeTriangleFace=TriangleFace;
                minIntersectionTime=IntersectionTime;
                memcpy(xMinTimeIntersection,xIntersectionTemp,3*sizeof(double));
              }
            }
          }
        }

      }
    }

    //determine the new node
    double xNodeExit[3],xFaceExitLocal[2];
    int nExitFace,iFace,jFace;

    if (PIC::RayTracing::GetBlockExitPoint(node->xmin,node->xmax,x,lRay,xNodeExit,xFaceExitLocal,nExitFace)==true) {
      iFace=(xFaceExitLocal[0]<0.5) ? 0 : 1;
      jFace=(xFaceExitLocal[1]<0.5) ? 0 : 1;

      node=node->GetNeibFace(nExitFace,iFace,jFace,PIC::Mesh::mesh);
    }
    else {
      PIC::RayTracing::GetBlockExitPoint(node->xmin,node->xmax,x,lRay,xNodeExit,xFaceExitLocal,nExitFace);
      node=PIC::Mesh::mesh->findTreeNode(xNodeExit,node);
    }

    memcpy(x,xNodeExit,3*sizeof(double));
  }

  //collect min intersection data from all processors
  if (ParallelCheck==true) {
    double IntersectionTimeTable[PIC::nTotalThreads];
    int iIntersectedFace,MinIntersectionThread=0;

    MPI_Gather(&minIntersectionTime,1,MPI_DOUBLE,IntersectionTimeTable,1,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);

    if (PIC::ThisThread==0) {
      for (int thread=0;thread<PIC::nTotalThreads;thread++) if ((minIntersectionTime<0.0) || ((IntersectionTimeTable[thread]>0.0)&&(IntersectionTimeTable[thread]<minIntersectionTime)) ) {
        minIntersectionTime=IntersectionTimeTable[thread],MinIntersectionThread=thread;
      }
    }

    MPI_Bcast(&MinIntersectionThread,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);
    MPI_Bcast(&minIntersectionTime,1,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);

    if (minIntersectionTime>0.0) {
      if (PIC::ThisThread==MinIntersectionThread) iIntersectedFace=minIntersectionTimeTriangleFace-PIC::Mesh::IrregularSurface::BoundaryTriangleFaces;

      MPI_Bcast(&iIntersectedFace,1,MPI_INT,MinIntersectionThread,MPI_GLOBAL_COMMUNICATOR);
      MPI_Bcast(xMinTimeIntersection,3,MPI_DOUBLE,MinIntersectionThread,MPI_GLOBAL_COMMUNICATOR);

      if (PIC::ThisThread!=MinIntersectionThread) minIntersectionTimeTriangleFace=iIntersectedFace+PIC::Mesh::IrregularSurface::BoundaryTriangleFaces;
    }
  }

  if (minIntersectionTimeTriangleFace!=NULL) {
    //intersection face is found -> return it
    memcpy(xIntersection,xMinTimeIntersection,3*sizeof(double));
    return minIntersectionTimeTriangleFace-PIC::Mesh::IrregularSurface::BoundaryTriangleFaces;
  }

  //no face of intersection was found -> return the default value
  return -1;
}

void PIC::RayTracing::FlushtCutCellShadowAttribute(int at) {
  for (int i=0;i<PIC::Mesh::IrregularSurface::nBoundaryTriangleFaces;i++) {
    auto face=PIC::Mesh::IrregularSurface::BoundaryTriangleFaces+i;

    PIC::Mesh::IrregularSurface::BoundaryTriangleFaces[i].pic__shadow_attribute=at;
  }
}  


void PIC::RayTracing::SetCutCellShadowAttribute(double *xLightSource,bool ParallelExecution) {
  double xTriangle[3],c;
  int i,idim,iStart,iFinish,nFaceThread;

  if (PIC::ThisThread==0) {
    printf("$PREFIX: Setting the cell shadow attribute.... ");
    fflush(stdout);
  }

  //init the ray traciing module if needed
  PIC::RayTracing::Init();

  if (ParallelExecution==true) {
    nFaceThread=PIC::Mesh::IrregularSurface::nBoundaryTriangleFaces/PIC::nTotalThreads;
    iStart=nFaceThread*PIC::ThisThread;
    iFinish=iStart+nFaceThread;
    if (PIC::ThisThread==PIC::nTotalThreads-1) iFinish=PIC::Mesh::IrregularSurface::nBoundaryTriangleFaces;
  }
  else iStart=0,iFinish=PIC::Mesh::IrregularSurface::nBoundaryTriangleFaces;


#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  #pragma omp parallel for schedule(dynamic,1) default (none) private (i,c,idim,xTriangle) shared (iStart,iFinish,PIC::Mesh::IrregularSurface::BoundaryTriangleFaces,xLightSource)
#endif
  for (i=iStart;i<iFinish;i++) {
    PIC::Mesh::IrregularSurface::BoundaryTriangleFaces[i].GetCenterPosition(xTriangle);
    for (c=0.0,idim=0;idim<3;idim++) c+=(xLightSource[idim]-xTriangle[idim])*PIC::Mesh::IrregularSurface::BoundaryTriangleFaces[i].ExternalNormal[idim];

    if (c<0.0) {
      //the point is in the shadow
      PIC::Mesh::IrregularSurface::BoundaryTriangleFaces[i].pic__shadow_attribute=_PIC__CUT_FACE_SHADOW_ATTRIBUTE__TRUE_;
    }
    else {
      //trace the ray to the light source
      if (TestDirectAccess(xTriangle,xLightSource,PIC::Mesh::IrregularSurface::BoundaryTriangleFaces+i)==false) {
        PIC::Mesh::IrregularSurface::BoundaryTriangleFaces[i].pic__shadow_attribute=_PIC__CUT_FACE_SHADOW_ATTRIBUTE__TRUE_;
      }
      else {
        PIC::Mesh::IrregularSurface::BoundaryTriangleFaces[i].pic__shadow_attribute=_PIC__CUT_FACE_SHADOW_ATTRIBUTE__FALSE_;
      }
    }
  }

  //collect the attributes when the search is performed in parallel model
  if (ParallelExecution==true) {
    char sendBuffer[2*nFaceThread];
    int thread,cnt;

    for (thread=0;thread<nTotalThreads;thread++) {
      iStart=nFaceThread*thread;
      iFinish=iStart+nFaceThread;
      if (thread==nTotalThreads-1) iFinish=PIC::Mesh::IrregularSurface::nBoundaryTriangleFaces;

      if (thread==PIC::ThisThread) {
        for (i=iStart,cnt=0;i<iFinish;i++,cnt++) sendBuffer[cnt]=PIC::Mesh::IrregularSurface::BoundaryTriangleFaces[i].pic__shadow_attribute;;
      }

      MPI_Bcast(sendBuffer,iFinish-iStart,MPI_CHAR,thread,MPI_GLOBAL_COMMUNICATOR);

      if (thread!=PIC::ThisThread) {
        for (i=iStart,cnt=0;i<iFinish;i++,cnt++) PIC::Mesh::IrregularSurface::BoundaryTriangleFaces[i].pic__shadow_attribute=sendBuffer[cnt];
      }
    }
  }


  //calculate the cosine of the face illumination angle
  for (i=0;i<PIC::Mesh::IrregularSurface::nBoundaryTriangleFaces;i++) {
    double xCenter[3],l=0.0,c=0.0,t;

    PIC::Mesh::IrregularSurface::BoundaryTriangleFaces[i].GetCenterPosition(xCenter);

    for (int idim=0;idim<3;idim++) {
      t=xLightSource[idim]-xCenter[idim];

      l+=t*t;
      c+=t*PIC::Mesh::IrregularSurface::BoundaryTriangleFaces[i].ExternalNormal[idim];
    }

    PIC::Mesh::IrregularSurface::BoundaryTriangleFaces[i].pic__cosine_illumination_angle=c/sqrt(l);
  }

  //the function return diagnostic message
  if (PIC::ThisThread==0) {
    printf("done\n");
    fflush(stdout);
  }
}

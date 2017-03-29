//$Id$

#include <stdio.h>
#include <stdlib.h>
#include <vector>
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

char *PIC::RayTracing::FaceAccessCounterTable=NULL;
char *PIC::RayTracing::FaceRayTracingOperationIDTable=NULL;

//init the ray tracking module
void PIC::RayTracing::Init() {
  int nThreadsOpenMP=1;

  static int AllocatedTableLength=-1;

  //get the number of the OpenMP threads


  if (AllocatedTableLength!=CutCell::nBoundaryTriangleFaces) {
    //reallocate the table

    if (FaceAccessCounterTable!=NULL) {
      delete [] FaceRayTracingOperationIDTable;
      delete [] FaceAccessCounterTable;
    }

    #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
    #pragma omp parallel shared(nThreadsOpenMP)
    {
      #pragma omp single
     {
       nThreadsOpenMP=omp_get_num_threads();
     }
   }
   #endif //_COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_

    AllocatedTableLength=CutCell::nBoundaryTriangleFaces;

    FaceAccessCounterTable=new char [nThreadsOpenMP];
    FaceRayTracingOperationIDTable=new char [nThreadsOpenMP*CutCell::nBoundaryTriangleFaces];

    int i,imax=nThreadsOpenMP*CutCell::nBoundaryTriangleFaces;
    for (i=0;i<imax;i++) FaceRayTracingOperationIDTable[i]=0;

    for (int iface=0;iface<CutCell::nBoundaryTriangleFaces;iface++) CutCell::BoundaryTriangleFaces[iface].pic__RayTracingOperationCounterTable=FaceRayTracingOperationIDTable+iface*nThreadsOpenMP;
  }
}

//reset the face access counter
void PIC::RayTracing::FlushFaceRayTracingOperationIDTable(int iThreadOpenMP) {
  int iface;

  for (iface=0;iface<CutCell::nBoundaryTriangleFaces;iface++) {
    CutCell::BoundaryTriangleFaces[iface].pic__RayTracingOperationCounterTable[iThreadOpenMP]=(char)0;
  }
}

//calculate the exist point of the ray from a given block
bool PIC::RayTracing::GetBlockExitPoint(double *xBlockMin,double *xBlockMax,double *x0Ray,double *lRay,double *xBlockExit, double *xFaceExitLocal, int &nExitFace) {
  double e0,e1,e2,dt,dtExit=-1.0;
  double xFace[2];

  e0=xBlockMax[0]-xBlockMin[0];
  e1=xBlockMax[1]-xBlockMin[1];
  e2=xBlockMax[2]-xBlockMin[2];

  nExitFace=-1;

  //calculate intersection of the ray with the faces
  //face 0
  if (lRay[0]<0.0) {
     dt=(xBlockMin[0]-x0Ray[0])/lRay[0];

     if (dt>dtExit) {
       xFace[0]=(x0Ray[1]+dt*lRay[1]-xBlockMin[1]); // /e1;
       xFace[1]=(x0Ray[2]+dt*lRay[2]-xBlockMin[2]); // /e2;

       //if ((xFace[0]>=0.0-PIC::Mesh::mesh.EPS)&&(xFace[0]<=1.0+PIC::Mesh::mesh.EPS)  && (xFace[1]>=0.0-PIC::Mesh::mesh.EPS)&&(xFace[1]<=1.0+PIC::Mesh::mesh.EPS)) {
       if ((xFace[0]>=0.0-PIC::Mesh::mesh.EPS)&&(xFace[0]<=e1+PIC::Mesh::mesh.EPS)  && (xFace[1]>=0.0-PIC::Mesh::mesh.EPS)&&(xFace[1]<=e2+PIC::Mesh::mesh.EPS)) {
         dtExit=dt,nExitFace=0;

         xFaceExitLocal[0]=xFace[0]/e1;
         xFaceExitLocal[1]=xFace[1]/e2;

         //memcpy(xFaceExitLocal,xFace,2*sizeof(double));
       }
    }
  }

  //face 1
  if (lRay[0]>0.0) {
    dt=(xBlockMax[0]-x0Ray[0])/lRay[0];

    if (dt>dtExit) {
      xFace[0]=(x0Ray[1]+dt*lRay[1]-xBlockMin[1]); // /e1;
      xFace[1]=(x0Ray[2]+dt*lRay[2]-xBlockMin[2]); // /e2;

      if ((xFace[0]>=0.0-PIC::Mesh::mesh.EPS)&&(xFace[0]<=e1+PIC::Mesh::mesh.EPS)  && (xFace[1]>=0.0-PIC::Mesh::mesh.EPS)&&(xFace[1]<=e2+PIC::Mesh::mesh.EPS)) {
        dtExit=dt,nExitFace=1;

        xFaceExitLocal[0]=xFace[0]/e1;
        xFaceExitLocal[1]=xFace[1]/e2;
        //memcpy(xFaceExitLocal,xFace,2*sizeof(double));
      }
    }
  }


  //face 2
  if (lRay[1]<0.0) {
    dt=(xBlockMin[1]-x0Ray[1])/lRay[1];

    if (dt>dtExit) {
      xFace[0]=(x0Ray[0]+dt*lRay[0]-xBlockMin[0]); // /e0;
      xFace[1]=(x0Ray[2]+dt*lRay[2]-xBlockMin[2]); // /e2;

      if ((xFace[0]>=0.0-PIC::Mesh::mesh.EPS)&&(xFace[0]<=e0+PIC::Mesh::mesh.EPS)  && (xFace[1]>=0.0-PIC::Mesh::mesh.EPS)&&(xFace[1]<=e2+PIC::Mesh::mesh.EPS)) {
        dtExit=dt,nExitFace=2;

        xFaceExitLocal[0]=xFace[0]/e0;
        xFaceExitLocal[1]=xFace[1]/e2;
        //memcpy(xFaceExitLocal,xFace,2*sizeof(double));
      }
    }
  }

  //face 3
  if (lRay[1]>0.0) {
    dt=(xBlockMax[1]-x0Ray[1])/lRay[1];

    if (dt>dtExit) {
      xFace[0]=(x0Ray[0]+dt*lRay[0]-xBlockMin[0]); // /e0;
      xFace[1]=(x0Ray[2]+dt*lRay[2]-xBlockMin[2]); // /e2;

      if ((xFace[0]>=0.0-PIC::Mesh::mesh.EPS)&&(xFace[0]<=e0+PIC::Mesh::mesh.EPS)  && (xFace[1]>=0.0-PIC::Mesh::mesh.EPS)&&(xFace[1]<=e2+PIC::Mesh::mesh.EPS)) {
        dtExit=dt,nExitFace=3;

        xFaceExitLocal[0]=xFace[0]/e0;
        xFaceExitLocal[1]=xFace[1]/e2;
        //memcpy(xFaceExitLocal,xFace,2*sizeof(double));
      }
    }
  }


  //face 4
  if (lRay[2]<0.0) {
    dt=(xBlockMin[2]-x0Ray[2])/lRay[2];

    if (dt>dtExit) {
      xFace[0]=(x0Ray[0]+dt*lRay[0]-xBlockMin[0]); // /e0;
      xFace[1]=(x0Ray[1]+dt*lRay[1]-xBlockMin[1]); // /e1;

      if ((xFace[0]>=0.0-PIC::Mesh::mesh.EPS)&&(xFace[0]<=e0+PIC::Mesh::mesh.EPS)  && (xFace[1]>=0.0-PIC::Mesh::mesh.EPS)&&(xFace[1]<=e1+PIC::Mesh::mesh.EPS)) {
        dtExit=dt,nExitFace=4;

        xFaceExitLocal[0]=xFace[0]/e0;
        xFaceExitLocal[1]=xFace[1]/e1;
        //memcpy(xFaceExitLocal,xFace,2*sizeof(double));
      }
    }
  }

  //face 5
  if (lRay[2]>0.0) {
    dt=(xBlockMax[2]-x0Ray[2])/lRay[2];

    if (dt>dtExit) {
      xFace[0]=(x0Ray[0]+dt*lRay[0]-xBlockMin[0]); // /e0;
      xFace[1]=(x0Ray[1]+dt*lRay[1]-xBlockMin[1]); // /e1;

      if ((xFace[0]>=0.0-PIC::Mesh::mesh.EPS)&&(xFace[0]<=e0+PIC::Mesh::mesh.EPS)  && (xFace[1]>=0.0-PIC::Mesh::mesh.EPS)&&(xFace[1]<=e1+PIC::Mesh::mesh.EPS)) {
        dtExit=dt,nExitFace=5;

        xFaceExitLocal[0]=xFace[0]/e0;
        xFaceExitLocal[1]=xFace[1]/e1;
        //memcpy(xFaceExitLocal,xFace,2*sizeof(double));
      }
    }
  }

  if (nExitFace==-1) return false;

  //calculate the exit point from the block



  xBlockExit[0]=x0Ray[0]+dtExit*lRay[0];
  xBlockExit[1]=x0Ray[1]+dtExit*lRay[1];
  xBlockExit[2]=x0Ray[2]+dtExit*lRay[2];

  return true;
}

bool PIC::RayTracing::TestDirectAccess(double *xStart,double *xTarget) {
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

	node=PIC::Mesh::mesh.findTreeNode(xStart);
	if (node==NULL) return false;

	//determine the thread number
	int thread=0;

  #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  thread=omp_get_thread_num();
  #endif

	//increment the call counter
  FaceAccessCounterTable[thread]++;

	if (FaceAccessCounterTable[thread]==(char)0xff) {
		//if the counter is 'owerflow' than re-init the counter value and the counters of the cut-faces
	  FaceAccessCounterTable[thread]=(char)1;
		FlushFaceRayTracingOperationIDTable(thread);
	}

	//trace the ray
	PIC::Mesh::IrregularSurface::cTriangleFaceDescriptor *t;

	while (node!=NULL) {
	  PIC::Mesh::IrregularSurface::cTriangleFace *TriangleFace;
		double IntersectionTime;
		bool code;

		if ((t=node->FirstTriangleCutFace)!=NULL) {
      for (;t!=NULL;t=t->next) {
        if ((TriangleFace=t->TriangleFace)->pic__RayTracingOperationCounterTable[thread]!=FaceAccessCounterTable[thread]) {
          code=TriangleFace->RayIntersection(x,l,IntersectionTime,PIC::Mesh::mesh.EPS);

          if ((code==true)&&(IntersectionTime>PIC::Mesh::mesh.EPS)) {
            //an intersection is found --> there is no direct access to the xTarget
            TriangleFace->RayIntersection(x,l,IntersectionTime,PIC::Mesh::mesh.EPS);

            return false;
          }
          else TriangleFace->pic__RayTracingOperationCounterTable[thread]=FaceAccessCounterTable[thread];
        }
      }
		}

		//determine the new node
		double xNodeExit[3],xFaceExitLocal[2];
		int nExitFace,iFace,jFace;

		if (PIC::RayTracing::GetBlockExitPoint(node->xmin,node->xmax,x,l,xNodeExit,xFaceExitLocal,nExitFace)==true) {
      iFace=(xFaceExitLocal[0]<0.5) ? 0 : 1;
		  jFace=(xFaceExitLocal[1]<0.5) ? 0 : 1;

		  node=node->GetNeibFace(nExitFace,iFace,jFace);
		}
		else {
		  PIC::RayTracing::GetBlockExitPoint(node->xmin,node->xmax,x,l,xNodeExit,xFaceExitLocal,nExitFace);
		  node=PIC::Mesh::mesh.findTreeNode(xNodeExit,node);
		}

		memcpy(x,xNodeExit,3*sizeof(double));
	}

	return true;
}

int PIC::RayTracing::CountFaceIntersectionNumber(double *xStart,double *xTarget,int MeshFileID,void* ExeptionFace) {
  double c=0.0,l[3],x[3];
  int idim,IntersectionCounter=0;


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

  node=PIC::Mesh::mesh.findTreeNode(xStart);
  if (node==NULL) return false;

  //determine the thread number
  int thread=0;

  #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  thread=omp_get_thread_num();
  #endif

  //increment the call counter
  FaceAccessCounterTable[thread]++;

  if (FaceAccessCounterTable[thread]==(char)0xff) {
    //if the counter is 'owerflow' than re-init the counter value and the counters of the cut-faces
    FaceAccessCounterTable[thread]=(char)1;
    FlushFaceRayTracingOperationIDTable(thread);
  }

  //trace the ray
  PIC::Mesh::IrregularSurface::cTriangleFaceDescriptor *t;

  while (node!=NULL) {
    PIC::Mesh::IrregularSurface::cTriangleFace *TriangleFace;
    double IntersectionTime;
    bool code;

    if ((t=node->FirstTriangleCutFace)!=NULL) {
      for (;t!=NULL;t=t->next) {
        if ((TriangleFace=t->TriangleFace)->pic__RayTracingOperationCounterTable[thread]!=FaceAccessCounterTable[thread]) {
          if ((MeshFileID<0)||(TriangleFace->MeshFileID==MeshFileID)) {
            code=TriangleFace->RayIntersection(x,l,IntersectionTime,PIC::Mesh::mesh.EPS);
            if ((TriangleFace!=(PIC::Mesh::IrregularSurface::cTriangleFace*)ExeptionFace) && (code==true)/*&&(IntersectionTime>PIC::Mesh::mesh.EPS)*/) IntersectionCounter++;
          }

          TriangleFace->pic__RayTracingOperationCounterTable[thread]=FaceAccessCounterTable[thread];
        }
      }
    }

    //determine the new node
    double xNodeExit[3],xFaceExitLocal[2];
    int nExitFace,iFace,jFace;

    if (PIC::RayTracing::GetBlockExitPoint(node->xmin,node->xmax,x,l,xNodeExit,xFaceExitLocal,nExitFace)==true) {
      iFace=(xFaceExitLocal[0]<0.5) ? 0 : 1;
      jFace=(xFaceExitLocal[1]<0.5) ? 0 : 1;

      node=node->GetNeibFace(nExitFace,iFace,jFace);
    }
    else {
      PIC::RayTracing::GetBlockExitPoint(node->xmin,node->xmax,x,l,xNodeExit,xFaceExitLocal,nExitFace);
      node=PIC::Mesh::mesh.findTreeNode(xNodeExit,node);
    }

    memcpy(x,xNodeExit,3*sizeof(double));
  }

  return IntersectionCounter;
}

int PIC::RayTracing::FindFistIntersectedFace(double *x0Ray,double *lRay,double *xIntersection,void* ExeptionFace) {
  double x[3],c=0.0;
  int idim;
  int IntersectedFaceNumber=-1;

  //determine the initial block
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=NULL;

  node=PIC::Mesh::mesh.findTreeNode(x0Ray);
  if (node==NULL) return false;

  //determine the thread number
  int thread=0;

  #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  thread=omp_get_thread_num();
  #endif

  //increment the call counter
  FaceAccessCounterTable[thread]++;

  if (FaceAccessCounterTable[thread]==(char)0xff) {
    //if the counter is 'owerflow' than re-init the counter value and the counters of the cut-faces
    FaceAccessCounterTable[thread]=(char)1;
    FlushFaceRayTracingOperationIDTable(thread);
  }

  memcpy(x,x0Ray,3*sizeof(double));

  //trace the ray
  PIC::Mesh::IrregularSurface::cTriangleFaceDescriptor *t;

  while (node!=NULL) {
    PIC::Mesh::IrregularSurface::cTriangleFace *TriangleFace;
    PIC::Mesh::IrregularSurface::cTriangleFace *minIntersectionTimeTriangleFace=NULL;
    double IntersectionTime,minIntersectionTime=-1.0;
    double xIntersectionTemp[3],xMinTimeIntersection[3];

    bool code;

    if ((t=node->FirstTriangleCutFace)!=NULL) {
      for (;t!=NULL;t=t->next) {
        if ((TriangleFace=t->TriangleFace)->pic__RayTracingOperationCounterTable[thread]!=FaceAccessCounterTable[thread])  {
          code=TriangleFace->RayIntersection(x0Ray,lRay,IntersectionTime,xIntersectionTemp,PIC::Mesh::mesh.EPS);

          if ((TriangleFace!=(PIC::Mesh::IrregularSurface::cTriangleFace*)ExeptionFace) && (code==true)) {
            if ((minIntersectionTimeTriangleFace==NULL) || (IntersectionTime>minIntersectionTime)) {
              minIntersectionTimeTriangleFace=TriangleFace;
              minIntersectionTime=IntersectionTime;
              memcpy(xMinTimeIntersection,xIntersectionTemp,3*sizeof(double));
            }
          }

          TriangleFace->pic__RayTracingOperationCounterTable[thread]=FaceAccessCounterTable[thread];
        }
      }
    }

    if (minIntersectionTimeTriangleFace!=NULL) {
      //intersection face is found -> return it
      memcpy(xIntersection,xMinTimeIntersection,3*sizeof(double));
      return minIntersectionTimeTriangleFace-PIC::Mesh::IrregularSurface::BoundaryTriangleFaces;
    }

    //determine the new node
    double xNodeExit[3],xFaceExitLocal[2];
    int nExitFace,iFace,jFace;

    if (PIC::RayTracing::GetBlockExitPoint(node->xmin,node->xmax,x,lRay,xNodeExit,xFaceExitLocal,nExitFace)==true) {
      iFace=(xFaceExitLocal[0]<0.5) ? 0 : 1;
      jFace=(xFaceExitLocal[1]<0.5) ? 0 : 1;

      node=node->GetNeibFace(nExitFace,iFace,jFace);
    }
    else {
      PIC::RayTracing::GetBlockExitPoint(node->xmin,node->xmax,x,lRay,xNodeExit,xFaceExitLocal,nExitFace);
      node=PIC::Mesh::mesh.findTreeNode(xNodeExit,node);
    }

    memcpy(x,xNodeExit,3*sizeof(double));
  }

  //no face of intersection was found -> return the default value
  return -1;
}


void PIC::RayTracing::SetCutCellShadowAttribute(double *xLightSource,bool ParallelExecution) {
  double xTriangle[3],c;
  int i,idim,iStart,iFinish,nFaceThread;

  //init the ray traciing module if needed
  PIC::RayTracing::Init();

  if (ParallelExecution==true) {
    nFaceThread=PIC::Mesh::IrregularSurface::nBoundaryTriangleFaces/PIC::nTotalThreads;
    iStart=nFaceThread*PIC::ThisThread;
    iFinish=iStart+nFaceThread;
    if (PIC::ThisThread==PIC::nTotalThreads-1) iFinish=PIC::Mesh::IrregularSurface::nBoundaryTriangleFaces;
  }
  else iStart=0,iFinish=PIC::Mesh::IrregularSurface::nBoundaryTriangleFaces;

  for (i=iStart;i<iFinish;i++) {
    PIC::Mesh::IrregularSurface::BoundaryTriangleFaces[i].GetCenterPosition(xTriangle);
    for (c=0.0,idim=0;idim<3;idim++) c+=(xLightSource[idim]-xTriangle[idim])*PIC::Mesh::IrregularSurface::BoundaryTriangleFaces[i].ExternalNormal[idim];

    if (c<0.0) {
      //the point is in the shadow
      PIC::Mesh::IrregularSurface::BoundaryTriangleFaces[i].pic__shadow_attribute=_PIC__CUT_FACE_SHADOW_ATTRIBUTE__TRUE_;
    }
    else {
      //trace the ray to the light source
      if (TestDirectAccess(xTriangle,xLightSource)==false) {
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
}

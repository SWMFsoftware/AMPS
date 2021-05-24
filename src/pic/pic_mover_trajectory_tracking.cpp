
#include "pic.h"


int PIC::Mover::TrajectoryTrackingMover(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode,CutCell::cTriangleFace* ExcludeCutTriangleFace) {
  namespace PB = PIC::ParticleBuffer;

  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *newNode=NULL;
  PIC::ParticleBuffer::byte *ParticleData;
  double vInit[3],xInit[3]={0.0,0.0,0.0};
  int idim,i,j,k,spec;

  double *vFinal=vInit,*xFinal=xInit;

  ParticleData=PB::GetParticleDataPointer(ptr);
  PB::GetV(vInit,ParticleData);
  PB::GetX(xInit,ParticleData);
  spec=PB::GetI(ParticleData);

  //the description of the boundaries of the block faces
  struct cExternalBoundaryFace {
    double norm[3];
    int nX0[3];
    double e0[3],e1[3],x0[3];
    double lE0,lE1;
  };

  static bool initExternalBoundaryFaceTable=false;

  static cExternalBoundaryFace ExternalBoundaryFaceTable[6]={
      {{-1.0,0.0,0.0}, {0,0,0}, {0,1,0},{0,0,1},{0,0,0}, 0.0,0.0}, {{1.0,0.0,0.0}, {1,0,0}, {0,1,0},{0,0,1},{0,0,0}, 0.0,0.0},
      {{0.0,-1.0,0.0}, {0,0,0}, {1,0,0},{0,0,1},{0,0,0}, 0.0,0.0}, {{0.0,1.0,0.0}, {0,1,0}, {1,0,0},{0,0,1},{0,0,0}, 0.0,0.0},
      {{0.0,0.0,-1.0}, {0,0,0}, {1,0,0},{0,1,0},{0,0,0}, 0.0,0.0}, {{0.0,0.0,1.0}, {0,0,1}, {1,0,0},{0,1,0},{0,0,0}, 0.0,0.0}
  };

  if (initExternalBoundaryFaceTable==false) {
    initExternalBoundaryFaceTable=true;

    for (int nface=0;nface<6;nface++) {
      double cE0=0.0,cE1=0.0;

      for (idim=0;idim<3;idim++) {
        ExternalBoundaryFaceTable[nface].x0[idim]=(ExternalBoundaryFaceTable[nface].nX0[idim]==0) ? PIC::Mesh::mesh->rootTree->xmin[idim] : PIC::Mesh::mesh->rootTree->xmax[idim];

        cE0+=pow(((ExternalBoundaryFaceTable[nface].e0[idim]+ExternalBoundaryFaceTable[nface].nX0[idim]<0.5) ? PIC::Mesh::mesh->rootTree->xmin[idim] : PIC::Mesh::mesh->rootTree->xmax[idim])-ExternalBoundaryFaceTable[nface].x0[idim],2);
        cE1+=pow(((ExternalBoundaryFaceTable[nface].e1[idim]+ExternalBoundaryFaceTable[nface].nX0[idim]<0.5) ? PIC::Mesh::mesh->rootTree->xmin[idim] : PIC::Mesh::mesh->rootTree->xmax[idim])-ExternalBoundaryFaceTable[nface].x0[idim],2);
      }

      ExternalBoundaryFaceTable[nface].lE0=sqrt(cE0);
      ExternalBoundaryFaceTable[nface].lE1=sqrt(cE1);
    }
  }

  static long int nLoop=0;
  static long int nCall=0;
  nCall++;



  //determine the flight time to the mearest boundary of the block
  auto GetBlockBoundaryFlightTime = [] (int& iIntersectedFace,int iFaceExclude,double *x,double *v,double *xmin,double *xmax) {
    int iface;
    double dt;

    iIntersectedFace=-1;

    for (iface=0;iface<6;iface++) if (iface!=iFaceExclude) {
      double dt=-1.0;
      int iOrthogonal1,iOrthogonal2;


      switch (iface) {
      case 0:
        iOrthogonal1=1,iOrthogonal2=2;
        if (v[0]!=0.0) dt=(xmin[0]-x[0])/v[0];
        break;

      case 1:
        iOrthogonal1=1,iOrthogonal2=2;
        if (v[0]!=0.0) dt=(xmax[0]-x[0])/v[0];
        break;

      case 2:
        iOrthogonal1=0,iOrthogonal2=2;
        if (v[1]!=0.0) dt=(xmin[1]-x[1])/v[1];
        break;

      case 3:
        iOrthogonal1=0,iOrthogonal2=2;
        if (v[1]!=0.0) dt=(xmax[1]-x[1])/v[1];

        break;
      case 4:
        iOrthogonal1=0,iOrthogonal2=1;
        if (v[2]!=0.0) dt=(xmin[2]-x[2])/v[2];

        break;
      case 5:
        iOrthogonal1=0,iOrthogonal2=1;
        if (v[2]!=0.0) dt=(xmax[2]-x[2])/v[2];

        break;
      }


      if (dt>0.0) {
        double t;

        t=x[iOrthogonal1]+v[iOrthogonal1]*dt;

        if ((xmin[iOrthogonal1]<=t)&&(t<=xmax[iOrthogonal1])) {
          t=x[iOrthogonal2]+dt*v[iOrthogonal2];

          if ((xmin[iOrthogonal2]<=t)&&(t<=xmax[iOrthogonal2])) {
              iIntersectedFace=iface;
              return dt;
          }
        }
      }
    }

    //check the distance of the point ot the boundary of the block: is the distance is less that EPS -> return dt=0
    for (int i=0;i<3;i++) {
      if (fabs(x[i]-xmin[i])<PIC::Mesh::mesh->EPS) {
        iIntersectedFace=2*i;
        return 0.0;
      }

      if (fabs(x[i]-xmax[i])<PIC::Mesh::mesh->EPS) {
        iIntersectedFace=1+2*i;
        return 0.0;
      }
    }

    return -1.0;
  };

  //determine the flight time to the neares cut surface
  auto GetCutSurfaceIntersectionTime = [] (CutCell::cTriangleFace* &CutTriangleFace,CutCell::cTriangleFace* ExcludeCutTriangleFace,double *x,double *v,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
    CutCell::cTriangleFaceDescriptor *t;
    CutCell::cTriangleFace *TriangleFace;
    double dt,FlightTime;

    CutTriangleFace=NULL,FlightTime=-1.0;

    for (t=node->FirstTriangleCutFace;t!=NULL;t=t->next) {
      TriangleFace=t->TriangleFace;

      if (TriangleFace!=ExcludeCutTriangleFace) {
        double xLocalIntersection[3],xIntersection[3];

        if (TriangleFace->RayIntersection(x,v,dt,xLocalIntersection,xIntersection,PIC::Mesh::mesh->EPS)==true) {
          if ((dt>0.0)&&(Vector3D::DotProduct(v,TriangleFace->ExternalNormal)<0.0)&&(dt*Vector3D::Length(v)>PIC::Mesh::mesh->EPS)) if ((CutTriangleFace==NULL)||(dt<FlightTime)) {
            CutTriangleFace=TriangleFace,FlightTime=dt;
          }
        }
      }
    }

    return FlightTime;
  };

  static const int iFaceExcludeTable[6]={1,0,3,2,5,4};
  int IntersectionMode,iIntersectedFace,iFaceExclude=-1;
  double FaceIntersectionFlightTime,dt;

  double CutTriangleIntersevtionFlightTime;
  CutCell::cTriangleFace* CutTriangleFace;

  const int _block_bounday=0;
  const int _cut_triangle=1;
  const int _undefined_boundary=2;


  while (dtTotal>0.0) {
    nLoop++;

    //determine the flight time to the boundary of the block
    //
    

int code,iIntersectedFace_debug=iIntersectedFace;
int iFaceExclude_debug=iFaceExclude;
double xInit_debug[3]={xInit[0],xInit[1],xInit[2]};
double vInit_debug[3]={vInit[0],vInit[1],vInit[2]};


    IntersectionMode=_undefined_boundary;


    FaceIntersectionFlightTime=GetBlockBoundaryFlightTime(iIntersectedFace,iFaceExclude,xInit,vInit,startNode->xmin,startNode->xmax);


    //determine the flight time to the nearest cut triangle
    if (startNode->FirstTriangleCutFace!=NULL) {
      CutTriangleIntersevtionFlightTime=GetCutSurfaceIntersectionTime(CutTriangleFace,ExcludeCutTriangleFace,xInit,vInit,startNode);


      bool cut_face_intersection=false;

      if (CutTriangleFace!=NULL) {
        if (CutTriangleIntersevtionFlightTime<=FaceIntersectionFlightTime) {
          cut_face_intersection=true;
        }
        else {
          double tt=CutTriangleIntersevtionFlightTime-FaceIntersectionFlightTime;

          if (Vector3D::DotProduct(vInit,vInit)*tt*tt<PIC::Mesh::mesh->EPS*PIC::Mesh::mesh->EPS) {
            cut_face_intersection=true;
          }
        }
      }

      if (cut_face_intersection==true) { //((CutTriangleIntersevtionFlightTime<FaceIntersectionFlightTime)&&(CutTriangleFace!=NULL)) {
        ExcludeCutTriangleFace=CutTriangleFace;
        dt=CutTriangleIntersevtionFlightTime;
        iFaceExclude=-1;
        IntersectionMode=_cut_triangle;
      }
      else if (iIntersectedFace!=-1) {
        ExcludeCutTriangleFace=NULL;
        dt=FaceIntersectionFlightTime;
        IntersectionMode=_block_bounday;
        iFaceExclude=iFaceExcludeTable[iIntersectedFace];
      }
    }
    else if (iIntersectedFace!=-1) {
      ExcludeCutTriangleFace=NULL;
      dt=FaceIntersectionFlightTime;
      IntersectionMode=_block_bounday;
      iFaceExclude=iFaceExcludeTable[iIntersectedFace];
    }


    if (IntersectionMode==_undefined_boundary) {
      exit(__LINE__,__FILE__,"Error: the boundary type is not defined");
    }


    //advance the particle location and velocity
    if (dtTotal<dt) {
      for (idim=0;idim<3;idim++) {
        xInit[idim]+=dtTotal*vInit[idim];
        vInit[idim]+=0.0;
      }

      dtTotal=0.0;
    }
    else {
      //the particle has intersected either with a cut-surface or the boundary of the block
      for (idim=0;idim<3;idim++) {
        xInit[idim]+=dt*vInit[idim];
        vInit[idim]+=0.0;
      }

      dtTotal-=dt;

      //determine the next block that particle is in
      double x_test[3],c_init;

      switch (IntersectionMode) {
      case _block_bounday:
        for (int i=0;i<3;i++) x_test[i]=xInit[i];

        switch (iIntersectedFace) {
        case 0:
          x_test[0]-=0.01*(startNode->xmax[0]-startNode->xmin[0]);
          break;
        case 1:
          x_test[0]+=0.01*(startNode->xmax[0]-startNode->xmin[0]);
          break;


        case 2:
          x_test[1]-=0.01*(startNode->xmax[1]-startNode->xmin[1]);
          break;
        case 3:
          x_test[1]+=0.01*(startNode->xmax[1]-startNode->xmin[1]);
          break;


        case 4:
          x_test[2]-=0.01*(startNode->xmax[2]-startNode->xmin[2]);
          break;
        case 5:
          x_test[2]+=0.01*(startNode->xmax[2]-startNode->xmin[2]);
          break;
        }

        startNode=PIC::Mesh::mesh->findTreeNode(x_test,startNode);
        if (startNode==NULL) {
          for (int i=0;i<3;i++) xInit[idim]=x_test[idim];
        }

        break;
      case _cut_triangle:
        c_init=Vector3D::DotProduct(vInit,CutTriangleFace->ExternalNormal);

        do {
          code=(ProcessTriangleCutFaceIntersection!=NULL) ? ProcessTriangleCutFaceIntersection(ptr,xInit,vInit,CutTriangleFace,startNode) : _PARTICLE_DELETED_ON_THE_FACE_;

          if (code==_PARTICLE_DELETED_ON_THE_FACE_) {
            PIC::ParticleBuffer::DeleteParticle(ptr);
            return _PARTICLE_LEFT_THE_DOMAIN_;
          }
        }
        while (c_init*Vector3D::DotProduct(vInit,CutTriangleFace->ExternalNormal)>=0.0); // (Vector3D::DotProduct(v,v_init)>=0.0);

        startNode=PIC::Mesh::mesh->findTreeNode(xInit,startNode);
        break;
      default:
        exit(__LINE__,__FILE__,"Error: the option is unknown");
      }

      if (startNode==NULL) {
        break;
      }
    }
  }


  newNode=PIC::Mesh::mesh->findTreeNode(xInit,startNode);

  //move the particle inside the block
  if (newNode!=NULL) for (int i=0;i<3;i++) {
    if (newNode->xmin[i]+PIC::Mesh::mesh->EPS>xInit[i]) xInit[i]=newNode->xmin[i]+PIC::Mesh::mesh->EPS;
    if (newNode->xmax[i]-PIC::Mesh::mesh->EPS<xInit[i]) xInit[i]=newNode->xmax[i]-PIC::Mesh::mesh->EPS; 
  }

  if (newNode==NULL) {
    //the particle left the computational domain
    int code=_PARTICLE_DELETED_ON_THE_FACE_;

#if _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE_ == _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__DELETE_
    //do nothing -> the particle deleting code is already set
#else
    //call the function that process particles that leaved the coputational domain
    //       if (ProcessOutsideDomainParticles!=NULL) {
    //determine through which face the particle left the domain

    int nface,nIntersectionFace;
    double tVelocityIncrement,cx,cv,r0[3],dt,vMiddle[3]={0.5*(vInit[0]+vFinal[0]),0.5*(vInit[1]+vFinal[1]),0.5*(vInit[2]+vFinal[2])},c,dtIntersection=-1.0;

    for (nface=0;nface<6;nface++) {
      for (idim=0,cx=0.0,cv=0.0;idim<3;idim++) {
        r0[idim]=xInit[idim]-ExternalBoundaryFaceTable[nface].x0[idim];
        cx+=r0[idim]*ExternalBoundaryFaceTable[nface].norm[idim];
        cv+=vMiddle[idim]*ExternalBoundaryFaceTable[nface].norm[idim];
      }

      if (cv>0.0) {
        dt=-cx/cv;

        if ((dtIntersection<0.0)||(dt<dtIntersection)&&(dt>0.0)) {
          double cE0=0.0,cE1=0.0;

          for (idim=0;idim<3;idim++) {
            c=r0[idim]+dt*vMiddle[idim];

            cE0+=c*ExternalBoundaryFaceTable[nface].e0[idim],cE1+=c*ExternalBoundaryFaceTable[nface].e1[idim];
          }

          if ((cE0<-PIC::Mesh::mesh->EPS)||(cE0>ExternalBoundaryFaceTable[nface].lE0+PIC::Mesh::mesh->EPS) || (cE1<-PIC::Mesh::mesh->EPS)||(cE1>ExternalBoundaryFaceTable[nface].lE1+PIC::Mesh::mesh->EPS)) continue;

          nIntersectionFace=nface,dtIntersection=dt;
        }
      }
    }

    if (nIntersectionFace==-1) exit(__LINE__,__FILE__,"Error: cannot find the face of the intersection");

    for (idim=0,tVelocityIncrement=((dtIntersection/dtTotal<1) ? dtIntersection/dtTotal : 1);idim<3;idim++) {
      xInit[idim]+=dtIntersection*vMiddle[idim]-ExternalBoundaryFaceTable[nIntersectionFace].norm[idim]*PIC::Mesh::mesh->EPS;
      vInit[idim]+=tVelocityIncrement*(vFinal[idim]-vInit[idim]);
    }

    newNode=PIC::Mesh::mesh->findTreeNode(xInit,startNode);

    if (newNode==NULL) {
      //the partcle is outside of the domain -> correct particle location and determine the newNode;
      double xmin[3],xmax[3];
      int ii;

      memcpy(xmin,PIC::Mesh::mesh->xGlobalMin,3*sizeof(double));
      memcpy(xmax,PIC::Mesh::mesh->xGlobalMax,3*sizeof(double));

      for (ii=0;ii<3;ii++) {
        if (xmin[ii]>=xInit[ii]) xInit[ii]=xmin[ii]+PIC::Mesh::mesh->EPS;
        if (xmax[ii]<=xInit[ii]) xInit[ii]=xmax[ii]-PIC::Mesh::mesh->EPS;
      }

      newNode=PIC::Mesh::mesh->findTreeNode(xInit,startNode);

      if (newNode==NULL) exit(__LINE__,__FILE__,"Error: cannot find the node");
    }

#if _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE_ == _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__USER_FUNCTION_
    code=ProcessOutsideDomainParticles(ptr,xInit,vInit,nIntersectionFace,newNode);
#elif _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE_ == _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__PERIODIC_CONDITION_
    exit(_LINE__,__FILE__,"Error: not implemented");
#elif _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE_ == _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__SPECULAR_REFLECTION_
    //reflect the particle back into the domain
    {
      double c=0.0;
      for (int idim=0;idim<3;idim++) c+=ExternalBoundaryFaceTable[nIntersectionFace].norm[idim]*vInit[idim];
      for (int idim=0;idim<3;idim++) vInit[idim]-=2.0*c*ExternalBoundaryFaceTable[nIntersectionFace].norm[idim];
    }

    code=_PARTICLE_REJECTED_ON_THE_FACE_;
#else
    exit(__LINE__,__FILE__,"Error: the option is unknown");
#endif




    //       }
#endif //_PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE_ == _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__DELETE

    //call the function that process particles that leaved the coputational domain
    switch (code) {
    case _PARTICLE_DELETED_ON_THE_FACE_:
      PIC::ParticleBuffer::DeleteParticle(ptr);
      return _PARTICLE_LEFT_THE_DOMAIN_;
    default:
      exit(__LINE__,__FILE__,"Error: not implemented");
    }
  }

  //save the trajectory point
#if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
  PIC::ParticleTracker::RecordTrajectoryPoint(xFinal,vFinal,spec,ParticleData,(void*)newNode);

#if _PIC_PARTICLE_TRACKER__TRACKING_CONDITION_MODE__DYNAMICS_ == _PIC_MODE_ON_
  PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(xFinal,vFinal,spec,ParticleData,(void*)newNode);
#endif
#endif

  //finish the trajectory integration procedure
  PIC::Mesh::cDataBlockAMR *block;

  if (PIC::Mesh::mesh->fingCellIndex(xFinal,i,j,k,newNode,false)==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located");

  if ((block=newNode->block)==NULL) {
    exit(__LINE__,__FILE__,"Error: the block is empty. Most probably hte tiime step is too long");
  }

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



  PIC::ParticleBuffer::SetV(vFinal,ParticleData);
  PIC::ParticleBuffer::SetX(xFinal,ParticleData);

  return _PARTICLE_MOTION_FINISHED_;
}






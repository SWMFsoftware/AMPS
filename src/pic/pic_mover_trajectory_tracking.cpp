
#include "pic.h"


int PIC::Mover::TrajectoryTrackingMover(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode,CutCell::cTriangleFace* ExcludeCutTriangleFace) {
  namespace PB = PIC::ParticleBuffer;

  PIC::ParticleBuffer::byte *ParticleData;
  double x[3],v[3];
  int spec;

  ParticleData=PB::GetParticleDataPointer(ptr);
  PB::GetV(v,ParticleData);
  PB::GetX(x,ParticleData);
  spec=PIC::ParticleBuffer::GetI(ParticleData);

  //determine the flight time to the mearest boundary of the block 
  auto GetBlockBoundaryFlightTime = [] (int& iIntersectedFace,double& FlightTime,int iFaceExclude,double *x,double *v,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* Node) {
    int iface;
    double *xmin=Node->xmin;
    double *xmax=Node->xmax;
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
      

      if (dt>=0.0) {
        double t;

        t=x[iOrthogonal1]+v[iOrthogonal1]*dt;

        if ((xmin[iOrthogonal1]<=t)&&(t<=xmax[iOrthogonal1])) {
          t=x[iOrthogonal2]+dt*v[iOrthogonal2];

          if ((xmin[iOrthogonal2]<=t)&&(t<=xmax[iOrthogonal2])) {
              iIntersectedFace=iface,FlightTime=dt;
              return;
          }
        }
      }
    }
  };
  

  //determine the flight time to the neares cut surface 
  auto GetCutSurfaceIntersectionTime = [] (CutCell::cTriangleFace* &CutTriangleFace,double& FlightTime,CutCell::cTriangleFace* ExcludeCutTriangleFace,double *x,double *v,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
    CutCell::cTriangleFaceDescriptor *t;
    CutCell::cTriangleFace *TriangleFace;
    double dt;

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
  };
                            

  int iFaceExclude=-1,iIntersectedFace;
  double FaceIntersectionFlightTime;
  double dt,accl[3];


  double CutTriangleIntersevtionFlightTime;
  CutCell::cTriangleFace* CutTriangleFace;

  const int _block_bounday=0;
  const int _cut_triangle=1; 
  const int _undefined_boundary=2; 

  int IntersectionMode=_undefined_boundary;

  while (dtTotal>0) {
	  //determine the flight time to the boundary of the block
	  GetBlockBoundaryFlightTime(iIntersectedFace,FaceIntersectionFlightTime,iFaceExclude,x,v,startNode);

	  //determine the flight time to the nearest cut triangle
	  if (startNode->FirstTriangleCutFace!=NULL) {
		  GetCutSurfaceIntersectionTime(CutTriangleFace,CutTriangleIntersevtionFlightTime,ExcludeCutTriangleFace,x,v,startNode);

		  if ((CutTriangleIntersevtionFlightTime<FaceIntersectionFlightTime)&&(CutTriangleFace!=NULL)) {
			  ExcludeCutTriangleFace=CutTriangleFace;
			  dt=CutTriangleIntersevtionFlightTime;
			  iFaceExclude=-1;
			  IntersectionMode=_cut_triangle;
		  }
		  else if (iIntersectedFace!=-1) {
			  ExcludeCutTriangleFace=NULL;
			  dt=FaceIntersectionFlightTime;
			  IntersectionMode=_block_bounday;
		  }
	  }
	  else if (iIntersectedFace!=-1) {
		  ExcludeCutTriangleFace=NULL;
		  dt=FaceIntersectionFlightTime;
		  IntersectionMode=_block_bounday;
	  }


	  if (IntersectionMode==_undefined_boundary) {
	    exit(__LINE__,__FILE__,"Error: the boundary type is not defined");
	  }

	  //advance the location of the particle
	  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node_tmp=PIC::Mesh::mesh->findTreeNode(x,startNode); 

    _PIC_PARTICLE_MOVER__TOTAL_PARTICLE_ACCELERATION_(accl,spec,ptr,x,v,node_tmp);


	  if (dtTotal<dt) {
		  //the particles did not reach the boundary
		  for (int idim=0;idim<3;idim++) {
		    x[idim]+=dtTotal*v[idim];
		    v[idim]+=dtTotal*accl[idim];
		  }

		  dtTotal=0.0;
	  }
	  else {
		  //the particle intersected the bounday
		  int code;
                  double x_test[3],c_init;
                  

		  switch (IntersectionMode) {
		  case _cut_triangle:
			  for (int idim=0;idim<3;idim++) {
			    x[idim]+=dt*v[idim];
			    v[idim]+=dt*accl[idim];
			  }

c_init=Vector3D::DotProduct(v,CutTriangleFace->ExternalNormal);


			  do {
			    code=(ProcessTriangleCutFaceIntersection!=NULL) ? ProcessTriangleCutFaceIntersection(ptr,x,v,CutTriangleFace,startNode) : _PARTICLE_DELETED_ON_THE_FACE_;

			    if (code==_PARTICLE_DELETED_ON_THE_FACE_) {
				    PIC::ParticleBuffer::DeleteParticle(ptr);
				    return _PARTICLE_LEFT_THE_DOMAIN_;
			    }
			  }
			  while (c_init*Vector3D::DotProduct(v,CutTriangleFace->ExternalNormal)>=0.0); // (Vector3D::DotProduct(v,v_init)>=0.0);

			  break;
		  case _block_bounday:
			  for (int idim=0;idim<3;idim++) {
			    x[idim]+=dt*v[idim];
			    v[idim]+=dt*accl[idim];
 
                            x_test[idim]=x[idim];
			  }

			  switch (iIntersectedFace) {
			  case 0:case 1:
				  if (iIntersectedFace==0) {
					  iFaceExclude=1;
					  if (x_test[0]>=startNode->xmin[0]) x_test[0]=startNode->xmin[0]-PIC::Mesh::mesh->EPS;
				  }
				  else {
					  iFaceExclude=0;
					  if (x_test[0]<startNode->xmax[0]) x_test[0]=startNode->xmax[0];
				  }

				  for (int iOrthogonal=1;iOrthogonal<3;iOrthogonal++) {
					  if (x_test[iOrthogonal]<startNode->xmin[iOrthogonal])  x_test[iOrthogonal]=startNode->xmin[iOrthogonal];
					  if (x_test[iOrthogonal]>=startNode->xmax[iOrthogonal]) x_test[iOrthogonal]=startNode->xmax[iOrthogonal]-PIC::Mesh::mesh->EPS;
				  }
				  break;

			  case 2:case 3:
				  if (iIntersectedFace==2) {
					  iFaceExclude=3;
					  if (x_test[1]>=startNode->xmin[1]) x_test[1]=startNode->xmin[1]-PIC::Mesh::mesh->EPS;
				  }
				  else {
					  iFaceExclude=2;
					  if (x_test[1]<startNode->xmax[1]) x_test[1]=startNode->xmax[1];
				  }

				  for (int iOrthogonal=0;iOrthogonal<3;iOrthogonal+=2) {
					  if (x_test[iOrthogonal]<startNode->xmin[iOrthogonal])  x_test[iOrthogonal]=startNode->xmin[iOrthogonal];
					  if (x_test[iOrthogonal]>=startNode->xmax[iOrthogonal]) x_test[iOrthogonal]=startNode->xmax[iOrthogonal]-PIC::Mesh::mesh->EPS;
				  }
				  break;

			  case 4:case 5:
				  if (iIntersectedFace==4) {
					  iFaceExclude=5;
					  if (x_test[2]>=startNode->xmin[2]) x_test[2]=startNode->xmin[2]-PIC::Mesh::mesh->EPS;
				  }
				  else {
					  iFaceExclude=4;
					  if (x_test[2]<startNode->xmax[2]) x_test[2]=startNode->xmax[2];
				  }

				  for (int iOrthogonal=0;iOrthogonal<2;iOrthogonal++) {
					  if (x_test[iOrthogonal]<startNode->xmin[iOrthogonal])  x_test[iOrthogonal]=startNode->xmin[iOrthogonal];
					  if (x_test[iOrthogonal]>=startNode->xmax[iOrthogonal]) x_test[iOrthogonal]=startNode->xmax[iOrthogonal]-PIC::Mesh::mesh->EPS;
				  }
				  break;
			  }

			  startNode=PIC::Mesh::mesh->findTreeNode(x_test,startNode);

			  if (startNode==NULL) {
				  //the partcle is outside of the domain
				  switch (_PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE_) {
				  case _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__USER_FUNCTION_:
					  code=ProcessOutsideDomainParticles(ptr,x,v,iFaceExclude,startNode);
					  break;

				  case _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__PERIODIC_CONDITION_:
					  exit(__LINE__,__FILE__,"Error: not implemented");
					  break;

				  case _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__SPECULAR_REFLECTION_:
					  exit(__LINE__,__FILE__,"Error: not implemented");
					  break;

				  case _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__DELETE_:
					  code=_PARTICLE_DELETED_ON_THE_FACE_;
					  break;
				  }


				  switch (code) {
				  case _PARTICLE_DELETED_ON_THE_FACE_:
					  PIC::ParticleBuffer::DeleteParticle(ptr);
					  return _PARTICLE_LEFT_THE_DOMAIN_;

				  case _PARTICLE_REJECTED_ON_THE_FACE_:
					  //apply a surface/particle interaction model
					  exit(__LINE__,__FILE__,"Error: not implemented");



					  break;
				  default:
					  exit(__LINE__,__FILE__,"Error: not implemented");
				  }
			  }
			  break;

			  default:
			    exit(__LINE__,__FILE__,"Error: the boundary time is undefied");
		  }

	    dtTotal-=dt;
	  }

	}

  if ((_PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_)&&(_PIC_DEBUGGER_MODE__VARIABLE_VALUE_RANGE_CHECK_ == _PIC_DEBUGGER_MODE__VARIABLE_VALUE_RANGE_CHECK_ON_)) {
    PIC::Debugger::CatchOutLimitValue(v,DIM,__LINE__,__FILE__);
    PIC::Debugger::CatchOutLimitValue(x,DIM,__LINE__,__FILE__);
  }

  //save the particle velocity and location, and update the particle list
  PB::SetV(v,ParticleData);
  PB::SetX(x,ParticleData);

  //get the cell where the particle is located
  int LocalCellNumber,i,j,k;
  
  if ((LocalCellNumber=PIC::Mesh::mesh->fingCellIndex(x,i,j,k,startNode,false))==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located4");

  #if _PIC_MOVER__MPI_MULTITHREAD_ == _PIC_MODE_ON_
  PIC::ParticleBuffer::SetPrev(-1,ParticleData);

  long int tempFirstCellParticle=atomic_exchange(startNode->block->tempParticleMovingListTable+i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k),ptr);
  PIC::ParticleBuffer::SetNext(tempFirstCellParticle,ParticleData);

  if (tempFirstCellParticle!=-1) PIC::ParticleBuffer::SetPrev(ptr,tempFirstCellParticle);

  #elif _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
  long int tempFirstCellParticle,*tempFirstCellParticlePtr;

  tempFirstCellParticlePtr=startNode->block->tempParticleMovingListTable+i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k);
  tempFirstCellParticle=(*tempFirstCellParticlePtr);

  PIC::ParticleBuffer::SetNext(tempFirstCellParticle,ParticleData);
  PIC::ParticleBuffer::SetPrev(-1,ParticleData);

  if (tempFirstCellParticle!=-1) PIC::ParticleBuffer::SetPrev(ptr,tempFirstCellParticle);
  *tempFirstCellParticlePtr=ptr;

  #elif _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  PIC::Mesh::cDataBlockAMR::cTempParticleMovingListMultiThreadTable* ThreadTempParticleMovingData=startNode->block->GetTempParticleMovingListMultiThreadTable(omp_get_thread_num(),i,j,k);

  PIC::ParticleBuffer::SetNext(ThreadTempParticleMovingData->first,ParticleData);
  PIC::ParticleBuffer::SetPrev(-1,ParticleData);

  if (ThreadTempParticleMovingData->last==-1) ThreadTempParticleMovingData->last=ptr;
  if (ThreadTempParticleMovingData->first!=-1) PIC::ParticleBuffer::SetPrev(ptr,ThreadTempParticleMovingData->first);
  ThreadTempParticleMovingData->first=ptr;

  #else
  #error The option is unknown
  #endif


  //save the trajectory point
  if (_PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_) { 
    PIC::ParticleTracker::RecordTrajectoryPoint(x,v,spec,ParticleData,(void*)startNode);
  }

  return _PARTICLE_MOTION_FINISHED_;
}



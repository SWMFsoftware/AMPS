//$Id$


#include "pic.h"
#include "Exosphere.dfn"
#include "Exosphere.h"

//define whether AVX will be used in Lapenta2017 functions
#define _AVX_LAPENTA_MOVER_  _ON_

#if _AVX_INSTRUCTIONS_USAGE_MODE_ == _AVX_INSTRUCTIONS_USAGE_MODE__OFF_
#undef _AVX_LAPENTA_MOVER_
#define _AVX_LAPENTA_MOVER_  _OFF_
#else
#if _AVX_PARTICLE_MOVER_ == _OFF_
#undef _AVX_LAPENTA_MOVER_
#define _AVX_LAPENTA_MOVER_  _OFF_
#endif
#endif


void PIC::Mover::BorisSplitAcceleration_default(double *accl, double *rotation, int spec,long int ptr,double *x,double *v,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  /* function finds acceleration and splits into a simple and gyroscopic parts
   * for Lorentz force (considered here)
   * a_{total} = (q/m) * E + (q/m) V\times B = a_{simple} + \Omega\times V
   *             \_simple_/ \_ gyroscopic _/
   * acceleration due to Lorentz Force only is considered 
   **********************************************************************/
  double accl_LOCAL[3]={0.0,0.0,0.0}, rotation_LOCAL[3]={0.0,0.0,0.0};

  //#if _FORCE_LORENTZ_MODE_ == _PIC_MODE_ON_
  // find electro-magnetic field
  double E[3]={0.0,0.0,0.0},B[3]={0.0,0.0,0.0};
  
  switch(_PIC_COUPLER_MODE_) {
  case _PIC_COUPLER_MODE__OFF_:
    // if no coupler is used -> use characteristic values
    //memcpy(E,Exosphere::swE_Typical,3*sizeof(double));
    //memcpy(B,Exosphere::swB_Typical,3*sizeof(double));
    break;
  default:
    // if coupler is used -> get values from it
    //......................................................................
    // find the cell based on particles' position x and block startNode
    // input: x, startNode; output: nd, i,j,k
    {
      long int nd;  // cell's number in the block
      int i,j,k;    // cell's coordinates in the block
  
      // fail-safe check: if the block doesn't exist => exit
      if (startNode->block==NULL) exit(__LINE__,__FILE__,"Error: the block is not initialized");
  
      // flag: true - exit if a point is not found in the block / false: don't
      nd = PIC::Mesh::mesh->FindCellIndex(x,i,j,k,startNode,false);
  
      // fail-safe check: if a point isn't found, try seacrhing in other blocks
      if (nd==-1) {
        // try to found the block the point is in;
        // starting point for search is block startNode
        startNode=PIC::Mesh::mesh->findTreeNode(x,startNode);
        nd=PIC::Mesh::mesh->FindCellIndex(x,i,j,k,startNode,false);
        // if still not found => exit
  
        if (nd==-1) exit(__LINE__,__FILE__,"Error: the cell is not found");
      }
    }
  
    //......................................................................
    // finally, get fields' values at the cell
    PIC::CPLR::InitInterpolationStencil(x,startNode);
    PIC::CPLR::GetBackgroundElectricField(E);
    PIC::CPLR::GetBackgroundMagneticField(B);
    break;
  }

  //......................................................................
  // calculate acceleraton due to Lorentz force
  double ElectricCharge=0.0,mass,Charge2Mass;

#if _PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_
  if ((_DUST_SPEC_<=spec) && (spec<_DUST_SPEC_+ElectricallyChargedDust::GrainVelocityGroup::nGroups)) {
    char ParticleData[PIC::ParticleBuffer::ParticleDataLength];

    memcpy((void*)ParticleData,(void*)PIC::ParticleBuffer::GetParticleDataPointer(ptr),PIC::ParticleBuffer::ParticleDataLength);

    ElectricCharge=ElectricallyChargedDust::GetGrainCharge((PIC::ParticleBuffer::byte*)ParticleData);
    mass=ElectricallyChargedDust::GetGrainMass((PIC::ParticleBuffer::byte*)ParticleData);
  }
  else {
    ElectricCharge=PIC::MolecularData::GetElectricCharge(spec);
    mass=PIC::MolecularData::GetMass(spec);
  }
#else
  ElectricCharge=PIC::MolecularData::GetElectricCharge(spec);
  mass=PIC::MolecularData::GetMass(spec);
#endif //_PIC_MODEL__DUST__MODE_ == _PIC_MODEL__DUST__MODE__ON_

  if (ElectricCharge!=0.0) {
    Charge2Mass=ElectricCharge/mass;

    for (int idim = 0; idim<DIM; idim++){
      accl_LOCAL[idim]    += Charge2Mass*E[idim];
      rotation_LOCAL[idim]-= Charge2Mass*B[idim];
    }
  }
  //#endif//_FORCE_LORENTZ_MODE_


  //calculate the acceleration due to gravity where applied
#if _TARGET_ID_(_TARGET_) != _TARGET_NONE__ID_
  double r2=x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
  double r=sqrt(r2);
  int idim;

  for (idim=0;idim<DIM;idim++) {
    accl_LOCAL[idim]-=GravityConstant*_MASS_(_TARGET_)/r2*x[idim]/r;
  }
#endif //_TARGET_ != _TARGET_NONE_


  memcpy(accl,    accl_LOCAL,    3*sizeof(double));
  memcpy(rotation,rotation_LOCAL,3*sizeof(double));
}


int PIC::Mover::Boris(long int ptr, double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *newNode=NULL;
  PIC::ParticleBuffer::byte *ParticleData;
  int idim,i,j,k,spec;

#if _AVX_INSTRUCTIONS_USAGE_MODE_ == _AVX_INSTRUCTIONS_USAGE_MODE__OFF_
  double vInit[3],xInit[3]={0.0,0.0,0.0},vFinal[3],xFinal[3],xminBlock[3],xmaxBlock[3];
  double u[3]={0.0,0.0,0.0}, U[3]={0.0,0.0,0.0};
  double acclInit[3],rotInit[3];
#else
  double xminBlock[3],xmaxBlock[3];

  union {__m256d vInit_v; double vInit[4];};
  union {__m256d vFinal_v; double vFinal[4];};

  union {__m256d xInit_v; double xInit[4];};
  union {__m256d xFinal_v; double xFinal[4];};

  union {__m256d u_v; double u[4];};
  union {__m256d U_v; double U[4];};

  union {__m256d acclInit_v; double acclInit[4];};
  union {__m256d rotInit_v; double rotInit[4];};
#endif

  ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  PIC::ParticleBuffer::GetV(vInit,ParticleData);
  PIC::ParticleBuffer::GetX(xInit,ParticleData);
  spec=PIC::ParticleBuffer::GetI(ParticleData);

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

  static long int nCall=0;
  nCall++;

  memcpy(xminBlock,startNode->xmin,DIM*sizeof(double));
  memcpy(xmaxBlock,startNode->xmax,DIM*sizeof(double));



  if (_PIC_PARTICLE_MOVER__FORCE_INTEGRTAION_MODE_ == _PIC_PARTICLE_MOVER__FORCE_INTEGRTAION_MODE__ON_) {
    _PIC_PARTICLE_MOVER__BORIS_SPLIT_ACCELERATION_(acclInit, rotInit,spec, ptr, xInit, vInit, startNode);
  }


  // Integrate the equations of motion
  /************************ Boris method: *************************
   * dV/dt    = A        + \Omega \times V                        *
   * X_1      = X_0      + dt * V_{+1/2}                          *
   * V_{+1/2} = U        + 0.5*dt*A_0                             *
   * U        = u        + (u + (u\times h))\times s              *
   * u        = V_{-1/2} + 0.5*dt*A_0                             *
   * h        =-0.5*dt*\Omega                                     *
   * s        = 2*h/(1+|h|^2)                                     *
   ****************************************************************/
  double dtTempOverTwo,dtTemp;

  switch (_PIC_PARTICLE_MOVER__BACKWARD_TIME_INTEGRATION_MODE_) {
  case  _PIC_PARTICLE_MOVER__BACKWARD_TIME_INTEGRATION_MODE__ENABLED_:
    switch (BackwardTimeIntegrationMode) {
    case _PIC_MODE_ON_ :
      dtTemp=-dtTotal,dtTempOverTwo=-dtTotal/2.0;
      break;
    default:
      dtTemp=dtTotal,dtTempOverTwo=dtTotal/2.0;
    }
    
    break;
  default:
    dtTemp=dtTotal,dtTempOverTwo=dtTotal/2.0;
  }

#if _AVX_INSTRUCTIONS_USAGE_MODE_ == _AVX_INSTRUCTIONS_USAGE_MODE__ON_
  __m256d h_v,h2_v,uh_v,omega_v;

  u_v=_mm256_fmadd_pd(_mm256_set1_pd(dtTempOverTwo),acclInit_v,vInit_v);
  h_v=_mm256_mul_pd(_mm256_set1_pd(-dtTempOverTwo),rotInit_v);

  h2_v=_mm256_mul_pd(h_v,h_v);
  uh_v=_mm256_mul_pd(u_v,h_v);

  double h2 = h2_v[0]+h2_v[1]+h2_v[2];
  double uh = uh_v[0]+uh_v[1]+uh_v[2];


  Vector3D::CrossProduct(omega_v,u_v,h_v);

  U_v=_mm256_add_pd(
      _mm256_mul_pd(_mm256_set1_pd((1.0-h2)/(1.0+h2)),u_v),
      _mm256_mul_pd(_mm256_set1_pd(2.0/(1.0+h2)),_mm256_fmadd_pd(_mm256_set1_pd(uh),h_v,omega_v))); 

  vFinal_v=_mm256_fmadd_pd(_mm256_set1_pd(dtTempOverTwo),acclInit_v,U_v);
  xFinal_v=_mm256_fmadd_pd(_mm256_set1_pd(dtTemp),vFinal_v,xInit_v);
#else

  u[0]=vInit[0]+dtTempOverTwo*acclInit[0];
  u[1]=vInit[1]+dtTempOverTwo*acclInit[1];
  u[2]=vInit[2]+dtTempOverTwo*acclInit[2];

  double h[3];

  h[0]=-dtTempOverTwo*rotInit[0];
  h[1]=-dtTempOverTwo*rotInit[1];
  h[2]=-dtTempOverTwo*rotInit[2];

  double h2 = h[0]*h[0]+h[1]*h[1]+h[2]*h[2];
  double uh = u[0]*h[0]+u[1]*h[1]+u[2]*h[2];

  U[0] = ( (1-h2)*u[0] + 2*(u[1]*h[2]-h[1]*u[2]+uh*h[0]) ) / (1+h2);
  U[1] = ( (1-h2)*u[1] + 2*(u[2]*h[0]-h[2]*u[0]+uh*h[1]) ) / (1+h2);
  U[2] = ( (1-h2)*u[2] + 2*(u[0]*h[1]-h[0]*u[1]+uh*h[2]) ) / (1+h2);


  vFinal[0]=U[0]+dtTempOverTwo*acclInit[0];
  vFinal[1]=U[1]+dtTempOverTwo*acclInit[1];
  vFinal[2]=U[2]+dtTempOverTwo*acclInit[2];

  xFinal[0]=xInit[0]+dtTemp*vFinal[0];
  xFinal[1]=xInit[1]+dtTemp*vFinal[1];
  xFinal[2]=xInit[2]+dtTemp*vFinal[2];
#endif


  //rotate the final position
  if (_PIC_SYMMETRY_MODE_ == _PIC_SYMMETRY_MODE__AXIAL_) {
    // rotate to the y=0 plane
    double cosPhi, sinPhi, vTmpX, vTmpY;
    double xNormFinal = pow(xFinal[0]*xFinal[0]+xFinal[1]*xFinal[1], 0.5);
    cosPhi = xFinal[0] / xNormFinal;
    sinPhi = xFinal[1] / xNormFinal;
    xFinal[0] = xNormFinal;
    xFinal[1] = 0.0;
    vTmpX = vFinal[0]; vTmpY = vFinal[1];
    vFinal[0] = vTmpX*cosPhi + vTmpY*sinPhi;
    vFinal[1] =-vTmpX*sinPhi + vTmpY*cosPhi;
  }

  //interaction with the faces of the block and internal surfaces
  //check whether the particle trajectory is intersected the spherical body
#if  _TARGET_ID_(_TARGET_) != _TARGET_NONE__ID_ && _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_ 
  double rFinal2;

  //if the particle is inside the sphere -> apply the boundary condition procedure
  if ((rFinal2=xFinal[0]*xFinal[0]+xFinal[1]*xFinal[1]+xFinal[2]*xFinal[2])<_RADIUS_(_TARGET_)*_RADIUS_(_TARGET_)) {
    double r=sqrt(rFinal2);
    int code;

    static cInternalSphericalData_UserDefined::fParticleSphereInteraction ParticleSphereInteraction=
        ((cInternalSphericalData*)(PIC::Mesh::mesh->InternalBoundaryList.front().BoundaryElement))->ParticleSphereInteraction;
    static void* BoundaryElement=PIC::Mesh::mesh->InternalBoundaryList.front().BoundaryElement;

    //move the particle location at the surface of the sphere
    for (idim=0;idim<DIM;idim++) xFinal[idim]*=_RADIUS_(_TARGET_)/r;

    //determine the block of the particle location
    newNode=PIC::Mesh::mesh->findTreeNode(xFinal,startNode);

    //apply the boundary condition
    code=ParticleSphereInteraction(spec,ptr,xFinal,vFinal,dtTotal,(void*)newNode,BoundaryElement);

    if (code==_PARTICLE_DELETED_ON_THE_FACE_) {
      PIC::ParticleBuffer::DeleteParticle(ptr);
      return _PARTICLE_LEFT_THE_DOMAIN_;
    }
  }
  else {
    newNode=PIC::Mesh::mesh->findTreeNode(xFinal,startNode);
  }
#else 
  newNode=PIC::Mesh::mesh->findTreeNode(xFinal,startNode);
#endif //_TARGET_ == _TARGET_NONE_


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

    switch(_PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE_) {
    case _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__USER_FUNCTION_:
      code=ProcessOutsideDomainParticles(ptr,xInit,vInit,nIntersectionFace,newNode);
      break;
    case _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__PERIODIC_CONDITION_:
      exit(__LINE__,__FILE__,"Error: not implemented");
      break;
    case _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__SPECULAR_REFLECTION_:
      //reflect the particle back into the domain
    {
      double c=0.0;
      for (int idim=0;idim<3;idim++) c+=ExternalBoundaryFaceTable[nIntersectionFace].norm[idim]*vInit[idim];
      for (int idim=0;idim<3;idim++) vInit[idim]-=2.0*c*ExternalBoundaryFaceTable[nIntersectionFace].norm[idim];
    }

    code=_PARTICLE_REJECTED_ON_THE_FACE_;
    break;
    default:
      exit(__LINE__,__FILE__,"Error: the option is unknown");
      break;
    }



    memcpy(vFinal,vInit,3*sizeof(double));
    memcpy(xFinal,xInit,3*sizeof(double));
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
  if (_PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_) {
    PIC::ParticleTracker::RecordTrajectoryPoint(xFinal,vFinal,spec,ParticleData,(void*)newNode);

    if (_PIC_PARTICLE_TRACKER__TRACKING_CONDITION_MODE__DYNAMICS_ == _PIC_MODE_ON_) {
      PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(xFinal,vFinal,spec,ParticleData,(void*)newNode);
    }
  }

  if (_PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ == _PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ON_) {   
    //model the generic particle transformation
    int GenericParticleTransformationReturnCode,specInit=spec;

    GenericParticleTransformationReturnCode=_PIC_PARTICLE_MOVER__GENERIC_TRANSFORMATION_PROCESSOR_(xInit,xFinal,vFinal,spec,ptr,ParticleData,dtTotal,startNode);   //xInit,xFinal,vFinal,spec,ptr,ParticleData,dtMin,startNode

    if (GenericParticleTransformationReturnCode==_GENERIC_PARTICLE_TRANSFORMATION_CODE__PARTICLE_REMOVED_) {
      PIC::ParticleBuffer::DeleteParticle(ptr);
      return _PARTICLE_LEFT_THE_DOMAIN_;
    }

    //adjust the value of the dtLeft to match the time step for the species 'spec'
    if (spec!=specInit) {
      if (_PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_) {
        if (_PIC_PARTICLE_TRACKER__TRACKING_CONDITION_MODE__CHEMISTRY_ == _PIC_MODE_ON_) {
          PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(xFinal,vFinal,spec,ParticleData,(void*)startNode);
        }
      }
    }
  } //_PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ == _PIC_GENERIC_PARTICLE_TRANSFORMATION_MODE_ON_




  //finish the trajectory integration procedure
  PIC::Mesh::cDataBlockAMR *block;

  if (PIC::Mesh::mesh->FindCellIndex(xFinal,i,j,k,newNode,false)==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located");

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


//particle mover for the energy conserving scheme (Stefano Markidis et al., 2010, Mathematics and Computers in Simulation 80 (2010) 1509–1519
int PIC::Mover::Markidis2010(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *newNode=NULL;
  PIC::ParticleBuffer::byte *ParticleData;
  double vInit[3],xInit[3]={0.0,0.0,0.0},vFinal[3],xFinal[3],B[3],E[3];
  int idim,i,j,k,spec;

  ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
  PIC::ParticleBuffer::GetV(vInit,ParticleData);
  PIC::ParticleBuffer::GetX(xInit,ParticleData);
  spec=PIC::ParticleBuffer::GetI(ParticleData);

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

  static long int nCall=0;
  nCall++;

  //interpolate the fields acting upon on the particle
  PIC::CPLR::InitInterpolationStencil(xInit,startNode);
  PIC::CPLR::GetBackgroundElectricField(E);
  PIC::CPLR::GetBackgroundMagneticField(B);

  double v_prime[3],QdT_over_m,QdT_over_2m;

  QdT_over_m=PIC::MolecularData::GetElectricCharge(spec)*dtTotal/PIC::MolecularData::GetMass(spec);
  QdT_over_2m=0.5*QdT_over_m;

  //Eq 22
  for (idim=0;idim<3;idim++) v_prime[idim]=vInit[idim]+QdT_over_m*E[idim];

  //calculate the new valuw od the particle velocity. Eq 23
  double Denominator=1.0/(1.0+QdT_over_2m*QdT_over_2m*(B[0]*B[0]+B[1]*B[1]+B[2]*B[2]));
  double n1[3],n2;

  Vector3D::CrossProduct(n1,v_prime,B);
  n2=QdT_over_2m*QdT_over_2m*(v_prime[0]*B[0]+v_prime[1]*B[1]+v_prime[2]*B[2]);

  for (idim=0;idim<3;idim++) {
    vFinal[idim]=Denominator*(v_prime[idim]+QdT_over_2m*n1[idim]+n2*B[idim]);
    xFinal[idim]=xInit[idim]+dtTotal*vFinal[idim];
  }

  //interaction with the faces of the block and internal surfaces
  //check whether the particle trajectory is intersected the spherical body
#if  _TARGET_ID_(_TARGET_) != _TARGET_NONE__ID_ && _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_   
  double rFinal2;

  //if the particle is inside the sphere -> apply the boundary condition procedure
  if ((rFinal2=xFinal[0]*xFinal[0]+xFinal[1]*xFinal[1]+xFinal[2]*xFinal[2])<_RADIUS_(_TARGET_)*_RADIUS_(_TARGET_)) {
    double r=sqrt(rFinal2);
    int code;

    static cInternalSphericalData_UserDefined::fParticleSphereInteraction ParticleSphereInteraction=
        ((cInternalSphericalData*)(PIC::Mesh::mesh->InternalBoundaryList.front().BoundaryElement))->ParticleSphereInteraction;
    static void* BoundaryElement=PIC::Mesh::mesh->InternalBoundaryList.front().BoundaryElement;

    //move the particle location at the surface of the sphere
    for (idim=0;idim<DIM;idim++) xFinal[idim]*=_RADIUS_(_TARGET_)/r;

    //determine the block of the particle location
    newNode=PIC::Mesh::mesh->findTreeNode(xFinal,startNode);

    //apply the boundary condition
    code=ParticleSphereInteraction(spec,ptr,xFinal,vFinal,dtTotal,(void*)newNode,BoundaryElement);

    if (code==_PARTICLE_DELETED_ON_THE_FACE_) {
      PIC::ParticleBuffer::DeleteParticle(ptr);
      return _PARTICLE_LEFT_THE_DOMAIN_;
    }
  }
  else {
    newNode=PIC::Mesh::mesh->findTreeNode(xFinal,startNode);
  }
#else
  newNode=PIC::Mesh::mesh->findTreeNode(xFinal,startNode);
#endif //_TARGET_ == _TARGET_NONE_


  if (newNode==NULL) {
    //the particle left the computational domain
    int code=_PARTICLE_DELETED_ON_THE_FACE_;

    switch (_PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE_) {
    case _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__DELETE_:
      code=_PARTICLE_DELETED_ON_THE_FACE_;
      break;

    default:
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

      switch(_PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE_){
      case _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__USER_FUNCTION_:
        code=ProcessOutsideDomainParticles(ptr,xInit,vInit,nIntersectionFace,newNode);
        break;
      case _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__PERIODIC_CONDITION_:
        exit(__LINE__,__FILE__,"Error: not implemented");
        break;
      case _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__SPECULAR_REFLECTION_:
        //reflect the particle back into the domain
      {
        double c=0.0;
        for (int idim=0;idim<3;idim++) c+=ExternalBoundaryFaceTable[nIntersectionFace].norm[idim]*vInit[idim];
        for (int idim=0;idim<3;idim++) vInit[idim]-=2.0*c*ExternalBoundaryFaceTable[nIntersectionFace].norm[idim];
      }

      code=_PARTICLE_REJECTED_ON_THE_FACE_;
      break;
      default:
        exit(__LINE__,__FILE__,"Error: the option is unknown");
        break;
      }



      memcpy(vFinal,vInit,3*sizeof(double));
      memcpy(xFinal,xInit,3*sizeof(double));
      //       }
    } //_PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE_ == _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__DELETE

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
  if (_PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_) {
    PIC::ParticleTracker::RecordTrajectoryPoint(xFinal,vFinal,spec,ParticleData,(void*)newNode);

    if (_PIC_PARTICLE_TRACKER__TRACKING_CONDITION_MODE__DYNAMICS_ == _PIC_MODE_ON_) {
      PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(xFinal,vFinal,spec,ParticleData,(void*)newNode);
    }
  }

  //finish the trajectory integration procedure
  PIC::Mesh::cDataBlockAMR *block;

  if (PIC::Mesh::mesh->FindCellIndex(xFinal,i,j,k,newNode,false)==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located");

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


//particle mover that is used in iPIC3D (G. Lapenta/JournalofComputationalPhysics334(2017)349–366)
int PIC::Mover::Lapenta2017(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  PIC::ParticleBuffer::byte *ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr); 

  int threadId=0;

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  threadId = omp_get_thread_num();
#endif

  double *E_Corner=(PIC::Mover::E_Corner!=NULL) ? PIC::Mover::E_Corner[threadId] : NULL; 
  double *B_Corner=(PIC::Mover::B_Corner!=NULL) ? PIC::Mover::B_Corner[threadId] : NULL;
  double *B_Center=(PIC::Mover::B_Center!=NULL) ? PIC::Mover::B_Center[threadId] : NULL;

  cLapentaInputData data;

  data.E_Corner=E_Corner;

  switch (_PIC_FIELD_SOLVER_B_MODE_) {
  case _PIC_FIELD_SOLVER_B_CENTER_BASED_:
    data.B_C=B_Center;
    break;
  case _PIC_FIELD_SOLVER_B_CORNER_BASED_:
    data.B_C=B_Corner;
    break;
  }


  data.MolMass=PIC::MolecularData::MolMass;
  data.ElectricChargeTable=PIC::MolecularData::ElectricChargeTable;
  data.TimeStepTable=PIC::ParticleWeightTimeStep::GlobalTimeStep;

  data.ParticleDataLength=PIC::ParticleBuffer::ParticleDataLength;
  data.ParticleDataBuffer=PIC::ParticleBuffer::ParticleDataBuffer;
  data.mesh=PIC::Mesh::mesh;
  data.node=startNode;

  return Lapenta2017(ParticleData,ptr,&data); 
}

int PIC::Mover::Lapenta2017(PIC::ParticleBuffer::byte *ParticleData,long int ptr,cLapentaInputData *data) {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=data->node,*newNode=NULL;
  int idim,i,j,k,spec;
  double dtTotal;

#if _AVX_LAPENTA_MOVER_ == _OFF_
  double vInit[3],xInit[3],vFinal[3],xFinal[3];
#else
  union {__m256d vInit_v; double vInit[4];};
  union {__m256d xInit_v; double xInit[4];};
  union {__m256d vFinal_v; double vFinal[4];};
  union {__m256d xFinal_v; double xFinal[4];};
#endif

  PIC::ParticleBuffer::GetV(vInit,ParticleData);
  PIC::ParticleBuffer::GetX(xInit,ParticleData);
  spec=PIC::ParticleBuffer::GetI(ParticleData);


  switch (_SIMULATION_TIME_STEP_MODE_) {
  case _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_:
    dtTotal=data->TimeStepTable[spec];
    break;
  case _SINGLE_GLOBAL_TIME_STEP_:
    dtTotal=data->TimeStepTable[0]; 
    break;
  default:
    dtTotal=startNode->block->GetLocalTimeStep(spec);
  }

  static long int nCall=0;
  nCall++;


  //interpolate the fields acting upon on the particle at the NEW location of the particle (Appendix D, Eq 2)

  PIC::InterpolationRoutines::CornerBased::cStencil ElectricFieldStencil;
  PIC::InterpolationRoutines::CellCentered::cStencil MagneticFieldStencil;
  int threadId=0;

#if _AVX_LAPENTA_MOVER_ == _OFF_
  double E[4]={0.0,0.0,0.0,0.0},B[4]={0.0,0.0,0.0,0.0};
#else
  union {__m256d B_v; double B[4];};
  union {__m256d E_v; double E[4];};

  B_v=_mm256_setzero_pd();
  E_v=_mm256_setzero_pd();
#endif

  int *LocalCellID,Length;
  double *Weight;

  switch (_PIC_FIELD_SOLVER_MODE_) {
  case _PIC_FIELD_SOLVER_MODE__OFF_:
    PIC::CPLR::InitInterpolationStencil(xInit,startNode);
    PIC::CPLR::GetBackgroundElectricField(E);
    PIC::CPLR::GetBackgroundMagneticField(B);

    break;

  default:
#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
    threadId = omp_get_thread_num();
#endif

    switch ( _PIC_FIELD_SOLVER_MODE_) {
    case _PIC_FIELD_SOLVER_MODE__ELECTROMAGNETIC__ECSIM_:
      //interpolate the elecric field (corner nodes)
      PIC::InterpolationRoutines::CornerBased::InitStencil(xInit,startNode,ElectricFieldStencil);

      Length=ElectricFieldStencil.Length;
      LocalCellID=ElectricFieldStencil.LocalCellID;
      Weight=ElectricFieldStencil.Weight;

      for (int iStencil=0;iStencil<Length;iStencil++) {
        double *tempE1=data->E_Corner+3*LocalCellID[iStencil];

#if  _PIC_FIELD_SOLVER_B_MODE_== _PIC_FIELD_SOLVER_B_CORNER_BASED_
        double *tempB1=data->B_C+3*LocalCellID[iStencil];
#endif

#if _AVX_LAPENTA_MOVER_ == _OFF_
#pragma ivdep
        for (idim=0;idim<3;idim++) {
          E[idim]+=Weight[iStencil]*tempE1[idim];

#if  _PIC_FIELD_SOLVER_B_MODE_== _PIC_FIELD_SOLVER_B_CORNER_BASED_
          B[idim]+=Weight[iStencil]*tempB1[idim];
#endif
        }

#else //_AVX_INSTRUCTIONS_USAGE_MODE_
        __m256d w=_mm256_set1_pd(Weight[iStencil]);

#if  _PIC_FIELD_SOLVER_B_MODE_== _PIC_FIELD_SOLVER_B_CORNER_BASED_
        B_v=_mm256_fmadd_pd(w,_mm256_loadu_pd(tempB1),B_v);
#endif

        E_v=_mm256_fmadd_pd(w,_mm256_loadu_pd(tempE1),E_v);
#endif //_AVX_INSTRUCTIONS_USAGE_MODE_
      }

#if  _PIC_FIELD_SOLVER_B_MODE_== _PIC_FIELD_SOLVER_B_CENTER_BASED_
      //interpolate the magnetic field (center nodes)
      PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(xInit,startNode,MagneticFieldStencil);

      Length=MagneticFieldStencil.Length;
      LocalCellID=MagneticFieldStencil.LocalCellID;
      Weight=MagneticFieldStencil.Weight;

      for (int iStencil=0;iStencil<Length;iStencil++) {
        double *tempB1 = data->B_C+3*LocalCellID[iStencil];

#if _AVX_LAPENTA_MOVER_ == _OFF_
#pragma ivdep
        for (idim=0;idim<3;idim++) {
          B[idim]+=Weight[iStencil]*tempB1[idim];
        }
#else //_AVX_INSTRUCTIONS_USAGE_MODE_

        B_v=_mm256_fmadd_pd(_mm256_set1_pd(Weight[iStencil]),_mm256_loadu_pd(tempB1),B_v);
#endif //_AVX_INSTRUCTIONS_USAGE_MODE_
      }
#endif

      break;
    default:
      exit(__LINE__,__FILE__,"Error: unknown value of _PIC_FIELD_SOLVER_MODE_");
    }
  }

#ifdef __CUDA_ARCH__
  __syncwarp;
#endif

  E[3]=0.0,B[3]=0.0;  //the line is important when AVX is used. The index [3] is correct since B and B are defined compatible with __m256d

  //advance the particle velocity
  double QdT_over_m,QdT_over_2m,alpha[3][3];
  double c0,QdT_over_2m_squared,mass,chargeQ;
  double mass_conv=1.0,charge_conv=1.0;

  if (_PIC_FIELD_SOLVER_INPUT_UNIT_== _PIC_FIELD_SOLVER_INPUT_UNIT_NORM_) {
    mass_conv =1.0/_AMU_;
    charge_conv=1.0/ElectronCharge;
  }

  chargeQ = data->ElectricChargeTable[spec]*charge_conv; 
  mass= data->MolMass[spec]*mass_conv;
  QdT_over_m=chargeQ*dtTotal/mass;
  QdT_over_2m=0.5*QdT_over_m;
  QdT_over_2m_squared=QdT_over_2m*QdT_over_2m;


#if _AVX_LAPENTA_MOVER_ == _OFF_
  double BB[3][3],P[3];

  for (i=0;i<3;i++) {
    P[i]=-QdT_over_2m*B[i];

#pragma ivdep
    for (j=0;j<=i;j++) {
      BB[i][j]=QdT_over_2m_squared*B[i]*B[j];
      BB[j][i]=BB[i][j];
    }
  }

  c0=1.0/(1.0+QdT_over_2m_squared*(B[0]*B[0]+B[1]*B[1]+B[2]*B[2]));

  //Eq. D.3
  alpha[0][0]=c0*(1.0+BB[0][0]);
  alpha[0][1]=c0*(-P[2]+BB[0][1]);
  alpha[0][2]=c0*(P[1]+BB[0][2]);

  alpha[1][0]=c0*(P[2]+BB[1][0]);
  alpha[1][1]=c0*(1.0+BB[1][1]);
  alpha[1][2]=c0*(-P[0]+BB[1][2]);

  alpha[2][0]=c0*(-P[1]+BB[2][0]);
  alpha[2][1]=c0*(P[0]+BB[2][1]);
  alpha[2][2]=c0*(1.0+BB[2][2]);

  //Eq. D.8
#pragma ivdep
  for (idim=0;idim<3;idim++) {
    double vp=0.0;

    for (j=0;j<3;j++) vp+=alpha[idim][j]*(vInit[j]+QdT_over_2m*E[j]);

    vFinal[idim]=2.0*vp-vInit[idim];
  }

  //advance the particle location
#pragma ivdep
  for (idim=0;idim<3;idim++) xFinal[idim]=xInit[idim]+dtTotal*vFinal[idim];
#else //_AVX_INSTRUCTIONS_USAGE_MODE_

  union {__m256d P_v; double P[4];};
  union {__m256d B2_v; double B2[4];};

  P_v=_mm256_mul_pd(_mm256_set1_pd(-QdT_over_2m),B_v);
  B2_v=_mm256_mul_pd(B_v,B_v);

  c0=1.0/(1.0+QdT_over_2m_squared*(B2[0]+B2[1]+B2[2]));

  //Eq. D.3
  __m256d alpha_0v,alpha_1v,alpha_2v;
  __m256d c0_v,t_v;

  c0_v=_mm256_set1_pd(c0);
  alpha_0v=_mm256_mul_pd(c0_v,_mm256_fmadd_pd(_mm256_set1_pd(QdT_over_2m_squared*B[0]),B_v,_mm256_set_pd(0.0,P[1],-P[2],1.0)));
  alpha_1v=_mm256_mul_pd(c0_v,_mm256_fmadd_pd(_mm256_set1_pd(QdT_over_2m_squared*B[1]),B_v,_mm256_set_pd(0.0,-P[0],1.0,P[2])));
  alpha_2v=_mm256_mul_pd(c0_v,_mm256_fmadd_pd(_mm256_set1_pd(QdT_over_2m_squared*B[2]),B_v,_mm256_set_pd(0.0,1.0,P[0],-P[1])));


  union {__m256d vp_v; double vp_v_ptr[4];};
  union {__m256d tt_v; double tt_v_ptr[4];};

  t_v=_mm256_fmadd_pd(_mm256_set1_pd(QdT_over_2m),E_v,vInit_v);

  tt_v=_mm256_mul_pd(alpha_0v,t_v);
  vp_v_ptr[0]=tt_v_ptr[0]+tt_v_ptr[1]+tt_v_ptr[2];

  tt_v=_mm256_mul_pd(alpha_1v,t_v);
  vp_v_ptr[1]=tt_v_ptr[0]+tt_v_ptr[1]+tt_v_ptr[2];


  tt_v=_mm256_mul_pd(alpha_2v,t_v);
  vp_v_ptr[2]=tt_v_ptr[0]+tt_v_ptr[1]+tt_v_ptr[2];

  vFinal_v=_mm256_fmsub_pd(_mm256_set1_pd(2.0),vp_v,vInit_v);
  xFinal_v=_mm256_fmadd_pd(_mm256_set1_pd(dtTotal),vFinal_v,xInit_v);
#endif //_AVX_INSTRUCTIONS_USAGE_MODE_

  newNode=data->mesh->findTreeNode(xFinal,startNode);

  //interaction with the faces of the block and internal surfaces
  //check whether the particle trajectory is intersected the spherical body
#if  _TARGET_ID_(_TARGET_) != _TARGET_NONE__ID_ && _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_ON_ 
  double rFinal2;

  //if the particle is inside the sphere -> apply the boundary condition procedure
  if ((rFinal2=xFinal[0]*xFinal[0]+xFinal[1]*xFinal[1]+xFinal[2]*xFinal[2])<_RADIUS_(_TARGET_)*_RADIUS_(_TARGET_)) {
    double r=sqrt(rFinal2);
    int code;

    static cInternalSphericalData_UserDefined::fParticleSphereInteraction ParticleSphereInteraction=
        ((cInternalSphericalData*)(PIC::Mesh::mesh->InternalBoundaryList.front().BoundaryElement))->ParticleSphereInteraction;
    static void* BoundaryElement=PIC::Mesh::mesh->InternalBoundaryList.front().BoundaryElement;

    //move the particle location at the surface of the sphere
    for (idim=0;idim<DIM;idim++) xFinal[idim]*=_RADIUS_(_TARGET_)/r;

    //determine the block of the particle location
    newNode=data->mesh->findTreeNode(xFinal,startNode);

    //apply the boundary condition
    code=ParticleSphereInteraction(spec,ptr,xFinal,vFinal,dtTotal,(void*)newNode,BoundaryElement);

    if (code==_PARTICLE_DELETED_ON_THE_FACE_) {
      PIC::ParticleBuffer::DeleteParticle(ptr);
      return _PARTICLE_LEFT_THE_DOMAIN_;
    }
  }
  else {
    newNode=data->mesh->findTreeNode(xFinal,startNode);
  }
#else
  newNode=data->mesh->findTreeNode(xFinal,startNode);
#endif //_TARGET_ == _TARGET_NONE_


  if (newNode==NULL) {
    //the particle left the computational domain
    int code=_PARTICLE_DELETED_ON_THE_FACE_;

    switch (_PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE_) {

    case _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__DELETE_:
      code=_PARTICLE_DELETED_ON_THE_FACE_;
      break;
    default:
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

      newNode=data->mesh->findTreeNode(xInit,startNode);

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

        newNode=data->mesh->findTreeNode(xInit,startNode);

        if (newNode==NULL) exit(__LINE__,__FILE__,"Error: cannot find the node");
      }

      switch(_PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE_) {
      case _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__USER_FUNCTION_:
        code = ProcessOutsideDomainParticles(ptr, xInit, vInit, nIntersectionFace, newNode);
        break;
      case _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__PERIODIC_CONDITION_:
        exit(__LINE__,__FILE__,"Error: not implemented");
        break;
      case _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__SPECULAR_REFLECTION_:
        //reflect the particle back into the domain
      {
        double c = 0.0;
        for (int idim = 0; idim < 3; idim++) c += ExternalBoundaryFaceTable[nIntersectionFace].norm[idim] * vInit[idim];
        for (int idim = 0; idim < 3; idim++) vInit[idim] -= 2.0 * c * ExternalBoundaryFaceTable[nIntersectionFace].norm[idim];
      }
      code = _PARTICLE_REJECTED_ON_THE_FACE_;
      break;
      default:
        exit(__LINE__,__FILE__,"Error: the option is unknown");
        break;
      }



      memcpy(vFinal,vInit,3*sizeof(double));
      memcpy(xFinal,xInit,3*sizeof(double));
      //       }
    }//_PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE_ == _PIC_PARTICLE_DOMAIN_BOUNDARY_INTERSECTION_PROCESSING_MODE__DELETE

    //call the function that process particles that leaved the coputational domain
    switch (code) {
    case _PARTICLE_DELETED_ON_THE_FACE_:
      PIC::ParticleBuffer::DeleteParticle(ptr);
      return _PARTICLE_LEFT_THE_DOMAIN_;
    default:
      exit(__LINE__,__FILE__,"Error: not implemented");
    }
  }
  else {
    //at this point newNode!=NULL
    if (newNode->IsUsedInCalculationFlag==false) {
      PIC::ParticleBuffer::DeleteParticle(ptr);
      return _PARTICLE_IN_NOT_IN_USE_NODE_;
    }
  }

  //save the trajectory point
  if (_PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_) {
    PIC::ParticleTracker::RecordTrajectoryPoint(xFinal,vFinal,spec,ParticleData,(void*)newNode);

    if (_PIC_PARTICLE_TRACKER__TRACKING_CONDITION_MODE__DYNAMICS_ == _PIC_MODE_ON_) {
      PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(xFinal,vFinal,spec,ParticleData,(void*)newNode);
    }
  }

  //finish the trajectory integration procedure
  PIC::Mesh::cDataBlockAMR *block;

  if (data->mesh->FindCellIndex(xFinal,i,j,k,newNode,false)==-1) exit(__LINE__,__FILE__,"Error: cannot find the cellwhere the particle is located");



  if ((block=newNode->block)==NULL) { // case newNode->IsUsedInCalculationFlag==false is already considered above
    if (_PIC_MOVER__UNKNOWN_ERROR_IN_PARTICLE_MOTION__STOP_EXECUTION_ == _PIC_MODE_OFF_) {
      double Rate;

      Rate=startNode->block->GetLocalParticleWeight(spec)*PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ptr)/
          startNode->block->GetLocalTimeStep(spec);

      PIC::Mover::Sampling::Errors::AddRemovedParticleData(Rate,spec,__LINE__,__FILE__);

      //remove the particle
      PIC::ParticleBuffer::DeleteParticle(ptr);
      return _PARTICLE_LEFT_THE_DOMAIN_;
    }
    else exit(__LINE__,__FILE__,"Error: the block is empty. Most probably the time step is too large");

  }

#ifdef __CUDA_ARCH__
  PIC::ParticleBuffer::SetPrev(-1,ParticleData);


  int tptr=ptr;
  int *source=(int*)(block->tempParticleMovingListTable+i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k));

  long int tempFirstCellParticle=atomicExch(source,tptr);

  if (sizeof(long int )>sizeof(int)) {
    *(source+1)=0;
  }

  PIC::ParticleBuffer::SetNext(tempFirstCellParticle,ParticleData);
  if (tempFirstCellParticle!=-1) PIC::ParticleBuffer::SetPrev(ptr,_GetParticleDataPointer(tempFirstCellParticle,data->ParticleDataLength,data->ParticleDataBuffer));



#elif _PIC_MOVER__MPI_MULTITHREAD_ == _PIC_MODE_ON_
  PIC::ParticleBuffer::SetPrev(-1,ParticleData); 

  long int tempFirstCellParticle=atomic_exchange(block->tempParticleMovingListTable+i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k),ptr); 

  PIC::ParticleBuffer::SetNext(tempFirstCellParticle,ParticleData);
  if (tempFirstCellParticle!=-1) PIC::ParticleBuffer::SetPrev(ptr,_GetParticleDataPointer(tempFirstCellParticle,data->ParticleDataLength,data->ParticleDataBuffer));

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


#if _PIC_TEMP_PARTICLE_LIST_MODE_ == _PIC_TEMP_PARTICLE_LIST_MODE__SHARED_
  PIC::Mesh::cDataCenterNode *CenterNode=block->GetCenterNode(i,j,k);

  while (CenterNode->lock_associated_data.test_and_set(std::memory_order_acquire)==true);
#endif


  PIC::ParticleBuffer::SetNext(ThreadTempParticleMovingData->first,ParticleData);
  PIC::ParticleBuffer::SetPrev(-1,ParticleData);

  if (ThreadTempParticleMovingData->last==-1) ThreadTempParticleMovingData->last=ptr;
  if (ThreadTempParticleMovingData->first!=-1) PIC::ParticleBuffer::SetPrev(ptr,ThreadTempParticleMovingData->first);
  ThreadTempParticleMovingData->first=ptr;

#if _PIC_TEMP_PARTICLE_LIST_MODE_ == _PIC_TEMP_PARTICLE_LIST_MODE__SHARED_
  CenterNode->lock_associated_data.clear(std::memory_order_release);
#endif


#else
#error The option is unknown
#endif



  PIC::ParticleBuffer::SetV(vFinal,ParticleData);
  PIC::ParticleBuffer::SetX(xFinal,ParticleData);



  return _PARTICLE_MOTION_FINISHED_;
}














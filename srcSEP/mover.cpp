/*
 * mover.cpp
 *
 *  Created on: May 16, 2020
 *      Author: vtenishe
 */

#include "sep.h"

int SEP::ParticleMover_default(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  double xInit[3];
  int res;

  PIC::ParticleBuffer::GetX(xInit,ptr);

  switch (SEP::ParticleTrajectoryCalculation) {
  case SEP::ParticleTrajectoryCalculation_RelativisticBoris:
    res=PIC::Mover::Relativistic::Boris(ptr,dtTotal,startNode);
    break;
  case ParticleTrajectoryCalculation_FieldLine:
    res=PIC::Mover::FieldLine::Mover_SecondOrder(ptr,dtTotal,startNode);
    break;
  }

  if (res==_PARTICLE_MOTION_FINISHED_) {
    //appply particle scatering model if needed 
    double x[3],v[3],speed; 
    double mean_free_path,l2,e,r;
    int idim,spec;

    PIC::ParticleBuffer::GetX(x,ptr);
    PIC::ParticleBuffer::GetV(v,ptr);
    spec=PIC::ParticleBuffer::GetI(ptr);
  
    speed=Vector3D::Length(v);
    e=Relativistic::Speed2E(speed,PIC::MolecularData::GetMass(spec));
    r=Vector3D::Length(x); 
  
    mean_free_path=3.4E9; 

    for (l2=0.0,idim=0;idim<3;idim++) {
      double t;

      t=xInit[idim]-x[idim];
      l2+=t*t;
    }

    if ((1.0-exp(-sqrt(l2)/mean_free_path))>rnd()) {
      //the scattering even occured
      Vector3D::Distribution::Uniform(v,speed);
   
      PIC::ParticleBuffer::SetV(v,ptr);
    }
  }

  return res;
} 

int SEP::ParticleMover__He_2019_AJL(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  namespace PB = PIC::ParticleBuffer; 

  double *x,*v,mu,SwU[3]={0.0,0.0,0.0},b[3]={0.0,0.0,0.0},p;
  PB::byte* ParticleData;
  int spec;   

  static int ncall=0;
  ncall++;
   
  
  ParticleData=PB::GetParticleDataPointer(ptr);
  v=PB::GetV(ParticleData);
  x=PB::GetX(ParticleData); 
  spec=PB::GetI(ParticleData);


  //get magnetic field, and plasma velocity at location of the particle
  PIC::InterpolationRoutines::CellCentered::cStencil CenterBasedStencil;
  char *AssociatedDataPointer;
  double weight;
  double *b_ptr,*u_ptr;
  int idim;

  PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(x,node,CenterBasedStencil);    

  for (int icell=0; icell<CenterBasedStencil.Length; icell++) {
    AssociatedDataPointer=CenterBasedStencil.cell[icell]->GetAssociatedDataBufferPointer();
    weight=CenterBasedStencil.Weight[icell];
  
    b_ptr=(double*)(AssociatedDataPointer+PIC::CPLR::SWMF::MagneticFieldOffset);
    u_ptr=(double*)(AssociatedDataPointer+PIC::CPLR::SWMF::BulkVelocityOffset); 
    
    for (idim=0;idim<3;idim++) {
      SwU[idim]+=weight*u_ptr[idim];
      b[idim]+=weight*b_ptr[idim];
    }
  } 

  //move to the frame of reference related to the solar wind
  for (idim=0;idim<3;idim++) v[idim]-=SwU[idim];

  //determine mu
  double AbsB=Vector3D::Normalize(b);
  mu=Vector3D::DotProduct(v,b)/Vector3D::Length(v);

  if (AbsB==0.0) {
   PIC::ParticleBuffer::DeleteParticle(ptr);
   return _PARTICLE_DELETED_ON_THE_FACE_;
  }

  //create a coordinate frame where 'x' and 'y' are orthogonal to 'b'
  double e0[3],e1[3];

  Vector3D::GetRandomNormFrame(e0,e1,b);

  //determine the spatial step for calculating the derivatives 
  double dxStencil,t; 

  dxStencil=(node->xmax[0]-node->xmin[0])/_BLOCK_CELLS_X_;
  dxStencil=min(dxStencil,node->xmax[1]-node->xmin[1])/_BLOCK_CELLS_Y_; 
  dxStencil=min(dxStencil,node->xmax[2]-node->xmin[2])/_BLOCK_CELLS_Z_;

   
  const int _Velocity=0;
  const int _MagneticField=1;
 


  auto GetShiftedSwSpeed = [&] (int iVelocityVectorIndex, int VariableId,int iShiftDimension,double dir) {
    //get location of the test point 
    double x_test[3];
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node_test; 

    switch (iShiftDimension) {
    case 2:
      for (idim=0;idim<3;idim++) x_test[idim]=x[idim]+dir*dxStencil*b[idim];
      break;
    case 0:
      for (idim=0;idim<3;idim++) x_test[idim]=x[idim]+dir*dxStencil*e0[idim];
      break;
    case 1: 
      for (idim=0;idim<3;idim++) x_test[idim]=x[idim]+dir*dxStencil*e1[idim];
      break;
    }

    node_test=PIC::Mesh::mesh->findTreeNode(x_test,node); 

    bool use_original_point=false;

    if (node_test==NULL) use_original_point=true;
    else if (node_test->block==NULL) use_original_point=true;

    if (use_original_point==true) {
      node_test=node;

      for (idim=0;idim<3;idim++) x_test[idim]=x[idim];
    }

    //get velocity of solar wind 
    double res=0.0,AbsB=0.0;

    char *AssociatedDataPointer;
    double weight,t;
    double *u_ptr;
    int idim,i,j,k;

    PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(x_test,node_test,CenterBasedStencil);

    for (int icell=0; icell<CenterBasedStencil.Length; icell++) {
      AssociatedDataPointer=CenterBasedStencil.cell[icell]->GetAssociatedDataBufferPointer();
      weight=CenterBasedStencil.Weight[icell];

      switch (VariableId) {
      case _Velocity:
        u_ptr=(double*)(AssociatedDataPointer+PIC::CPLR::SWMF::BulkVelocityOffset);

        switch (iVelocityVectorIndex) {
        case 0:
          res+=weight*Vector3D::DotProduct(u_ptr,e0);
          break;
        case 1:
          res+=weight*Vector3D::DotProduct(u_ptr,e1);
          break; 
        case 2:
          res+=weight*Vector3D::DotProduct(u_ptr,b);
          break;
        }

        break;
   
      case _MagneticField:
        t=0.0;

        for (int idim=0;idim<3;idim++) {
          double t0=((double*)(AssociatedDataPointer+PIC::CPLR::SWMF::MagneticFieldOffset))[idim];
          t+=t0*t0;
        }

        res+=weight*sqrt(t);
      }
    }

    return res;
  };

  //calculate the derivatives
  double dUx_dx=(GetShiftedSwSpeed(0,_Velocity,0,1.0)-GetShiftedSwSpeed(0,_Velocity,0,-1.0))/(2.0*dxStencil);   
  double dUy_dy=(GetShiftedSwSpeed(1,_Velocity,1,1.0)-GetShiftedSwSpeed(1,_Velocity,1,-1.0))/(2.0*dxStencil);
  double dUz_dz=(GetShiftedSwSpeed(2,_Velocity,2,1.0)-GetShiftedSwSpeed(2,_Velocity,2,-1.0))/(2.0*dxStencil);

  double dB_dz=(GetShiftedSwSpeed(0,_MagneticField,2,1.0)-GetShiftedSwSpeed(0,_MagneticField,2,-1.0))/(2.0*dxStencil);
  
  //calculate the updated value of 'p'
  double p_v[3];
  double m0=PIC::MolecularData::GetMass(spec);
  double mu2=mu*mu;
  double AbsV=Vector3D::Length(v);

  switch (_PIC_PARTICLE_MOVER__RELATIVITY_MODE_) {
  case _PIC_MODE_OFF_:
    p=AbsV*m0;
    break;
  case _PIC_MODE_ON_:
    p=Relativistic::Speed2Momentum(AbsV,m0);
    break;
  }

  p-=dtTotal*p*( (1.0-mu2)/2.0*(dUx_dx+dUy_dy) + mu2*dUz_dz); 
  mu+=dtTotal*(1.0-mu2)*0.5*( -AbsV/AbsB*dB_dz+ mu*(dUx_dx+dUy_dy-2.0*dUz_dz) ); 

  if (mu>1.0) mu=1.0;
  if (mu<-1.0) mu=-1.0;

  //determine the final location of the particle
  for (idim=0;idim<3;idim++) {
    x[idim]+=(AbsV*b[idim]+SwU[idim])*dtTotal; 
  }

  //determine the final velocity of the particle in the frame related to the Sun 
  switch (_PIC_PARTICLE_MOVER__RELATIVITY_MODE_) {
  case _PIC_MODE_OFF_:
    AbsV=p/m0;
    break;
  case _PIC_MODE_ON_:
    AbsV=Relativistic::Momentum2Speed(p,m0);
    break;
  }

  double sin_mu=sqrt(1.0-mu*mu);

  for (idim=0;idim<3;idim++) {
    v[idim]=AbsV*(b[idim]*mu+e0[idim]*sin_mu)+SwU[idim];  
  }

  //check whether the particle is still in the domain
  node=PIC::Mesh::mesh->findTreeNode(x,node); 

  if (node==NULL) {
    //the particle left the computational domain
    PIC::ParticleBuffer::DeleteParticle(ptr);
    return _PARTICLE_DELETED_ON_THE_FACE_;
  }
  else if (node->IsUsedInCalculationFlag==false) {
    //the particle left the computational domain
    PIC::ParticleBuffer::DeleteParticle(ptr);
    return _PARTICLE_DELETED_ON_THE_FACE_;
  }  
  else {
    if (node->block==NULL) exit(__LINE__,__FILE__,"Error: node->block==NULL -> the time step is too large");
  }

  //check that the new particle location is outside of the Sun and Earth
  if (x[0]*x[0]+x[1]*x[1]+x[2]*x[2]<_RADIUS_(_SUN_)*_RADIUS_(_SUN_)) {
    //the particle is inside the Sun -> remove it
    PIC::ParticleBuffer::DeleteParticle(ptr);
    return _PARTICLE_LEFT_THE_DOMAIN_;
  }
 
  
  //save the trajectory point
  if (_PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_) { 
    PIC::ParticleTracker::RecordTrajectoryPoint(x,v,spec,ParticleData,(void*)node);

    if (_PIC_PARTICLE_TRACKER__TRACKING_CONDITION_MODE__DYNAMICS_ == _PIC_MODE_ON_) { 
      PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(x,v,spec,ParticleData,(void*)node);
    }
  }

  //determine the cell where the particle is located
  int i,j,k;

  PIC::Mesh::mesh->FindCellIndex(x,i,j,k,node,false);

  //update particle lists
  #if _PIC_MOVER__MPI_MULTITHREAD_ == _PIC_MODE_ON_
  PIC::ParticleBuffer::SetPrev(-1,ParticleData);

  long int tempFirstCellParticle=atomic_exchange(block->tempParticleMovingListTable+i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k),ptr);

  PIC::ParticleBuffer::SetNext(tempFirstCellParticle,ParticleData);
  if (tempFirstCellParticle!=-1) PIC::ParticleBuffer::SetPrev(ptr,tempFirstCellParticle);

  #elif _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
  long int tempFirstCellParticle,*tempFirstCellParticlePtr;

  tempFirstCellParticlePtr=node->block->tempParticleMovingListTable+i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k);
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

   
  return _PARTICLE_MOTION_FINISHED_;
}






    
 
int SEP::ParticleMover_Kartavykh_2016_AJ(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  namespace PB = PIC::ParticleBuffer;


  if (PIC::CPLR::SWMF::CouplingTime_last<0.0) {
    return ParticleMover__He_2019_AJL(ptr,dtTotal,node);
  }
  else if (_PIC_SWMF_COUPLER__SAVE_TWO_DATA_SETS_!=_PIC_MODE_ON_) {
    exit(__LINE__,__FILE__,"For this mover _PIC_SWMF_COUPLER__SAVE_TWO_DATA_SETS_ should be _PIC_MODE_ON_");
  } 

  double *x,*v,mu,mu2,p,m0,AbsV;
  PB::byte* ParticleData;
  int spec,idim;

  static int ncall=0;
  ncall++;

  ParticleData=PB::GetParticleDataPointer(ptr);
  v=PB::GetV(ParticleData);
  x=PB::GetX(ParticleData);
  spec=PB::GetI(ParticleData);
  m0=PIC::MolecularData::GetMass(spec); 

  //get the interpolation stencils that will be used for calculating the derivatives
  PIC::InterpolationRoutines::CellCentered::cStencil xp_stencil;
  PIC::InterpolationRoutines::CellCentered::cStencil x0_plus_stencil,x0_minus_stencil;
  PIC::InterpolationRoutines::CellCentered::cStencil x1_plus_stencil,x1_minus_stencil;
  PIC::InterpolationRoutines::CellCentered::cStencil x2_plus_stencil,x2_minus_stencil;


  const int _dx_plus=0;
  const int _dx_minus=1;

  PIC::InterpolationRoutines::CellCentered::cStencil stencil_table[3][2];

  //pointing vectors the would be used for calcualting the derivatives
  double b[3],e0[3],e1[3];

  //determine the spatial step for calculating the derivatives 
  double dxStencil;
  
  dxStencil=(node->xmax[0]-node->xmin[0])/_BLOCK_CELLS_X_;
  dxStencil=min(dxStencil,node->xmax[1]-node->xmin[1])/_BLOCK_CELLS_Y_;
  dxStencil=min(dxStencil,node->xmax[2]-node->xmin[2])/_BLOCK_CELLS_Z_;
  
  PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(x,node,xp_stencil);

  //prepare the stencil table
  for (int idim=0;idim<3;idim++) for (int idir=-1;idir<=1;idir+=2) {
    double x_test[3];
    int i;
    PIC::InterpolationRoutines::CellCentered::cStencil stencil_test;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node_test;

    for (i=0;i<3;i++) x_test[i]=x[i];

    x_test[idim]+=idir*dxStencil;

    node_test=PIC::Mesh::mesh->findTreeNode(x_test,node);

    if (node_test->IsUsedInCalculationFlag==false) {
      stencil_test=xp_stencil;   
    }
    else {
      PIC::InterpolationRoutines::CellCentered::Linear::InitStencil(x_test,node_test,stencil_test);
    }
     
    switch (idir) {
    case 1:
      stencil_table[idim][_dx_plus]=stencil_test;
      break;
    case -1:
      stencil_table[idim][_dx_minus]=stencil_test;
      break; 
    }
  }

  //calculate the interpoalted value
  auto CalculateInterpolatedVector = [&] (double* res, int length,int offset,PIC::InterpolationRoutines::CellCentered::cStencil& Stencil) {
    char *AssociatedDataPointer;
    double weight,*ptr;
    int i;

    for (i=0;i<length;i++) res[i]=0.0; 

    for (int icell=0; icell<Stencil.Length; icell++) {
      AssociatedDataPointer=Stencil.cell[icell]->GetAssociatedDataBufferPointer();
      weight=Stencil.Weight[icell];
      ptr=(double*)(offset+AssociatedDataPointer);

      for (i=0;i<length;i++) res[i]+=weight*ptr[i];
    }
  }; 

  auto CalculateInterpolatedValue = [&] (int offset,PIC::InterpolationRoutines::CellCentered::cStencil& Stencil) {
    double res;

    CalculateInterpolatedVector(&res,1,offset,Stencil);
    return res;
  };


  //generate a random frame of reference 
  double Usw[3],AbsB;
  
  CalculateInterpolatedVector(b,3,PIC::CPLR::SWMF::MagneticFieldOffset,xp_stencil);
  AbsB=Vector3D::Normalize(b);

  if (AbsB==0.0) {
    PIC::ParticleBuffer::DeleteParticle(ptr);
    return _PARTICLE_DELETED_ON_THE_FACE_;
  }

  //create a coordinate frame where 'x' and 'y' are orthogonal to 'b'
  Vector3D::GetRandomNormFrame(e0,e1,b);
  
  //all further calculations are done in the frame (e0,e1,b):
  double dU0_dx0,dU1_dx1,dU2_dx2,dU_dt[3],U[3],B_shifted[3],AbsB_shifted,U_shifted[3]; 

  //time detivative
  double U_last[3];

  CalculateInterpolatedVector(U,3,PIC::CPLR::SWMF::BulkVelocityOffset,xp_stencil);
  CalculateInterpolatedVector(U_last,3,PIC::CPLR::SWMF::BulkVelocityOffset_last,xp_stencil);

  //move to the frame of reference moving with the ambient plasma
  for (idim=0;idim<3;idim++) v[idim]-=U[idim];

  AbsV=Vector3D::Length(v);
  mu=Vector3D::DotProduct(v,b)/AbsV;
  mu2=mu*mu;

  switch (_PIC_PARTICLE_MOVER__RELATIVITY_MODE_) {
  case _PIC_MODE_OFF_:
    p=AbsV*m0;
    break;
  case _PIC_MODE_ON_:
    p=Relativistic::Speed2Momentum(AbsV,m0);
    break;
  }


  //parameters 
  double stencil_multiplyer=1.0/(2.0*dxStencil);

  double db_i__dx_i=0.0; //\frac{\partial b}{\partial b}
  double dU_i__dx_i=0.0; //\frac{\partial U_i}{\partial x_i}
  double b_i__b_j__dU_i__dx_j=0.0; //b_i*b_j\frac{\partial U_i}{\partial x_j}
  double b_i__du_i__dt=0.0; //b_i* \frac{\partial U_i}{\partial t} 
  double b_i__U_j__dU_i__dx_j=0.0; //b_i*  U_j\frac{\partial U_i}{\partial x_j} 

  double B_shifted_plus[3],B_shifted_minus[3],AbsB_shifted_plus,AbsB_shifted_minus;
  double U_shifted_plus[3],U_shifted_minus[3];

  double swmf_coupling_time_interval=PIC::CPLR::SWMF::CouplingTime-PIC::CPLR::SWMF::CouplingTime_last;

  for (int j=0;j<3;j++) {
    //magnetic field 
    CalculateInterpolatedVector(B_shifted_plus,3,PIC::CPLR::SWMF::MagneticFieldOffset,stencil_table[j][_dx_plus]);        
    CalculateInterpolatedVector(B_shifted_minus,3,PIC::CPLR::SWMF::MagneticFieldOffset,stencil_table[j][_dx_minus]); 

    AbsB_shifted_plus=Vector3D::Normalize(B_shifted_plus);
    AbsB_shifted_minus=Vector3D::Normalize(B_shifted_minus);

    if ((AbsB_shifted_plus==0.0)||(AbsB_shifted_minus==0.0)) {
      PIC::ParticleBuffer::DeleteParticle(ptr);
      return _PARTICLE_DELETED_ON_THE_FACE_;
    }      

    db_i__dx_i+=stencil_multiplyer*(B_shifted_plus[j]-B_shifted_minus[j]);

    //plasma bulk velocity
    CalculateInterpolatedVector(U_shifted_plus,3,PIC::CPLR::SWMF::BulkVelocityOffset,stencil_table[j][_dx_plus]); 
    CalculateInterpolatedVector(U_shifted_minus,3,PIC::CPLR::SWMF::BulkVelocityOffset,stencil_table[j][_dx_minus]);

    dU_i__dx_i+=stencil_multiplyer*(U_shifted_plus[j]-U_shifted_minus[j]); 

    b_i__b_j__dU_i__dx_j+=stencil_multiplyer*b[j]*(
      (U_shifted_plus[0]-U_shifted_minus[0])*b[0]+
      (U_shifted_plus[1]-U_shifted_minus[1])*b[1]+
      (U_shifted_plus[2]-U_shifted_minus[2])*b[2]); 

    b_i__U_j__dU_i__dx_j+=stencil_multiplyer*U[j]*(
      (U_shifted_plus[0]-U_shifted_minus[0])*b[0]+
      (U_shifted_plus[1]-U_shifted_minus[1])*b[1]+
      (U_shifted_plus[2]-U_shifted_minus[2])*b[2]);  

    //time derivative: p3 
    b_i__du_i__dt+=(b[j]*(U[j]-U_last[j]))/swmf_coupling_time_interval;
  }


  double dmu_dt,dp_dt;

  dmu_dt=(1.0-mu2)/2.0*(
    AbsV*db_i__dx_i+
    mu*dU_i__dx_i-3.0*mu*b_i__b_j__dU_i__dx_j- 
    2.0/AbsV*(b_i__du_i__dt+b_i__U_j__dU_i__dx_j));

  dp_dt=p*(
    (1.0-3.0*mu2)/2.0*b_i__b_j__dU_i__dx_j-
    (1.0-mu2)/2.0*dU_i__dx_i-
    mu/AbsV*(b_i__du_i__dt+b_i__U_j__dU_i__dx_j));

  //calculate the updated value of 'mu' and 'p'
  p+=dtTotal*dp_dt;
  mu+=dtTotal*dmu_dt;

  if (mu<-1.0) mu=-1.0;
  if (mu>1.0) mu=1.0;

  //determine the final location of the particle
  if (_SEP_MOVER_DRIFT_==_PIC_MODE_OFF_) { 
    for (idim=0;idim<3;idim++) {
      x[idim]+=(AbsV*b[idim]+U[idim])*dtTotal;
    }
  }
  else {
    double v_drift[3],t[3];
    double v_parallel,v_perp;
    double ElectricCharge=PIC::MolecularData::GetElectricCharge(spec);

    v_parallel=Vector3D::DotProduct(v,b);

    memcpy(t,v,3*sizeof(double));
    Vector3D::Orthogonalize(b,t); 
    v_perp=Vector3D::Length(t); 

    GetDriftVelocity(v_drift,x,v_parallel,v_perp,ElectricCharge,m0,node);

    for (idim=0;idim<3;idim++) {
      x[idim]+=(AbsV*b[idim]+U[idim]+v_drift[idim])*dtTotal;
    }
  }


  //determine the final velocity of the particle in the frame related to the Sun
  switch (_PIC_PARTICLE_MOVER__RELATIVITY_MODE_) {
  case _PIC_MODE_OFF_:
    AbsV=p/m0;
    break;
  case _PIC_MODE_ON_:
    AbsV=Relativistic::Momentum2Speed(p,m0);
    break;
  }

  double sin_mu=sqrt(1.0-mu*mu);

  for (idim=0;idim<3;idim++) {
    v[idim]=AbsV*(b[idim]*mu+e0[idim]*sin_mu)+U[idim];
  }


  //check whether the particle is still in the domain
  node=PIC::Mesh::mesh->findTreeNode(x,node);

  if (node==NULL) {
    //the particle left the computational domain
    PIC::ParticleBuffer::DeleteParticle(ptr);
    return _PARTICLE_DELETED_ON_THE_FACE_;
  }
  else if (node->IsUsedInCalculationFlag==false) {
    //the particle left the computational domain
    PIC::ParticleBuffer::DeleteParticle(ptr);
    return _PARTICLE_DELETED_ON_THE_FACE_;
  }
  else {
    if (node->block==NULL) exit(__LINE__,__FILE__,"Error: node->block==NULL -> the time step is too large");
  }

  //check that the new particle location is outside of the Sun and Earth
  if (x[0]*x[0]+x[1]*x[1]+x[2]*x[2]<_RADIUS_(_SUN_)*_RADIUS_(_SUN_)) {
    //the particle is inside the Sun -> remove it
    PIC::ParticleBuffer::DeleteParticle(ptr);
    return _PARTICLE_LEFT_THE_DOMAIN_;
  }

  //save the trajectory point
  if (_PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_) {
    PIC::ParticleTracker::RecordTrajectoryPoint(x,v,spec,ParticleData,(void*)node);

    if (_PIC_PARTICLE_TRACKER__TRACKING_CONDITION_MODE__DYNAMICS_ == _PIC_MODE_ON_) {
      PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(x,v,spec,ParticleData,(void*)node);
    }
  }


  //determine the cell where the particle is located
  int i,j,k;

  PIC::Mesh::mesh->FindCellIndex(x,i,j,k,node,false);

  //update particle lists
  #if _PIC_MOVER__MPI_MULTITHREAD_ == _PIC_MODE_ON_
  PIC::ParticleBuffer::SetPrev(-1,ParticleData);

  long int tempFirstCellParticle=atomic_exchange(block->tempParticleMovingListTable+i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k),ptr);

  PIC::ParticleBuffer::SetNext(tempFirstCellParticle,ParticleData);
  if (tempFirstCellParticle!=-1) PIC::ParticleBuffer::SetPrev(ptr,tempFirstCellParticle);

  #elif _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
  long int tempFirstCellParticle,*tempFirstCellParticlePtr;

  tempFirstCellParticlePtr=node->block->tempParticleMovingListTable+i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k);
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


  return _PARTICLE_MOTION_FINISHED_;
}




int SEP::ParticleMover_Droge_2009_AJ(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  namespace PB = PIC::ParticleBuffer;
  namespace FL = PIC::FieldLine;

  PIC::ParticleBuffer::byte *ParticleData;
  double mu,AbsB,L,vParallel,vNormal,v,DivAbsB;
  double FieldLineCoord;
  int iFieldLine,spec;
  FL::cFieldLineSegment *Segment; 

  ParticleData=PB::GetParticleDataPointer(ptr);

  FieldLineCoord=PB::GetFieldLineCoord(ParticleData);
  iFieldLine=PB::GetFieldLineId(ParticleData); 
  spec=PB::GetI(ParticleData);

  vParallel=PB::GetVParallel(ParticleData);
  vNormal=PB::GetVNormal(ParticleData);
  v=sqrt(vParallel*vParallel+vNormal*vNormal);

  //determine the segment of the particle location 
  Segment=FL::FieldLinesAll[iFieldLine].GetSegment(FieldLineCoord); 

  //calcualte mu of the particle
  mu=vParallel/v;

  //calculate B and L
  double B[3],B0[3],B1[3], AbsBDeriv;

  FL::FieldLinesAll[iFieldLine].GetMagneticField(B0, (int)FieldLineCoord);
  FL::FieldLinesAll[iFieldLine].GetMagneticField(B,       FieldLineCoord);
  FL::FieldLinesAll[iFieldLine].GetMagneticField(B1, (int)FieldLineCoord+1-1E-7);
  AbsB   = pow(B[0]*B[0] + B[1]*B[1] + B[2]*B[2], 0.5);

  AbsBDeriv = (pow(B1[0]*B1[0] + B1[1]*B1[1] + B1[2]*B1[2], 0.5) -
    pow(B0[0]*B0[0] + B0[1]*B0[1] + B0[2]*B0[2], 0.5)) /  FL::FieldLinesAll[iFieldLine].GetSegmentLength(FieldLineCoord);

  L=-Vector3D::Length(B)/AbsBDeriv;

  //move the particle along the magnetic field line 
  double FieldLineCoord_init=FieldLineCoord;
  FieldLineCoord=FL::FieldLinesAll[iFieldLine].move(FieldLineCoord,dtTotal*vParallel);

  //get the new value of 'mu'
  double D,dD_dmu;

  double mu_init=mu;
  double time_counter=0.0;
  double dt=dtTotal;
  double dmu=0.0;

  while (time_counter<dtTotal) { 
   if (time_counter+dt>dtTotal) dt=dtTotal-time_counter;

   dmu=0.0;

  switch (_SEP_DIFFUSION_MODEL_) {
  case _DIFFUSION_NONE_:
    //no diffution model is used -> do nothing
    break; 

  case _DIFFUSION_ROUX2004AJ_:
    SEP::Diffusion::Roux2004AJ::GetPitchAngleDiffusionCoefficient(D,dD_dmu,mu,vParallel,vNormal,spec,FieldLineCoord_init,Segment); 

    dmu=sqrt(2.0*D*dt)*Vector3D::Distribution::Normal();
    dmu+=dD_dmu*dt;
    break;
   
  case _DIFFUSION_BOROVIKOV_2019_ARXIV_:
    SEP::Diffusion::Borovokov_2019_ARXIV::GetPitchAngleDiffusionCoefficient(D,dD_dmu,mu,vParallel,vNormal,spec,FieldLineCoord_init,Segment);

    dmu=sqrt(2.0*D*dt)*Vector3D::Distribution::Normal();
    dmu+=dD_dmu*dt;
    break;

  case _DIFFUSION_JOKIPII1966AJ_:
    SEP::Diffusion::Jokopii1966AJ::GetPitchAngleDiffusionCoefficient(D,dD_dmu,mu,vParallel,vNormal,spec,FieldLineCoord_init,Segment);

    dmu=sqrt(2.0*D*dt)*Vector3D::Distribution::Normal();
    dmu+=dD_dmu*dt;
    break;


  default: 
    exit(__LINE__,__FILE__,"Error: the option is unknown");
  }

  dmu+=(1.0-mu*mu)/(2.0*L)*v*dt;

  if (fabs(dmu)>0.1) {
    dt/=2.0;
    time_counter=0.0;
    mu=mu_init;
    continue;
  }

  mu+=dmu;
  time_counter+=dt;

  if (mu<-1.0) mu=-1.0;
  if (mu>1.0) mu=1.0;
}


  //get the segment of the new particle location 
  if ((Segment=FL::FieldLinesAll[iFieldLine].GetSegment(FieldLineCoord))==NULL) {
    //the particle left the computational domain
    int code=_PARTICLE_DELETED_ON_THE_FACE_;
    
    //call the function that process particles that leaved the coputational domain
    switch (code) {
    case _PARTICLE_DELETED_ON_THE_FACE_:
      PIC::ParticleBuffer::DeleteParticle(ptr);
      return _PARTICLE_LEFT_THE_DOMAIN_;

    default:
      exit(__LINE__,__FILE__,"Error: not implemented");
    }
  }

  //set the newvalues of the normal and parallel particle velocities 
  vParallel=mu*v;
  vNormal=sqrt(1.0-mu*mu)*v;

  PB::SetVParallel(vParallel,ParticleData);
  PB::SetVNormal(vNormal,ParticleData);

  //set the new particle coordinate 
  PB::SetFieldLineCoord(FieldLineCoord,ParticleData);

  //attach the particle to the temporaty list
  switch (_PIC_PARTICLE_LIST_ATTACHING_) {
  case  _PIC_PARTICLE_LIST_ATTACHING_NODE_:
    exit(__LINE__,__FILE__,"Error: the function was developed for the case _PIC_PARTICLE_LIST_ATTACHING_==_PIC_PARTICLE_LIST_ATTACHING_FL_SEGMENT_");
    break;
  case _PIC_PARTICLE_LIST_ATTACHING_FL_SEGMENT_:
    PIC::ParticleBuffer::SetNext(Segment->tempFirstParticleIndex,ParticleData);
    PIC::ParticleBuffer::SetPrev(-1,ParticleData);

    if (Segment->tempFirstParticleIndex!=-1) PIC::ParticleBuffer::SetPrev(ptr,Segment->tempFirstParticleIndex);
    Segment->tempFirstParticleIndex=ptr;
    break;
  default:
    exit(__LINE__,__FILE__,"Error: the option is unknown");
  }

  return _PARTICLE_MOTION_FINISHED_;
}

//neutron physics for the Earth's model

#include "pic.h"
#include "Earth.h" 


void Earth::EnergeticParticlesPhysics(long int ptr,long int& FirstParticleCell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  int spec;

  spec=PIC::ParticleBuffer::GetI(ptr);

  switch(spec) {
  case _ELECTRON_SPEC_:
    ElectronPhysics(ptr,FirstParticleCell,node);
    break;
  case _H_PLUS_SPEC_:
    ProtonPhysics(ptr,FirstParticleCell,node);
    break;
  case _NEUTRON_SPEC_:
    NeutronPhysics(ptr,FirstParticleCell,node);
    break;
  }
}
   
void Earth::ElectronPhysics(long int ptr,long int& FirstParticleCell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  double dt=node->block->GetLocalTimeStep(_ELECTRON_SPEC_);
  double *v,n,*x,freq;

  const double CrossSection=1.2E-19; //m^2 J. W. Daiber and H. F. Waldron Scattering Cross Sections of Argon and Atomic Oxygen to Thermal Electrons 
    
  x=PIC::ParticleBuffer::GetX(ptr);
  v=PIC::ParticleBuffer::GetV(ptr);

  n=GetAtmosphereTotalNumberDensity(x); 
  freq=n*CrossSection*Vector3D::Length(v); 

  if (-log(rnd())/freq<dt) {
    //collsion of the electron with the background atmosphere had happed -> the electron is removed from the simulation
    PIC::ParticleBuffer::DeleteParticle(ptr);
  }
  else {
    //the electon survived the time step
    PIC::ParticleBuffer::SetNext(FirstParticleCell,ptr);
    PIC::ParticleBuffer::SetPrev(-1,ptr);

    if (FirstParticleCell!=-1) PIC::ParticleBuffer::SetPrev(ptr,FirstParticleCell);
    FirstParticleCell=ptr;
  }
}   
   
void Earth::ProtonPhysics(long int ptr,long int& FirstParticleCell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  double dt=node->block->GetLocalTimeStep(_ELECTRON_SPEC_);
  double *v,n,*x,freq;
    
  const double CrossSection=1.5E-19;        

  x=PIC::ParticleBuffer::GetX(ptr);
  v=PIC::ParticleBuffer::GetV(ptr);

  n=GetAtmosphereTotalNumberDensity(x);
  freq=n*CrossSection*Vector3D::Length(v);

  if (-log(rnd())/freq<dt) {
    //collsion of the electron with the background atmosphere had happed -> the electron is removed from the simulation
    PIC::ParticleBuffer::DeleteParticle(ptr);
  }
  else {
    //the electon survived the time step
    PIC::ParticleBuffer::SetNext(FirstParticleCell,ptr);
    PIC::ParticleBuffer::SetPrev(-1,ptr);
    
    if (FirstParticleCell!=-1) PIC::ParticleBuffer::SetPrev(ptr,FirstParticleCell);
    FirstParticleCell=ptr;
  }
}

void Earth::NeutronPhysics(long int ptr,long int& FirstParticleCell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  double dt=node->block->GetLocalTimeStep(_ELECTRON_SPEC_);
  double *v,n,*x,freq_total,freq_collision;

  const double CrossSection=1.5E-19; 
  const double BetaDecayLifetime=879.6;  
   
  //two parallel processes can occur: 1. interaction with the atmosphere, and 2. beta-decay
  const double freq_decay=1.0/BetaDecayLifetime;
  
  x=PIC::ParticleBuffer::GetX(ptr);
  v=PIC::ParticleBuffer::GetV(ptr);

  n=GetAtmosphereTotalNumberDensity(x);
  freq_collision=n*CrossSection*Vector3D::Length(v);

  freq_total=freq_collision+freq_decay;
  if (-log(rnd())/freq_total<dt) { 
    //the neutron has decayed through one of the two channels 
   
    if (rnd()<freq_decay/freq_total) {
      //the neutron decayed via bete-decay
      long int ptr_electron,ptr_proton;

      ptr_electron=PIC::ParticleBuffer::GetNewParticle();
      ptr_proton=PIC::ParticleBuffer::GetNewParticle();

      PIC::ParticleBuffer::SetNext(FirstParticleCell,ptr_electron);
      PIC::ParticleBuffer::SetPrev(-1,ptr_electron);
      if (FirstParticleCell!=-1) PIC::ParticleBuffer::SetPrev(ptr_electron,FirstParticleCell);
      FirstParticleCell=ptr_electron;

      PIC::ParticleBuffer::SetNext(FirstParticleCell,ptr_proton);
      PIC::ParticleBuffer::SetPrev(-1,ptr_proton);
      if (FirstParticleCell!=-1) PIC::ParticleBuffer::SetPrev(ptr_proton,FirstParticleCell);
      FirstParticleCell=ptr_proton;

      PIC::ParticleBuffer::SetI(_ELECTRON_SPEC_,ptr_electron); 
      PIC::ParticleBuffer::SetI(_H_PLUS_SPEC_,ptr_proton);

      PIC::ParticleBuffer::SetX(x,ptr_electron);
      PIC::ParticleBuffer::SetX(x,ptr_proton);

      //determine velocities of the new particles 
      double e_total=5.5*MeV2J;
      double e_electron=0.78*MeV2J;
      double e_proton=e_total-e_electron; 
      double v_new[3],speed;


      speed=Relativistic::E2Speed(e_electron,PIC::MolecularData::GetMass(_ELECTRON_SPEC_));
      Vector3D::Distribution::Uniform(v_new,speed);

      for (int idim=0;idim<3;idim++) v_new[idim]+=v[idim];
      PIC::ParticleBuffer::SetV(v_new,ptr_electron);  

      speed=Relativistic::E2Speed(e_proton,PIC::MolecularData::GetMass(_H_PLUS_SPEC_));
      Vector3D::Distribution::Uniform(v_new,speed);

      for (int idim=0;idim<3;idim++) v_new[idim]+=v[idim];
      PIC::ParticleBuffer::SetV(v_new,ptr_proton);
    }

    //delete the original neutron
    PIC::ParticleBuffer::DeleteParticle(ptr);
  }
  else {
    //the neutron survived the time step    
    PIC::ParticleBuffer::SetNext(FirstParticleCell,ptr);
    PIC::ParticleBuffer::SetPrev(-1,ptr);

    if (FirstParticleCell!=-1) PIC::ParticleBuffer::SetPrev(ptr,FirstParticleCell);
    FirstParticleCell=ptr;
  }
}











    

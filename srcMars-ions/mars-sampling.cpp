//$Id$
//sampling function for the Mars-ion model

#include "mars-ions.h"
#include "pic.h"


/*
 * mars-sampling.cpp
 *
 *  Created on: May 4, 2017
 *      Author: vtenishe
 */



void MarsIon::Sampling::PrintManager(int nDataSet) {
  char fname[300];
  int iShell,spec;

  for (iShell=0;iShell<SphericalShells::nSphericalShells;iShell++) {
    for (spec=0;spec<PIC::nTotalSpecies;spec++) {
      sprintf(fname,"%s/mars.shell=%i.spec=%i.out=%i.dat",PIC::OutputDataFileDirectory,iShell,spec,nDataSet);
      SphericalShells::SamplingSphericlaShell[iShell].PrintSurfaceData(fname,spec,true);
    }

    SphericalShells::SamplingSphericlaShell[iShell].flush();
  }
}

//empty fucntion. needed for compartibility with the core
void MarsIon::Sampling::SamplingManager() {}


void MarsIon::Sampling::SphericalShells::Output::PrintTitle(FILE* fout) {
  fprintf(fout,"TITLE=\"Sampled particle ion output flux\"");
}

void MarsIon::Sampling::SphericalShells::Output::PrintVariableList(FILE* fout) {
  fprintf(fout,", \"Total Flux Up\", \"Total Flux Down\"");

  for (int iEnergyLevel=0;iEnergyLevel<MarsIon::Sampling::SphericalShells::nSampledLevels;iEnergyLevel++) {
    double emin=MarsIon::Sampling::SphericalShells::minSampledEnergy*pow(10.0,iEnergyLevel*MarsIon::Sampling::SphericalShells::dLogEnergy);
    double emax=emin*pow(10.0,MarsIon::Sampling::SphericalShells::dLogEnergy);

    emin*=J2eV;
    emax*=J2eV;

    fprintf(fout,", \"Ion Flux Up (%e<E<%e eV)\", \" Ion Flux Down (%e<E<%e) eV\"",emin,emax, emin,emax);
  }
}

void MarsIon::Sampling::SphericalShells::Output::PrintDataStateVector(FILE* fout,long int nZenithPoint,long int nAzimuthalPoint,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength,cInternalSphericalData *Sphere,
    int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads) {

  double Flux=0.0,Rigidity=-1.0,ParticleFluxUp=0.0,ParticleFluxDown=0.0;
  int i,el;
  double elArea,TotalStencilArea=0.0;

  //energy distribution of the particle fluency
  static double *ParticleFluenceUp=NULL;
  static double *ParticleFluenceDown=NULL;

  if (ParticleFluenceUp==NULL) {
    ParticleFluenceUp=new double [MarsIon::Sampling::SphericalShells::nSampledLevels];
    ParticleFluenceDown=new double [MarsIon::Sampling::SphericalShells::nSampledLevels];
  }

  for (i=0;i<MarsIon::Sampling::SphericalShells::nSampledLevels;i++) ParticleFluenceDown[i]=0.0,ParticleFluenceUp[i]=0.0;

  //get averaged on a given processor
  for (i=0;i<SurfaceElementsInterpolationListLength;i++) {
    el=SurfaceElementsInterpolationList[i];

    elArea=Sphere->GetSurfaceElementArea(el);
    TotalStencilArea+=elArea;

    Flux+=Sphere->Flux[spec][el];
    ParticleFluxUp+=Sphere->ParticleFluxUp[spec][el];
    ParticleFluxDown+=Sphere->ParticleFluxDown[spec][el];

    for (int iLevel=0;iLevel<MarsIon::Sampling::SphericalShells::nSampledLevels;iLevel++) {
      ParticleFluenceDown[iLevel]+=Sphere->ParticleFluencyDown[spec][el][iLevel];
      ParticleFluenceUp[iLevel]+=Sphere->ParticleFluencyUp[spec][el][iLevel];
    }

    if (Sphere->minRigidity[spec][el]>0.0) {
      if ((Rigidity<0.0) || (Sphere->minRigidity[spec][el]<Rigidity)) Rigidity=Sphere->minRigidity[spec][el];
    }
  }

  double SamplingTime;

  #if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  SamplingTime=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec]*PIC::LastSampleLength;
  #else
  exit(__LINE__,__FILE__,"Error: not implemented for a local times step model");
  #endif

  Flux/=TotalStencilArea*SamplingTime;

  ParticleFluxUp/=TotalStencilArea*SamplingTime;
  ParticleFluxDown/=TotalStencilArea*SamplingTime;
  for (i=0;i<MarsIon::Sampling::SphericalShells::nSampledLevels;i++) ParticleFluenceDown[i]/=TotalStencilArea*SamplingTime,ParticleFluenceUp[i]/=TotalStencilArea*SamplingTime;

  //send data to the root processor, and print the data
  if (PIC::ThisThread==0) {
    double t;
    int thread;

    for (thread=1;thread<PIC::nTotalThreads;thread++) {
      Flux+=pipe->recv<double>(thread);

      ParticleFluxUp+=pipe->recv<double>(thread);
      ParticleFluxDown+=pipe->recv<double>(thread);

      for (i=0;i<MarsIon::Sampling::SphericalShells::nSampledLevels;i++) {
        ParticleFluenceDown[i]+=pipe->recv<double>(thread);
        ParticleFluenceUp[i]+=pipe->recv<double>(thread);
      }

      t=pipe->recv<double>(thread);
      if ((t>0.0) && ((t<Rigidity)||(Rigidity<0.0)) ) Rigidity=t;
    }

    fprintf(fout," %e  %e  ",ParticleFluxUp,ParticleFluxDown);
    for (i=0;i<MarsIon::Sampling::SphericalShells::nSampledLevels;i++) fprintf(fout," %e  %e ", ParticleFluenceUp[i],ParticleFluenceDown[i]);
  }
  else {
    pipe->send(Flux);

    pipe->send(ParticleFluxUp);
    pipe->send(ParticleFluxDown);

    for (i=0;i<MarsIon::Sampling::SphericalShells::nSampledLevels;i++) {
      pipe->send(ParticleFluenceDown[i]);
      pipe->send(ParticleFluenceUp[i]);
    }

    pipe->send(Rigidity);
  }
}

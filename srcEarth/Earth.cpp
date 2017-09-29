/*
 * Europa.cpp
 *
 *  Created on: Feb 13, 2012
 *      Author: vtenishe
 */

//$Id$



#include "pic.h"
#include "Earth.h"



char Earth::Mesh::sign[_MAX_STRING_LENGTH_PIC_]="";
char Exosphere::ObjectName[_MAX_STRING_LENGTH_PIC_]="";
char Exosphere::SO_FRAME[_MAX_STRING_LENGTH_PIC_]="";

//composition of the GCRs
cCompositionGroupTable *Earth::CompositionGroupTable=NULL;
int *Earth::CompositionGroupTableIndex=NULL;
int Earth::nCompositionGroups;

//sampling of the particle flux at different altitudes
cInternalSphericalData Earth::Sampling::SamplingSphericlaShell[Earth::Sampling::nSphericalShells];

//sampling of the paericle fluency
int Earth::Sampling::Fluency::nSampledLevels=10;
double Earth::Sampling::Fluency::dLogEnergy=0.0;
double Earth::Sampling::Fluency::minSampledEnergy=1.0*MeV2J;
double Earth::Sampling::Fluency::maxSampledEnergy=100.0*MeV2J;

int Earth_Sampling_Fluency_nSampledLevels=10;

void Earth::Sampling::PrintVariableList(FILE* fout) {
  fprintf(fout,", \"TotalFlux\", \"min Rigidity\", \"FluxUp\", \"FluxDown\"");

  for (int iEnergyLevel=0;iEnergyLevel<Earth::Sampling::Fluency::nSampledLevels;iEnergyLevel++) {
    double emin=Fluency::minSampledEnergy*pow(10.0,iEnergyLevel*Fluency::dLogEnergy);
    double emax=emin*pow(10.0,Fluency::dLogEnergy);

    emin*=J2eV;
    emax*=J2eV;

    fprintf(fout,", \"TotalFluency(%e<E<%e)\", \"FluencyUp(%e<E<%e)\", \"FluencyDown(%e<E<%e)\"",emin,emax, emin,emax, emin,emax);
  }
}

void Earth::Sampling::PrintTitle(FILE* fout) {
  fprintf(fout,"TITLE=\"Sampled particle flux and max rigidity\"");
}

void Earth::Sampling::PrintDataStateVector(FILE* fout,long int nZenithPoint,long int nAzimuthalPoint,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength,cInternalSphericalData *Sphere,
    int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads) {

  double Flux=0.0,Rigidity=-1.0,ParticleFluxUp=0.0,ParticleFluxDown=0.0;
  int i,el;
  double elArea,TotalStencilArea=0.0;

  //energy distribution of the particle fluency
  static double *ParticleFluenceUp=NULL;
  static double *ParticleFluenceDown=NULL;

  if (ParticleFluenceUp==NULL) {
    ParticleFluenceUp=new double [Fluency::nSampledLevels];
    ParticleFluenceDown=new double [Fluency::nSampledLevels];
  }

  for (i=0;i<Fluency::nSampledLevels;i++) ParticleFluenceDown[i]=0.0,ParticleFluenceUp[i]=0.0;

  //get averaged on a given processor
  for (i=0;i<SurfaceElementsInterpolationListLength;i++) {
    el=SurfaceElementsInterpolationList[i];

    elArea=Sphere->GetSurfaceElementArea(el);
    TotalStencilArea+=elArea;

    Flux+=Sphere->Flux[spec][el];
    ParticleFluxUp+=Sphere->ParticleFluxUp[spec][el];
    ParticleFluxDown+=Sphere->ParticleFluxDown[spec][el];

    for (int iLevel=0;iLevel<Fluency::nSampledLevels;iLevel++) {
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
  for (i=0;i<Fluency::nSampledLevels;i++) ParticleFluenceDown[i]/=TotalStencilArea*SamplingTime,ParticleFluenceUp[i]/=TotalStencilArea*SamplingTime;

  //send data to the root processor, and print the data
  if (PIC::ThisThread==0) {
    double t;
    int thread;

    for (thread=1;thread<PIC::nTotalThreads;thread++) {
      Flux+=pipe->recv<double>(thread);

      ParticleFluxUp+=pipe->recv<double>(thread);
      ParticleFluxDown+=pipe->recv<double>(thread);

      for (i=0;i<Fluency::nSampledLevels;i++) {
        ParticleFluenceDown[i]+=pipe->recv<double>(thread);
        ParticleFluenceUp[i]+=pipe->recv<double>(thread);
      }

      t=pipe->recv<double>(thread);
      if ((t>0.0) && ((t<Rigidity)||(Rigidity<0.0)) ) Rigidity=t;
    }

    fprintf(fout," %e  %e  %e  %e ", ParticleFluxUp+ParticleFluxDown,((Rigidity>0.0) ? Rigidity : 0.0),ParticleFluxUp,ParticleFluxDown);
    for (i=0;i<Fluency::nSampledLevels;i++) fprintf(fout," %e  %e  %e ", ParticleFluenceUp[i]+ParticleFluenceDown[i],ParticleFluenceUp[i],ParticleFluenceDown[i]);
  }
  else {
    pipe->send(Flux);

    pipe->send(ParticleFluxUp);
    pipe->send(ParticleFluxDown);

    for (i=0;i<Fluency::nSampledLevels;i++) {
      pipe->send(ParticleFluenceDown[i]);
      pipe->send(ParticleFluenceUp[i]);
    }

    pipe->send(Rigidity);
  }
}

void Earth::Sampling::PrintManager(int nDataSet) {
  char fname[300];
  int iShell,spec;

  for (iShell=0;iShell<nSphericalShells;iShell++) {
    for (spec=0;spec<PIC::nTotalSpecies;spec++) {
      sprintf(fname,"%s/earth.shell=%i.spec=%i.out=%i.dat",PIC::OutputDataFileDirectory,iShell,spec,nDataSet);
      SamplingSphericlaShell[iShell].PrintSurfaceData(fname,spec,true);
    }

    SamplingSphericlaShell[iShell].flush();
  }
}

//empty fucntion. needed for compartibility with the core
void Earth::Sampling::SamplingManager() {}

//particle mover: call relativistic Boris, ans sample particle flux
int Earth::ParticleMover(long int ptr,double dtTotal,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* startNode) {
  double xInit[3],xFinal[3];
  int res,iShell;

  PIC::ParticleBuffer::GetX(xInit,ptr);



  switch (PIC::ParticleBuffer::GetI(ptr)) {
  case _ELECTRON_SPEC_:
//    res=PIC::Mover::GuidingCenter::Mover_SecondOrder(ptr,dtTotal,startNode);
res=PIC::Mover::Relativistic::Boris(ptr,dtTotal,startNode);

    break;
  default:
    res=PIC::Mover::Relativistic::Boris(ptr,dtTotal,startNode);
  }

  if ((Sampling::SamplingMode==true)&&(res==_PARTICLE_MOTION_FINISHED_)) {
    //trajectory integration was successfully completed
    double r0,r1;

    PIC::ParticleBuffer::GetX(xFinal,ptr);
    r0=sqrt(xInit[0]*xInit[0]+xInit[1]*xInit[1]+xInit[2]*xInit[2]);
    r1=sqrt(xFinal[0]*xFinal[0]+xFinal[1]*xFinal[1]+xFinal[2]*xFinal[2]);

    //Important: now r0<=r1
    //search for a shell that was intersected by the trajectory segment
    for (iShell=0;iShell<Sampling::nSphericalShells;iShell++) if ((r0-Sampling::SampleSphereRadii[iShell])*(r1-Sampling::SampleSphereRadii[iShell])<=0.0) {
      //the trajectory segment has intersected a sampling shell
      //get the point of intersection
      double l[3],t=-1.0,a=0.0,b=0.0,c=0.0,d4,t1,t2;
      int idim;

      for (idim=0;idim<3;idim++) {
        a+=pow(xFinal[idim]-xInit[idim],2);
        b+=xInit[idim]*(xFinal[idim]-xInit[idim]);
        c+=pow(xInit[idim],2);
      }

      c-=pow(Sampling::SampleSphereRadii[iShell],2);
      d4=sqrt(b*b-a*c);

      if ((isfinite(d4)==true)&&(a!=0.0)) {
        //the determinent is real
        t1=(-b+d4)/a;
        t2=(-b-d4)/a;

        if ((0.0<t1)&&(t1<1.0)) t=t1;
        if ((0.0<t2)&&(t2<1.0)) if ((t<0.0)||(t2<t)) t=t2;

        if (t>0.0) {
          //the point of intersection is found. Sample flux and rigidity at that location
          double xIntersection[3];
          long int iSurfaceElementNumber,iZenithElement,iAzimuthalElement;

          for (idim=0;idim<3;idim++) xIntersection[idim]=xInit[idim]+t*(xFinal[idim]-xInit[idim]);

          Sampling::SamplingSphericlaShell[iShell].GetSurfaceElementProjectionIndex(xIntersection,iZenithElement,iAzimuthalElement);
          iSurfaceElementNumber=Sampling::SamplingSphericlaShell[iShell].GetLocalSurfaceElementNumber(iZenithElement,iAzimuthalElement);

          //increment sampling counters
          double ParticleWeight,Rigidity;
          int spec;


          spec=PIC::ParticleBuffer::GetI(ptr);
          ParticleWeight=startNode->block->GetLocalParticleWeight(spec)*PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ptr);


          #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
          double *t;
          t=Sampling::SamplingSphericlaShell[iShell].Flux[spec]+iSurfaceElementNumber;

          #pragma omp atomic
          *t+=ParticleWeight;
          #else
          Sampling::SamplingSphericlaShell[iShell].Flux[spec][iSurfaceElementNumber]+=ParticleWeight;
          #endif

          //Fluence Sampling Energy level
          int iEnergyLevel=-1;
          double Energy,Speed,v[3];

          PIC::ParticleBuffer::GetV(v,ptr);
          Speed=Vector3D::Length(v);
          Energy=Relativistic::Speed2E(Speed,PIC::MolecularData::GetMass(spec));

          if ((Sampling::Fluency::minSampledEnergy<Energy)&&(Energy<Sampling::Fluency::maxSampledEnergy)) {
            iEnergyLevel=(int)(log10(Energy/Sampling::Fluency::minSampledEnergy)/Sampling::Fluency::dLogEnergy);

            if ((iEnergyLevel<0)||(iEnergyLevel>=Sampling::Fluency::nSampledLevels)) iEnergyLevel=-1;
          }
          else iEnergyLevel=-1;

          if (r0>Sampling::SampleSphereRadii[iShell]) {
            //the particle moves toward hte Earth
            #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
            t=Sampling::SamplingSphericlaShell[iShell].ParticleFluxDown[spec]+iSurfaceElementNumber;

            #pragma omp atomic
            *t+=ParticleWeight;

            if (iEnergyLevel!=-1) {
              t=Sampling::SamplingSphericlaShell[iShell].ParticleFluencyDown[spec][iSurfaceElementNumber]+iEnergyLevel;

              #pragma omp atomic
              *t+=ParticleWeight;

              //Sampling::SamplingSphericlaShell[iShell].ParticleFluencyDown[spec][iSurfaceElementNumber][iEnergyLevel]+=ParticleWeight;
            }


            #else //_COMPILATION_MODE__HYBRID_
            Sampling::SamplingSphericlaShell[iShell].ParticleFluxDown[spec][iSurfaceElementNumber]+=ParticleWeight;
            if (iEnergyLevel!=-1) Sampling::SamplingSphericlaShell[iShell].ParticleFluencyDown[spec][iSurfaceElementNumber][iEnergyLevel]+=ParticleWeight;
            #endif //_COMPILATION_MODE__HYBRID_

          }
          else {
            //the particle moves outward the Earth
            #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
            t=Sampling::SamplingSphericlaShell[iShell].ParticleFluxUp[spec]+iSurfaceElementNumber;

            #pragma omp atomic
            *t+=ParticleWeight;

            if (iEnergyLevel!=-1) {
              t=Sampling::SamplingSphericlaShell[iShell].ParticleFluencyUp[spec][iSurfaceElementNumber]+iEnergyLevel;

              #pragma omp atomic
              *t+=ParticleWeight;

              //Sampling::SamplingSphericlaShell[i].ParticleFluencyUp[spec][iSurfaceElementNumber][iEnergyLevel]+=ParticleWeight;
            }

            #else //_COMPILATION_MODE__HYBRID_
            Sampling::SamplingSphericlaShell[iShell].ParticleFluxUp[spec][iSurfaceElementNumber]+=ParticleWeight;
            if (iEnergyLevel!=-1) Sampling::SamplingSphericlaShell[iShell].ParticleFluencyUp[spec][iSurfaceElementNumber][iEnergyLevel]+=ParticleWeight;
            #endif //_COMPILATION_MODE__HYBRID_
          }

          //calculate particle rigidity
          double ElectricCharge;

          ElectricCharge=fabs(PIC::MolecularData::GetElectricCharge(spec));

          if (ElectricCharge>0.0) {
            Rigidity=Relativistic::Speed2Momentum(Speed,PIC::MolecularData::GetMass(spec))/PIC::MolecularData::GetElectricCharge(spec);

            #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
            #pragma omp critical
            #endif
            {
              if ((Sampling::SamplingSphericlaShell[iShell].minRigidity[spec][iSurfaceElementNumber]<0.0) || (Rigidity<Sampling::SamplingSphericlaShell[iShell].minRigidity[spec][iSurfaceElementNumber])) {
                Sampling::SamplingSphericlaShell[iShell].minRigidity[spec][iSurfaceElementNumber]=Rigidity;
              }
            }

          }
        }
      }

    }
  }
}


void Earth::Sampling::Init() {
  int iShell;


  if (Earth::Sampling::SamplingMode==true) {
    //set the user-function for output of the data files in the core
    PIC::Sampling::ExternalSamplingLocalVariables::RegisterSamplingRoutine(SamplingManager,PrintManager);

    //set parameters of the fluency sampling
    Earth_Sampling_Fluency_nSampledLevels=Fluency::nSampledLevels;
    Fluency::dLogEnergy=log10(Fluency::maxSampledEnergy/Fluency::minSampledEnergy)/Fluency::nSampledLevels;

    //init the shells
     for (iShell=0;iShell<nSphericalShells;iShell++) {
      double sx0[3]={0.0,0.0,0.0};

      SamplingSphericlaShell[iShell].PrintDataStateVector=PrintDataStateVector;
      SamplingSphericlaShell[iShell].PrintTitle=PrintTitle;
      SamplingSphericlaShell[iShell].PrintVariableList=PrintVariableList;

      SamplingSphericlaShell[iShell].SetSphereGeometricalParameters(sx0,SampleSphereRadii[iShell]);
      SamplingSphericlaShell[iShell].Allocate<cInternalSphericalData>(PIC::nTotalSpecies,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber,_EXOSPHERE__SOURCE_MAX_ID_VALUE_,SamplingSphericlaShell+iShell);
    }
  }
}

//manager of the data recovering
void Earth::DataRecoveryManager(list<pair<string,list<int> > >& SampledDataRecoveryTable ,int MinOutputFileNumber,int MaxOutputFilenumber) {
  int iOutputFile;
  char fname[_MAX_STRING_LENGTH_PIC_];
  pair<string,list<int> > DataRecoveryTableElement;

  DataRecoveryTableElement.second.push_back(0);

  for (iOutputFile=MinOutputFileNumber;iOutputFile<=MaxOutputFilenumber;iOutputFile++) {
    sprintf(fname,"pic.SamplingDataRestart.out=%i.dat",iOutputFile);

    DataRecoveryTableElement.first=fname;
    SampledDataRecoveryTable.push_back(DataRecoveryTableElement);
  }
}

//calcualte the true anomaly angle
double Exosphere::OrbitalMotion::GetTAA(SpiceDouble et) {
  return 0.0;
}

int Exosphere::ColumnIntegral::GetVariableList(char *vlist) {
  int spec,nVariables=0;
  return nVariables;
}

void Exosphere::ColumnIntegral::ProcessColumnIntegrationVector(double *res,int resLength) {
//do nothing
}

double Exosphere::SurfaceInteraction::StickingProbability(int spec,double& ReemissionParticleFraction,double Temp) {
  double res=1.0;
  ReemissionParticleFraction=0.0;

  if (spec==_O2_SPEC_) res=0.0;

  return res;
}

double Exosphere::GetSurfaceTemeprature(double cosSubsolarAngle,double *x) {
  return 300.0;
}

void Exosphere::ColumnIntegral::CoulumnDensityIntegrant(double *res,int resLength,double* x,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  //do nothing
}

//particle/sphere interactions
int Earth::BC::ParticleSphereInteraction(int spec,long int ptr,double *x,double *v, double &dtTotal, void *NodeDataPonter,void *SphereDataPointer)  {
  return _PARTICLE_DELETED_ON_THE_FACE_;
}

//the total injection rate from the Earth
double Earth::BC::sphereInjectionRate(int spec,void *SphereDataPointer) {
  double res=0.0;
  return res;
}

//init the Earth magnetosphere model
void Earth::Init() {
  //init the composition gourp tables
  //!!!!!!!!!!!!! For now only hydrogen is considered !!!!!!!!!!!!!!!!!

  //composition of the GCRs
  nCompositionGroups=1;
  CompositionGroupTable=new cCompositionGroupTable[nCompositionGroups];
  CompositionGroupTableIndex=new int[PIC::nTotalSpecies];

  for (int spec=0;spec<PIC::nTotalSpecies;spec++) CompositionGroupTableIndex[spec]=0; //all simulated model species are hydrogen
  CompositionGroupTable[0].FistGroupSpeciesNumber=0;
  CompositionGroupTable[0].nModelSpeciesGroup=PIC::nTotalSpecies;

  CompositionGroupTable[0].minVelocity=Relativistic::E2Speed(Earth::BoundingBoxInjection::minEnergy,PIC::MolecularData::GetMass(0));
  CompositionGroupTable[0].maxVelocity=Relativistic::E2Speed(Earth::BoundingBoxInjection::maxEnergy,PIC::MolecularData::GetMass(0));

  CompositionGroupTable[0].GroupVelocityStep=(CompositionGroupTable[0].maxVelocity-CompositionGroupTable[0].minVelocity)/CompositionGroupTable[0].nModelSpeciesGroup;

  if (PIC::ThisThread==0) {
    cout << "$PREFIX: Composition Group Velocity and Energy Characteristics:\nspec\tmin Velocity [m/s]\tmax Velovity[m/s]\t min Energy[eV]\tmax Energy[eV]" << endl;

    for (int s=0;s<PIC::nTotalSpecies;s++) {
      double minV,maxV,minE,maxE,mass;

      mass=PIC::MolecularData::GetMass(s);

      minV=Earth::CompositionGroupTable[0].GetMinVelocity(s);
      maxV=Earth::CompositionGroupTable[0].GetMaxVelocity(s);

      //convert velocity into energy and distribute energy of a new particles
      minE=Relativistic::Speed2E(minV,mass);
      maxE=Relativistic::Speed2E(maxV,mass);

      cout << s << "\t" << minV << "\t" << maxV << "\t" << minE*J2eV << "\t" <<  maxE*J2eV << endl;
    }
  }

  //init source models of SEP and GCR
  if (_PIC_EARTH_SEP__MODE_==_PIC_MODE_ON_) BoundingBoxInjection::SEP::Init();
  if (_PIC_EARTH_GCR__MODE_==_PIC_MODE_ON_) BoundingBoxInjection::GCR::Init();
  if (_PIC_EARTH_ELECTRON__MODE_==_PIC_MODE_ON_) BoundingBoxInjection::Electrons::Init();
}



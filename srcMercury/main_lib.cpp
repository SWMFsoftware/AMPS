

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <list>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <iostream>
#include <fstream>
#include <time.h>

#include <sys/time.h>
#include <sys/resource.h>

//$Id$


//#include "vt_user.h"
//#include <VT.h>

//the particle class
#include "pic.h"
#include "constants.h"
#include "Mercury.h"


//the modeling case
#define _MERCURY_MODEL_MODE__NEAR_PLANET_    0
#define _MERCURY_MODEL_MODE__TAIL_MODEL_     1

#define _MERCURY_MODEL_MODE_ _MERCURY_MODEL_MODE__TAIL_MODEL_

//Check the mesh consistency
//#define _CHECK_MESH_CONSISTENCY_
//#define _ICES_CREATE_COORDINATE_LIST_
#define _ICES_LOAD_DATA_


//defiend the modes for included models
#define _MERCURY_MODE_ON_    0
#define _MERCURY_MODE_OFF_   1

#define _MERCURY_IMPACT_VAPORIZATION_MODE_ _MERCURY_MODE_ON_
#define _MERCURY_PSD_MODE_ _MERCURY_MODE_ON_
#define _MERCURY_THERMAL_DESORPTION_MODE_ _MERCURY_MODE_ON_

//sample particles that enter FIPS
#define _MERCURY_FIPS_SAMPLING_  _MERCURY_MODE_ON_

//the parameters of the domain and the sphere

const double DebugRunMultiplier=4.0;


const double rSphere=_RADIUS_(_TARGET_);

//const double xMaxDomain=400.0; //modeling of the tail
const double xMaxDomain=10; //modeling the vicinity of the planet
const double yMaxDomain=15; //the minimum size of the domain in the direction perpendicular to the direction to the sun


const double dxMinGlobal=DebugRunMultiplier*2.0,dxMaxGlobal=DebugRunMultiplier*10.0;
const double dxMinSphere=DebugRunMultiplier*4.0*1.0/100,dxMaxSphere=DebugRunMultiplier*2.0/10.0;



//the species
int NA=0;
int NAPLUS=1;


//the total acceleration of the particles
#include "Na.h"

//photoionization lifetime of sodium
double sodiumPhotoionizationLifeTime(double *x,int spec,int ptr,bool &PhotolyticReactionAllowedFlag) {

  static const double LifeTime=3600.0*5.8;

  double res,r2=x[1]*x[1]+x[2]*x[2];

  if ((r2>rSphere*rSphere)||(x[0]<0.0)) {
    res=LifeTime,PhotolyticReactionAllowedFlag=true;
  }
  else {
    res=-1.0,PhotolyticReactionAllowedFlag=false;
  }

  return res;
}

int sodiumPhotoionizationReactionProcessor(double *xInit,double *xFinal,int ptr,int &spec,PIC::ParticleBuffer::byte *ParticleData) {
  spec=NAPLUS;

  PIC::ParticleBuffer::SetI(spec,ParticleData);
//  return _PHOTOLYTIC_REACTIONS_PARTICLE_SPECIE_CHANGED_;

  return _PHOTOLYTIC_REACTIONS_PARTICLE_REMOVED_;
}


//the codes for the soruces processes
#define _ALL_SOURCE_PROCESSES_                -1
#define _IMPACT_VAPORIZATION_SOURCE_PROCESS_   0
#define _PSD_SOURCE_PROCESS_                   1

//sodium surface production
double sodiumTotalProductionRate(int SourceProcessCode=-1) {
  double res=0.0;

#if _MERCURY_IMPACT_VAPORIZATION_MODE_ == _MERCURY_MODE_ON_
  if ((SourceProcessCode==-1)||(SourceProcessCode==_IMPACT_VAPORIZATION_SOURCE_PROCESS_)) { //the total impact vaporization flux
    res+=2.6E23;
  }
#endif

#if _MERCURY_PSD_MODE_ == _MERCURY_MODE_ON_
  if ((SourceProcessCode==-1)||(SourceProcessCode==_PSD_SOURCE_PROCESS_)) { //the photon stimulated desorption flux
    res+=3.6E24;
  }
#endif

  return res;
}


bool sodiumDistributeParticleParameters(double *x,double *v, cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode,double &ParticleWeightCorrection) {

  int const nIntervalsLoockupTable=10000;

  //determine throught which process the particle will be generated
  static bool initflag=false;


  static double injectionProbabilityIV_LimitMin=-1.0,injectionProbabilityIV_LimitMax=-1.0;
  static double injectionProbabilityPSD_LimitMin=-1.0,injectionProbabilityPSD_LimitMax=-1.0;

  static double EnergyDistributionPSDTable[nIntervalsLoockupTable];

  ParticleWeightCorrection=1.0;


  if (initflag==false) {
    initflag=true;

    //init the probability limits
    double upperProbabilityUpperLimit=0.0,totalPtoductionRate=sodiumTotalProductionRate();

#if _MERCURY_IMPACT_VAPORIZATION_MODE_ == _MERCURY_MODE_ON_
    injectionProbabilityIV_LimitMin=upperProbabilityUpperLimit;
    upperProbabilityUpperLimit+=sodiumTotalProductionRate(_IMPACT_VAPORIZATION_SOURCE_PROCESS_)/totalPtoductionRate;
    injectionProbabilityIV_LimitMax=upperProbabilityUpperLimit;
#endif

#if _MERCURY_PSD_MODE_ == _MERCURY_MODE_ON_
    injectionProbabilityPSD_LimitMin=upperProbabilityUpperLimit;
    upperProbabilityUpperLimit+=sodiumTotalProductionRate(_PSD_SOURCE_PROCESS_)/totalPtoductionRate;
    injectionProbabilityPSD_LimitMax=upperProbabilityUpperLimit;


    //generte the loock-up table for the generation of the energy of injected particle
    double emin=0.0,emax=0.5*PIC::MolecularData::GetMass(NA)*pow(2000.0,2);

    const int iIntegratedIntervals=10.0*nIntervalsLoockupTable;
    double e,de=(emax-emin)/iIntegratedIntervals,A=0.0,dA;
    int i,iTable;

    const double beta=0.7,U=0.052*ElectronCharge;

    for (i=0;i<iIntegratedIntervals;i++) {
      e=emin+((double)i+0.5)*de;
      A+=e/pow(e+U,2.0+beta);
    }

    dA=A/nIntervalsLoockupTable;
    A=0.0;
    EnergyDistributionPSDTable[0]=emin;

    for (i=0,iTable=0;i<iIntegratedIntervals;i++) {
      e=emin+((double)i+0.5)*de;
      A+=e/pow(e+U,2.0+beta);

      if (A>dA) {
        A-=dA;
        ++iTable;

        if (iTable==nIntervalsLoockupTable) break;

        EnergyDistributionPSDTable[iTable]=e;
      }
    }

#endif
  }


  //when the process is determines: defined the position and the velocty of the newly generated particle
  double r=rnd();

  if  ((injectionProbabilityIV_LimitMin<=r)&&(r<=injectionProbabilityIV_LimitMax)) { //injection though the impact vaporization process
    static const double InjectionTemperature=2000.0;  //the temeprature of the injected sodium atoms

    int idim;
    double ExternalNormal[3];
    double vbulk[3]={0.0,0.0,0.0};
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;

    for (r=0.0,idim=0;idim<DIM;idim++) {
      ExternalNormal[idim]=sqrt(-2.0*log(rnd()))*cos(PiTimes2*rnd());
      r+=ExternalNormal[idim]*ExternalNormal[idim];
    }


//================   DEBUG ===================
    /*
    ExternalNormal[0]=1,ExternalNormal[1]=0.0,ExternalNormal[2]=0.0;
    r=1;
    */
//================  END DEBUG =================

    r=-sqrt(r);

    for (idim=0;idim<DIM;idim++) {
      ExternalNormal[idim]/=r;
      x[idim]=-rSphere*ExternalNormal[idim];
    }

    node=PIC::Mesh::mesh->findTreeNode(x,startNode);
    if (node->Thread!=PIC::Mesh::mesh->ThisThread) return false;
    startNode=node;

    //generate the particle velocity
    static const double MinParticleWeightCorrectionFactor=1.0E-5;
    bool flag;

    do {
      ParticleWeightCorrection=PIC::Distribution::InjectMaxwellianDistribution(v,vbulk,InjectionTemperature,ExternalNormal,NA,_PIC_DISTRIBUTION_WEIGHT_CORRECTION_MODE__INDIVIDUAL_PARTICLE_WEIGHT_);

      if (ParticleWeightCorrection>MinParticleWeightCorrectionFactor) flag=true;
      else {
        double p; // <= the probability to accept the particle

        p=ParticleWeightCorrection/MinParticleWeightCorrectionFactor;

        if (p>rnd()) flag=true,ParticleWeightCorrection=MinParticleWeightCorrectionFactor;
        else flag=false;
      }
    }
    while (flag==false);
  }
  else if ((injectionProbabilityPSD_LimitMin<=r)&&(r<=injectionProbabilityPSD_LimitMax)) { //injection through the PSD process

    int idim;
    double ExternalNormal[3]={0.0,0.0,0.0},p,l;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;

    //1. determine the position on the sphere where the new particle is generated
    do {  //calculate the probability of choosing the positions

      do { //generate the random position on the day side
        for (r=0.0,idim=0;idim<DIM;idim++) {
          ExternalNormal[idim]=sqrt(-2.0*log(rnd()))*cos(PiTimes2*rnd());
          r+=ExternalNormal[idim]*ExternalNormal[idim];
        }
      }
      while (ExternalNormal[0]>0.0);

      r=sqrt(r);

      for (idim=0;idim<DIM;idim++) {
        ExternalNormal[idim]/=r;
        x[idim]=-(rSphere+2.0*PIC::Mesh::mesh->EPS)*ExternalNormal[idim];
      }

      //r -> the random position in teh dayside of the sphere
      //p -> the probability to generate particle at the location 'r'
      //calculate the Longitude and the Latitude of the position (Lon, Lat)

      p=x[0]/rSphere;
    }
    while (p<rnd());

    node=PIC::Mesh::mesh->findTreeNode(x,startNode);
    if (node->Thread!=PIC::Mesh::mesh->ThisThread) return false;
    startNode=node;

    //2. determine the velocity  vector of the newly generated particle
    double e,speed,emin,emax;
    int eLevel;

    eLevel=(int)(rnd()*nIntervalsLoockupTable);
    emin=EnergyDistributionPSDTable[eLevel];
    emax=(eLevel<nIntervalsLoockupTable-1) ? EnergyDistributionPSDTable[eLevel+1] : EnergyDistributionPSDTable[eLevel];

    e=emin+rnd()*(emax-emin);
    speed=sqrt(2.0*e/PIC::MolecularData::GetMass(NA));

    //destribute the direction of the particle velocity vector
    do {
      l=0.0,p=0.0;

      for (idim=0;idim<3;idim++) {
        v[idim]=sqrt(-2.0*log(rnd()))*cos(PiTimes2*rnd());

        l+=v[idim]*v[idim];
        p+=v[idim]*ExternalNormal[idim];
      }
    }
    while (p>0.0);

    for (l=speed/sqrt(l),idim=0;idim<3;idim++) v[idim]*=l;




  }
  else exit(__LINE__,__FILE__,"Error: the injection process is not determined");


  return true;
}





//the mesh resolution
double localSphericalSurfaceResolution(double *x) {
  double res,r,l[3];
  int idim;
  double SubsolarAngle;

  for (r=0.0,idim=0;idim<3;idim++) r+=pow(x[idim],2);
  for (r=sqrt(r),idim=0;idim<3;idim++) l[idim]=x[idim]/r;

  SubsolarAngle=acos(l[0]);


  SubsolarAngle=0.0;

  res=dxMinSphere+(dxMaxSphere-dxMinSphere)/Pi*SubsolarAngle;

  #if _PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_
  res*=2*2.2;
  #elif _MERCURY_MESH_RESOLUTION_MODE_ == _MERCURY_MESH_RESOLUTION_MODE__FULL_
  res/=2.2;
  #else
  exit(__LINE__,__FILE__,"Error: the option is unknown");
  #endif

  return rSphere*res;
}

double localResolution(double *x) {
  int idim;
  double lnR,res,r=0.0;

  for (idim=0;idim<DIM;idim++) r+=pow(x[idim],2);

  r=sqrt(r);

  if (r>dxMinGlobal*rSphere) {
    lnR=log(r);
    res=dxMinGlobal+(dxMaxGlobal-dxMinGlobal)/log(xMaxDomain*rSphere)*lnR;
  }
  else res=dxMinGlobal;


  #if _PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_
  res*=4*2.2;
  #elif _MERCURY_MESH_RESOLUTION_MODE_ == _MERCURY_MESH_RESOLUTION_MODE__FULL_
  //do nothing
  #else
  exit(__LINE__,__FILE__,"Error: the option is unknown");
  #endif


  return rSphere*res;
}

//set up the local time step

double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  double CellSize,CharacteristicSpeed;

  switch (spec) {
  case _NA_SPEC_:
    CharacteristicSpeed=5.0E3;
    break;
  case _NA_PLUS_SPEC_:
    CharacteristicSpeed=400.0E3;
    break;
  case _H_PLUS_SPEC_: case _HE_2PLUS_SPEC_:
    CharacteristicSpeed=1500.0E3;
    break;
  default:
    exit(__LINE__,__FILE__,"Error: unknown species");
  }


  CellSize=startNode->GetCharacteristicCellSize();
  return 0.3*CellSize/CharacteristicSpeed;
}

double localParticleInjectionRate(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {

  /*
  bool ExternalFaces[6];
  double res=0.0,ExternalNormal[3],BlockSurfaceArea,ModelParticlesInjectionRate;
  int nface;

  static double vNA[3]={0.0,000.0,000.0},nNA=5.0E6,tempNA=8.0E4;
*/

  return 0.0;

  /*
  if (spec!=NA) return 0.0;

  if (PIC::Mesh::mesh->ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
      startNode->GetExternalNormal(ExternalNormal,nface);
      BlockSurfaceArea=startNode->GetBlockFaceSurfaceArea(nface);

      ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(nNA,tempNA,vNA,ExternalNormal,NA);

      res+=ModelParticlesInjectionRate*BlockSurfaceArea;
    }
  }

  return res;
  */
}

bool BoundingBoxParticleInjectionIndicator(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  bool ExternalFaces[6];
  double ExternalNormal[3],ModelParticlesInjectionRate;
  int nface;

  static double vNA[3]={0.0,000.0,000.0},nNA=5.0E6,tempNA=8.0E4;

  if (PIC::Mesh::mesh->ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
      startNode->GetExternalNormal(ExternalNormal,nface);
      ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(nNA,tempNA,vNA,ExternalNormal,NA);

      if (ModelParticlesInjectionRate>0.0) return true;
    }
  }

  return false;
}



//injection of model particles through the faces of the bounding box
int  BoundingBoxInjection(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  bool ExternalFaces[6];
  double ParticleWeight,LocalTimeStep,TimeCounter,ExternalNormal[3],x[3],x0[3],e0[3],e1[3],c0,c1;
  int nface,idim;
//  int nInjectedParticles;
  int newParticle;
  PIC::ParticleBuffer::byte *newParticleData;
  int nInjectedParticles=0;

  if (spec!=NA) return 0; //inject only spec=0

  static double vNA[3]={0.0,000.0,000.0},nNA=5.0E6,tempNA=8.0E4;
  double v[3];


  double ModelParticlesInjectionRate;

  if (PIC::Mesh::mesh->ExternalBoundaryBlock(startNode,ExternalFaces)==_EXTERNAL_BOUNDARY_BLOCK_) {
    ParticleWeight=startNode->block->GetLocalParticleWeight(spec);
    LocalTimeStep=startNode->block->GetLocalTimeStep(spec);


    for (nface=0;nface<2*DIM;nface++) if (ExternalFaces[nface]==true) {
      startNode->GetExternalNormal(ExternalNormal,nface);
      TimeCounter=0.0;

      ModelParticlesInjectionRate=PIC::BC::CalculateInjectionRate_MaxwellianDistribution(nNA,tempNA,vNA,ExternalNormal,NA);


      if (ModelParticlesInjectionRate>0.0) {
        ModelParticlesInjectionRate*=startNode->GetBlockFaceSurfaceArea(nface)/ParticleWeight;

        PIC::Mesh::mesh->GetBlockFaceCoordinateFrame_3D(x0,e0,e1,nface,startNode);

        while ((TimeCounter+=-log(rnd())/ModelParticlesInjectionRate)<LocalTimeStep) {
          //generate the new particle position on the face
          for (idim=0,c0=rnd(),c1=rnd();idim<DIM;idim++) x[idim]=x0[idim]+c0*e0[idim]+c1*e1[idim];

          //generate a particle
          newParticle=PIC::ParticleBuffer::GetNewParticle();
          newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
          nInjectedParticles++;

          PIC::BC::CalculateInjectionRate_MaxwellianDistribution(nNA,tempNA,vNA,ExternalNormal,NA);

          PIC::ParticleBuffer::SetX(x,newParticleData);
          PIC::ParticleBuffer::SetV(v,newParticleData);
          PIC::ParticleBuffer::SetI(spec,newParticleData);

          PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,newParticleData);


          //inject the particle into the system
//          PIC::Mover::MoveParticleTimeStep[spec](newParticle,LocalTimeStep-TimeCounter,startNode);

        }
      }


    }
  }

  return nInjectedParticles;
}

int BoundingBoxInjection(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  int nInjectedParticles=0;



//  for (int s=0;s<PIC::nTotalSpecies;s++) nInjectedParticles+=BoundingBoxInjection(s,startNode);

  return nInjectedParticles;
}






double InitLoadMeasure(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  double res=1.0;

 // for (int idim=0;idim<DIM;idim++) res*=(node->xmax[idim]-node->xmin[idim]);

  return res;
}

int ParticleSphereInteraction(int spec,int ptr,double *x,double *v,double &dtTotal,void *NodeDataPonter,void *SphereDataPointer)  {
   double radiusSphere,*x0Sphere,l[3],r,vNorm,c;
   cInternalSphericalData *Sphere;
   cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode;
   int idim;

//   int newParticle;
//   PIC::ParticleBuffer::byte *newParticleData;
//   double ParticleStatWeight,WeightCorrection;


   Sphere=(cInternalSphericalData*)SphereDataPointer;
   startNode=(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*)NodeDataPonter;

   Sphere->GetSphereGeometricalParameters(x0Sphere,radiusSphere);

   for (r=0.0,idim=0;idim<DIM;idim++) {
     l[idim]=x[idim]-x0Sphere[idim];
     r+=pow(l[idim],2);
   }

   for (r=sqrt(r),vNorm=0.0,idim=0;idim<DIM;idim++) vNorm+=v[idim]*l[idim]/r;
   if (vNorm<0.0) for (c=2.0*vNorm/r,idim=0;idim<DIM;idim++) v[idim]-=c*l[idim];

   //sample the particle data
   double *SampleData;
   int nSurfaceElement,nZenithElement,nAzimuthalElement;

   Sphere->GetSurfaceElementProjectionIndex(x,nZenithElement,nAzimuthalElement);
   nSurfaceElement=Sphere->GetLocalSurfaceElementNumber(nZenithElement,nAzimuthalElement);
   SampleData=Sphere->SamplingBuffer+PIC::BC::InternalBoundary::Sphere::collectingSpecieSamplingDataOffset(spec,nSurfaceElement);


   SampleData[PIC::BC::InternalBoundary::Sphere::sampledFluxDownRelativeOffset]+=startNode->block->GetLocalParticleWeight(spec)/startNode->block->GetLocalTimeStep(spec)/Sphere->GetSurfaceElementArea(nZenithElement,nAzimuthalElement);


//   if (x[0]*x[0]+x[1]*x[1]+x[2]*x[2]<0.9*rSphere*rSphere) exit(__LINE__,__FILE__,"Particle inside the sphere");


   r=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);


   //particle-surface interaction
   if (false) { /////(spec!=NA) { //surface reactiona
     exit(__LINE__,__FILE__,"no BC for the space is implemented");

		/*
     ParticleStatWeight=startNode->block->GetLocalParticleWeight(0);
     ParticleStatWeight*=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ptr);

     //model the rejected H+ and neutralized H
     //1. model rejected SW protons (SPEC=1)
     newParticle=PIC::ParticleBuffer::GetNewParticle();
     newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);

     PIC::ParticleBuffer::CloneParticle(newParticle,ptr);

     PIC::ParticleBuffer::SetV(v,newParticleData);
     PIC::ParticleBuffer::SetX(x,newParticleData);
     PIC::ParticleBuffer::SetI(1,newParticleData);


     //Set the correction of the individual particle weight
     WeightCorrection=0.01*ParticleStatWeight/startNode->block->GetLocalParticleWeight(1);
     PIC::ParticleBuffer::SetIndividualStatWeightCorrection(WeightCorrection,newParticleData);

     PIC::Mover::MoveParticleTimeStep[1](newParticle,dtTotal,startNode);

     //1. model rejected SW NEUTRALIZED protons (SPEC=2)
     newParticle=PIC::ParticleBuffer::GetNewParticle();
     newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);

     PIC::ParticleBuffer::CloneParticle(newParticle,ptr);

     PIC::ParticleBuffer::SetV(v,newParticleData);
     PIC::ParticleBuffer::SetX(x,newParticleData);
     PIC::ParticleBuffer::SetI(2,newParticleData);


     //Set the correction of the individual particle weight
     WeightCorrection=0.2*ParticleStatWeight/startNode->block->GetLocalParticleWeight(2);
     PIC::ParticleBuffer::SetIndividualStatWeightCorrection(WeightCorrection,newParticleData);

     PIC::Mover::MoveParticleTimeStep[2](newParticle,dtTotal,startNode);

     //remove the oroginal particle (s=0)
     PIC::ParticleBuffer::DeleteParticle(ptr);
     return _PARTICLE_DELETED_ON_THE_FACE_;
     */
   }


   //delete all particles that was not reflected on the surface
   PIC::ParticleBuffer::DeleteParticle(ptr);
   return _PARTICLE_DELETED_ON_THE_FACE_;
}


/*
void prePopulateSWprotons(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  double x[3],v[3],aNpart,NodeMeasure=1.0;
  int nPart,idim,i,j,k;

  static const double nSW=5.0E6; //solar wind number denstiy
  static const double tempSW=8.0E4;
  static double swVel[3]={4.0E5,0.0,0.0};

  int newParticle,nd;
  PIC::ParticleBuffer::byte *newParticleData;

  static int nTotalGeneratedParticles=0,nTotalProcessorBlocks=0;
  static double GlobalParticleWeight=0.0,aNpartTotal=0.0,TotalDomainVolume=0.0;


  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    if (startNode->Thread==PIC::Mesh::mesh->ThisThread) {
      nTotalProcessorBlocks++;


      //place particles in the blocks
      for (idim=0;idim<DIM;idim++) NodeMeasure*=(startNode->xmax[idim]-startNode->xmin[idim]);
      aNpart=NodeMeasure*nSW/startNode->block->GetLocalParticleWeight(SW);

      aNpartTotal+=aNpart;
      TotalDomainVolume+=NodeMeasure;

      GlobalParticleWeight=startNode->block->GetLocalParticleWeight(SW);

      nPart=(int)aNpart;
      aNpart-=nPart;
      if (aNpart>rnd()) nPart++;


      //generate particles
      for (;nPart>0;nPart--) {

        //generate a random particle's position
        for (idim=0;idim<DIM;idim++) x[idim]=startNode->xmin[idim]+(startNode->xmax[idim]-startNode->xmin[idim])*rnd();
        if (x[0]*x[0]+x[1]*x[1]+x[2]*x[2]<rSphere*rSphere) continue;

        nTotalGeneratedParticles++;

        PIC::Mesh::mesh->fingCellIndex(x,i,j,k,startNode);
        nd=startNode->block->getCenterNodeLocalNumber(i,j,k);

        newParticle=PIC::ParticleBuffer::GetNewParticle(startNode->block->GetCenterNode(nd)->FirstCellParticle);
        newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);

        PIC::Distribution::MaxwellianVelocityDistribution(v,swVel,tempSW,SW);

        PIC::ParticleBuffer::SetV(v,newParticleData);
        PIC::ParticleBuffer::SetX(x,newParticleData);
        PIC::ParticleBuffer::SetI(SW,newParticleData);
        PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,newParticleData);

      }


    }

  }
  else for (int nDownNode=0;nDownNode<8;nDownNode++) prePopulateSWprotons(startNode->downNode[nDownNode]);


  if (startNode==PIC::Mesh::mesh->rootTree) {
    int *GeneratedParticle=new int [PIC::Mesh::mesh->nTotalThreads];
    int *GeneratedNodes=new int [PIC::Mesh::mesh->nTotalThreads];
    double *anpart=new double [PIC::Mesh::mesh->nTotalThreads];
    double *volume=new double [PIC::Mesh::mesh->nTotalThreads];

    MPI_Gather(&nTotalGeneratedParticles,1,MPI_LONG,GeneratedParticle,1,MPI_LONG,0,MPI_GLOBAL_COMMUNICATOR);
    MPI_Gather(&nTotalProcessorBlocks,1,MPI_LONG,GeneratedNodes,1,MPI_LONG,0,MPI_GLOBAL_COMMUNICATOR);
    MPI_Gather(&aNpartTotal,1,MPI_DOUBLE,anpart,1,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);
    MPI_Gather(&TotalDomainVolume,1,MPI_DOUBLE,volume,1,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);


    if (PIC::Mesh::mesh->ThisThread==0) {
      cout << "Pre Population of the domain:\n Thread\tGenerated Particles\tDomain Block Number\t aNpartTotal\tSubDomain Volume" << endl;

      for (int t=0;t<PIC::Mesh::mesh->nTotalThreads;t++) cout << t << "\t" << GeneratedParticle[t] << "\t" << GeneratedNodes[t] << "\t" << anpart[t] << "\t" << volume[t] << endl;

      cout << "Global Particle's weight=" << GlobalParticleWeight << endl;
    }

    delete [] GeneratedNodes;
    delete [] GeneratedParticle;
    delete [] anpart;
    delete [] volume;
  }

}
*/

double sphereInjectionRate(int spec,void *SphereDataPointer) {
  double res=0.0;

  if (spec==NA) res=sodiumTotalProductionRate();
  else if (spec==NAPLUS) res=0.0;
  else exit(__LINE__,__FILE__,"Error: the source rate for the species is not determined");


  return res;
}

int sphereParticleInjection(void *SphereDataPointer) {
  cInternalSphericalData *Sphere;
  double ParticleWeight,LocalTimeStep,/*ExternalNormal[3],*/x[3],v[3],/*r,*/*sphereX0,sphereRadius;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=NULL;
  int newParticle,nInjectedParticles=0;
  PIC::ParticleBuffer::byte *newParticleData;
//  int idim;

  double ParticleWeightCorrection;

//  static const double Temp=200.0;
//  double vbulk[3]={0.0,0.0,0.0};


//  return 0;

//====================  DEBUG ===========================
//  static bool FirstPArticleGenerated=false;
//====================  END DEBUG ===================================


  Sphere=(cInternalSphericalData*)SphereDataPointer;
  Sphere->GetSphereGeometricalParameters(sphereX0,sphereRadius);

#if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
  ParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[NA];
#else
  exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
#endif


#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  LocalTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[NA];
#elif _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
  LocalTimeStep=Sphere->maxIntersectedNodeTimeStep[NA];
#else
  exit(__LINE__,__FILE__,"Error: the time step node is not defined");
#endif

  /*
  double TimeCounter=0.0;
  double ModelParticlesInjectionRate=sphereInjectionRate(NA,SphereDataPointer)/ParticleWeight;

  while ((TimeCounter+=-log(rnd())/ModelParticlesInjectionRate)<LocalTimeStep) {

*/

  static double InjectParticles=0.0;

  bool flag;

  InjectParticles+=sphereInjectionRate(NA,SphereDataPointer)*LocalTimeStep;

  while (InjectParticles>0.0) {

    /*
    //generate the particle position
    for (r=0.0,idim=0;idim<DIM;idim++) {
      ExternalNormal[idim]=sqrt(-2.0*log(rnd()))*cos(PiTimes2*rnd());
      r+=ExternalNormal[idim]*ExternalNormal[idim];
    }

    r=-sqrt(r);

    for (idim=0;idim<DIM;idim++) {
      ExternalNormal[idim]/=r;
      x[idim]=sphereX0[idim]-sphereRadius*ExternalNormal[idim];
    }

    startNode=PIC::Mesh::mesh->findTreeNode(x,startNode);
    if (startNode->Thread!=PIC::Mesh::mesh->ThisThread) continue;

    //generate the particle velocity
    PIC::Distribution::InjectMaxwellianDistribution(v,vbulk,Temp,ExternalNormal,NA);
    */

//====================  DEBUG ===========================
/*
    if (FirstPArticleGenerated==true) break;
    FirstPArticleGenerated=true;
    */
//====================  END DEBUG ===================================


    startNode=NULL;


    //if (sodiumDistributeParticleParameters(x,v,startNode,ParticleWeightCorrection)==false) continue;


    flag=sodiumDistributeParticleParameters(x,v,startNode,ParticleWeightCorrection);
    InjectParticles-=ParticleWeight*ParticleWeightCorrection;
    if (flag==false) continue;

//====================  DEBUG ===========================
    {
static double InjectionRadialVelocity=0.0,InjectionTangentionalSpeed=0.0;
static int nTotalInjectedParticles=0;

double l[3],r=0.0,v0=0.0,v1=0.0;
int idim;

for (idim=0;idim<3;idim++) r+=pow(x[idim],2);
r=sqrt(r);
for (idim=0;idim<3;idim++) {
  l[idim]=x[idim]/r;

  v0+=v[idim]*l[idim];
}

for (idim=0;idim<3;idim++) v1+=pow(v[idim]-v0*l[idim],2);

nTotalInjectedParticles++;
InjectionRadialVelocity+=v0;
InjectionTangentionalSpeed+=sqrt(v1);
    }
//====================  END DEBUG ===================================




    if (startNode->block->GetLocalTimeStep(NA)/LocalTimeStep<rnd()) continue;

    //generate a particle
    newParticle=PIC::ParticleBuffer::GetNewParticle();
    newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
    nInjectedParticles++;

    PIC::ParticleBuffer::SetX(x,newParticleData);
    PIC::ParticleBuffer::SetV(v,newParticleData);
    PIC::ParticleBuffer::SetI(NA,newParticleData);

    PIC::ParticleBuffer::SetIndividualStatWeightCorrection(ParticleWeightCorrection,newParticleData);


    //inject the particle into the system
    _PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_(newParticle,startNode->block->GetLocalTimeStep(NA)*rnd() /*LocalTimeStep-TimeCounter*/,startNode,true);
//    PIC::Mover::MoveParticleBoundaryInjection[NA](newParticle,startNode->block->GetLocalTimeStep(NA)*rnd() /*LocalTimeStep-TimeCounter*/,startNode,true);
  }

  return nInjectedParticles;
}

void amps_init_mesh() {
//  MPI_Init(&argc,&argv);
  PIC::InitMPI();


  //SetUp the alarm
//  PIC::Alarm::SetAlarm(2000);



//  VT_OFF();
//  VT_traceoff();


#ifdef _MAIN_PROFILE_
  initSaturn ("");
#endif


  rnd_seed();



  char inputFile[]="mercury.input";



  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);


  //init the Mercury model
  Mercury::Init_BeforeParser();

  //init the particle solver

  PIC::Init_BeforeParser();
//  PIC::Parser::Run(inputFile);


  PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber=400; //*10;
  PIC::RequiredSampleLength=200; //0; //0;

  Mercury::OrbitalMotion::nOrbitalPositionOutputMultiplier=10;




  //output the PDS enerfy distribution function
/*  if (PIC::ThisThread==0) {
    Mercury::SourceProcesses::PhotonStimulatedDesorption::EnergyDistribution.fPrintCumulativeDistributionFunction("CumulativeEnergyDistribution-PSD.dat");
    Mercury::SourceProcesses::PhotonStimulatedDesorption::EnergyDistribution.fPrintDistributionFunction("EnergyDistribution-PSD.dat");

    cout << Mercury::SourceProcesses::PhotonStimulatedDesorption::EnergyDistribution.DistributeVariable() << endl;
  }*/

  if (PIC::ThisThread==0) {
    for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
      char fname[200];

      sprintf(fname,"%s/CumulativeEnergyDistribution-PSD.nspec=%i.%s.dat",PIC::OutputDataFileDirectory,spec,PIC::MolecularData::GetChemSymbol(spec));
      Mercury::SourceProcesses::PhotonStimulatedDesorption::EnergyDistribution[spec].fPrintCumulativeDistributionFunction(fname);

      sprintf(fname,"%s/EnergyDistribution-PSD.nspec=%i.%s.dat",PIC::OutputDataFileDirectory,spec,PIC::MolecularData::GetChemSymbol(spec));
      Mercury::SourceProcesses::PhotonStimulatedDesorption::EnergyDistribution[spec].fPrintDistributionFunction(fname,&spec);
    }

//    cout << Moon::SourceProcesses::PhotonStimulatedDesorption::EnergyDistribution.DistributeVariable() << endl;
  }



  //register the sphere
  {
    double sx0[3]={0.0,0.0,0.0};
    cInternalBoundaryConditionsDescriptor SphereDescriptor;
    cInternalSphericalData *Sphere;


    //reserve memory for sampling of the surface balance of sticking species
    int ReserveSamplingSpace[PIC::nTotalSpecies];

    for (int s=0;s<PIC::nTotalSpecies;s++) ReserveSamplingSpace[s]=_OBJECT_SURFACE_SAMPLING__TOTAL_SAMPLED_VARIABLES_;


    cInternalSphericalData::SetGeneralSurfaceMeshParameters(60,100);



    PIC::BC::InternalBoundary::Sphere::Init(ReserveSamplingSpace,NULL);
    SphereDescriptor=PIC::BC::InternalBoundary::Sphere::RegisterInternalSphere();
    Sphere=(cInternalSphericalData*) SphereDescriptor.BoundaryElement;
    Sphere->SetSphereGeometricalParameters(sx0,rSphere);


    //init the object for distribution of the injection surface elements
/*    Mercury::SourceProcesses::PhotonStimulatedDesorption::SurfaceInjectionDistribution.SetLimits(0,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber-1,
        Mercury::SourceProcesses::PhotonStimulatedDesorption::GetSurfaceElementSodiumProductionRate);


  #if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__ON_
    Mercury::SourceProcesses::ThermalDesorption::SurfaceInjectionDistribution.SetLimits(0,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber-1,
        Mercury::SourceProcesses::ThermalDesorption::GetSurfaceElementSodiumProductionRate);
  #endif

  #if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
    Mercury::SourceProcesses::SolarWindSputtering::SurfaceInjectionDistribution.SetLimits(0,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber-1,
        Mercury::SourceProcesses::SolarWindSputtering::GetSurfaceElementSodiumProductionRate);
  #endif



    Sphere->PrintSurfaceMesh("Sphere.dat");
    Sphere->PrintSurfaceData("SpheraData.dat",0);*/


    for (int spec=0;spec<PIC::nTotalSpecies;spec++) {
    	Mercury::SourceProcesses::PhotonStimulatedDesorption::SurfaceInjectionDistribution[spec].SetLimits(0,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber-1,
    			Mercury::SourceProcesses::PhotonStimulatedDesorption::GetSurfaceElementProductionRate);


    #if _EXOSPHERE_SOURCE__THERMAL_DESORPTION_ == _EXOSPHERE_SOURCE__ON_
    	Mercury::SourceProcesses::ThermalDesorption::SurfaceInjectionDistribution[spec].SetLimits(0,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber-1,
    			Mercury::SourceProcesses::ThermalDesorption::GetSurfaceElementProductionRate);
    #endif

    #if _EXOSPHERE_SOURCE__SOLAR_WIND_SPUTTERING_ == _EXOSPHERE_SOURCE__ON_
    	Mercury::SourceProcesses::SolarWindSputtering::SurfaceInjectionDistribution[spec].SetLimits(0,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber-1,
    			Mercury::SourceProcesses::SolarWindSputtering::GetSurfaceElementProductionRate);
    #endif
    }


    char fname[_MAX_STRING_LENGTH_PIC_];

    sprintf(fname,"%s/Sphere.dat",PIC::OutputDataFileDirectory);
    Sphere->PrintSurfaceMesh(fname);

    sprintf(fname,"%s/SpheraData.dat",PIC::OutputDataFileDirectory);
    Sphere->PrintSurfaceData(fname,0);


    Sphere->Radius=_RADIUS_(_TARGET_);
    Sphere->localResolution=localSphericalSurfaceResolution;
    Sphere->InjectionRate=Mercury::SourceProcesses::totalProductionRate;
    Sphere->faceat=0;
    Sphere->ParticleSphereInteraction=Mercury::SurfaceInteraction::ParticleSphereInteraction_SurfaceAccomodation;
    Sphere->InjectionBoundaryCondition=Mercury::SourceProcesses::InjectionBoundaryModel; ///sphereParticleInjection;

    Sphere->PrintTitle=Mercury::Sampling::OutputSurfaceDataFile::PrintTitle;
    Sphere->PrintVariableList=Mercury::Sampling::OutputSurfaceDataFile::PrintVariableList;
    Sphere->PrintDataStateVector=Mercury::Sampling::OutputSurfaceDataFile::PrintDataStateVector;

    //set up the planet pointer in Mercury model
    Mercury::Planet=Sphere;
    Sphere->Allocate<cInternalSphericalData>(PIC::nTotalSpecies,PIC::BC::InternalBoundary::Sphere::TotalSurfaceElementNumber,_EXOSPHERE__SOURCE_MAX_ID_VALUE_,Sphere);



  }

//=============   DEBUG ==============

  /*
  double xProbeMin[3]={-2172234,-1.236913E-10,2000742};
  double dx=-2172234-(-2229398),v=1.0,t;
  double xProbeMax[3];
  int IntersectionStatus;

  for (int i=0;i<3;i++) {
    xProbeMax[i]=xProbeMin[i]+dx;
    v*=dx;
  }

  t=Sphere->GetRemainedBlockVolume(xProbeMin,xProbeMax,0.0,&IntersectionStatus);
*/

//============= END DEBUG =============

  //init the solver
  PIC::Mesh::initCellSamplingDataBuffer();

  //init the mesh
  cout << "Init the mesh" << endl;

  int maxBlockCellsnumber,minBlockCellsnumber,idim;

  maxBlockCellsnumber=_BLOCK_CELLS_X_;
  if (DIM>1) maxBlockCellsnumber=max(maxBlockCellsnumber,_BLOCK_CELLS_Y_);
  if (DIM>2) maxBlockCellsnumber=max(maxBlockCellsnumber,_BLOCK_CELLS_Z_);

  minBlockCellsnumber=_BLOCK_CELLS_X_;
  if (DIM>1) minBlockCellsnumber=min(minBlockCellsnumber,_BLOCK_CELLS_Y_);
  if (DIM>2) minBlockCellsnumber=min(minBlockCellsnumber,_BLOCK_CELLS_Z_);

  double DomainLength[3],DomainCenterOffset[3],xmax[3]={0.0,0.0,0.0},xmin[3]={0.0,0.0,0.0};

  if (maxBlockCellsnumber==minBlockCellsnumber) {
    for (idim=0;idim<DIM;idim++) {
      DomainLength[idim]=2.0*xMaxDomain*rSphere;
      DomainCenterOffset[idim]=-xMaxDomain*rSphere;
    }
  }
  else {
    if (maxBlockCellsnumber!=_BLOCK_CELLS_X_) exit(__LINE__,__FILE__);
    if (minBlockCellsnumber!=_BLOCK_CELLS_Y_) exit(__LINE__,__FILE__);
    if (minBlockCellsnumber!=_BLOCK_CELLS_Z_) exit(__LINE__,__FILE__);

    DomainLength[0]=xMaxDomain*rSphere*(1.0+double(_BLOCK_CELLS_Y_)/_BLOCK_CELLS_X_);
    DomainLength[1]=DomainLength[0]*double(_BLOCK_CELLS_Y_)/_BLOCK_CELLS_X_;
    DomainLength[2]=DomainLength[0]*double(_BLOCK_CELLS_Z_)/_BLOCK_CELLS_X_;

    if (DomainLength[1]<2.01*yMaxDomain*_RADIUS_(_TARGET_)) {
      double r;

      printf("Size of the domain is smaller that the radius of the body: the size of the domain is readjusted\n");
      r=2.01*yMaxDomain*_RADIUS_(_TARGET_)/DomainLength[1];

      for (idim=0;idim<DIM;idim++) DomainLength[idim]*=r;
    }

    DomainCenterOffset[0]=-yMaxDomain*rSphere;////-xMaxDomain*rSphere*double(_BLOCK_CELLS_Y_)/_BLOCK_CELLS_X_;
    DomainCenterOffset[1]=-DomainLength[1]/2.0;
    DomainCenterOffset[2]=-DomainLength[2]/2.0;
  }

  /*
  for (idim=0;idim<DIM;idim++) {
    xmin[idim]=DomainCenterOffset[idim];
    xmax[idim]=DomainLength[idim]+DomainCenterOffset[idim];
  }
  */


#if _PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_
  for (idim=0;idim<DIM;idim++) {
    xmax[idim]=-DomainCenterOffset[idim];
    xmin[idim]=+DomainCenterOffset[idim];
  }
#else
  for (idim=0;idim<DIM;idim++) {
    xmax[idim]=-DomainCenterOffset[idim];
    xmin[idim]=-(DomainLength[idim]+DomainCenterOffset[idim]);
  }
#endif //_PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_

//  double xmin[3]={-xMaxDomain*rSphere*_BLOCK_CELLS_X_/double(maxBlockCellsnumber),-xMaxDomain*rSphere*_BLOCK_CELLS_Y_/double(maxBlockCellsnumber),-xMaxDomain*rSphere*_BLOCK_CELLS_Z_/double(maxBlockCellsnumber)};
//  double xmax[3]={xMaxDomain*rSphere*_BLOCK_CELLS_X_/double(maxBlockCellsnumber),xMaxDomain*rSphere*_BLOCK_CELLS_Y_/double(maxBlockCellsnumber),xMaxDomain*rSphere*_BLOCK_CELLS_Z_/double(maxBlockCellsnumber)};


  //generate only the tree
  PIC::Mesh::mesh->AllowBlockAllocation=false;
  PIC::Mesh::mesh->init(xmin,xmax,localResolution);
  PIC::Mesh::mesh->memoryAllocationReport();


//  VT_ON();
//  VT_USER_START("name");
//  VT_ON();

//  {
//    VT_TRACER("name");

  bool NewMeshGeneratedFlag=false;
  char mesh[200]="amr.bin";

  FILE *fmesh=fopen(mesh,"r");

  if (fmesh!=NULL) {
    fclose(fmesh);
    PIC::Mesh::mesh->readMeshFile(mesh);
  }
  else {
    NewMeshGeneratedFlag=true;

    if (PIC::Mesh::mesh->ThisThread==0) {
       PIC::Mesh::mesh->buildMesh();
       PIC::Mesh::mesh->saveMeshFile("mesh.msh");
       MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    }
    else {
       MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
       PIC::Mesh::mesh->readMeshFile("mesh.msh");
    }
  }

  if (NewMeshGeneratedFlag==true) PIC::Mesh::mesh->outputMeshTECPLOT("mesh.dat");

  PIC::Mesh::mesh->memoryAllocationReport();
  PIC::Mesh::mesh->GetMeshTreeStatistics();

#ifdef _CHECK_MESH_CONSISTENCY_
  PIC::Mesh::mesh->checkMeshConsistency(PIC::Mesh::mesh->rootTree);
#endif

  PIC::Mesh::mesh->SetParallelLoadMeasure(InitLoadMeasure);
  PIC::Mesh::mesh->CreateNewParallelDistributionLists();

  //initialize the blocks
  PIC::Mesh::mesh->AllowBlockAllocation=true;
  PIC::Mesh::mesh->AllocateTreeBlocks();

  PIC::Mesh::mesh->memoryAllocationReport();
  PIC::Mesh::mesh->GetMeshTreeStatistics();

#ifdef _CHECK_MESH_CONSISTENCY_
  PIC::Mesh::mesh->checkMeshConsistency(PIC::Mesh::mesh->rootTree);
#endif

  //get the mesh signature
  unsigned int MeshSignature;
  MeshSignature=PIC::Mesh::mesh->getMeshSignature(); 

  //init the volume of the cells'
  PIC::Mesh::mesh->InitCellMeasure();

  //if the new mesh was generated => rename created mesh.msh into amr.sig=0x%lx.mesh.bin
  if (NewMeshGeneratedFlag==true) {
    unsigned int MeshSignature=PIC::Mesh::mesh->getMeshSignature();

    if (PIC::Mesh::mesh->ThisThread==0) {
      char command[300];

      sprintf(command,"mv mesh.msh amr.sig=0x%lx.mesh.bin",MeshSignature);
      system(command);
    }
  }

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);


  if (PIC::ThisThread==0) cout << "AMPS' Initialization is complete" << endl;

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
}

void amps_init() {
  //init the boundary injection procedure
  PIC::BC::BlockInjectionBCindicatior=BoundingBoxParticleInjectionIndicator;
  PIC::BC::userDefinedBoundingBlockInjectionFunction=BoundingBoxInjection;
  PIC::BC::InitBoundingBoxInjectionBlockList();


  //load the background plasma data
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__DATAFILE_
#if _PIC_COUPLER_DATAFILE_READER_MODE_ == _PIC_COUPLER_DATAFILE_READER_MODE__TECPLOT_

  if (PIC::CPLR::DATAFILE::BinaryFileExists("MERCURY-BATSRUS")==true)  {
    PIC::CPLR::DATAFILE::LoadBinaryFile("MERCURY-BATSRUS");
  }
  else {

    double xminTECPLOT[3]={-32,-32,-32},xmaxTECPLOT[3]={16,32,32};

    PIC::CPLR::DATAFILE::TECPLOT::UnitLength=_RADIUS_(_TARGET_);
    PIC::CPLR::DATAFILE::TECPLOT::SetDomainLimitsXYZ(xminTECPLOT,xmaxTECPLOT);
    PIC::CPLR::DATAFILE::TECPLOT::SetDomainLimitsSPHERICAL(1.0,50.0);

    PIC::CPLR::DATAFILE::TECPLOT::DataMode=PIC::CPLR::DATAFILE::TECPLOT::DataMode_SPHERICAL;
    PIC::CPLR::DATAFILE::TECPLOT::SetLoadedVelocityVariableData(5,1.0E3);
    PIC::CPLR::DATAFILE::TECPLOT::SetLoadedIonPressureVariableData(14,1.0E-9);
    PIC::CPLR::DATAFILE::TECPLOT::SetLoadedMagneticFieldVariableData(11,1.0E-9);
    PIC::CPLR::DATAFILE::TECPLOT::SetLoadedDensityVariableData(4,1.0E6);
    PIC::CPLR::DATAFILE::TECPLOT::nTotalVarlablesTECPLOT=21;
    PIC::CPLR::DATAFILE::TECPLOT::ImportData("3d__var_7_t00000200_n0300072.plt"); //data/input/Mercury/040915-Jia/3d__var_7_t00000200_n0300072.plt

    PIC::CPLR::DATAFILE::SaveBinaryFile("MERCURY-BATSRUS");
  }
#else
  exit(__LINE__,__FILE__,"ERROR: unrecognized datafile reader mode");

#endif //_PIC_COUPLER_DATAFILE_READER_MODE_
#endif //_PIC_COUPLER_MODE_

  Mercury::Init_AfterParser();

//init the PIC solver
  PIC::Init_AfterParser ();
	PIC::Mover::Init();
//	PIC::Mover::TotalParticleAcceleration=TotalParticleAcceleration;

//	for (int s=0;s<PIC::nTotalSpecies;s++) {
//	  PIC::Mover::MoveParticleTimeStep[s]=PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory_SecondOrder; ///UniformWeight_UniformTimeStep_SecondOrder;
//	  PIC::Mover::MoveParticleBoundaryInjection[s]=PIC::Mover::UniformWeight_UniformTimeStep_noForce_TraceTrajectory_BoundaryInjection_SecondOrder;
//	}


  //set up the time step
  PIC::ParticleWeightTimeStep::LocalTimeStep=localTimeStep;
  PIC::ParticleWeightTimeStep::initTimeStep();

  //set up the particle weight
  PIC::ParticleWeightTimeStep::UserDefinedExtraSourceRate=Exosphere::SourceProcesses::SolarWindSputtering::TypicalIonFluxSputteringRate;

  PIC::ParticleWeightTimeStep::LocalBlockInjectionRate=localParticleInjectionRate;
  PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(_NA_SPEC_);

  //copy the weight and time step from Na neutra to Na ions
  if (_NA_PLUS_SPEC_>0) {
    PIC::ParticleWeightTimeStep::copyLocalParticleWeightDistribution(_NA_PLUS_SPEC_,_NA_SPEC_,5.0E3/800.0E3);
    PIC::ParticleWeightTimeStep::copyLocalTimeStepDistribution(_NA_PLUS_SPEC_,_NA_SPEC_,5.0E3/800.0E3);
  }

  if (_H_PLUS_SPEC_>0) {
    PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(_H_PLUS_SPEC_);
  }

  if (_HE_2PLUS_SPEC_>0) {
    PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(_HE_2PLUS_SPEC_);
  }

  //set photolytic reactions
/*  PIC::ChemicalReactions::PhotolyticReactions::SetReactionProcessor(sodiumPhotoionizationReactionProcessor,NA);
  PIC::ChemicalReactions::PhotolyticReactions::SetSpeciesTotalPhotolyticLifeTime(sodiumPhotoionizationLifeTime,NA);*/


//  PIC::Mesh::mesh->outputMeshTECPLOT("mesh.dat");
//  PIC::Mesh::mesh->outputMeshDataTECPLOT("mesh.data.dat",0);

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  if (PIC::Mesh::mesh->ThisThread==0) cout << "The mesh is generated" << endl;


/*  PIC::Mesh::mesh->outputMeshDataTECPLOT("loaded.data.dat",0);


  PIC::CPLR::DATAFILE::SaveBinaryFile("SavedCellData.bin");
  PIC::CPLR::DATAFILE::LoadBinaryFile("SavedCellData.bin");
  PIC::Mesh::mesh->outputMeshDataTECPLOT("loaded.SavedCellData.dat",0);



  PIC::CPLR::SaveCenterNodeAssociatedData("SavedAssociatedData.bin");
  PIC::CPLR::LoadCenterNodeAssociatedData("SavedAssociatedData.bin");
  PIC::Mesh::mesh->outputMeshDataTECPLOT("loaded.SavedAssociatedData.dat",0);*/





	//output final data
//  PIC::Mesh::mesh->outputMeshDataTECPLOT("final.data.dat",0);

  //create the list of mesh nodes where the injection boundary conditinos are applied
/*  PIC::BC::BlockInjectionBCindicatior=BoundingBoxParticleInjectionIndicator;
  PIC::BC::userDefinedBoundingBlockInjectionFunction=BoundingBoxInjection;
  PIC::BC::InitBoundingBoxInjectionBlockList();*/


  //init the particle buffer
  PIC::ParticleBuffer::Init(1000000);
//  double TimeCounter=time(NULL);
//  int LastDataOutputFileNumber=-1;


  //init the sampling of the particls' distribution functions
/*
  const int nSamplePoints=3;
  double SampleLocations[nSamplePoints][DIM]={{7.6E5,6.7E5,0.0}, {2.8E5,5.6E5,0.0}, {-2.3E5,3.0E5,0.0}};

  PIC::DistributionFunctionSample::vMin=-40.0E3;
  PIC::DistributionFunctionSample::vMax=40.0E3;
  PIC::DistributionFunctionSample::nSampledFunctionPoints=500;

  PIC::DistributionFunctionSample::Init(SampleLocations,nSamplePoints);
*/


#if _MERCURY_FIPS_SAMPLING_ == _MERCURY_MODE_ON_
  //init sampling points along the s/c trajectory
  const int nFlybySamplePasses=6;
  const double FlybySamplingInterval=20.0*60.0,FlybySamplingIntervalStep=60.0; //in seconds
  const int nSampleSteps=(int)(FlybySamplingInterval/FlybySamplingIntervalStep);

  const char *FlybySamplePassesUTC[nFlybySamplePasses]={"2011-04-13T16:15:00","2011-04-14T16:30:00","2011-04-16T04:35:00",
      "2011-04-13T04:40:00","2011-04-15T04:55:00","2011-04-21T17:50:00"};

  SpiceDouble et,lt;
  SpiceDouble state[6];
  int nFlybyPass,n;

  /* FIPS POINTING
      INS-236720_FOV_FRAME       = 'MSGR_EPPS_FIPS'
      INS-236720_FOV_SHAPE       = 'CIRCLE'
      INS-236720_BORESIGHT       = ( 0.0, 0.0, 1.0 )
   */

  SpiceDouble pointing[3],bsight[3],bsight_INIT[3]={0.0,0.0,1.0};
  SpiceDouble rotate[3][3];

  const SpiceInt lenout = 35;
  SpiceChar utcstr[lenout+2];


  int nFluxSamplePoint=0;

  int nTotalFluxSamplePoints=nFlybySamplePasses*nSampleSteps;
  double FluxSampleLocations[nTotalFluxSamplePoints][3];
  double FluxSampleDirections[nTotalFluxSamplePoints][3];



/*


  for (nFlybyPass=0;nFlybyPass<nFlybySamplePasses;nFlybyPass++) {
    utc2et_c(FlybySamplePassesUTC[nFlybyPass],&et);

    if (PIC::ThisThread==0) {
      cout << "S/C Flyby Sampling: Pass=" << nFlybyPass << ":" << endl;
      cout << "Flux Sample Point\tUTS\t\t\t x[km]\t\ty[km]\t\tz[km]\t\t\t lx\t\tly\t\tlz\t" << endl;
    }

    for (n=0;n<nSampleSteps;n++) {
      //position of the s/c
      spkezr_c("MESSENGER",et,"MSGR_MSO","NONE","MERCURY",state,&lt);


      //get the pointing vector in the 'SO' frame
      memcpy(bsight,bsight_INIT,3*sizeof(double));

      pxform_c ("MSGR_EPPS_FIPS","MSGR_MSO",et,rotate);
      mxv_c(rotate,bsight,pointing);

      //print the pointing information
      if (PIC::ThisThread==0) {
        et2utc_c(et,"C",0,lenout,utcstr);
        printf("%i\t\t\t%s\t",nFluxSamplePoint,utcstr);
        for (idim=0;idim<3;idim++) printf("%e\t",state[idim]);

        cout << "\t";

        for (idim=0;idim<3;idim++) printf("%e\t",pointing[idim]);
        cout << endl;
      }

      //save the samlpe pointing information
      for (idim=0;idim<3;idim++) {
        FluxSampleLocations[nFluxSamplePoint][idim]=state[idim]*1.0E3;
        FluxSampleDirections[nFluxSamplePoint][idim]=pointing[idim];
      }

      //increment the flyby time
      et+=FlybySamplingIntervalStep;
      ++nFluxSamplePoint;
    }

    if (PIC::ThisThread==0) cout << endl;
  }

  PIC::ParticleFluxDistributionSample::Init(FluxSampleLocations,FluxSampleDirections,30.0/180.0*Pi,nTotalFluxSamplePoints);
*/

#elif _MERCURY_FIPS_SAMPLING_ == _MERCURY_MODE_OFF_
  printf("No FIPS sampling\n");
#else
  exit(__LINE__,__FILE__,"Error: the option is not recognized");
#endif




/*
#if _PIC_COUPLER_MODE_ ==       _PIC_COUPLER_MODE__ICES_
  //#ifdef _ICES_LOAD_DATA_
//  PIC::Mesh::mesh->outputMeshDataTECPLOT("ices.data.dat",0);

  //output the solar wind ion flux at the palnet's surface
  PIC::CPLR::ICES::PrintSphereSurfaceIonFlux("SurfaceIonFlux.dat",1.05*_RADIUS_(_TARGET_));
  PIC::CPLR::ICES::EvaluateSurfaceIonFlux(1.05);

  PIC::BC::InternalBoundary::Sphere::InternalSpheres.GetEntryPointer(0)->PrintSurfaceData("Surface.test.dat",0);

  Exosphere::SourceProcesses::Init();
  PIC::BC::InternalBoundary::Sphere::InternalSpheres.GetEntryPointer(0)->PrintSurfaceData("Surface.test-1.dat",0);
#endif
*/

  //init the source model
  Exosphere::SourceProcesses::Init();
  PIC::BC::InternalBoundary::Sphere::InternalSpheres.GetEntryPointer(0)->PrintSurfaceData("Surface.test.dat",0);

}







  //time step
// for (int niter=0;niter<100000001;niter++) {
void amps_time_step(){

  int idim;

    //determine the parameters of the orbital motion of Mercury
    SpiceDouble StateBegin[6],StateEnd[6],lt,StateSun[6],StateMiddle[6];
    double lBegin[3],rBegin,lEnd[3],rEnd,vTangentialBegin=0.0,vTangentialEnd=0.0,c0=0.0,c1=0.0;

    SpiceDouble HCI_to_MSO_TransformationMartix[6][6];

    spkezr_c("Mercury",Mercury::OrbitalMotion::et,"MSGR_HCI","none","SUN",StateBegin,&lt);
    spkezr_c("SUN",Mercury::OrbitalMotion::et,"MSGR_MSO","none","Mercury",StateSun,&lt);

    //calculate Mercury's velocity in an itertial frame, which have dirtectional vectors that coinsides with that of MSO
    sxform_c("MSGR_HCI","MSGR_MSO",Mercury::OrbitalMotion::et+0.5*PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_],HCI_to_MSO_TransformationMartix);
    spkezr_c("Mercury",Mercury::OrbitalMotion::et+0.5*PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_],"MSGR_HCI","none","SUN",StateMiddle,&lt);


    Mercury::OrbitalMotion::et+=PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_];
    spkezr_c("Mercury",Mercury::OrbitalMotion::et,"MSGR_HCI","none","SUN",StateEnd,&lt);


    for (rBegin=0.0,rEnd=0.0,idim=0;idim<3;idim++) {
      StateBegin[idim]*=1.0E3,StateBegin[3+idim]*=1.0E3;
      StateEnd[idim]*=1.0E3,StateEnd[3+idim]*=1.0E3;
      StateMiddle[idim]*=1.0E3,StateMiddle[3+idim]*=1.0E3;

      rBegin+=pow(StateBegin[idim],2);
      rEnd+=pow(StateEnd[idim],2);

      Mercury::xObject_HCI[idim]=StateBegin[idim];
      Mercury::vObject_HCI[idim]=StateBegin[3+idim];

      Mercury::xSun_SO[idim]=1.0E3*StateSun[idim];
      Mercury::vSun_SO[idim]=1.0E3*StateSun[3+idim];
    }

    //calculate parameters of SO_FROZEN
    //velocity of the coordinate frame
    for (idim=0;idim<3;idim++) {
      Mercury::vObject_SO_FROZEN[idim]=
          HCI_to_MSO_TransformationMartix[idim][0]*StateMiddle[3+0]+
          HCI_to_MSO_TransformationMartix[idim][1]*StateMiddle[3+1]+
          HCI_to_MSO_TransformationMartix[idim][2]*StateMiddle[3+2];
    }

    //the axis of rotation of the MSO fraim in MSO_FROZEN during the next time step
    //get pointing direction to the Sun at the end of the current iteration in MSO_FROZEN
    SpiceDouble fmatrix[6][6];
    double SunPointingDirectionEnd[3],SunPointingDirectionEnd_MSO_FROZEN[3];

    //calculate Sun pointing at the end of the iteration in HCI frame (et is already incremented!!!!!!)
    sxform_c("MSGR_MSO","MSGR_HCI",Mercury::OrbitalMotion::et,fmatrix);

    SunPointingDirectionEnd[0]=fmatrix[0][0];
    SunPointingDirectionEnd[1]=fmatrix[1][0];
    SunPointingDirectionEnd[2]=fmatrix[2][0];

    //convert the pointing direction vector into MSO_FROZEN frame
    sxform_c("MSGR_HCI","MSGR_MSO",Mercury::OrbitalMotion::et-PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_],fmatrix);

    for (idim=0;idim<3;idim++) {
      SunPointingDirectionEnd_MSO_FROZEN[idim]=
          fmatrix[idim][0]*SunPointingDirectionEnd[0]+
          fmatrix[idim][1]*SunPointingDirectionEnd[1]+
          fmatrix[idim][2]*SunPointingDirectionEnd[2];
    }

    //calculate the rate of rotation in MSO_FROZEN
    Mercury::RotationRate_SO_FROZEN=acos(SunPointingDirectionEnd_MSO_FROZEN[0])/PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_];


    //calculate the direction of rotation
    double c=sqrt(pow(SunPointingDirectionEnd_MSO_FROZEN[1],2)+pow(SunPointingDirectionEnd_MSO_FROZEN[2],2));

    if (c>0.0) {
      Mercury::RotationVector_SO_FROZEN[0]=0.0;
      Mercury::RotationVector_SO_FROZEN[1]=-SunPointingDirectionEnd_MSO_FROZEN[2]/c*Mercury::RotationRate_SO_FROZEN;
      Mercury::RotationVector_SO_FROZEN[2]=SunPointingDirectionEnd_MSO_FROZEN[1]/c*Mercury::RotationRate_SO_FROZEN;
    }
    else {
      Mercury::RotationVector_SO_FROZEN[0]=0.0;
      Mercury::RotationVector_SO_FROZEN[1]=0.0;
      Mercury::RotationVector_SO_FROZEN[2]=0.0;
    }


    //RECALCUALTE THE ROTATION VECTOR USING THE TRANSOFRMATON MARTICX FROM MSO_FROSEN at the time step (n) to the MSO_FROZEN at the time step (n+1)
    //the rotation vector is the eigrnvector of the transformation matrix
    //Zhuravlev, Osnovy teoreticheskoi mehaniki, Chapter 2, paragraph 6.2 (sposoby zadaniya orientacii tverdogo tela)

    //get the transformation matrix T(LSO[n]->LSO[n+1])=T1(LSO[n]->MSGR_HCI)*T2(MSGR_HCI->LSO[n+1])
    SpiceDouble T1[6][6],T2[6][6];
    double T[3][3];
    int i,j,k;


    /*
    sxform_c("MSGR_MSO","MSGR_HCI",Mercury::OrbitalMotion::et-PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_],T1);
    sxform_c("MSGR_HCI","MSGR_MSO",Mercury::OrbitalMotion::et,T2);
*/

    sxform_c("MSGR_MSO","MSGR_HCI",Mercury::OrbitalMotion::et,T1);
    sxform_c("MSGR_HCI","MSGR_MSO",Mercury::OrbitalMotion::et-PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_],T2);


/*
    furnsh_c("/Users/vtenishe/SPICE/Kernels/OTHER/Moon.LSO.tf");
    sxform_c("LSO","MSGR_HCI",Mercury::OrbitalMotion::et,T1);
    sxform_c("MSGR_HCI","LSO",Mercury::OrbitalMotion::et-PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_],T2);
*/


    for (i=0;i<3;i++) for (j=0;j<3;j++) {
      T[i][j]=0.0;

      for (k=0;k<3;k++) T[i][j]+=T2[i][k]*T1[k][j];
    }

    //determine the rate and the vectrot of the rotation
    double RotationAngle,t,RotationVector[3],RotationRate;

    RotationAngle=acos((T[0][0]+T[1][1]+T[2][2]-1.0)/2.0);

    t=2.0*sin(RotationAngle);
    RotationVector[0]=(T[2][1]-T[1][2])/t;
    RotationVector[1]=(T[0][2]-T[2][0])/t;
    RotationVector[2]=(T[1][0]-T[0][1])/t;

    RotationRate=RotationAngle/PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_];

//    t=RotationRate/sqrt(RotationVector[0]*RotationVector[0]+RotationVector[1]*RotationVector[1]+RotationVector[2]*RotationVector[2]);
//    RotationVector[0]*=t,RotationVector[1]*=t,RotationVector[2]*=t;

t=1.0;

    //TEST THE ROTATION RATE AND THE ROTATION VECTOR
    double testRoptationMatrix[3][3]={{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
    double cosRotationAngle,sinRotationAngle;

    cosRotationAngle=cos(RotationAngle);
    sinRotationAngle=sin(RotationAngle);

    for (i=0;i<3;i++) testRoptationMatrix[i][i]+=cosRotationAngle;
    for (i=0;i<3;i++) for (j=0;j<3;j++) testRoptationMatrix[i][j]+=(1.0-cosRotationAngle)*RotationVector[i]*RotationVector[j]/pow(t,2);

    testRoptationMatrix[0][1]-=sinRotationAngle*RotationVector[2]/t,testRoptationMatrix[0][2]+=sinRotationAngle*RotationVector[1]/t;
    testRoptationMatrix[1][0]+=sinRotationAngle*RotationVector[2]/t,testRoptationMatrix[1][2]-=sinRotationAngle*RotationVector[0]/t;
    testRoptationMatrix[2][0]-=sinRotationAngle*RotationVector[1]/t,testRoptationMatrix[2][1]+=sinRotationAngle*RotationVector[0]/t;


    //CALCULATE THE EFECT OF THE TRANSFORMATION AND TRANSTER THE TRANSFORMED VEWCTROS TO HCI
    double T3[3][3];

    sxform_c("MSGR_MSO","MSGR_HCI",Mercury::OrbitalMotion::et-PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_],T2);

    for (i=0;i<3;i++) for (j=0;j<3;j++) {
      T3[i][j]=0.0;

      for (k=0;k<3;k++) T3[i][j]+=T2[i][k]*testRoptationMatrix[k][j];
    }


    //GET THE ROTATION MATRIX FROM SPICE
//    SpiceDouble rVect1[3],rVect2[3],rVect[3],rot[3][3],xform[6][6];

//    sxform_c (  "MSGR_HCI","MSGR_MSO", et, tsipm ) ;

    //RECALCULATE THE MATRIX AGIN
    double newRotationVector[3],newRate;

    newRate=Exosphere::OrbitalMotion::FrameRotation::GetRotationVector(newRotationVector,"MSGR_MSO",Mercury::OrbitalMotion::et-PIC::ParticleWeightTimeStep::GlobalTimeStep[_NA_SPEC_],Mercury::OrbitalMotion::et);


    //END OF TRANSFORMATIUON TESTS ------------------

    rBegin=sqrt(rBegin);
    rEnd=sqrt(rEnd);

    for (idim=0;idim<3;idim++) {
      lBegin[idim]=StateBegin[idim]/rBegin;
      lEnd[idim]=StateEnd[idim]/rEnd;

      c0+=StateBegin[3+idim]*lBegin[idim];
      c1+=StateEnd[3+idim]*lEnd[idim];
    }

    Mercury::xObjectRadial=0.5*(rBegin+rEnd);
    Mercury::vObjectRadial=0.5*(c0+c1);

    //calculate TAA
    Mercury::OrbitalMotion::TAA=Mercury::OrbitalMotion::GetTAA(Mercury::OrbitalMotion::et);

    for (idim=0;idim<3;idim++) {
      vTangentialBegin+=pow(StateBegin[3+idim]-c0*lBegin[idim],2);
      vTangentialEnd+=pow(StateEnd[3+idim]-c1*lEnd[idim],2);
    }

    vTangentialBegin=sqrt(vTangentialBegin);
    vTangentialEnd=sqrt(vTangentialEnd);

    Mercury::OrbitalMotion::CoordinateFrameRotationRate=0.5*(vTangentialBegin/rBegin+vTangentialEnd/rEnd);


    //determine direction to the Sun and rotation angle in the coordiname frame related to Mercury
    SpiceDouble state[6],l=0.0;

    spkezr_c("SUN",Mercury::OrbitalMotion::et,"IAU_MERCURY","none","MERCURY",state,&lt);

    for (idim=0;idim<3;idim++) l+=pow(state[idim],2);

    for (l=sqrt(l),idim=0;idim<3;idim++) {
      Mercury::OrbitalMotion::SunDirection_IAU_OBJECT[idim]=state[idim]/l;
    }

    //matrixes for tranformation MSO->IAU and IAU->MSO coordinate frames
    sxform_c("MSGR_MSO","IAU_MERCURY",Mercury::OrbitalMotion::et,Mercury::OrbitalMotion::SO_to_IAU_TransformationMartix);
    sxform_c("IAU_MERCURY","MSGR_MSO",Mercury::OrbitalMotion::et,Mercury::OrbitalMotion::IAU_to_SO_TransformationMartix);




    //make the time advance
     PIC::TimeStep();

     //increase the sample length
     static int LastDataOutputFileNumber=-1;

     if ((PIC::DataOutputFileNumber!=0)&&(PIC::DataOutputFileNumber!=LastDataOutputFileNumber)) {
       PIC::RequiredSampleLength*=2;
       if (PIC::RequiredSampleLength>20000) PIC::RequiredSampleLength=20000;

       LastDataOutputFileNumber=PIC::DataOutputFileNumber;
       if (PIC::Mesh::mesh->ThisThread==0) cout << "The new sample length is " << PIC::RequiredSampleLength << endl;
     }
}




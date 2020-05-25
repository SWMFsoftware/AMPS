/*
 * inner_boundary_injection.cpp
 *
 *  Created on: May 16, 2020
 *      Author: vtenishe
 */

#include "sep.h"


double SEP::ParticleSource::InnerBoundary::sphereInjectionRate(int spec,int BoundaryElementType,void *BoundaryElement) {

  double res=1.0E20;

  return res;
}

long int SEP::ParticleSource::InnerBoundary::sphereParticleInjection(int spec,int BoundaryElementType,void *SphereDataPointer) {
  cInternalSphericalData *Sphere;
  double ParticleWeight,LocalTimeStep,/*ExternalNormal[3],*/x[3],v[3],/*r,*/*sphereX0,sphereRadius;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode=NULL;
  long int newParticle,nInjectedParticles=0;
  PIC::ParticleBuffer::byte *newParticleData;
//  int idim;

  double ParticleWeightCorrection=1.0;

//  static const double Temp=200.0;
//  double vbulk[3]={0.0,0.0,0.0};


//  return 0;

//====================  DEBUG ===========================
//  static bool FirstPArticleGenerated=false;
//====================  END DEBUG ===================================


  Sphere=(cInternalSphericalData*)SphereDataPointer;
  Sphere->GetSphereGeometricalParameters(sphereX0,sphereRadius);

#if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
  ParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec];
#else
  exit(__LINE__,__FILE__,"Error: the weight mode is node defined");
#endif


#if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  LocalTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
#elif _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
  LocalTimeStep=Sphere->maxIntersectedNodeTimeStep[spec];
#else
  exit(__LINE__,__FILE__,"Error: the time step node is not defined");
#endif

  
  double TimeCounter=0.0;
  double ModelParticlesInjectionRate=sphereInjectionRate(spec,BoundaryElementType,SphereDataPointer)/ParticleWeight;
  int idim;

  double vbulk[3]={0.0,0.0,0.0};
  double Temp=10.0E6; 
  double ExternalNormal[3];


  double emin=0.1*MeV2J;
  double emax=100.0*MeV2J;

  double s=4.0;
  double q=3.0*s/(3-1.0);

  double p,pmin,pmax,speed,pvect[3]; 

  double cMin=pow(pmin,-q);

  speed=Relativistic::E2Speed(emin,PIC::MolecularData::GetMass(spec));  
  pmin=Relativistic::Speed2Momentum(speed,PIC::MolecularData::GetMass(spec)); 

  speed=Relativistic::E2Speed(emax,PIC::MolecularData::GetMass(spec));
  pmax=Relativistic::Speed2Momentum(speed,PIC::MolecularData::GetMass(spec));

  double A0=pow(pmin,-q+1.0);
  double A=pow(pmax,-q+1.0)-A0; 

  double WeightNorm=pow(pmin,-q);

  int iFieldLine;
  
  while ((TimeCounter+=-log(rnd())/ModelParticlesInjectionRate)<LocalTimeStep) {

    switch (SEP::ParticleTrajectoryCalculation) {
    case SEP::ParticleTrajectoryCalculation_RelativisticBoris:
      Vector3D::Distribution::Uniform(x,sphereRadius);   
      break;
    case ParticleTrajectoryCalculation_FieldLine: case  ParticleTrajectoryCalculation_IgorFieldLine:
      iFieldLine=(int)(PIC::FieldLine::nFieldLine*rnd());
      PIC::FieldLine::FieldLinesAll[iFieldLine].GetFirstVertex()->GetX(x);
      break;
    default:
      exit(__LINE__,__FILE__,"Error: the option is unknown");
    }
    
    startNode=PIC::Mesh::mesh.findTreeNode(x,startNode);

    if (startNode->Thread!=PIC::Mesh::mesh.ThisThread) continue; 
    if (startNode->block->GetLocalTimeStep(spec)/LocalTimeStep<rnd()) continue;

    //generate the particle velocity
    for (idim=0;idim<3;idim++) ExternalNormal[idim]=-x[idim]/sphereRadius; 
    PIC::Distribution::InjectMaxwellianDistribution(v,vbulk,Temp,ExternalNormal,spec);

//    p=pow(A0+rnd()*A,1.0/(-q+1.0));
//    speed=Relativistic::Momentum2Speed(p,PIC::MolecularData::GetMass(spec));

    p=pmin+rnd()*(pmax-pmin);
    ParticleWeightCorrection=pow(p,-q)/WeightNorm;
    speed=Relativistic::Momentum2Speed(p,PIC::MolecularData::GetMass(spec));
    nInjectedParticles++;

    switch (SEP::ParticleTrajectoryCalculation) {
    case SEP::ParticleTrajectoryCalculation_RelativisticBoris:
      Vector3D::Distribution::Uniform(v,speed);
      if (Vector3D::DotProduct(x,v)<0.0) for (int i=0;i<3;i++) v[i]=-v[i];

      //generate a particle
      newParticle=PIC::ParticleBuffer::GetNewParticle();
      newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);

      PIC::ParticleBuffer::SetX(x,newParticleData);
      PIC::ParticleBuffer::SetV(v,newParticleData);
      PIC::ParticleBuffer::SetI(spec,newParticleData);

      PIC::ParticleBuffer::SetIndividualStatWeightCorrection(ParticleWeightCorrection,newParticleData);
      break;
    case ParticleTrajectoryCalculation_FieldLine: case  ParticleTrajectoryCalculation_IgorFieldLine:
      Vector3D::Distribution::Uniform(pvect,p);
      if (Vector3D::DotProduct(x,pvect)<0.0) for (int i=0;i<3;i++) pvect[i]=-pvect[i];

      //add a new aprticle inthe system
      PIC::FieldLine::InjectParticle_default(spec,pvect,ParticleWeightCorrection,iFieldLine,0);
      break;
    default:
      exit(__LINE__,__FILE__,"Error: the option is unknown");
    } 


    //inject the particle into the system
    _PIC_PARTICLE_MOVER__MOVE_PARTICLE_TIME_STEP_(newParticle,startNode->block->GetLocalTimeStep(spec)*rnd(),startNode);
  }

  return nInjectedParticles;

/*
  static double InjectParticles=0.0;

  bool flag;
  int idim;
  double r,ExternalNormal[3];

  InjectParticles+=sphereInjectionRate(spec,BoundaryElementType,SphereDataPointer)*LocalTimeStep;

  while (InjectParticles>0.0) {

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

    InjectParticles-=ParticleWeight;

    startNode=PIC::Mesh::mesh.findTreeNode(x,startNode);
    if (startNode->Thread!=PIC::Mesh::mesh.ThisThread) continue;

    //generate the particle velocity
//    PIC::Distribution::InjectMaxwellianDistribution(v,vbulk,Temp,ExternalNormal,NA);


    for (idim=0;idim<3;idim++) v[idim]=-ExternalNormal[idim]*3.0e3;


//    InjectParticles-=ParticleWeight*ParticleWeightCorrection;
//    if (flag==false) continue;

//====================  DEBUG ===========================
    {
static double InjectionRadialVelocity=0.0,InjectionTangentionalSpeed=0.0;
static long int nTotalInjectedParticles=0;

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




    if (startNode->block->GetLocalTimeStep(spec)/LocalTimeStep<rnd()) continue;

    //generate a particle
    newParticle=PIC::ParticleBuffer::GetNewParticle();
    newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);
    nInjectedParticles++;

    PIC::ParticleBuffer::SetX(x,newParticleData);
    PIC::ParticleBuffer::SetV(v,newParticleData);
    PIC::ParticleBuffer::SetI(spec,newParticleData);

    PIC::ParticleBuffer::SetIndividualStatWeightCorrection(ParticleWeightCorrection,newParticleData);


    //inject the particle into the system
    _PIC_PARTICLE_MOVER__MOVE_PARTICLE_BOUNDARY_INJECTION_(newParticle,startNode->block->GetLocalTimeStep(spec)*rnd() / *LocalTimeStep-TimeCounter* /,startNode,true);
//    PIC::Mover::MoveParticleBoundaryInjection[NA](newParticle,startNode->block->GetLocalTimeStep(NA)*rnd() / *LocalTimeStep-TimeCounter* /,startNode,true);
  }

  return nInjectedParticles;
*/
}

long int SEP::ParticleSource::InnerBoundary::sphereParticleInjection(int BoundaryElementType,void *BoundaryElement) {
  long int spec,res=0;

  for (spec=0;spec<PIC::nTotalSpecies;spec++) res+=sphereParticleInjection(spec,BoundaryElementType,BoundaryElement);

  return res;
}

long int SEP::ParticleSource::InnerBoundary::sphereParticleInjection(void *SphereDataPointer)  {
  long int res=0.0;
  int spec;

  for (spec=0;spec<PIC::nTotalSpecies;spec++) res+=sphereParticleInjection(spec,SphereDataPointer);
  return res;
}

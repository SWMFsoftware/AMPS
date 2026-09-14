/*
 * inner_boundary_injection.cpp
 *
 *  Created on: May 16, 2020
 *      Author: vtenishe
 */

#include "sep.h"
#include "transport_common.h"


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


  switch (_SIMULATION_TIME_STEP_MODE_) {
  case _SINGLE_GLOBAL_TIME_STEP_:
    LocalTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[0];
    break;
  case _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_:
    LocalTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
    break;
  case   _SPECIES_DEPENDENT_LOCAL_TIME_STEP_:
    LocalTimeStep=Sphere->maxIntersectedNodeTimeStep[spec];
    break;
  default:
    exit(__LINE__,__FILE__,"Error: the time step node is not defined");
  }

  double TimeCounter=0.0;
  double ModelParticlesInjectionRate=sphereInjectionRate(spec,BoundaryElementType,SphereDataPointer)/ParticleWeight;
  int idim;

  double vbulk[3]={0.0,0.0,0.0};
  double Temp=10.0E6;
  double ExternalNormal[3];


  double emin=0.1*MeV2J;
  double emax=100.0*MeV2J;

  double s=4.0;
  double q=3.0*s/(s-1.0);

  double p,pmin,pmax,speed,pvect[3];
  double mass=PIC::MolecularData::GetMass(spec);

  pmin=Relativistic::Energy2Momentum(emin,mass);
  pmax=Relativistic::Energy2Momentum(emax,mass);

  double cMin=pow(pmin,-q);

  speed=Relativistic::E2Speed(emin,PIC::MolecularData::GetMass(spec));
  pmin=Relativistic::Speed2Momentum(speed,mass);

  speed=Relativistic::E2Speed(emax,PIC::MolecularData::GetMass(spec));
  pmax=Relativistic::Speed2Momentum(speed,mass);

  double A0=pow(pmin,-q+1.0);
  double A=pow(pmax,-q+1.0)-A0;

  double WeightNorm=pow(pmin,-q);

  int iFieldLine;

  while ((TimeCounter+=-log(rnd())/ModelParticlesInjectionRate)<LocalTimeStep) {

    // Inner-boundary particles are born on an existing magnetic field line.
    // Its first vertex supplies a 3-D position for ownership checks, but that
    // position is not an independent Cartesian particle degree of freedom.
    iFieldLine=(int)(PIC::FieldLine::nFieldLine*rnd());
    PIC::FieldLine::FieldLinesAll[iFieldLine].GetFirstVertex()->GetX(x);

    startNode=PIC::Mesh::mesh->findTreeNode(x,startNode);

    if (startNode->Thread!=PIC::Mesh::mesh->ThisThread) continue;
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

    Vector3D::Distribution::Uniform(pvect,p);
    if (Vector3D::DotProduct(x,pvect)<0.0) {
      for (int i=0;i<3;i++) pvect[i]=-pvect[i];
    }

    // Inject through the field-line API so line id, coordinate, momentum
    // components, and segment-list membership are initialized as one contract.
    newParticle=PIC::FieldLine::InjectParticle_default(
        spec,pvect,ParticleWeightCorrection,iFieldLine,0);
    if (newParticle != -1)
      SEP::Transport::PICAdapter::InitializeParticleTransportState(
          newParticle, UINT64_C(3),
          (static_cast<std::uint64_t>(iFieldLine) << 32) |
              static_cast<std::uint64_t>(nInjectedParticles),
          pvect);
    newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);


    //inject the particle into the system
    int code;

    code=_PIC_PARTICLE_MOVER__MOVE_PARTICLE_TIME_STEP_(newParticle,startNode->block->GetLocalTimeStep(spec)*rnd(),startNode);

    //apply condition of tracking the particle
    if ((_PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_)&&(code==_PARTICLE_MOTION_FINISHED_)) {
      PIC::ParticleTracker::InitParticleID(newParticleData);
      PIC::ParticleTracker::ApplyTrajectoryTrackingCondition(x,v,spec,newParticleData,(void*)startNode);
    }
  }

  return nInjectedParticles;

}

long int SEP::ParticleSource::InnerBoundary::sphereParticleInjection(int BoundaryElementType,void *BoundaryElement) {
  (void)BoundaryElementType;
  (void)BoundaryElement;

  // AMR-surface injection existed only for Cartesian movers.  Field-line
  // sources use SEP::FieldLine::InjectParticles, registered in main_lib.cpp.
  return 0;
}

long int SEP::ParticleSource::InnerBoundary::sphereParticleInjection(void *SphereDataPointer)  {
  (void)SphereDataPointer;

  // Retain the AMPS callback signature for linkage, but do not manufacture a
  // Cartesian particle record in the field-line-only srcSEP application.
  return 0;
}

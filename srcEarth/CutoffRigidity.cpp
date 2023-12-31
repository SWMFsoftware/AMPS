
//$Id$
//functionality for calculating of the cutogg rigidity

#include "pic.h"
#include "Earth.h"
#include "specfunc.h"

/*
 * CutoffRigidity.cpp
 *
 *  Created on: Dec 25, 2016
 *      Author: vtenishe
 */

bool Earth::CutoffRigidity::SampleRigidityMode=false;
long int Earth::CutoffRigidity::InitialRigidityOffset=-1;
long int Earth::CutoffRigidity::InitialLocationOffset=-1;
long int Earth::CutoffRigidity::IntegratedPathLengthOffset=-1;
long int Earth::CutoffRigidity::IntegratedTimeOffset=-1;
array_2d<double> Earth::CutoffRigidity::CutoffRigidityTable;
array_2d<double> Earth::CutoffRigidity::SampledFluxTable;
array_2d<int> Earth::CutoffRigidity::InjectedParticleMap;
array_2d<double> Earth::CutoffRigidity::MaxEnergyInjectedParticles;

//enable/disable the particle injection procedure
bool Earth::CutoffRigidity::ParticleInjector::ParticleInjectionMode=true;

//the total number of iteration used for calcualting the cutoff rigidity
int Earth::CutoffRigidity::nTotalIterations=-1;

//the number of the locations where the rigidity and particle flux is simulated
int Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength=0;
double** Earth::CutoffRigidity::IndividualLocations::xTestLocationTable=NULL;
double Earth::CutoffRigidity::IndividualLocations::MaxEnergyLimit=0.0;
double Earth::CutoffRigidity::IndividualLocations::MinEnergyLimit=0.0;
array_2d<double>  Earth::CutoffRigidity::IndividualLocations::CutoffRigidityTable;
array_2d<double>  Earth::CutoffRigidity::IndividualLocations::SampledFluxTable;

double Earth::CutoffRigidity::IndividualLocations::MinInjectionRigidityLimit=0.0;
double Earth::CutoffRigidity::IndividualLocations::MaxInjectionRigidityLimit=0.0;
int Earth::CutoffRigidity::IndividualLocations::nRigiditySearchIntervals=1;
int Earth::CutoffRigidity::IndividualLocations::InjectionMode=Earth::CutoffRigidity::IndividualLocations::_energy_injection; 

double CutoffRigidityTestLocationTable[][3]={{0.0,0.0,0.0}};
int CutoffRigidityTestLocationTableLength=0;

long int Earth::CutoffRigidity::ParticleDataOffset::OriginLocationIndex=-1;
long int Earth::CutoffRigidity::ParticleDataOffset::OriginalSpeed=-1;
long int Earth::CutoffRigidity::ParticleDataOffset::OriginalVelocityDirectionID=-1;


int Earth::CutoffRigidity::IndividualLocations::nTotalTestParticlesPerLocations=0;  //the total number of model particles ejected from a test location
int Earth::CutoffRigidity::IndividualLocations::nParticleInjectionIterations=1;

//sample velocity dependent rigidity cutoff
cInternalSphericalData *Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable=NULL;

void Earth::CutoffRigidity::Init_BeforeParser() {
  if ((SampleRigidityMode==true)&&(InitialRigidityOffset==-1)) {
    //request data for sampling of the cutoff rigidity in the particle state vector
    PIC::ParticleBuffer::RequestDataStorage(InitialRigidityOffset,sizeof(double));
    PIC::ParticleBuffer::RequestDataStorage(InitialLocationOffset,3*sizeof(double));
    PIC::ParticleBuffer::RequestDataStorage(IntegratedPathLengthOffset,sizeof(double));
    PIC::ParticleBuffer::RequestDataStorage(IntegratedTimeOffset,sizeof(double));
  }

  //allocate 'IndividualLocations::xTestLocationTable'
  if ((CutoffRigidityTestLocationTableLength!=0)||(IndividualLocations::xTestLocationTableLength!=0)) {
    //de allocate previously allocated buffers
    /*
    if (IndividualLocations::xTestLocationTable!=NULL) {
      delete [] IndividualLocations::xTestLocationTable[0];
      delete [] IndividualLocations::xTestLocationTable;
    }
    */


    //allocate new sampling buffers
    if (IndividualLocations::xTestLocationTableLength==0) {
      IndividualLocations::xTestLocationTableLength=CutoffRigidityTestLocationTableLength;

      if (IndividualLocations::xTestLocationTable==NULL) {
        IndividualLocations::xTestLocationTable=new double* [CutoffRigidityTestLocationTableLength];
        IndividualLocations::xTestLocationTable[0]=new double [3*CutoffRigidityTestLocationTableLength];

        for (int i=1;i<CutoffRigidityTestLocationTableLength;i++) IndividualLocations::xTestLocationTable[i]=IndividualLocations::xTestLocationTable[i-1]+3;

        for (int i=0;i<IndividualLocations::xTestLocationTableLength;i++) {
          for (int j=0;j<3;j++) IndividualLocations::xTestLocationTable[i][j]=CutoffRigidityTestLocationTable[i][j];
        }
      }
    }

    if (IndividualLocations::CutoffRigidityTable.IsAllocated()==false) {
      IndividualLocations::CutoffRigidityTable.init(PIC::nTotalSpecies,IndividualLocations::xTestLocationTableLength);
      IndividualLocations::SampledFluxTable.init(PIC::nTotalSpecies,IndividualLocations::xTestLocationTableLength);
    }

    IndividualLocations::CutoffRigidityTable=-1.0;
    IndividualLocations::SampledFluxTable=0.0;


    //request particle storage if needed
    if (Earth::CutoffRigidity::ParticleDataOffset::OriginLocationIndex==-1) {
      PIC::ParticleBuffer::RequestDataStorage(Earth::CutoffRigidity::ParticleDataOffset::OriginLocationIndex,sizeof(int));
      PIC::ParticleBuffer::RequestDataStorage(Earth::CutoffRigidity::ParticleDataOffset::OriginalSpeed,sizeof(double));
      PIC::ParticleBuffer::RequestDataStorage(Earth::CutoffRigidity::ParticleDataOffset::OriginalVelocityDirectionID,sizeof(int));
    }
  }
}

void Earth::CutoffRigidity::AllocateCutoffRigidityTable() {
  int i,j,k,offset;

  if ((CutoffRigidityTestLocationTableLength!=0)&&(CutoffRigidityTable.IsAllocated()==false)) {
    //allocate the cutoff rigidity table
    //access pattern CutoffRigidityTable[spec][iZenith][iAzimuthal]
    CutoffRigidityTable.init(PIC::nTotalSpecies,Earth::Planet->nZenithSurfaceElements*Earth::Planet->nAzimuthalSurfaceElements);
    SampledFluxTable.init(PIC::nTotalSpecies,Earth::Planet->nZenithSurfaceElements*Earth::Planet->nAzimuthalSurfaceElements);

    CutoffRigidityTable=-1.0;
    SampledFluxTable=0.0;
  }
}


//process model particles that leaves the computational domain
int Earth::CutoffRigidity::ProcessOutsideDomainParticles(long int ptr,double* xInit,double* vInit,int nIntersectionFace,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  double *x,Rigidity,v[3],ParticleEnergy;
  long int iAzimuth,iZenith;
  int spec;
  PIC::ParticleBuffer::byte *ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);

  if (Earth::ModelMode!=Earth::CutoffRigidityMode) {
    return _PARTICLE_DELETED_ON_THE_FACE_;   
  }

  spec=PIC::ParticleBuffer::GetI(ParticleData);
  PIC::ParticleBuffer::GetV(v,ParticleData);
  ParticleEnergy=PIC::MolecularData::GetMass(spec)*Vector3D::DotProduct(v,v)/2.0;

  //register the particle velocity vector
  int iOriginIndex=0;

  iOriginIndex=*((int*)(ParticleData+Earth::CutoffRigidity::ParticleDataOffset::OriginLocationIndex));
  if ((DomainBoundaryParticleProperty::SampleDomainBoundaryParticleProperty==true)&&(nIntersectionFace!=-1)) DomainBoundaryParticleProperty::RegisterParticleProperties(PIC::ParticleBuffer::GetI(ptr),xInit,vInit,iOriginIndex,nIntersectionFace);

  int VelocityDirectionId=-1;
  double charge,rigidity,momentum,energy;

  if (Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable!=NULL) {
    VelocityDirectionId=*((int*)(ParticleData+Earth::CutoffRigidity::ParticleDataOffset::OriginalVelocityDirectionID));
    rigidity=*((double*)(ParticleData+InitialRigidityOffset));

//if (VelocityDirectionId==22499) {
//DebuggerTrap();
//}

    charge=fabs(PIC::MolecularData::GetElectricCharge(spec));
    momentum=rigidity*1.0E9*charge/SpeedOfLight;
    energy=Relativistic::Momentum2Speed(momentum,PIC::MolecularData::GetMass(spec)); 

    if (Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable[iOriginIndex].minRigidity[spec][VelocityDirectionId]<0.0) {
      Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable[iOriginIndex].minRigidity[spec][VelocityDirectionId]=rigidity;
    }
    else if (Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable[iOriginIndex].minRigidity[spec][VelocityDirectionId]>rigidity) {
      Earth::CutoffRigidity::IndividualLocations::SamplingSphereTable[iOriginIndex].minRigidity[spec][VelocityDirectionId]=rigidity;
    }  
  }

  //update the rigidity data
  if (SampleRigidityMode==true) {
    x=(double*)(ParticleData+InitialLocationOffset);

    Rigidity=*((double*)(ParticleData+InitialRigidityOffset));

    //get coordinates of the point of origin of the particle
    if ((Earth::Planet!=NULL)&&(CutoffRigidityTable.IsAllocated()==true)) {
      Earth::Planet->GetSurfaceElementProjectionIndex(x,iZenith,iAzimuth);

      if ((iZenith<0)||(iZenith>=PIC::nTotalSpecies*Earth::Planet->nZenithSurfaceElements)||(iAzimuth<0)||(iAzimuth>=Earth::Planet->nAzimuthalSurfaceElements)) exit(__LINE__,__FILE__,"Error: out of range");

/*
      double *RigidityTableElement=CutoffRigidityTable.GetPtr(spec,iOriginIndex);

      #pragma omp critical
      if ((*RigidityTableElement<0.0)||(*RigidityTableElement>Rigidity)) *RigidityTableElement=Rigidity;
*/

      if ((CutoffRigidityTable(spec,iOriginIndex)<0.0)||(CutoffRigidityTable(spec,iOriginIndex)>Rigidity)) CutoffRigidityTable(spec,iOriginIndex)=Rigidity; 
      SampledFluxTable(spec,iOriginIndex)+=Earth::OutsideParticleFlux(ParticleEnergy);

    }
  }

  //cutoff rigidity at individual locations
  if ((Earth::CutoffRigidity::IndividualLocations::xTestLocationTableLength!=0)&&(IndividualLocations::CutoffRigidityTable.IsAllocated()==true)) {
    if (SampleRigidityMode==true) {
      //sample the cutoff rigidity
      double *RigidityTableElement=IndividualLocations::CutoffRigidityTable.GetPtr(spec,iOriginIndex);

      if ((*RigidityTableElement<0.0)||(Rigidity<*RigidityTableElement)) *RigidityTableElement=Rigidity;
      IndividualLocations::SampledFluxTable(spec,iOriginIndex)+=Earth::OutsideParticleFlux(ParticleEnergy);
    }
  }


  return _PARTICLE_DELETED_ON_THE_FACE_;
}

//output sampled cutoff rigidity map
void Earth::CutoffRigidity::OutputDataFile::PrintVariableList(FILE* fout) {
  if (Earth::ModelMode!=Earth::CutoffRigidityMode) return;

  fprintf(fout,", \"Cutoff Rigidity\", \"Injected Particle Number\"");
}

void Earth::CutoffRigidity::OutputDataFile::PrintDataStateVector(FILE* fout,long int nZenithPoint,long int nAzimuthalPoint,long int *SurfaceElementsInterpolationList,long int SurfaceElementsInterpolationListLength,cInternalSphericalData *Sphere,int spec,CMPI_channel* pipe,int ThisThread,int nTotalThreads) {
  int nInterpolationElement,nSurfaceElement,iZenith,iAzimuth;
  double InterpolationNormalization=0.0,InterpolationCoefficient;

  double CutoffRigidity=0.0,SurfaceElementCutoffRigidity;
  double InterpolatedInjectedParticleNumber=0.0,normInterpolatedInjectedParticleNumber=0.0;
  int InjectedParticleNumber;

  if (Earth::ModelMode!=Earth::CutoffRigidityMode) return;

  for (nInterpolationElement=0;nInterpolationElement<SurfaceElementsInterpolationListLength;nInterpolationElement++) {
    nSurfaceElement=SurfaceElementsInterpolationList[nInterpolationElement];
    InterpolationCoefficient=Sphere->GetSurfaceElementArea(nSurfaceElement);

    Sphere->GetSurfaceElementIndex(iZenith,iAzimuth,nSurfaceElement);

    InjectedParticleNumber=InjectedParticleMap(spec,nSurfaceElement); 

    if (PIC::ThisThread!=0) {
      pipe->send(InjectedParticleNumber);
      pipe->send(CutoffRigidityTable(spec,nSurfaceElement));
    }
    else {
      double t;
      int thread;

      SurfaceElementCutoffRigidity=CutoffRigidityTable(spec,nSurfaceElement);

      for (thread=1;thread<PIC::nTotalThreads;thread++) {
        InjectedParticleNumber+=pipe->recv<int>(thread);        

        t=pipe->recv<double>(thread);
        if ((SurfaceElementCutoffRigidity<0.0)||(SurfaceElementCutoffRigidity>t)) SurfaceElementCutoffRigidity=t;
      }

      InterpolatedInjectedParticleNumber+=InjectedParticleNumber*InterpolationCoefficient;
      normInterpolatedInjectedParticleNumber+=InterpolationCoefficient; 

      if (SurfaceElementCutoffRigidity>0.0) {
        CutoffRigidity+=SurfaceElementCutoffRigidity*InterpolationCoefficient;
        InterpolationNormalization+=InterpolationCoefficient;
      }
    }
  }

  if (PIC::ThisThread==0) fprintf(fout," %e  %e ",((InterpolationNormalization>0.0) ? CutoffRigidity/InterpolationNormalization : -1),InterpolatedInjectedParticleNumber/normInterpolatedInjectedParticleNumber);
}

//injection rate of the test particles when calculate the cut-off rigidity
double Earth::CutoffRigidity::ParticleInjector::GetTotalProductionRate(int spec,int BoundaryElementType,void *SphereDataPointer) {return 1.0E14;}

//generate a new particle
bool Earth::CutoffRigidity::ParticleInjector::GenerateParticleProperties(int spec,PIC::ParticleBuffer::byte* tempParticleData,double *x_SO_OBJECT,
  double *x_IAU_OBJECT,double *v_SO_OBJECT,double *v_IAU_OBJECT,double *sphereX0,
  double sphereRadius,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* &startNode,int BoundaryElementType,void *BoundaryElement) {

  int idim;
  double ExternalNormal[3];

  //if the injection model is disabled => exit from the procedure
  if (ParticleInjectionMode==false) return false;

  //Generate new particle position
  Vector3D::Distribution::Uniform(ExternalNormal);

  for (idim=0;idim<DIM;idim++) {
    x_IAU_OBJECT[idim]=-RigidityTestRadiusVector*ExternalNormal[idim];
    x_SO_OBJECT[idim]=x_IAU_OBJECT[idim];
  }

  //determine if the particle belongs to this processor
  startNode=PIC::Mesh::mesh->findTreeNode(x_SO_OBJECT,startNode);
  if (startNode->Thread!=PIC::Mesh::mesh->ThisThread) return false;

  //generate velocity of the injected particle
  //IMPORTANT! The velocity vector is directed inside the sphere; and the trajectory is traced backward

  //direction of the new particle velocity
  do {
    Vector3D::Distribution::Uniform(v_IAU_OBJECT);
  }
  while (Vector3D::DotProduct(v_IAU_OBJECT,ExternalNormal)<=0.0);

  //energy and rigidity of the new particles
  double mass,speed,energy,rigidity,momentum,charge;

  mass=PIC::MolecularData::GetMass(spec);
  charge=fabs(PIC::MolecularData::GetElectricCharge(spec));

  static double logRigidityTestMinEnergy=log(RigidityTestMinEnergy);
  static double logRigidityTestMaxEnergy=log(RigidityTestMaxEnergy);

  energy=exp(logRigidityTestMinEnergy+rnd()*(logRigidityTestMaxEnergy-logRigidityTestMinEnergy));

  speed=Relativistic::E2Speed(energy,mass);
  momentum=Relativistic::Speed2Momentum(speed,mass);

  rigidity=(charge>0.0) ? momentum/charge : 0.0;

  for (idim=0;idim<3;idim++) {
    v_IAU_OBJECT[idim]*=speed;
    v_SO_OBJECT[idim]=v_IAU_OBJECT[idim];
  }

  //save the initial location, and rigidity of the particle
  if (SampleRigidityMode==true) {
    *((double*)(tempParticleData+InitialRigidityOffset))=rigidity;
    memcpy((tempParticleData+InitialLocationOffset),v_IAU_OBJECT,3*sizeof(double));
  }

  return true;
}


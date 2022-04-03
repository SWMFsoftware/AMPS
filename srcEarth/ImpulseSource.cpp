//$Id$
//functionality of the impulse source of the energetic particles in the magnetosphere

/*
 * ImpulseSource.cpp
 *
 *  Created on: Dec 26, 2016
 *      Author: vtenishe
 */

#include "pic.h"
#include "Earth.h"

Earth::ImpulseSource::cImpulseSourceData Earth::ImpulseSource::ImpulseSourceData[]={0,0.0,false,{0.0,0.0,0.0},0.0,0.0};
int Earth::ImpulseSource::nTotalSourceLocations=0;
double Earth::ImpulseSource::TimeCounter=0.0;
double Earth::ImpulseSource::EnergySpectrum::Constant::e=1.0*MeV2J;
int Earth::ImpulseSource::EnergySpectrum::Mode=Earth::ImpulseSource::EnergySpectrum::Mode_Constatant;
bool Earth::ImpulseSource::Mode=false;


void amps_init();
void amps_init_mesh();
void amps_time_step();

//inject the energetic particles
long int Earth::ImpulseSource::InjectParticles() {
  int nTotalInjectedParticles=0;
  int idim,spec,iSource;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode;
  double mass,ElectricCharge,a;
  long int newParticle;
  PIC::ParticleBuffer::byte *newParticleData;

  //sampling buffer of the gyro-frequency and gyro-radii
  double *GyroFrequencySample=new double[PIC::nTotalThreadsOpenMP];
  double *GyroRadiiSample=new double [PIC::nTotalThreadsOpenMP];
  int *SampleCounter=new int [PIC::nTotalThreadsOpenMP];

  for (int thread=0;thread<PIC::nTotalThreadsOpenMP;thread++) {
    SampleCounter[thread]=0;
    GyroFrequencySample[thread]=0.0,GyroRadiiSample[thread]=0.0;
  }

  double SourceLocationB[3];
  TimeCounter+=PIC::ParticleWeightTimeStep::GetGlobalTimeStep(0);

  for (iSource=0;iSource<nTotalSourceLocations;iSource++) {
    if ((TimeCounter>=ImpulseSourceData[iSource].time)&&(ImpulseSourceData[iSource].ProcessedFlag==false)) {
      ImpulseSourceData[iSource].ProcessedFlag=true;

      //inject energetic particles
      startNode=PIC::Mesh::mesh->findTreeNode(ImpulseSourceData[iSource].x);
      if (startNode->Thread!=PIC::Mesh::mesh->ThisThread) continue;

      //determine the number of the injected particles
      spec=ImpulseSourceData[iSource].spec;
      a=ImpulseSourceData[iSource].Source/startNode->block->GetLocalParticleWeight(spec);
      nTotalInjectedParticles=(int)a;
      a-=nTotalInjectedParticles;
      if (rnd()<a) nTotalInjectedParticles++;

      mass=PIC::MolecularData::GetMass(spec);
      ElectricCharge=PIC::MolecularData::GetElectricCharge(spec);

      //get magnetic field at the source location
      PIC::CPLR::InitInterpolationStencil(ImpulseSourceData[iSource].x,startNode);
      PIC::CPLR::GetBackgroundMagneticField(SourceLocationB);

      //reset the sampling buffers
      for (int thread=0;thread<PIC::nTotalThreadsOpenMP;thread++) {
        GyroFrequencySample[thread]=0.0,GyroRadiiSample[thread]=0.0;
        SampleCounter[thread]=0;
      }


      long int *ParticleTable=new long int [nTotalInjectedParticles];

      
      for (int iPart=0;iPart<nTotalInjectedParticles;iPart++) {
        //generate particles' velocity
        ParticleTable[iPart]=PIC::ParticleBuffer::GetNewParticle(true);
      }

       #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
       #pragma omp parallel for firstprivate (newParticle) private (idim,newParticleData)  \
       firstprivate (SourceLocationB,GyroFrequencySample,GyroRadiiSample,SampleCounter,TimeCounter,iSource,nTotalInjectedParticles,startNode,spec,mass,ElectricCharge,EnergySpectrum::Mode,EnergySpectrum::Constant::e,ImpulseSourceData)
       #endif
      for (int iPart=0;iPart<nTotalInjectedParticles;iPart++) {
        //generate particles' velocity
        newParticle=ParticleTable[iPart]; 
//PIC::ParticleBuffer::GetNewParticle(true);

        //generate new particle velocity
        double v[3],speed;

        switch (EnergySpectrum::Mode) {
        case EnergySpectrum::Mode_Constatant:
          speed=Relativistic::E2Speed(EnergySpectrum::Constant::e,mass);
          break;
        default:
          exit(__LINE__,__FILE__,"Error: the option is not found");
        }

        Vector3D::Distribution::Uniform(v);
        for (idim=0;idim<3;idim++) v[idim]*=speed;

        //generate particles' velocity
        newParticleData=PIC::ParticleBuffer::GetParticleDataPointer(newParticle);

        PIC::ParticleBuffer::SetX(ImpulseSourceData[iSource].x,newParticleData);
        PIC::ParticleBuffer::SetV(v,newParticleData);
        PIC::ParticleBuffer::SetI(spec,newParticleData);
        PIC::ParticleBuffer::SetIndividualStatWeightCorrection(1.0,newParticleData);

        //sample the gyrofequency and gyroradius of the source
        #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
        int thread=omp_get_thread_num();
        #else
        int thread=0;
        #endif

        GyroFrequencySample[thread]+=Relativistic::GetGyroFrequency(v,mass,ElectricCharge,SourceLocationB);
        GyroRadiiSample[thread]+=Relativistic::GetGyroRadius(v,mass,ElectricCharge,SourceLocationB);
        SampleCounter[thread]++;

        //inject the particle into the system
        _PIC_PARTICLE_MOVER__MOVE_PARTICLE_TIME_STEP_(newParticle,TimeCounter-ImpulseSourceData[iSource].time,startNode);
      }

      //end of the particle injetion loop

      //output sampled information
      for (int thread=1;thread<PIC::nTotalThreadsOpenMP;thread++) {
        GyroFrequencySample[0]+=GyroFrequencySample[thread];
        GyroRadiiSample[0]+=GyroRadiiSample[thread];
        SampleCounter[0]+=SampleCounter[thread];
      }

      delete [] ParticleTable;

      if (SampleCounter[0]!=0) GyroRadiiSample[0]/=SampleCounter[0],GyroFrequencySample[0]/=SampleCounter[0];
      printf("$PREFIX: Impulse source location %i:  spec=%i\n",iSource,spec);
      printf("$PREFIX: Impulse source location %i:  Mean Gyro Radii=%e\n",iSource,GyroRadiiSample[0]);
      printf("$PREFIX: Impulse source location %i:  Mean Gyro Frequency=%e\n",iSource,GyroFrequencySample[0]);
    }

    //increment the species dependent time counter

  }


  //deallocate the sampling buffers
  delete [] GyroFrequencySample;
  delete [] GyroRadiiSample;
  delete [] SampleCounter;

  return nTotalInjectedParticles;
}

//set weights of the model particles
void Earth::ImpulseSource::InitParticleWeight() {
  int spec,iTotalSourceLocations;
  double WeightTable[PIC::nTotalSpecies],TotalSourceTable[PIC::nTotalSpecies],TotalRequestedModelParticlesTable[PIC::nTotalSpecies];
  double MinParticleWeight=-1.0;

  for (spec=0;spec<PIC::nTotalSpecies;spec++) WeightTable[spec]=-1.0,TotalSourceTable[spec]=0.0,TotalRequestedModelParticlesTable[spec]=0.0;

  //loop through all injection events and collect the injection data
  for (iTotalSourceLocations=0;iTotalSourceLocations<nTotalSourceLocations;iTotalSourceLocations++) {
    TotalSourceTable[ImpulseSourceData[iTotalSourceLocations].spec]+=ImpulseSourceData[iTotalSourceLocations].Source;
    TotalRequestedModelParticlesTable[ImpulseSourceData[iTotalSourceLocations].spec]+=ImpulseSourceData[iTotalSourceLocations].NumberInjectedParticles;
  }

  //calculate the weight
  for (spec=0;spec<PIC::nTotalSpecies;spec++) if (TotalRequestedModelParticlesTable[spec]>0.0) {
    WeightTable[spec]=TotalSourceTable[spec]/TotalRequestedModelParticlesTable[spec];
    if ((MinParticleWeight<0.0)||(WeightTable[spec]<MinParticleWeight)) MinParticleWeight=WeightTable[spec];
  }

  //1) check is weights os all model spaces are defined. if not -> used the min weight
  //2) set the particle weight in the system
  for (spec=0;spec<PIC::nTotalSpecies;spec++) {
    if (WeightTable[spec]<0.0) WeightTable[spec]=MinParticleWeight;
    PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(spec,WeightTable[spec]);
  }
}

//run the impulse source model
void RunImpulseSource() {
  Earth::ImpulseSource::Mode==true; 

  amps_init_mesh();
  amps_init();


  int nTotalIterations=1000000;
  int LastDataOutputFileNumber=0;

  if (_PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_) {
    nTotalIterations=100;
    PIC::RequiredSampleLength=50;
  } 
  
  for (long int niter=0;niter<nTotalIterations;niter++) {
    amps_time_step();

    if (PIC::Mesh::mesh->ThisThread==0) {
      time_t TimeValue=time(NULL);
      tm *ct=localtime(&TimeValue);
      printf(": (%i/%i %i:%i:%i), Iteration: %ld  (current sample length:%ld, %ld interations to the next output)\n",
             ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec,niter,
             PIC::RequiredSampleLength,
             PIC::RequiredSampleLength-PIC::CollectingSampleCounter);
    }

     if ((PIC::DataOutputFileNumber!=0)&&(PIC::DataOutputFileNumber!=LastDataOutputFileNumber)) {
//       PIC::RequiredSampleLength*=2;
       if (PIC::RequiredSampleLength>50000) PIC::RequiredSampleLength=50000;


       LastDataOutputFileNumber=PIC::DataOutputFileNumber;
       if (PIC::Mesh::mesh->ThisThread==0) cout << "The new sample length is " << PIC::RequiredSampleLength << endl;
     }
  }

  char fname[400];

  sprintf(fname,"%s/test_Earth-Impulse-Neutron.dat",PIC::OutputDataFileDirectory);
  PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname);


}

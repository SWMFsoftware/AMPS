//====================================================
//$Id$
//====================================================
//the general functions for the pic solver

#include "pic.h"


//sampling variables

long int PIC::LastSampleLength=0,PIC::CollectingSampleCounter=0,PIC::RequiredSampleLength=100,PIC::DataOutputFileNumber=0;
int PIC::SamplingMode=_RESTART_SAMPLING_MODE_;

//====================================================
//perform one time step
void PIC::TimeStep() {
   double ParticleExchangeTime,IterationExecutionTime,SamplingTime,StartTime=MPI_Wtime();

  //sampling of the particle data
  PIC::Sampling();
  SamplingTime=MPI_Wtime()-StartTime;

  //move existing particles
  PIC::Mover::MoveParticles();

  //injection boundary conditions
  PIC::BC::InjectionBoundaryConditions();

  IterationExecutionTime=MPI_Wtime()-StartTime;

  //syncrinie processors and exchnge particle data
  ParticleExchangeTime=MPI_Wtime();
  PIC::Parallel::ExchangeParticleData();
  ParticleExchangeTime=MPI_Wtime()-ParticleExchangeTime;


  //Collect and exchange the run's statictic information
  static const int nRunStatisticExchangeIterationsMin=5,nRunStatisticExchangeIterationsMax=500,nRunStatisticExchangeTime=120;
  static long int nTotalIterations=0,nInteractionsAfterRunStatisticExchange=0;
  static int nExchangeStatisticsIterationNumberSteps=10;

  struct cExchangeStatisticData {
    double TotalInterationRunTime;
    double IterationExecutionTime;
    long int TotalParticlesNumber;
    double ParticleExchangeTime;
    double SamplingTime;
    double Latency;
    long int recvParticleCounter,sendParticleCounter;
    long int nInjectedParticles;
  };

  nTotalIterations++;
  nInteractionsAfterRunStatisticExchange++;

  if (nInteractionsAfterRunStatisticExchange==nExchangeStatisticsIterationNumberSteps) { //collect and exchenge the statistical data of the run
    cExchangeStatisticData localRunStatisticData;
    int thread;
    double Latency=MPI_Wtime();
    cExchangeStatisticData *ExchangeBuffer=NULL;

    nInteractionsAfterRunStatisticExchange=0;

    MPI_Barrier(MPI_COMM_WORLD);
    Latency=MPI_Wtime()-Latency;


    localRunStatisticData.TotalInterationRunTime=MPI_Wtime()-StartTime;
    localRunStatisticData.IterationExecutionTime=IterationExecutionTime;
    localRunStatisticData.TotalParticlesNumber=PIC::ParticleBuffer::NAllPart;
    localRunStatisticData.ParticleExchangeTime=ParticleExchangeTime;
    localRunStatisticData.SamplingTime=SamplingTime;
    localRunStatisticData.Latency=Latency;
    localRunStatisticData.recvParticleCounter=PIC::Parallel::recvParticleCounter;
    localRunStatisticData.sendParticleCounter=PIC::Parallel::sendParticleCounter;
    localRunStatisticData.nInjectedParticles=PIC::BC::nInjectedParticles;

    if (PIC::Mesh::mesh.ThisThread==0) {
      time_t TimeValue=time(NULL);
      tm *ct=localtime(&TimeValue);

      printf("\nPIC: (%i/%i %i:%i:%i), Iteration: %ld  (%ld interations to the next output)\n",ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec,nTotalIterations,PIC::RequiredSampleLength-PIC::CollectingSampleCounter);



      ExchangeBuffer=new cExchangeStatisticData[PIC::Mesh::mesh.nTotalThreads];
      MPI_Gather((char*)&localRunStatisticData,sizeof(cExchangeStatisticData),MPI_CHAR,(char*)ExchangeBuffer,sizeof(cExchangeStatisticData),MPI_CHAR,0,MPI_COMM_WORLD);


      //output the data
      long int nTotalModelParticles=0,nTotalInjectedParticels=0;
      double MinExecutionTime=localRunStatisticData.IterationExecutionTime,MaxExecutionTime=localRunStatisticData.IterationExecutionTime,MaxLatency=0.0;

      printf("Thread, Total Particle's number, Total Interation Time, Iteration Execution Time, Particle Exchange Time, Latency, Send Particles, Recv Particles, nInjected Particls\n");

      for (thread=0;thread<PIC::Mesh::mesh.nTotalThreads;thread++) {
        printf("%i\t %ld\t %e\t %e\t %e\t %e\t %ld\t %ld\t %ld\n",thread,ExchangeBuffer[thread].TotalParticlesNumber,ExchangeBuffer[thread].TotalInterationRunTime,
            ExchangeBuffer[thread].IterationExecutionTime,ExchangeBuffer[thread].ParticleExchangeTime,ExchangeBuffer[thread].Latency,ExchangeBuffer[thread].sendParticleCounter,
            ExchangeBuffer[thread].recvParticleCounter,ExchangeBuffer[thread].nInjectedParticles);

        nTotalModelParticles+=ExchangeBuffer[thread].TotalParticlesNumber;
        nTotalInjectedParticels+=ExchangeBuffer[thread].nInjectedParticles;


        if (MinExecutionTime>ExchangeBuffer[thread].IterationExecutionTime) MinExecutionTime=ExchangeBuffer[thread].IterationExecutionTime;
        if (MaxExecutionTime<ExchangeBuffer[thread].IterationExecutionTime) MaxExecutionTime=ExchangeBuffer[thread].IterationExecutionTime;
        if (MaxLatency<ExchangeBuffer[thread].Latency) MaxLatency=ExchangeBuffer[thread].Latency;
      }

      printf("Total number of particles: %ld\n",nTotalModelParticles);
      printf("Total number of injected particles: %ld\n",nTotalInjectedParticels);
      printf("Iteration Execution Time: min=%e, max=%e\n",MinExecutionTime,MaxExecutionTime);
      printf("Latency: max=%e\n",MaxLatency);

      //determine the new value for the interation number between exchanging of the run stat data
      if (localRunStatisticData.TotalInterationRunTime>0.0) {
        nExchangeStatisticsIterationNumberSteps=(int)(nRunStatisticExchangeTime/localRunStatisticData.TotalInterationRunTime);
        if (nExchangeStatisticsIterationNumberSteps<nRunStatisticExchangeIterationsMin) nExchangeStatisticsIterationNumberSteps=nRunStatisticExchangeIterationsMin;
        if (nExchangeStatisticsIterationNumberSteps>nRunStatisticExchangeIterationsMax) nExchangeStatisticsIterationNumberSteps=nRunStatisticExchangeIterationsMax;
      }
      else nExchangeStatisticsIterationNumberSteps=nRunStatisticExchangeIterationsMin;

      MPI_Bcast(&nExchangeStatisticsIterationNumberSteps,1,MPI_INT,0,MPI_COMM_WORLD);


      delete [] ExchangeBuffer;
    }
    else {
      MPI_Gather((char*)&localRunStatisticData,sizeof(cExchangeStatisticData),MPI_CHAR,(char*)ExchangeBuffer,sizeof(cExchangeStatisticData),MPI_CHAR,0,MPI_COMM_WORLD);
      MPI_Bcast(&nExchangeStatisticsIterationNumberSteps,1,MPI_INT,0,MPI_COMM_WORLD);
    }


    MPI_Barrier(MPI_COMM_WORLD);
  }

}
//====================================================
//the general sampling procedure
void PIC::Sampling() {
  int s,i,j,k,idim;
  long int LocalCellNumber,ptr;

  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];
  PIC::ParticleBuffer::byte *ParticleData;
  PIC::Mesh::cDataCenterNode *cell;
  PIC::Mesh::cDataBlockAMR *block;
  char *SamplingData;
  double *v,LocalParticleWeight;

  //parallel efficientcy measure
#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_PARTICLE_NUMBER_
  long int TreeNodeTotalParticleNumber;
#elif _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
  double TreeNodeProcessingTime;
#endif

  //go through the 'local nodes'
  while (node!=NULL) {
    block=node->block;

#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_PARTICLE_NUMBER_
    TreeNodeTotalParticleNumber=0;
#elif _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
    TreeNodeProcessingTime=MPI_Wtime();
#endif

    //sample the distrivution function
#if _SAMPLING_DISTRIBUTION_FUNCTION_MODE_ == _SAMPLING_DISTRIBUTION_FUNCTION_ON_
    if (PIC::DistributionFunctionSample::SamplingInitializedFlag==true) PIC::DistributionFunctionSample::SampleDistributionFnction(node);
#endif

    for (k=0;k<_BLOCK_CELLS_Z_;k++) {
       for (j=0;j<_BLOCK_CELLS_Y_;j++) {
          for (i=0;i<_BLOCK_CELLS_X_;i++) {
            LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
            cell=block->GetCenterNode(LocalCellNumber);
            SamplingData=cell->GetAssociatedDataBufferPointer()+PIC::Mesh::collectingCellSampleDataPointerOffset;

            ptr=cell->FirstCellParticle;

            while (ptr!=-1) {
#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_PARTICLE_NUMBER_
              TreeNodeTotalParticleNumber++;
#endif

              ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
              s=PIC::ParticleBuffer::GetI(ParticleData);
              v=PIC::ParticleBuffer::GetV(ParticleData);



//=====================  DEBUG =========================

              /*
double *x;
x=PIC::ParticleBuffer::GetX(ParticleData);
if (sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])<1500.0E3) {
  cout << __FILE__ << __LINE__ << endl;
}

*/
//===================== END DEBUG ==================




              LocalParticleWeight=block->GetLocalParticleWeight(s);
              LocalParticleWeight*=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);


              *(s+(double*)(SamplingData+PIC::Mesh::sampledParticleWeghtRelativeOffset))+=LocalParticleWeight;
              *(s+(double*)(SamplingData+PIC::Mesh::sampledParticleNumberRelativeOffset))+=1;
              *(s+(double*)(SamplingData+PIC::Mesh::sampledParticleNumberDensityRelativeOffset))+=LocalParticleWeight/cell->Measure;

              for (idim=0;idim<3;idim++) {
                *(3*s+idim+(double*)(SamplingData+PIC::Mesh::sampledParticleVelocityRelativeOffset))+=v[idim]*LocalParticleWeight;
                *(3*s+idim+(double*)(SamplingData+PIC::Mesh::sampledParticleVelocity2RelativeOffset))+=v[idim]*v[idim]*LocalParticleWeight;
              }

              ptr=PIC::ParticleBuffer::GetNext(ParticleData);
            }
          }
       }
    }


    //Sample the parallel load: the total particle number
#if _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_PARTICLE_NUMBER_
    node->ParallelLoadMeasure+=TreeNodeTotalParticleNumber;
#elif _PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_
    node->ParallelLoadMeasure+=MPI_Wtime()-TreeNodeProcessingTime;
#endif


    node=node->nextNodeThisThread;
  }

  //Increment the sample length
  CollectingSampleCounter++;

  //Output the data flow file if needed
  if (CollectingSampleCounter==RequiredSampleLength) {
    //exchnge the sampling data
    PIC::Mesh::mesh.ParallelBlockDataExchange();

    //check different sampling modes
    if (SamplingMode==_RESTART_SAMPLING_MODE_) {
      PIC::Mesh::switchSamplingBuffers();
      LastSampleLength=CollectingSampleCounter;
      CollectingSampleCounter=0;

      for (node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
        block=node->block;

        for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
          LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
          PIC::Mesh::flushCollectingSamplingBuffer(block->GetCenterNode(LocalCellNumber));
        }
      }

      //flush sampling buffers in the internal surfaces installed into the mesh
  #if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_OFF_
      //do nothing
  #elif _INTERNAL_BOUNDARY_MODE_ ==  _INTERNAL_BOUNDARY_MODE_ON_
      long int iSphericalSurface,nTotalSphericalSurfaces=PIC::BC::InternalBoundary::Sphere::InternalSpheres.usedElements();
//      PIC::BC::InternalBoundary::Sphere::cSurfaceDataSphere* Sphere;

      cInternalSphericalData *Sphere;

      for (iSphericalSurface=0;iSphericalSurface<nTotalSphericalSurfaces;iSphericalSurface++) {
//        Sphere=(PIC::BC::InternalBoundary::Sphere::cSurfaceDataSphere*)PIC::BC::InternalBoundary::Sphere::InternalSpheres.GetEntryPointer(iSphericalSurface)->GetSurfaceDataPointer();

        Sphere=PIC::BC::InternalBoundary::Sphere::InternalSpheres.GetEntryPointer(iSphericalSurface);

        PIC::BC::InternalBoundary::Sphere::flushCollectingSamplingBuffer(Sphere);
      }


  #else
      exit(__LINE__,__FILE__,"Error: unknown option");
  #endif


    }
    else if (SamplingMode==_ACCUMULATE_SAMPLING_MODE_) {
      exit(__LINE__,__FILE__,"Error: the sampling mode '_ACCUMULATE_SAMPLING_MODE_' is not implemented");
    }
    else exit(__LINE__,__FILE__,"Error: the sampling mode is not defined");

    //print output file
    char fname[_MAX_STRING_LENGTH_PIC_],ChemSymbol[_MAX_STRING_LENGTH_PIC_];

    for (s=0;s<PIC::nTotalSpecies;s++) {
      PIC::MolecularData::GetChemSymbol(ChemSymbol,s);
      sprintf(fname,"pic.%s.s=%i.out=%ld.dat",ChemSymbol,s,DataOutputFileNumber);

      if (PIC::Mesh::mesh.ThisThread==0) {
        printf("printing output file: %s.........",fname);
        fflush(stdout);
      }

      PIC::Mesh::mesh.outputMeshDataTECPLOT(fname,s);

      if (PIC::Mesh::mesh.ThisThread==0) {
        printf("done.\n");
        fflush(stdout);
      }

      //print the sampled distribution function into a file
#if _SAMPLING_DISTRIBUTION_FUNCTION_MODE_ == _SAMPLING_DISTRIBUTION_FUNCTION_ON_
      if (PIC::DistributionFunctionSample::SamplingInitializedFlag==true) {
        sprintf(fname,"pic.distribution.%s.s=%i.out=%ld",ChemSymbol,s,DataOutputFileNumber);
        PIC::DistributionFunctionSample::printDistributionFunction(fname,s);
      }
#endif

      //print into the file the sampled data of the internal surfaces installed into the mesh
#if _INTERNAL_BOUNDARY_MODE_ == _INTERNAL_BOUNDARY_MODE_OFF_
      //do nothing
#elif _INTERNAL_BOUNDARY_MODE_ ==  _INTERNAL_BOUNDARY_MODE_ON_
      long int iSphericalSurface,nTotalSphericalSurfaces=PIC::BC::InternalBoundary::Sphere::InternalSpheres.usedElements();

      for (iSphericalSurface=0;iSphericalSurface<nTotalSphericalSurfaces;iSphericalSurface++) {
        sprintf(fname,"pic.Sphere=%ld.%s.s=%i.out=%ld.dat",iSphericalSurface,ChemSymbol,s,DataOutputFileNumber);

        if (PIC::Mesh::mesh.ThisThread==0) {
          printf("printing output file: %s.........",fname);
          fflush(stdout);
        }

        PIC::BC::InternalBoundary::Sphere::InternalSpheres.GetEntryPointer(iSphericalSurface)->PrintSurfaceData(fname,s);

        if (PIC::Mesh::mesh.ThisThread==0) {
          printf("done.\n");
          fflush(stdout);
        }
      }


#else
      exit(__LINE__,__FILE__,"Error: unknown option");
#endif
    }

    //print the statistic information for the run
    CMPI_channel pipe(1000000);
    long int nTotalSimulatedParticles;

    nTotalSimulatedParticles=PIC::ParticleBuffer::NAllPart;

    if (PIC::Mesh::mesh.ThisThread!=0) {
      pipe.openSend(0);

      pipe.send(nTotalSimulatedParticles);

      pipe.closeSend();
    }
    else {
      pipe.openRecvAll();

      for (int thread=1;thread<PIC::Mesh::mesh.nTotalThreads;thread++) {
        nTotalSimulatedParticles+=pipe.recv<long int>(thread);
      }

      printf("Model run statistics:\n");
      printf("The total number of model particles: %ld\n",nTotalSimulatedParticles);

      pipe.closeRecvAll();
    }

#if _SAMPLING_DISTRIBUTION_FUNCTION_MODE_ == _SAMPLING_DISTRIBUTION_FUNCTION_ON_
    if (PIC::DistributionFunctionSample::SamplingInitializedFlag==true) PIC::DistributionFunctionSample::flushSamplingBuffers();
#endif

    //increment the output file number
    DataOutputFileNumber++;

    //redistribute the processor load
    PIC::Mesh::mesh.CreateNewParallelDistributionLists(_PIC_DYNAMIC_BALANCE_SEND_RECV_MESH_NODE_EXCHANGE_TAG_);
  }





}

//====================================================
//init the particle solver
void PIC::Init() {

  //init MPI variables
  MPI_Comm_rank(MPI_COMM_WORLD,&ThisThread);
  MPI_Comm_size(MPI_COMM_WORLD,&nTotalThreads);
}



//====================================================
//set up the particle weight and time step
void PIC::Mesh::cDataBlockAMR::SetLocalParticleWeight(double weight, int spec) {
  #if _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
  *(spec+(double *)(associatedDataPointer+LocalParticleWeightOffset))=weight;
  #elif _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_

  if (PIC::ParticleWeightTimeStep::GlobalParticleWeight==NULL) {
    PIC::ParticleWeightTimeStep::GlobalParticleWeight=new double [PIC::nTotalSpecies];
    for (int s=0;s<PIC::nTotalSpecies;s++) PIC::ParticleWeightTimeStep::GlobalParticleWeight[s]=-1.0;
  }

  if ((PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec]<0.0)||(PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec]>weight)) {
    PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec]=weight;
  }

  #else
  exit(__LINE__,__FILE__,"not implemented");
  #endif
}

double PIC::Mesh::cDataBlockAMR::GetLocalParticleWeight(int spec) {
  #if _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_LOCAL_PARTICLE_WEIGHT_
  double *res;
  res=spec+(double *)(associatedDataPointer+LocalParticleWeightOffset);
  return *res;
  #elif _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_

  return PIC::ParticleWeightTimeStep::GlobalParticleWeight[spec];
  #else
  exit(__LINE__,__FILE__,"not implemented");
  #endif
}

//set up the particle time step
void PIC::Mesh::cDataBlockAMR::SetLocalTimeStep(double dt, int spec) {
  #if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
  *(spec+(double *)(associatedDataPointer+LocalTimeStepOffset))=dt;
  #elif _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_

  if (PIC::ParticleWeightTimeStep::GlobalTimeStep==NULL) {
    PIC::ParticleWeightTimeStep::GlobalTimeStep=new double [PIC::nTotalSpecies];
    for (int s=0;s<PIC::nTotalSpecies;s++) PIC::ParticleWeightTimeStep::GlobalTimeStep[s]=-1.0;
  }

  if ((PIC::ParticleWeightTimeStep::GlobalTimeStep[spec]<0.0)||(PIC::ParticleWeightTimeStep::GlobalTimeStep[spec]>dt)) {
    PIC::ParticleWeightTimeStep::GlobalTimeStep[spec]=dt;
  }


  #else
  exit(__LINE__,__FILE__,"not implemented");
  #endif
}

double PIC::Mesh::cDataBlockAMR::GetLocalTimeStep(int spec) {

  #if _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_LOCAL_TIME_STEP_
  double *res;
  res=spec+(double *)(associatedDataPointer+LocalTimeStepOffset);
  return *res;
  #elif _SIMULATION_TIME_STEP_MODE_ == _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_
  return  PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
  #else
  exit(__LINE__,__FILE__,"not implemented");
  #endif


}





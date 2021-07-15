//the file contains parts of PIC::TimeStep

#include "pic.h"

void PIC::TimeStepInternal::PrintTimeStep() {
  using namespace PIC;
  
  if (ThisThread==0) {
    static int InteractionCouinter=0;
    time_t TimeValue=time(NULL);
    tm *ct=localtime(&TimeValue);

    printf("$PREFIX: (%i/%i %i:%i:%i), Iteration: %i  (current sample length:%ld, %ld interations to the next output)\n",ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec,InteractionCouinter,RequiredSampleLength,RequiredSampleLength-CollectingSampleCounter);
    InteractionCouinter++;
  }
}


//===============================================================================================
//recover the sampling data from the sampling data restart file, print the TECPLOT files and quit
void PIC::TimeStepInternal::RecoverSamplingDataRestart() {
  static bool RestartFileReadFlag=false;

  if (RestartFileReadFlag==false) {
    RestartFileReadFlag=true;

    //prepare the list of the files/speces which data fille be revobered and saved
    list<pair<string,list<int> > > SampledDataRecoveryTable;

    if (Restart::SamplingData::DataRecoveryManager==NULL) {
      //no user-defined fucntion to create the 'SampledDataRecoveringTable' is provided
      pair<string,list<int> > Table;

      for (int s=0;s<PIC::nTotalSpecies;s++) Table.second.push_back(s);

      Table.first=Restart::SamplingData::RestartFileName;
      SampledDataRecoveryTable.push_back(Table);
    }
    else Restart::SamplingData::DataRecoveryManager(SampledDataRecoveryTable,Restart::SamplingData::minReadFileNumber,Restart::SamplingData::maxReadFileNumber);

    //iterate through SampledDataRecoveryTable
    list<pair<string,list<int> > >::iterator RecoveryEntry;

    for (RecoveryEntry=SampledDataRecoveryTable.begin();RecoveryEntry!=SampledDataRecoveryTable.end();RecoveryEntry++) {
      Restart::SamplingData::Read(RecoveryEntry->first.c_str());
      MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

      for (list<int>::iterator s=RecoveryEntry->second.begin();s!=RecoveryEntry->second.end();s++) {
        char fname[_MAX_STRING_LENGTH_PIC_],ChemSymbol[_MAX_STRING_LENGTH_PIC_];

        PIC::MolecularData::GetChemSymbol(ChemSymbol,*s);
        sprintf(fname,"RECOVERED.%s.%s.s=%i.dat",RecoveryEntry->first.c_str(),ChemSymbol,*s);
        PIC::Mesh::mesh->outputMeshDataTECPLOT(fname,*s);

        //preplot the recovered file if needed
        if (Restart::SamplingData::PreplotRecoveredData==true) {
          char cmd[_MAX_STRING_LENGTH_PIC_];

          sprintf(cmd,"preplot %s",fname);
          system(cmd);
        }
      }
    }

    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

    if (_PIC_RECOVER_SAMPLING_DATA_RESTART_FILE__EXECUTION_MODE_ == _PIC_RECOVER_SAMPLING_DATA_RESTART_FILE__EXECUTION_MODE__STOP_) {
      if (PIC::ThisThread==0) {
        printf("The sampled data from restart file/files\n");
        for (RecoveryEntry=SampledDataRecoveryTable.begin();RecoveryEntry!=SampledDataRecoveryTable.end();RecoveryEntry++) printf("%s\n",RecoveryEntry->first.c_str());
        printf("is successfully recoved. Execution is complete. See you later :-)\n");
      }

      MPI_Finalize();
      exit(EXIT_SUCCESS);
    }
  }
}


//===============================================================================================
//recover the particle data restart file
void PIC::TimeStepInternal::ReadParticleDataRestartFile() {
  static bool RestartFileReadFlag=false;

  using namespace PIC;

  if (RestartFileReadFlag==false) {
    RestartFileReadFlag=true;

    Restart::ReadParticleData(Restart::recoverParticleDataRestartFileName);
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  }
}

//===============================================================================================
//save the particle data restart file
void PIC::TimeStepInternal::SaveParticleRestartFile() {
  static long int IterationCounter=0;
  static int saveRestartFileCounter=0;
  
  using namespace PIC;

  IterationCounter++;

  if (IterationCounter%Restart::ParticleRestartAutosaveIterationInterval==0) {
    char fname[_MAX_STRING_LENGTH_PIC_];

    if (Restart::ParticleDataRestartFileOverwriteMode==true) sprintf(fname,"%s",Restart::saveParticleDataRestartFileName);
    else sprintf(fname,"%s.restart=%i",Restart::saveParticleDataRestartFileName,saveRestartFileCounter);

    Restart::SaveParticleData(fname);
    saveRestartFileCounter++;

    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  }
}

//===============================================================================================
//simulate particle collisions
void PIC::TimeStepInternal::ParticleCollisions(double &ParticleCollisionTime) {
  SetExitErrorCode(__LINE__,_PIC__EXIT_CODE__LAST_FUNCTION__PIC_TimeStep_);
  ParticleCollisionTime=MPI_Wtime();
  
  switch (_PIC__PARTICLE_COLLISION_MODEL_) {
  case _PIC__PARTICLE_COLLISION_MODEL__NTC_:
    PIC::MolecularCollisions::ParticleCollisionModel::ntc();
    break;
    
  case _PIC__PARTICLE_COLLISION_MODEL__MF_:
    PIC::MolecularCollisions::ParticleCollisionModel::mf();
    break;
    
  case _PIC__PARTICLE_COLLISION_MODEL__USER_DEFINED_:
    exit(__LINE__,__FILE__,"Error: the option is not implemented");
    break;
    
  default:  
    exit(__LINE__,__FILE__,"Error: the option is not implemented");
  }
  
  ParticleCollisionTime=MPI_Wtime()-ParticleCollisionTime;
}



//===============================================================================================
//injection boundary conditions
void PIC::TimeStepInternal::ParticleInjectionBC(double &InjectionBoundaryTime) {
  InjectionBoundaryTime=MPI_Wtime();

  //inject particle through the domain's boundaries
  SetExitErrorCode(__LINE__,_PIC__EXIT_CODE__LAST_FUNCTION__PIC_TimeStep_);

  PIC::BC::InjectionBoundaryConditions();

  //inject particles into the volume of the domain
  #if _PIC_VOLUME_PARTICLE_INJECTION_MODE_ == _PIC_VOLUME_PARTICLE_INJECTION_MODE__ON_
  if (PIC::VolumeParticleInjection::nRegistratedInjectionProcesses!=0) PIC::BC::nTotalInjectedParticles+=PIC::VolumeParticleInjection::InjectParticle();
  #endif

  //call a user-defined injection function
  if (PIC::BC::UserDefinedParticleInjectionFunction!=NULL) {
    SetExitErrorCode(__LINE__,_PIC__EXIT_CODE__LAST_FUNCTION__PIC_TimeStep_);
    PIC::BC::nTotalInjectedParticles+=PIC::BC::UserDefinedParticleInjectionFunction();
  }

  //the extra injection process by the exosphere model (src/models/exosphere)
  if (BC::ExosphereModelExtraInjectionFunction!=NULL) {
    SetExitErrorCode(__LINE__,_PIC__EXIT_CODE__LAST_FUNCTION__PIC_TimeStep_);
    PIC::BC::nTotalInjectedParticles+=BC::ExosphereModelExtraInjectionFunction();
  }

  InjectionBoundaryTime=MPI_Wtime()-InjectionBoundaryTime;
}


//===============================================================================================
//initialization at the first call of PIC::TimeStep()
void PIC::TimeStepInternal::Init() {
  if (_SIMULATION_TIME_STEP_MODE_ != _SINGLE_GLOBAL_TIME_STEP_) exit(__LINE__,__FILE__,"Error: the option is only implemented for single global time step");

  if (PIC::Sampling::SampleTimeInterval<0) exit(__LINE__,__FILE__,"Error: the sample time interval is negative");

  if (PIC::ParticleWeightTimeStep::GetGlobalTimeStep(0)<0) {
    exit(__LINE__,__FILE__,"Error: the global time step is negative");
  }
  else {
    PIC::RequiredSampleLength=(long int)(PIC::Sampling::SampleTimeInterval/PIC::ParticleWeightTimeStep::GetGlobalTimeStep(0)+0.5);

    if (PIC::RequiredSampleLength==0) {
      char msg[600];
      sprintf(msg,"Error: the required sample length is 0, PIC::Sampling::SampleTimeInterval:%e, PIC::ParticleWeightTimeStep::GetGlobalTimeStep(0):%e", PIC::Sampling::SampleTimeInterval,PIC::ParticleWeightTimeStep::GetGlobalTimeStep(0));
      exit(__LINE__,__FILE__,msg);
    }

    if (PIC::ThisThread==0) printf("PIC::RequiredSampleLength is set to %ld \n", PIC::RequiredSampleLength);
  } 
}

//===============================================================================================
//sampling of the particle properties and creating an output file
void PIC::TimeStepInternal::Sampling(double &SamplingTime) {
  SetExitErrorCode(__LINE__,_PIC__EXIT_CODE__LAST_FUNCTION__PIC_TimeStep_);
  SamplingTime=MPI_Wtime();

  if ((_PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_) && 
      (_PIC_DEBUGGER_MODE__SAMPLING_BUFFER_VALUE_RANGE_CHECK_ == _PIC_DEBUGGER_MODE__VARIABLE_VALUE_RANGE_CHECK_ON_)) {
    Sampling::CatchOutLimitSampledValue();
  }

  timing_start("PT::Sampling");
  PIC::Sampling::Sampling();
  timing_stop("PT::Sampling");


  if ((_PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_) && 
      (_PIC_DEBUGGER_MODE__SAMPLING_BUFFER_VALUE_RANGE_CHECK_ == _PIC_DEBUGGER_MODE__VARIABLE_VALUE_RANGE_CHECK_ON_)) {
    Sampling::CatchOutLimitSampledValue();
  }

  SamplingTime=MPI_Wtime()-SamplingTime;
}





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
void PIC::TimeStepInternal::ParticleCollisions() {
  
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
}













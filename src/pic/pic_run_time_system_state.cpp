//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//$Id$
//calculate integrated parameters that describes the "state" of the entire simulation

#include "pic.h"

//get the check sum of all particles
void PIC::RunTimeSystemState::GetIndividualParticleFieldCheckSum_CallCounter(void *ParticleData,const char *msg) {
  CRC32 CheckSum;
  static int CallCounter=0;

  CheckSum.add<char>((char*)ParticleData,PIC::ParticleBuffer::ParticleDataLength);

  printf("$PREFIX: Particle Field Check Sum ");
  if (msg!=NULL) printf("(message=\"%s\")",msg);
  printf(":  0x%lx [Counter=%ld]\n",CheckSum.checksum(),CallCounter);
  CallCounter++;
}

void PIC::RunTimeSystemState::GetIndividualParticleFieldCheckSum_CallCounter(long int ptr,const char *msg) {
  GetIndividualParticleFieldCheckSum_CallCounter((void*)PIC::ParticleBuffer::GetParticleDataPointer(ptr),msg);
}

void PIC::RunTimeSystemState::GetIndividualParticleFieldCheckSum_CallCounter(void *ParticleData,long int nline,const char *fname,const char *msg) {
  char buffer[_MAX_STRING_LENGTH_PIC_];

  sprintf(buffer,"file=%s, line=%ld",fname,nline);
  GetIndividualParticleFieldCheckSum_CallCounter(ParticleData,buffer);
}

void PIC::RunTimeSystemState::GetIndividualParticleFieldCheckSum_CallCounter(long int ptr,long int nline,const char *fname,const char *msg) {
  GetIndividualParticleFieldCheckSum_CallCounter((void*)PIC::ParticleBuffer::GetParticleDataPointer(ptr),nline,fname,msg);
}

void PIC::RunTimeSystemState::GetParticleFieldCheckSum(const char *msg) {
  CRC32 CheckSum;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];
  PIC::ParticleBuffer::byte *ParticleData;
  int i,j,k,ParticleCounter=0;
  long int ptr;

  //calcualte the check sum of the particles with in the domain
  while (node!=NULL) {
    for (k=0;k<_BLOCK_CELLS_Z_;k++) {
      for (j=0;j<_BLOCK_CELLS_Y_;j++)  {
        for (i=0;i<_BLOCK_CELLS_X_;i++)  {
          for (int npass=0;npass<2;npass++) {
            ptr=(npass==0) ? node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)] : node->block->tempParticleMovingListTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

            while (ptr!=-1) {
              ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
              CheckSum.add<char>((char*)ParticleData,PIC::ParticleBuffer::ParticleDataLength);
              ParticleCounter++;

              ptr=PIC::ParticleBuffer::GetNext(ParticleData);
            }

          }
        }
      }
    }

    node=node->nextNodeThisThread;
  }

  //calcualte the check sum of particles in the 'ghost'cells
  for (int To=0;To<PIC::Mesh::mesh.nTotalThreads;To++) if ((PIC::ThisThread!=To)&&(PIC::Mesh::mesh.ParallelSendRecvMap[PIC::ThisThread][To]==true)) {
    for (node=PIC::Mesh::mesh.DomainBoundaryLayerNodesList[To];node!=NULL;node=node->nextNodeThisThread) {
      for (k=0;k<_BLOCK_CELLS_Z_;k++) {
        for (j=0;j<_BLOCK_CELLS_Y_;j++)  {
          for (i=0;i<_BLOCK_CELLS_X_;i++)  {
            for (int npass=0;npass<2;npass++) {
              ptr=(npass==0) ? node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)] : node->block->tempParticleMovingListTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

              while (ptr!=-1) {
                ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
                CheckSum.add<char>((char*)ParticleData,PIC::ParticleBuffer::ParticleDataLength);
                ParticleCounter++;

                ptr=PIC::ParticleBuffer::GetNext(ParticleData);
              }

            }
          }
        }
      }

    }
  }

  unsigned long CheckSumBuffer[PIC::nTotalThreads],t;
  int ParticleCounterBuffer[PIC::nTotalThreads];

  t=CheckSum.checksum();
  MPI_Gather(&t,1,MPI_LONG,CheckSumBuffer,1,MPI_LONG,0,MPI_GLOBAL_COMMUNICATOR);
  MPI_Gather(&ParticleCounter,1,MPI_INT,ParticleCounterBuffer,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);

  if (PIC::ThisThread==0) {
     printf("$PREFIX: Particle Field Check Sum:");
     if (msg!=NULL) printf("(message=\"%s\") \n",msg);

     for (int thread=0;thread<PIC::nTotalThreads;thread++) printf("  0x%lx",CheckSumBuffer[thread]);
     printf("\n");

     for (int thread=0;thread<PIC::nTotalThreads;thread++) printf("  %ld",ParticleCounterBuffer[thread]);
     printf("\n");
  }
}

void PIC::RunTimeSystemState::GetParticleFieldCheckSum(long int nline,const char *fname) {
  char msg[_MAX_STRING_LENGTH_PIC_];

  sprintf(msg,"file=%s, line=%ld",fname,nline);
  GetParticleFieldCheckSum(msg);
}

void PIC::RunTimeSystemState::GetParticleFieldCheckSum_CallCounter(long int nline,const char *fname) {
  char msg[_MAX_STRING_LENGTH_PIC_];
  static long int CallCounter=0;

  sprintf(msg,"file=%s, line=%ld [Counter=%i]",fname,nline,CallCounter);
  GetParticleFieldCheckSum(msg);
  CallCounter++;
}

void PIC::RunTimeSystemState::GetParticleFieldCheckSum_CallCounter(const char *msg) {
  char buffer[_MAX_STRING_LENGTH_PIC_];
  static long int CallCounter=0;

  if (msg!=NULL) {
    sprintf(buffer,"%s [Counter=%i]",msg,CallCounter);
  }
  else {
    sprintf(buffer,"[Counter=%i]",CallCounter);
  }

  GetParticleFieldCheckSum(msg);
  CallCounter++;
}


//compare the domain decomposition: calcualte the chech sum of all block's TempID belong to the currect processor
void PIC::RunTimeSystemState::GetDomainDecompositionCheckSum(const char *msg) {
  CRC32 CheckSum;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];

  while (node!=NULL) {
    CheckSum.add(node->Temp_ID);
    node=node->nextNodeThisThread;
  }

  unsigned long CheckSumBuffer[PIC::nTotalThreads],t;

  t=CheckSum.checksum();
  MPI_Gather(&t,1,MPI_LONG,CheckSumBuffer,1,MPI_LONG,0,MPI_GLOBAL_COMMUNICATOR);

  if (PIC::ThisThread==0) {
    printf("$PREFIX: Domain Decomposition Check Sum:");
    if (msg!=NULL) printf("(message=\"%s\") ",msg);

    for (int thread=0;thread<PIC::nTotalThreads;thread++) printf("  0x%lx",CheckSumBuffer[thread]);
    printf("\n");
  }
}

void PIC::RunTimeSystemState::GetDomainDecompositionCheckSum(long int nline,const char *fname) {
  char msg[_MAX_STRING_LENGTH_PIC_];

  sprintf(msg,"file=%s, line=%ld",fname,nline);
  GetDomainDecompositionCheckSum(msg);
}

void PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(FILE* fout,const char *msg) {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];
  PIC::ParticleBuffer::byte *ParticleData;
  int i,j,k,idim;
  long int ptr;

  int s;
  double *v;
  double StatWeight,TotalStatWeight[PIC::nTotalSpecies],MeanSpeed[PIC::nTotalSpecies],MeanVelocity[3*PIC::nTotalSpecies];

  for (s=0;s<PIC::nTotalSpecies;s++) {
    MeanSpeed[s]=0.0;
    TotalStatWeight[s]=0.0;

    for (idim=0;idim<3;idim++) MeanVelocity[idim+3*s]=0.0;
  }

  //sample the mean values
  for (int thread=0;thread<PIC::nTotalThreads;thread++) {
    for (node=((thread==PIC::ThisThread) ? PIC::Mesh::mesh.ParallelNodesDistributionList[thread] : PIC::Mesh::mesh.DomainBoundaryLayerNodesList[thread]);node!=NULL;node=node->nextNodeThisThread) {

      for (k=0;k<_BLOCK_CELLS_Z_;k++) {
        for (j=0;j<_BLOCK_CELLS_Y_;j++)  {
          for (i=0;i<_BLOCK_CELLS_X_;i++)  {
            ptr=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

            while (ptr!=-1) {
              ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
              s=PIC::ParticleBuffer::GetI(ParticleData);

              StatWeight=node->block->GetLocalParticleWeight(s)*PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData);
              v=PIC::ParticleBuffer::GetV(ParticleData);

              TotalStatWeight[s]+=StatWeight;
              MeanSpeed[s]+=StatWeight*sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
              for (idim=0;idim<3;idim++) MeanVelocity[idim+3*s]+=StatWeight*v[idim];

              ptr=PIC::ParticleBuffer::GetNext(ParticleData);
            }

          }
        }
      }

    }
  }

  double TempBuffer[3*PIC::nTotalSpecies];

  MPI_Reduce(TotalStatWeight,TempBuffer,PIC::nTotalSpecies,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
  memcpy(TotalStatWeight,TempBuffer,PIC::nTotalSpecies*sizeof(double));

  MPI_Reduce(MeanSpeed,TempBuffer,PIC::nTotalSpecies,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
  memcpy(MeanSpeed,TempBuffer,PIC::nTotalSpecies*sizeof(double));

  MPI_Reduce(MeanVelocity,TempBuffer,3*PIC::nTotalSpecies,MPI_DOUBLE,MPI_SUM,0,MPI_GLOBAL_COMMUNICATOR);
  memcpy(MeanVelocity,TempBuffer,3*PIC::nTotalSpecies*sizeof(double));

  if (PIC::ThisThread==0) {
    fprintf(fout,"$PREFIX: Averaged particles microscopic aprameters:");
    if (msg!=NULL) fprintf(fout,"(message=\"%s\") ",msg);
    fprintf(fout,"\n$PREFIX: spec\t<Speed>\t\t<v[0]>\t\t<v[1]>\t\t<v[2]>\n");

    for (s=0;s<PIC::nTotalSpecies;s++) {
      if (TotalStatWeight[s]>0.0) {
        MeanSpeed[s]/=TotalStatWeight[s];

        for (idim=0;idim<3;idim++) MeanVelocity[3*s+idim]/=TotalStatWeight[s];
      }

      fprintf(fout,"$PREFIX: %i\t%e\t%e\t%e\t%e\t\n",s,MeanSpeed[s],MeanVelocity[0+3*s],MeanVelocity[1+3*s],MeanVelocity[2+3*s]);
    }

    fprintf(fout,"\n");
  }

  fflush(fout);
}

void PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(FILE* fout,long int nline,const char *fname) {
  char msg[_MAX_STRING_LENGTH_PIC_];

  sprintf(msg,"file=%s, line=%ld",fname,nline);
  GetMeanParticleMicroscopicParameters(fout,msg);
}

void PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(const char *fname) {
  FILE *fout=NULL;

  if (PIC::ThisThread==0) fout=fopen(fname,"w");
  GetMeanParticleMicroscopicParameters(fout,NULL);
  if (PIC::ThisThread==0) fclose(fout);
}


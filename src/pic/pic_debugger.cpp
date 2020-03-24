/*
 * pic_debugger.cpp
 *
 *  Created on: Feb 19, 2014
 *      Author: vtenishe
 */
//$Id$
//contains functions used for debugging AMPS

#include "pic.h"


//Save particle data into a debugger data stream
void PIC::Debugger::SaveParticleDataIntoDebuggerDataStream(void* data,int length,int nline,const char* fname) {
  char msg[_MAX_STRING_LENGTH_PIC_];

  sprintf(msg,"%s, line %i",fname,nline);
  PIC::Debugger::SaveParticleDataIntoDebuggerDataStream(data,length,msg);
}


void PIC::Debugger::SaveParticleDataIntoDebuggerDataStream(void* data,int length,const char* msg) {

  struct cStreamBuffer {
    int CollCounter;
    unsigned long CheckSum;
    char CallPoint[200];
  };

  const int CheckSumBufferLength=500;
  static cStreamBuffer StreamBuffer[CheckSumBufferLength];

  static int BufferPointer=0;
  static int CallCounter=0;
  static CRC32 CheckSum;

  //create a new copy of the dbuffer stream file at the first call of this function
  if (CallCounter==0) {
    FILE *fout;
    char fn[200];

    sprintf(fn,"DebuggerStream.thread=%i.dbg",PIC::ThisThread);
    fout=fopen(fn,"w");
    fclose(fout);
  }

  //increment the call counter
  CallCounter++;

  //update the check sum
  CheckSum.add((char*)data,length);

  //save the checksum in the buffer
  StreamBuffer[BufferPointer].CollCounter=CallCounter;
  StreamBuffer[BufferPointer].CheckSum=CheckSum.checksum();
  sprintf(StreamBuffer[BufferPointer].CallPoint,"%s",msg);
  BufferPointer++;

  if (BufferPointer>=CheckSumBufferLength) {
    //save the accumulated checksums into a file
    FILE *fout;
    char fn[200];

    sprintf(fn,"DebuggerStream.thread=%i.dbg",PIC::ThisThread);
    fout=fopen(fn,"a");

    for (int i=0;i<BufferPointer;i++) fprintf(fout,"%i: 0x%lx\t%s\n",StreamBuffer[i].CollCounter,StreamBuffer[i].CheckSum,StreamBuffer[i].CallPoint);

    BufferPointer=0;
    fclose(fout);
  }
}


//InfiniteLoop==false ==> no problem found; InfiniteLoop==true => the actual number of particles does not consider the that in teh particle buffer
bool PIC::Debugger::InfiniteLoop(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  int nDownNode,i,j,k;
  bool res=false;
  long int ptr;
  static long int nAllCountedParticles=0;

  if (startNode==NULL) startNode=PIC::Mesh::mesh.rootTree;
  if (startNode==PIC::Mesh::mesh.rootTree) nAllCountedParticles=0;


  for (nDownNode=0;nDownNode<(1<<DIM);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) InfiniteLoop(startNode->downNode[nDownNode]);

  if (startNode->block!=NULL) {
    //the block is allocated; check the particle lists associated with the block
    int nTotalThreads_OpenMP=1,thread_OpenMP;

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
    nTotalThreads_OpenMP=omp_get_thread_num();
#endif


    for (thread_OpenMP=0;thread_OpenMP<nTotalThreads_OpenMP;thread_OpenMP++) for (i=0;i<_BLOCK_CELLS_X_;i++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (k=0;k<_BLOCK_CELLS_Z_;k++) {

#if _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
      ptr=startNode->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
#elif _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
      ptr=startNode->block->GetTempParticleMovingListMultiThreadTable(thread_OpenMP,i,j,k)->first;
#else
#error The option is unknown
#endif

      while (ptr!=-1) {
        ptr=PIC::ParticleBuffer::GetNext(ptr);

        if (++nAllCountedParticles>PIC::ParticleBuffer::NAllPart) {
          FindDoubleReferencedParticle();

          exit(__LINE__,__FILE__,"The counted particle number exeeds the number of particles stored in the particle buffer");
        }
      }
    }
  }

  if (startNode==PIC::Mesh::mesh.rootTree) {
    //collect the information from all processors;
    int Flag,FlagTable[PIC::nTotalThreads];

    Flag=(nAllCountedParticles==PIC::ParticleBuffer::NAllPart) ? true : false;
    MPI_Gather(FlagTable,1,MPI_INT,&Flag,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);

    if (PIC::ThisThread==0) {
      for (int thread=0;thread<PIC::nTotalThreads;thread++) if (FlagTable[thread]==false) {
        printf("$PREFIX: the actual number of particles does not coniside with that in the particle buffer (Thread=%i,line=%i,file=%s)\n",thread,__LINE__,__FILE__);
        res=true;
      }
    }

    MPI_Bcast(&Flag,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);
    nAllCountedParticles=0;
  }

  return res;
}



void PIC::Debugger::FindDoubleReferencedParticle(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  int nDownNode,i,j,k;
  long int ptr,ptrNext;
  static long int nAllCountedParticles=0;

  static char* ParticleAllocationTable=NULL;
  static long int* ptrPrevTable=NULL;

  if (startNode==NULL) startNode=PIC::Mesh::mesh.rootTree;

  //create table of allocated/not-allocated flag table
  if (startNode==PIC::Mesh::mesh.rootTree) {
    long int cnt=0;

    nAllCountedParticles=0;
    ParticleAllocationTable=new char [PIC::ParticleBuffer::MaxNPart];
    ptrPrevTable=new long int [PIC::ParticleBuffer::MaxNPart];

    //the definition of the bits of ParticleAllocationTable[]
    //ParticleAllocationTable[] & 1: == 0 --> particle is not allocated; == 1 --> the particle is allocated
    //ParticleAllocationTable[] & 2: == 0 --> the particle is not lonked yet; == 2 --> the particle is already linked

    for (ptr=0;ptr<PIC::ParticleBuffer::MaxNPart;ptr++) ptrPrevTable[ptr]=-1,ParticleAllocationTable[ptr]=1; //particle is allocated

    if (PIC::ParticleBuffer::FirstPBufferParticle!=-1) ParticleAllocationTable[PIC::ParticleBuffer::FirstPBufferParticle]=2; //particle is links and not allocated

    for (cnt=0,ptr=PIC::ParticleBuffer::FirstPBufferParticle;ptr!=-1;ptr=ptrNext) {
      ptrNext=PIC::ParticleBuffer::GetNext(ptr);

      if (ptrNext!=-1) {
        if ((ParticleAllocationTable[ptrNext]&2)==2) {
          printf("Error: have found double-referenced particle in the list of un-allocated particles\n%ld --> %ld\n%ld --> %ld\n",ptr,ptrNext,ptrPrevTable[ptrNext],ptrNext);

          exit(__LINE__,__FILE__,"Error: an un-allocated particle is double-referenced");
        }
        else ParticleAllocationTable[ptrNext]=2,ptrPrevTable[ptrNext]=ptr;
      }

      if (++cnt>PIC::ParticleBuffer::MaxNPart) {
        exit(__LINE__,__FILE__,"The counted particle number exeeds the number of particles stored in the particle buffer");
      }

      ptr=ptrNext;
    }
  }

  //check the particle lists
  for (nDownNode=0;nDownNode<(1<<DIM);nDownNode++) if (startNode->downNode[nDownNode]!=NULL) FindDoubleReferencedParticle(startNode->downNode[nDownNode]);

  if (startNode->block!=NULL) {
    //the block is allocated; check the particle lists associated with the block

    int nTotalThreads_OpenMP=1,thread_OpenMP;

  #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
    nTotalThreads_OpenMP=omp_get_thread_num();
  #endif

    for (i=0;i<_BLOCK_CELLS_X_;i++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (k=0;k<_BLOCK_CELLS_Z_;k++) for (int npass=0;npass<2;npass++)  for (thread_OpenMP=0;thread_OpenMP<((npass==0) ? 1 : nTotalThreads_OpenMP);thread_OpenMP++) {
      if (npass==0) ptr=startNode->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
      else {
#if _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
        ptr=startNode->block->tempParticleMovingListTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
#elif _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
        ptr=startNode->block->GetTempParticleMovingListMultiThreadTable(thread_OpenMP,i,j,k)->first;
#endif
      }

//      ptr=(npass==0) ? startNode->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)] : startNode->block->tempParticleMovingListTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

      if (ptr!=-1) {
        if ((ParticleAllocationTable[ptr]&1)==0) {
          exit(__LINE__,__FILE__,"Error: un-allocated particles is found in the particle list");
        }

        if ((ParticleAllocationTable[ptr]&2)==2) {
          printf("Error: the first particle in the list is referenced: %ld --> %ld\n%ld --> %ld\n",ptr,ptrNext,ptrPrevTable[ptr],ptr);
          exit(__LINE__,__FILE__,"Error: the first particle in the list is referenced");
        }
        else ParticleAllocationTable[ptr]|=2;

        if (PIC::ParticleBuffer::GetPrev(ptr)>=0) {
          exit(__LINE__,__FILE__,"Error: the first particle in the list has non-negarive value of GetPrev");
        }
      }


      while (ptr!=-1) {
        if ((ParticleAllocationTable[ptr]&1)==0) {
          exit(__LINE__,__FILE__,"Error: un-allocated particles is found in the particle list");
        }

        ptrNext=PIC::ParticleBuffer::GetNext(ptr);

        if (ptrNext!=-1) {
          if ((ParticleAllocationTable[ptrNext]&1)==0) {
            exit(__LINE__,__FILE__,"Error: reference to an un-allocated particles is found in the particle list");
          }

          if ((ParticleAllocationTable[ptrNext]&2)==2) {
            printf("Error: have found double-referenced particle in the list of un-allocated particles\n%ld --> %ld\n%ld --> %ld\n",ptr,ptrNext,ptrPrevTable[ptrNext],ptrNext);
            exit(__LINE__,__FILE__,"Error: have found double-referenced particle in the list");
          }

          if (PIC::ParticleBuffer::GetPrev(ptrNext)!=ptr) {
            exit(__LINE__,__FILE__,"Error: PIC::ParticleBuffer::GetPrev(ptrNext) != ptr");
          }

          ParticleAllocationTable[ptrNext]|=2;
          ptrPrevTable[ptrNext]=ptr;
        }

        ptr=ptrNext;

        if (++nAllCountedParticles>PIC::ParticleBuffer::NAllPart) {
          exit(__LINE__,__FILE__,"The counted particle number exeeds the number of particles stored in the particle buffer");
        }
      }
    }
  }

  if (startNode==PIC::Mesh::mesh.rootTree) {
    for (ptr=0;ptr<PIC::ParticleBuffer::MaxNPart;ptr++) {
      if ((ParticleAllocationTable[ptr]&2)==0) {
        exit(__LINE__,__FILE__,"Error: found particles that is not referenced at all");
      }
    }

    delete [] ParticleAllocationTable;
    delete [] ptrPrevTable;
  }
}


//catch the out of limit value in the sample buffer (check only the base quantaty)
void PIC::Sampling::CatchOutLimitSampledValue() {
#if _PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_
#if _PIC_DEBUGGER_MODE__SAMPLING_BUFFER_VALUE_RANGE_CHECK_ == _PIC_DEBUGGER_MODE__VARIABLE_VALUE_RANGE_CHECK_ON_

  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  int s,i,j,k;
  PIC::Mesh::cDataCenterNode *cell;
  PIC::Mesh::cDataBlockAMR *block;
  long int LocalCellNumber;
  char *SamplingBuffer;



  for (node=PIC::Mesh::mesh.ParallelNodesDistributionList[PIC::Mesh::mesh.ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    block=node->block;

    for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++)  for (i=0;i<_BLOCK_CELLS_X_;i++) {
      LocalCellNumber=_getCenterNodeLocalNumber(i,j,k);
      cell=block->GetCenterNode(LocalCellNumber);

      SamplingBuffer=cell->GetAssociatedDataBufferPointer()+PIC::Mesh::collectingCellSampleDataPointerOffset;


      for (s=0;s<PIC::nTotalSpecies;s++) {
        if (PIC::Mesh::DatumParticleWeight.is_active()==true) {
          PIC::Debugger::CatchOutLimitValue((s+(double*)(SamplingBuffer+PIC::Mesh::DatumParticleWeight.offset)),1,__LINE__,__FILE__);
        }

        if (PIC::Mesh::DatumParticleNumber.is_active()==true) {
          PIC::Debugger::CatchOutLimitValue((s+(double*)(SamplingBuffer+PIC::Mesh::DatumParticleNumber.offset)),1,__LINE__,__FILE__);
        }

        if (PIC::Mesh::DatumNumberDensity.is_active()==true) {
          PIC::Debugger::CatchOutLimitValue((s+(double*)(SamplingBuffer+PIC::Mesh::DatumNumberDensity.offset)),1,__LINE__,__FILE__);
        }

  
        if (PIC::Mesh::DatumParticleVelocity.is_active()==true) {
          PIC::Debugger::CatchOutLimitValue((3*s+(double*)(SamplingBuffer+PIC::Mesh::DatumParticleVelocity.offset)),DIM,__LINE__,__FILE__);
        }


        if (PIC::Mesh::DatumParticleVelocity2.is_active()==true) {
          PIC::Debugger::CatchOutLimitValue((3*s+(double*)(SamplingBuffer+PIC::Mesh::DatumParticleVelocity2.offset)),DIM,__LINE__,__FILE__);
        }


        if (PIC::Mesh::DatumParticleSpeed.is_active()==true) {
          PIC::Debugger::CatchOutLimitValue((s+(double*)(SamplingBuffer+PIC::Mesh::DatumParticleSpeed.offset)),1,__LINE__,__FILE__);
        }

      }

    }
  }

#endif
#endif
}


//==========================================================================================
//get checksum of the corner and center node associated data
unsigned long int PIC::Debugger::SaveCornerNodeAssociatedDataSignature(long int nline,const char* fnameSource,const char* fnameOutput,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  return PIC::Debugger::SaveCornerNodeAssociatedDataSignature(0,PIC::Mesh::cDataCornerNode::totalAssociatedDataLength,nline,fnameSource,fnameOutput,startNode);
}


unsigned long int PIC::Debugger::SaveCornerNodeAssociatedDataSignature(int SampleVectorOffset,int SampleVectorLength,long int nline,const char* fnameSource,const char* fnameOutput,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  int i,j,k;
  PIC::Mesh::cDataCornerNode *CornerNode;
  PIC::Mesh::cDataBlockAMR *block;

  static int EntryCounter;  //used to ensure the same allocation of the blocks and nodes
  static CRC32 CheckSum,SingleVectorCheckSum;
  static CMPI_channel pipe(1000000);
  static int initflag=false;
  static FILE *fout=NULL;

  //coundate of the fucntion calls
  static int nCallCounter=0;

  if (startNode==NULL) {
    startNode=PIC::Mesh::mesh.rootTree;
    EntryCounter=0;
    CheckSum.clear();

    nCallCounter++;

    if ((fnameOutput!=NULL)&&(PIC::ThisThread==0)) fout=fopen(fnameOutput,"w");

    if (initflag==false) {
      initflag=true;

      if (PIC::ThisThread==0) {
        pipe.openRecvAll();
      }
      else {
        pipe.openSend(0);
      }
    }
  }

  //add the associated node data
  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    EntryCounter++;
    CheckSum.add(EntryCounter);

    if ((startNode->Thread==PIC::ThisThread)||(PIC::ThisThread==0)) {
      block=startNode->block;

      for (k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_+1;k++) {
        for (j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_+1;j++) {
          for (i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_+1;i++) {
            EntryCounter++;
            CheckSum.add(EntryCounter);
            SingleVectorCheckSum.clear();

            if (startNode->Thread==0) {
              //the associated data is located the the root
              if (block!=NULL) if ((CornerNode=block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k)))!=NULL) {
                CheckSum.add(CornerNode->GetAssociatedDataBufferPointer()+SampleVectorOffset,SampleVectorLength);

                if (fnameOutput!=NULL) {
                  SingleVectorCheckSum.add(CornerNode->GetAssociatedDataBufferPointer()+SampleVectorOffset,SampleVectorLength);
                }
              }

              if (fnameOutput!=NULL) {
                fprintf(fout,"node: id=%ld, i=%i, j=%i, k=%i, CheckSum=0x%lx\n",startNode->Temp_ID,i,j,k,SingleVectorCheckSum.checksum());
              }
            }
            else {
              if (PIC::ThisThread==0) {
                MPI_Status status;

                MPI_Send(&CheckSum.crc_accum,1,MPI_UNSIGNED_LONG,startNode->Thread,0,MPI_GLOBAL_COMMUNICATOR);
                MPI_Recv(&CheckSum.crc_accum,1,MPI_UNSIGNED_LONG,startNode->Thread,0,MPI_GLOBAL_COMMUNICATOR,&status);

                if (fnameOutput!=NULL) {
                  MPI_Recv(&SingleVectorCheckSum.crc_accum,1,MPI_UNSIGNED_LONG,startNode->Thread,0,MPI_GLOBAL_COMMUNICATOR,&status);

                  fprintf(fout,"node: id=%ld, i=%i, j=%i, k=%i, CheckSum=0x%lx\n",startNode->Temp_ID,i,j,k,SingleVectorCheckSum.checksum());
                }
              }
              else {
                MPI_Status status;
                MPI_Recv(&CheckSum.crc_accum,1,MPI_UNSIGNED_LONG,0,0,MPI_GLOBAL_COMMUNICATOR,&status);

                if (block!=NULL) if ((CornerNode=block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k)))!=NULL) {
                  CheckSum.add(CornerNode->GetAssociatedDataBufferPointer()+SampleVectorOffset,SampleVectorLength);

                  if (fnameOutput!=NULL) {
                    SingleVectorCheckSum.add(CornerNode->GetAssociatedDataBufferPointer()+SampleVectorOffset,SampleVectorLength);
                  }
                }

                MPI_Send(&CheckSum.crc_accum,1,MPI_UNSIGNED_LONG,0,0,MPI_GLOBAL_COMMUNICATOR);

                if (fnameOutput!=NULL) {
                  MPI_Send(&SingleVectorCheckSum.crc_accum,1,MPI_UNSIGNED_LONG,0,0,MPI_GLOBAL_COMMUNICATOR);
                }

              }
            }

          }

        }
      }
    }
  }
  else {
    //add daugher blocks
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *downNode;

    for (i=0;i<(1<<DIM);i++) if ((downNode=startNode->downNode[i])!=NULL) SaveCornerNodeAssociatedDataSignature(SampleVectorOffset,SampleVectorLength,nline,fnameSource,fnameOutput,downNode);
  }

  //output the checksum
  if (startNode==PIC::Mesh::mesh.rootTree) {
    pipe.flush();

    if (PIC::ThisThread==0) {
      char msg[500];

      sprintf(msg," line=%ld, file=%s (Call Counter=%i)",nline,fnameSource,nCallCounter);
      CheckSum.PrintChecksumSingleThread(msg);

      if (fnameOutput!=NULL) {
        fclose(fout);
        fout=NULL;
      }
    }
  }

  return CheckSum.checksum();
}

unsigned long int PIC::Debugger::GetCornerNodeAssociatedDataSignature(long int nline,const char* fname,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  return SaveCornerNodeAssociatedDataSignature(nline,fname,NULL,startNode);
}

unsigned long int PIC::Debugger::GetCenterNodeAssociatedDataSignature(long int nline,const char* fname,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  return SaveCenterNodeAssociatedDataSignature(nline,fname,NULL,startNode);
}

unsigned long int PIC::Debugger::SaveCenterNodeAssociatedDataSignature(long int nline,const char* fnameSource,const char* fnameOutput,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  return SaveCenterNodeAssociatedDataSignature(0,PIC::Mesh::cDataCenterNode::totalAssociatedDataLength,nline,fnameSource,fnameOutput,startNode);
}


unsigned long int PIC::Debugger::SaveCenterNodeAssociatedDataSignature(int SampleVectorOffset, int SampleVectorLength, long int nline,const char* fnameSource,const char* fnameOutput,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *startNode) {
  int i,j,k;
  PIC::Mesh::cDataCenterNode *CenterNode;
  PIC::Mesh::cDataBlockAMR *block;

  static int EntryCounter;  //used to ensure the same allocation of the blocks and nodes
  static CRC32 CheckSum,SingleVectorCheckSum;
  static CMPI_channel pipe(1000000);
  static int initflag=false;
  static FILE *fout=NULL;

  //coundate of the fucntion calls
  static int nCallCounter=0;

  if (startNode==NULL) {
    startNode=PIC::Mesh::mesh.rootTree;
    EntryCounter=0;
    CheckSum.clear();

    if ((fnameOutput!=NULL)&&(PIC::ThisThread==0)) fout=fopen(fnameOutput,"w");

    nCallCounter++;

    if (initflag==false) {
      initflag=true;

      if (PIC::ThisThread==0) {
        pipe.openRecvAll();
      }
      else {
        pipe.openSend(0);
      }
    }
  }

  //add the associated node data
  if (startNode->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    EntryCounter++;
    CheckSum.add(EntryCounter);

    if ((startNode->Thread==PIC::ThisThread)||(PIC::ThisThread==0)) {
      block=startNode->block;

      for (k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;k++) {
        for (j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;j++) {
          for (i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;i++) {
            EntryCounter++;
            CheckSum.add(EntryCounter);
            SingleVectorCheckSum.clear();

            if (startNode->Thread==0) {
              //the associated data is located the the root
              if (block!=NULL) if ((CenterNode=block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k)))!=NULL) {
                CheckSum.add(CenterNode->GetAssociatedDataBufferPointer()+SampleVectorOffset,SampleVectorLength);

                if (fnameOutput!=NULL) {
                  SingleVectorCheckSum.add(CenterNode->GetAssociatedDataBufferPointer()+SampleVectorOffset,SampleVectorLength);
                }
              }

              if (fnameOutput!=NULL) {
                fprintf(fout,"node: id=%ld, i=%i, j=%i, k=%i, CheckSum=0x%lx\n",startNode->Temp_ID,i,j,k,SingleVectorCheckSum.checksum());
              }

            }
            else {
              if (PIC::ThisThread==0) {
                MPI_Status status;

                MPI_Send(&CheckSum.crc_accum,1,MPI_UNSIGNED_LONG,startNode->Thread,0,MPI_GLOBAL_COMMUNICATOR);
                MPI_Recv(&CheckSum.crc_accum,1,MPI_UNSIGNED_LONG,startNode->Thread,0,MPI_GLOBAL_COMMUNICATOR,&status);


                if (fnameOutput!=NULL) {
                  MPI_Recv(&SingleVectorCheckSum.crc_accum,1,MPI_UNSIGNED_LONG,startNode->Thread,0,MPI_GLOBAL_COMMUNICATOR,&status);

                  fprintf(fout,"node: id=%ld, i=%i, j=%i, k=%i, CheckSum=0x%lx\n",startNode->Temp_ID,i,j,k,SingleVectorCheckSum.checksum());
                }
              }
              else {
                MPI_Status status;
                MPI_Recv(&CheckSum.crc_accum,1,MPI_UNSIGNED_LONG,0,0,MPI_GLOBAL_COMMUNICATOR,&status);

                if (block!=NULL) if ((CenterNode=block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k)))!=NULL) {
                  CheckSum.add(CenterNode->GetAssociatedDataBufferPointer()+SampleVectorOffset,SampleVectorLength);

                  if (fnameOutput!=NULL) {
                    SingleVectorCheckSum.add(CenterNode->GetAssociatedDataBufferPointer()+SampleVectorOffset,SampleVectorLength);
                  }
                }

                MPI_Send(&CheckSum.crc_accum,1,MPI_UNSIGNED_LONG,0,0,MPI_GLOBAL_COMMUNICATOR);

                if (fnameOutput!=NULL) {
                  MPI_Send(&SingleVectorCheckSum.crc_accum,1,MPI_UNSIGNED_LONG,0,0,MPI_GLOBAL_COMMUNICATOR);
                }

              }
            }

          }

        }
      }
    }
  }
  else {
    //add daugher blocks
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *downNode;

    for (i=0;i<(1<<DIM);i++) if ((downNode=startNode->downNode[i])!=NULL) SaveCenterNodeAssociatedDataSignature(SampleVectorOffset,SampleVectorLength,nline,fnameSource,fnameOutput,downNode);
  }

  //output the checksum
  if (startNode==PIC::Mesh::mesh.rootTree) {
    pipe.flush();

    if (PIC::ThisThread==0) {
      char msg[500];

      sprintf(msg," line=%ld, file=%s (Call Counter=%i)",nline,fnameSource,nCallCounter);
      CheckSum.PrintChecksumSingleThread(msg);

      if (fnameOutput!=NULL) {
        fclose(fout);
        fout=NULL;
      }
    }
  }

  return CheckSum.checksum();
}


//=====================================================================================
//get signature describe the particle population
unsigned long int PIC::Debugger::GetParticlePopulationSignature(long int nline,const char* fname) {
  CRC32 Checksum;
  PIC::ParticleBuffer::byte *ParticleDataPtr,ParticleBuffer[PIC::ParticleBuffer::ParticleDataLength];
  int i,j,k,ptr;

  //init the particle buffer
  for (i=0;i<PIC::ParticleBuffer::ParticleDataLength;i++) ParticleBuffer[i]=0;

  //loop through all blocks
  for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=PIC::Mesh::mesh.BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) {
    if ((node->Thread==PIC::ThisThread)||(PIC::ThisThread==0)) {

      if (node->Thread!=0) {
        //recieve the checksum object
        if (PIC::ThisThread==0) {
          MPI_Send(&Checksum,sizeof(Checksum),MPI_CHAR,node->Thread,0,MPI_GLOBAL_COMMUNICATOR);
        }
        else {
          MPI_Status status;

          MPI_Recv(&Checksum,sizeof(Checksum),MPI_CHAR,0,0,MPI_GLOBAL_COMMUNICATOR,&status);
        }
      }

      if (node->block!=NULL) for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
        ptr=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

        //collect signature
        while (ptr!=-1) {
          //copy the state vector of the particle without 'next' and 'prev'
          ParticleDataPtr=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
          PIC::ParticleBuffer::CloneParticle(ParticleBuffer,ParticleDataPtr);

          //add signature of the particle
          Checksum.add(ParticleBuffer,PIC::ParticleBuffer::ParticleDataLength);
          ptr=PIC::ParticleBuffer::GetNext(ptr);
        }
      }

      //send the checksum object back to the root
      if (node->Thread!=0) {
        //recieve the checksum object
        if (PIC::ThisThread!=0) {
          MPI_Send(&Checksum,sizeof(Checksum),MPI_CHAR,0,0,MPI_GLOBAL_COMMUNICATOR);
        }
        else {
          MPI_Status status;

          MPI_Recv(&Checksum,sizeof(Checksum),MPI_CHAR,node->Thread,0,MPI_GLOBAL_COMMUNICATOR,&status);
        }
      }
    }
  }


  //output the checksum
  if (PIC::ThisThread==0) {
    char msg[500];

    sprintf(msg," line=%ld, file=%s",nline,fname);
    Checksum.PrintChecksumSingleThread(msg);
  }

  return Checksum.checksum();
}


unsigned long int PIC::Debugger::GetParticlePopulationStateVectorSignature(int offset,int length,long int nline,const char* fname) {
  CRC32 Checksum;
  PIC::ParticleBuffer::byte *ParticleDataPtr,ParticleBuffer[PIC::ParticleBuffer::ParticleDataLength];
  int i,j,k,ptr;

  CMPI_channel pipe;
  pipe.init(1000000);

  //init the particle buffer
  for (i=0;i<PIC::ParticleBuffer::ParticleDataLength;i++) ParticleBuffer[i]=0;

  const int CommunicationCompleted_SIGNAL=0;
  const int ParticleDataSend_SIGNAL=1;
  const int BockDataStarted_SIGNAL=2;

  if (PIC::ThisThread==0) pipe.openRecvAll();
  else pipe.openSend(0);

  //loop through all blocks
  for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=PIC::Mesh::mesh.BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) {
    if ((node->Thread==0)&&(PIC::ThisThread==0)) {
      //the block belongs to the root

      if (node->block!=NULL) for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
        ptr=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

        //collect signature
        while (ptr!=-1) {
          //copy the state vector of the particle without 'next' and 'prev'
          ParticleDataPtr=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
          PIC::ParticleBuffer::CloneParticle(ParticleBuffer,ParticleDataPtr);

          //add signature of the particle
          Checksum.add(ParticleBuffer+offset,length);
          ptr=PIC::ParticleBuffer::GetNext(ptr);
        }
      }
    }
    else if (PIC::ThisThread==0) {
      //this is the root BUT the block belongs to another MPI process
      int Signal;

      pipe.recv(Signal,node->Thread);

      switch (Signal) {
      case BockDataStarted_SIGNAL:
        for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
          pipe.recv(Signal,node->Thread);

          //collect signatures
          while (Signal!=CommunicationCompleted_SIGNAL) {
            pipe.recv(ParticleBuffer,PIC::ParticleBuffer::ParticleDataLength,node->Thread);
            Checksum.add(ParticleBuffer+offset,length);

            pipe.recv(Signal,node->Thread);
          }
        }

        break;
      case CommunicationCompleted_SIGNAL:
        break;
      default:
        exit(__LINE__,__FILE__,"Error: the sigmal is not recognized");
      }


    }
    else if (node->Thread==PIC::ThisThread) {
      //this is NOT the root BUT the block belongs to the current MPI process
      //loop through all cells and particles

      if (node->block!=NULL) {
        pipe.send(BockDataStarted_SIGNAL);

        for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
          ptr=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

          while (ptr!=-1) {
            pipe.send(ParticleDataSend_SIGNAL);

            //copy the state vector of the particle without 'next' and 'prev'
            ParticleDataPtr=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
            PIC::ParticleBuffer::CloneParticle(ParticleBuffer,ParticleDataPtr);
            pipe.send(ParticleBuffer,PIC::ParticleBuffer::ParticleDataLength);

            ptr=PIC::ParticleBuffer::GetNext(ptr);
          }

          pipe.send(CommunicationCompleted_SIGNAL);
        }

      }
      else {
        pipe.send(CommunicationCompleted_SIGNAL);
      }

    }
  }

  //output the checksum
  if (PIC::ThisThread==0) {
    char msg[500];

    pipe.closeRecvAll();

    sprintf(msg," line=%ld, file=%s",nline,fname);
    Checksum.PrintChecksumSingleThread(msg);
  }
  else pipe.closeSend();

  return Checksum.checksum();
}

unsigned long int PIC::Debugger::GetParticlePopulationLocationSignature(long int nline,const char* fname) {
  return GetParticlePopulationStateVectorSignature(_PIC_PARTICLE_DATA__POSITION_OFFSET_,DIM*sizeof(double),nline,fname);
}
unsigned long int PIC::Debugger::GetParticlePopulationVelocitySignature(long int nline,const char* fname) {
  return GetParticlePopulationStateVectorSignature(_PIC_PARTICLE_DATA__VELOCITY_OFFSET_,3*sizeof(double),nline,fname);
}

//=========================================================================================================
//save the map of the domain decomposition
void PIC::Debugger::SaveDomainDecompositionMap(long int nline,const char* fname,int Index) {
  FILE *fout;
  char FullFileName[100];
  int id,i,j,iface,iedge,icorner;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *neibNode;

  sprintf(FullFileName,"DomainDecompositionMap.thread=%i(line=%ld,file=%s,Index=%i).dat",PIC::ThisThread,nline,fname,Index);
  fout=fopen(FullFileName,"w");

  fprintf(fout,"VARIABLES=\"NodeTempID\", \"Thread\", \"Face Neib\", \"Edge Neib\", \"Corner Neib\"\n");

  for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=PIC::Mesh::mesh.BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) {
    fprintf(fout,"node->Temp_ID=%ld, thread=%i\n",node->Temp_ID,node->Thread);

    //face neib
    for (iface=0;iface<6;iface++) for (i=0;i<2;i++) for (j=0;j<2;j++) {
      if ((neibNode=node->GetNeibFace(iface,i,j))!=NULL) id=neibNode->Temp_ID;
      else id=-1;

      fprintf(fout,"iface=%i,i=%i,j=%i,neib=%i\n",iface,i,j,id);
    }

    //edge neib
    for (iedge=0;iedge<12;iedge++) for (i=0;i<2;i++) {
      if ((neibNode=node->GetNeibEdge(iedge,i))!=NULL) id=neibNode->Temp_ID;
      else id=-1;

      fprintf(fout,"iedge=%i,i=%i,neib=%i\n",iface,i,id);
    }

    //corner
    for (icorner=0;icorner<8;icorner++) {
      if ((neibNode=node->GetNeibCorner(icorner))!=NULL) id=neibNode->Temp_ID;
      else id=-1;

      fprintf(fout,"icorner=%i,neib=%i\n",icorner,id);
    }

    //end the line
    fprintf(fout,"\n");
  }

  //close the file
  fclose(fout);
}


//====================================================================================================
//output the debug debug particle data
list<PIC::Debugger::ParticleDebugData::cDebugData> PIC::Debugger::ParticleDebugData::DebugParticleData;

//method for sorting the debug particle list
bool PIC::Debugger::ParticleDebugData::CompareParticleDebugData(const PIC::Debugger::ParticleDebugData::cDebugData& first, const PIC::Debugger::ParticleDebugData::cDebugData& second) {
  return (first.initCheckSum < second.initCheckSum);
}

//accumulate the partilce debug data
void PIC::Debugger::ParticleDebugData::AddParticleDebugData(long int ptr,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node,bool InitChckSumMode) {
  int i,j,k,di,dj,dk;
  double *x;
  cDebugData p;

  //add the particle data
  CRC32 CheckSum;

  if (InitChckSumMode==true) {
    p.initCheckSum=PIC::ParticleBuffer::GetParticleSignature(ptr);

    x=PIC::ParticleBuffer::GetX(ptr);
    PIC::Mesh::mesh.fingCellIndex(x,i,j,k,node,false);

    p.i=i;
    p.j=j;
    p.k=k;
    p.nodeid=node->Temp_ID;
    p.ptr=ptr;
  }
  else {
    //look for the particle
    list<cDebugData>::iterator pp;

    for (pp=DebugParticleData.begin();pp!=DebugParticleData.end();pp++) {
      if (pp->ptr==ptr) {
        //the particle has been found
        pp->finalCheckSum=PIC::ParticleBuffer::GetParticleSignature(ptr);
        return;
      }
    }

    //when come to this point --> the particle has not beed found
    exit(__LINE__,__FILE__,"the particle has not been found");
  }

  //add checksum of the corner and center nodes
  for (dk=0;dk<2;dk++) {
    for (dj=0;dj<2;dj++)  {
      for (di=0;di<2;di++) {
        PIC::Mesh::cDataCornerNode *CornerNode;

        if (node->Thread==PIC::ThisThread) {
          if (node->block!=NULL) {
            CheckSum.clear();

            CornerNode=node->block->GetCornerNode(_getCornerNodeLocalNumber(p.i+di,p.j+dj,p.k+dk));
            CheckSum.add(CornerNode->GetAssociatedDataBufferPointer(),PIC::Mesh::cDataCornerNode::totalAssociatedDataLength);

            p.CornerNodeChecksum[di][dj][dk]=CheckSum.checksum();
          }
        }
      }
    }
  }


  for (dk=-1;dk<1;dk++) {
    for (dj=-1;dj<1;dj++)  {
      for (di=-1;di<1;di++) {
        PIC::Mesh::cDataCenterNode *CenterNode;

        if (node->Thread==PIC::ThisThread) {
          if (node->block!=NULL) {
            CheckSum.clear();

            CenterNode=node->block->GetCenterNode(_getCenterNodeLocalNumber(p.i+di,p.j+dj,p.k+dk));
            CheckSum.add(CenterNode->GetAssociatedDataBufferPointer(),PIC::Mesh::cDataCenterNode::totalAssociatedDataLength);

            p.CenterNodeChecksum[1+di][1+dj][1+dk]=CheckSum.checksum();
          }
        }
      }
    }
  }


  //add the particle data to the list
  DebugParticleData.push_front(p);
}

//============================================================================================
//output particle debug data
void PIC::Debugger::ParticleDebugData::OutputParticleDebugData(int nLineSource,const char *fNameSource,int Index) {
  int size;
  cDebugData p;
  list<cDebugData>::iterator pp;

  CMPI_channel pipe;
  pipe.init(1000000);

  if (PIC::ThisThread==0) {
      pipe.openRecvAll();
  }
  else pipe.openSend(0);

  //collect the debug data from all MPI processes
  if (PIC::ThisThread==0) {
    for (int thread=1;thread<PIC::nTotalThreads;thread++) {
      pipe.recv(size,thread);

      for (int i=0;i<size;i++) {
        pipe.recv(p,thread);
        DebugParticleData.push_front(p);
      }
    }


    //sort and output the particle data
    DebugParticleData.sort(CompareParticleDebugData);

    //output the debug data
    FILE *fout;
    char fname[100];
    int di,dj,dk;

    sprintf(fname,"ParticleDebugData(nline=%i,file=%s).Index=%i.dat",nLineSource,fNameSource,Index);
    fout=fopen(fname,"w");

    for (pp=DebugParticleData.begin();pp!=DebugParticleData.end();pp++) {
      fprintf(fout,"particle(init)=0x%lx, particle(final)=0x%lx, nodeid=%i, i=%i, j=%i, k=%i\n",
        pp->initCheckSum,pp->finalCheckSum,pp->nodeid,pp->i,pp->j,pp->k);

      fprintf(fout,"CenterNodeChecksum:\n");

      for (di=0;di<2;di++) for (dj=0;dj<2;dj++) for (dk=0;dk<2;dk++) {
        fprintf(fout,"di=%i, dj=%i, dk=%i, CenterNodeChecksum=0x%lx\n",di-1,dj-1,dk-1,pp->CenterNodeChecksum[di][dj][dk]);
      }


      fprintf(fout,"CornerNodeChecksum:\n");

      for (di=0;di<2;di++) for (dj=0;dj<2;dj++) for (dk=0;dk<2;dk++) {
        fprintf(fout,"di=%i, dj=%i, dk=%i, CornerNodeChecksum=0x%lx\n",di,dj,dk,pp->CornerNodeChecksum[di][dj][dk]);
      }

      fprintf(fout,"\n");
    }

    fclose(fout);
  }
  else {
    size=DebugParticleData.size();
    pipe.send(size);

    for (pp=DebugParticleData.begin();pp!=DebugParticleData.end();pp++) {
      p=*pp;
      pipe.send(p);
    }
  }


   //cloe the pipe
   if (PIC::ThisThread==0) {
     pipe.closeRecvAll();
   }
   else pipe.closeSend();

   //remove the particle debug data list
   DebugParticleData.clear();
}

//========================================================================================================================
//save signatures of the nodes
void PIC::Debugger::SaveNodeSignature(int nline,const char *fname) {
  int i,j,k;
  static FILE *fout=NULL;
  static int ncall=0;

  ncall++;

  if ((fout==NULL)&&(PIC::ThisThread==0))  {
    fout=fopen("SaveNodeSignature.dat","w");
  }

  unsigned long s=PIC::Debugger::GetCornerNodeAssociatedDataSignature(nline,fname);
  unsigned long int p=PIC::Debugger::GetParticlePopulationSignature(nline,fname);

  //'corner' data
  for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=PIC::Mesh::mesh.BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) {
    for (k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_+1;k++) {
      for (j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_+1;j++)  {
        for (i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_+1;i++) {
          unsigned long int cs;
          CRC32 CheckSum;
          PIC::Mesh::cDataCornerNode *CornerNode;

          if (node->Thread==PIC::ThisThread) {
            if (node->block!=NULL) {
              CornerNode=node->block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k));
              CheckSum.add(CornerNode->GetAssociatedDataBufferPointer(),PIC::Mesh::cDataCornerNode::totalAssociatedDataLength);
            }

            cs=CheckSum.checksum();

            if (PIC::ThisThread!=0) {
              //send the checksum to the root
              MPI_Send(&cs,1,MPI_UNSIGNED_LONG,0,0,MPI_GLOBAL_COMMUNICATOR);
            }
          }
          else if (PIC::ThisThread==0) {
            //recieve the checksum
            MPI_Status status;

            MPI_Recv(&cs,1,MPI_UNSIGNED_LONG,node->Thread,0,MPI_GLOBAL_COMMUNICATOR,&status);
          }

          //print the checksum to a file
          if (PIC::ThisThread==0) {
            fprintf(fout,"Corner CheckSum=0x%lx, i=%i, j=%i, k=%i,id=%ld, ncall=%i, line=%i,file=%s \n",cs,i,j,k,node->Temp_ID,ncall,nline,fname);
          }
        }
      }
    }
  }

  //'center' data
  for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=PIC::Mesh::mesh.BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) {
    for (k=-_GHOST_CELLS_Z_;k<_BLOCK_CELLS_Z_+_GHOST_CELLS_Z_;k++) {
      for (j=-_GHOST_CELLS_Y_;j<_BLOCK_CELLS_Y_+_GHOST_CELLS_Y_;j++)  {
         for (i=-_GHOST_CELLS_X_;i<_BLOCK_CELLS_X_+_GHOST_CELLS_X_;i++) {
           unsigned long int cs;
           CRC32 CheckSum;
           PIC::Mesh::cDataCenterNode *CenterNode;

           if (node->Thread==PIC::ThisThread) {
             if (node->block!=NULL) {
               CenterNode=node->block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k));
               CheckSum.add(CenterNode->GetAssociatedDataBufferPointer(),PIC::Mesh::cDataCenterNode::totalAssociatedDataLength);
             }

             cs=CheckSum.checksum();

             if (PIC::ThisThread!=0) {
               //send the checksum to the root
               MPI_Send(&cs,1,MPI_UNSIGNED_LONG,0,0,MPI_GLOBAL_COMMUNICATOR);
             }
           }
           else if (PIC::ThisThread==0) {
             //recieve the checksum
             MPI_Status status;

             MPI_Recv(&cs,1,MPI_UNSIGNED_LONG,node->Thread,0,MPI_GLOBAL_COMMUNICATOR,&status);
           }

           //print the checksum to a file
           if (PIC::ThisThread==0) {
             fprintf(fout,"Center CheckSum=0x%lx, i=%i, j=%i, k=%i,id=%ld, ncall=%i, line=%i,file=%s \n",cs,i,j,k,node->Temp_ID,ncall,nline,fname);
           }
        }
      }
    }
  }

  if (fout!=NULL) {
    fflush(fout);
  }
}


double PIC::Debugger::read_mem_usage() {
  // This function returns the resident set size (RSS) of                                                                                                                             
  // this processor in unit MB.                                                                                                                                                       

  // From wiki:                                                                                                                                                                       
  // RSS is the portion of memory occupied by a process that is                                                                                                                       
  // held in main memory (RAM).                                                                                                                                                       

  double rssMB = 0.0;

  ifstream stat_stream("/proc/self/stat", ios_base::in);

  if (!stat_stream.fail()) {
    // Dummy vars for leading entries in stat that we don't care about                                                                                                                
    string pid, comm, state, ppid, pgrp, session, tty_nr;
    string tpgid, flags, minflt, cminflt, majflt, cmajflt;
    string utime, stime, cutime, cstime, priority, nice;
    string O, itrealvalue, starttime;

    // Two values we want                                                                                                                                                             
    unsigned long vsize;
    unsigned long rss;

    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr >>
      tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt >> utime >>
      stime >> cutime >> cstime >> priority >> nice >> O >> itrealvalue >>
      starttime >> vsize >> rss; // Ignore the rest                                                                                                                                   
    stat_stream.close();

    rssMB = rss*sysconf(_SC_PAGE_SIZE)/1024.0/1024.0;
  }

  return rssMB;
}


#if defined(__linux__)
struct mallinfo PIC::Debugger::MemoryLeakCatch::Baseline;
bool PIC::Debugger::MemoryLeakCatch::Active=false;

void PIC::Debugger::MemoryLeakCatch::SetBaseline() {Baseline=mallinfo();} 

void PIC::Debugger::MemoryLeakCatch::SetActive(bool flag) {
  Active=flag;
  if (Active==true) SetBaseline();
}

void PIC::Debugger::MemoryLeakCatch::Trap(int nline,const char* fname) {
   char msg[500];
   struct mallinfo t=mallinfo();

   sprintf(msg,"The amount of memory allocated with malloc have increased by %e MB (thread=%i,line=%i,file=%s)\n",(t.uordblks-Baseline.uordblks)/1048576.0,PIC::ThisThread,nline,fname);
   printf(msg); //here the memory leack can be intersepted in a debugger
}  

bool PIC::Debugger::MemoryLeakCatch::Test(int nline,const char* fname) {
  bool res=false;

  if (Active==true) {
    struct mallinfo t=mallinfo();

    if (t.uordblks>Baseline.uordblks) {
      Trap(nline,fname);
      SetBaseline();
     
      res=true; 
    }
  }

  return res;
}
#else  //defined(__linux__) 
double PIC::Debugger::MemoryLeakCatch::Baseline=0.0;
bool PIC::Debugger::MemoryLeakCatch::Active=false;

void PIC::Debugger::MemoryLeakCatch::SetBaseline() {Baseline=PIC::Debugger::read_mem_usage();}

void PIC::Debugger::MemoryLeakCatch::SetActive(bool flag) {
  Active=flag;
  if (Active==true) SetBaseline();
}

void PIC::Debugger::MemoryLeakCatch::Trap(int nline,const char* fname) {
   char msg[500];

   sprintf(msg,"The size of the allocated  heap have increased by %e MB (thread=%i,line=%i,file=%s)\n",PIC::Debugger::read_mem_usage()-Baseline,PIC::ThisThread,nline,fname);
   printf(msg); //here the memory leack can be intersepted in a debugger
}

bool PIC::Debugger::MemoryLeakCatch::Test(int nline,const char* fname) {
  bool res=false;

  if (Active==true) {
    if (PIC::Debugger::read_mem_usage()>Baseline) {
      Trap(nline,fname);
      SetBaseline();

      res=true;
    }
  }

  return res;
}
#endif //defined(__linux__)


void PIC::Debugger::check_max_mem_usage(string tag) {
  double memLocal = read_mem_usage();
    
  cout << "$PREFIX: " << tag << " Maximum memory usage = " << memLocal << "Mb(MB?) on rank = " << PIC::ThisThread << endl;
}

void PIC::Debugger::GetMemoryUsageStatus(long int nline,const char *fname,bool ShowUsagePerProcessFlag) {
  double LocalMemoryUsage,GlobalMemoryUsage;
  double *MemoryUsageTable=NULL;

  if (PIC::ThisThread==0) MemoryUsageTable=new double [PIC::nTotalThreads];

  //collect momery usage information
  LocalMemoryUsage=read_mem_usage();

  //gather the memory usage table
  MPI_Gather(&LocalMemoryUsage,1,MPI_DOUBLE,MemoryUsageTable,1,MPI_DOUBLE,0,MPI_GLOBAL_COMMUNICATOR);

  //output the memory usage status
  if (PIC::ThisThread==0) {
    int thread;

    if (ShowUsagePerProcessFlag==true) {
      printf("$PREFIX: Memory Usage Status (file=%s,line=%ld)\nThread\tUsed Memory (MB)\n",fname,nline);

      for (thread=0,GlobalMemoryUsage=0.0;thread<PIC::nTotalThreads;thread++) {
        GlobalMemoryUsage+=MemoryUsageTable[thread];
        printf("$PREFIX: %i\t%e [MB]\n",thread,MemoryUsageTable[thread]);
      }

      printf("$PREFIX: Total=%e [MB], %e [MB per MPI Process]\n",GlobalMemoryUsage,GlobalMemoryUsage/PIC::nTotalThreads);
    }
    else {
      for (thread=0,GlobalMemoryUsage=0.0;thread<PIC::nTotalThreads;thread++) GlobalMemoryUsage+=MemoryUsageTable[thread];

      printf("$PREFIX: Memory Usage Status (file=%s,line=%ld): Total Memory Used=%e [MB], %e [MB per MPI Process]\n",fname,nline,GlobalMemoryUsage,GlobalMemoryUsage/PIC::nTotalThreads);
    }

    delete [] MemoryUsageTable;
  }
}

//=======================================================================================
//verify that the number of particles in the lists is the same as the number of used particles in the buffer

int PIC::Debugger::GetParticleNumberInLists(bool CurrentThreadOnly) {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  int ptr,i,j,k,nTotalParticles=0;

  for (node=PIC::Mesh::mesh.BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) if ((node->block!=NULL) && ((CurrentThreadOnly==false)||(node->Thread==PIC::ThisThread)) ) {
    for (k=0;k<_BLOCK_CELLS_Z_;k++) {
      for (j=0;j<_BLOCK_CELLS_Y_;j++) {
        for (i=0;i<_BLOCK_CELLS_X_;i++) {
          ptr=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

          while (ptr!=-1) {
            nTotalParticles++;
            ptr=PIC::ParticleBuffer::GetNext(ptr);
          }
        }
      }
    }
  }

  return nTotalParticles;
}


void PIC::Debugger::VerifyTotalParticleNumber(int line,const char* fname,bool CurrentThreadOnly) {
  int nTotalParticles=GetParticleNumberInLists(CurrentThreadOnly);


  if (nTotalParticles!=PIC::ParticleBuffer::GetAllPartNum()) {
    char msg[1000];

    sprintf(msg,"Error: the particle number is inconsistent (particles in the lists: %i, particle buffer: %ld",nTotalParticles,PIC::ParticleBuffer::GetAllPartNum());
    exit(line,fname,msg);
  }
}


//=======================================================================================
//ger check summs of the corner abd center nodes in all blocks excluding the ghost blocks

void PIC::Debugger:: GetBlockAssociatedDataSignature_no_ghost_blocks(long int nline,const char* fname) {
  CRC32 CenterNodeCheckSum,CornerNodeCheckSum;


  //reserve/release the flags
  int periodic_bc_pair_real_block=-1;
  int periodic_bc_pair_ghost_block=-1;

  auto ReleasePeriodicBCFlags = [&] () {
    PIC::Mesh::mesh.rootTree->ReleaseFlag(periodic_bc_pair_real_block);
    PIC::Mesh::mesh.rootTree->ReleaseFlag(periodic_bc_pair_ghost_block);
  };


  auto ReservePeriodicBCFlags = [&] () {
    periodic_bc_pair_real_block=PIC::Mesh::mesh.rootTree->CheckoutFlag();
    periodic_bc_pair_ghost_block=PIC::Mesh::mesh.rootTree->CheckoutFlag();

    if ((periodic_bc_pair_real_block==-1)||(periodic_bc_pair_ghost_block==-1)) exit(__LINE__,__FILE__,"Error: cannot reserve a flag");
  };


  //set the flag
  std::function<void(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*)> ResetPeriodicBCFlags;

  ResetPeriodicBCFlags = [&] (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) -> void {
    startNode->SetFlag(false,periodic_bc_pair_real_block);
    startNode->SetFlag(false,periodic_bc_pair_ghost_block);

    int i;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *downNode;

    for (i=0;i<(1<<DIM);i++) if ((downNode=startNode->downNode[i])!=NULL) {
      ResetPeriodicBCFlags(downNode);
    }

    if (startNode==PIC::Mesh::mesh.rootTree) {
      //loop throught the list of the block pairs used to impose the periodic BC
      int iBlockPair,RealBlockThread,GhostBlockThread;
      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *RealBlock,*GhostBlock;

      //loop through all block pairs
      for (iBlockPair=0;iBlockPair<PIC::BC::ExternalBoundary::Periodic::BlockPairTableLength;iBlockPair++) {
        GhostBlock=PIC::BC::ExternalBoundary::Periodic::BlockPairTable[iBlockPair].GhostBlock;
        RealBlock=PIC::BC::ExternalBoundary::Periodic::BlockPairTable[iBlockPair].RealBlock;

        GhostBlock->SetFlag(true,periodic_bc_pair_ghost_block);
        RealBlock->SetFlag(true,periodic_bc_pair_real_block);
      }
    }
  };


  //get the checksum of a blocks
  auto CenterNodeBlockCheckSum = [&] (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
    int i,j,k;
    PIC::Mesh::cDataBlockAMR *block=node->block;

    if (node->Thread!=0) {
      if (PIC::ThisThread==0) {
        MPI_Send(&CenterNodeCheckSum,sizeof(CenterNodeCheckSum),MPI_CHAR,node->Thread,0,MPI_GLOBAL_COMMUNICATOR);
      }
      else if (PIC::ThisThread==node->Thread) {
        MPI_Status status;

        MPI_Recv(&CenterNodeCheckSum,sizeof(CenterNodeCheckSum),MPI_CHAR,0,0,MPI_GLOBAL_COMMUNICATOR,&status);
      }
    }


    if ((node->Thread==PIC::ThisThread)&&(block!=NULL)) {
      for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
        PIC::Mesh::cDataCenterNode *CenterNode;

        CenterNode=block->GetCenterNode(_getCenterNodeLocalNumber(i,j,k));
        if (CenterNode!=NULL) CenterNodeCheckSum.add(CenterNode->GetAssociatedDataBufferPointer(),PIC::Mesh::cDataCenterNode::totalAssociatedDataLength);
      }
    }

    if (node->Thread!=0) {
      if (PIC::ThisThread==node->Thread) {
        MPI_Send(&CenterNodeCheckSum,sizeof(CenterNodeCheckSum),MPI_CHAR,0,0,MPI_GLOBAL_COMMUNICATOR);
      }
      else if (PIC::ThisThread==0) {
        MPI_Status status;

        MPI_Recv(&CenterNodeCheckSum,sizeof(CenterNodeCheckSum),MPI_CHAR,node->Thread,0,MPI_GLOBAL_COMMUNICATOR,&status);
      }
    }
  };

  auto CornerNodeBlockCheckSum = [&] (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
    int i,j,k;
    PIC::Mesh::cDataBlockAMR *block=node->block;

    if (node->Thread!=0) {
      if (PIC::ThisThread==0) {
        MPI_Send(&CornerNodeCheckSum,sizeof(CornerNodeCheckSum),MPI_CHAR,node->Thread,0,MPI_GLOBAL_COMMUNICATOR);
      }
      else if (PIC::ThisThread==node->Thread) {
        MPI_Status status;

        MPI_Recv(&CornerNodeCheckSum,sizeof(CornerNodeCheckSum),MPI_CHAR,0,0,MPI_GLOBAL_COMMUNICATOR,&status);
      }
    }


    if ((node->Thread==PIC::ThisThread)&&(block!=NULL)) {
      for (k=0;k<_BLOCK_CELLS_Z_+1;k++) for (j=0;j<_BLOCK_CELLS_Y_+1;j++) for (i=0;i<_BLOCK_CELLS_X_+1;i++) {
        PIC::Mesh::cDataCornerNode *CornerNode;

        CornerNode=block->GetCornerNode(_getCornerNodeLocalNumber(i,j,k));
        if (CornerNode!=NULL) CornerNodeCheckSum.add(CornerNode->GetAssociatedDataBufferPointer(),PIC::Mesh::cDataCornerNode::totalAssociatedDataLength);
      }
    }

    if (node->Thread!=0) {
      if (PIC::ThisThread==node->Thread)  {
        MPI_Send(&CornerNodeCheckSum,sizeof(CornerNodeCheckSum),MPI_CHAR,0,0,MPI_GLOBAL_COMMUNICATOR);
      }
      else if (PIC::ThisThread==0) {
        MPI_Status status;

        MPI_Recv(&CornerNodeCheckSum,sizeof(CornerNodeCheckSum),MPI_CHAR,node->Thread,0,MPI_GLOBAL_COMMUNICATOR,&status);
      }
    }
  };




  //loop through all blocks
  ReservePeriodicBCFlags();
  ResetPeriodicBCFlags(PIC::Mesh::mesh.rootTree);


  for (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=PIC::Mesh::mesh.BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) {
    if (node->TestFlag(periodic_bc_pair_ghost_block)==false) {
      CenterNodeBlockCheckSum(node);
      CornerNodeBlockCheckSum(node);
    }
  }


  ReleasePeriodicBCFlags();


  //output the checksum
  if (PIC::ThisThread==0) {
    char msg[1000];

    static int cnt=0;

    sprintf(msg,"CenterNode CheckSum (%s@%ld, cnt=%i)",fname,nline,cnt);
    CenterNodeCheckSum.PrintChecksumSingleThread(msg);

    sprintf(msg,"CornerNode CheckSum (%s@%ld, cnt=%i)",fname,nline,cnt);
    CornerNodeCheckSum.PrintChecksumSingleThread(msg);

    cnt++;
  }
}

//=======================================================================================================================================
//order particle lists
void PIC::Debugger::OrderParticleList(long int &first_particle) {
  long int ptr=first_particle;

  class cLocalParticleData {
  public:
    long int ptr;
    unsigned long int checksum;

    bool operator < (const cLocalParticleData& __y) const {return (checksum < __y.checksum);}
  };

  list <cLocalParticleData> ParticleData;
  cLocalParticleData p_data;

  CRC32 p_checksum;
  PIC::ParticleBuffer::byte ParticleInternalData[PIC::ParticleBuffer::ParticleDataLength];
  PIC::ParticleBuffer::byte *p_data_ptr;

  if (first_particle==-1) return;

  for (int i=0;i<PIC::ParticleBuffer::ParticleDataLength;i++) ParticleInternalData[i]=0;

  //create the list
  while (ptr!=-1) {
    p_data_ptr=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
    p_checksum.clear();


    PIC::ParticleBuffer::CloneParticle(ParticleInternalData,p_data_ptr);
    p_checksum.add(ParticleInternalData,PIC::ParticleBuffer::ParticleDataLength);


    p_data.ptr=ptr;
    p_data.checksum=p_checksum.crc_accum;

    ParticleData.push_front(p_data);
    ptr=PIC::ParticleBuffer::GetNext(ptr);
  }

  //sort the list
  ParticleData.sort();

  //create the ordered particle list
  first_particle=-1;

  for (list <cLocalParticleData>::iterator it=ParticleData.begin();it!=ParticleData.end();it++) {
    ptr=it->ptr;

    PIC::ParticleBuffer::SetNext(first_particle,ptr);
    PIC::ParticleBuffer::SetPrev(-1,ptr);

    if (first_particle>=0) PIC::ParticleBuffer::SetPrev(ptr,first_particle);

    first_particle=ptr;
  }
}


void PIC::Debugger::OrderParticleLists() {
  std::function<void(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*)> ProcessBlock;

  ProcessBlock=[&] (cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *Node) -> void {
    if (Node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
      PIC::Mesh::cDataBlockAMR *block=Node->block;

      if (block!=NULL) {
         for (int k=0;k<_BLOCK_CELLS_Z_;k++) {
           for (int j=0;j<_BLOCK_CELLS_Y_;j++) {
             for (int i=0;i<_BLOCK_CELLS_X_;i++) {
               if (block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]!=-1) {
                 OrderParticleList(block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]);
               }
             }
           }
         }
      }
    }
    else {
      //add daugher blocks
      cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>  *downNode;

      for (int i=0;i<(1<<DIM);i++) if ((downNode=Node->downNode[i])!=NULL) ProcessBlock(downNode);
    }
  };

  ProcessBlock(PIC::Mesh::mesh.rootTree);
}















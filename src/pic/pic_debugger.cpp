/*
 * pic_debugger.cpp
 *
 *  Created on: Feb 19, 2014
 *      Author: vtenishe
 */
//$Id$
//contains functions used for debugging AMPS

#include "pic.h"

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
      ptr=*(startNode->block->GetTempParticleMovingListTableThread(thread_OpenMP,i,j,k));
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
        ptr=(*startNode->block->GetTempParticleMovingListTableThread(thread_OpenMP,i,j,k));
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
      LocalCellNumber=PIC::Mesh::mesh.getCenterNodeLocalNumber(i,j,k);
      cell=block->GetCenterNode(LocalCellNumber);

      SamplingBuffer=cell->GetAssociatedDataBufferPointer()+PIC::Mesh::collectingCellSampleDataPointerOffset;


      for (s=0;s<PIC::nTotalSpecies;s++) {
        PIC::Debugger::CatchOutLimitValue((s+(double*)(SamplingBuffer+PIC::Mesh::sampledParticleWeghtRelativeOffset)),1,__LINE__,__FILE__);
        PIC::Debugger::CatchOutLimitValue((s+(double*)(SamplingBuffer+PIC::Mesh::sampledParticleNumberRelativeOffset)),1,__LINE__,__FILE__);
        PIC::Debugger::CatchOutLimitValue((s+(double*)(SamplingBuffer+PIC::Mesh::sampledParticleNumberDensityRelativeOffset)),1,__LINE__,__FILE__);

        PIC::Debugger::CatchOutLimitValue((3*s+(double*)(SamplingBuffer+PIC::Mesh::sampledParticleVelocityRelativeOffset)),DIM,__LINE__,__FILE__);
        PIC::Debugger::CatchOutLimitValue((3*s+(double*)(SamplingBuffer+PIC::Mesh::sampledParticleVelocity2RelativeOffset)),DIM,__LINE__,__FILE__);
        PIC::Debugger::CatchOutLimitValue((s+(double*)(SamplingBuffer+PIC::Mesh::sampledParticleSpeedRelativeOffset)),1,__LINE__,__FILE__);
      }

    }
  }

#endif
#endif
}



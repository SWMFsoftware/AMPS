//random number generator


#include "pic.h"

/*
 * pic_rnd.cpp
 *
 *  Created on: Jul 18, 2019
 *      Author: vtenishe
 */

int PIC::Rnd::CenterNode::Offset=-1;

int PIC::Rnd::CenterNode::RequestDataBuffer(int OffsetIn) {
  Offset=OffsetIn;

  #if _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
  int nThreadsOpenMP=1;
  #else
  int nThreadsOpenMP=omp_get_num_threads();
  #endif

  return nThreadsOpenMP*sizeof(unsigned long int);
}

void PIC::Rnd::CenterNode::Init() {
  //reserve memory to store the seed in the center node state vector
  PIC::IndividualModelSampling::RequestStaticCellData.push_back(RequestDataBuffer);
}

void PIC::Rnd::CenterNode::Seed(int i,int j,int k,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  int iThreadOpneMP,nTotalThreadsOpenMP;
  char* CenterNodeAssociatedDataBufferPointer;
  PIC::Mesh::cDataCenterNode *cell;
  int LocalCellNumber;
  PIC::Mesh::cDataBlockAMR *block;
  unsigned long int Seed=0;

  LocalCellNumber=_getCenterNodeLocalNumber(i,j,k);
  block=node->block;

  if (block!=NULL) {
    LocalCellNumber=_getCenterNodeLocalNumber(i,j,k);
    cell=block->GetCenterNode(LocalCellNumber);

    if (cell!=NULL) {
      CenterNodeAssociatedDataBufferPointer=cell->GetAssociatedDataBufferPointer();

      if (CenterNodeAssociatedDataBufferPointer!=NULL) {
        //copy the seed
        int LengthBlockId=node->AMRnodeID.Length(),iByteSeed=0,iByteBlockId=0,iByteSeedZero=-1;
        unsigned char *SeedPtr,*BlockIdPtr;

        SeedPtr=(unsigned char*)(&Seed);
        BlockIdPtr=(unsigned char*)(node->AMRnodeID.id);

        for (iByteBlockId=0,iByteSeed=0;iByteBlockId<LengthBlockId;iByteBlockId++,iByteSeed++) {
          if (iByteSeed>=sizeof(unsigned long)) iByteSeed=0;

          SeedPtr[iByteSeed]^=BlockIdPtr[iByteBlockId];
          if (SeedPtr[iByteSeed]==0) iByteSeedZero=iByteSeed;
        }

        //Increment seed
        if (iByteSeedZero==-1) iByteSeedZero=0;

        for (int i=0;(i<sizeof(int))&&(iByteSeedZero<sizeof(unsigned long));iByteSeedZero++,i++) {
          unsigned char t;

          t=((unsigned char*)&LocalCellNumber)[i];
          SeedPtr[iByteSeedZero]&=t;
        }

        #if _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
        nTotalThreadsOpenMP=1;
        #else
        nTotalThreadsOpenMP=omp_get_num_threads();
        #endif

        for (iThreadOpneMP=0;iThreadOpneMP<nTotalThreadsOpenMP;iThreadOpneMP++) {
          *(iThreadOpneMP+(unsigned long int*)(CenterNodeAssociatedDataBufferPointer+Offset))=Seed;
          Seed=Seed>>1;
        }
      }
    }
  }
}

void PIC::Rnd::CenterNode::Seed(cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  if (node==NULL) node=PIC::Mesh::mesh.rootTree;

  //loop throught the tree
  if (node->lastBranchFlag()==_BOTTOM_BRANCH_TREE_) {
    PIC::Mesh::cDataBlockAMR *block=node->block;

    if (block!=NULL) {
      int i,j,k;

      for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
        Seed(i,j,k,node);
      }
    }
  }
  else {
    int i;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *downNode;

    for (i=0;i<(1<<DIM);i++) if ((downNode=node->downNode[i])!=NULL) {
       Seed(downNode);
    }
  }
}

unsigned long int *PIC::Rnd::CenterNode::GetSeed(char* CenterNodeAssociatedDataBufferPointer) {
  unsigned long int* seed;

  #if _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
  seed=(unsigned long int*)(CenterNodeAssociatedDataBufferPointer+Offset);
  #elif _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  int thread=omp_get_thread_num();
  seed=thread+(unsigned long int*)(CenterNodeAssociatedDataBufferPointer+Offset);
  #else
  #error Unknown option
  #endif //_COMPILATION_MODE_

  return seed;
}

unsigned long int *PIC::Rnd::CenterNode::GetSeed(int i,int j,int k,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
  int LocalCellNumber;
  PIC::Mesh::cDataCenterNode *cell;
  char *CenterNodeAssociatedDataBufferPointer;

  LocalCellNumber=_getCenterNodeLocalNumber(i,j,k);
  cell=node->block->GetCenterNode(LocalCellNumber);
  CenterNodeAssociatedDataBufferPointer=cell->GetAssociatedDataBufferPointer();

  return GetSeed(CenterNodeAssociatedDataBufferPointer);
}


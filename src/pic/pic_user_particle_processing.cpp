//$Id$
//the functionality of the model particle processing

/*
 * pic_user_particle_processing.cpp
 *
 *  Created on: Mar 14, 2016
 *      Author: vtenishe
 */

#include "pic.h"

//the default function for the particle processing (the default --> do nothing)
void PIC::UserParticleProcessing::Processing_default(int ptr,int& FirstParticleCell,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node) {
  PIC::ParticleBuffer::SetNext(FirstParticleCell,ptr);
  PIC::ParticleBuffer::SetPrev(-1,ptr);

  if (FirstParticleCell!=-1) PIC::ParticleBuffer::SetPrev(ptr,FirstParticleCell);
  FirstParticleCell=ptr;
}

//call the particle processing function
void PIC::UserParticleProcessing::Processing() {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  int i,j,k;
  int oldFirstCellParticle,newFirstCellParticle,p,pnext;

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
#pragma omp parallel default(none)  shared(PIC::Mesh::mesh) \
  private (node,i,j,k,oldFirstCellParticle,newFirstCellParticle,p,pnext)
   {
#pragma omp single
     {
#endif //_COMPILATION_MODE_
  for (node=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::Mesh::mesh->ThisThread];node!=NULL;node=node->nextNodeThisThread) {
    if (node->block!=NULL) {

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
#pragma omp task default (none) firstprivate (node) private (i,j,k,oldFirstCellParticle,newFirstCellParticle,p,pnext)
    {
#endif
      for (k=0;k<_BLOCK_CELLS_Z_;k++) {
         for (j=0;j<_BLOCK_CELLS_Y_;j++) {
            for (i=0;i<_BLOCK_CELLS_X_;i++) {
              oldFirstCellParticle=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
              newFirstCellParticle=-1;

              if (oldFirstCellParticle!=-1) {
                p=oldFirstCellParticle;

                while (p!=-1) {
                  pnext=PIC::ParticleBuffer::GetNext(p);
                  _PIC_USER_PARTICLE_PROCESSING__FUNCTION_(p,newFirstCellParticle,node);
                  p=pnext;
                }

                node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]=newFirstCellParticle;
              }
           }
         }
      }

      //end of the block processing
#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  }
#endif
    }
  }

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
}}
#endif

}



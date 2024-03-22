//the model of electron impact  reactions

#include "pic.h"


void PIC::ChemicalReactions::ElectronImpactIonizationReactions::ExecuteElectronImpactIonizationModel() {
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  double ReactionRate,ParticleWeight;
  int spec,i,j,k,ProductSpeciesIndex;
  long int ptr,ptrnext;


#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
#pragma omp parallel for schedule(dynamic,1) default (none)  \
    private (k,j,i,node,oldFirstCellParticle,newFirstCellParticle,p,pnext) \
    shared (DomainBlockDecomposition::nLocalBlocks,PIC::DomainBlockDecomposition::BlockTable)
#endif

for (int nLocalNode=0;nLocalNode<DomainBlockDecomposition::nLocalBlocks;nLocalNode++) for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++)  for (i=0;i<_BLOCK_CELLS_X_;i++) {
    double StartTime=MPI_Wtime();


    node=DomainBlockDecomposition::BlockTable[nLocalNode];
    if (node->block!=NULL) {
       auto block=node->block;
       PIC::ParticleBuffer::byte *ParticleData;

       ptr=block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

       if (ptr!=-1) {
         while (ptr!=-1) {
           ParticleData=PIC::ParticleBuffer::GetParticleDataPointer(ptr);
           ptrnext=PIC::ParticleBuffer::GetNext(ParticleData);
	   spec=PIC::ParticleBuffer::GetI(ParticleData);

           //check if the reaction occured 
//	   ParticleWeight=block->GetLocalParticleWeight(spec)*PIC::ParticleBuffer::GetIndividualStatWeightCorrection(ParticleData); 
	   ReactionRate=_PIC_ELECTRON_IMPACT_IONIZATION_RATE_(ParticleData,ProductSpeciesIndex,node); ///ParticleWeight; 

	   if (ReactionRate>0.0) {
             if (block->GetLocalTimeStep(spec)>-log(rnd())/ReactionRate) {
                //electron impact reaction occured
		
	        if (ProductSpeciesIndex<0) {
	          //the ion species is not exist in the simulation -> delete the particle 
                  PIC::ParticleBuffer::DeleteParticle(ptr,node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]); 
                }
	        else {
                  if (block->GetLocalParticleWeight(spec)!=block->GetLocalParticleWeight(spec)) {
                    exit(__LINE__,__FILE__,"Error: the case when stat weight of the parent and daughter species are difference is not implemented yet -- need to implement it, or use the same weight for parent and daougher species");
                  }

                  PIC::ParticleBuffer::SetI(ProductSpeciesIndex,ParticleData);
	        }
             }
	   }

	   ptr=ptrnext;
         }
      }
    }

    if (_PIC_DYNAMIC_LOAD_BALANCING_MODE_ == _PIC_DYNAMIC_LOAD_BALANCING_EXECUTION_TIME_) {
      node->ParallelLoadMeasure+=MPI_Wtime()-StartTime;
    }
  }
}

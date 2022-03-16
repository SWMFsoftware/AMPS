
#include "pic.h"

//functions for splitting/merging particles 

void PIC::ParticleSplitting::Split::Scatter(int particle_num_limit) {
  int inode,i,j,k,nParticles;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  long int p;

  for (inode=0;inode<PIC::DomainBlockDecomposition::nLocalBlocks;inode++) {
    node=PIC::DomainBlockDecomposition::BlockTable[inode];

    if (node->block!=NULL) {
      nParticles=0;

      for (i=0;i<_BLOCK_CELLS_X_;i++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (k=0;k<_BLOCK_CELLS_Z_;k++) {
        p=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

        while (p!=-1) {
          nParticles++;
          p=PIC::ParticleBuffer::GetNext(p);
        }
      } 

      if (nParticles<8*particle_num_limit) {
        //the particle number below the limit -> split particles 
 
        class cParticleId {
        public:
          long int ptr;
          int nd;
        };

        //collect particles in the block
        list <cParticleId> ParticleList;
        cParticleId t;
        double ParticleWeight;

        for (i=0;i<_BLOCK_CELLS_X_;i++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (k=0;k<_BLOCK_CELLS_Z_;k++) {
          p=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
          t.nd=_getCenterNodeLocalNumber(i,j,k);


          while (p!=-1) {
            t.ptr=p; 
            ParticleList.push_back(t);

            ParticleWeight=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(p);
            ParticleWeight/=8.0;
            PIC::ParticleBuffer::SetIndividualStatWeightCorrection(ParticleWeight,p);

            p=PIC::ParticleBuffer::GetNext(p);
          }
        }
        
        //scatter the particles  
        double dx[3],*xmin,x[3];
        int nd,idim;

        xmin=node->xmin;

        dx[0]=(node->xmax[0]-xmin[0])/_BLOCK_CELLS_X_;
        dx[1]=(node->xmax[1]-xmin[1])/_BLOCK_CELLS_Y_;
        dx[2]=(node->xmax[2]-xmin[2])/_BLOCK_CELLS_Z_;

        for (i=0;i<_BLOCK_CELLS_X_;i++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (k=0;k<_BLOCK_CELLS_Z_;k++) {
          nd=_getCenterNodeLocalNumber(i,j,k); 

          for (auto& it : ParticleList) {
            if (it.nd!=nd) {
              p=PIC::ParticleBuffer::GetNewParticle(node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]); 
              PIC::ParticleBuffer::CloneParticle(p,it.ptr);
 
              x[0]=xmin[0]+(i+rnd())*dx[0];
              x[1]=xmin[1]+(j+rnd())*dx[1];
              x[2]=xmin[2]+(k+rnd())*dx[2]; 

              PIC::ParticleBuffer::SetX(x,p);
            }
          }
        }
      }
    }
  }
}






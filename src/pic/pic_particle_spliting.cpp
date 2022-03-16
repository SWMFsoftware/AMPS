
#include "pic.h"

//functions for splitting/merging particles 

void PIC::ParticleSplitting::Split::Scatter(int particle_num_limit) {
  int inode,i,j,k,i0,j0,k0,di,dj,dk,nParticles;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  long int p;

  if ((_BLOCK_CELLS_X_%2!=0)||(_BLOCK_CELLS_Y_%2!=0)||(_BLOCK_CELLS_Z_%2!=0)) {
    exit(__LINE__,__FILE__,"Error: the number of cells in a block has to be even in each direction. Change the input file"); 
  }

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

      if (nParticles<particle_num_limit*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_) {
        //the particle number below the limit -> split particles 
 
        class cParticleId {
        public:
          long int ptr;
          int nd;
        };

        double dx[3],*xmin,x[3];
        int nd,idim;

        xmin=node->xmin;

        dx[0]=(node->xmax[0]-xmin[0])/_BLOCK_CELLS_X_;
        dx[1]=(node->xmax[1]-xmin[1])/_BLOCK_CELLS_Y_;
        dx[2]=(node->xmax[2]-xmin[2])/_BLOCK_CELLS_Z_;

        //collect particles in the block
        list <cParticleId> ParticleList;
        cParticleId t;
        double ParticleWeight;

        for (i0=0;i0+1<_BLOCK_CELLS_X_;i0+=2) for (j0=0;j0+1<_BLOCK_CELLS_Y_;j0+=2) for (k0=0;k0+1<_BLOCK_CELLS_Z_;k0+=2) {
          ParticleList.clear(); 

          for (di=0;di<2;di++) for (dj=0;dj<2;dj++) for (dk=0;dk<2;dk++) { 
            i=i0+di;   
            j=j0+dj;
            k=k0+dk; 

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
          for (di=0;di<2;di++) for (dj=0;dj<2;dj++) for (dk=0;dk<2;dk++) {       
            i=i0+di;
            j=j0+dj;
            k=k0+dk;

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
}






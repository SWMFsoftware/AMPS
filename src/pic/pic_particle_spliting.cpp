
#include "pic.h"

//functions for splitting/merging particles 

void PIC::ParticleSplitting::Split::Scatter(int particle_num_limit_min,int particle_num_limit_max) {
  int inode,i,j,k,i0,j0,k0,di,dj,dk,nParticles;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  long int p;

  if ((_BLOCK_CELLS_X_%2!=0)||(_BLOCK_CELLS_Y_%2!=0)||(_BLOCK_CELLS_Z_%2!=0)) {
    exit(__LINE__,__FILE__,"Error: the number of cells in a block has to be even in each direction. Change the input file"); 
  }

  //get the number of model particle in the cell block 2x2x2 started from position i0,j0,k0
  auto GetParticleNumber = [] (int i,int j,int k,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node)  {
    int ParticleNumber=0;
    long int p;

    p=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

    while (p!=-1) {
      ParticleNumber++;
      p=PIC::ParticleBuffer::GetNext(p);
    }

    return ParticleNumber;
  };
  
  auto GetBlockParticleNumber = [&] (int i0,int j0,int k0,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node)  {
    int di,dj,dk,i,j,k;
    int ParticleNumber=0;
    
    for (di=0;di<2;di++) for (dj=0;dj<2;dj++) for (dk=0;dk<2;dk++) { 
      i=i0+di;   
      j=j0+dj;
      k=k0+dk; 
      
      ParticleNumber+=GetParticleNumber(i,j,k,node);
    }
    
    return ParticleNumber;
  };  
  
  class cParticleListElement {
  public:
    long int ptr;
    double w;
  };
  
  
  class cParticleId {
  public:
    long int ptr;
    int nd;
  };
  
  vector<cParticleListElement>  ParticleList;
  
  //reduec the number of model particles 
  auto ReduceParicleNumber = [&] (int i,int j,int k,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,vector<cParticleListElement>& ParticleList,int nRequestedParticles)  { 
    int nTotalParticles=0;
    long int p;
    double SummWeightCorrection=0.0,WeightCorrectionMax=0.0,w,RemovedWeightCorrectionFactor=0.0;
    cParticleListElement t;
    
    ParticleList.clear();
    p=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
    
    while (p!=-1) {
      w=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(p);
      
      t.ptr=p;
      t.w=w;
      SummWeightCorrection+=w;
      
      if (WeightCorrectionMax<w) WeightCorrectionMax=w;
      
      ParticleList.push_back(t);
      nTotalParticles++;
      
      p=PIC::ParticleBuffer::GetNext(p);      
    }
  
    //remove particles 
    int nRemovedParticles=nTotalParticles-nRequestedParticles;
    
    for (int i=0;i<nRemovedParticles;i++) {
      int ii;
      bool flag=false;

      do {
        ii=(int)(rnd()*nTotalParticles);
      
        p=ParticleList[ii].ptr;
      
        if (ParticleList[ii].w/WeightCorrectionMax<rnd()) { //particles with higher weight has more chanse to remain in the system
          RemovedWeightCorrectionFactor+=ParticleList[ii].w;
        
          PIC::ParticleBuffer::DeleteParticle(p);
          ParticleList[ii]=ParticleList[nTotalParticles-1];
          nTotalParticles--;
          flag=true;
        }
      }
      while (flag==false);
    }
    
    //rebuild the particle list 
    double WeightRatio=SummWeightCorrection/(SummWeightCorrection-RemovedWeightCorrectionFactor);
    
    for (int i=0;i<nTotalParticles;i++) {
      w=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(p);
      w*=WeightRatio;
      PIC::ParticleBuffer::SetIndividualStatWeightCorrection(w,p);
      
      if (i==0) {
        PIC::ParticleBuffer::SetPrev(-1,ParticleList[0].ptr);
      }
      else  {
        PIC::ParticleBuffer::SetPrev(ParticleList[i-1].ptr,ParticleList[i].ptr);
      }
      
      if (i==nTotalParticles-1) {
        PIC::ParticleBuffer::SetNext(-1,ParticleList[i].ptr);
      }
      else {
        PIC::ParticleBuffer::SetNext(ParticleList[i+1].ptr,ParticleList[i].ptr);
      }
    }
    
    //replace the starting point of the list 
    node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]=ParticleList[0].ptr;
  };
  
  
  //collect particles in the cell block 2x2x2
  auto CollectParticleBlock = [&] (int i0,int j0,int k0,list<cParticleId>& ParticleList) {
    cParticleId t;
    long int p;
    double ParticleWeight;
    
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
  };
  
  //scatter particles 
  auto ScatterParticles = [&] (int i,int j, int k,double *xmin,double *dx,list<cParticleId>& ParticleList) {
    int nd;
    double x[3];

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
  };
  
  
  for (inode=0;inode<PIC::DomainBlockDecomposition::nLocalBlocks;inode++) {
    node=PIC::DomainBlockDecomposition::BlockTable[inode];

    if (node->block!=NULL) {
      for (i0=0;i0+1<_BLOCK_CELLS_X_;i0+=2) for (j0=0;j0+1<_BLOCK_CELLS_Y_;j0+=2) for (k0=0;k0+1<_BLOCK_CELLS_Z_;k0+=2) {
        nParticles=GetBlockParticleNumber(i0,j0,k0,node);
        
        if (nParticles<particle_num_limit_min*8) { // 8<- the number of cell in the cell block 2x2x2
          //the number of particles is below the limit: scatter the particles 
          list<cParticleId> ParticleList;
          double dx[3],*xmin,x[3];

          xmin=node->xmin;

          dx[0]=(node->xmax[0]-xmin[0])/_BLOCK_CELLS_X_;
          dx[1]=(node->xmax[1]-xmin[1])/_BLOCK_CELLS_Y_;
          dx[2]=(node->xmax[2]-xmin[2])/_BLOCK_CELLS_Z_;
          
          CollectParticleBlock(i0,j0,k0,ParticleList);
      
          
          for (di=0;di<2;di++) for (dj=0;dj<2;dj++) for (dk=0;dk<2;dk++) { 
            i=i0+di;   
            j=j0+dj;
            k=k0+dk; 
          
            ScatterParticles(i,j,k,xmin,dx,ParticleList);
          }
        }
        
        
        vector<cParticleListElement> ParticleList;
      
        for (di=0;di<2;di++) for (dj=0;dj<2;dj++) for (dk=0;dk<2;dk++) { 
          i=i0+di;   
          j=j0+dj;
          k=k0+dk; 
          
          if (GetParticleNumber(i,j,k,node)>particle_num_limit_max) {
            ReduceParicleNumber(i,j,k,node,ParticleList,particle_num_limit_max);
          }
        }
      }
    }
  }
}






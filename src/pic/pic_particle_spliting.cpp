
#include "pic.h"

//functions for splitting/merging particles 


void PIC::ParticleSplitting::Split::SplitWithVelocityShift(double shift_max,int particle_num_limit_min,int particle_num_limit_max) {
  int i,j,k,inode;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;

  auto GetParticleNumber = [&] (int *ParticleNumberTable,int i,int j,int k,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
    long int p;
    int spec;
    
    for (spec=0;spec<PIC::nTotalSpecies;spec++) ParticleNumberTable[spec]=0;
    
    p=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
    
    while (p!=-1) {
      ParticleNumberTable[PIC::ParticleBuffer::GetI(p)]++;
      p=PIC::ParticleBuffer::GetNext(p);
    }    
  };
    
  //Reduce particle number  
  class cParticleDescriptor {
  public:
    long int p;
    double w;
  };

  auto ReduceParticleNumber = [&] (int spec,int i,int j,int k,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node)  {
    vector<cParticleDescriptor> ParticleList;  
    cParticleDescriptor t;
    double *v,w_max=0.0,SummedWeight=0.0,w;
    int nModelParticles=0;
    long int p;

    p=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

    while (p!=-1) {
      if (spec==PIC::ParticleBuffer::GetI(p)) { 
        t.p=p;
        t.w=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(p);
        ParticleList.push_back(t);

        SummedWeight+=t.w;
        nModelParticles++;

        if (w_max<t.w) w_max=t.w;
      }

      p=PIC::ParticleBuffer::GetNext(p);
    }

    //sort the particle list  
    std::sort(ParticleList.begin(),ParticleList.end(),
        [](const cParticleDescriptor& a,const cParticleDescriptor& b) {return a.w>b.w;});
    
    //delete particles 
    int nDeleteParticles=nModelParticles-particle_num_limit_max,ip; 
    double RemovedParticleWeight=0.0;
    int ii;
    
    for (ii=0;ii<nDeleteParticles;ii++) {
      RemovedParticleWeight+=ParticleList[ii].w;
      PIC::ParticleBuffer::DeleteParticle(ParticleList[ii].p,node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]);
    }
    

    //distribute the removed weight among particles that left in the system 
    for (;ii<nModelParticles;ii++) {
      w=ParticleList[ii].w+RemovedParticleWeight/(SummedWeight-RemovedParticleWeight);    
      p=ParticleList[ii].p;

      PIC::ParticleBuffer::SetIndividualStatWeightCorrection(w,p);
    }
  }; 

  auto IncreseParticleNumber = [&] (int spec,int i,int j,int k,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
    double MeanV[3]={0.0,0.0,0.0},*v,w;
    double MeanV2[3]={0.0,0.0,0.0};
    
    vector<cParticleDescriptor> ParticleList;
    cParticleDescriptor t;
    double w_max=0.0,SummedWeight=0.0,ThermalSpeed=0.0;
    int idim,nModelParticles=0;
    long int p;

    p=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

    while (p!=-1) {
      if (spec==PIC::ParticleBuffer::GetI(p)) {
        t.p=p;
        t.w=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(p);
        ParticleList.push_back(t); 

        SummedWeight+=t.w;
        if (w_max<t.w) w_max=t.w;

        v=PIC::ParticleBuffer::GetV(p); 

        for (int idim=0;idim<3;idim++) {
          MeanV[idim]+=t.w*v[idim];
          MeanV2[idim]+=t.w*v[idim]*v[idim];
        }

        nModelParticles++;
      }

      p=PIC::ParticleBuffer::GetNext(p);
    }

    for (idim=0;idim<3;idim++) ThermalSpeed+=MeanV2[idim]/SummedWeight-MeanV[idim]*MeanV[idim]/(SummedWeight*SummedWeight); 

    ThermalSpeed=sqrt(fabs(ThermalSpeed));

    //add new particles in the system
    int nNewParticles=particle_num_limit_min-nModelParticles;
    int ip;
    long int pnew;
    double l[3],*vnew;
    bool flag;

    for (int ii=0;ii<nNewParticles;ii++) {
      flag=true;
 
      while (flag==true) {
        ip=(int)(rnd()*nModelParticles);
  
        if (ParticleList[ip].w/w_max>rnd()) { //split particles with larger weight  
          pnew=PIC::ParticleBuffer::GetNewParticle(node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]); 
          PIC::ParticleBuffer::CloneParticle(pnew,ParticleList[ip].p);
  
          //set new weight for both particles 
          ParticleList[ip].w/=2.0;
          t=ParticleList[ip];
          t.p=pnew;
          ParticleList.push_back(t);
          nModelParticles++;
          flag=false;
  
          PIC::ParticleBuffer::SetIndividualStatWeightCorrection(ParticleList[ip].w,ParticleList[ip].p);
          PIC::ParticleBuffer::SetIndividualStatWeightCorrection(ParticleList[ip].w,pnew);
  
          //perturb a new particle velocity;
          double dSpeed=shift_max*rnd()*ThermalSpeed; 
  
          Vector3D::Distribution::Uniform(l); 
  
          //get velocity for the new particle 
          vnew=PIC::ParticleBuffer::GetV(pnew); 
  
          for (idim=0;idim<3;idim++) vnew[idim]+=dSpeed*l[idim];
        }
      }
    }
  };

  //loop through the nodes
  int ParticleNumberTable[PIC::nTotalSpecies],spec; 

  for (inode=0;inode<PIC::DomainBlockDecomposition::nLocalBlocks;inode++) {
    node=PIC::DomainBlockDecomposition::BlockTable[inode];

    if (node->block!=NULL) {
      for (i=0;i<_BLOCK_CELLS_X_;i++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (k=0;k<_BLOCK_CELLS_Z_;k++) {
        GetParticleNumber(ParticleNumberTable,i,j,k,node);
        
        for (spec=0;spec<PIC::nTotalSpecies;spec++) {
          if (ParticleNumberTable[spec]>particle_num_limit_max) {
            ReduceParticleNumber(spec,i,j,k,node);
          }
          else if ((0<ParticleNumberTable[spec])&&(ParticleNumberTable[spec]<particle_num_limit_min)) {
            IncreseParticleNumber(spec,i,j,k,node);
          }
        }
      }
    }      
  }
  
}




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






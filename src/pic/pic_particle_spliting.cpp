
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
    double *v,w_max=0.0,w_min=-1.0,SummedWeight=0.0,w;
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
        if ((w_min<0.0)||(w_min>t.w)) w_min=t.w; 
      }

      p=PIC::ParticleBuffer::GetNext(p);
    }

    //sort the particle list  
//    std::sort(ParticleList.begin(),ParticleList.end(),
//        [](const cParticleDescriptor& a,const cParticleDescriptor& b) {return a.w>b.w;});
    
    //delete particles 
    int nDeleteParticles=nModelParticles-particle_num_limit_max,ip; 
    double RemovedParticleWeight=0.0;
    int ii;
    
    for (ii=0;ii<nDeleteParticles;ii++) {
      const bool _search_continue=true;
      const bool _search_completed=false; 

      bool flag=_search_continue;
      int ip, cnt=0;

      while ((++cnt<1000000)&&(flag==_search_continue)) {
        ip=(int)(rnd()*nModelParticles);

        if (w_min/ParticleList[ip].w>rnd()) { // delete particles with smaller weight  
          RemovedParticleWeight+=ParticleList[ip].w;
          PIC::ParticleBuffer::DeleteParticle(ParticleList[ip].p,node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]);

          ParticleList[ip]=ParticleList[nModelParticles-1];
          nModelParticles--;
          flag=_search_completed; 
        }
      } 
    }
    

    //distribute the removed weight among particles that left in the system 
    for (ii=0;ii<nModelParticles;ii++) {
      w=ParticleList[ii].w*(1.0+RemovedParticleWeight/(SummedWeight-RemovedParticleWeight));    
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
  auto GetParticleNumber = [] (int spec,int i,int j,int k,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node)  {
    int ParticleNumber=0;
    long int p;

    p=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

    while (p!=-1) {
      if (PIC::ParticleBuffer::GetI(p)==spec) ParticleNumber++;
      p=PIC::ParticleBuffer::GetNext(p);
    }

    return ParticleNumber;
  };
  
  auto GetBlockParticleNumber = [&] (int spec,int i0,int j0,int k0,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node)  {
    int di,dj,dk,i,j,k;
    int ParticleNumber=0;
    
    for (di=0;di<2;di++) for (dj=0;dj<2;dj++) for (dk=0;dk<2;dk++) { 
      i=i0+di;   
      j=j0+dj;
      k=k0+dk; 
      
      ParticleNumber+=GetParticleNumber(spec,i,j,k,node);
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
  
  class cParticleDescriptor {
  public:
    long int p;
    double w;
  };
  
  vector<cParticleListElement>  ParticleList;
  
  //reduce the number of model particles 
  
auto GetVelocityRange = [&] (int spec,double *vmin,double *vmax,int i,int j,int k,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node)  {
  long int p;
  int res=0,idim;
  bool initflag=false;
  double v[3];

    p=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

    while (p!=-1) {
      if (spec==PIC::ParticleBuffer::GetI(p)) {

      res++;
      PIC::ParticleBuffer::GetV(v,p);

      if (initflag==false) {
        for (idim=0;idim<3;idim++) vmin[idim]=v[idim],vmax[idim]=v[idim];
        
        initflag=true;
      }
      else {
        for (idim=0;idim<3;idim++) {
          if (vmin[idim]>v[idim]) vmin[idim]=v[idim];
          if (vmax[idim]<v[idim]) vmax[idim]=v[idim]; 
        }
      }
    }

    p=PIC::ParticleBuffer::GetNext(p);
  }

  return res;
};


/*
auto GetParticleNumber = [&] (int spec,double *vmin,double *vmax,int i,int j,int k,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node)  {
  long int p;
  int res=0;

    p=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

    while (p!=-1) {
      if (spec==PIC::ParticleBuffer::GetI(p)) {
        res++;
      }

      p=PIC::ParticleBuffer::GetNext(p);
    }

  return res;
};
*/
 
class cVelocitySpaceElement {
public:
  int iv,jv,kv;
  vector<long int> p;
}; 

list<cVelocitySpaceElement> VelocitySpaceTable;

//create the initial velovity space table
auto InitVelocitySpaceTable = [&] (int spec,double *vmin,double *dv,int i, int j,int k,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node)  {
  long int p;
  int idim,iv[3];
  double v[3];
  bool foundflag;

  VelocitySpaceTable.clear();

    p=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

    while (p!=-1) {
      if (spec==PIC::ParticleBuffer::GetI(p)) {
        foundflag=false;

      PIC::ParticleBuffer::GetV(v,p);
      

      for (idim=0;idim<3;idim++) iv[idim]=(int)((v[idim]-vmin[idim])/dv[idim]); 

      for (list<cVelocitySpaceElement>::iterator it=VelocitySpaceTable.begin();it!=VelocitySpaceTable.end();it++) {
        if ((it->iv==iv[0])&&(it->jv==iv[1])&&(it->kv==iv[2])) {
          it->p.push_back(p);
          foundflag=true;
          break;
        }
      }

      if (foundflag==false) {
        cVelocitySpaceElement t;

        t.iv=iv[0],t.jv=iv[1],t.kv=iv[2]; 
         
        VelocitySpaceTable.push_front(t);
        VelocitySpaceTable.begin()->p.push_back(p);
      }
    }

    p=PIC::ParticleBuffer::GetNext(p);
  }
}; 

auto RemoveOneParticle= [&] (int nMaxParticlesReduce,int i, int j, int k,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,bool& more_then_one_particle_left) { 
  int cnt=0,npart;
  long int p;
  double w_remove;

  more_then_one_particle_left=false;

  if (nMaxParticlesReduce==0) return cnt;

  for (list<cVelocitySpaceElement>::iterator it=VelocitySpaceTable.begin();it!=VelocitySpaceTable.end();it++) {
    if ((npart=it->p.size())>1) {
      double w_tot=0.0,w;

      p=it->p[0];
      w_remove=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(p);

      //delete particles
      PIC::ParticleBuffer::DeleteParticle(p,node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]);
      it->p.erase(it->p.begin());

      if (npart-1>1) more_then_one_particle_left=true;

      //distribute the weights
      for (auto& t : it->p) w_tot+=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(t); 
      
      for (auto& t : it->p) {
        w=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(t);

        w=w*(1.0+w_remove/w_tot);
        PIC::ParticleBuffer::SetIndividualStatWeightCorrection(w,t);
      } 

      if (++cnt==nMaxParticlesReduce) break;
    }
  } 

  return cnt;
}; 

auto RebuildVelocitySpaceTable = [&] () {
  list<cVelocitySpaceElement> NewVelocitySpaceTable;

  for (auto& it : VelocitySpaceTable) {
    bool found=false;

    for (auto& itNew : NewVelocitySpaceTable) {
      if ((itNew.iv==it.iv/2)&&(itNew.jv==it.jv/2)&&(itNew.kv==it.kv/2)) {
        itNew.p.insert(itNew.p.begin(),it.p.begin(),it.p.end());
        found=true;
        break;
      }
    }   

    if (found==false) {
      NewVelocitySpaceTable.push_front(it);

      NewVelocitySpaceTable.begin()->iv/=2;
      NewVelocitySpaceTable.begin()->jv/=2; 
      NewVelocitySpaceTable.begin()->kv/=2;

//      NewVelocitySpaceTable.begin()->p.insert(NewVelocitySpaceTable.begin()->p.begin(),it.p.begin(),it.p.end());
    }
  }

  VelocitySpaceTable=NewVelocitySpaceTable;
}; 


auto ReduceParticleNumber1 = [&] (int spec,int nModelParticles,int i,int j,int k,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,int nRequestedParticles)  {
  double vmin[3],vmax[3],dv[3]; 
  bool more_then_one_particle_left;

  const int nTestLevels=6;

  if (nModelParticles>nRequestedParticles) {
    int nDeleteParticles=nModelParticles-nRequestedParticles;  
    int nRemovedParticles=0;

    GetVelocityRange(spec,vmin,vmax,i,j,k,node);

    for (int idim=0;idim<3;idim++) dv[idim]=(vmax[idim]-vmin[idim])/(1<<nTestLevels);

    InitVelocitySpaceTable(spec,vmin,dv,i,j,k,node);

    while (nRemovedParticles<nDeleteParticles) {
      nRemovedParticles+=RemoveOneParticle(nDeleteParticles-nRemovedParticles,i,j,k,node,more_then_one_particle_left); 

      if ((nRemovedParticles<nDeleteParticles)&&(more_then_one_particle_left==false)) {
        RebuildVelocitySpaceTable();
      } 
    }
  }
};

auto ReduceParticleNumber = [&] (int spec,int i,int j,int k,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,int nRequestedParticles)  {
    vector<cParticleDescriptor> ParticleList;  
    cParticleDescriptor t;
    double *v,w_max=0.0,w_min=-1.0,SummedWeight=0.0,w;
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
        if ((w_min<0.0)||(w_min>t.w)) w_min=t.w; 
      }

      p=PIC::ParticleBuffer::GetNext(p);
    }

    //delete particles 
    int nDeleteParticles=nModelParticles-nRequestedParticles,ip; 
    double RemovedParticleWeight=0.0;
    int ii;
    
    for (ii=0;ii<nDeleteParticles;ii++) {
      const bool _search_continue=true;
      const bool _search_completed=false; 

      bool flag=_search_continue;
      int ip, cnt=0;

      while ((++cnt<1000000)&&(flag==_search_continue)) {
        ip=(int)(rnd()*nModelParticles);

        if (w_min/ParticleList[ip].w>rnd()) { // delete particles with smaller weight  
          RemovedParticleWeight+=ParticleList[ip].w;
          PIC::ParticleBuffer::DeleteParticle(ParticleList[ip].p,node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]);

          ParticleList[ip]=ParticleList[nModelParticles-1];
          nModelParticles--;
          flag=_search_completed; 
        }
      } 
    }
    

    //distribute the removed weight among particles that left in the system 
    for (ii=0;ii<nModelParticles;ii++) {
      w=ParticleList[ii].w*(1.0+RemovedParticleWeight/(SummedWeight-RemovedParticleWeight));    
      p=ParticleList[ii].p;

      PIC::ParticleBuffer::SetIndividualStatWeightCorrection(w,p);
    }
  }; 
  
    
  //collect particles in the cell block 2x2x2
  auto CollectParticleBlock = [&] (int spec,int i0,int j0,int k0,list<cParticleId>& ParticleList) {
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
        if (PIC::ParticleBuffer::GetI(p)==spec) {
          t.ptr=p; 
          ParticleList.push_back(t);
  
          ParticleWeight=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(p);
          ParticleWeight/=8.0;
          PIC::ParticleBuffer::SetIndividualStatWeightCorrection(ParticleWeight,p);
        }

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
      for (int s=0;s<PIC::nTotalSpecies;s++) for (i0=0;i0+1<_BLOCK_CELLS_X_;i0+=2) for (j0=0;j0+1<_BLOCK_CELLS_Y_;j0+=2) for (k0=0;k0+1<_BLOCK_CELLS_Z_;k0+=2) {
        nParticles=GetBlockParticleNumber(s,i0,j0,k0,node);
        
        if ((nParticles!=0)&&(nParticles<particle_num_limit_min*8)) { // 8<- the number of cell in the cell block 2x2x2
          //the number of particles is below the limit: scatter the particles 
          list<cParticleId> ParticleList;
          double dx[3],*xmin,x[3];

          xmin=node->xmin;

          dx[0]=(node->xmax[0]-xmin[0])/_BLOCK_CELLS_X_;
          dx[1]=(node->xmax[1]-xmin[1])/_BLOCK_CELLS_Y_;
          dx[2]=(node->xmax[2]-xmin[2])/_BLOCK_CELLS_Z_;
          
          CollectParticleBlock(s,i0,j0,k0,ParticleList);
      
          
          for (di=0;di<2;di++) for (dj=0;dj<2;dj++) for (dk=0;dk<2;dk++) { 
            i=i0+di;   
            j=j0+dj;
            k=k0+dk; 
          
            ScatterParticles(i,j,k,xmin,dx,ParticleList);
          }
        }
              
        for (di=0;di<2;di++) for (dj=0;dj<2;dj++) for (dk=0;dk<2;dk++) { 
          i=i0+di;   
          j=j0+dj;
          k=k0+dk; 

          int npart;
          
          if ((npart=GetParticleNumber(s,i,j,k,node))>particle_num_limit_max) {
            ReduceParticleNumber1(s,npart,i,j,k,node,particle_num_limit_min+0.85*(particle_num_limit_max-particle_num_limit_min));
          }
        }
      }
    }
  }
}  






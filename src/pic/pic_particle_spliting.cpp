


#include <list>
#include <algorithm> 

#include "pic.h"


int PIC::ParticleSplitting::particle_num_limit_min=100;
int PIC::ParticleSplitting::particle_num_limit_max=150;
double PIC::ParticleSplitting::v_shift_max=0.01;
double PIC::ParticleSplitting::x_shift_max=0.2;
bool PIC::ParticleSplitting::apply_non_uniform_x_shift=false;
int PIC::ParticleSplitting::Mode=PIC::ParticleSplitting::_disactivated;

//functions for splitting/merging particles 
void PIC::ParticleSplitting::Split::SplitWithVelocityShift(int particle_num_limit_min,int particle_num_limit_max) {
  int inode,i,j,k,i0,j0,k0,di,dj,dk,nParticles;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  long int p;

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


  class cVelocitySpaceElement {
  public:
    int iv,jv,kv;
    vector<long int> p;
  }; 

  cVelocitySpaceElement **VelocitySpaceTable=NULL;
  const int nSerachLevels=6;

  int nVelCells1d=1<<nSerachLevels;
  int nVelCells=nVelCells1d*nVelCells1d*nVelCells1d;


  //create the initial velovity space table
  auto InitVelocitySpaceTable = [&] (int spec,double *vmin,double *dv,int i, int j,int k,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node)  {
    long int p;
    int idim,iv[3],iVelCell;
    double v[3];
    bool foundflag;

    VelocitySpaceTable=new cVelocitySpaceElement* [nVelCells];
    for (int i=0;i<nVelCells;i++) VelocitySpaceTable[i]=NULL;

    p=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

    while (p!=-1) {
      if (spec==PIC::ParticleBuffer::GetI(p)) {
        PIC::ParticleBuffer::GetV(v,p);

        for (idim=0;idim<3;idim++) iv[idim]=(int)((v[idim]-vmin[idim])/dv[idim]); 

        iVelCell=iv[0]+nVelCells1d*(iv[1]+nVelCells1d*iv[2]);  

        if (VelocitySpaceTable[iVelCell]==NULL) {
          VelocitySpaceTable[iVelCell]=new cVelocitySpaceElement [1];

          VelocitySpaceTable[iVelCell]->iv=iv[0];
          VelocitySpaceTable[iVelCell]->jv=iv[1];
          VelocitySpaceTable[iVelCell]->kv=iv[2];
        } 

        VelocitySpaceTable[iVelCell]->p.push_back(p);
      }

      p=PIC::ParticleBuffer::GetNext(p);
    }
  }; 

  auto RemoveOneParticle= [&] (int nCurrentResolutionLevel,int nMaxParticlesReduce,int i, int j, int k,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,bool& more_then_one_particle_left) { 
    int cnt=0,npart;
    long int p;
    double w_remove;

    more_then_one_particle_left=false;

    if (nMaxParticlesReduce==0) return cnt;


    int nCells1d=1<<nCurrentResolutionLevel;
    int nCells=nCells1d*nCells1d*nCells1d;


    for (int ii=0;ii<nCells;ii++) if (VelocitySpaceTable[ii]!=NULL) if ((npart=VelocitySpaceTable[ii]->p.size())>1) { 
      double w_tot=0.0,w,wmin;
      auto itmin=VelocitySpaceTable[ii]->p.begin(); 


      //find particle with smallest weight
      p=-1,wmin=-1.0;
    
  
      for (auto it=VelocitySpaceTable[ii]->p.begin();it!=VelocitySpaceTable[ii]->p.end();it++) {
        w=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(*it);

        if ((wmin<0.0)||(w<wmin)) p=*it,wmin=w,itmin=it;
      }


/*
      for (auto& t : VelocitySpaceTable[ii]->p) {
        w=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(t); 

        if ((wmin<0.0)||(w<wmin)) p=t,wmin=w;
      }
*/

      //find particles that closest to 'p' in the velocity space 
      long int p1;
      double d2,d2min=-1.0;
      double v1[3],v2[3];

      PIC::ParticleBuffer::GetV(v1,p);

      for (auto& t : VelocitySpaceTable[ii]->p) if (t!=p) {
        PIC::ParticleBuffer::GetV(v2,t);
        d2=0.0;
        
        for (int idim=0;idim<3;idim++) {
          double tt=v1[idim]-v2[idim];

          d2+=tt*tt;
        }

        if ((d2min<0.0)||(d2<d2min)) d2min=d2,p1=t;
      } 
 
      //merge particle p1 and p: weight of 'p' is added to 'p1'
      w=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(p1); 

      PIC::ParticleBuffer::GetV(v1,p1);
      PIC::ParticleBuffer::GetV(v2,p);

      for (int idim=0;idim<3;idim++) v1[idim]=(w*v1[idim]+wmin*v2[idim])/(w+wmin);

      PIC::ParticleBuffer::SetV(v1,p1);  
      PIC::ParticleBuffer::SetIndividualStatWeightCorrection(w+wmin,p1);
      
      //delete particles
      PIC::ParticleBuffer::DeleteParticle(p,node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]);
      VelocitySpaceTable[ii]->p.erase(itmin);

      if (npart-1>1) more_then_one_particle_left=true;


/*
      p=VelocitySpaceTable[ii]->p[0];
      w_remove=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(p);

      //delete particles
      PIC::ParticleBuffer::DeleteParticle(p,node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]);
      VelocitySpaceTable[ii]->p.erase(VelocitySpaceTable[ii]->p.begin());

      if (npart-1>1) more_then_one_particle_left=true;

      //distribute the weights
      for (auto& t : VelocitySpaceTable[ii]->p) w_tot+=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(t); 

      for (auto& t : VelocitySpaceTable[ii]->p) {
        w=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(t);

        w=w*(1.0+w_remove/w_tot);
        PIC::ParticleBuffer::SetIndividualStatWeightCorrection(w,t);
      } 
*/

      if (++cnt==nMaxParticlesReduce) break;
    }

    return cnt;
  }; 

  auto RebuildVelocitySpaceTable = [&] (int CurrentResolutionLevel) {
    cVelocitySpaceElement **NewVelocitySpaceTable;
    int iVelCell,iNewVelCell,iv,jv,kv;


    int NewResolutionLevel=CurrentResolutionLevel-1;

    int nVelCells1d=1<<CurrentResolutionLevel;
    int nNewVelCells1d=nVelCells1d/2;

    int nVelCells=nVelCells1d*nVelCells1d*nVelCells1d;
    int nNewVelCells=nNewVelCells1d*nNewVelCells1d*nNewVelCells1d;


    NewVelocitySpaceTable=new cVelocitySpaceElement* [nNewVelCells];
    for (int i=0;i<nNewVelCells;i++) NewVelocitySpaceTable[i]=NULL;

    for (int i=0;i<nVelCells;i++) {
      if (VelocitySpaceTable[i]!=NULL) {
        iv=VelocitySpaceTable[i]->iv;
        jv=VelocitySpaceTable[i]->jv;
        kv=VelocitySpaceTable[i]->kv;

        iVelCell=iv+nVelCells1d*(jv+nVelCells1d*kv);  

        iv/=2,jv/=2,kv/=2;
        iNewVelCell=iv+nNewVelCells1d*(jv+nNewVelCells1d*kv);

        if (NewVelocitySpaceTable[iNewVelCell]==NULL) {
          NewVelocitySpaceTable[iNewVelCell]=new cVelocitySpaceElement [1];

          *NewVelocitySpaceTable[iNewVelCell]=*VelocitySpaceTable[i];

          NewVelocitySpaceTable[iNewVelCell]->iv/=2;
          NewVelocitySpaceTable[iNewVelCell]->jv/=2;
          NewVelocitySpaceTable[iNewVelCell]->kv/=2;

        }
        else {
          NewVelocitySpaceTable[iNewVelCell]->p.insert(NewVelocitySpaceTable[iNewVelCell]->p.begin(),
              VelocitySpaceTable[i]->p.begin(),VelocitySpaceTable[i]->p.end());
        }
      }
    }

    //delete the previous table 
    for (int i=0;i<nVelCells;i++) if (VelocitySpaceTable[i]!=NULL) delete [] VelocitySpaceTable[i];

    delete [] VelocitySpaceTable;
    VelocitySpaceTable=NewVelocitySpaceTable;
  }; 


  auto DeleteVelocitySpaceTable = [&] (int CurrentResolutionLevel) {
    int nVelCells1d=1<<CurrentResolutionLevel;
    int nVelCells=nVelCells1d*nVelCells1d*nVelCells1d;

    for (int i=0;i<nVelCells;i++) if (VelocitySpaceTable[i]!=NULL) delete [] VelocitySpaceTable[i];

    delete [] VelocitySpaceTable;
  };

  auto ReduceParticleNumber1 = [&] (int spec,int nModelParticles,int i,int j,int k,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,int nRequestedParticles)  {
    double vmin[3],vmax[3],dv[3]; 
    bool more_then_one_particle_left;


if (false) {
    class cParticleListElement {
    public:
      long int p;
      double w;
    }; 

    list<cParticleListElement> p_list;
    long int p;

    p=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

    while (p!=-1) {
      if (spec==PIC::ParticleBuffer::GetI(p)) {
        cParticleListElement t;

        t.p=p;
        t.w=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(p);

        p_list.push_back(t);
      }

      p=PIC::ParticleBuffer::GetNext(p);
    }

    while ((nModelParticles>nRequestedParticles)&&(p_list.size()!=0)) {
      //1. find a particle with a smallest weight
      list<cParticleListElement>::iterator p_wmin_it=p_list.begin();
      list<cParticleListElement>::iterator it;
      double wmin=p_list.begin()->w;

      for (it=p_list.begin();it!=p_list.end();it++) { 
        if (wmin>it->w) {
          wmin=it->w;
          p_wmin_it=it;
        } 
      }

      //2. find a particle that is closest to p_wmin_it in the velocity space 
      double r2,r2min=-1.0;
      list<cParticleListElement>::iterator p_r2min_it;
      double d,v0[3],v1[3];
      int idim;

      PIC::ParticleBuffer::GetV(v0,p_wmin_it->p);   

      for (it=p_list.begin();it!=p_list.end();it++) if (it!=p_wmin_it) {
        PIC::ParticleBuffer::GetV(v1,it->p);
        r2=0.0;

        for (idim=0;idim<3;idim++) {
          d=v0[idim]-v1[idim]; 
          r2+=d*d;
        }

        if ((r2min<0.0)||(r2min>r2)) {
          r2min=r2,p_r2min_it=it;
        }
      }

      //merge particles if the velocity difference to not too large 
      if ((r2min>0.0)&&(r2min<0.0000001*Vector3D::DotProduct(v0,v0))) {
        double w0,w1;

        PIC::ParticleBuffer::GetV(v1,p_r2min_it->p); 

        w0=p_wmin_it->w;
        w1=p_r2min_it->w;

        for (idim=0;idim<3;idim++) v1[idim]=(w0*v0[idim]+w1*v1[idim])/(w0+w1); 

        PIC::ParticleBuffer::SetV(v1,p_r2min_it->p);
        PIC::ParticleBuffer::SetIndividualStatWeightCorrection(w0+w1,p_r2min_it->p); 

        p_r2min_it->w=w0+w1;

        PIC::ParticleBuffer::DeleteParticle(p_wmin_it->p,node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]);
        nModelParticles--;
      }

      p_list.erase(p_wmin_it); 
    }

   return;
}

    if (nModelParticles>nRequestedParticles) {
      int nDeleteParticles=nModelParticles-nRequestedParticles;  
      int nRemovedParticles=0;


static int ncall=0;
ncall++;



if (ncall==10792) {
  cout << "sfpsodfposiapfoipo"<< endl;
}




      GetVelocityRange(spec,vmin,vmax,i,j,k,node);

      for (int idim=0;idim<3;idim++) {
        double d=vmax[idim]-vmin[idim]; 

        vmin[idim]-=0.05*d; 
        vmax[idim]+=0.05*d;

        dv[idim]=(vmax[idim]-vmin[idim])/nVelCells1d;
      }


      InitVelocitySpaceTable(spec,vmin,dv,i,j,k,node);

      int CurrentResolutionLevel=nSerachLevels;

      while (nRemovedParticles<nDeleteParticles) {
        nRemovedParticles+=RemoveOneParticle(CurrentResolutionLevel,nDeleteParticles-nRemovedParticles,i,j,k,node,more_then_one_particle_left); 

        if ((nRemovedParticles<nDeleteParticles)&&(more_then_one_particle_left==false)) {
          RebuildVelocitySpaceTable(CurrentResolutionLevel);
          --CurrentResolutionLevel;
        } 
      }

      DeleteVelocitySpaceTable(CurrentResolutionLevel);
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

  auto IncreseParticleNumber = [&] (int spec,int i,int j,int k,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
    double MeanV[3]={0.0,0.0,0.0},v[3],w;
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

        PIC::ParticleBuffer::GetV(v,p); 

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
    int nNewParticles=particle_num_limit_max-nModelParticles;
    int ip;
    long int pnew;
    double l[3],vnew[3],vold[3];
    bool flag;

    for (int ii=0;ii<nNewParticles;ii++) {
      flag=true;

      while (flag==true) {
        ip=(int)(rnd()*nModelParticles);

        int ip_wmax=-1;
        w_max=-1.0;
 
        for (int i=0;i<nModelParticles;i++) {
          if (w_max<ParticleList[i].w) w_max=ParticleList[i].w,ip_wmax=i;
        } 

        ip=ip_wmax;

        if (true) { // (ParticleList[ip].w/w_max>rnd()) { //split particles with larger weight  
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
          Vector3D::Distribution::Uniform(l,0.5*v_shift_max*ThermalSpeed); 

          //get velocity for the new particle 
          PIC::ParticleBuffer::GetV(vnew,pnew); 

          for (idim=0;idim<3;idim++) {
            vold[idim]=vnew[idim]-l[idim];
            vnew[idim]+=l[idim];
          }

          PIC::ParticleBuffer::SetV(vnew,pnew);
          PIC::ParticleBuffer::SetV(vold,ParticleList[ip].p);

          double dx[3],*xmin,*xmax,x[3],x0[3];

          xmin=node->xmin;
          

          dx[0]=(node->xmax[0]-xmin[0])/_BLOCK_CELLS_X_;
          dx[1]=(node->xmax[1]-xmin[1])/_BLOCK_CELLS_Y_;
          dx[2]=(node->xmax[2]-xmin[2])/_BLOCK_CELLS_Z_;




 



double r;
bool flag=true;
do {
  flag=true;


    r=pow(rnd(),1.0/3.0);

  
  if (apply_non_uniform_x_shift==true) if (r>rnd()) {
    flag=false;
    continue;
  }

  PIC::ParticleBuffer::GetX(x0,pnew);

  double ll[3];

  Vector3D::Distribution::Uniform(ll,x_shift_max*r*sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]));

  x0[0]+=ll[0];
  if (x0[0]<xmin[0]+i*dx[0]) {flag=false;continue;}
  if (x0[0]>xmin[0]+(i+1)*dx[0]) {flag=false;continue;} 

  x0[1]+=ll[1];
  if (x0[1]<xmin[1]+j*dx[1]) {flag=false;continue;}
  if (x0[1]>xmin[1]+(j+1)*dx[1]) {flag=false;continue;}

  x0[2]+=ll[2];
  if (x0[2]<xmin[2]+k*dx[2]) {flag=false;continue;}
  if (x0[2]>xmin[2]+(k+1)*dx[2]) {flag=false;continue;}
}
while (flag==false);



PIC::ParticleBuffer::SetX(x0,pnew);



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
            ReduceParticleNumber1(spec,ParticleNumberTable[spec],i,j,k,node,(int)(particle_num_limit_min+0.9*(particle_num_limit_max-particle_num_limit_min)));
          }
          else if ((1<ParticleNumberTable[spec])&&(ParticleNumberTable[spec]<particle_num_limit_min)) {
            IncreseParticleNumber(spec,i,j,k,node);
          }
        }
      }
    }      
  }

//  PIC::Parallel::ExchangeParticleData();
}  


//=======================================================================================
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


  class cVelocitySpaceElement {
  public:
    int iv,jv,kv;
    vector<long int> p;
  }; 

  cVelocitySpaceElement **VelocitySpaceTable=NULL;
  const int nSerachLevels=6;

  int nVelCells1d=1<<nSerachLevels;
  int nVelCells=nVelCells1d*nVelCells1d*nVelCells1d;


  //create the initial velovity space table
  auto InitVelocitySpaceTable = [&] (int spec,double *vmin,double *dv,int i, int j,int k,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node)  {
    long int p;
    int idim,iv[3],iVelCell;
    double v[3];
    bool foundflag;

    VelocitySpaceTable=new cVelocitySpaceElement* [nVelCells];
    for (int i=0;i<nVelCells;i++) VelocitySpaceTable[i]=NULL;

    p=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

    while (p!=-1) {
      if (spec==PIC::ParticleBuffer::GetI(p)) {
        PIC::ParticleBuffer::GetV(v,p);

        for (idim=0;idim<3;idim++) iv[idim]=(int)((v[idim]-vmin[idim])/dv[idim]); 

        iVelCell=iv[0]+nVelCells1d*(iv[1]+nVelCells1d*iv[2]);  

        if (VelocitySpaceTable[iVelCell]==NULL) {
          VelocitySpaceTable[iVelCell]=new cVelocitySpaceElement [1];

          VelocitySpaceTable[iVelCell]->iv=iv[0];
          VelocitySpaceTable[iVelCell]->jv=iv[1];
          VelocitySpaceTable[iVelCell]->kv=iv[2];
        } 

        VelocitySpaceTable[iVelCell]->p.push_back(p);
      }

      p=PIC::ParticleBuffer::GetNext(p);
    }
  }; 

  auto RemoveOneParticle= [&] (int nCurrentResolutionLevel,int nMaxParticlesReduce,int i, int j, int k,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,bool& more_then_one_particle_left) { 
    int cnt=0,npart;
    long int p;
    double w_remove;

    more_then_one_particle_left=false;

    if (nMaxParticlesReduce==0) return cnt;


    int nCells1d=1<<nCurrentResolutionLevel;
    int nCells=nCells1d*nCells1d*nCells1d;



    //create a list of the non-zero cells in the velocity space
    class cVelListEl {
    public:
      cVelocitySpaceElement* SpaceTableEl;
      int nTotalParticles;
    }; 

    vector <cVelListEl> SpaceTableElList;

    for (int ii=0;ii<nCells;ii++) if (VelocitySpaceTable[ii]!=NULL) if ((npart=VelocitySpaceTable[ii]->p.size())>1) {
      cVelListEl t;
  
      t.SpaceTableEl=VelocitySpaceTable[ii];
      t.nTotalParticles=npart;

      SpaceTableElList.push_back(t);
    } 


    auto cmp = [] (const cVelListEl& a, const cVelListEl &b) {
      return a.nTotalParticles>b.nTotalParticles;
    }; 
    
    sort(SpaceTableElList.begin(),SpaceTableElList.end(),cmp ); 
//      sort(SpaceTableElList.begin(),SpaceTableElList.end(), [] (const cVelListEl& a, const cVelListEl &b) {return (a.nTotalParticles)>(b.nTotalParticles);});  


    auto RemoveParticle = [&] (vector <cVelListEl>::iterator& it) {
      if (it->nTotalParticles==1) return;
      it->nTotalParticles--;   

      cVelocitySpaceElement* VelocitySpaceTableEl=it->SpaceTableEl;  

      double w_tot=0.0,w;

      p=VelocitySpaceTableEl->p[0];
      w_remove=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(p);

      //delete particles
      PIC::ParticleBuffer::DeleteParticle(p,node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)]);
      VelocitySpaceTableEl->p.erase(VelocitySpaceTableEl->p.begin());
      
      //distribute the weights
      for (auto& t : VelocitySpaceTableEl->p) w_tot+=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(t);
      
      for (auto& t : VelocitySpaceTableEl->p) {
        w=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(t);
      
        w=w*(1.0+w_remove/w_tot);
        PIC::ParticleBuffer::SetIndividualStatWeightCorrection(w,t);
      }
    }; 

    auto it_last=SpaceTableElList.begin();

    if (SpaceTableElList.size()==0) return 0;


    for (it_last=SpaceTableElList.begin();it_last!=SpaceTableElList.end();it_last++) {
      if (SpaceTableElList.begin()->nTotalParticles!=it_last->nTotalParticles) {
        break;
      }
    }

while (cnt<nMaxParticlesReduce) {
    for (auto it=SpaceTableElList.begin();it!=it_last;it++) {
      if (it->nTotalParticles!=1) {
        RemoveParticle(it);
        if (++cnt==nMaxParticlesReduce) return cnt;
      }
    }

    if (it_last!=SpaceTableElList.end()) if (it_last->nTotalParticles==SpaceTableElList.begin()->nTotalParticles) {
      for (;it_last!=SpaceTableElList.end();it_last++) {
        if (it_last->nTotalParticles!=SpaceTableElList.begin()->nTotalParticles) {
          break;
        }
      }

      auto t=it_last;

      if (++t==SpaceTableElList.end()) it_last=SpaceTableElList.end();
    }  

    if (SpaceTableElList.begin()->nTotalParticles==1) {
      break;
    }
}
     

    return cnt;
  }; 

  auto RebuildVelocitySpaceTable = [&] (int CurrentResolutionLevel) {
    cVelocitySpaceElement **NewVelocitySpaceTable;
    int iVelCell,iNewVelCell,iv,jv,kv;


    int NewResolutionLevel=CurrentResolutionLevel-1;

    int nVelCells1d=1<<CurrentResolutionLevel;
    int nNewVelCells1d=nVelCells1d/2;

    int nVelCells=nVelCells1d*nVelCells1d*nVelCells1d;
    int nNewVelCells=nNewVelCells1d*nNewVelCells1d*nNewVelCells1d;


    NewVelocitySpaceTable=new cVelocitySpaceElement* [nNewVelCells];
    for (int i=0;i<nNewVelCells;i++) NewVelocitySpaceTable[i]=NULL;

    for (int i=0;i<nVelCells;i++) {
      if (VelocitySpaceTable[i]!=NULL) {
        iv=VelocitySpaceTable[i]->iv;
        jv=VelocitySpaceTable[i]->jv;
        kv=VelocitySpaceTable[i]->kv;

        iVelCell=iv+nVelCells1d*(jv+nVelCells1d*kv);  

        iv/=2,jv/=2,kv/=2;
        iNewVelCell=iv+nNewVelCells1d*(jv+nNewVelCells1d*kv);

        if (NewVelocitySpaceTable[iNewVelCell]==NULL) {
          NewVelocitySpaceTable[iNewVelCell]=new cVelocitySpaceElement [1];

          *NewVelocitySpaceTable[iNewVelCell]=*VelocitySpaceTable[i];

          NewVelocitySpaceTable[iNewVelCell]->iv/=2;
          NewVelocitySpaceTable[iNewVelCell]->jv/=2;
          NewVelocitySpaceTable[iNewVelCell]->kv/=2;

        }
        else {
          NewVelocitySpaceTable[iNewVelCell]->p.insert(NewVelocitySpaceTable[iNewVelCell]->p.begin(),
              VelocitySpaceTable[i]->p.begin(),VelocitySpaceTable[i]->p.end());
        }
      }
    }

    //delete the previous table 
    for (int i=0;i<nVelCells;i++) if (VelocitySpaceTable[i]!=NULL) delete [] VelocitySpaceTable[i];

    delete [] VelocitySpaceTable;
    VelocitySpaceTable=NewVelocitySpaceTable;
  }; 

  auto DeleteVelocitySpaceTable = [&] (int CurrentResolutionLevel) {
    int nVelCells1d=1<<CurrentResolutionLevel;
    int nVelCells=nVelCells1d*nVelCells1d*nVelCells1d;

    for (int i=0;i<nVelCells;i++) if (VelocitySpaceTable[i]!=NULL) delete [] VelocitySpaceTable[i];

    delete [] VelocitySpaceTable;
  };

  auto ReduceParticleNumber1 = [&] (int spec,int nModelParticles,int i,int j,int k,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node,int nRequestedParticles)  {
    double vmin[3],vmax[3],dv[3]; 
    bool more_then_one_particle_left;

    if (nModelParticles>nRequestedParticles) {
      int nDeleteParticles=nModelParticles-nRequestedParticles;  
      int nRemovedParticles=0;

      GetVelocityRange(spec,vmin,vmax,i,j,k,node);

      for (int idim=0;idim<3;idim++) {
        double d=vmax[idim]-vmin[idim]; 

        vmin[idim]-=0.05*d; 
        vmax[idim]+=0.05*d;

        dv[idim]=(vmax[idim]-vmin[idim])/nVelCells1d;
      }


      InitVelocitySpaceTable(spec,vmin,dv,i,j,k,node);

      int CurrentResolutionLevel=nSerachLevels;

      while (nRemovedParticles<nDeleteParticles) {
        nRemovedParticles+=RemoveOneParticle(CurrentResolutionLevel,nDeleteParticles-nRemovedParticles,i,j,k,node,more_then_one_particle_left); 

        if ((nRemovedParticles<nDeleteParticles)&&(more_then_one_particle_left==false)) {
          RebuildVelocitySpaceTable(CurrentResolutionLevel);
          --CurrentResolutionLevel;
        } 
      }

      DeleteVelocitySpaceTable(CurrentResolutionLevel);    
    }
  };

  auto PrintStat = [&] (int spec,int i,int j,int k,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node) {
    double w_tot=0.0,v[3],v_tot[3]={0.0,0.0,0.0},w;
    int nModelParticles=0;
    long int p;

    p=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

    while (p!=-1) {
      if (spec==PIC::ParticleBuffer::GetI(p)) {
        PIC::ParticleBuffer::GetV(v,p);
        w=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(p);

        for (int i=0;i<3;i++) v_tot[i]+=w*v[i];


        w_tot+=w;
        nModelParticles++;
      }

      p=PIC::ParticleBuffer::GetNext(p);
    }

    if (nModelParticles!=0) for (int i=0;i<3;i++) v_tot[i]/nModelParticles;

//    cout << "w_tot=" << w_tot << ", npart=" << nModelParticles << ", v=" << v_tot[0] << ", " << v_tot[1] << ", " << v_tot[2] << ", i,j,k=" << i << ", " << j << ", " << k << endl;
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
          int CurrentResolutionLevel=nSerachLevels;

          if ((npart=GetParticleNumber(s,i,j,k,node))>particle_num_limit_max) {

if (PIC::ThisThread==0) PrintStat(s,i,j,k,node);

            ReduceParticleNumber1(s,npart,i,j,k,node,particle_num_limit_min+0.9*(particle_num_limit_max-particle_num_limit_min));

if (PIC::ThisThread==0) PrintStat(s,i,j,k,node);
          }
        }
      }
    }
  }
}  

///////////////////////////////////////////////////////////////////////////////////////////////////////
void PIC::ParticleSplitting::Split::SplitWithVelocityShift_FL(int particle_num_limit_min,int particle_num_limit_max,double WeightSplittingLimit) {
  int inode,i,j,k,i0,j0,k0,di,dj,dk,nParticles;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  long int p;
  
  namespace FL = PIC::FieldLine;
  namespace PB = PIC::ParticleBuffer;

  auto GetParticleNumber = [&] (int *ParticleNumberTable,FL::cFieldLineSegment *Segment) {
    long int p;
    int spec;

    for (spec=0;spec<PIC::nTotalSpecies;spec++) ParticleNumberTable[spec]=0;

    p=Segment->FirstParticleIndex;

    while (p!=-1) {
      ParticleNumberTable[PIC::ParticleBuffer::GetI(p)]++;
      p=PIC::ParticleBuffer::GetNext(p);
    }
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
  auto GetVelocityRange = [&] (int spec,double *vmin,double *vmax,FL::cFieldLineSegment *Segment)  {
    long int p;
    int res=0,idim;
    bool initflag=false;
    double vNorm,vParallel,v[3];

    p=Segment->FirstParticleIndex;

    while (p!=-1) {
      if (spec==PIC::ParticleBuffer::GetI(p)) {

        res++;
        v[0]=PB::GetVParallel(p);
        v[1]=PB::GetVNormal(p);
        v[2]=0.0;

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


  class cVelocitySpaceElement {
  public:
    int iv,jv,kv;
    vector<long int> p;
  }; 

  cVelocitySpaceElement **VelocitySpaceTable=NULL;
  const int nSerachLevels=6;

  int nVelCells1d=1<<nSerachLevels;
  int nVelCells=nVelCells1d*nVelCells1d;


  //create the initial velovity space table
  auto InitVelocitySpaceTable = [&] (int spec,double *vmin,double *dv,FL::cFieldLineSegment *Segment)  {
    long int p;
    int idim,iv[3],iVelCell;
    double v[3];
    bool foundflag;

    VelocitySpaceTable=new cVelocitySpaceElement* [nVelCells];
    for (int i=0;i<nVelCells;i++) VelocitySpaceTable[i]=NULL;

    p=Segment->FirstParticleIndex;

    while (p!=-1) {
      if (spec==PIC::ParticleBuffer::GetI(p)) {
        v[0]=PB::GetVParallel(p);
        v[1]=PB::GetVNormal(p);
        v[2]=0.0;

        for (idim=0;idim<2;idim++) {
          iv[idim]=(int)((v[idim]-vmin[idim])/dv[idim]); 

          if (iv[idim]<0) iv[idim]=0;
          if (iv[idim]>=nVelCells1d) iv[idim]=nVelCells1d-1;   
        } 

        iVelCell=iv[0]+nVelCells1d*iv[1];  

        if (VelocitySpaceTable[iVelCell]==NULL) {
          VelocitySpaceTable[iVelCell]=new cVelocitySpaceElement [1];

          VelocitySpaceTable[iVelCell]->iv=iv[0];
          VelocitySpaceTable[iVelCell]->jv=iv[1];
          VelocitySpaceTable[iVelCell]->kv=iv[2];
        } 

        VelocitySpaceTable[iVelCell]->p.push_back(p);
      }

      p=PIC::ParticleBuffer::GetNext(p);
    }
  }; 

  auto RemoveOneParticle= [&] (int nCurrentResolutionLevel,int nMaxParticlesReduce,FL::cFieldLineSegment *Segment,bool& more_then_one_particle_left) { 
    int cnt=0,npart;
    long int p;
    double w_remove;

    more_then_one_particle_left=false;

    if (nMaxParticlesReduce==0) return cnt;


    int nCells1d=1<<nCurrentResolutionLevel;
    int nCells=nCells1d*nCells1d;


    for (int ii=0;ii<nCells;ii++) if (VelocitySpaceTable[ii]!=NULL) if ((npart=VelocitySpaceTable[ii]->p.size())>1) { 
      double w_tot=0.0,w,wmin;
      auto itmin=VelocitySpaceTable[ii]->p.begin(); 


      //find particle with smallest weight
      p=-1,wmin=-1.0;
    
  
      for (auto it=VelocitySpaceTable[ii]->p.begin();it!=VelocitySpaceTable[ii]->p.end();it++) {
        w=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(*it);

        if ((wmin<0.0)||(w<wmin)) p=*it,wmin=w,itmin=it;
      }


      //find particles that closest to 'p' in the velocity space 
      long int p1;
      double d2,d2min=-1.0;
      double v1[3],v2[3];

    //  PIC::ParticleBuffer::GetV(v1,p);
      
      v1[0]=PB::GetVParallel(p);
      v1[1]=PB::GetVNormal(p);
      v1[2]=0.0;

      for (auto& t : VelocitySpaceTable[ii]->p) if (t!=p) {
        //PIC::ParticleBuffer::GetV(v2,p);
        v2[0]=PB::GetVParallel(t);
        v2[1]=PB::GetVNormal(t);
        v2[2]=0.0;


        d2=0.0;
        
        for (int idim=0;idim<3;idim++) {
          double tt=v1[idim]-v2[idim];

          d2+=tt*tt;
        }

        if ((d2min<0.0)||(d2<d2min)) d2min=d2,p1=t;
      } 
 
      //merge particle p1 and p: weight of 'p' is added to 'p1'
      w=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(p1); 

      //PIC::ParticleBuffer::GetV(v1,p1);
      v1[0]=PB::GetVParallel(p1);
      v1[1]=PB::GetVNormal(p1);
      v1[2]=0.0;
      
      
      
      //PIC::ParticleBuffer::GetV(v2,p);
      v2[0]=PB::GetVParallel(p);
      v2[1]=PB::GetVNormal(p);
      v2[2]=0.0;

      for (int idim=0;idim<3;idim++) v1[idim]=(w*v1[idim]+wmin*v2[idim])/(w+wmin);

      //PIC::ParticleBuffer::SetV(v1,p1);  
      PB::SetVParallel(v1[0],p1);
      PB::SetVNormal(v1[1],p1);
      
      PIC::ParticleBuffer::SetIndividualStatWeightCorrection(w+wmin,p1);
      
      //delete particles
      PIC::ParticleBuffer::DeleteParticle(p,Segment->FirstParticleIndex);
      VelocitySpaceTable[ii]->p.erase(itmin);

      if (npart-1>1) more_then_one_particle_left=true;

      if (++cnt==nMaxParticlesReduce) break;
    }

    return cnt;
  }; 

  auto RebuildVelocitySpaceTable = [&] (int CurrentResolutionLevel) {
    cVelocitySpaceElement **NewVelocitySpaceTable;
    int iVelCell,iNewVelCell,iv,jv;


    int NewResolutionLevel=CurrentResolutionLevel-1;

    int nVelCells1d=1<<CurrentResolutionLevel;
    int nNewVelCells1d=nVelCells1d/2;

    int nVelCells=nVelCells1d*nVelCells1d;
    int nNewVelCells=nNewVelCells1d*nNewVelCells1d;


    NewVelocitySpaceTable=new cVelocitySpaceElement* [nNewVelCells];
    for (int i=0;i<nNewVelCells;i++) NewVelocitySpaceTable[i]=NULL;

    for (int i=0;i<nVelCells;i++) {
      if (VelocitySpaceTable[i]!=NULL) {
        iv=VelocitySpaceTable[i]->iv;
        jv=VelocitySpaceTable[i]->jv;

        iVelCell=iv+nVelCells1d*jv;  

        iv/=2,jv/=2;
        iNewVelCell=iv+nNewVelCells1d*jv;

        if (NewVelocitySpaceTable[iNewVelCell]==NULL) {
          NewVelocitySpaceTable[iNewVelCell]=new cVelocitySpaceElement [1];

          *NewVelocitySpaceTable[iNewVelCell]=*VelocitySpaceTable[i];

          NewVelocitySpaceTable[iNewVelCell]->iv/=2;
          NewVelocitySpaceTable[iNewVelCell]->jv/=2;
        }
        else {
          NewVelocitySpaceTable[iNewVelCell]->p.insert(NewVelocitySpaceTable[iNewVelCell]->p.begin(),
              VelocitySpaceTable[i]->p.begin(),VelocitySpaceTable[i]->p.end());
        }
      }
    }

    //delete the previous table 
    for (int i=0;i<nVelCells;i++) if (VelocitySpaceTable[i]!=NULL) delete [] VelocitySpaceTable[i];

    delete [] VelocitySpaceTable;
    VelocitySpaceTable=NewVelocitySpaceTable;
  }; 


  auto DeleteVelocitySpaceTable = [&] (int CurrentResolutionLevel) {
    int nVelCells1d=1<<CurrentResolutionLevel;
    int nVelCells=nVelCells1d*nVelCells1d;

    for (int i=0;i<nVelCells;i++) if (VelocitySpaceTable[i]!=NULL) delete [] VelocitySpaceTable[i];

    delete [] VelocitySpaceTable;
  };

  auto ReduceParticleNumber1 = [&] (int spec,int nModelParticles,FL::cFieldLineSegment *Segment,int nRequestedParticles)  {
    double vmin[3],vmax[3],dv[3]; 
    bool more_then_one_particle_left;


    if (nModelParticles>nRequestedParticles) {
      int nDeleteParticles=nModelParticles-nRequestedParticles;  
      int nRemovedParticles=0;


      static int ncall=0;
      ncall++;

      GetVelocityRange(spec,vmin,vmax,Segment);

      for (int idim=0;idim<3;idim++) {
        double d=vmax[idim]-vmin[idim]; 

        if (d==0.0) d=1.0;

        vmin[idim]-=0.05*d; 
        vmax[idim]+=0.05*d;

        dv[idim]=(vmax[idim]-vmin[idim])/nVelCells1d;
      }


      InitVelocitySpaceTable(spec,vmin,dv,Segment);

      int CurrentResolutionLevel=nSerachLevels;

      while (nRemovedParticles<nDeleteParticles) {
        nRemovedParticles+=RemoveOneParticle(CurrentResolutionLevel,nDeleteParticles-nRemovedParticles,Segment,more_then_one_particle_left); 

        if ((nRemovedParticles<nDeleteParticles)&&(more_then_one_particle_left==false)) {
          RebuildVelocitySpaceTable(CurrentResolutionLevel);
          --CurrentResolutionLevel;
        } 
      }

      DeleteVelocitySpaceTable(CurrentResolutionLevel);
    }
  };

  auto ReduceParticleNumber = [&] (int spec,FL::cFieldLineSegment *Segment,int nRequestedParticles)  {
    vector<cParticleDescriptor> ParticleList;  
    cParticleDescriptor t;
    double *v,w_max=0.0,w_min=-1.0,SummedWeight=0.0,w;
    int nModelParticles=0;
    long int p;

    p=Segment->FirstParticleIndex;

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
          PIC::ParticleBuffer::DeleteParticle(ParticleList[ip].p,Segment->FirstParticleIndex);

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




  auto IncreseParticleNumber = [&] (int spec,FL::cFieldLineSegment *Segment) {
    double MeanV[3]={0.0,0.0,0.0},v[3],w;
    double MeanV2[3]={0.0,0.0,0.0};

    vector<cParticleDescriptor> ParticleList;
    cParticleDescriptor t;
    double w_max=0.0,SummedWeight=0.0,ThermalSpeed=0.0;
    int idim,nModelParticles=0;
    long int p;

    p=Segment->FirstParticleIndex;

    while (p!=-1) {
      if (spec==PIC::ParticleBuffer::GetI(p)) {
        t.p=p;
        t.w=PIC::ParticleBuffer::GetIndividualStatWeightCorrection(p);

        if (t.w>WeightSplittingLimit) {
          ParticleList.push_back(t); 

          SummedWeight+=t.w;
          if (w_max<t.w) w_max=t.w;

          //PIC::ParticleBuffer::GetV(v,p); 
          v[0]=PB::GetVParallel(p);
          v[1]=PB::GetVNormal(p);
          v[2]=0.0;

          for (int idim=0;idim<3;idim++) {
            MeanV[idim]+=t.w*v[idim];
            MeanV2[idim]+=t.w*v[idim]*v[idim];
          }

          nModelParticles++;
        }
      }

      p=PIC::ParticleBuffer::GetNext(p);
    }

    if (nModelParticles<2) return;

    for (idim=0;idim<3;idim++) ThermalSpeed+=MeanV2[idim]/SummedWeight-MeanV[idim]*MeanV[idim]/(SummedWeight*SummedWeight); 

    ThermalSpeed=sqrt(fabs(ThermalSpeed));

    //add new particles in the system
    int nNewParticles=particle_num_limit_max-nModelParticles;
    int ip;
    long int pnew;
    double l[3],vnew[3],vold[3];
    bool flag;

    for (int ii=0;ii<nNewParticles;ii++) {
      flag=true;

      while (flag==true) {
        ip=(int)(rnd()*nModelParticles);

        int ip_wmax=-1;
        w_max=-1.0;
 
        for (int i=0;i<nModelParticles;i++) {
          if (w_max<ParticleList[i].w) w_max=ParticleList[i].w,ip_wmax=i;
        } 

        ip=ip_wmax;

        if (true) { // (ParticleList[ip].w/w_max>rnd()) { //split particles with larger weight  
          pnew=PIC::ParticleBuffer::GetNewParticle(Segment->FirstParticleIndex); 
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
          Vector3D::Distribution::Uniform(l,0.5*v_shift_max*ThermalSpeed); 

          //get velocity for the new particle 
       //   PIC::ParticleBuffer::GetV(vnew,pnew); 
        vnew[0]=PB::GetVParallel(pnew);
        vnew[1]=PB::GetVNormal(pnew);
        vnew[2]=0.0;

          for (idim=0;idim<3;idim++) {
            vold[idim]=vnew[idim]-l[idim];
            vnew[idim]+=l[idim];
          }

          //PIC::ParticleBuffer::SetV(vnew,pnew);
          PB::SetVParallel(vnew[0],pnew);
          PB::SetVNormal(vnew[1],pnew);

       //   PIC::ParticleBuffer::SetV(vold,ParticleList[ip].p);
          PB::SetVParallel(vold[0],ParticleList[ip].p);
          PB::SetVNormal(vold[1],ParticleList[ip].p);


/*
          double dx[3],*xmin,*xmax,x[3],x0[3];

          xmin=node->xmin;
          

          dx[0]=(node->xmax[0]-xmin[0])/_BLOCK_CELLS_X_;
          dx[1]=(node->xmax[1]-xmin[1])/_BLOCK_CELLS_Y_;
          dx[2]=(node->xmax[2]-xmin[2])/_BLOCK_CELLS_Z_;
*/

/*

 



double r;
bool flag=true;
do {
  flag=true;


    r=pow(rnd(),1.0/3.0);

  
  if (apply_non_uniform_x_shift==true) if (r>rnd()) {
    flag=false;
    continue;
  }

  PIC::ParticleBuffer::GetX(x0,pnew);

  double ll[3];

  Vector3D::Distribution::Uniform(ll,x_shift_max*r*sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]));

  x0[0]+=ll[0];
  if (x0[0]<xmin[0]+i*dx[0]) {flag=false;continue;}
  if (x0[0]>xmin[0]+(i+1)*dx[0]) {flag=false;continue;} 

  x0[1]+=ll[1];
  if (x0[1]<xmin[1]+j*dx[1]) {flag=false;continue;}
  if (x0[1]>xmin[1]+(j+1)*dx[1]) {flag=false;continue;}

  x0[2]+=ll[2];
  if (x0[2]<xmin[2]+k*dx[2]) {flag=false;continue;}
  if (x0[2]>xmin[2]+(k+1)*dx[2]) {flag=false;continue;}
}
while (flag==false);



PIC::ParticleBuffer::SetX(x0,pnew);*/



        }
      }
    }
  };


  //loop through the nodes
  int ParticleNumberTable[PIC::nTotalSpecies],spec; 
  
  
  
  for (int iFieldLine=0; iFieldLine<FL::nFieldLine; iFieldLine++) {  
    FL::cFieldLineSegment* Segment;
    int iSegment;
    
    for (iSegment=0,Segment=FL::FieldLinesAll[iFieldLine].GetFirstSegment(); iSegment<FL::FieldLinesAll[iFieldLine].GetTotalSegmentNumber(); iSegment++,Segment=Segment->GetNext()) {
      GetParticleNumber(ParticleNumberTable,Segment);

      for (spec=0;spec<PIC::nTotalSpecies;spec++) {
        if (ParticleNumberTable[spec]>particle_num_limit_max) {
          ReduceParticleNumber1(spec,ParticleNumberTable[spec],Segment,particle_num_limit_min+0.9*(particle_num_limit_max-particle_num_limit_min));
        }
        else if ((1<ParticleNumberTable[spec])&&(ParticleNumberTable[spec]<particle_num_limit_min)) {
          IncreseParticleNumber(spec,Segment);
        }
      }
    }
    
  }
  

/*  for (inode=0;inode<PIC::DomainBlockDecomposition::nLocalBlocks;inode++) {
    node=PIC::DomainBlockDecomposition::BlockTable[inode];

    if (node->block!=NULL) {
      for (i=0;i<_BLOCK_CELLS_X_;i++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (k=0;k<_BLOCK_CELLS_Z_;k++) {
        GetParticleNumber(ParticleNumberTable,i,j,k,node);

        for (spec=0;spec<PIC::nTotalSpecies;spec++) {
          if (ParticleNumberTable[spec]>particle_num_limit_max) {
            ReduceParticleNumber1(spec,ParticleNumberTable[spec],i,j,k,node,particle_num_limit_min+0.9*(particle_num_limit_max-particle_num_limit_min));
          }
          else if ((1<ParticleNumberTable[spec])&&(ParticleNumberTable[spec]<particle_num_limit_min)) {
            IncreseParticleNumber(spec,i,j,k,node);
          }
        }
      }
    }      
  }*/

//  PIC::Parallel::ExchangeParticleData();
} 





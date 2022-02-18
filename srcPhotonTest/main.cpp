/*
 * main.cpp
 *
 *  Created on: Jun 21, 2012
 *      Author: fougere and vtenishe
 */

//$Id$


#include "pic.h"
#include "constants.h"

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <list>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <algorithm>
#include <sys/time.h>
#include <sys/resource.h>
#include <ctime>

#include "meshAMRcutcell.h"
#include "cCutBlockSet.h"
#include "meshAMRgeneric.h"

#include "radiation.h"

#include "../../srcInterface/LinearSystemCornerNode.h"
#include "linear_solver_wrapper_c.h"

#include "PeriodicBCTest.dfn"


//for lapenta mover

#include "pic.h"
#include "Exosphere.dfn"
#include "Exosphere.h"


int nVars=3; //number of variables in center associated data
double Background[3]={100.0,-20.0,10.0};


//#define _UNIFORM_MESH_ 1
//#define _NONUNIFORM_MESH_ 2

#ifndef _TEST_MESH_MODE_
#define _TEST_MESH_MODE_ _UNIFORM_MESH_
#endif


double xmin[3]={-16.0,-8.0,-4.0};
double xmax[3]={16.0,8.0,4.0};

int CurrentCenterNodeOffset=-1,NextCenterNodeOffset=-1;
int CurrentCornerNodeOffset=-1,NextCornerNodeOffset=-1;

int iCase;

void CleanParticles(){
  
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;

  for (node=PIC::Mesh::mesh->BranchBottomNodeList;node!=NULL;node=node->nextBranchBottomNode) if (node->block!=NULL) {
   
     long int *  FirstCellParticleTable=node->block->FirstCellParticleTable;
     if (FirstCellParticleTable==NULL) continue;
     for (int k=0;k<_BLOCK_CELLS_Z_;k++) {
       for (int j=0;j<_BLOCK_CELLS_Y_;j++)  {
	 for (int i=0;i<_BLOCK_CELLS_X_;i++) {
	   long int * ptr=FirstCellParticleTable+(i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k));
	   while ((*ptr)!=-1) PIC::ParticleBuffer::DeleteParticle(*ptr,*ptr);

//////
/*
long int next,ptr=FirstCellParticleTable[(i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k))]; 

while (ptr!=-1) {

  next=PIC::ParticleBuffer::GetNext(ptr);
  PIC::ParticleBuffer::DeleteParticle(ptr);
  ptr=next; 
}

FirstCellParticleTable[(i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k))]=-1;
*/

/////	   
	 }   
       }
     }
     
  }

}


long int PrepopulateDomain() {
  int iCell,jCell,kCell;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
  PIC::Mesh::cDataCenterNode *cell;
  long int nd,nGlobalInjectedParticles,nLocalInjectedParticles=0;
  double Velocity[3];
  /*
  //local copy of the block's cells
  int cellListLength=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread]->block->GetCenterNodeListLength();
  PIC::Mesh::cDataCenterNode *cellList[cellListLength];
  */
  //particle ejection parameters
  double ParticleWeight;//beta=PIC::MolecularData::GetMass(spec)/(2*Kbol*Temperature);
  double waveNumber[3]={0.0,0.0,0.0};
  double lambda=32.0;
 
  waveNumber[0]=2*Pi/lambda;
  
  int nBlock[3]={_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};

  //the boundaries of the block and middle point of the cell
  double *xminBlock,*xmaxBlock;
  double v[3],anpart;
  int npart;
  char * offset=NULL;
  int ionSpec=0, electronSpec=1;
  double ionMass = PIC::MolecularData::GetMass(ionSpec)/_AMU_;
  double electronMass = PIC::MolecularData::GetMass(electronSpec)/_AMU_;
  //  for (node=PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::Mesh::mesh->ThisThread];node!=NULL;node=node->nextNodeThisThread) {
  //  {
        for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
      
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> * node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];
    if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_) {
      bool BoundaryBlock=false;
      
      for (int iface=0;iface<6;iface++) if (node->GetNeibFace(iface,0,0,PIC::Mesh::mesh)==NULL) {
	  //the block is at the domain boundary, and thresefor it is a 'ghost' block that is used to impose the periodic boundary conditions
	  BoundaryBlock=true;
	  break;
	}
      
      if (BoundaryBlock==true) continue;
    }

    if (node->Thread!=PIC::ThisThread) continue;


    // }

    // PIC::Mesh::cDataCenterNode *cellList[cellListLength];
  
    //memcpy(cellList,node->block->GetCenterNodeList(),cellListLength*sizeof(PIC::Mesh::cDataCenterNode*));

    xminBlock=node->xmin,xmaxBlock=node->xmax;
    double dx[3];
    double CellVolume=1;
    for (int idim=0;idim<3;idim++) {
      dx[idim]=(xmaxBlock[idim]-xminBlock[idim])/nBlock[idim];
      CellVolume *= dx[idim];
    }
    //particle stat weight
#ifndef _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
#error ERROR: _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_ is used but not defined
#endif
#ifndef _SIMULATION_PARTICLE_WEIGHT_MODE_
#error ERROR: _SIMULATION_PARTICLE_WEIGHT_MODE_ is used but not defined
#endif

    //assume ion and electron have the same particle weight
    #if  _SIMULATION_PARTICLE_WEIGHT_MODE_ == _SPECIES_DEPENDENT_GLOBAL_PARTICLE_WEIGHT_
    ParticleWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[ionSpec];
    #else
    ParticleWeight=node->block->GetLocalParticleWeight(ionSpec);
    #endif

    for (kCell=0;kCell<nBlock[2];kCell++) for (jCell=0;jCell<nBlock[1];jCell++) for (iCell=0;iCell<nBlock[0];iCell++) {
	  //      nd=PIC::Mesh::mesh->getCenterNodeLocalNumber(iCell,jCell,kCell);

      // cell=cellList[nd];
      //  xMiddle=cell->GetX();
      //offset = cell->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
	  int ind[3]={iCell,jCell,kCell};
	  double x[3];
	  for (int idim=0; idim<3; idim++) x[idim]=xminBlock[idim]+(ind[idim]+0.5)*dx[idim];
	  
          
          //double waveShape =sin(waveNumber[0]*(x[0]-xmin[0])+waveNumber[1]*(x[1]-xmin[1])+waveNumber[2]*(x[2]-xmin[2]));
          double waveShape=1.0;
          double rho = 1, ux=0.0, p=4.5e-4, Ppar=4.5e-4;
          double Ppar1 = 4.5e-5*waveShape;
          double p1 = 7.5e-5*waveShape;
          double rho1 = 0.1*waveShape;
          double rho_conv=0.0795774715459477;
          double p_conv = 0.0795774715459477;
          
          rho *=rho_conv;
          rho1 *= rho_conv;
          p *= p_conv;
          p1 *= p_conv;
          
          Ppar *=p_conv;
          Ppar1 *=p_conv;
          

          double rho_i = (rho+rho1)*ionMass/(ionMass+electronMass);
          double rho_e = (rho+rho1)*electronMass/(ionMass+electronMass); 
          double ux1= 0.005*waveShape;



          double NumberDensity=(rho+rho1)/(ionMass+electronMass);
          double kTemp_par = (Ppar+Ppar1)/NumberDensity*0.5;
          double kTemp = (p+p1)/NumberDensity*0.5;
          double kTemp_perp = 3*kTemp-2*kTemp_par;

          
          double pi_par=(Ppar+Ppar1)*0.5;
          double pe_par=pi_par;
          double pi_perp = 3*(p+p1)*0.5-2*pi_par;
          double pe_perp = pi_perp;
          double ionBulkVelocity[3]={0,0,0};
          double electronBulkVelocity[3]={0,0,0};
          
          ionBulkVelocity[0] = ux+ux1;
          electronBulkVelocity[0] = ux+ux1;

          //inject particles into the cell
          anpart=NumberDensity*CellVolume/ParticleWeight;
          //std::cout<<"CellLoc:"<<x[0]<<" "<<x[1]<<" "<<x[2]<<" NumberDensity: "<<NumberDensity<<"cell volume: "<<CellVolume<<"anpart: "<<anpart<<std::endl;
          npart=(int)(anpart);
          //if (rnd()<anpart-npart) npart++;
          nLocalInjectedParticles+=npart*2;
          //std::cout<<"need to inject npart: "<<npart<<std::endl;
          
          while (npart-->0) {
            double xPar[3];
            xPar[0]=x[0]+dx[0]*(rnd()-0.5);
            xPar[1]=x[1]+dx[1]*(rnd()-0.5);
            
            // xPar[0]=x[0];
            // xPar[1]=x[1];
            xPar[2]=x[2];

            
            double electronVelocity[3],ionVelocity[3];
            for (int idim=0;idim<3;idim++) { 
              double uth_e = idim!=1?sqrt(pe_perp/rho_e):sqrt(pe_par/rho_e);
              double uth_i = idim!=1?sqrt(pi_perp/rho_i):sqrt(pi_par/rho_i);
              
              electronVelocity[idim]=uth_e* sqrt(-2.0 * log(1.0 - .999999999 * rnd()))*cos(2*Pi*rnd())+electronBulkVelocity[idim];
              ionVelocity[idim]=uth_i*sqrt(-2.0 * log(1.0 - .999999999 * rnd()))*cos(2*Pi*rnd())+ionBulkVelocity[idim];   
            }
            
            /*  
            for (int idim=0;idim<3;idim++) {
              //in this test case B field is in y-direction
              double ElectronTemp= idim!=1?kTemp_perp/electronMass:kTemp_par/electronMass; 
              double IonTemp= idim!=1?kTemp_perp/ionMass:kTemp_par/ionMass; 
              
              electronVelocity[idim]=cos(2*Pi*rnd())*sqrt(-log(rnd())*(2*ElectronTemp))+electronBulkVelocity[idim];
              ionVelocity[idim]=cos(2*Pi*rnd())*sqrt(-log(rnd())*(2*IonTemp))+ionBulkVelocity[idim];       
              }
            */      
            //initiate the new particle
            PIC::ParticleBuffer::InitiateParticle(xPar, electronVelocity,NULL,&electronSpec,NULL,_PIC_INIT_PARTICLE_MODE__ADD2LIST_,(void*)node);
            PIC::ParticleBuffer::InitiateParticle(xPar, ionVelocity,NULL,&ionSpec,NULL,_PIC_INIT_PARTICLE_MODE__ADD2LIST_,(void*)node);
            
          }
      //end of the particle injection block
      //std::cout<<"finished injecting npart: "<<npart<<std::endl;
        }
        }

  MPI_Allreduce(&nLocalInjectedParticles,&nGlobalInjectedParticles,1,MPI_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  printf("particles prepopulated!\n");
  return nGlobalInjectedParticles;
}


/*
void SetIC() {
  
    int i,j,k;
    char *offset;
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *node;
    double cPi = 3.14159265;
    double waveNumber[3]={0.0,0.0,0.0};
    double lambda=32.0;
   
    waveNumber[0]=2*cPi/lambda;
  
    double x[3];
   
    using namespace PIC::FieldSolver::Electromagnetic::ECSIM;
    int nBreak=0;

    printf("User Set IC called\n");

    int nBlocks[3] ={_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};
    for (int nLocalNode=0;nLocalNode<PIC::DomainBlockDecomposition::nLocalBlocks;nLocalNode++) {
      node=PIC::DomainBlockDecomposition::BlockTable[nLocalNode];
      
      if (_PIC_BC__PERIODIC_MODE_==_PIC_BC__PERIODIC_MODE_ON_) {
	bool BoundaryBlock=false;
	
	for (int iface=0;iface<6;iface++) if (node->GetNeibFace(iface,0,0)==NULL) {
	    //the block is at the domain boundary, and thresefor it is a 'ghost' block that is used to impose the periodic boundary conditions
	    BoundaryBlock=true;
	    break;
	  }
	
	if (BoundaryBlock==true) continue;
      }
      
      double dx[3];
      double *xminBlock= node->xmin, *xmaxBlock= node->xmax;

      for (int idim=0;idim<3;idim++) dx[idim]=(xmaxBlock[0]-xminBlock[0])/nBlocks[idim];
      
     
      for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
	    
            int ind[3]={i,j,k};
            
	    PIC::Mesh::cDataCornerNode *CornerNode= node->block->GetCornerNode(PIC::Mesh::mesh->getCornerNodeLocalNumber(i,j,k));
	    if (CornerNode!=NULL){
	      offset=CornerNode->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;
	     
              
	      for (int idim=0; idim<3; idim++) x[idim]=xminBlock[idim]+ind[idim]*dx[idim];
	     
              ((double*)(offset+CurrentEOffset))[ExOffsetIndex]=0.0;
              ((double*)(offset+CurrentEOffset))[EyOffsetIndex]=0.0;
              ((double*)(offset+CurrentEOffset))[EzOffsetIndex]=0.0;
	  
              ((double*)(offset+OffsetE_HalfTimeStep))[ExOffsetIndex]=0.0;
              ((double*)(offset+OffsetE_HalfTimeStep))[ExOffsetIndex]=0.0;
              ((double*)(offset+OffsetE_HalfTimeStep))[ExOffsetIndex]=0.0;
          

	    // ((double*)(offset+CurrentCornerNodeOffset))[EzOffsetIndex]=i+j*_BLOCK_CELLS_X_+k*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_+nLocalNode*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;


	    
	    //((double*)(offset+NextCornerNodeOffset))[EzOffsetIndex]=i+j*_BLOCK_CELLS_X_+k*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_+nLocalNode*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;
	    }//
	  }//for (k=0;k<_BLOCK_CELLS_Z_+1;k++) for (j=0;j<_BLOCK_CELLS_Y_+1;j++) for (i=0;i<_BLOCK_CELLS_X_+1;i++) 
      for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) {
            int ind[3]={i,j,k};
	    PIC::Mesh::cDataCenterNode *CenterNode= node->block->GetCenterNode(PIC::Mesh::mesh->getCenterNodeLocalNumber(i,j,k));
	    if (CenterNode!=NULL){
	      offset=node->block->GetCenterNode(PIC::Mesh::mesh->getCenterNodeLocalNumber(i,j,k))->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;

              for (int idim=0; idim<3; idim++) x[idim]=xminBlock[idim]+(ind[idim]+0.5)*dx[idim];

	      //double B = 0.004*sin(waveNumber[0]*(x[0]-xmin[0])+waveNumber[1]*(x[1]-xmin[1])+waveNumber[2]*(x[2]-xmin[2]));
              //double By = 0.04;
              double B=0.0, By=0.0;
              double Bfactor =1.120998243279586E-3*892.062058076386;
              ((double*)(offset+CurrentBOffset))[BxOffsetIndex]=0.0;
              ((double*)(offset+CurrentBOffset))[ByOffsetIndex]=(B+By)*Bfactor;
              ((double*)(offset+CurrentBOffset))[BzOffsetIndex]=0.0;
		
		
              ((double*)(offset+PrevBOffset))[BxOffsetIndex]=0.0;
              ((double*)(offset+PrevBOffset))[ByOffsetIndex]=(B+By)*Bfactor;
              ((double*)(offset+PrevBOffset))[BzOffsetIndex]=0.0;
            }// if (CenterNode!=NULL)
	  }//for (k=0;k<_BLOCK_CELLS_Z_;k++) for (j=0;j<_BLOCK_CELLS_Y_;j++) for (i=0;i<_BLOCK_CELLS_X_;i++) 
    }
   
    switch (_PIC_BC__PERIODIC_MODE_) {
    case _PIC_BC__PERIODIC_MODE_OFF_:
      PIC::Mesh::mesh->ParallelBlockDataExchange();
      break;
      
    case _PIC_BC__PERIODIC_MODE_ON_:
      PIC::BC::ExternalBoundary::UpdateData();
      break;
    }
}
*/

double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
  double res,CellSize;

  CellSize=startNode->GetCharacteristicCellSize();
  res=0.01*CellSize/Radiation::SpeedOfLight_cm;


  return (res>1.0) ? 1.0 : res; 
}


double BulletLocalResolution(double *x) {                                                                                           
  double res;

  res=0.01;
  return res;
}
                       

int main(int argc,char **argv) {
  
  time_t TimeValue=time(NULL);
  tm *ct=localtime(&TimeValue);
  
  printf("start: (%i/%i %i:%i:%i)\n",ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec);

  Radiation::Init();

  PIC::InitMPI();
  PIC::Init_BeforeParser();
  
  //seed the random number generator
  rnd_seed(100);

  double xmin[3]={0.0,0.0,0.0},xmax[3]={0.2,0.05,0.05};

  //generate mesh or read from file
  PIC::Mesh::mesh->AllowBlockAllocation=false;

  PIC::Mesh::mesh->init(xmin,xmax,BulletLocalResolution);
  PIC::Mesh::mesh->memoryAllocationReport();

  //generate mesh or read from file
  if (PIC::Mesh::mesh->ThisThread==0) {
     PIC::Mesh::mesh->buildMesh();
     PIC::Mesh::mesh->saveMeshFile("mesh.msh");
     MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  }
  else {
     MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
     PIC::Mesh::mesh->readMeshFile("mesh.msh");
   }
  
  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);

  PIC::Mesh::initCellSamplingDataBuffer();

  PIC::Mesh::mesh->CreateNewParallelDistributionLists();

  PIC::Mesh::mesh->AllowBlockAllocation=true;
  PIC::Mesh::mesh->AllocateTreeBlocks();
  PIC::Mesh::mesh->InitCellMeasure();

  //create the list of mesh nodes where the injection boundary conditinos are applied
  PIC::BC::BlockInjectionBCindicatior=Radiation::Injection::BoundingBoxParticleInjectionIndicator;
  PIC::BC::userDefinedBoundingBlockInjectionFunction=Radiation::Injection::BoundingBoxInjection;
  PIC::BC::InitBoundingBoxInjectionBlockList();

  PIC::Init_AfterParser();
  PIC::Mover::Init();

  //set up the time step
  PIC::ParticleWeightTimeStep::LocalTimeStep=localTimeStep;
  PIC::ParticleWeightTimeStep::initTimeStep();

  if (PIC::ThisThread==0) printf("test1\n");
  PIC::Mesh::mesh->outputMeshTECPLOT("mesh_test.dat");
  
  if(_PIC_BC__PERIODIC_MODE_== _PIC_BC__PERIODIC_MODE_ON_){
  PIC::BC::ExternalBoundary::Periodic::InitBlockPairTable();
  }

  int s,i,j,k;
  if (PIC::ThisThread==0) printf("test2\n");
 
  PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(0,1e-2*0.0795774715459477*10*0.001*0.01*0.001*0.1);

  PIC::DomainBlockDecomposition::UpdateBlockTable();

  Radiation::IC::Set();
  PIC::Mesh::mesh->outputMeshDataTECPLOT("ic.dat",0);
  
  int totalIter=300000;

  if (_PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_) totalIter=100;


  Radiation::ClearCellCounters();

  for (int niter=0;niter<totalIter;niter++) {

if (PIC::ThisThread==0) cout << "t=" << niter*PIC::ParticleWeightTimeStep::GlobalTimeStep[0] << "[ns], pass=" << niter*PIC::ParticleWeightTimeStep::GlobalTimeStep[0]*Radiation::SpeedOfLight_cm << "[cm]"<<endl;  

    PIC::TimeStep();
    Radiation::Emission();
  }
  
  PIC::RunTimeSystemState::CumulativeTiming::Print();
  MPI_Finalize();
  
  TimeValue=time(NULL);
  ct=localtime(&TimeValue);
  
  printf("end: (%i/%i %i:%i:%i)\n",ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec);

  cout << "End of the run" << endl;
  return EXIT_SUCCESS;
}

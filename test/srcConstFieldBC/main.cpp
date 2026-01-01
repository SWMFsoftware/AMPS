/*
 * main.cpp
 *
 *  Created on: Jun 21, 2012
 *      Author: fougere and vtenishe
 */

/***************************************************************************************
 * Test driver for PIC with/without particles and constant field boundary conditions
 *
 * OVERVIEW
 * --------
 * This test can run in four modes controlled via the command line:
 *
 *   1) -particles
 *        • PIC with particles enabled.
 *        • Field boundary conditions are held constant.
 *        • Particle boundary condition is *specular reflection* at domain faces.
 *        • Initial fields: by default E = 0, B = 0 unless you also set a background.
 *
 *   2) -no-particles
 *        • Field-only (no particles).
 *        • Field boundary conditions are held constant.
 *        • Initial fields default to E = 0, B = 0 (unless -B/-E are provided).
 *
 *   3) -B [Bx By Bz]
 *        • Field-only (implies -no-particles).
 *        • Initialize B to a uniform vector; E = 0.
 *        • If the three numbers are omitted, defaults to B = (0, 1, 0).
 *
 *   4) -E [Ex Ey Ez]
 *        • Field-only (implies -no-particles).
 *        • Initialize E to a uniform vector; B = 0.
 *        • If the three numbers are omitted, defaults to E = (1, 0, 0).
 *
 * Additionally:
 *   -stencil-order=N
 *        • Sets the finite-difference stencil order used by divergence/curl/Poisson
 *          operators in this test (typical choices: 2,4,6,8). Default is 2.
 *
 * BOUNDARY CONDITIONS
 * -------------------
 * • Fields: “constant BC” — the field values at the domain boundary are held fixed to
 *   the values chosen by the mode (e.g., the uniform E or B you set). In practice this
 *   can be enforced via your project’s boundary callback or by reapplying the uniform
 *   values to halo layers after each step in field-only runs.
 * • Particles: “specular reflection” — when a particle hits a domain face, its velocity
 *   component normal to that face is flipped (tangential components are preserved).
 *
 * HOW SetIC() WORKS NOW
 * ---------------------
 * SetIC(cfg) zeros all fields, then applies the requested initial condition:
 *   • -particles: leaves E and B at your default (often zeros) unless your test adds a
 *     background; particles are (pre)populated later as usual.
 *   • -no-particles: leaves both fields zero.
 *   • -B: sets cell-centered B uniformly to (Bx,By,Bz), keeps corner E = 0.
 *   • -E: sets corner E uniformly to (Ex,Ey,Ez), keeps cell-centered B = 0.
 *
 * EXAMPLES
 * --------
 *   # PIC with particles, specular reflection, 4th-order stencil
 *   ./amps_test -particles -stencil-order=4
 *
 *   # Field-only, uniform B=(0,1,0), default 2nd-order
 *   ./amps_test -B
 *
 *   # Field-only, E=(0.2,0,0), 8th-order stencil
 *   ./amps_test -E 0.2 0.0 0.0 -stencil-order 8
 *
 *   # Field-only, all-zero fields, 6th-order stencil
 *   ./amps_test -no-particles -stencil-order=6
 *
 ***************************************************************************************/


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

#include "../../srcInterface/LinearSystemCornerNode.h"
#include "linear_solver_wrapper_c.h"

//#include "PeriodicBCTest.dfn"

#if _CUDA_MODE_ == _ON_
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#endif

//for lapenta mover

#include "pic.h"
#include "Exosphere.dfn"
#include "Exosphere.h"


int nVars=3; //number of variables in center associated data
double Background[3]={100.0,-20.0,10.0};



double xmin[3]={-16.0,-8.0,-4.0};
double xmax[3]={16.0,8.0,4.0};

int CurrentCenterNodeOffset=-1,NextCenterNodeOffset=-1;
int CurrentCornerNodeOffset=-1,NextCornerNodeOffset=-1;

int iCase;


#include "main_lib.h"

// Global configuration for this test (used by SetIC/PrepopulateDomain in main_lib.cpp)
TestConfig cfg;

// Optional: export stencil order for helper operators in this test module
int g_TestStencilOrder = 2;

int main(int argc,char **argv) {
  
   time_t TimeValue=time(NULL);
   tm *ct=localtime(&TimeValue);
  
   printf("start: (%i/%i %i:%i:%i)\n",ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec);

  ConfigureTestFromArgsWithInput(cfg,argc,argv);
  g_TestStencilOrder = cfg.stencilOrder;

// If -L was provided, redefine the domain to be centered at (0,0,0).
// NOTE: xmin/xmax are used for mesh initialization and (if periodic) for periodic BC setup.
if (cfg.use_domain_L) {
  for (int d=0; d<3; ++d) {
    xmin[d] = -0.5*cfg.domain_L[d];
    xmax[d] =  0.5*cfg.domain_L[d];
  }
}



  // Dirichlet on all 6 faces
  PIC::FieldSolver::Electromagnetic::DomainBC.SetAll(PIC::Mesh::BCTypeDirichlet);

  PIC::InitMPI();
  PIC::Init_BeforeParser();

  
  int RelativeOffset=0;
  
#ifndef _NONUNIFORM_MESH_
#error ERROR: _NONUNIFORM_MESH_ is used but not defined
#endif
#ifndef _TEST_MESH_MODE_
#error ERROR: _TEST_MESH_MODE_ is used but not defined
#endif
#if _TEST_MESH_MODE_==_NONUNIFORM_MESH_
  printf("non-uniform mesh!\n");
#endif
#ifndef _UNIFORM_MESH_
#error ERROR: _UNIFORM_MESH_ is used but not defined
#endif
#ifndef _TEST_MESH_MODE_
#error ERROR: _TEST_MESH_MODE_ is used but not defined
#endif
#if _TEST_MESH_MODE_==_UNIFORM_MESH_
  printf("uniform mesh!\n");
#endif


#ifndef _PIC_MODE_ON_
#error ERROR: _PIC_MODE_ON_ is used but not defined
#endif
#ifndef _CURRENT_MODE_
#error ERROR: _CURRENT_MODE_ is used but not defined
#endif
#if _CURRENT_MODE_==_PIC_MODE_ON_
  printf("current on!\n");
#endif
#ifndef _PIC_MODE_OFF_
#error ERROR: _PIC_MODE_OFF_ is used but not defined
#endif
#ifndef _CURRENT_MODE_
#error ERROR: _CURRENT_MODE_ is used but not defined
#endif
#if _CURRENT_MODE_==_PIC_MODE_OFF_
  printf("current mode off!\n");
#endif



  //seed the random number generator
  rnd_seed(100);

  //generate mesh or read from file
  char mesh[_MAX_STRING_LENGTH_PIC_]="none";  ///"amr.sig=0xd7058cc2a680a3a2.mesh.bin";
  sprintf(mesh,"amr.sig=%s.mesh.bin","test_mesh");

  PIC::Mesh::mesh->AllowBlockAllocation=false;
  if(_PIC_BC__PERIODIC_MODE_== _PIC_BC__PERIODIC_MODE_ON_){
  PIC::BC::ExternalBoundary::Periodic::Init(xmin,xmax,BulletLocalResolution);
  }else{
    PIC::Mesh::mesh->init(xmin,xmax,BulletLocalResolution);
  }
  PIC::Mesh::mesh->memoryAllocationReport();

  //generate mesh or read from file
  bool NewMeshGeneratedFlag=false;

  char fullname[STRING_LENGTH];
  sprintf(fullname,"%s/%s",PIC::UserModelInputDataPath,mesh);

  FILE *fmesh=NULL;

  fmesh=fopen(fullname,"r");

  if (fmesh!=NULL) {
    fclose(fmesh);
    PIC::Mesh::mesh->readMeshFile(fullname);
  }
  else {
    NewMeshGeneratedFlag=true;

    if (PIC::Mesh::mesh->ThisThread==0) {
       PIC::Mesh::mesh->buildMesh();
       PIC::Mesh::mesh->saveMeshFile("mesh.msh");
       MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    }
    else {
       MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
       PIC::Mesh::mesh->readMeshFile("mesh.msh");
    }
  }


  //if the new mesh was generated => rename created mesh.msh into amr.sig=0x%lx.mesh.bin
  if (NewMeshGeneratedFlag==true) {
    unsigned long MeshSignature=PIC::Mesh::mesh->getMeshSignature();

    if (PIC::Mesh::mesh->ThisThread==0) {
      char command[300];

      sprintf(command,"mv mesh.msh amr.sig=0x%lx.mesh.bin",MeshSignature);
      system(command);
    }
  }

  MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
  

  PIC::Mesh::initCellSamplingDataBuffer();

  PIC::Mesh::mesh->CreateNewParallelDistributionLists();

  PIC::Mesh::mesh->AllowBlockAllocation=true;
  PIC::Mesh::mesh->AllocateTreeBlocks();
  PIC::Mesh::mesh->InitCellMeasure();

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
  //-387.99e2
  int s,i,j,k;


  if (PIC::ThisThread==0) printf("test2\n");
 
  // PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(0);
  //PIC::ParticleWeightTimeStep::initParticleWeight_ConstantWeight(1);

    // Initialize global particle weight to target ~ppc/spec for the uniform particle IC.
  InitGlobalParticleWeight_TargetPPC(cfg);


  //solve the transport equation
  //set the initial conditions for the transport equation
  //  TransportEquation::SetIC(3);
 

  switch (_PIC_BC__PERIODIC_MODE_) {
  case _PIC_BC__PERIODIC_MODE_OFF_:
      PIC::Mesh::mesh->ParallelBlockDataExchange();
      break;
      
  case _PIC_BC__PERIODIC_MODE_ON_:
    PIC::BC::ExternalBoundary::UpdateData();
      break;
  }
  //PIC::FieldSolver::Init(); 
  PIC::FieldSolver::Electromagnetic::ECSIM::SetIC=SetIC;
    
  int  totalIter,CaseNumber;
  //PIC::FieldSolver::Init();
  PIC::FieldSolver::Electromagnetic::ECSIM::Init_IC();
  PIC::CPLR::FLUID::EFieldTol = 1.0e-8;

  totalIter=60;
     
    switch (_PIC_BC__PERIODIC_MODE_) {
    case _PIC_BC__PERIODIC_MODE_OFF_:
      PIC::Mesh::mesh->ParallelBlockDataExchange();
      break;
      
    case _PIC_BC__PERIODIC_MODE_ON_:
      PIC::BC::ExternalBoundary::UpdateData();
      break;
    }
    PIC::Mesh::mesh->outputMeshDataTECPLOT("ic.dat",0);
  

      int LocalParticleNumber=PIC::ParticleBuffer::GetAllPartNum();
      int GlobalParticleNumber;
      MPI_Allreduce(&LocalParticleNumber,&GlobalParticleNumber,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
      printf("Before cleaning, LocalParticleNumber,GlobalParticleNumber,iThread:%d,%d,%d\n",LocalParticleNumber,GlobalParticleNumber,PIC::ThisThread);
      std::cout<<"LocalParticleNumber: "<<LocalParticleNumber<<" GlobalParticleNumber:"<<GlobalParticleNumber<<std::endl;

      CleanParticles();
      LocalParticleNumber=PIC::ParticleBuffer::GetAllPartNum();
      MPI_Allreduce(&LocalParticleNumber,&GlobalParticleNumber,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
      printf("After cleaning, LocalParticleNumber,GlobalParticleNumber,iThread:%d,%d,%d\n",LocalParticleNumber,GlobalParticleNumber,PIC::ThisThread);

      if (cfg.mode==TestConfig::Mode::WithParticles) {
        PrepopulateDomain();
      }

      LocalParticleNumber=PIC::ParticleBuffer::GetAllPartNum();
      MPI_Allreduce(&LocalParticleNumber,&GlobalParticleNumber,1,MPI_INT,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
      printf("After prepopulating, LocalParticleNumber,GlobalParticleNumber,iThread:%d,%d,%d\n",LocalParticleNumber,GlobalParticleNumber,PIC::ThisThread);
      std::cout<<"LocalParticleNumber: "<<LocalParticleNumber<<" GlobalParticleNumber:"<<GlobalParticleNumber<<std::endl;
   
    switch (_PIC_BC__PERIODIC_MODE_) {
    case _PIC_BC__PERIODIC_MODE_OFF_:
      PIC::Mesh::mesh->ParallelBlockDataExchange();
      break;
      
    case _PIC_BC__PERIODIC_MODE_ON_:
      PIC::BC::ExternalBoundary::UpdateData();
      break;
    }
    
    PIC::Sampling::Sampling();

    for (int niter=0;niter<totalIter;niter++) {

if (_CUDA_MODE_ == 12384) { ////_ON_
#if _CUDA_MODE_ == _ON_

      int *ParticlePopulationNumberTable=NULL;

      amps_malloc_managed<int>(ParticlePopulationNumberTable,PIC::DomainBlockDecomposition::nLocalBlocks*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_);

      auto CreateParticlePopulationNumberTable = [=] _TARGET_HOST_ _TARGET_DEVICE_ (int *ParticleNumberTable,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> **BlockTable) { 
        int TableLength=PIC::DomainBlockDecomposition::nLocalBlocks*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;

        //get the thread global id
        #ifdef __CUDA_ARCH__ 
        int id=blockIdx.x*blockDim.x+threadIdx.x;
        int increment=gridDim.x*blockDim.x;
        int  SearchIndexLimit=warpSize*(1+TableLength/warpSize); 
        #else 
        int id=0,increment=1;
        int SearchIndexLimit=TableLength;
        #endif


        for (int icell=id;icell<SearchIndexLimit;icell+=increment) {
          int nLocalNode,ii=icell;
          int i,j,k;
          long int ptr;

          if (icell<TableLength) {
            nLocalNode=ii/(_BLOCK_CELLS_Z_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_);
            ii-=nLocalNode*_BLOCK_CELLS_Z_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_;

            k=ii/(_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_);
            ii-=k*_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_;

            j=ii/_BLOCK_CELLS_X_;
            ii-=j*_BLOCK_CELLS_X_;

            i=ii;

            cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> * node=BlockTable[nLocalNode];
            ParticleNumberTable[icell]=0;

            if (node->block!=NULL) {
              ptr=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];

               while (ptr!=-1) {
                 ParticleNumberTable[icell]++;
                 ptr=PIC::ParticleBuffer::GetNext(ptr);
               }
            }
          }

          #ifdef __CUDA_ARCH__
          __syncwarp();	
          #endif
       }
     }; 


      auto CreateParticlePopulationTable = [=] _TARGET_HOST_ _TARGET_DEVICE_ (long int *ParticlePopulationTable,int *ParticleOffsetTable,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> **BlockTable) {
        int TableLength=PIC::DomainBlockDecomposition::nLocalBlocks*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;

        #ifdef __CUDA_ARCH__ 
        int id=blockIdx.x*blockDim.x+threadIdx.x;
        int increment=gridDim.x*blockDim.x;
        int  SearchIndexLimit=warpSize*(1+TableLength/warpSize);
        #else
        int id=0,increment=1;
        int SearchIndexLimit=TableLength;
        #endif


        for (int icell=id;icell<SearchIndexLimit;icell+=increment) {
          int nLocalNode,ii=icell;
          int i,j,k,offset;
          long int ptr;

          if (icell<TableLength) {
            nLocalNode=ii/(_BLOCK_CELLS_Z_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_);
            ii-=nLocalNode*_BLOCK_CELLS_Z_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_;

            k=ii/(_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_);
            ii-=k*_BLOCK_CELLS_Y_*_BLOCK_CELLS_X_;

            j=ii/_BLOCK_CELLS_X_;
            ii-=j*_BLOCK_CELLS_X_;

            i=ii;

            cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> * node=BlockTable[nLocalNode];

            if (node->block!=NULL) {
              ptr=node->block->FirstCellParticleTable[i+_BLOCK_CELLS_X_*(j+_BLOCK_CELLS_Y_*k)];
              offset=ParticleOffsetTable[icell];

               while (ptr!=-1) {
                 ParticlePopulationTable[offset++]=ptr;
                 ptr=PIC::ParticleBuffer::GetNext(ptr);
               }
            }
          }

          #ifdef __CUDA_ARCH__
          __syncwarp();
          #endif
       }
     };

     kernel_2<<<3,128>>>(CreateParticlePopulationNumberTable,ParticlePopulationNumberTable,PIC::DomainBlockDecomposition::BlockTable);
     cudaDeviceSynchronize();

     int total_number=0;
     
     for (int i=0;i<PIC::DomainBlockDecomposition::nLocalBlocks*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;i++)total_number+=ParticlePopulationNumberTable[i]; 

     if (total_number!=PIC::ParticleBuffer::NAllPart) exit(__LINE__,__FILE__,"Error: the particle number is not consistent");


      long int *ParticlePopulationTable=NULL;
      int *ParticleOffsetNumber;

      amps_malloc_managed<long int>(ParticlePopulationTable,PIC::ParticleBuffer::NAllPart);
      amps_malloc_managed<int>(ParticleOffsetNumber,PIC::DomainBlockDecomposition::nLocalBlocks*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_);
    

      total_number=0;

      for (int i=0;i<PIC::DomainBlockDecomposition::nLocalBlocks*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;i++) {
        ParticleOffsetNumber[i]=total_number;
        total_number+=ParticlePopulationNumberTable[i];
      }

      kernel_3<<<3,128>>>(CreateParticlePopulationTable,ParticlePopulationTable,ParticleOffsetNumber,PIC::DomainBlockDecomposition::BlockTable); 
      cudaDeviceSynchronize();





//      CreateParticlePopulationnumberTable(ParticlePopulationNumberTable,PIC::DomainBlockDecomposition::BlockTable); 

total_number=0;
for (int i=0;i<PIC::DomainBlockDecomposition::nLocalBlocks*_BLOCK_CELLS_X_*_BLOCK_CELLS_Y_*_BLOCK_CELLS_Z_;i++)total_number+=ParticlePopulationNumberTable[i];



     amps_free_managed(ParticlePopulationNumberTable);
#endif
}

    
      //PIC::Mesh::mesh->outputMeshDataTECPLOT("1.dat",0);
    
      //TransportEquation::TimeStep();
  
      PIC::TimeStep();
      //PIC::FieldSolver::Electromagnetic::ECSIM::TimeStep();

      //PIC::Mesh::mesh->outputMeshDataTECPLOT("2.dat",0);


      switch (_PIC_BC__PERIODIC_MODE_) {
      case _PIC_BC__PERIODIC_MODE_OFF_:
	PIC::Mesh::mesh->ParallelBlockDataExchange();
	break;

      case _PIC_BC__PERIODIC_MODE_ON_:
	exit(__LINE__,__FILE__,"error: the test is intended to test non-periodic BC");
	PIC::BC::ExternalBoundary::UpdateData();
	break;
      }

    }
  
  PIC::RunTimeSystemState::CumulativeTiming::Print();
  MPI_Finalize();
  
  TimeValue=time(NULL);
  ct=localtime(&TimeValue);
  
  printf("end: (%i/%i %i:%i:%i)\n",ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec);

  cout << "End of the run" << endl;
  return EXIT_SUCCESS;


}

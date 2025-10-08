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

#include "PeriodicBCTest.dfn"

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

// ---------------- CLI configuration ----------------
struct TestConfig {
  enum class Mode { WithParticles, NoParticles, FieldOnlyB, FieldOnlyE } mode = Mode::WithParticles;
  bool userB = false, userE = false;
  double B0[3] = {0.0, 0.0, 0.0};
  double E0[3] = {0.0, 0.0, 0.0};
  int stencilOrder = 2;
};

static bool TryRead3(int& i, int argc, char** argv, double v[3]) {
  if (i + 3 >= argc) return false;
  char* end=nullptr;
  for (int k=0; k<3; ++k) {
    end = nullptr;
    v[k] = std::strtod(argv[i+1+k], &end);
    if (end==argv[i+1+k] || !end) return false;
  }
  i += 3;
  return true;
}

static bool ParseIntAfterEqOrNext(int& i, int argc, char** argv, const char* /*opt*/, int& outVal) {
  std::string a(argv[i]);
  auto pos = a.find('=');
  if (pos != std::string::npos) {
    const char* s = a.c_str() + pos + 1;
    char* end = nullptr;
    long v = std::strtol(s, &end, 10);
    if (end && end != s) { outVal = static_cast<int>(v); return true; }
    return false;
  } else {
    if (i + 1 >= argc) return false;
    char* end = nullptr;
    long v = std::strtol(argv[i+1], &end, 10);
    if (end && end != argv[i+1]) { outVal = static_cast<int>(v); ++i; return true; }
    return false;
  }
}

static TestConfig ConfigureTestFromArgs(int argc, char** argv) {
  TestConfig cfg;
  for (int i=1; i<argc; ++i) {
    std::string a(argv[i]);
    if (a=="-particles") {
      cfg.mode = TestConfig::Mode::WithParticles;
    } else if (a=="-no-particles") {
      cfg.mode = TestConfig::Mode::NoParticles;
    } else if (a=="-B") {
      cfg.mode = TestConfig::Mode::FieldOnlyB;
      cfg.userB = TryRead3(i, argc, argv, cfg.B0);
      if (!cfg.userB) { cfg.B0[0]=0.0; cfg.B0[1]=1.0; cfg.B0[2]=0.0; }
    } else if (a=="-E") {
      cfg.mode = TestConfig::Mode::FieldOnlyE;
      cfg.userE = TryRead3(i, argc, argv, cfg.E0);
      if (!cfg.userE) { cfg.E0[0]=1.0; cfg.E0[1]=0.0; cfg.E0[2]=0.0; }
    } else if (a.rfind("-stencil-order",0)==0) {
      int val = cfg.stencilOrder;
      if (!ParseIntAfterEqOrNext(i, argc, argv, "-stencil-order", val)) {
        std::printf("Invalid -stencil-order value; keeping default %d\n", cfg.stencilOrder);
      } else {
        cfg.stencilOrder = val;
      }
    } else if (a=="-h" || a=="--help") {
      std::printf(
        "Usage:\n"
        "  %s [-particles | -no-particles | -B [Bx By Bz] | -E [Ex Ey Ez]] [-stencil-order=N]\n"
        "Modes (choose one; -B/-E imply -no-particles):\n"
        "  -particles           : PIC with particles; constant field BC; specular reflection.\n"
        "  -no-particles        : Field-only (no particles).\n"
        "  -B [Bx By Bz]        : Field-only with uniform B (E=0). Default: (0,1,0).\n"
        "  -E [Ex Ey Ez]        : Field-only with uniform E (B=0). Default: (1,0,0).\n"
        "Options:\n"
        "  -stencil-order=N     : FD stencil order for this test (default 2).\n",
        argv[0]);
      std::exit(0);
    } else {
      std::printf("Unknown option: %s (use -h for help)\n", argv[i]);
    }
  }
  return cfg;
}

// Global that operators may read (only if your code lacks a setter)
int g_TestStencilOrder = 2;

// Uniform setters used by SetIC()
// NOTE: These use the same buffer/offset accessors that your legacy SetIC() uses.
void SetUniformCornerE(const double E0[3]) {
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;

  for (auto* node = PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];
       node; node = node->nextNodeThisThread) {
    if (!node->block) continue;

    char *offset;
    for (int k=0; k<=_BLOCK_CELLS_Z_; ++k)
    for (int j=0; j<=_BLOCK_CELLS_Y_; ++j)
    for (int i=0; i<=_BLOCK_CELLS_X_; ++i) {
      auto* cn = node->block->GetCornerNode(PIC::Mesh::mesh->getCornerNodeLocalNumber(i,j,k));
      if (!cn) continue;
      offset = cn->GetAssociatedDataBufferPointer() + PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset;

      ((double*)(offset + CurrentEOffset))[ExOffsetIndex] = E0[0];
      ((double*)(offset + CurrentEOffset))[EyOffsetIndex] = E0[1];
      ((double*)(offset + CurrentEOffset))[EzOffsetIndex] = E0[2];

      ((double*)(offset + OffsetE_HalfTimeStep))[ExOffsetIndex] = E0[0];
      ((double*)(offset + OffsetE_HalfTimeStep))[EyOffsetIndex] = E0[1];
      ((double*)(offset + OffsetE_HalfTimeStep))[EzOffsetIndex] = E0[2];
    }
  }
}

void SetUniformCenterB(const double B0[3]) {
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;

  for (auto* node = PIC::Mesh::mesh->ParallelNodesDistributionList[PIC::ThisThread];
       node; node = node->nextNodeThisThread) {
    if (!node->block) continue;

    char *offset;
    for (int k=0; k<_BLOCK_CELLS_Z_; ++k)
    for (int j=0; j<_BLOCK_CELLS_Y_; ++j)
    for (int i=0; i<_BLOCK_CELLS_X_; ++i) {
      auto* cc = node->block->GetCenterNode(PIC::Mesh::mesh->getCenterNodeLocalNumber(i,j,k));
      if (!cc) continue;
      offset = cc->GetAssociatedDataBufferPointer() + PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;

      ((double*)(offset + CurrentBOffset))[BxOffsetIndex] = B0[0];
      ((double*)(offset + CurrentBOffset))[ByOffsetIndex] = B0[1];
      ((double*)(offset + CurrentBOffset))[BzOffsetIndex] = B0[2];

      ((double*)(offset + PrevBOffset))[BxOffsetIndex] = B0[0];
      ((double*)(offset + PrevBOffset))[ByOffsetIndex] = B0[1];
      ((double*)(offset + PrevBOffset))[BzOffsetIndex] = B0[2];
    }
  }
}

// New SetIC based on CLI
const TestConfig cfg;
void SetIC() {
  const double Z[3] = {0.0,0.0,0.0};
  // Zero fields
  SetUniformCornerE(Z);
  SetUniformCenterB(Z);

  // Apply chosen mode
  switch (cfg.mode) {
    case TestConfig::Mode::FieldOnlyB:
      SetUniformCenterB(cfg.B0);
      break;
    case TestConfig::Mode::FieldOnlyE:
      SetUniformCornerE(cfg.E0);
      break;
    case TestConfig::Mode::WithParticles:
    case TestConfig::Mode::NoParticles:
    default:
      // keep zero fields (or your legacy defaults)
      break;
  }

  PIC::Mesh::mesh->ParallelBlockDataExchange();
}

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
	 }   
       }
     }
     
  }

}


long int PrepopulateDomain() {
  using namespace PIC::FieldSolver::Electromagnetic::ECSIM;
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

  double *ParticleDataTable=NULL,*ParticleDataTable_dev=NULL;
  int ParticleDataTableIndex=0,ParticleDataTableLength=0;

  
  int nBlock[3]={_BLOCK_CELLS_X_,_BLOCK_CELLS_Y_,_BLOCK_CELLS_Z_};

  //the boundaries of the block and middle point of the cell
  double *xminBlock,*xmaxBlock;
  double v[3],anpart;
  int npart;
  char * offset=NULL;
  int ionSpec=0, electronSpec=1;
  double ionMass = PIC::MolecularData::GetMass(ionSpec)/_AMU_;
  double electronMass = PIC::MolecularData::GetMass(electronSpec)/_AMU_;

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

//    double *ParticleDataTable=NULL,*ParticleDataTable_dev=NULL;
//    int ParticleDataTableIndex=0,ParticleDataTableLength=0;

    for (kCell=0;kCell<nBlock[2];kCell++) for (jCell=0;jCell<nBlock[1];jCell++) for (iCell=0;iCell<nBlock[0];iCell++) {
	  //      nd=PIC::Mesh::mesh->getCenterNodeLocalNumber(iCell,jCell,kCell);

      // cell=cellList[nd];
      //  xMiddle=cell->GetX();
      //offset = cell->GetAssociatedDataBufferPointer()+PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset;
	  int ind[3]={iCell,jCell,kCell};
	  double x[3];
	  for (int idim=0; idim<3; idim++) x[idim]=xminBlock[idim]+(ind[idim]+0.5)*dx[idim];
	  
          
          double waveShape =sin(waveNumber[0]*(x[0]-xmin[0])+waveNumber[1]*(x[1]-xmin[1])+waveNumber[2]*(x[2]-xmin[2]));
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
          
          #if _CUDA_MODE_ == _ON_
          if (ParticleDataTableLength<npart) {
            if (ParticleDataTable!=NULL) {
              delete [] ParticleDataTable;
              cudaFree(ParticleDataTable_dev);
            }

            ParticleDataTable=new double [9*npart];
            cudaMalloc(&ParticleDataTable_dev,9*npart*sizeof(double));
          }
          
          ParticleDataTableIndex=0;
          #endif 

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
            
            #if _CUDA_MODE_ == _OFF_ 
            PIC::ParticleBuffer::InitiateParticle(xPar, electronVelocity,NULL,&electronSpec,NULL,_PIC_INIT_PARTICLE_MODE__ADD2LIST_,(void*)node);
            PIC::ParticleBuffer::InitiateParticle(xPar, ionVelocity,NULL,&ionSpec,NULL,_PIC_INIT_PARTICLE_MODE__ADD2LIST_,(void*)node);
            #else 
           
            memcpy(ParticleDataTable+0+9*ParticleDataTableIndex,xPar,3*sizeof(double));
            memcpy(ParticleDataTable+3+9*ParticleDataTableIndex,electronVelocity,3*sizeof(double));
            memcpy(ParticleDataTable+6+9*ParticleDataTableIndex,ionVelocity,3*sizeof(double));

            ParticleDataTableIndex++;
            #endif
            
          }
      //end of the particle injection block
      //std::cout<<"finished injecting npart: "<<npart<<std::endl;
      
     #if _CUDA_MODE_ == _ON_ 
          auto InitParticle = [=] _TARGET_DEVICE_ (double *ParticleDataTable, int ParticleDataTableIndex,int electronSpec, int ionSpec, void *node) {
            int id=blockIdx.x*blockDim.x+threadIdx.x;
            int increment=gridDim.x*blockDim.x;

            for (int i=id;i<ParticleDataTableIndex;i+=increment) {
              PIC::ParticleBuffer::InitiateParticle(ParticleDataTable+0+9*i,ParticleDataTable+3+9*i,NULL,&electronSpec,NULL,_PIC_INIT_PARTICLE_MODE__ADD2LIST_,(void*)node);
              PIC::ParticleBuffer::InitiateParticle(ParticleDataTable+0+9*i,ParticleDataTable+6+9*i,NULL,&ionSpec,NULL,_PIC_INIT_PARTICLE_MODE__ADD2LIST_,(void*)node);
            }
          };


          cudaMemcpy(ParticleDataTable_dev,ParticleDataTable,9*ParticleDataTableIndex*sizeof(double),cudaMemcpyHostToDevice);

          kernel_5<<<1,1>>>(InitParticle,ParticleDataTable_dev,ParticleDataTableIndex,electronSpec,ionSpec,node);
          cudaDeviceSynchronize();
      #endif
      
      
        }
        }

    #if _CUDA_MODE_ == _ON_
    delete [] ParticleDataTable;
    cudaFree(ParticleDataTable_dev);
    #endif

  MPI_Allreduce(&nLocalInjectedParticles,&nGlobalInjectedParticles,1,MPI_LONG,MPI_SUM,MPI_GLOBAL_COMMUNICATOR);
  printf("particles prepopulated!\n");
  return nGlobalInjectedParticles;
}



double localTimeStep(int spec,cTreeNodeAMR<PIC::Mesh::cDataBlockAMR> *startNode) {
    double CellSize;
    double CharacteristicSpeed;
    double dt;


    CellSize=startNode->GetCharacteristicCellSize();
    //return 0.3*CellSize/CharacteristicSpeed;

    //return 0.05;
    return 1;
}


double BulletLocalResolution(double *x) {                                                                                           
  double dist = xmax[0]-xmin[0];

#ifndef _UNIFORM_MESH_
#error ERROR: _UNIFORM_MESH_ is used but not defined
#endif
#ifndef _TEST_MESH_MODE_
#error ERROR: _TEST_MESH_MODE_ is used but not defined
#endif
#if _TEST_MESH_MODE_==_UNIFORM_MESH_  
  double res = 3;
#endif

#ifndef _NONUNIFORM_MESH_
#error ERROR: _NONUNIFORM_MESH_ is used but not defined
#endif
#ifndef _TEST_MESH_MODE_
#error ERROR: _TEST_MESH_MODE_ is used but not defined
#endif
#if _TEST_MESH_MODE_==_NONUNIFORM_MESH_
  double highRes = dist/32.0, lowRes= dist/2.0;     
  double res =(5-1)/dist*(x[0]-xmin[0])+1;  
#endif

  res=sqrt(3)+0.1;
  return res;
}
                       

int main(int argc,char **argv) {
  
   time_t TimeValue=time(NULL);
   tm *ct=localtime(&TimeValue);
  
   printf("start: (%i/%i %i:%i:%i)\n",ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec);


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
  
  PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(0,1e-2*0.0795774715459477);
  PIC::ParticleWeightTimeStep::SetGlobalParticleWeight(1,1e-2*0.0795774715459477);

  PIC::DomainBlockDecomposition::UpdateBlockTable();

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

      PrepopulateDomain();

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

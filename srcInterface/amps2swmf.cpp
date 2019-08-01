//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//$Id$
//the interface between AMPS and SWMF

#include "mpi.h"

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string.h>
#include <list>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <iostream>
#include <fstream>
#include <signal.h>
#include <sstream>
#include <string>

#include <sys/time.h>
#include <sys/resource.h>

#include "pic.h"
#include "amps2swmf.h"
#include "FluidPicInterface.h"

using namespace std;


void amps_init();
void amps_init_mesh();
void amps_time_step();


extern "C" { 
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__FLUID_
  void amps_from_gm_init(int *ParamInt, double *ParamReal, char *NameVar);
  void amps_dynamic_allocate_blocks();
  void amps_dynamic_init_blocks();
#endif
  void amps_timestep_(double* TimeSimulation, double* TimeSimulationLimit);
  int initamps_();
  void amps_impose_global_time_step_(double *Dt);
  void amps_setmpicommunicator_(signed int* iComm,signed int* iProc,signed int* nProc, signed int* nThread);
  void amps_finalize_();

  //import magnetic field from GM onto the 'center' nodes
  void amps_get_center_point_number(int*);
  void amps_get_center_point_coordinates(double*);

  //return the number of the AMPS' mesh rebalancing operations
  void amps_mesh_id_(int* id) {
    //*id=PIC::Mesh::mesh.nParallelListRedistributions;
    PIC::Mesh::mesh.SyncMeshID();
    *id=PIC::Mesh::mesh.GetMeshID();
    //*id=aa;
  }

  void amps_setmpicommunicator_(signed int* iComm,signed int* iProc,signed int* nProc, signed int* nThread) {
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
    omp_set_num_threads(*nThread);
    PIC::CPLR::SWMF::ConvertMpiCommunicatorFortran2C(iComm,iProc,nProc);
    //initialize the coupler and AMPS
    PIC::CPLR::SWMF::init();
    amps_init_mesh();
#elif _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__FLUID_
    omp_set_num_threads(*nThread);
    PIC::CPLR::FLUID::ConvertMpiCommunicatorFortran2C(iComm,iProc,nProc);
    //initialize the coupler and AMPS
    //PIC::CPLR::FLUID::init();
#endif
    // amps_init_mesh();
  }


  void amps_get_center_point_number_(int *nCenterPoints) {

#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
    PIC::CPLR::SWMF::GetCenterPointNumber(nCenterPoints);
#elif _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__FLUID_
    //PIC::CPLR::FLUID::GetCenterPointNumber(nCenterPoints); 
#endif 

  }

  void amps_get_center_point_coordinates_(double *x) {

#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
    PIC::CPLR::SWMF::GetCenterPointCoordinates(x);
#elif _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__FLUID_
    //PIC::CPLR::FLUID::GetCenterPointCoordinates(x);
#endif 

  }


  void amps_recieve_batsrus2amps_center_point_data_(char *NameVar, int *nVar, double *data,int *index) {
    
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
    PIC::CPLR::SWMF::RecieveCenterPointData(NameVar,*nVar,data,index);
#elif _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__FLUID_
    //PIC::CPLR::FLUID::RecieveCenterPointData(NameVar,*nVar,data,index);
#endif 
  }

  
  void amps_get_corner_point_number_(int *nCornerPoints) {
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__FLUID_
    PIC::CPLR::FLUID::GetCornerPointNumber(nCornerPoints);
#endif
  }
  
  void amps_get_corner_point_coordinates_(double *x) {
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__FLUID_
    PIC::CPLR::FLUID::GetCornerPointCoordinates(x);
#endif
  }
  
  
  void amps_recieve_batsrus2amps_corner_point_data_(char *NameVar, int *nVar, double *data,int *index) {
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__FLUID_    
    PIC::CPLR::FLUID::ReceiveCornerPointData(NameVar,*nVar,data,index);
#endif  
  }


  void amps_impose_global_time_step_(double *Dt) {
    for(int spec=0; spec<PIC::nTotalSpecies; spec++)
      if(PIC::ParticleWeightTimeStep::GlobalTimeStep[spec] != *Dt){
	if(PIC::ThisThread == 0)
	  std::cout<<"AMPS: global time step has been changed to "<<
	    *Dt<<" sec for species "<< spec << std::endl;
	PIC::ParticleWeightTimeStep::GlobalTimeStep[spec] = *Dt;
      }
  }

  void amps_send_batsrus2amps_center_point_data_(char *NameVar, int *nVarIn, int *nDimIn, int *nPoint, double *Xyz_DI, double *Data_VI) {
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_ 
    list<PIC::CPLR::SWMF::fSendCenterPointData>::iterator f;
#elif  _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__FLUID_ 
    list<PIC::CPLR::FLUID::fSendCenterPointData>::iterator f;
#endif

#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__FLUID_
    static bool initFlag=false;
    if (!initFlag){
      amps_init();
      initFlag = true; 
    }
#endif

#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_ 
    for (f=PIC::CPLR::SWMF::SendCenterPointData.begin();f!=PIC::CPLR::SWMF::SendCenterPointData.end();f++) {
      (*f)(NameVar,nVarIn,nDimIn,nPoint,Xyz_DI,Data_VI);
    }
#elif  _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__FLUID_ 
    for (f=PIC::CPLR::FLUID::SendCenterPointData.begin();f!=PIC::CPLR::FLUID::SendCenterPointData.end();f++) {
      (*f)(NameVar,nVarIn,nDimIn,nPoint,Xyz_DI,Data_VI);
    }
#endif

  }

  void amps_timestep_(double* TimeSimulation, double* TimeSimulationLimit) {
    static bool InitFlag=false;
    static double swmfTimeSimulation=-1.0;

    if (swmfTimeSimulation<0.0) swmfTimeSimulation=*TimeSimulation;

#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_    
    //call AMPS only after the first coupling has occured
    if (PIC::CPLR::SWMF::FirstCouplingOccured==false) {
      *TimeSimulation=*TimeSimulationLimit;
      return;
    }
#elif _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__FLUID_    
    //call AMPS only after the first coupling has occured
    if (PIC::CPLR::FLUID::FirstCouplingOccured==false) {
      *TimeSimulation=*TimeSimulationLimit;
      return;
    }
#endif

#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
    //init AMPS
    if (InitFlag==false) {
      //initamps_();

      amps_init();
      InitFlag=true;

      //print the output file on each iteration
      //PIC::RequiredSampleLength=1;
    }
#endif

    //determine whether to proceed with the current iteraction
    if (swmfTimeSimulation+PIC::ParticleWeightTimeStep::GlobalTimeStep[0]*PIC::CPLR::FLUID::FluidInterface.getNo2SiT()>*TimeSimulationLimit) {
      *TimeSimulation=*TimeSimulationLimit;
      return;
    }
    else {
      swmfTimeSimulation+=PIC::ParticleWeightTimeStep::GlobalTimeStep[0]*PIC::CPLR::FLUID::FluidInterface.getNo2SiT();
      *TimeSimulation=swmfTimeSimulation;
    }

    //call AMPS
    static long int counter=0;
    counter++;

    amps_time_step();

/*
   //GetMeanParticleMicroscopicParameters() is called ub amps_finalize_()
    if (PIC::ModelTestRun::mode==true) if (counter==PIC::ModelTestRun::nTotalIteraction) {
      char fname[400];

      sprintf(fname,"%s/amps.dat",PIC::OutputDataFileDirectory);
      PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname);

      exit(0);
    }
*/

  }

  void amps_finalize_() {
    char fname[_MAX_STRING_LENGTH_PIC_];

    //output the test run particle data
    sprintf(fname,"%s/amps.dat",PIC::OutputDataFileDirectory);
    
    #if _PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_ 
    PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname);
    #endif

    //save particle trajectory file
    #if _PIC_PARTICLE_TRACKER_MODE_ == _PIC_MODE_ON_
    sprintf(fname,"%s/amps.TrajectoryTracking.out=Final",PIC::OutputDataFileDirectory);
    PIC::ParticleTracker::OutputTrajectory(fname);
   #endif
  }


  int amps_read_param_(char *param, int *nlines, int *ncharline, int *iProc){
    // convert character array to string stream object                                                                                                      
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
    std::stringstream ss;

    AMPS2SWMF::PARAMIN::char_to_stringstream(param, *nlines, *ncharline,&ss);
    AMPS2SWMF::PARAMIN::read_paramin(&ss);
#endif

#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__FLUID_
    std::string paramString; 
    // Convert char* to a string with the function below in ReadParam.h
    char_to_string(paramString, param, (*nlines)*(*ncharline), (*ncharline));
    // Use a string to initialize readParam.
    PIC::CPLR::FLUID::FluidInterface.readParam = paramString; 
#endif

    return 0;
  }

  //find the thread number that point 'x' is located at
  void amps_get_point_thread_number_(int *thread,double *x) {
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
    static cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=NULL;

    node=PIC::Mesh::mesh.findTreeNode(x,node);
    *thread=(node!=NULL) ? node->Thread : -1;
#endif

#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__FLUID_
    static cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=NULL;    
    double pic_D[3]; 

    PIC::CPLR::FLUID::FluidInterface.mhd_to_Pic_Vec(x, pic_D);
    node=PIC::Mesh::mesh.findTreeNode(pic_D,node);
    *thread=(node!=NULL) ? node->Thread : -2;
    if (node->IsUsedInCalculationFlag==false) *thread=-1;
   
    if (*thread==-2) printf("cannot find x:%e,%e,%e\n",x[0],x[1],x[2]);
#endif
  }

  #if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__FLUID_
  void amps_from_gm_init_(int *ParamInt, double *ParamReal, char *NameVar){
    string nameFunc = "PC: amps_from_gm_init";
    std::stringstream *ss;
    ss = new std::stringstream;
    (*ss) << NameVar;

    int nPIC = 1, iPIC = 0, nParamRegion = 21; 


    PIC::CPLR::FLUID::FluidInterface.set_myrank(PIC::ThisThread);
    PIC::CPLR::FLUID::FluidInterface.set_nProcs(PIC::nTotalThreads);

    PIC::CPLR::FLUID::FluidInterface.ReadFromGMinit(ParamInt, 
						 &ParamReal[iPIC*nParamRegion], 
						 &ParamReal[nPIC*nParamRegion], 
						 ss);
    
    // The domain size and resolution is in the FluidInterface now. 
    amps_init_mesh();
    
    PIC::CPLR::FLUID::read_param();

    PIC::CPLR::FLUID::FluidInterface.PrintFluidPicInterface();
    
    //amps_init();
    delete ss;
  }    

  void amps_dynamic_allocate_blocks_(){
    if (PIC::FieldSolver::Electromagnetic::ECSIM::dynamicAllocateBlocks)
      PIC::FieldSolver::Electromagnetic::ECSIM::dynamicAllocateBlocks();

  }

  void amps_dynamic_init_blocks_(){
    if (PIC::FieldSolver::Electromagnetic::ECSIM::initNewBlocks)
    PIC::FieldSolver::Electromagnetic::ECSIM::initNewBlocks();
  }
  #endif

}

/*
int main () {
  return 1;
*/

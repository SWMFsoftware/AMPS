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


//the location of the Earth as calcualted with the SWMF. Used for heliophysics modeling  
double AMPS2SWMF::xEarthHgi[3]={0.0,0.0,0.0}; 

char AMPS2SWMF::ComponentName[10]="";
int AMPS2SWMF::ComponentID=_AMPS_SWMF_UNDEFINED_; 

//the namespace containds variables used in heliosphere simulations
double AMPS2SWMF::Heliosphere::rMin=-1.0;  

//parameters of the current SWMF session
int AMPS2SWMF::iSession=-1;
double AMPS2SWMF::swmfTimeSimulation=-1.0;
bool AMPS2SWMF::swmfTimeAccurate=true;

//amps_init_flag
bool AMPS2SWMF::amps_init_flag=false;

//amps execution timer 
PIC::Debugger::cTimer AMPS2SWMF::ExecutionTimer(_PIC_TIMER_MODE_HRES_);  

//hook that AMPS applications can use so a user-defined function is called at the end of the SWMF simulation
AMPS2SWMF::fUserFinalizeSimulation AMPS2SWMF::UserFinalizeSimulation=NULL;

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
  void amps_save_restart_();
  void amps_finalize_();

  void amps_init_session_(int* iSession,double *TimeSimulation,int* swmfTimeAccurateMode);

  //set the component name 
  void amps_set_component_name_(char* l) {
//    sprintf(AMPS2SWMF::ComponentName,"%s",l);

    for (int i=0;i<2;i++) AMPS2SWMF::ComponentName[i]=l[i];
    AMPS2SWMF::ComponentName[2]=0;


    if (strcmp("PC",AMPS2SWMF::ComponentName)==0) AMPS2SWMF::ComponentID=_AMPS_SWMF_PC_;
    else if (strcmp("PT",AMPS2SWMF::ComponentName)==0) AMPS2SWMF::ComponentID=_AMPS_SWMF_PT_;
    else exit(__LINE__,__FILE__,"Error: the component code is umnknown");  
  } 

  //determine whether a point belongs to a fomain of a particlelar SWMF component 
  bool IsDomainSC(double *x) {
    double r2=x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
    bool res=true;

    if (AMPS2SWMF::Heliosphere::rMin>0.0) if (r2<AMPS2SWMF::Heliosphere::rMin*AMPS2SWMF::Heliosphere::rMin) res=false;  

    double t=23.0*_SUN__RADIUS_;
    if (r2>t*t) res=false;

    return res;
  }


  bool IsDomainIH(double *x) {return x[0]*x[0]+x[1]*x[1]+x[2]*x[2]>=23.0*23.0*_SUN__RADIUS_*_SUN__RADIUS_;} 

  //set the location of the Earth. The function is caled by the coupler 
  void set_earth_locaton_hgi_(double *x) {for (int idim=0;idim<3;idim++) AMPS2SWMF::xEarthHgi[idim]=x[idim];}

  //init coupling with SWMF
  void amps_from_oh_init_(int *nIonFluids,int *OhCouplingCode,int *IhCouplingCode); 
  
  //get fluid number from other SWMF component
  void amps_get_fluid_number_(int * nVarIn);

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

  void amps_get_fluid_number_(int * fluidNumberIn){
    
    PIC::CPLR::SWMF::nCommunicatedIonFluids= (*fluidNumberIn-3)/5;
    //for test
    //printf("coupler fluid number is %d\n", PIC::CPLR::SWMF::nFluid );
  }

  void amps_setmpicommunicator_(signed int* iComm,signed int* iProc,signed int* nProc, signed int* nThread) {
    #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
    omp_set_num_threads(*nThread);
    #endif

    PIC::CPLR::SWMF::ConvertMpiCommunicatorFortran2C(iComm,iProc,nProc);
  }

  void amps_init_session_(int* iSession,double *TimeSimulation,int* swmfTimeAccurateMode) {
    AMPS2SWMF::iSession=*iSession;
    AMPS2SWMF::swmfTimeSimulation=*TimeSimulation;

    AMPS2SWMF::swmfTimeAccurate=(*swmfTimeAccurateMode==1) ? true : false;
  }
  
  void amps_save_restart_(){
    printf("amps_save_restart start\n");
    
    switch (AMPS2SWMF::ComponentID) {
    case _AMPS_SWMF_PC_:
      system("mkdir -p PC");
      system("mkdir -p PC/restartOUT"); 

      PIC::Restart::SamplingData::Save("PC/restartOUT/restart_field.dat");
      PIC::Restart::SaveParticleData("PC/restartOUT/restart_particle.dat");
      break;

    case _AMPS_SWMF_PT_:
      system("mkdir -p PT");
      system("mkdir -p PT/restartOUT");

      PIC::Restart::SamplingData::Save("PT/restartOUT/restart_field.dat");
      PIC::Restart::SaveParticleData("PT/restartOUT/restart_particle.dat");
      break;

    default:
      exit(__LINE__,__FILE__,"Error: the option is unlnown"); 
    }

 
    printf("amps_save_restart end\n");
  }

  void amps_get_center_point_number_(int *nCenterPoints) {

#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
    PIC::CPLR::SWMF::GetCenterPointNumber(nCenterPoints);
#elif _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__FLUID_
    //PIC::CPLR::FLUID::GetCenterPointNumber(nCenterPoints); 
#endif 

  }

  void amps_get_center_point_number_sc_(int *nCenterPoints) {
    PIC::CPLR::SWMF::GetCenterPointNumber(nCenterPoints,IsDomainSC);
  }

  void amps_get_center_point_number_ih_(int *nCenterPoints) {
    PIC::CPLR::SWMF::GetCenterPointNumber(nCenterPoints,IsDomainIH);
  }

  void amps_get_center_point_coordinates_(double *x) {

#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
    PIC::CPLR::SWMF::GetCenterPointCoordinates(x);
#elif _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__FLUID_
    //PIC::CPLR::FLUID::GetCenterPointCoordinates(x);
#endif 

  }

  void amps_get_center_point_coordinates_sc_(double *x) {
    PIC::CPLR::SWMF::GetCenterPointCoordinates(x,IsDomainSC);
  }

  void amps_get_center_point_coordinates_ih_(double *x) {
    PIC::CPLR::SWMF::GetCenterPointCoordinates(x,IsDomainIH);
  }


  void amps_recieve_batsrus2amps_center_point_data_(char *NameVar, int *nVar, double *data,int *index) {
    
#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
    PIC::CPLR::SWMF::RecieveCenterPointData(NameVar,*nVar,data,index);
#elif _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__FLUID_
    //PIC::CPLR::FLUID::RecieveCenterPointData(NameVar,*nVar,data,index);
#endif 
  }

  void amps_recieve_batsrus2amps_center_point_data_sc_(char *NameVar, int *nVar, double *data,int *index) {
    PIC::CPLR::SWMF::RecieveCenterPointData(NameVar,*nVar,data,index,IsDomainSC);
  }

  void amps_recieve_batsrus2amps_center_point_data_ih_(char *NameVar, int *nVar, double *data,int *index) {
    PIC::CPLR::SWMF::RecieveCenterPointData(NameVar,*nVar,data,index,IsDomainIH);
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

    if (AMPS2SWMF::amps_init_flag==false) {
      amps_init();
      AMPS2SWMF::amps_init_flag=true; 
    }

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
    using namespace AMPS2SWMF;

    if (swmfTimeSimulation<0.0) swmfTimeSimulation=*TimeSimulation;
    
    auto coupler_swmf = [&] () {
      bool res=true;

      //call AMPS only after the first coupling has occured
      if (PIC::CPLR::SWMF::FirstCouplingOccured==false) {
        *TimeSimulation=*TimeSimulationLimit;
        res=false;
      }
      else { 
        //init AMPS
        if (AMPS2SWMF::amps_init_flag==false) {
          amps_init();
          AMPS2SWMF::amps_init_flag=true;
        }
      
        //determine whether to proceed with the current iteraction
        if (swmfTimeSimulation+PIC::ParticleWeightTimeStep::GlobalTimeStep[0]>*TimeSimulationLimit) {
          *TimeSimulation=*TimeSimulationLimit;
          res=false;
        }
        else {
          swmfTimeSimulation+=PIC::ParticleWeightTimeStep::GlobalTimeStep[0];
          *TimeSimulation=swmfTimeSimulation;
          res=true;
        }
      }

      return res;
    };
    
    auto coupler_fluid = [&] () {
      bool res=true;

      //call AMPS only after the first coupling has occured
      if (PIC::CPLR::FLUID::FirstCouplingOccured==false) {
        *TimeSimulation=*TimeSimulationLimit;
        res=false;
      }
      else { 
        //determine whether to proceed with the current iteraction
        if (swmfTimeSimulation+PIC::ParticleWeightTimeStep::GlobalTimeStep[0]*PIC::CPLR::FLUID::FluidInterface.getNo2SiT()>*TimeSimulationLimit) {
          *TimeSimulation=*TimeSimulationLimit;
          res=false;
        }
        else {
          swmfTimeSimulation+=PIC::ParticleWeightTimeStep::GlobalTimeStep[0]*PIC::CPLR::FLUID::FluidInterface.getNo2SiT();
          *TimeSimulation=swmfTimeSimulation;
          res=true;
        }
      }

      return res;
    };

    bool call_amps_flag=true;

do {

    switch (_PIC_COUPLER_MODE_) {
    case _PIC_COUPLER_MODE__SWMF_:
      call_amps_flag=coupler_swmf();
      break;
    case _PIC_COUPLER_MODE__FLUID_:
      call_amps_flag=coupler_fluid();
      break;
    }

    //call AMPS
    static long int counter=0;
    counter++;

    if (call_amps_flag==true) {
      ExecutionTimer.Start();
      amps_time_step();
      ExecutionTimer.UpdateTimer();

      PIC::Restart::LoadRestartSWMF=false; //in case the AMPS was set to read a restart file  
    }
}
while (false); // ((swmfTimeAccurate==true)&&(call_amps_flag==true));


  }

  void  amps_from_oh_init_(int *nIonFluids,int *OhCouplingCode,int *IhCouplingCode) {
    static bool init_flag=false;

    if (init_flag==true) return;

    init_flag=true;

    //the the coupling flags 
    if (*OhCouplingCode==1) PIC::CPLR::SWMF::OhCouplingFlag=true;
    if (*IhCouplingCode==1) PIC::CPLR::SWMF::IhCouplingFlag=true; 

    //set the total number of the ion fluids
    PIC::CPLR::SWMF::nCommunicatedIonFluids=*nIonFluids; 

    if (PIC::CPLR::SWMF::nCommunicatedIonFluids<=0) exit(__LINE__,__FILE__,"Error: the number of communicated ion fluids has to be positive");

    //initialize the coupler and AMPS
    PIC::CPLR::SWMF::init();
    amps_init_mesh();
  }

  void amps_finalize_() {
    char fname[_MAX_STRING_LENGTH_PIC_];

    //print the executed time 
     AMPS2SWMF::ExecutionTimer.PrintMeanMPI("AMPS execution time avaraged over all MPI processes involved in the simulation");

    //output the test run particle data
    sprintf(fname,"%s/amps.dat",PIC::OutputDataFileDirectory);
    
    #if _PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_ 
    PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname);
    #endif

    //call a user-defined function to finalize the application
    if (AMPS2SWMF::UserFinalizeSimulation!=NULL) AMPS2SWMF::UserFinalizeSimulation();

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

    if (node){
      if(node->IsUsedInCalculationFlag==false) *thread=-1;
    }
  
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

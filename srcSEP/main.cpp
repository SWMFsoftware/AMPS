

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
#include <iostream>
#include <fstream>
#include <time.h>

#include <sys/time.h>
#include <sys/resource.h>

//$Id$


//#include "vt_user.h"
//#include <VT.h>

//the particle class
#include "constants.h"
#include "sep.h"

#include "tests.h"

void amps_init();
void amps_init_mesh();
void amps_time_step();


int main(int argc,char **argv) {
  //      MPI_Init(&argc,&argv);


  //read post-compile input file  
  if (PIC::PostCompileInputFileName!="") {
     SEP::Parser::ReadFile(PIC::PostCompileInputFileName);
  }
  

  //setup datum to store the segment's data for the Alfven turbulence model 
  PIC::FieldLine::cFieldLineSegment::AddDatumStored(&SEP::AlfvenTurbulence_Kolmogorov::WaveEnergyDensity); 
  PIC::FieldLine::cFieldLineSegment::AddDatumStored(&SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::S);
  PIC::FieldLine::cFieldLineSegment::AddDatumStored(&SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::S_pm);
	  
 
  amps_init_mesh();
  amps_init();

  //init the Alfven turbulence IC
  SEP::AlfvenTurbulence_Kolmogorov::ModelInit::Init(); 

  if (_PIC_FIELD_LINE_MODE_==_PIC_MODE_ON_) { 
    TestManager();
  }

  int TotalIterations=(_PIC_NIGHTLY_TEST_MODE_==_PIC_MODE_ON_) ? PIC::RequiredSampleLength+10 : 100000001;  

  SEP::Diffusion::GetPitchAngleDiffusionCoefficient=SEP::Diffusion::Jokopii1966AJ::GetPitchAngleDiffusionCoefficient;


  //init MPI exchange of field line segent's associated data 
  PIC::FieldLine::Parallel::InitializeDatumStoredAtEdgeMPI();

  //time step
  for (long int niter=0;niter<TotalIterations;niter++) {
    //SEP::InitDriftVelData();
    amps_time_step();

    //reduce S
    PIC::FieldLine::Parallel::MPIAllReduceDatumStoredAtEdge(SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::S);

    //calculate S+/-


    //scatter S+/-  
    
    
    //update Alfven turbulence energy density due to particle interaction  


    //scatter Alfven turbulence energy density 


    //convect Alfven turbument energy density 



    //PIC::ParticleSplitting::Split::SplitWithVelocityShift_FL(10,200);
    //
    //PIC::ParticleSplitting::FledLine::WeightedParticleMerging(20,20,20,500,800); 
    //PIC::ParticleSplitting::FledLine::WeightedParticleSplitting(20,20,20,500,800);    
  }


  char fname[400];

  sprintf(fname,"%s/test_SEP.dat",PIC::OutputDataFileDirectory);
  PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname);

  cout << "End of the run:" << PIC::nTotalSpecies << endl;

  MPI_Finalize();
  return EXIT_SUCCESS;
}


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


//#include <saturn.h>

//#define _MAIN_PROFILE_



//#include "vt_user.h"
//#include <VT.h>

//the particle class
#include "pic.h"
#include "constants.h"
#include "Europa.h"


void amps_init();
void amps_init_mesh();
void amps_time_step();





int main(int argc,char **argv) {
//	MPI_Init(&argc,&argv);

  amps_init_mesh();
  amps_init();

  //init the location of the Sun and Jupiter in respect to Europa
  Europa::OrbitalMotion::UpdateSunJupiterLocation();


PIC::RequiredSampleLength=3000;
//PIC::RequiredSampleLength=500;



//time step
int nTotalIterations=(_PIC_NIGHTLY_TEST_MODE_==_PIC_MODE_OFF_) ? 1000000 : 150;
//int nTotalIterations=(_PIC_NIGHTLY_TEST_MODE_==_PIC_MODE_OFF_) ? 100 : 150;

	for (long int niter=0;niter<nTotalIterations;niter++) {

	  amps_time_step();

    if (PIC::Mesh::mesh->ThisThread==0) {
      time_t TimeValue=time(NULL);
      tm *ct=localtime(&TimeValue);

      printf(": (%i/%i %i:%i:%i), Iteration: %ld  (currect sample length:%ld, %ld interations to the next output)\n",ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec,niter,PIC::RequiredSampleLength,PIC::RequiredSampleLength-PIC::CollectingSampleCounter);
    }

	}


  //output the particle statistics of the test run
  #if _PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_
	char fname[300];

  sprintf(fname,"%s/test_Europa.dat",PIC::OutputDataFileDirectory);
  PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname);
  #endif


  //finish the run
  MPI_Finalize();
  cout << "End of the run:" << PIC::nTotalSpecies << endl;

  return EXIT_SUCCESS;
}

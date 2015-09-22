
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


void amps_init();
void amps_time_step();

char Exosphere::IAU_FRAME[_MAX_STRING_LENGTH_PIC_]="IAU_MERCURY";



int main(int argc,char **argv) {
//	MPI_Init(&argc,&argv);

  amps_init();




	//time step
  //determine the total number of the iterations to perform
  //in the test-mode run 100 iterations and than output the particle data statistics
  int nIterations,nTotalIterations=100000001;

  if (_PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_) nTotalIterations=100;

	for (long int niter=0;niter<nTotalIterations;niter++) {

    //print the iteration number
    if (PIC::Mesh::mesh.ThisThread==0) {
      time_t TimeValue=time(NULL);
      tm *ct=localtime(&TimeValue);

      printf(": (%i/%i %i:%i:%i), Iteration: %ld  (currect sample length:%ld, %ld interations to the next output)\n",ct->tm_mon+1,ct->tm_mday,ct->tm_hour,ct->tm_min,ct->tm_sec,niter,PIC::RequiredSampleLength,PIC::RequiredSampleLength-PIC::CollectingSampleCounter);
    }

	  amps_time_step();


    //check whether particle tracing if finished
    static int LastDataOutputFileNumber=-1;

    if ((PIC::DataOutputFileNumber!=0)&&(PIC::DataOutputFileNumber!=LastDataOutputFileNumber)) {
     int spec,TrajectoryAccumulationFinished=true;

      for (spec=0;spec<PIC::nTotalSpecies;spec++) if (PIC::ParticleTracker::maxSampledTrajectoryNumber!=0) {
        if (PIC::ParticleTracker::totalSampledTrajectoryNumber[spec]<PIC::ParticleTracker::maxSampledTrajectoryNumber) {
          TrajectoryAccumulationFinished=false;
          break;
        }
      }

      MPI_Bcast(&TrajectoryAccumulationFinished,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);

      if (TrajectoryAccumulationFinished==true) {
        //the requested number of the trajectories is reached
        //1. stop injection of new particles


        //2. check whether any particle left in the system
        int thread,flag,ParticlePresent[PIC::nTotalThreads];

        flag=(PIC::ParticleBuffer::NAllPart!=0) ? true : false;
        MPI_Allgather(&flag,1,MPI_INT,ParticlePresent,1,MPI_INT,MPI_GLOBAL_COMMUNICATOR);
        for (thread=0;thread<PIC::nTotalThreads;thread++) if (ParticlePresent[thread]==true) flag=true;

        if (flag==false) {
          //there is no particles in the system any more -> the particle tracing procedure is finished
          break;
        }
      }

      LastDataOutputFileNumber=PIC::DataOutputFileNumber;
    }
	}

  //output the particle statistics for the nightly tests
  if (_PIC_NIGHTLY_TEST_MODE_ == _PIC_MODE_ON_) {
    char fname[400];

    sprintf(fname,"%s/test_CCMC.dat",PIC::OutputDataFileDirectory);
    PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname);
  }

	if (PIC::ThisThread==0) cout << "End of the run:" << PIC::nTotalSpecies << endl;
	MPI_Finalize();

	return EXIT_SUCCESS;
}

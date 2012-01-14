//=================================================================
//$Id$
//=================================================================
//define the random generator

#include "mpi.h"

#include "rnd.h"

int RandomNumberGenerator::rndLastSeed=0;

void rnd_seed(int seed) {
  int thread;
  MPI_Comm_rank(MPI_COMM_WORLD,&thread);

  if (seed==-1) seed=thread;

  RandomNumberGenerator::rndLastSeed=seed;
}

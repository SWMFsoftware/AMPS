//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//=================================================================
//$Id$
//=================================================================
//define the random generator

#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#ifndef _RND_
#define _RND_

#ifndef _DO_NOT_LOAD_GLOBAL_H_
#include "global.h"
#endif 

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
#include <omp.h>
#endif //_PIC_COMPILATION_MODE_ == _PIC_COMPILATION_MODE__HYBRID_



namespace RandomNumberGenerator {
  extern unsigned long int rndLastSeed;
  extern unsigned long int *rndLastSeedArray;

}


struct cRndSeedContainer {
  unsigned long int Seed;
};

void rnd_seed(int seed=-1);

inline double rnd(cRndSeedContainer *SeedIn) {
  double res;
  unsigned long int Seed=SeedIn->Seed;

  //mask part of the value of the Seed to coarsen the set of the generated random numbers
  //the practical mask values should be 0xf -> mask lower 4 birs, 0xff -> lover 8 bits, 0xfff -> lower 12 bits, etc
  const int RndCoarseningBitMask=0;

  Seed*=48828125;
  Seed&=2147483647; // pow(2,31) - 1
  if (Seed==0) Seed=1;

  res=double((Seed&(~RndCoarseningBitMask))/2147483648.0); //(pow(2,31) - 1) + 1
  SeedIn->Seed=Seed;

  return res;
}

inline double rnd() {
  double res;
  cRndSeedContainer SeedContainer;

  #if _COMPILATION_MODE_ == _COMPILATION_MODE__MPI_
  SeedContainer.Seed=RandomNumberGenerator::rndLastSeed;
  res=rnd(&SeedContainer);
  RandomNumberGenerator::rndLastSeed=SeedContainer.Seed;
  #elif _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
  int thread=omp_get_thread_num();

  SeedContainer.Seed=RandomNumberGenerator::rndLastSeedArray[thread];
  res=rnd(&SeedContainer);
  RandomNumberGenerator::rndLastSeedArray[thread]=SeedContainer.Seed;
  #else
  #error Unknown option
  #endif //_COMPILATION_MODE_

  return res;
}

#endif

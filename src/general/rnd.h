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

#include "global.h"

namespace RandomNumberGenerator {
  extern unsigned long int rndLastSeed;
}

void rnd_seed(int seed=-1);

inline double rnd() {
  RandomNumberGenerator::rndLastSeed*=48828125;
  RandomNumberGenerator::rndLastSeed&=2147483647; // pow(2,31) - 1
  if (RandomNumberGenerator::rndLastSeed==0) RandomNumberGenerator::rndLastSeed=1;

  return double(RandomNumberGenerator::rndLastSeed/2147483648.0); //(pow(2,31) - 1) + 1
}


#endif

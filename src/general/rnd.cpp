//  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
//  For more information, see http://csem.engin.umich.edu/tools/swmf
//=================================================================
//$Id$
//=================================================================
//define the random generator

#include "mpi.h"

#include "rnd.h"
#include "specfunc.h"

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
#include <omp.h>
#endif //_PIC_COMPILATION_MODE_ == _PIC_COMPILATION_MODE__HYBRID_

_CUDA_MANAGED_ _TARGET_DEVICE_ unsigned long int RandomNumberGenerator::rndLastSeed=0;
_CUDA_MANAGED_ _TARGET_DEVICE_ unsigned long int *RandomNumberGenerator::rndLastSeedArray=NULL;
thread_local std::mt19937 RandomNumberGenerator::gen(42);
thread_local std::uniform_real_distribution<double> RandomNumberGenerator::dist(0.0, 1.0);
thread_local bool RandomNumberGenerator::NormalSpareAvailable=false;
thread_local double RandomNumberGenerator::NormalSpare=0.0;

#if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_ && \
    _CUDA_MODE_ != _ON_
namespace {
// Host-only capacity for the OpenMP seed table. CUDA builds use the managed
// allocation made by rnd_seedGPU() and therefore must not delete it here.
int RndHostSeedArrayCapacity=0;
}
#endif

#if _CUDA_MODE_ == _ON_
void rnd_seedGPU(int seed) {
  int i,nThreads=_CUDA_BLOCKS_*_CUDA_THREADS_;

  cudaMallocManaged(&RandomNumberGenerator::rndLastSeedArray,nThreads*sizeof(unsigned long int));
  for (i=0;i<nThreads;i++) RandomNumberGenerator::rndLastSeedArray[i]=abs(seed)+i;
}
#endif

void rnd_seed(int seed) {
  int thread;
  MPI_Comm_rank(MPI_GLOBAL_COMMUNICATOR,&thread);

  if (seed==-1) seed=1+thread;

  #if _CUDA_MODE_ == _ON_
  rnd_seedGPU(seed);
  #endif

  // In a hybrid build every worker owns both a uniform seed and thread_local
  // derived-distribution state.  Reset them inside one OpenMP team so an odd
  // number of pre-reseed Gaussian draws on a worker cannot leak a cached
  // Box--Muller variate into the restarted stream.  The same rule also seeds
  // every worker's thread_local Mersenne engine instead of only the caller's.
  #if _COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_
    int nThreadsOpenMP=1;
    #pragma omp parallel shared(nThreadsOpenMP)
    {
      #pragma omp single
      {
        nThreadsOpenMP=omp_get_num_threads();
#if _CUDA_MODE_ != _ON_
        if (RandomNumberGenerator::rndLastSeedArray==NULL ||
            RndHostSeedArrayCapacity<nThreadsOpenMP) {
          delete [] RandomNumberGenerator::rndLastSeedArray;
          RandomNumberGenerator::rndLastSeedArray=
              new unsigned long int[nThreadsOpenMP];
          RndHostSeedArrayCapacity=nThreadsOpenMP;
        }
#endif
      }
      #pragma omp barrier

      const int ompThread=omp_get_thread_num();
      const unsigned long int magnitude=(seed<0) ?
          0UL-static_cast<unsigned long int>(seed) :
          static_cast<unsigned long int>(seed);
      const unsigned long int streamSeed=(ompThread==0) ?
          static_cast<unsigned long int>(seed) :
          magnitude+static_cast<unsigned long int>(ompThread)+
          static_cast<unsigned long int>(thread)*
              static_cast<unsigned long int>(nThreadsOpenMP);

      RandomNumberGenerator::rndLastSeedArray[ompThread]=streamSeed;
      RandomNumberGenerator::NormalSpareAvailable=false;
      RandomNumberGenerator::NormalSpare=0.0;

      switch (_RND_MODE_) {
      case _RND_MODE_DEFAULT_:
        break;
      case _RND_MODE_MERSENNE_TWISTER_:
        RandomNumberGenerator::gen.seed(streamSeed);
        break;
      default:
        exit(__LINE__,__FILE__,"Error: unknown option");
      }
    }
  #else
    RandomNumberGenerator::NormalSpareAvailable=false;
    RandomNumberGenerator::NormalSpare=0.0;

    switch (_RND_MODE_) {
    case _RND_MODE_DEFAULT_:
      RandomNumberGenerator::rndLastSeed=seed;
      break;
    case _RND_MODE_MERSENNE_TWISTER_:
      RandomNumberGenerator::gen.seed(seed);
      break;
    default:
      exit(__LINE__,__FILE__,"Error: unknown option");
    }
  #endif //_COMPILATION_MODE_ == _COMPILATION_MODE__HYBRID_

}

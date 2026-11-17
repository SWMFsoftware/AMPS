//$Id$
//Deterministic work partition of the OpenMP parallel loops over the blocks (cells) of an MPI process.
//
//The particle mover and the ECSIM UpdateJMassMatrix() accumulate the contributions of the particles into shared
//lists and corner nodes. With 'schedule(dynamic)' or 'schedule(guided)' the order of these contributions depends on
//the thread timing, and so do the floating-point sums (run-to-run differences of the hybrid MPI+OpenMP builds). With
//the partition below each thread processes one contiguous range of the work items, the ranges follow the thread
//order, and every thread processes its range in increasing order: the contributions can then be merged in the same
//order as the MPI-only code processes them (see PIC::Mover::MoveParticles() and ECSIM::UpdateJMassMatrix()).
//
//The ranges are balanced by the given weights (e.g. the number of particles in a block). The weights must not
//depend on the thread timing, so that the partition, and with it everything that depends on the thread that
//processes a work item (e.g. the per-thread random number streams), is the same in every run.

#ifndef _PIC_OPENMP_PARTITION_H_
#define _PIC_OPENMP_PARTITION_H_

#include <vector>

namespace PIC {
namespace OpenMP {

  //split the work items [0,weight.size()) into nThreads contiguous ranges: thread t gets [begin[t],begin[t+1]),
  //begin[0]=0, begin[nThreads]=weight.size(). The item i goes to the thread t for which the middle of the item's
  //weight interval, (W(i)+W(i+1))/2 with W(i)=weight[0]+...+weight[i-1], falls into [t*W/nThreads,(t+1)*W/nThreads).
  //Integer arithmetic only, so the result does not depend on the compiler's floating-point code.
  inline void ContiguousPartition(const std::vector<long int>& weight,int nThreads,std::vector<int>& begin) {
    const int n=(int)weight.size();
    long int total=0;

    if (nThreads<1) nThreads=1;
    begin.assign(nThreads+1,n);
    begin[0]=0;

    for (int i=0;i<n;i++) total+=(weight[i]>0) ? weight[i] : 0;

    if (total==0) {
      //no weight: equal numbers of items
      for (int t=1;t<nThreads;t++) begin[t]=(int)(((long int)n*t)/nThreads);
      return;
    }

    long int w=0;
    int t=0;

    for (int i=0;i<n;i++) {
      long int wi=(weight[i]>0) ? weight[i] : 0;

      //2*nThreads*(W(i)+wi/2) compared with 2*total*(t+1): the item goes to the first thread whose upper limit is above
      //the middle of the item
      while ((t<nThreads-1)&&(nThreads*(2*w+wi)>=2*total*(t+1))) begin[++t]=i;
      w+=wi;
    }

    while (t<nThreads-1) begin[++t]=n;
  }

}
}

#endif

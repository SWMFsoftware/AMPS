// Compile/link probe: which spellings of fmod have a device implementation in nvc++?
// (Build only.  Do not run the binary: it would initialize a GPU.)
//
//   source /data/vtenishe/openacc/env.sh
//   nvc++ -std=c++17 -O2 -acc=gpu -gpu=cc86,cuda12.4       fmod_probe.cpp -o /dev/null
//   nvc++ -std=c++17 -O2 -acc=gpu -gpu=cc86,cuda12.4,nordc fmod_probe.cpp -o /dev/null
//   nvc++ ... '-DFN(a,b)=std::remainder(a,b)' fmod_probe.cpp -o /dev/null
//
// nvc++ 24.5 results (2026-10-21):
//   std::fmod, ::fmod, __builtin_fmod, fmodf (float)
//       rdc  : nvlink error   : Undefined reference to 'fmod' (resp. 'fmodf')
//       nordc: ptxas fatal    : Unresolved extern function 'fmod' (resp. 'fmodf')
//   std::remainder, std::remquo, std::trunc, std::modf : compile and link
#include <cmath>
#include <cstdio>
#ifndef FN
#define FN(a,b) std::fmod(a,b)
#endif
int main() {
  const int n=1024;
  double* x=new double[n];
  for (int i=0;i<n;++i) x[i]=i*0.37;
#pragma acc parallel loop copy(x[0:n])
  for (int i=0;i<n;++i) x[i]=FN(x[i],6.283185307179586);
  double s=0;
  for (int i=0;i<n;++i) s+=x[i];
  std::printf("%.17g\n",s);
  delete[] x;
}

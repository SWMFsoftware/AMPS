#ifndef AMPS_COMPLETE_STENCIL_H
#define AMPS_COMPLETE_STENCIL_H

// Append an already constructed eight-point pointer stencil without repeating
// generic AddCell checks at every insertion. This is a guarded fast path, not
// a new interpolant: the caller supplies the SAME node order and floating-point
// weights as the general path. No weight normalization or arithmetic regrouping
// occurs here. A rejected candidate leaves the destination unchanged, allowing
// the full AMR/missing-node/domain-clipping algorithm to handle it normally.
#if defined(__CUDACC__)
#define AMPS_STENCIL_HD __host__ __device__
#else
#define AMPS_STENCIL_HD
#endif
namespace amps_complete_stencil {
template<class Node,class Stencil>
AMPS_STENCIL_HD inline bool Assign(Node* const nodes[8],const int ids[8],
    const double weights[8],bool periodic,const double* lower,const double* upper,
    int dimensions,Stencil& target) {
  static_assert(sizeof(target.Weight)/sizeof(target.Weight[0])>=8,"eight-point stencil capacity required");
  for (int i=0;i<8;i++) {
    if (nodes[i]==nullptr || weights[i]<0.0) return false;
    if (!periodic) for (int d=0;d<dimensions;d++) {
      const double x=nodes[i]->GetX()[d];
      if (x<lower[d] || x>upper[d]) return false;
    }
  }
  for (int i=0;i<8;i++) {
    target.Weight[i]=weights[i]; target.cell[i]=nodes[i]; target.LocalCellID[i]=ids[i];
  }
  target.Length=8;
  return true;
}
}
#undef AMPS_STENCIL_HD
#endif

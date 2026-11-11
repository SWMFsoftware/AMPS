#ifndef AMPS_ORDERED_SPARSE_MATVEC_H
#define AMPS_ORDERED_SPARSE_MATVEC_H

// Compact, nonowning-RHS view of a fixed sparse operator. Coefficients belong
// to one explicitly bounded solve snapshot; unknown-vector VALUES remain live
// through pointers supplied by the caller. Build again after coefficients,
// row topology or RHS storage change. No matrix, MPI or PIC dependency exists.
// Row evaluation retains insertion order and scalar += operations, rather
// than reassociating a floating-point reduction or enabling a relaxed solver.
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace amps_ordered_sparse {
struct Entry { double coefficient; const double* rhs; };
class View {
  std::vector<std::size_t> offsets_;
  std::vector<double> coefficients_;
  std::vector<const double*> rhs_;
public:
  template<class Count,class Element>
  void Build(std::size_t rows,Count count,Element element) {
    offsets_.resize(rows+1); offsets_[0]=0;
    for (std::size_t row=0;row<rows;row++) {
      const auto size=count(row);
      if (size<0 || static_cast<std::size_t>(size)>std::numeric_limits<std::size_t>::max()-offsets_[row])
        throw std::invalid_argument("invalid compact sparse row length");
      offsets_[row+1]=offsets_[row]+static_cast<std::size_t>(size);
    }
    coefficients_.resize(offsets_.back()); rhs_.resize(offsets_.back());
    for (std::size_t row=0;row<rows;row++)
      for (std::size_t i=offsets_[row];i<offsets_[row+1];i++) {
        const Entry entry=element(row,i-offsets_[row]);
        if (entry.rhs==nullptr) throw std::invalid_argument("null compact sparse RHS");
        coefficients_[i]=entry.coefficient; rhs_[i]=entry.rhs;
      }
  }
  double Row(std::size_t row) const {
    if (offsets_.empty() || row>=offsets_.size()-1) throw std::out_of_range("compact sparse row");
    double result=0.0;
    for (std::size_t i=offsets_[row];i<offsets_[row+1];i++)
      result+=coefficients_[i]*(*rhs_[i]);
    return result;
  }
};
}
#endif

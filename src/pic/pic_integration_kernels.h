#ifndef AMPS_PIC_INTEGRATION_KERNELS_H
#define AMPS_PIC_INTEGRATION_KERNELS_H

/* Small deterministic kernels used by column integration and tabulated
 * stopping-power providers.  Keeping the numerical contract independent of
 * mesh/MPI state permits exact constant/affine and knot/limit tests. */

#include <cmath>
#include <cstddef>

namespace PIC {
namespace IntegrationKernels {

inline bool TrapezoidAccumulate(double* integral, const double* left,
                                const double* right, std::size_t count,
                                double path_increment) {
  if (integral==0 || left==0 || right==0 ||
      !std::isfinite(path_increment) || path_increment<0.0) return false;
  for (std::size_t i=0;i<count;++i) {
    if (!std::isfinite(left[i]) || !std::isfinite(right[i])) return false;
    integral[i]+=0.5*(left[i]+right[i])*path_increment;
  }
  return true;
}

enum class TableStatus {
  Success,
  BelowValidityRange,
  AboveValidityRange,
  InvalidTable,
  InvalidArgument
};

/* Linear interpolation without silent extrapolation.  Abscissae must be
 * strictly increasing and all values finite.  End knots are accepted exactly;
 * values outside the validity range return a typed status and leave output
 * untouched.  Physical stopping-power tables may wrap this kernel with a
 * documented log transform, but their validity policy remains explicit. */
inline TableStatus InterpolateTable(const double* x, const double* y,
                                    std::size_t count, double query,
                                    double& value) {
  if (x==0 || y==0 || count<2 || !std::isfinite(query))
    return TableStatus::InvalidArgument;
  for (std::size_t i=0;i<count;++i) {
    if (!std::isfinite(x[i]) || !std::isfinite(y[i]) ||
        (i!=0 && !(x[i]>x[i-1]))) return TableStatus::InvalidTable;
  }
  if (query<x[0]) return TableStatus::BelowValidityRange;
  if (query>x[count-1]) return TableStatus::AboveValidityRange;
  if (query==x[count-1]) { value=y[count-1]; return TableStatus::Success; }
  std::size_t low=0,high=count-1;
  while (high-low>1) {
    const std::size_t middle=low+(high-low)/2;
    if (query<x[middle]) high=middle;
    else low=middle;
  }
  const double fraction=(query-x[low])/(x[high]-x[low]);
  value=y[low]+fraction*(y[high]-y[low]);
  return TableStatus::Success;
}

}  // namespace IntegrationKernels
}  // namespace PIC

#endif

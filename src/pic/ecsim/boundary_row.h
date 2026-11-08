#ifndef AMPS_PIC_ECSIM_BOUNDARY_ROW_H
#define AMPS_PIC_ECSIM_BOUNDARY_ROW_H

/* Dependency-light ECSIM increment-row policy.
 *
 * The solver unknown is delta-E.  A fixed-field (Dirichlet) face therefore
 * uses delta-E_b=0, while homogeneous Neumann uses
 * delta-E_b - average(delta-E_inward)=0 at edges/corners.  Periodic faces do
 * not create a physical constraint row.  A perfect-conductor vector closure
 * requires component/normal information absent from the current scalar face
 * policy and is rejected explicitly instead of being treated as Dirichlet.
 */

#include <cmath>

namespace PIC { namespace FieldSolver { namespace Electromagnetic {
namespace ECSIM { namespace BoundaryRow {

enum class Status { InteriorOrPeriodic, Dirichlet, Neumann, Unsupported, Invalid };

struct Policy {
  Status status;
  int entries;
  double diagonal;
  double donor_weight;
};

inline Policy Resolve(bool periodic, int bc_type, int donor_count) {
  if (periodic) return {Status::InteriorOrPeriodic,0,0.0,0.0};
  if (bc_type==1) return {Status::Dirichlet,1,1.0,0.0};
  if (bc_type==0) {
    if (donor_count<1 || donor_count>3)
      return {Status::Invalid,0,0.0,0.0};
    return {Status::Neumann,1+donor_count,1.0,
            -1.0/static_cast<double>(donor_count)};
  }
  // bc_type==2 is reserved for a future component-aware perfect-conductor
  // closure; all other nonnegative codes are likewise unsupported.
  if (bc_type>=2) return {Status::Unsupported,0,0.0,0.0};
  return {Status::InteriorOrPeriodic,0,0.0,0.0};
}

}}}}}  // namespace PIC::FieldSolver::Electromagnetic::ECSIM::BoundaryRow

#endif

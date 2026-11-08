#ifndef _AMPS_PIC_FIELD_SOLVER_KERNELS_H_
#define _AMPS_PIC_FIELD_SOLVER_KERNELS_H_

/*
 * Small, dependency-light reference kernels used by PIC field modules and
 * their qualification tests.  SI factors are intentionally explicit at the
 * API boundary; the discrete equations below are otherwise nondimensional.
 *
 * SolvePoissonDirichlet1D solves -d2(phi)/dx2=rhs on [0,L], phi=0 at both
 * ends, with the symmetric second-order three-point matrix and unpreconditioned
 * conjugate gradients.  ResidualHistory stores ||b-Ax||_2, including the
 * initial residual.  Breakdown and nonfinite input are typed failures.
 *
 * AdvanceVacuumYee1D advances the periodic one-dimensional Maxwell pair
 *   B^(n+1/2)=B^(n-1/2)-c dt D+ E^n,
 *   E^(n+1)=E^n-c dt D- B^(n+1/2).
 * E lives at integer nodes and B at half nodes.  Stability requires
 * |c dt/dx|<=1.  For a Fourier mode, the independently testable dispersion is
 * sin(omega dt/2)=C sin(k dx/2).  The staggered invariant is
 * 1/2 sum(E_n^2 + B_(n+1/2) B_(n-1/2)); it is not the instantaneous SI field
 * energy, but converges to it as dt->0 and is exactly conserved by the update.
 */

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace PIC { namespace FieldSolverKernels {

enum class Status { Success,InvalidArgument,NoConvergence,NumericalBreakdown };

struct IterativeResult {
  Status status=Status::InvalidArgument;
  int iterations=0;
  std::vector<double> residual_history;
};

inline double Dot(const std::vector<double>& a,const std::vector<double>& b) {
  double value=0.0;
  for (std::size_t i=0;i<a.size();++i) value+=a[i]*b[i];
  return value;
}

inline IterativeResult SolvePoissonDirichlet1D(
    const std::vector<double>& rhs,double dx,double relative_tolerance,
    int maximum_iterations,std::vector<double>& solution) {
  IterativeResult result;
  if (rhs.empty() || !std::isfinite(dx) || dx<=0.0 ||
      !std::isfinite(relative_tolerance) || relative_tolerance<=0.0 ||
      maximum_iterations<=0) return result;
  for (double value:rhs) if (!std::isfinite(value)) return result;

  const std::size_t n=rhs.size();
  solution.assign(n,0.0);
  std::vector<double> residual=rhs,direction=rhs,product(n,0.0);
  const double inverse_dx2=1.0/(dx*dx);
  auto multiply=[&](const std::vector<double>& in,std::vector<double>& out) {
    for (std::size_t i=0;i<n;++i) {
      const double left=i==0 ? 0.0 : in[i-1];
      const double right=i+1==n ? 0.0 : in[i+1];
      out[i]=(2.0*in[i]-left-right)*inverse_dx2;
    }
  };
  double rr=Dot(residual,residual);
  const double initial=std::sqrt(rr);
  result.residual_history.push_back(initial);
  if (initial==0.0) { result.status=Status::Success; return result; }
  const double threshold=relative_tolerance*initial;

  for (int iteration=1;iteration<=maximum_iterations;++iteration) {
    multiply(direction,product);
    const double denominator=Dot(direction,product);
    if (!std::isfinite(denominator) || denominator<=0.0) {
      result.status=Status::NumericalBreakdown; return result;
    }
    const double alpha=rr/denominator;
    for (std::size_t i=0;i<n;++i) {
      solution[i]+=alpha*direction[i];
      residual[i]-=alpha*product[i];
    }
    const double rr_new=Dot(residual,residual);
    const double norm=std::sqrt(rr_new);
    result.residual_history.push_back(norm);
    result.iterations=iteration;
    if (!std::isfinite(norm)) {
      result.status=Status::NumericalBreakdown; return result;
    }
    if (norm<=threshold) { result.status=Status::Success; return result; }
    const double beta=rr_new/rr;
    for (std::size_t i=0;i<n;++i) direction[i]=residual[i]+beta*direction[i];
    rr=rr_new;
  }
  result.status=Status::NoConvergence;
  return result;
}

inline Status AdvanceVacuumYee1D(std::vector<double>& electric,
                                 std::vector<double>& magnetic_half,
                                 double courant) {
  if (electric.size()<3 || electric.size()!=magnetic_half.size() ||
      !std::isfinite(courant) || std::abs(courant)>1.0)
    return Status::InvalidArgument;
  for (double value:electric) if (!std::isfinite(value)) return Status::InvalidArgument;
  for (double value:magnetic_half) if (!std::isfinite(value)) return Status::InvalidArgument;
  const std::size_t n=electric.size();
  for (std::size_t i=0;i<n;++i)
    magnetic_half[i]-=courant*(electric[(i+1)%n]-electric[i]);
  std::vector<double> next=electric;
  for (std::size_t i=0;i<n;++i)
    next[i]-=courant*(magnetic_half[i]-magnetic_half[(i+n-1)%n]);
  electric.swap(next);
  return Status::Success;
}

inline double YeeInvariant(const std::vector<double>& electric,
                           const std::vector<double>& magnetic_new_half,
                           const std::vector<double>& magnetic_old_half) {
  if (electric.size()!=magnetic_new_half.size() ||
      electric.size()!=magnetic_old_half.size())
    return std::numeric_limits<double>::quiet_NaN();
  double value=0.0;
  for (std::size_t i=0;i<electric.size();++i)
    value+=0.5*(electric[i]*electric[i]+
                magnetic_new_half[i]*magnetic_old_half[i]);
  return value;
}

inline double ParticleFieldEnergy(double electric_energy,double magnetic_energy,
                                  const std::vector<double>& mass,
                                  const std::vector<double>& speed_squared) {
  if (!std::isfinite(electric_energy) || !std::isfinite(magnetic_energy) ||
      electric_energy<0.0 || magnetic_energy<0.0 ||
      mass.size()!=speed_squared.size())
    return std::numeric_limits<double>::quiet_NaN();
  double total=electric_energy+magnetic_energy;
  for (std::size_t i=0;i<mass.size();++i) {
    if (!std::isfinite(mass[i]) || !std::isfinite(speed_squared[i]) ||
        mass[i]<0.0 || speed_squared[i]<0.0)
      return std::numeric_limits<double>::quiet_NaN();
    total+=0.5*mass[i]*speed_squared[i];
  }
  return total;
}

}} // namespace PIC::FieldSolverKernels
#endif

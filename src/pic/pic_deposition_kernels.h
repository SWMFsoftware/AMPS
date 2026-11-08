#ifndef _AMPS_PIC_DEPOSITION_KERNELS_H_
#define _AMPS_PIC_DEPOSITION_KERNELS_H_

/*
 * Dependency-light particle-to-grid deposition kernels.
 *
 * Coordinates xi are dimensionless cell coordinates in [0,1].  Corner
 * number c=ix+2*iy+4*iz follows the AMPS corner-node convention.  The shape
 * function is the tensor product cloud-in-cell (CIC) basis
 *
 *   S_c(xi,eta,zeta) = sx(ix) sy(iy) sz(iz),
 *   sx(0)=1-xi, sx(1)=xi,
 *
 * so its partition of unity preserves the represented charge q and the
 * velocity-current sum q*v.  These routines do not divide by a control
 * volume: their outputs are nodal charge [C] and charge flux [C m s^-1].
 * Callers that store densities must apply their mesh dual volume exactly
 * once.
 *
 * DepositChargeConservingCurrent1D implements the finite-volume identity
 * for a CIC particle segment on a uniform one-dimensional node lattice:
 *
 *   (Q_i^{n+1}-Q_i^n)/dt + (I_{i+1/2}-I_{i-1/2}) = 0.
 *
 * Here I is charge per time [C s^-1].  The prefix construction is algebraic,
 * supports crossings of any number of cells, and avoids a tolerance-based
 * geometric case split.  It is a segment-integrated current, not the
 * point-current q*v used by the implicit ECSIM susceptibility deposition.
 * The two must not be silently interchanged.
 */

#include <cmath>
#include <cstddef>
#include <vector>

namespace PIC { namespace DepositionKernels {

enum class Status {
  Success,
  InvalidArgument,
  OutsideDomain
};

inline Status TrilinearWeights(const double xi[3],double weight[8],
                               double tolerance=32.0e-15) {
  if (xi==nullptr || weight==nullptr || !std::isfinite(tolerance) ||
      tolerance<0.0) return Status::InvalidArgument;

  double u[3];
  for (int d=0;d<3;++d) {
    if (!std::isfinite(xi[d])) return Status::InvalidArgument;
    if (xi[d]<-tolerance || xi[d]>1.0+tolerance)
      return Status::OutsideDomain;
    // Roundoff at a cell face must not create a negative nodal charge.
    u[d]=xi[d]<0.0 ? 0.0 : (xi[d]>1.0 ? 1.0 : xi[d]);
  }

  const double x0=1.0-u[0],x1=u[0];
  const double y0=1.0-u[1],y1=u[1];
  const double z0=1.0-u[2],z1=u[2];
  weight[0]=x0*y0*z0;
  weight[1]=x1*y0*z0;
  weight[2]=x0*y1*z0;
  weight[3]=x1*y1*z0;
  weight[4]=x0*y0*z1;
  weight[5]=x1*y0*z1;
  weight[6]=x0*y1*z1;
  weight[7]=x1*y1*z1;
  return Status::Success;
}

inline Status DepositChargeAndVelocityCurrent(
    const double xi[3],double charge,const double velocity[3],
    double nodal_charge[8],double nodal_qv[8][3]) {
  if (velocity==nullptr || nodal_charge==nullptr || nodal_qv==nullptr ||
      !std::isfinite(charge) || !std::isfinite(velocity[0]) ||
      !std::isfinite(velocity[1]) || !std::isfinite(velocity[2]))
    return Status::InvalidArgument;

  double weight[8];
  const Status status=TrilinearWeights(xi,weight);
  if (status!=Status::Success) return status;
  for (int c=0;c<8;++c) {
    nodal_charge[c]=charge*weight[c];
    for (int d=0;d<3;++d) nodal_qv[c][d]=nodal_charge[c]*velocity[d];
  }
  return Status::Success;
}

inline Status DepositCICCharge1D(double x,double x_min,double dx,
                                 std::size_t node_count,double charge,
                                 std::vector<double>& nodal_charge) {
  if (!std::isfinite(x) || !std::isfinite(x_min) || !std::isfinite(dx) ||
      !std::isfinite(charge) || dx<=0.0 || node_count<2)
    return Status::InvalidArgument;
  const double coordinate=(x-x_min)/dx;
  if (coordinate<0.0 || coordinate>static_cast<double>(node_count-1))
    return Status::OutsideDomain;

  nodal_charge.assign(node_count,0.0);
  std::size_t left=static_cast<std::size_t>(std::floor(coordinate));
  double fraction=coordinate-static_cast<double>(left);
  if (left==node_count-1) {
    left=node_count-2;
    fraction=1.0;
  }
  nodal_charge[left]=charge*(1.0-fraction);
  nodal_charge[left+1]=charge*fraction;
  return Status::Success;
}

inline Status DepositChargeConservingCurrent1D(
    double x_old,double x_new,double x_min,double dx,std::size_t node_count,
    double charge,double dt,std::vector<double>& old_charge,
    std::vector<double>& new_charge,std::vector<double>& edge_current) {
  if (!std::isfinite(dt) || dt<=0.0) return Status::InvalidArgument;
  Status status=DepositCICCharge1D(
      x_old,x_min,dx,node_count,charge,old_charge);
  if (status!=Status::Success) return status;
  status=DepositCICCharge1D(x_new,x_min,dx,node_count,charge,new_charge);
  if (status!=Status::Success) return status;

  edge_current.assign(node_count-1,0.0);
  double cumulative_change=0.0;
  for (std::size_t edge=0;edge+1<node_count;++edge) {
    cumulative_change+=new_charge[edge]-old_charge[edge];
    edge_current[edge]=-cumulative_change/dt;
  }
  return Status::Success;
}

}} // namespace PIC::DepositionKernels

#endif

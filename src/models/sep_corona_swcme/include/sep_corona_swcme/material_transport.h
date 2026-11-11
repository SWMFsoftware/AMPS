#ifndef SEP_CORONA_SWCME_MATERIAL_TRANSPORT_H
#define SEP_CORONA_SWCME_MATERIAL_TRANSPORT_H

#include "sep_corona_swcme/cross_ray_transport.h"

#include <array>

namespace SEP { namespace CoronaSwcme { namespace Material {

using Matrix3=std::array<double,9>;

struct MapJet {
  // A SINGLE prescribed map X(a,t), with a in reference HCI metres. These
  // derivatives must come from that same map; this pointwise kernel cannot
  // prove global spatial integrability, coverage, topology or map ownership.
  CoronalCME::Vec3 positionHciM,velocityHciMPerS,accelerationHciMPerS2;
  Matrix3 deformation{};           // F_ij=dX_i/da_j, dimensionless, row major
  Matrix3 deformationRatePerS{};    // dF/dt at fixed material label, 1/s
};

struct TransportedState {
  CrossRay::Primitive primitive;
  double jacobian=0;
  Matrix3 velocityGradientPerS{};   // grad_X v = Fdot * inverse(F), 1/s
  double densityMaterialRateKgM3S=0,pressureMaterialRatePaPerS=0;
  CoronalCME::Vec3 magneticMaterialRateTPerS;
};

// Exact local mass, ideal frozen-flux and ADIABATIC transport for a regular
// orientation-preserving material map:
// rho=rho0/J, B=F*B0/J, p=p0*J^(-gamma), J=det F.
// Jdot/J=tr(Fdot*F^-1), DB/Dt=(Fdot*B0)/J-B*Jdot/J.
// All vectors use HCI Cartesian components and SI. No density/field floors,
// diagonal-only flux scaling or passive radial B approximation is applied.
//
// This is the material-transport building block for the user-selected
// analytical/semi-empirical closure, NOT an MHD runtime or a choice of the
// missing velocity/material map. It neither supplies a computed shock nor
// qualifies momentum, heating, interfaces, inverse maps or native storage.
// Nonadiabatic work/heating requires a separately specified closure; it must
// not be hidden by changing this adiabatic pressure law.
//
// Scaled determinant/cofactor algebra checks numerical rank before inverse
// construction. Tiny positive J is allowed if representable and resolved;
// physical minimum-J/force/work gates remain separate and unchanged. Negative
// or numerically unresolved maps reject instead of regularizing a fold.
Core::Result<TransportedState> PushForwardAdiabatic(
    const CrossRay::Primitive& reference,const MapJet& map,
    const CrossRay::Thermodynamics& gas);

} } } // namespace SEP::CoronaSwcme::Material

#endif

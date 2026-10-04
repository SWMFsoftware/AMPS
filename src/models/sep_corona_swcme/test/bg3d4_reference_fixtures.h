#ifndef SEP_CORONA_SWCME_BG3D4_REFERENCE_FIXTURES_H
#define SEP_CORONA_SWCME_BG3D4_REFERENCE_FIXTURES_H

#include "sep_coronal_cme/mhd_jump_solver.h"

#include <cmath>

namespace BG3D4Reference {

using SEP::CoronalCME::Dot;
using SEP::CoronalCME::MhdPrimitiveState;
using SEP::CoronalCME::Vec3;

constexpr double kPi=3.141592653589793238462643383279502884;
constexpr double kMu0=1.25663706212e-6;

// Exact finite-inventory planar shock/piston solution.  At t=0 the downstream
// slab -H<=a<=0 already exists.  Its left boundary is a material piston/contact
// Xc=-H+U2*t, while the shock Xs=Vs*t admits new material.  Existing labels map
// as x=a+U2*t and shock-born labels as x=Vs*tau+U2*(t-tau).  The selected
// gamma=5/3, rho1=2, p1=1.2, w1=3 solution has compression 3, w2=1 and
// p2=13.2.  A uniform normal magnetic field is nonzero and divergence-free;
// it cancels from the normal RH momentum/energy jumps and remains continuous.
// Gravity and external forcing are absent.  This is a manufactured local MHD
// verification problem, not a global CME model.
struct FinitePlanarShockPiston {
  double gamma=5.0/3.0;
  double areaM2=7.0;
  double initialThicknessM=5.0;
  double shockSpeedMPerS=4.0;
  MhdPrimitiveState upstream{2.0,1.2,{1.0,0.0,0.0},{2e-9,0.0,0.0}};
  MhdPrimitiveState downstream{6.0,13.2,{3.0,0.0,0.0},{2e-9,0.0,0.0}};
  // A transverse-field ideal contact can use this ejecta-side state: U, p and
  // B are continuous while density/entropy differ.  It therefore supplies an
  // independently specified two-sided boundary without asserting flux-rope
  // topology.
  MhdPrimitiveState pistonSide{4.0,13.2,{3.0,0.0,0.0},{2e-9,0.0,0.0}};

  double UpstreamInflow() const {
    return shockSpeedMPerS-upstream.velocityMPerS.x;
  }
  double DownstreamInflow() const {
    return shockSpeedMPerS-downstream.velocityMPerS.x;
  }
  double ContactPosition(double timeS) const {
    return -initialThicknessM+downstream.velocityMPerS.x*timeS;
  }
  double ShockPosition(double timeS) const { return shockSpeedMPerS*timeS; }
  double InitialMap(double labelM,double timeS) const {
    return labelM+downstream.velocityMPerS.x*timeS;
  }
  double AdmittedMap(double crossingTimeS,double timeS) const {
    return shockSpeedMPerS*crossingTimeS+
        downstream.velocityMPerS.x*(timeS-crossingTimeS);
  }
  double InitialMassKg() const {
    return downstream.massDensityKgM3*areaM2*initialThicknessM;
  }
  double AdmittedMassKg(double timeS) const {
    return upstream.massDensityKgM3*UpstreamInflow()*areaM2*timeS;
  }
  double VolumeM3(double timeS) const {
    return areaM2*(initialThicknessM+
        (shockSpeedMPerS-downstream.velocityMPerS.x)*timeS);
  }
  double TotalMassKg(double timeS) const {
    return downstream.massDensityKgM3*VolumeM3(timeS);
  }
  double EnergyDensity(const MhdPrimitiveState& state) const {
    return 0.5*state.massDensityKgM3*Dot(state.velocityMPerS,
        state.velocityMPerS)+state.pressurePa/(gamma-1)+
        Dot(state.magneticFieldT,state.magneticFieldT)/(2*kMu0);
  }
  double EnergyFluxX(const MhdPrimitiveState& state) const {
    const double energy=EnergyDensity(state);
    const double magnetic=Dot(state.magneticFieldT,state.magneticFieldT)/(2*kMu0);
    return (energy+state.pressurePa+magnetic)*state.velocityMPerS.x-
        Dot(state.velocityMPerS,state.magneticFieldT)*state.magneticFieldT.x/kMu0;
  }
  double EnergyRateW() const {
    return EnergyDensity(downstream)*areaM2*
        (shockSpeedMPerS-downstream.velocityMPerS.x);
  }
  double MovingBoundaryEnergyPowerW() const {
    const double energy=EnergyDensity(downstream);
    const double flux=EnergyFluxX(downstream);
    // -integral (F-E*V_boundary).n dA at the right shock (n=+x)
    // and left material piston (n=-x).  Keeping both terms prevents pressure
    // work at the contact from being silently omitted.
    return areaM2*((energy*shockSpeedMPerS-flux)+
        (flux-energy*downstream.velocityMPerS.x));
  }
};

// Exact curved material shell.  Reference labels are a=s*e_r with
// Rc0<=s<=Rs0 and X(a,t)=C+lambda(t)a, lambda=1+alpha*t.  Thus F=lambda I,
// J=lambda^3, U=(alpha/lambda)(X-C), rho=rho0/J,
// p=p0*J^-gamma and B=B0/lambda^2.  B0 is a uniform Cartesian field: it is
// nonzero, divergence-free, and threads both material boundaries.  The
// interface topology is the declared transverse-field ideal contact, not a
// closed flux rope.  Because lambda''=0 and all scalar/vector fields are
// spatially uniform except U, inertia, grad(p), curl(B), and manufactured
// body forcing are identically zero.  Expansion energy changes only through
// the analytically known isotropic-pressure and Maxwell boundary work.
struct HomologousCurvedShell {
  double innerRadiusM=2.0;
  double outerRadiusM=3.5;
  double alphaPerS=0.02;
  double rho0KgM3=5.0;
  double pressure0Pa=11.0;
  double gamma=5.0/3.0;
  Vec3 magnetic0T{2e-3,-3e-3,4e-3};

  double Lambda(double timeS) const { return 1+alphaPerS*timeS; }
  Vec3 Map(Vec3 referenceM,double timeS) const {
    return Lambda(timeS)*referenceM;
  }
  Vec3 Velocity(Vec3 referenceM,double) const { return alphaPerS*referenceM; }
  double Jacobian(double timeS) const {
    return std::pow(Lambda(timeS),3);
  }
  double Density(double timeS) const { return rho0KgM3/Jacobian(timeS); }
  double Pressure(double timeS) const {
    return pressure0Pa*std::pow(Jacobian(timeS),-gamma);
  }
  Vec3 MagneticField(double timeS) const {
    const double lambda=Lambda(timeS);
    return magnetic0T/(lambda*lambda);
  }
  double ReferenceVolumeM3() const {
    return 4*kPi*(std::pow(outerRadiusM,3)-std::pow(innerRadiusM,3))/3;
  }
  double VolumeM3(double timeS) const {
    return Jacobian(timeS)*ReferenceVolumeM3();
  }
  double UnweightedInnerAreaColumnM3(double timeS) const {
    const double lambda=Lambda(timeS);
    const double rc=lambda*innerRadiusM,rs=lambda*outerRadiusM;
    return 4*kPi*rc*rc*(rs-rc);
  }
  double InternalEnergyJ(double timeS) const {
    return Pressure(timeS)*VolumeM3(timeS)/(gamma-1);
  }
  double MagneticEnergyJ(double timeS) const {
    const Vec3 field=MagneticField(timeS);
    return Dot(field,field)*VolumeM3(timeS)/(2*kMu0);
  }
  double KineticEnergyJ() const {
    // Integral .5*rho0*alpha^2*s^2 dV0 over the reference shell.
    return 2*kPi*rho0KgM3*alphaPerS*alphaPerS*
        (std::pow(outerRadiusM,5)-std::pow(innerRadiusM,5))/5;
  }
  double TotalEnergyJ(double timeS) const {
    return KineticEnergyJ()+InternalEnergyJ(timeS)+MagneticEnergyJ(timeS);
  }
  double BoundaryWorkW(double timeS) const {
    const double expansion=alphaPerS/Lambda(timeS);
    return -expansion*VolumeM3(timeS)*(3*Pressure(timeS)+
        Dot(MagneticField(timeS),MagneticField(timeS))/(2*kMu0));
  }
};

} // namespace BG3D4Reference

#endif

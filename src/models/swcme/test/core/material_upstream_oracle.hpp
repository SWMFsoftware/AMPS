#ifndef SWCME_TEST_MATERIAL_UPSTREAM_ORACLE_HPP
#define SWCME_TEST_MATERIAL_UPSTREAM_ORACLE_HPP

// Independent D10-1 references for the Stage12 prescribed upstream contract.
// This header is TEST ONLY. It never replaces a production evaluator and does
// not import the production material profile, its gradients, force helpers,
// map, shock residuals, or cached samples. SI literals below are the declared
// verification inputs and adopted constants, not fitted candidate output.
// Long-double arithmetic reduces reference rounding when checking binary64
// primitive values and cancellation-sensitive vector balances.
#include <array>
#include <cmath>
#include <limits>

namespace swcme_material_oracle {
using Scalar=long double;
using Vector=std::array<Scalar,3>;
using Tensor=std::array<Scalar,9>;
constexpr Scalar AU=149597870700.0L;
constexpr Scalar ProtonMass=1.67262192595e-27L;
constexpr Scalar Boltzmann=1.380649e-23L;
constexpr Scalar Mu0=1.25663706127e-6L;
constexpr Scalar Gravity=1.32712440018e20L;
constexpr Scalar Speed=400000.0L;
constexpr Scalar ElectronDensity=6000000.0L;
constexpr Scalar ElectronTemperature=120000.0L;
constexpr Scalar ProtonTemperature=120000.0L;
constexpr Scalar FieldMagnitude=5.0e-9L;
constexpr Scalar Rotation=2.86533e-6L;
constexpr Scalar Gamma=5.0L/3.0L;
// Frozen BEFORE inspecting candidate results. These are the approved D10-1
// budget and inherited short-arithmetic/core derivative criteria; no per-case
// residual is allowed to select or enlarge a threshold after execution.
constexpr Scalar ExactTolerance=512.0L*std::numeric_limits<double>::epsilon();
constexpr Scalar BudgetTolerance=1.0e-10L;
constexpr Scalar MinimumDerivativeOrder=1.90L;
constexpr Scalar MaximumDerivativeOrder=2.10L;
constexpr Scalar ShockMassTolerance=1.0e-9L;
constexpr Scalar ShockNormalBTolerance=1.0e-10L;
constexpr Scalar ShockElectricTolerance=1.0e-8L;
constexpr Scalar ShockMomentumTolerance=1.0e-8L;
constexpr Scalar ShockEnergyTolerance=1.0e-8L;

inline Scalar dot(const Vector& a,const Vector& b) {
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
inline Scalar norm(const Vector& a) {return std::sqrt(dot(a,a));}
inline Vector cross(const Vector& a,const Vector& b) {
  return {{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}};
}
inline Vector scaled(Scalar factor,const Vector& a) {
  return {{factor*a[0],factor*a[1],factor*a[2]}};
}
inline Vector subtract(const Vector& a,const Vector& b) {
  return {{a[0]-b[0],a[1]-b[1],a[2]-b[2]}};
}
inline Vector product(const Tensor& a,const Vector& x) {
  Vector y{{0,0,0}};
  for(unsigned i=0;i<3;++i)for(unsigned j=0;j<3;++j)y[i]+=a[3*i+j]*x[j];
  return y;
}
struct State {
  Vector position,radial,velocity,magnetic,pressureGradient,lorentz,force;
  Tensor velocityGradient,magneticGradient;
  Scalar radius,electrons,protons,density,electronTemperature,protonTemperature;
  Scalar pressure,soundSpeed,radialField,work,heat,divergence;
};

inline State evaluate(const Vector& x,int polarity,const Vector& unitAxis) {
  State s{};s.position=x;s.radius=norm(x);s.radial=scaled(1/s.radius,x);
  const Scalar ratio=AU/s.radius;
  s.electrons=ElectronDensity*ratio*ratio;s.protons=s.electrons;
  s.density=ProtonMass*s.protons;
  const Scalar temperatureFactor=std::pow(ratio,4.0L/3.0L);
  s.electronTemperature=ElectronTemperature*temperatureFactor;
  s.protonTemperature=ProtonTemperature*temperatureFactor;
  s.pressure=Boltzmann*(s.protons*s.protonTemperature+s.electrons*s.electronTemperature);
  s.soundSpeed=std::sqrt(Gamma*s.pressure/s.density);
  s.velocity=scaled(Speed,s.radial);s.divergence=2*Speed/s.radius;

  // The 5nT magnitude is anchored ONCE at the declared reference equator.
  // Changing query latitude must not renormalize the field on each ray. The
  // Cartesian expression remains regular at the rotation axis, where e_phi is
  // undefined but both the physical field and all derivatives remain finite.
  const Scalar radialAtReference=polarity*FieldMagnitude/
      std::sqrt(1+std::pow(Rotation*AU/Speed,2));
  const Scalar a=radialAtReference*AU*AU;
  const Scalar c=a*Rotation/Speed;
  const Vector winding=cross(unitAxis,x);
  s.radialField=a/(s.radius*s.radius);
  for(unsigned i=0;i<3;++i)
    s.magnetic[i]=a*x[i]/std::pow(s.radius,3)-c*winding[i]/std::pow(s.radius,2);
  const Tensor rotationMatrix{{0,-unitAxis[2],unitAxis[1],unitAxis[2],0,-unitAxis[0],-unitAxis[1],unitAxis[0],0}};
  for(unsigned i=0;i<3;++i)for(unsigned j=0;j<3;++j) {
    const Scalar identity=i==j?1:0;
    s.velocityGradient[3*i+j]=Speed/s.radius*(identity-s.radial[i]*s.radial[j]);
    s.magneticGradient[3*i+j]=a*(identity/std::pow(s.radius,3)-3*x[i]*x[j]/std::pow(s.radius,5))-
        c*(rotationMatrix[3*i+j]/std::pow(s.radius,2)-2*winding[i]*x[j]/std::pow(s.radius,4));
  }

  // A constant radial speed has zero material acceleration, not zero force.
  // Closed curl(B) gives a transverse Lorentz force away from equator/poles.
  // Compute it independently of the gradient tensor above so an omitted or
  // transposed derivative cannot validate a force error against itself.
  const Scalar latitudeProjection=dot(unitAxis,s.radial);
  const Scalar lorentzFactor=2*c*c*latitudeProjection/(Mu0*std::pow(s.radius,3));
  for(unsigned i=0;i<3;++i) {
    s.lorentz[i]=lorentzFactor*(unitAxis[i]-latitudeProjection*s.radial[i]);
    s.pressureGradient[i]=-2*Gamma*s.pressure/s.radius*s.radial[i];
    s.force[i]=s.pressureGradient[i]+s.density*Gravity/(s.radius*s.radius)*s.radial[i]-s.lorentz[i];
  }
  // Lorentz is transverse to the radial wind and therefore does no work in
  // this prescribed realization. Mechanical sustaining work remains separate
  // from the exactly zero thermal source required by adiabatic entropy flow.
  s.work=Speed*(s.density*Gravity/(s.radius*s.radius)-2*Gamma*s.pressure/s.radius);
  s.heat=0;
  return s;
}

struct MaterialMap {Vector position;Tensor deformation;Scalar jacobian;};
inline MaterialMap move(const Vector& initial,Scalar elapsed) {
  MaterialMap m{};const Scalar a=norm(initial);const Vector e=scaled(1/a,initial);
  const Scalar r=a+Speed*elapsed,transverseStretch=r/a;
  m.position=scaled(r,e);m.jacobian=transverseStretch*transverseStretch;
  // The map leaves radial separations unchanged and stretches BOTH transverse
  // directions by r/a. This eigenprojector is the independent full Cartesian
  // deformation; J is not inferred from any production density ratio.
  for(unsigned i=0;i<3;++i)for(unsigned j=0;j<3;++j)
    m.deformation[3*i+j]=e[i]*e[j]+transverseStretch*((i==j?1:0)-e[i]*e[j]);
  return m;
}

inline Scalar massPerSolidAngle(Scalar inner,Scalar outer) {
  return ProtonMass*ElectronDensity*AU*AU*(outer-inner);
}
inline Scalar massFluxPerSolidAngle() {return ProtonMass*ElectronDensity*AU*AU*Speed;}
inline Scalar entropy(const State& s) {return s.pressure/std::pow(s.density,Gamma);}
inline Scalar energyFluxPerSolidAngle(const State& s) {
  const Scalar magnetic2=dot(s.magnetic,s.magnetic);
  const Scalar radialMagnetic=dot(s.magnetic,s.radial);
  const Scalar kinetic=0.5L*s.density*Speed*Speed*Speed;
  const Scalar enthalpy=Gamma/(Gamma-1)*s.pressure*Speed;
  const Scalar poynting=Speed*(magnetic2-radialMagnetic*radialMagnetic)/Mu0;
  return s.radius*s.radius*(kinetic+enthalpy+poynting);
}
inline Scalar shellGravityPowerPerSolidAngle(Scalar inner,Scalar outer) {
  return massFluxPerSolidAngle()*Gravity*(1/outer-1/inner);
}
inline Scalar shellPressureWorkPerSolidAngle(Scalar inner,Scalar outer) {
  const Scalar referencePressure=Boltzmann*ElectronDensity*(ElectronTemperature+ProtonTemperature);
  const Scalar exponent=2-2*Gamma;
  return Gamma/(Gamma-1)*Speed*referencePressure*AU*AU*
      (std::pow(outer/AU,exponent)-std::pow(inner/AU,exponent));
}
} // namespace swcme_material_oracle
#endif

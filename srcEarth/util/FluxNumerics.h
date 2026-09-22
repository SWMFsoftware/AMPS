#ifndef _SRC_EARTH_UTIL_FLUX_NUMERICS_H_
#define _SRC_EARTH_UTIL_FLUX_NUMERICS_H_

//======================================================================================
// FluxNumerics.h
//======================================================================================
// Shared, field-backend-independent numerical kernels for cutoff-derived particle
// products.  Both gridless and Mode3D/SWMF density solvers include this file.
//
// Why this module exists
// ----------------------
// The two solvers used to carry separate copies of the relativistic conversions,
// energy-grid builder, angular grid, channel clipping, and trapezoidal quadrature.  A
// copy in Mode3D also contained a quadratic factor in the LINEAR energy grid.  Keeping
// these operations here makes a backend comparison test a comparison of field/tracing
// physics rather than a comparison of two slightly different post-processing codes.
//
// Unit contract
// -------------
//   kinetic energy passed to spectra/integrators : joule [J]
//   public energy-grid coordinates               : megaelectronvolt [MeV]
//   rigidity                                     : gigavolt [GV]
//   mass                                         : kilogram [kg]
//   charge                                       : absolute coulomb [C]
//   differential intensity J(E)                  : m^-2 s^-1 sr^-1 J^-1
//   omnidirectional integral flux                : m^-2 s^-1
//   number density                               : m^-3
//
// Unresolved-trajectory contract
// ------------------------------
// Numerical limits and invalid-field/mover failures are not physical FORBIDDEN
// classifications.  AccessAccumulator therefore retains sampled, resolved, allowed,
// retried, and termination counts.  ResolveAccess() reports three transmissivities:
//
//   nominal = allowed weighted sum / resolved count
//   lower   = allowed weighted sum / sampled count
//   upper   = (allowed weighted sum + unresolved*maximum weight) / sampled count
//
// The bounds mean "all unresolved are forbidden" and "all unresolved are maximally
// allowed."  If no trajectory resolves, nominal is NaN--never a physical zero--while
// the conservative [lower,upper] interval remains available.  Folding each curve with
// the same non-negative boundary spectrum produces corresponding density/flux bounds.
//======================================================================================

#include "TrajectoryTermination.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace Earth {
namespace FluxNumerics {

static constexpr double kPi = 3.141592653589793238462643383279502884;
static constexpr double kSpeedOfLight_m_s = 299792458.0;
static constexpr double kElementaryCharge_C = 1.602176634e-19;
static constexpr double kAtomicMassUnit_kg = 1.66053906660e-27;
static constexpr double kMeVToJ = 1.0e6*kElementaryCharge_C;

inline double RigidityFromEnergyGV(double kineticEnergy_J,
                                   double absoluteCharge_C,
                                   double restMass_kg) {
  if (!(kineticEnergy_J >= 0.0) || !(absoluteCharge_C > 0.0) ||
      !(restMass_kg > 0.0)) return std::numeric_limits<double>::quiet_NaN();
  const double mc2 = restMass_kg*kSpeedOfLight_m_s*kSpeedOfLight_m_s;
  // Algebraically pc=sqrt((K+mc^2)^2-(mc^2)^2).  The factored K*(K+2mc^2)
  // form avoids subtracting nearly equal rest-energy terms at low kinetic energy.
  const double pc_J = std::sqrt(std::max(0.0,kineticEnergy_J*(kineticEnergy_J+2.0*mc2)));
  return (pc_J/absoluteCharge_C)*1.0e-9;
}

inline double EnergyFromRigidityMeV(double rigidity_GV,
                                    double absoluteCharge_C,
                                    double restMass_kg) {
  if (!(rigidity_GV >= 0.0) || !(absoluteCharge_C > 0.0) ||
      !(restMass_kg > 0.0)) return std::numeric_limits<double>::quiet_NaN();
  const double mc2 = restMass_kg*kSpeedOfLight_m_s*kSpeedOfLight_m_s;
  const double pc_J = rigidity_GV*1.0e9*absoluteCharge_C;
  const double totalEnergy_J=std::sqrt(pc_J*pc_J+mc2*mc2);
  // K=(pc)^2/(E_total+mc^2) is the cancellation-free form of E_total-mc^2.
  return (pc_J*pc_J/(totalEnergy_J+mc2))/kMeVToJ;
}

inline double RelativisticSpeed(double kineticEnergy_J,double restMass_kg) {
  if (!(kineticEnergy_J >= 0.0) || !(restMass_kg > 0.0))
    return std::numeric_limits<double>::quiet_NaN();
  const double mc2 = restMass_kg*kSpeedOfLight_m_s*kSpeedOfLight_m_s;
  const double gamma = 1.0+kineticEnergy_J/mc2;
  if (gamma <= 1.0) return 0.0;
  return kSpeedOfLight_m_s*std::sqrt(std::max(0.0,1.0-1.0/(gamma*gamma)));
}

enum class EnergySpacing { Linear, Log };

// Build the fixed energy coordinate used by both output and quadrature.  In rigidity
// scan mode, equal log-R nodes are converted back to kinetic energy because magnetic
// access is naturally organized in rigidity but spectra are functions of energy.
inline std::vector<double> BuildEnergyGridMeV(double Emin_MeV,double Emax_MeV,
                                               int legacyPointCount,
                                               EnergySpacing spacing,
                                               bool rigidityScan,
                                               int rigidityScanPointCount,
                                               int maximumPointCount,
                                               double absoluteCharge_C,
                                               double restMass_kg) {
  if (!(Emin_MeV > 0.0) || !(Emax_MeV > Emin_MeV) || legacyPointCount < 2)
    throw std::invalid_argument("BuildEnergyGridMeV: invalid energy range or point count");

  int n = legacyPointCount;
  if (rigidityScan) {
    if (rigidityScanPointCount > 0) n=rigidityScanPointCount;
    if (maximumPointCount > 0) n=std::min(n,maximumPointCount);
    n=std::max(2,n);
  }

  std::vector<double> E(static_cast<std::size_t>(n),0.0);
  if (rigidityScan) {
    const double Rmin=RigidityFromEnergyGV(Emin_MeV*kMeVToJ,absoluteCharge_C,restMass_kg);
    const double Rmax=RigidityFromEnergyGV(Emax_MeV*kMeVToJ,absoluteCharge_C,restMass_kg);
    if (!(Rmin > 0.0) || !(Rmax > Rmin))
      throw std::invalid_argument("BuildEnergyGridMeV: invalid charge, mass, or rigidity range");
    const double logRmin=std::log(Rmin), logRmax=std::log(Rmax);
    for (int i=0;i<n;++i) {
      const double a=static_cast<double>(i)/static_cast<double>(n-1);
      const double R=std::exp(logRmin+a*(logRmax-logRmin));
      E[static_cast<std::size_t>(i)]=EnergyFromRigidityMeV(R,absoluteCharge_C,restMass_kg);
    }
  }
  else if (spacing==EnergySpacing::Log) {
    const double logMin=std::log(Emin_MeV), logMax=std::log(Emax_MeV);
    for (int i=0;i<n;++i) {
      const double a=static_cast<double>(i)/static_cast<double>(n-1);
      E[static_cast<std::size_t>(i)]=std::exp(logMin+a*(logMax-logMin));
    }
  }
  else {
    // Deliberately one factor of a.  The former Mode3D copy used a*a here and
    // silently produced a quadratic, not LINEAR, energy coordinate.
    for (int i=0;i<n;++i) {
      const double a=static_cast<double>(i)/static_cast<double>(n-1);
      E[static_cast<std::size_t>(i)]=Emin_MeV+a*(Emax_MeV-Emin_MeV);
    }
  }
  // Pin endpoints to input values so channel clipping and backend regression files are
  // bitwise-consistent despite inverse conversion roundoff in rigidity-scan mode.
  E.front()=Emin_MeV;
  E.back()=Emax_MeV;
  return E;
}

struct DirectionSample {
  double x{0.0},y{0.0},z{0.0};
  double solidAngleWeight_sr{0.0};
  DirectionSample() = default;
  DirectionSample(double xIn,double yIn,double zIn,double weightIn)
      : x(xIn),y(yIn),z(zIn),solidAngleWeight_sr(weightIn) {}
};

// Midpoint sampling uniform in mu=cos(theta) and phi.  Every cell has the same solid
// angle, 4*pi/(N_mu*N_phi), and the weights sum to 4*pi to roundoff.
inline std::vector<DirectionSample> BuildEqualSolidAngleDirections(int nMu,int nPhi) {
  if (nMu < 1 || nPhi < 1)
    throw std::invalid_argument("BuildEqualSolidAngleDirections: dimensions must be positive");
  std::vector<DirectionSample> result;
  result.reserve(static_cast<std::size_t>(nMu)*static_cast<std::size_t>(nPhi));
  const double w=4.0*kPi/static_cast<double>(nMu*nPhi);
  for (int i=0;i<nMu;++i) {
    const double mu=-1.0+2.0*(static_cast<double>(i)+0.5)/static_cast<double>(nMu);
    const double st=std::sqrt(std::max(0.0,1.0-mu*mu));
    for (int j=0;j<nPhi;++j) {
      const double phi=2.0*kPi*(static_cast<double>(j)+0.5)/static_cast<double>(nPhi);
      result.push_back(DirectionSample{st*std::cos(phi),st*std::sin(phi),mu,w});
    }
  }
  return result;
}

template<class T>
inline std::vector<T> SelectDeterministic(const std::vector<T>& full,int nUse) {
  if (full.empty()) return std::vector<T>();
  nUse=std::max(1,nUse);
  if (nUse >= static_cast<int>(full.size())) return full;
  std::vector<T> result;
  result.reserve(static_cast<std::size_t>(nUse));
  if (nUse==1) {
    result.push_back(full[full.size()/2]);
    return result;
  }
  for (int k=0;k<nUse;++k) {
    const double a=static_cast<double>(k)/static_cast<double>(nUse-1);
    const int index=static_cast<int>(std::floor(a*static_cast<double>(full.size()-1)+0.5));
    result.push_back(full[static_cast<std::size_t>(index)]);
  }
  return result;
}

static constexpr int kTerminationCount=
    static_cast<int>(GridlessMode::TrajectoryTermination::Count);

struct AccessAccumulator {
  double weightSum{0.0};
  int sampled{0};
  int resolved{0};
  int allowed{0};
  int retried{0};
  std::array<int,kTerminationCount> terminationCounts{};

  int unresolved() const { return sampled-resolved; }
  double transmission() const {
    return resolved>0 ? weightSum/static_cast<double>(resolved)
                      : std::numeric_limits<double>::quiet_NaN();
  }

  void Record(GridlessMode::TrajectoryTermination termination,
              double allowedTrajectoryWeight,bool wasRetried=false) {
    ++sampled;
    if (wasRetried) ++retried;
    const int code=static_cast<int>(termination);
    if (code>=0 && code<kTerminationCount)
      ++terminationCounts[static_cast<std::size_t>(code)];
    if (!GridlessMode::IsResolvedTermination(termination)) return;
    ++resolved;
    if (!GridlessMode::IsAllowedTermination(termination)) return;
    ++allowed;
    weightSum+=allowedTrajectoryWeight;
  }

  void Merge(const AccessAccumulator& other) {
    weightSum+=other.weightSum;
    sampled+=other.sampled;
    resolved+=other.resolved;
    allowed+=other.allowed;
    retried+=other.retried;
    for (int i=0;i<kTerminationCount;++i)
      terminationCounts[static_cast<std::size_t>(i)]+=
          other.terminationCounts[static_cast<std::size_t>(i)];
  }
};

struct AccessEstimate {
  double nominal{std::numeric_limits<double>::quiet_NaN()};
  double lower{0.0};
  double upper{1.0};
  double unresolvedFraction{1.0};
  bool hasResolved{false};
};

inline AccessEstimate ResolveAccess(const AccessAccumulator& a,double maximumAllowedWeight=1.0) {
  AccessEstimate r;
  maximumAllowedWeight=std::max(0.0,maximumAllowedWeight);
  if (a.sampled <= 0) return r;
  r.hasResolved=(a.resolved>0);
  if (r.hasResolved) r.nominal=a.weightSum/static_cast<double>(a.resolved);
  r.lower=a.weightSum/static_cast<double>(a.sampled);
  r.upper=(a.weightSum+static_cast<double>(a.unresolved())*maximumAllowedWeight)/
          static_cast<double>(a.sampled);
  r.unresolvedFraction=static_cast<double>(a.unresolved())/static_cast<double>(a.sampled);
  // Clamp only roundoff excursions.  Anisotropic weights may legitimately exceed one.
  r.lower=std::max(0.0,std::min(maximumAllowedWeight,r.lower));
  r.upper=std::max(r.lower,std::min(maximumAllowedWeight,r.upper));
  return r;
}

inline bool ExceedsUnresolvedTolerance(const AccessAccumulator& a,double tolerance) {
  return ResolveAccess(a).unresolvedFraction>tolerance;
}

inline double Trapezoid(const std::vector<double>& x,const std::vector<double>& y) {
  if (x.size()!=y.size() || x.size()<2) return 0.0;
  double result=0.0;
  for (std::size_t i=0;i+1<x.size();++i)
    result+=0.5*(y[i]+y[i+1])*(x[i+1]-x[i]);
  return result;
}

inline double InterpolateLinear(const std::vector<double>& x,
                                const std::vector<double>& y,double value) {
  if (x.size()!=y.size() || x.empty()) return std::numeric_limits<double>::quiet_NaN();
  std::vector<double>::const_iterator it=std::lower_bound(x.begin(),x.end(),value);
  if (it==x.begin()) return y.front();
  if (it==x.end()) return y.back();
  const std::size_t hi=static_cast<std::size_t>(it-x.begin()), lo=hi-1;
  if (!(x[hi]>x[lo])) return y[lo];
  return y[lo]+(y[hi]-y[lo])*(value-x[lo])/(x[hi]-x[lo]);
}

// Spectrum is any callable double(double energy_J).  Channel endpoints are clipped to
// the solver grid, inserted explicitly, and T is linearly interpolated there.  The
// spectrum itself is evaluated at the inserted energy rather than interpolated.
template<class Spectrum>
inline double IntegrateFlux(const std::vector<double>& E_MeV,
                            const std::vector<double>& T,
                            double lower_MeV,double upper_MeV,
                            Spectrum spectrumPerJ) {
  if (E_MeV.size()!=T.size() || E_MeV.size()<2 || !(lower_MeV<upper_MeV)) return 0.0;
  const double lo=std::max(lower_MeV,E_MeV.front());
  const double hi=std::min(upper_MeV,E_MeV.back());
  if (!(lo<hi)) return 0.0;
  std::vector<double> nodes;
  nodes.reserve(E_MeV.size()+2);
  nodes.push_back(lo);
  for (std::size_t i=0;i<E_MeV.size();++i)
    if (E_MeV[i]>lo && E_MeV[i]<hi) nodes.push_back(E_MeV[i]);
  nodes.push_back(hi);
  double result=0.0;
  for (std::size_t i=0;i+1<nodes.size();++i) {
    const double E0_J=nodes[i]*kMeVToJ, E1_J=nodes[i+1]*kMeVToJ;
    const double f0=InterpolateLinear(E_MeV,T,nodes[i])*spectrumPerJ(E0_J);
    const double f1=InterpolateLinear(E_MeV,T,nodes[i+1])*spectrumPerJ(E1_J);
    result+=0.5*(f0+f1)*(E1_J-E0_J);
  }
  return 4.0*kPi*result;
}

template<class Spectrum>
inline double IntegrateDensity(const std::vector<double>& E_MeV,
                               const std::vector<double>& T,
                               double restMass_kg,Spectrum spectrumPerJ) {
  if (E_MeV.size()!=T.size() || E_MeV.size()<2) return 0.0;
  std::vector<double> E_J(E_MeV.size(),0.0),integrand(E_MeV.size(),0.0);
  for (std::size_t i=0;i<E_MeV.size();++i) {
    E_J[i]=E_MeV[i]*kMeVToJ;
    const double speed=RelativisticSpeed(E_J[i],restMass_kg);
    integrand[i]=(speed>0.0) ? 4.0*kPi*T[i]*spectrumPerJ(E_J[i])/speed : 0.0;
  }
  return Trapezoid(E_J,integrand);
}

struct TransmissionDiagnostics {
  double RcLower_GV{0.0};
  double RcEffective_GV{0.0};
  double RcUpper_GV{0.0};
  double PenumbraWidth_GV{0.0};
  double THigh{0.0};
};

inline TransmissionDiagnostics ComputeTransmissionDiagnostics(
    const std::vector<double>& E_MeV,const std::vector<double>& T,
    double absoluteCharge_C,double restMass_kg) {
  TransmissionDiagnostics d;
  if (E_MeV.size()!=T.size() || E_MeV.size()<2) return d;
  for (std::size_t i=0;i<T.size();++i)
    if (!std::isfinite(T[i])) return d;
  std::vector<double> R(E_MeV.size(),0.0);
  for (std::size_t i=0;i<E_MeV.size();++i)
    R[i]=RigidityFromEnergyGV(E_MeV[i]*kMeVToJ,absoluteCharge_C,restMass_kg);
  d.THigh=std::max(0.0,T.back());
  const double tiny=1.0e-12, tolerance=1.0e-3;
  if (!(d.THigh>tiny)) return d;
  int first=-1,lastBelow=-1;
  for (int i=0;i<static_cast<int>(T.size());++i) {
    if (first<0 && T[static_cast<std::size_t>(i)]>tiny*d.THigh) first=i;
    if (T[static_cast<std::size_t>(i)]<(1.0-tolerance)*d.THigh) lastBelow=i;
  }
  if (first>=0) d.RcLower_GV=R[static_cast<std::size_t>(first)];
  d.RcUpper_GV=(lastBelow<0) ? R.front() :
      ((lastBelow+1<static_cast<int>(R.size())) ? R[static_cast<std::size_t>(lastBelow+1)] : R.back());
  double blockedArea=0.0;
  for (std::size_t i=0;i+1<R.size();++i) {
    const double f0=1.0-std::max(0.0,std::min(1.0,T[i]/d.THigh));
    const double f1=1.0-std::max(0.0,std::min(1.0,T[i+1]/d.THigh));
    blockedArea+=0.5*(f0+f1)*(R[i+1]-R[i]);
  }
  d.RcEffective_GV=std::max(R.front(),std::min(R.back(),R.front()+blockedArea));
  d.PenumbraWidth_GV=std::max(0.0,d.RcUpper_GV-d.RcLower_GV);
  return d;
}

} // namespace FluxNumerics
} // namespace Earth

#endif

#ifndef _SRC_EARTH_UTIL_BOUNDARY_PRODUCTS_H_
#define _SRC_EARTH_UTIL_BOUNDARY_PRODUCTS_H_

//======================================================================================
// BoundaryProducts.h -- Roadmap Step 6
//======================================================================================
// Field-backend-independent boundary-distribution and product kernels.
//
// This header is deliberately free of AMPS, MPI, SPICE, Geopack, and SWMF types.  A
// trajectory backend supplies access estimates A(E,Omega); a spectrum supplies the
// boundary differential intensity; this module performs the remaining physics and
// quadrature exactly once for gridless, standalone Mode3D, and coupled Mode3D callers.
//
// Unit contract
// -------------
//   energy coordinate in public tables             MeV or MeV/nucleon (declared)
//   kinetic energy passed to spectrum callables     J per particle
//   differential intensity passed to integrators    m^-2 s^-1 sr^-1 J^-1
//   differential spectra stored in ProductSet       m^-2 s^-1 sr^-1 MeV^-1
//   omnidirectional integral flux                    m^-2 s^-1
//   one-way planar integral flux                     m^-2 s^-1
//   number density                                   m^-3
//   detector response                                effective area*solid-angle [m^2 sr]
//   detector folded rate                            s^-1
//
// Static magnetic snapshots conserve kinetic energy, so J_local=A*J_boundary.  The
// GeneralPhaseSpace mapping is provided for a later validated electric/time-dependent
// characteristic and applies Liouville's j/p^2 invariant:
//
//   J_local = A * (p_local^2/p_boundary^2) * J_boundary.
//
// Production Step-6 callers select StaticMagnetic.  Merely enabling an electric field
// must not silently select the general mode; the trajectory contract continues to fail
// fast until the physical backward electromagnetic mover is released.
//======================================================================================

#include "FluxNumerics.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace Earth {
namespace BoundaryProducts {

using Earth::FluxNumerics::kMeVToJ;
using Earth::FluxNumerics::kPi;

enum class EnergyBasis { PerParticle, PerNucleon };

struct SpectrumUnits {
  EnergyBasis energyBasis{EnergyBasis::PerParticle};
  // Number of nucleons represented by one physical particle.  It must be one for a
  // per-particle proton/electron coordinate and positive for a per-nucleon ion grid.
  double massNumber{1.0};
  std::string intensityLabel{"m^-2 s^-1 sr^-1 MeV^-1"};

  double ParticleEnergyJ(double coordinateEnergy_MeV) const {
    if (!(coordinateEnergy_MeV >= 0.0) || !(massNumber > 0.0))
      throw std::invalid_argument("SpectrumUnits: invalid energy or mass number");
    const double factor=(energyBasis==EnergyBasis::PerNucleon) ? massNumber : 1.0;
    return coordinateEnergy_MeV*factor*kMeVToJ;
  }

  // Convert dJ/d(coordinate MeV) to dJ/d(particle joule).  For a MeV/nucleon
  // coordinate, dE_particle=A*dE_n, hence the additional 1/A factor.
  double PerCoordinateMeVToPerParticleJ(double intensity) const {
    if (!(massNumber > 0.0))
      throw std::invalid_argument("SpectrumUnits: mass number must be positive");
    const double factor=(energyBasis==EnergyBasis::PerNucleon) ? massNumber : 1.0;
    return intensity/(factor*kMeVToJ);
  }

  double PerParticleJToPerCoordinateMeV(double intensity) const {
    if (!(massNumber > 0.0))
      throw std::invalid_argument("SpectrumUnits: mass number must be positive");
    const double factor=(energyBasis==EnergyBasis::PerNucleon) ? massNumber : 1.0;
    return intensity*factor*kMeVToJ;
  }
};

struct Bounds {
  double nominal{0.0};
  double lower{0.0};
  double upper{0.0};

  Bounds() = default;
  Bounds(double n,double lo,double hi) : nominal(n),lower(lo),upper(hi) {}
};

inline Bounds OrderedBounds(double nominal,double lower,double upper) {
  Bounds r;
  r.lower=std::max(0.0,std::min(lower,upper));
  r.upper=std::max(r.lower,std::max(lower,upper));
  r.nominal=std::isfinite(nominal)
      ? std::max(r.lower,std::min(nominal,r.upper)) : nominal;
  return r;
}

inline Bounds MultiplyNonNegative(const Bounds& a,const Bounds& b) {
  return OrderedBounds(a.nominal*b.nominal,a.lower*b.lower,a.upper*b.upper);
}

enum class CharacteristicMapping { StaticMagnetic, GeneralPhaseSpace };

inline Bounds MapBoundaryIntensity(const Bounds& access,
                                   const Bounds& boundaryIntensity,
                                   CharacteristicMapping mapping,
                                   double localMomentum_kg_m_s=1.0,
                                   double boundaryMomentum_kg_m_s=1.0) {
  double ratio=1.0;
  if (mapping==CharacteristicMapping::GeneralPhaseSpace) {
    if (!(localMomentum_kg_m_s > 0.0) || !(boundaryMomentum_kg_m_s > 0.0) ||
        !std::isfinite(localMomentum_kg_m_s) ||
        !std::isfinite(boundaryMomentum_kg_m_s)) {
      throw std::invalid_argument(
          "GeneralPhaseSpace mapping requires finite positive local/boundary momenta");
    }
    const double q=localMomentum_kg_m_s/boundaryMomentum_kg_m_s;
    ratio=q*q;
  }
  const Bounds mapped=MultiplyNonNegative(access,boundaryIntensity);
  return Bounds(mapped.nominal*ratio,mapped.lower*ratio,mapped.upper*ratio);
}

//--------------------------------------------------------------------------------------
// PAD and spatial factors
//--------------------------------------------------------------------------------------
enum class NormalizationMode { Raw, UnitMean };
enum class PadModel { Isotropic, SinAlphaN, CosAlphaN, Bidirectional };
enum class SpatialModel { Uniform, DaysideNightside };

inline double PadRawMean(PadModel model,double exponent) {
  if (!(exponent >= 0.0) || !std::isfinite(exponent))
    throw std::invalid_argument("PAD exponent must be finite and non-negative");
  if (model==PadModel::Isotropic) return 1.0;
  if (model==PadModel::CosAlphaN || model==PadModel::Bidirectional)
    return 1.0/(exponent+1.0);
  // <sin(alpha)^n> over 4pi = 1/2 integral[-1,1](1-mu^2)^(n/2)dmu.
  return std::sqrt(kPi)*std::tgamma(0.5*exponent+1.0) /
         (2.0*std::tgamma(0.5*exponent+1.5));
}

inline double PadWeight(PadModel model,double cosAlpha,double exponent,
                        NormalizationMode normalization) {
  if (!(exponent >= 0.0) || !std::isfinite(exponent) || !std::isfinite(cosAlpha))
    throw std::invalid_argument("PAD input must be finite and exponent non-negative");
  const double mu=std::max(-1.0,std::min(1.0,cosAlpha));
  double raw=1.0;
  if (model==PadModel::SinAlphaN)
    raw=std::pow(std::max(0.0,1.0-mu*mu),0.5*exponent);
  else if (model==PadModel::CosAlphaN || model==PadModel::Bidirectional)
    raw=std::pow(std::fabs(mu),exponent);
  if (normalization==NormalizationMode::UnitMean) raw/=PadRawMean(model,exponent);
  return raw;
}

inline double SpatialWeight(SpatialModel model,double gsmX_m,
                            double daysideFactor,double nightsideFactor,
                            NormalizationMode normalization) {
  if (!(daysideFactor >= 0.0) || !(nightsideFactor >= 0.0) ||
      !std::isfinite(daysideFactor) || !std::isfinite(nightsideFactor))
    throw std::invalid_argument("Spatial factors must be finite and non-negative");
  if (model==SpatialModel::Uniform) return 1.0;
  double value=(gsmX_m>0.0) ? daysideFactor : nightsideFactor;
  if (normalization==NormalizationMode::UnitMean) {
    // The reference boundary measure gives the two x hemispheres equal solid angle.
    const double mean=0.5*(daysideFactor+nightsideFactor);
    if (!(mean>0.0))
      throw std::invalid_argument("Cannot normalize zero dayside+nightside distribution");
    value/=mean;
  }
  return value;
}

inline double MaximumPadWeight(PadModel model,double exponent,
                               NormalizationMode normalization) {
  const double rawMax=1.0;
  return normalization==NormalizationMode::Raw
      ? rawMax : rawMax/PadRawMean(model,exponent);
}

inline double MaximumSpatialWeight(SpatialModel model,double daysideFactor,
                                   double nightsideFactor,
                                   NormalizationMode normalization) {
  if (model==SpatialModel::Uniform) return 1.0;
  const double rawMax=std::max(daysideFactor,nightsideFactor);
  if (normalization==NormalizationMode::Raw) return rawMax;
  const double mean=0.5*(daysideFactor+nightsideFactor);
  if (!(mean>0.0))
    throw std::invalid_argument("Cannot normalize zero dayside+nightside distribution");
  return rawMax/mean;
}

//--------------------------------------------------------------------------------------
// Deterministic log-intensity interpolation for time-dependent spectrum tables
//--------------------------------------------------------------------------------------
enum class OutOfRangePolicy { Clamp, Zero, Fail };
enum class GapPolicy { InterpolateAndFlag, HoldNearest, Fail };
enum class TemporalStatus {
  Exact,
  Interpolated,
  GapInterpolated,
  GapHeld,
  ClampedBefore,
  ClampedAfter,
  ZeroBefore,
  ZeroAfter
};

struct TemporalSpectrumRow {
  double time_s{0.0};
  std::string epochUTC;
  std::vector<double> intensity;
};

struct TemporalSpectrumSelection {
  std::vector<double> intensity;
  TemporalStatus status{TemporalStatus::Exact};
  std::size_t leftIndex{0};
  std::size_t rightIndex{0};
  double interpolationFraction{0.0};
  bool gapFlag{false};
};

inline const char* TemporalStatusName(TemporalStatus status) {
  switch (status) {
    case TemporalStatus::Exact: return "EXACT";
    case TemporalStatus::Interpolated: return "INTERPOLATED";
    case TemporalStatus::GapInterpolated: return "GAP_INTERPOLATED";
    case TemporalStatus::GapHeld: return "GAP_HELD";
    case TemporalStatus::ClampedBefore: return "CLAMPED_BEFORE";
    case TemporalStatus::ClampedAfter: return "CLAMPED_AFTER";
    case TemporalStatus::ZeroBefore: return "ZERO_BEFORE";
    case TemporalStatus::ZeroAfter: return "ZERO_AFTER";
  }
  return "UNKNOWN";
}

inline TemporalSpectrumSelection SelectTemporalSpectrum(
    const std::vector<TemporalSpectrumRow>& rows,double evaluationTime_s,
    double maximumGap_s,OutOfRangePolicy outOfRangePolicy,
    GapPolicy gapPolicy) {
  if (rows.empty() || !std::isfinite(evaluationTime_s))
    throw std::invalid_argument("SelectTemporalSpectrum: empty rows or invalid time");
  const std::size_t nValue=rows.front().intensity.size();
  if (nValue==0) throw std::invalid_argument("SelectTemporalSpectrum: empty intensity row");
  for (std::size_t i=0;i<rows.size();++i) {
    if (!std::isfinite(rows[i].time_s) || rows[i].intensity.size()!=nValue ||
        (i>0 && !(rows[i].time_s>rows[i-1].time_s)))
      throw std::invalid_argument("SelectTemporalSpectrum: rows must be ordered and rectangular");
    for (double v:rows[i].intensity)
      if (!(v>0.0) || !std::isfinite(v))
        throw std::invalid_argument("SelectTemporalSpectrum: log interpolation requires positive finite intensity");
  }

  TemporalSpectrumSelection out;
  const auto useEndpoint=[&](std::size_t index,TemporalStatus status) {
    out.intensity=rows[index].intensity;
    out.status=status;
    out.leftIndex=out.rightIndex=index;
    out.interpolationFraction=0.0;
  };
  const auto useZero=[&](TemporalStatus status) {
    out.intensity.assign(nValue,0.0);
    out.status=status;
    out.leftIndex=out.rightIndex=(status==TemporalStatus::ZeroBefore ? 0 : rows.size()-1);
  };

  if (evaluationTime_s < rows.front().time_s) {
    if (outOfRangePolicy==OutOfRangePolicy::Fail)
      throw std::out_of_range("Spectrum evaluation time precedes table");
    if (outOfRangePolicy==OutOfRangePolicy::Zero) useZero(TemporalStatus::ZeroBefore);
    else useEndpoint(0,TemporalStatus::ClampedBefore);
    return out;
  }
  if (evaluationTime_s > rows.back().time_s) {
    if (outOfRangePolicy==OutOfRangePolicy::Fail)
      throw std::out_of_range("Spectrum evaluation time follows table");
    if (outOfRangePolicy==OutOfRangePolicy::Zero) useZero(TemporalStatus::ZeroAfter);
    else useEndpoint(rows.size()-1,TemporalStatus::ClampedAfter);
    return out;
  }

  const auto it=std::lower_bound(rows.begin(),rows.end(),evaluationTime_s,
      [](const TemporalSpectrumRow& row,double t) { return row.time_s<t; });
  const std::size_t right=static_cast<std::size_t>(it-rows.begin());
  if (right<rows.size() && rows[right].time_s==evaluationTime_s) {
    useEndpoint(right,TemporalStatus::Exact);
    return out;
  }
  const std::size_t left=right-1;
  const double dt=rows[right].time_s-rows[left].time_s;
  const bool gap=(maximumGap_s>0.0 && dt>maximumGap_s);
  out.leftIndex=left;
  out.rightIndex=right;
  out.gapFlag=gap;

  if (gap && gapPolicy==GapPolicy::Fail)
    throw std::runtime_error("Spectrum interpolation interval exceeds maximum gap");
  if (gap && gapPolicy==GapPolicy::HoldNearest) {
    const std::size_t nearest=(evaluationTime_s-rows[left].time_s <=
                               rows[right].time_s-evaluationTime_s) ? left : right;
    useEndpoint(nearest,TemporalStatus::GapHeld);
    out.gapFlag=true;
    return out;
  }

  const double fraction=(evaluationTime_s-rows[left].time_s)/dt;
  out.interpolationFraction=fraction;
  out.status=gap ? TemporalStatus::GapInterpolated : TemporalStatus::Interpolated;
  out.intensity.resize(nValue);
  for (std::size_t i=0;i<nValue;++i) {
    const double logValue=(1.0-fraction)*std::log(rows[left].intensity[i])+
                          fraction*std::log(rows[right].intensity[i]);
    out.intensity[i]=std::exp(logValue);
  }
  return out;
}

//--------------------------------------------------------------------------------------
// Directional differential products at one energy
//--------------------------------------------------------------------------------------
struct DirectionalAccessSample {
  double x{0.0},y{0.0},z{1.0}; // arrival velocity unit vector at the observer
  double solidAngleWeight_sr{0.0};
  Bounds access{0.0,0.0,0.0};
  double boundaryFactor{1.0};  // PAD*spatial multiplier at the boundary exit
  double localMomentum_kg_m_s{1.0};
  double boundaryMomentum_kg_m_s{1.0};
};

struct DirectionalDifferentialProduct {
  Bounds angularMeanPerMeV;
  Bounds omnidirectionalPerMeV;
  Bounds oneWayPlanarPerMeV;
};

inline DirectionalDifferentialProduct FoldDirectionalDifferential(
    const std::vector<DirectionalAccessSample>& samples,const Bounds& boundaryPerMeV,
    CharacteristicMapping mapping,double detectorNormalX=0.0,
    double detectorNormalY=0.0,double detectorNormalZ=1.0) {
  DirectionalDifferentialProduct out;
  const double nn=std::sqrt(detectorNormalX*detectorNormalX+
                            detectorNormalY*detectorNormalY+
                            detectorNormalZ*detectorNormalZ);
  if (!(nn>0.0)) throw std::invalid_argument("Detector normal must be non-zero");
  const double nx=detectorNormalX/nn,ny=detectorNormalY/nn,nz=detectorNormalZ/nn;
  double totalWeight=0.0;
  for (const DirectionalAccessSample& s:samples) {
    if (!(s.solidAngleWeight_sr>=0.0) || !(s.boundaryFactor>=0.0))
      throw std::invalid_argument("Directional weights and boundary factors must be non-negative");
    Bounds b(boundaryPerMeV.nominal*s.boundaryFactor,
             boundaryPerMeV.lower*s.boundaryFactor,
             boundaryPerMeV.upper*s.boundaryFactor);
    const Bounds local=MapBoundaryIntensity(s.access,b,mapping,
        s.localMomentum_kg_m_s,s.boundaryMomentum_kg_m_s);
    const double w=s.solidAngleWeight_sr;
    out.omnidirectionalPerMeV.nominal+=w*local.nominal;
    out.omnidirectionalPerMeV.lower+=w*local.lower;
    out.omnidirectionalPerMeV.upper+=w*local.upper;
    // Arrival velocity into a surface with outward normal n has Omega.n < 0.
    const double projected=std::max(0.0,-(s.x*nx+s.y*ny+s.z*nz));
    out.oneWayPlanarPerMeV.nominal+=w*projected*local.nominal;
    out.oneWayPlanarPerMeV.lower+=w*projected*local.lower;
    out.oneWayPlanarPerMeV.upper+=w*projected*local.upper;
    totalWeight+=w;
  }
  if (!(totalWeight>0.0))
    throw std::invalid_argument("Directional product requires positive solid-angle support");
  out.angularMeanPerMeV.nominal=out.omnidirectionalPerMeV.nominal/totalWeight;
  out.angularMeanPerMeV.lower=out.omnidirectionalPerMeV.lower/totalWeight;
  out.angularMeanPerMeV.upper=out.omnidirectionalPerMeV.upper/totalWeight;
  return out;
}

//--------------------------------------------------------------------------------------
// Energy-integrated products for a direction-averaged access curve
//--------------------------------------------------------------------------------------
struct EnergyChannel {
  std::string name;
  double lower_MeV{0.0};
  double upper_MeV{0.0};
  EnergyChannel() = default;
  EnergyChannel(const std::string& nameIn,double lowerIn,double upperIn)
      : name(nameIn),lower_MeV(lowerIn),upper_MeV(upperIn) {}
};

struct DetectorResponse {
  std::string name;
  std::vector<double> energy_MeV;
  std::vector<double> relativeResponse;
  double geometricFactor_m2_sr{0.0};

  void Validate() const {
    if (name.empty() || energy_MeV.size()<2 ||
        energy_MeV.size()!=relativeResponse.size() ||
        !(geometricFactor_m2_sr>=0.0) || !std::isfinite(geometricFactor_m2_sr))
      throw std::invalid_argument("Invalid detector response definition");
    for (std::size_t i=0;i<energy_MeV.size();++i) {
      if (!(energy_MeV[i]>0.0) || !std::isfinite(energy_MeV[i]) ||
          !(relativeResponse[i]>=0.0) || !std::isfinite(relativeResponse[i]) ||
          (i>0 && !(energy_MeV[i]>energy_MeV[i-1])))
        throw std::invalid_argument("Detector response nodes must be finite, non-negative, and ordered");
    }
  }

  double Evaluate(double energy) const {
    if (energy<energy_MeV.front() || energy>energy_MeV.back()) return 0.0;
    const auto it=std::lower_bound(energy_MeV.begin(),energy_MeV.end(),energy);
    if (it==energy_MeV.begin()) return relativeResponse.front();
    if (it==energy_MeV.end()) return relativeResponse.back();
    const std::size_t i1=static_cast<std::size_t>(it-energy_MeV.begin());
    if (*it==energy) return relativeResponse[i1];
    const std::size_t i0=i1-1;
    const double a=(energy-energy_MeV[i0])/(energy_MeV[i1]-energy_MeV[i0]);
    return (1.0-a)*relativeResponse[i0]+a*relativeResponse[i1];
  }
};

struct SpectrumPoint {
  double energy_MeV{0.0};
  Bounds boundaryPerMeV;
  Bounds localPerMeV;
  Bounds omnidirectionalPerMeV;
  Bounds oneWayPlanarPerMeV;
};

struct ProductSet {
  std::vector<SpectrumPoint> spectrum;
  Bounds numberDensity_m3;
  Bounds omnidirectionalFlux_m2_s;
  Bounds oneWayPlanarFlux_m2_s;
  std::vector<Bounds> channelFlux_m2_s;
  std::vector<Bounds> detectorRate_s;
};

inline double InterpolateLinear(const std::vector<double>& x,const std::vector<double>& y,
                                double q) {
  if (x.size()!=y.size() || x.empty())
    throw std::invalid_argument("InterpolateLinear: incompatible arrays");
  if (q<=x.front()) return y.front();
  if (q>=x.back()) return y.back();
  const auto it=std::upper_bound(x.begin(),x.end(),q);
  const std::size_t i1=static_cast<std::size_t>(it-x.begin()),i0=i1-1;
  const double a=(q-x[i0])/(x[i1]-x[i0]);
  return (1.0-a)*y[i0]+a*y[i1];
}

template<class SpectrumEvaluator>
inline Bounds IntegrateDetectorRate(const std::vector<double>& energy_MeV,
                                    const std::vector<double>& nominalAccess,
                                    const std::vector<double>& lowerAccess,
                                    const std::vector<double>& upperAccess,
                                    const DetectorResponse& response,
                                    SpectrumEvaluator boundarySpectrumPerJ) {
  response.Validate();
  if (energy_MeV.size()<2 || energy_MeV.size()!=nominalAccess.size() ||
      energy_MeV.size()!=lowerAccess.size() || energy_MeV.size()!=upperAccess.size())
    throw std::invalid_argument("IntegrateDetectorRate: incompatible access arrays");
  const double lo=std::max(energy_MeV.front(),response.energy_MeV.front());
  const double hi=std::min(energy_MeV.back(),response.energy_MeV.back());
  if (!(hi>lo)) return Bounds();

  // Add both solver and response nodes.  This is essential for top-hat/triangular
  // channel edges: evaluating only the solver grid can miss a narrow response.
  std::vector<double> nodes{lo,hi};
  for (double e:energy_MeV) if (e>lo && e<hi) nodes.push_back(e);
  for (double e:response.energy_MeV) if (e>lo && e<hi) nodes.push_back(e);
  std::sort(nodes.begin(),nodes.end());
  nodes.erase(std::unique(nodes.begin(),nodes.end()),nodes.end());

  Bounds out;
  for (std::size_t i=0;i+1<nodes.size();++i) {
    const double e0=nodes[i],e1=nodes[i+1];
    const double j0=boundarySpectrumPerJ(e0*kMeVToJ);
    const double j1=boundarySpectrumPerJ(e1*kMeVToJ);
    const double r0=response.Evaluate(e0),r1=response.Evaluate(e1);
    const double dE=(e1-e0)*kMeVToJ;
    const double tn0=InterpolateLinear(energy_MeV,nominalAccess,e0);
    const double tn1=InterpolateLinear(energy_MeV,nominalAccess,e1);
    const double tl0=InterpolateLinear(energy_MeV,lowerAccess,e0);
    const double tl1=InterpolateLinear(energy_MeV,lowerAccess,e1);
    const double tu0=InterpolateLinear(energy_MeV,upperAccess,e0);
    const double tu1=InterpolateLinear(energy_MeV,upperAccess,e1);
    const double scale=0.5*dE*response.geometricFactor_m2_sr;
    out.nominal+=scale*(tn0*j0*r0+tn1*j1*r1);
    out.lower+=scale*(tl0*j0*r0+tl1*j1*r1);
    out.upper+=scale*(tu0*j0*r0+tu1*j1*r1);
  }
  return out;
}

template<class SpectrumEvaluator>
inline double IntegrateDensityWithUnits(const std::vector<double>& energy_MeV,
                                        const std::vector<double>& access,
                                        double particleMass_kg,
                                        SpectrumEvaluator boundarySpectrumPerCoordinateJ,
                                        const SpectrumUnits& units) {
  if (energy_MeV.size()<2 || energy_MeV.size()!=access.size() ||
      !(particleMass_kg>0.0))
    throw std::invalid_argument("IntegrateDensityWithUnits: invalid input");
  double result=0.0;
  for (std::size_t i=0;i+1<energy_MeV.size();++i) {
    if (!(energy_MeV[i+1]>energy_MeV[i]))
      throw std::invalid_argument("IntegrateDensityWithUnits: energy grid is not increasing");
    const double e0CoordinateJ=energy_MeV[i]*kMeVToJ;
    const double e1CoordinateJ=energy_MeV[i+1]*kMeVToJ;
    const double v0=Earth::FluxNumerics::RelativisticSpeed(
        units.ParticleEnergyJ(energy_MeV[i]),particleMass_kg);
    const double v1=Earth::FluxNumerics::RelativisticSpeed(
        units.ParticleEnergyJ(energy_MeV[i+1]),particleMass_kg);
    if (!(v0>0.0) || !(v1>0.0)) continue;
    const double f0=4.0*kPi*access[i]*boundarySpectrumPerCoordinateJ(e0CoordinateJ)/v0;
    const double f1=4.0*kPi*access[i+1]*boundarySpectrumPerCoordinateJ(e1CoordinateJ)/v1;
    result+=0.5*(f0+f1)*(e1CoordinateJ-e0CoordinateJ);
  }
  return result;
}

template<class SpectrumEvaluator>
inline ProductSet EvaluateIsotropicProducts(
    const std::vector<double>& energy_MeV,
    const std::vector<double>& nominalAccess,
    const std::vector<double>& lowerAccess,
    const std::vector<double>& upperAccess,
    double particleMass_kg,SpectrumEvaluator boundarySpectrumPerJ,
    const std::vector<EnergyChannel>& channels=std::vector<EnergyChannel>(),
    const std::vector<DetectorResponse>& responses=std::vector<DetectorResponse>(),
    double boundaryRelativeUncertainty=0.0,
    const SpectrumUnits& units=SpectrumUnits()) {
  const std::size_t n=energy_MeV.size();
  if (n<2 || nominalAccess.size()!=n || lowerAccess.size()!=n || upperAccess.size()!=n ||
      !(particleMass_kg>0.0) || !(boundaryRelativeUncertainty>=0.0) ||
      !std::isfinite(boundaryRelativeUncertainty))
    throw std::invalid_argument("EvaluateIsotropicProducts: invalid input arrays or metadata");

  ProductSet out;
  const double spectrumLowerScale=std::max(0.0,1.0-boundaryRelativeUncertainty);
  const double spectrumUpperScale=1.0+boundaryRelativeUncertainty;
  out.spectrum.reserve(n);
  for (std::size_t i=0;i<n;++i) {
    const double boundaryJ=boundarySpectrumPerJ(energy_MeV[i]*kMeVToJ);
    const double boundaryPerMeV=boundaryJ*kMeVToJ;
    const Bounds boundary(boundaryPerMeV,
        boundaryPerMeV*std::max(0.0,1.0-boundaryRelativeUncertainty),
        boundaryPerMeV*(1.0+boundaryRelativeUncertainty));
    const Bounds access=OrderedBounds(nominalAccess[i],lowerAccess[i],upperAccess[i]);
    const Bounds local=MapBoundaryIntensity(access,boundary,
                                            CharacteristicMapping::StaticMagnetic);
    SpectrumPoint p;
    p.energy_MeV=energy_MeV[i];
    p.boundaryPerMeV=boundary;
    p.localPerMeV=local;
    p.omnidirectionalPerMeV=Bounds(4.0*kPi*local.nominal,
                                   4.0*kPi*local.lower,
                                   4.0*kPi*local.upper);
    p.oneWayPlanarPerMeV=Bounds(kPi*local.nominal,kPi*local.lower,kPi*local.upper);
    out.spectrum.push_back(p);
  }

  out.numberDensity_m3=Bounds(
      IntegrateDensityWithUnits(energy_MeV,nominalAccess,particleMass_kg,boundarySpectrumPerJ,units),
      spectrumLowerScale*IntegrateDensityWithUnits(energy_MeV,lowerAccess,particleMass_kg,boundarySpectrumPerJ,units),
      spectrumUpperScale*IntegrateDensityWithUnits(energy_MeV,upperAccess,particleMass_kg,boundarySpectrumPerJ,units));
  out.omnidirectionalFlux_m2_s=Bounds(
      Earth::FluxNumerics::IntegrateFlux(energy_MeV,nominalAccess,energy_MeV.front(),energy_MeV.back(),boundarySpectrumPerJ),
      spectrumLowerScale*Earth::FluxNumerics::IntegrateFlux(energy_MeV,lowerAccess,energy_MeV.front(),energy_MeV.back(),boundarySpectrumPerJ),
      spectrumUpperScale*Earth::FluxNumerics::IntegrateFlux(energy_MeV,upperAccess,energy_MeV.front(),energy_MeV.back(),boundarySpectrumPerJ));
  out.oneWayPlanarFlux_m2_s=Bounds(
      0.25*out.omnidirectionalFlux_m2_s.nominal,
      0.25*out.omnidirectionalFlux_m2_s.lower,
      0.25*out.omnidirectionalFlux_m2_s.upper);

  out.channelFlux_m2_s.reserve(channels.size());
  for (const EnergyChannel& channel:channels) {
    if (!(channel.lower_MeV>0.0) || !(channel.upper_MeV>channel.lower_MeV))
      throw std::invalid_argument("Energy channel requires 0 < lower < upper");
    out.channelFlux_m2_s.push_back(Bounds(
        Earth::FluxNumerics::IntegrateFlux(energy_MeV,nominalAccess,channel.lower_MeV,channel.upper_MeV,boundarySpectrumPerJ),
        spectrumLowerScale*Earth::FluxNumerics::IntegrateFlux(energy_MeV,lowerAccess,channel.lower_MeV,channel.upper_MeV,boundarySpectrumPerJ),
        spectrumUpperScale*Earth::FluxNumerics::IntegrateFlux(energy_MeV,upperAccess,channel.lower_MeV,channel.upper_MeV,boundarySpectrumPerJ)));
  }
  out.detectorRate_s.reserve(responses.size());
  for (const DetectorResponse& response:responses) {
    Bounds rate=IntegrateDetectorRate(
        energy_MeV,nominalAccess,lowerAccess,upperAccess,response,boundarySpectrumPerJ);
    rate.lower*=spectrumLowerScale;
    rate.upper*=spectrumUpperScale;
    out.detectorRate_s.push_back(rate);
  }
  return out;
}

} // namespace BoundaryProducts
} // namespace Earth

#endif

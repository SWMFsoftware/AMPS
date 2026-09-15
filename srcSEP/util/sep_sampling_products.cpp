#include "sep_sampling_products.h"

#include <algorithm>
#include <cmath>

namespace SEP {
namespace SamplingCore {
namespace {

Transport::Status Error(const std::string& text) {
  return Transport::Status::Error(Transport::StatusCode::InvalidArgument,text);
}

bool ValidGrid(const BinGrid& grid) {
  if (grid.edges.size()<2) return false;
  for (std::size_t i=0;i<grid.edges.size();++i)
    if (!std::isfinite(grid.edges[i]) ||
        (i && grid.edges[i]<=grid.edges[i-1])) return false;
  return true;
}

}  // namespace

Kinematics EvaluateKinematics(double vp,double vt,double mass,double charge,
                              double magneticField,double c) {
  Kinematics result;
  if (!std::isfinite(vp) || !std::isfinite(vt) || !std::isfinite(mass) ||
      !std::isfinite(charge) || !std::isfinite(magneticField) ||
      !std::isfinite(c) || mass<=0.0 || charge==0.0 || magneticField<=0.0 || c<=0.0) {
    result.status=Error("sampling kinematics require finite velocity and physical species/background");
    return result;
  }
  result.speedMPerS=std::hypot(vp,vt);
  if (!(result.speedMPerS>0.0) || result.speedMPerS>=c) {
    result.status=Transport::Status::Error(Transport::StatusCode::InvalidParticleState,
        "sampling excludes zero-speed and luminal/superluminal particles");
    return result;
  }
  result.mu=vp/result.speedMPerS;
  const double beta=result.speedMPerS/c;
  const double gamma=1.0/std::sqrt(1.0-beta*beta);
  result.momentumKgMPerS=gamma*mass*result.speedMPerS;
  result.kineticEnergyJ=(gamma-1.0)*mass*c*c;
  const double perpendicularMomentum=gamma*mass*std::fabs(vt);
  result.larmorRadiusM=perpendicularMomentum/(std::fabs(charge)*magneticField);
  result.status=std::isfinite(result.larmorRadiusM) &&
      std::isfinite(result.kineticEnergyJ) ? Transport::Status::Ok()
      : Error("sampling kinematic conversion overflowed");
  return result;
}

BinResult Classify(const BinGrid& grid,double value) {
  BinResult result;
  if (!ValidGrid(grid) || !std::isfinite(value)) return result;
  if (value<grid.edges.front()) {
    result.classification=BinClassification::Underflow; return result;
  }
  if (value>=grid.edges.back()) {
    result.classification=BinClassification::Overflow; return result;
  }
  result.classification=BinClassification::InRange;
  result.index=static_cast<std::size_t>(
      std::upper_bound(grid.edges.begin(),grid.edges.end(),value)-
      grid.edges.begin()-1);
  return result;
}

Transport::Status Initialize(const BinGrid& grid,WeightedBins* bins) {
  if (!bins || !ValidGrid(grid)) return Error("sampling grid is invalid");
  bins->grid=grid;
  bins->sumWeight.assign(grid.edges.size()-1,0.0);
  bins->sumWeightSquared.assign(grid.edges.size()-1,0.0);
  bins->counters=InvalidCounters();
  return Transport::Status::Ok();
}

Transport::Status Accumulate(WeightedBins* bins,double value,double weight) {
  if (!bins || bins->sumWeight.size()+1!=bins->grid.edges.size() ||
      !std::isfinite(weight) || weight<0.0) return Error("sampling accumulation is invalid");
  const BinResult bin=Classify(bins->grid,value);
  if (bin.classification==BinClassification::Underflow) ++bins->counters.underflow;
  else if (bin.classification==BinClassification::Overflow) ++bins->counters.overflow;
  else if (bin.classification==BinClassification::Invalid) ++bins->counters.invalidKinematics;
  else {
    bins->sumWeight[bin.index]+=weight;
    bins->sumWeightSquared[bin.index]+=weight*weight;
    ++bins->counters.accepted;
  }
  return Transport::Status::Ok();
}

double EffectiveSampleSize(const WeightedBins& bins,std::size_t i) {
  if (i>=bins.sumWeight.size() || bins.sumWeightSquared[i]<=0.0) return 0.0;
  return bins.sumWeight[i]*bins.sumWeight[i]/bins.sumWeightSquared[i];
}

double StandardErrorOfWeightedCount(const WeightedBins& bins,std::size_t i) {
  if (i>=bins.sumWeightSquared.size()) return 0.0;
  // For independent weighted events, sum(w_i^2) is the Poisson variance
  // estimator of the weighted count.  Product-specific geometry/bin-width
  // factors are applied to both count and error by the output adapter.
  return std::sqrt(std::max(0.0,bins.sumWeightSquared[i]));
}

ProductValues BuildProduct(const WeightedBins& bins,Product product,
                           const NormalizationContext& context) {
  ProductValues result;
  if (!ValidGrid(bins.grid) ||
      bins.sumWeight.size()+1!=bins.grid.edges.size() ||
      bins.sumWeightSquared.size()!=bins.sumWeight.size()) {
    result.status=Error("sampling product has inconsistent bins");
    return result;
  }
  double denominatorScale=0.0;
  if (product==Product::NumberDensityPerEnergy ||
      product==Product::OmnidirectionalDifferentialIntensity) {
    if (!(context.volumeM3>0.0) || !std::isfinite(context.volumeM3)) {
      result.status=Error("density/intensity product requires positive volume");
      return result;
    }
    denominatorScale=context.volumeM3;
    result.units=product==Product::NumberDensityPerEnergy ? "m^-3 J^-1"
                                                          : "m^-2 s^-1 sr^-1 J^-1";
  }
  else if (product==Product::DirectionalCrossingFlux) {
    if (!(context.areaM2>0.0) || !(context.durationS>0.0) ||
        !std::isfinite(context.areaM2) || !std::isfinite(context.durationS)) {
      result.status=Error("crossing flux requires positive area and duration");
      return result;
    }
    denominatorScale=context.areaM2*context.durationS;
    result.units="m^-2 s^-1 J^-1";
  }
  else {
    long double total=0.0L;
    for (std::size_t i=0;i<bins.sumWeight.size();++i) total+=bins.sumWeight[i];
    if (!(total>0.0L)) {
      result.status=Error("normalized shape requires positive sampled weight");
      return result;
    }
    denominatorScale=static_cast<double>(total);
    result.units="J^-1";
  }
  if (product==Product::OmnidirectionalDifferentialIntensity &&
      (!(context.representativeSpeedMPerS>0.0) ||
       !std::isfinite(context.representativeSpeedMPerS))) {
    result.status=Error("intensity requires positive representative speed");
    return result;
  }
  result.value.resize(bins.sumWeight.size());
  result.standardError.resize(bins.sumWeight.size());
  const double intensityFactor=product==Product::OmnidirectionalDifferentialIntensity
      ? context.representativeSpeedMPerS/(4.0*3.14159265358979323846) : 1.0;
  for (std::size_t i=0;i<bins.sumWeight.size();++i) {
    const double width=bins.grid.edges[i+1]-bins.grid.edges[i];
    result.value[i]=intensityFactor*bins.sumWeight[i]/(denominatorScale*width);
    result.standardError[i]=intensityFactor*
        StandardErrorOfWeightedCount(bins,i)/(denominatorScale*width);
  }
  result.status=Transport::Status::Ok();
  return result;
}

std::string ProductName(Product p) {
  switch (p) {
    case Product::NumberDensityPerEnergy: return "number-density-per-energy";
    case Product::OmnidirectionalDifferentialIntensity:
      return "omnidirectional-differential-intensity";
    case Product::DirectionalCrossingFlux: return "directional-crossing-flux";
    case Product::NormalizedShape: return "normalized-shape";
  }
  return "unknown";
}

Transport::Status ValidateMetadata(const ProductMetadata& m) {
  if (m.schemaVersion!="srcsep-sampling-v2" || m.units.empty() ||
      m.geometryFingerprint.empty() || m.backgroundFingerprint.empty() ||
      m.coefficientFingerprint.empty() || m.runFingerprint.empty() ||
      m.species.empty() || !std::isfinite(m.snapshotEpochS))
    return Error("sampling product metadata is incomplete");
  return Transport::Status::Ok();
}

}  // namespace SamplingCore
}  // namespace SEP

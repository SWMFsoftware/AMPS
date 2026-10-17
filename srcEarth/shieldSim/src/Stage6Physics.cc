#include "Stage6Physics.hh"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace shield_sim {

LogarithmicLETBinning::LogarithmicLETBinning(
    double minimum,double maximum,int binsPerDecade){
  if(!std::isfinite(minimum) || !std::isfinite(maximum) || minimum<=0.0 ||
     maximum<=minimum || binsPerDecade<=0 || binsPerDecade>1000)
    throw std::invalid_argument("invalid logarithmic LET binning");
  const double decades=std::log10(maximum/minimum);
  const double exactCount=decades*binsPerDecade;
  const long count=std::lround(exactCount);
  if(count<=0 || count>100000 ||
     std::abs(exactCount-count)>1.0e-10*std::max(1.0,std::abs(exactCount)))
    throw std::invalid_argument(
        "LET range must span an integer number of bins at bins_per_decade");
  edges_.reserve(static_cast<std::size_t>(count)+1);
  for(long index=0;index<=count;++index)
    edges_.push_back(minimum*std::pow(10.0,
        static_cast<double>(index)/binsPerDecade));
  // Preserve the user's exact endpoints against accumulated pow/log roundoff.
  edges_.front()=minimum;
  edges_.back()=maximum;
}

long LogarithmicLETBinning::Locate(double value) const {
  if(edges_.size()<2 || !std::isfinite(value) || value<edges_.front()) return -1;
  if(value>edges_.back()) return static_cast<long>(BinCount());
  if(value==edges_.back()) return static_cast<long>(BinCount()-1);
  const auto upper=std::upper_bound(edges_.begin(),edges_.end(),value);
  return static_cast<long>(upper-edges_.begin()-1);
}

double LogarithmicLETBinning::Lower(std::size_t index) const {
  if(index>=BinCount()) throw std::out_of_range("LET bin index");
  return edges_[index];
}
double LogarithmicLETBinning::Upper(std::size_t index) const {
  if(index>=BinCount()) throw std::out_of_range("LET bin index");
  return edges_[index+1];
}
double LogarithmicLETBinning::Width(std::size_t index) const {
  return Upper(index)-Lower(index);
}
double LogarithmicLETBinning::GeometricCenter(std::size_t index) const {
  return std::sqrt(Lower(index)*Upper(index));
}

double IntegratePowerLawAbove(const PowerLawFluenceBin& bin,double threshold){
  if(!std::isfinite(bin.lowerMeV) || !std::isfinite(bin.upperMeV) ||
     !std::isfinite(bin.referenceEnergyMeV) ||
     !std::isfinite(bin.valueAtReferencePerCm2PerMeV) ||
     !std::isfinite(bin.spectralIndex) || bin.lowerMeV<=0.0 ||
     bin.upperMeV<=bin.lowerMeV || bin.referenceEnergyMeV<=0.0 ||
     bin.valueAtReferencePerCm2PerMeV<0.0 || !std::isfinite(threshold))
    throw std::invalid_argument("invalid power-law fluence bin or threshold");
  if(bin.valueAtReferencePerCm2PerMeV==0.0 || threshold>=bin.upperMeV)
    return 0.0;
  const double lower=std::max(bin.lowerMeV,threshold);
  const double upper=bin.upperMeV;
  // Integrate C*(E/Eref)^p exactly, including the p=-1 logarithmic limit.
  if(std::abs(bin.spectralIndex+1.0)<1.0e-12)
    return bin.valueAtReferencePerCm2PerMeV*bin.referenceEnergyMeV*
        std::log(upper/lower);
  const double power=bin.spectralIndex+1.0;
  return bin.valueAtReferencePerCm2PerMeV*bin.referenceEnergyMeV/power*
      (std::pow(upper/bin.referenceEnergyMeV,power)-
       std::pow(lower/bin.referenceEnergyMeV,power));
}

HardnessFoldResult FoldH100Over10(const std::vector<PowerLawFluenceBin>& bins){
  HardnessFoldResult result;
  for(const auto& bin:bins){
    result.denominatorFluencePerCm2+=IntegratePowerLawAbove(bin,10.0);
    result.numeratorFluencePerCm2+=IntegratePowerLawAbove(bin,100.0);
  }
  if(result.denominatorFluencePerCm2<=0.0){
    result.reason="zero fluence above the strict 10 MeV threshold";
    return result;
  }
  result.defined=true;
  result.value=result.numeratorFluencePerCm2/result.denominatorFluencePerCm2;
  result.reason="exact partial-bin integration at strict 10 and 100 MeV thresholds";
  return result;
}

std::string RatioStatusName(RatioStatus status){
  switch(status){
    case RatioStatus::Available: return "available";
    case RatioStatus::ZeroDenominator: return "zero_denominator";
    case RatioStatus::StatisticallyUnresolvedDenominator:
      return "statistically_unresolved_denominator";
    case RatioStatus::InvalidMoments: return "invalid_moments";
  }
  return "invalid_moments";
}

RatioEstimate EstimateRatioOfMeans(const RatioMoments& moments,double z){
  RatioEstimate result;
  if(moments.histories<=0 || !std::isfinite(z) || z<=0.0){
    result.reason="ratio moments require histories and a positive finite quantile";
    return result;
  }
  const double n=static_cast<double>(moments.histories);
  const double values[]={moments.numeratorSum,moments.denominatorSum,
      moments.numeratorSumSquares,moments.denominatorSumSquares,
      moments.numeratorDenominatorSumProducts};
  for(double value:values) if(!std::isfinite(value)){
    result.reason="ratio moments contain a non-finite value";
    return result;
  }
  result.numeratorMean=moments.numeratorSum/n;
  result.denominatorMean=moments.denominatorSum/n;
  if(result.denominatorMean==0.0){
    result.status=RatioStatus::ZeroDenominator;
    result.reason="the estimated denominator is exactly zero";
    return result;
  }
  if(result.denominatorMean<0.0 || moments.numeratorSumSquares<0.0 ||
     moments.denominatorSumSquares<0.0){
    result.reason="ratio moments violate non-negative fluence/LET sums";
    return result;
  }
  if(moments.histories<2){
    result.status=RatioStatus::StatisticallyUnresolvedDenominator;
    result.reason="at least two event histories are required for a finite interval";
    return result;
  }
  const double divisor=n-1.0;
  const double sampleVarN=std::max(0.0,
      (moments.numeratorSumSquares-moments.numeratorSum*moments.numeratorSum/n)/divisor);
  const double sampleVarD=std::max(0.0,
      (moments.denominatorSumSquares-moments.denominatorSum*moments.denominatorSum/n)/divisor);
  const double sampleCov=(moments.numeratorDenominatorSumProducts-
      moments.numeratorSum*moments.denominatorSum/n)/divisor;
  result.numeratorVarianceOfMean=sampleVarN/n;
  result.denominatorVarianceOfMean=sampleVarD/n;
  result.covarianceOfMeans=sampleCov/n;
  result.value=result.numeratorMean/result.denominatorMean;
  const double deltaVariance=(result.numeratorVarianceOfMean+
      result.value*result.value*result.denominatorVarianceOfMean-
      2.0*result.value*result.covarianceOfMeans)/
      (result.denominatorMean*result.denominatorMean);
  result.standardError=std::sqrt(std::max(0.0,deltaVariance));

  // Fieller confidence set for mean(N)/mean(D): A r^2 + B r + C <= 0.
  const double z2=z*z;
  const double A=result.denominatorMean*result.denominatorMean-
      z2*result.denominatorVarianceOfMean;
  const double B=-2.0*(result.numeratorMean*result.denominatorMean-
      z2*result.covarianceOfMeans);
  const double C=result.numeratorMean*result.numeratorMean-
      z2*result.numeratorVarianceOfMean;
  const double discriminant=B*B-4.0*A*C;
  if(A<=0.0 || discriminant<0.0){
    result.status=RatioStatus::StatisticallyUnresolvedDenominator;
    result.reason="the 95% Fieller denominator confidence set includes zero";
    return result;
  }
  const double root=std::sqrt(std::max(0.0,discriminant));
  result.interval95Low=(-B-root)/(2.0*A);
  result.interval95High=(-B+root)/(2.0*A);
  if(result.interval95Low>result.interval95High)
    std::swap(result.interval95Low,result.interval95High);
  result.status=RatioStatus::Available;
  result.reason="finite two-sided 95% Fieller interval from event-level covariance";
  return result;
}

} // namespace shield_sim

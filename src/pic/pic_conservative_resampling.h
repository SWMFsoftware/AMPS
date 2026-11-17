#ifndef AMPS_PIC_CONSERVATIVE_RESAMPLING_H
#define AMPS_PIC_CONSERVATIVE_RESAMPLING_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace PIC { namespace ConservativeResampling {

enum class Status { Success,InvalidInput,NoNullDirection,NumericalRejection };
struct Endpoints {
  Status status=Status::InvalidInput;
  std::string message;
  std::vector<double> positive,negative;
  double positiveProbability=0.0;
  double maximumRelativeResidual=0.0;
  double maximumUnbiasedRelativeError=0.0;
  std::size_t rank=0;
};

// Construct the two conservative, positive-weight endpoints of Gonoskov's
// agnostic downsampling (CPC 271 (2022)108200, arXiv:1607.03755, eqs.5-9).
// Rows are physical per-particle quantities; row zero MUST be all ones. No
// phase-space value is reconstructed. Both endpoints conserve all supplied
// moments, and their probabilities preserve each input weight in expectation.
// These guarantees are checked again after rounding to actual double weights.
inline Endpoints Prepare(const std::vector<double>& weights,
                         const std::vector<std::vector<double>>& rows,
                         double tolerance) {
  Endpoints result;
  const std::size_t n=weights.size(),m=rows.size();
  if(n<2||m==0||!std::isfinite(tolerance)||tolerance<0.0) {
    result.message="resampling requires two weights, constraints and a finite nonnegative tolerance";
    return result;
  }
  for(double weight:weights)if(!std::isfinite(weight)||weight<=0.0) {
    result.message="resampling weights must be finite and positive";return result;
  }
  for(const auto& row:rows) {
    if(row.size()!=n) {result.message="resampling constraint dimension differs from weights";return result;}
    for(double value:row)if(!std::isfinite(value)) {
      result.message="resampling constraints must be finite";return result;
    }
  }
  for(double value:rows.front())if(value!=1.0) {
    result.message="first resampling constraint must preserve represented number";return result;
  }

  // Column-pivoted Householder QR of A^T reveals the constraint row space.
  // Centering any row by a multiple of row0 leaves its nullspace unchanged;
  // scaling each centered row avoids comparisons between SI joules/momenta
  // and dimensionless quantities. Never form A*A^T or invert normal equations.
  std::vector<std::vector<long double>> matrix(n,std::vector<long double>(m));
  for(std::size_t j=0;j<m;++j) {
    const long double reference=j?rows[j][0]:0.0L;
    long double scale=0.0L;
    for(std::size_t i=0;i<n;++i)
      scale=std::max(scale,std::fabs(static_cast<long double>(rows[j][i])-reference));
    for(std::size_t i=0;i<n;++i)matrix[i][j]=scale?
        (static_cast<long double>(rows[j][i])-reference)/scale:0.0L;
  }
  std::vector<std::vector<long double>> reflectors;
  const long double threshold=64.0L*std::numeric_limits<long double>::epsilon()*
      std::max(n,m)*std::sqrt(static_cast<long double>(n));
  for(std::size_t k=0;k<std::min(n,m);++k) {
    std::size_t pivot=k;
    long double maximumNorm2=0.0L;
    for(std::size_t j=k;j<m;++j) {
      long double norm2=0.0L;
      for(std::size_t i=k;i<n;++i)norm2+=matrix[i][j]*matrix[i][j];
      if(norm2>maximumNorm2) {maximumNorm2=norm2;pivot=j;}
    }
    if(std::sqrt(maximumNorm2)<=threshold)break;
    if(pivot!=k)for(std::size_t i=0;i<n;++i)std::swap(matrix[i][pivot],matrix[i][k]);
    std::vector<long double> reflector(n,0.0L);
    for(std::size_t i=k;i<n;++i)reflector[i]=matrix[i][k];
    reflector[k]+=std::copysign(std::sqrt(maximumNorm2),reflector[k]);
    long double norm2=0.0L;
    for(std::size_t i=k;i<n;++i)norm2+=reflector[i]*reflector[i];
    const long double norm=std::sqrt(norm2);
    for(std::size_t i=k;i<n;++i)reflector[i]/=norm;
    for(std::size_t j=k;j<m;++j) {
      long double dot=0.0L;
      for(std::size_t i=k;i<n;++i)dot+=reflector[i]*matrix[i][j];
      for(std::size_t i=k;i<n;++i)matrix[i][j]-=2.0L*reflector[i]*dot;
    }
    reflectors.push_back(reflector);
  }
  result.rank=reflectors.size();
  if(result.rank==n) {
    result.status=Status::NoNullDirection;
    result.message="no reduction: candidate count equals independent constraint rank";
    return result;
  }
  // A column of Q orthogonal to its constraint columns is a null vector of A.
  std::vector<long double> direction(n,0.0L);
  direction[result.rank]=1.0L;
  for(std::size_t k=reflectors.size();k>0;--k) {
    const auto& reflector=reflectors[k-1];
    long double dot=0.0L;
    for(std::size_t i=0;i<n;++i)dot+=reflector[i]*direction[i];
    for(std::size_t i=0;i<n;++i)direction[i]-=2.0L*reflector[i]*dot;
  }
  long double amplitude=0.0L;
  for(long double value:direction)amplitude=std::max(amplitude,std::fabs(value));
  if(!(amplitude>0.0L)) {
    result.status=Status::NumericalRejection;result.message="QR produced an empty null direction";return result;
  }
  for(auto& value:direction)value/=amplitude;
  // Canonical sign prevents implementation/pivot sign choices from swapping
  // outcome labels (and thus a campaign's physical random-bit interpretation).
  for(long double value:direction)if(value!=0.0L) {
    if(value<0.0L)for(auto& entry:direction)entry=-entry;
    break;
  }
  long double plus=std::numeric_limits<long double>::infinity();
  long double minus=plus;
  for(std::size_t i=0;i<n;++i) {
    if(direction[i]<0.0L)plus=std::min(plus,-weights[i]/direction[i]);
    if(direction[i]>0.0L)minus=std::min(minus,weights[i]/direction[i]);
  }
  if(!std::isfinite(plus)||!std::isfinite(minus)||plus<=0.0L||minus<=0.0L) {
    result.status=Status::NumericalRejection;result.message="null direction does not have both weight-conserving signs";return result;
  }
  const long double probability=plus>minus?
      (minus/plus)/(1.0L+minus/plus):1.0L/(1.0L+plus/minus);
  result.positiveProbability=static_cast<double>(probability);
  if(!(result.positiveProbability>0.0&&result.positiveProbability<1.0)) {
    result.status=Status::NumericalRejection;result.message="endpoint probability is outside representable open interval";return result;
  }
  result.positive.resize(n);result.negative.resize(n);
  unsigned positiveZeros=0,negativeZeros=0;
  for(std::size_t i=0;i<n;++i) {
    // Ratio form makes every minimizing component exactly zero by construction
    // and keeps other decreasing weights nonnegative; this is not clipping.
    const long double positive=direction[i]<0.0L?
        weights[i]*(1.0L-plus/(-weights[i]/direction[i])):
        weights[i]+plus*direction[i];
    const long double negative=direction[i]>0.0L?
        weights[i]*(1.0L-minus/(weights[i]/direction[i])):
        weights[i]-minus*direction[i];
    result.positive[i]=static_cast<double>(positive);
    result.negative[i]=static_cast<double>(negative);
    if(!std::isfinite(result.positive[i])||!std::isfinite(result.negative[i])||
       result.positive[i]<0.0||result.negative[i]<0.0) {
      result.status=Status::NumericalRejection;result.message="rounded endpoint has invalid weight";return result;
    }
    positiveZeros+=result.positive[i]==0.0;negativeZeros+=result.negative[i]==0.0;
    const long double expected=result.positiveProbability*static_cast<long double>(result.positive[i])+
        (1.0L-result.positiveProbability)*result.negative[i];
    result.maximumUnbiasedRelativeError=std::max(result.maximumUnbiasedRelativeError,
        static_cast<double>(std::fabs(expected-weights[i])/weights[i]));
  }
  if(!positiveZeros||!negativeZeros) {
    result.status=Status::NumericalRejection;result.message="endpoint did not eliminate a representative";return result;
  }
  for(const auto& row:rows) {
    long double before=0.0L,scale=0.0L,afterPlus=0.0L,afterMinus=0.0L;
    for(std::size_t i=0;i<n;++i) {
      const long double moment=static_cast<long double>(weights[i])*row[i];
      before+=moment;scale+=std::fabs(moment);
      afterPlus+=static_cast<long double>(result.positive[i])*row[i];
      afterMinus+=static_cast<long double>(result.negative[i])*row[i];
    }
    const long double error=std::max(std::fabs(afterPlus-before),std::fabs(afterMinus-before));
    const double relative=scale?static_cast<double>(error/scale):
        (error==0.0L?0.0:std::numeric_limits<double>::infinity());
    result.maximumRelativeResidual=std::max(result.maximumRelativeResidual,relative);
  }
  if(result.maximumRelativeResidual>tolerance||result.maximumUnbiasedRelativeError>tolerance) {
    result.status=Status::NumericalRejection;
    result.message="rounded endpoints failed a declared SI invariant or conditional-unbiasedness gate";
    return result;
  }
  result.status=Status::Success;result.message="positive conservative endpoints validated";
  return result;
}
}} // namespace PIC::ConservativeResampling
#endif

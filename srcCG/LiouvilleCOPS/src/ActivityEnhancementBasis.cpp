#include "ActivityEnhancementBasis.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

namespace ActivityEnhancementBasis {
namespace {

constexpr double Pi=3.141592653589793238462643383279502884;

void ValidateDegrees(int minDegree,int maxDegree) {
  if ((minDegree<MinSupportedDegree)||(minDegree>maxDegree)||(maxDegree>MaxSupportedDegree)) {
    throw std::invalid_argument("activity-enhancement degree range must satisfy 0 <= min <= max <= 8");
  }
}

// Fully normalized real-harmonic factor.  P_l^m below already contains the
// Condon--Shortley (-1)^m phase; sqrt(2) creates separately normalized sine
// and cosine branches for nonzero order.  lgamma evaluates the small factorial
// ratio without depending on an integer-factorial implementation.
double Normalization(int degree,int absOrder) {
  //Only 45 (l,|m|) pairs exist through degree 8.  C++11 guarantees
  //thread-safe initialization of this immutable table, avoiding millions of
  //repeated lgamma calls when a 100,000-facet mesh is initialized.
  using Table=std::array<std::array<double,MaxSupportedDegree+1>,MaxSupportedDegree+1>;
  static const Table values=[]() {
    Table result={};
    for (int l=0;l<=MaxSupportedDegree;l++) for (int m=0;m<=l;m++) {
      const double ratio=std::exp(std::lgamma(l-m+1.0)-std::lgamma(l+m+1.0));
      const double base=std::sqrt((2.0*l+1.0)/(4.0*Pi)*ratio);
      result[l][m]=(m==0) ? base : std::sqrt(2.0)*base;
    }
    return result;
  }();
  return values[degree][absOrder];
}

} // anonymous namespace

std::size_t CoefficientCount(int minDegree,int maxDegree) {
  ValidateDegrees(minDegree,maxDegree);
  return static_cast<std::size_t>((maxDegree+1)*(maxDegree+1)-minDegree*minDegree);
}

std::vector<Mode> Modes(int minDegree,int maxDegree) {
  ValidateDegrees(minDegree,maxDegree);
  std::vector<Mode> result;
  result.reserve(CoefficientCount(minDegree,maxDegree));
  for (int degree=minDegree;degree<=maxDegree;degree++) {
    for (int order=-degree;order<=degree;order++) result.push_back({degree,order});
  }
  return result;
}

std::size_t CoefficientIndex(int degree,int order,int minDegree,int maxDegree) {
  ValidateDegrees(minDegree,maxDegree);
  if ((degree<minDegree)||(degree>maxDegree)||(order<-degree)||(order>degree)) {
    throw std::invalid_argument("activity-enhancement mode is outside the requested degree range");
  }
  return static_cast<std::size_t>(degree*degree+degree+order-minDegree*minDegree);
}

void EvaluateCartesian(double x,double y,double z,int minDegree,int maxDegree,
                       double* output,std::size_t outputSize) {
  ValidateDegrees(minDegree,maxDegree);
  const std::size_t count=CoefficientCount(minDegree,maxDegree);
  if ((output==nullptr)||(outputSize!=count)) {
    throw std::invalid_argument("activity-enhancement output has the wrong size");
  }
  if ((!std::isfinite(x))||(!std::isfinite(y))||(!std::isfinite(z))) {
    throw std::invalid_argument("activity-enhancement Cartesian coordinates must be finite");
  }

  // Nested hypot is robust against avoidable overflow and underflow.  The
  // centroid direction, not its distance from the origin, defines the basis.
  const double rho=std::hypot(x,y),radius=std::hypot(rho,z);
  if ((!std::isfinite(radius))||(radius<=0.0)) {
    throw std::invalid_argument("activity-enhancement Cartesian radius must be finite and positive");
  }
  const double mu=std::max(-1.0,std::min(1.0,z/radius));
  const double phi=(rho==0.0) ? 0.0 : std::atan2(y,x);
  const double root=std::sqrt(std::max(0.0,1.0-mu*mu));

  auto Store=[&](int degree,int absOrder,double legendre,double sine,double cosine) {
    if (degree<minDegree) return;
    const double factor=Normalization(degree,absOrder);
    //The enclosing call validated the entire degree range.  Use the normative
    //closed-form index directly in this hot per-mode path rather than repeating
    //argument validation for every stored value.
    const auto index=[&](int order) {
      return static_cast<std::size_t>(degree*degree+degree+order-minDegree*minDegree);
    };
    if (absOrder==0) {
      output[index(0)]=factor*legendre;
    }
    else {
      output[index(-absOrder)]=factor*legendre*sine;
      output[index(+absOrder)]=factor*legendre*cosine;
    }
  };

  // Explicit Condon--Shortley associated-Legendre recurrence.  pMM advances
  // diagonally between orders; pLm2/pLm1 advance in degree for fixed |m|.
  // Only scalar state is needed for one point.
  double pMM=1.0;
  for (int absOrder=0;absOrder<=maxDegree;absOrder++) {
    if (absOrder>0) pMM=-(2*absOrder-1)*root*pMM;
    const double sine=(absOrder==0) ? 0.0 : std::sin(absOrder*phi);
    const double cosine=(absOrder==0) ? 1.0 : std::cos(absOrder*phi);
    Store(absOrder,absOrder,pMM,sine,cosine);
    if (absOrder==maxDegree) continue;

    double pLm2=pMM;
    double pLm1=(2*absOrder+1)*mu*pMM;
    Store(absOrder+1,absOrder,pLm1,sine,cosine);
    for (int degree=absOrder+2;degree<=maxDegree;degree++) {
      const double pLm=((2*degree-1)*mu*pLm1-(degree+absOrder-1)*pLm2)/(degree-absOrder);
      Store(degree,absOrder,pLm,sine,cosine);
      pLm2=pLm1,pLm1=pLm;
    }
  }
}

std::vector<double> EvaluateCartesian(double x,double y,double z,int minDegree,int maxDegree) {
  std::vector<double> result(CoefficientCount(minDegree,maxDegree));
  EvaluateCartesian(x,y,z,minDegree,maxDegree,result.data(),result.size());
  return result;
}

} // namespace ActivityEnhancementBasis

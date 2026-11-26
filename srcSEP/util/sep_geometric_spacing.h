#ifndef SRCSEP_UTIL_SEP_GEOMETRIC_SPACING_H
#define SRCSEP_UTIL_SEP_GEOMETRIC_SPACING_H
// Dependency-free geometric arc spacing for the reduced material line
// (mesh.field_line.first_segment_length_m).  The reduced line's inner end
// moves with the plasma and is re-extended along B by ONE first-segment arc,
// so that arc bounds how far the line can dip below its start radius.  Fine
// arcs near the Sun keep that excursion small without refining the whole line.
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace SEP { namespace ReducedShock {

// Arcs a q^k (k=0..N-2), q>1, with sum exactly lengthM.  q is the unique root
// of S(q)=a(q^(N-1)-1)/(q-1)=L: S is increasing in q and S(1)=a(N-1)<L.
inline bool GeometricArcs(double length,std::size_t count,double first,
    std::vector<double>* arcs,std::string* error) {
  if(!arcs||!error)return false;
  if(!std::isfinite(length)||length<=0||count<3||!std::isfinite(first)||first<=0||
      !(first<length/static_cast<double>(count-1))) {
    *error="geometric line spacing requires 0<first arc<length/(point_count-1) and at least three points";
    return false;
  }
  const auto total=[&](double q) {
    double sum=0,term=first;
    for(std::size_t k=0;k+1<count;++k){sum+=term;term*=q;}
    return sum;
  };
  // Bisection in log q to roundoff. The upper bracket starts above the
  // continuous estimate and doubles until it encloses the root.
  double low=0,high=std::log(length/first)/static_cast<double>(count-2)+1e-3;
  for(int i=0;i<64&&total(std::exp(high))<length;++i)high*=2;
  if(!(total(std::exp(high))>=length)){*error="geometric line spacing cannot be bracketed";return false;}
  for(int i=0;i<400&&high-low>4*std::numeric_limits<double>::epsilon()*high;++i) {
    const double middle=0.5*(low+high);
    if(total(std::exp(middle))<length)low=middle;else high=middle;
  }
  const double q=std::exp(0.5*(low+high));
  std::vector<double> out;out.reserve(count-1);
  double term=first,sum=0;
  for(std::size_t k=0;k+1<count;++k){out.push_back(term);sum+=term;term*=q;}
  // Remove only the residual bisection roundoff; a larger mismatch is a
  // solver failure rather than a silently rescaled mesh.
  const double scale=length/sum;
  if(!std::isfinite(scale)||std::abs(scale-1)>1e-9){*error="geometric line spacing did not converge";return false;}
  for(double& arc:out)arc*=scale;
  arcs->swap(out);
  return true;
}

} } // namespace SEP::ReducedShock
#endif

#ifndef SRCMOON_SPECIES_TRAJECTORY_INTERVALS_H
#define SRCMOON_SPECIES_TRAJECTORY_INTERVALS_H
// Outward-rounded arithmetic for finite polynomial geometry, not a physical
// uncertainty model. Inputs are the exact represented binary floating values.
// Every elementary operation widens to adjacent representable long doubles;
// no compiler fast-math/reassociation mode is admissible for this service.
#ifdef __FAST_MATH__
#error "Moon trajectory interval geometry requires IEEE finite arithmetic without fast-math"
#endif
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>
namespace Moon { namespace Species { namespace TrajectoryIntervals {
struct Interval {long double lo=0,hi=0;};
using Polynomial=std::vector<Interval>;
inline Interval Point(long double x){return {x,x};}
inline bool Zero(const Interval&a){return a.lo==0&&a.hi==0;}
inline long double Down(long double x){return std::nextafter(x,-std::numeric_limits<long double>::infinity());}
inline long double Up(long double x){return std::nextafter(x,std::numeric_limits<long double>::infinity());}
inline bool Finite(const Interval&a){return std::isfinite(a.lo)&&std::isfinite(a.hi)&&a.lo<=a.hi;}
inline Interval Add(const Interval&a,const Interval&b){if(Zero(a))return b;if(Zero(b))return a;return {Down(a.lo+b.lo),Up(a.hi+b.hi)};}
inline Interval Negate(const Interval&a){return {-a.hi,-a.lo};}
inline Interval Subtract(const Interval&a,const Interval&b){if(a.lo==a.hi&&b.lo==b.hi&&a.lo==b.lo)return Point(0);return Add(a,Negate(b));}
inline Interval Multiply(const Interval&a,const Interval&b){if(Zero(a)||Zero(b))return Point(0);
  const long double p[]={a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi};
  return {Down(*std::min_element(p,p+4)),Up(*std::max_element(p,p+4))};}
inline Interval Divide(const Interval&a,const Interval&b){
  if(b.lo<=0&&b.hi>=0)return {-std::numeric_limits<long double>::infinity(),std::numeric_limits<long double>::infinity()};
  return Multiply(a,{Down(1/b.hi),Up(1/b.lo)});
}
inline Interval Value(const Polynomial&p,const Interval&x){Interval value;for(auto i=p.rbegin();i!=p.rend();++i)value=Add(Multiply(value,x),*i);return value;}
inline Polynomial Product(const Polynomial&a,const Polynomial&b){Polynomial c(a.size()+b.size()-1);for(std::size_t i=0;i<a.size();++i)for(std::size_t j=0;j<b.size();++j)c[i+j]=Add(c[i+j],Multiply(a[i],b[j]));return c;}
inline void Accumulate(Polynomial&a,const Polynomial&b,const Interval&factor=Point(1)){if(a.size()<b.size())a.resize(b.size());for(std::size_t i=0;i<b.size();++i)a[i]=Add(a[i],Multiply(b[i],factor));}

// Cover EVERY possible sign-changing root of the exact represented geometry.
// A complementary interval is excluded only when outward Horner evaluation
// proves a strict sign. Near tangencies/coefficient cancellation therefore
// enlarge real uncertainty brackets instead of disappearing behind a root
// tolerance or an epsilon. Identically zero polynomials cannot change sign
// and need no transition bracket. Work exhaustion leaves output unchanged.
inline bool Isolate(const Polynomial&p,long double timeTolerance,
    std::size_t maximumVisits,std::size_t* visits,std::vector<Interval>*output,
    std::string*error){
  if(!visits||!output||!error||!(timeTolerance>0)||!std::isfinite(timeTolerance)||maximumVisits==0){if(error)*error="invalid explicit interval-geometry work/time controls";return false;}
  bool identicallyZero=true;for(const auto&c:p){if(!Finite(c)){*error="nonfinite interval geometry coefficients";return false;}identicallyZero=identicallyZero&&Zero(c);}
  if(identicallyZero)return true;
  std::vector<Interval> work{{0,1}},uncertain;
  while(!work.empty()){
    if(*visits>=maximumVisits){*error="interval shadow geometry exhausted maximumGeometryIntervals; reduce physical substep, increase the explicit geometry work bound, or supply better-conditioned geometry. No root or exposure was silently dropped.";return false;}
    ++*visits;const auto x=work.back();work.pop_back();const auto range=Value(p,x);
    if(!Finite(range)){*error="interval shadow polynomial overflow; rescale geometry or reduce the physical interval";return false;}
    if(range.lo>0||range.hi<0)continue;
    const long double width=Up(x.hi-x.lo),mid=x.lo+(x.hi-x.lo)/2;
    if(width<=timeTolerance||!(mid>x.lo&&mid<x.hi)){uncertain.push_back(x);continue;}
    work.push_back({mid,x.hi});work.push_back({x.lo,mid});
  }
  output->insert(output->end(),uncertain.begin(),uncertain.end());return true;
}
}}}
#endif

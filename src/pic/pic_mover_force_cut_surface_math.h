// Numerical primitives used by BOTH production force-aware cut-cell tracers.
// No mesh, target force law or observational reference belongs in this file.
#ifndef AMPS_FORCE_CUT_SURFACE_MATH_H
#define AMPS_FORCE_CUT_SURFACE_MATH_H
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <vector>

namespace PIC { namespace Mover { namespace ForceCutSurface { namespace Math {
struct State {double x[3],v[3];};

// A contact callback receives the STORED velocity and normal, not rounded
// D=v*h. Require entry both in wider arithmetic and in the ordinary callback
// dot-product order. No empirical angular tolerance is introduced. A leaving
// or unresolved exactly tangent state must never be reflected into a solid.
inline bool IncomingVelocity(const double* velocity,const double* normal) {
  double ordinary=0;long double wide=0;
  for(int d=0;d<3;++d){ordinary+=velocity[d]*normal[d];wide+=static_cast<long double>(velocity[d])*normal[d];}
  return ordinary<0&&wide<0;
}

// Constant-acceleration position extension of one force quadrature interval.
// Frozen-force velocity is globally first order for general a(x,v,t).
// Evaluating acceleration at the ballistic midpoint gives second order,
// including velocity-dependent Lorentz and rotating-frame acceleration.
template<class Force> State Step(const State& start,double time,double h,
    int order,Force& force,double* startAcceleration=nullptr) {
  double a[3],midAcceleration[3];force(start,time,a);
  State finish;
  for (int d=0;d<3;++d) {
    finish.x[d]=start.x[d]+h*start.v[d]+0.5*h*h*a[d];
    if (startAcceleration) startAcceleration[d]=a[d];
  }
  if (order==2) {
    State middle;
    for (int d=0;d<3;++d) {
      middle.x[d]=start.x[d]+0.5*h*start.v[d]+0.125*h*h*a[d];
      middle.v[d]=start.v[d]+0.5*h*a[d];
    }
    force(middle,time+0.5*h,midAcceleration);
  }
  for (int d=0;d<3;++d)
    finish.v[d]=start.v[d]+h*(order==2?midAcceleration[d]:a[d]);
  return finish;
}

inline long double Value(const std::vector<long double>& c,long double t) {
  long double value=0;for (auto p=c.rbegin();p!=c.rend();++p)value=value*t+*p;
  return value;
}

// All real roots of a polynomial in the CLOSED unit interval. A quadratic
// uses the cancellation-safe q formula. Cubic/quartic sphere intersections
// are bracketed between derivative roots, so a turn-back crossing cannot be
// missed just because both interval endpoints lie on the gas side.
inline std::vector<long double> Roots01(std::vector<long double> c,
    long double tolerance=1.0e-15L) {
  while (c.size()>1 && c.back()==0)c.pop_back();
  std::vector<long double> roots;
  auto add=[&](long double t) {
    if (t>=-tolerance && t<=1+tolerance) {
      t=std::max(0.0L,std::min(1.0L,t));
      if (roots.empty() || std::fabs(t-roots.back())>4*tolerance)roots.push_back(t);
    }
  };
  if (c.size()==1)return roots; // Identically on a plane is not a crossing.
  if (c.size()==2) {add(-c[0]/c[1]);return roots;}
  if (c.size()==3) {
    const long double b=c[1],a=c[2],z=c[0];
    long double discriminant=b*b-4*a*z;
    const long double bound=64*LDBL_EPSILON*(b*b+std::fabs(4*a*z));
    if (discriminant<0 && discriminant>=-bound)discriminant=0;
    if (discriminant<0)return roots;
    const long double q=-0.5L*(b+std::copysign(std::sqrt(discriminant),b));
    std::vector<long double> candidates;
    if (q==0)candidates.push_back(-b/(2*a));
    else {candidates.push_back(q/a);candidates.push_back(z/q);}
    std::sort(candidates.begin(),candidates.end());
    for (auto t:candidates)add(t);
    return roots;
  }
  std::vector<long double> derivative(c.size()-1);
  for (std::size_t i=1;i<c.size();++i)derivative[i-1]=i*c[i];
  std::vector<long double> limits=Roots01(derivative,tolerance);
  limits.insert(limits.begin(),0);limits.push_back(1);
  long double scale=0;for(auto value:c)scale+=std::fabs(value);
  const long double zero=128*LDBL_EPSILON*scale;
  for(std::size_t i=0;i+1<limits.size();++i) {
    long double lo=limits[i],hi=limits[i+1],flo=Value(c,lo),fhi=Value(c,hi);
    if(std::fabs(flo)<=zero)add(lo);
    if((flo<0 && fhi>0)||(flo>0 && fhi<0)) {
      for(int iteration=0;iteration<100 && hi-lo>tolerance;++iteration) {
        const long double mid=(lo+hi)/2, fm=Value(c,mid);
        if(fm==0){lo=hi=mid;break;}
        if((flo<0 && fm>0)||(flo>0 && fm<0)){hi=mid;fhi=fm;}
        else {lo=mid;flo=fm;}
      }
      add((lo+hi)/2);
    }
    if(i+2==limits.size() && std::fabs(fhi)<=zero)add(hi);
  }
  std::sort(roots.begin(),roots.end());
  roots.erase(std::unique(roots.begin(),roots.end(),[&](long double a,long double b){return std::fabs(a-b)<=4*tolerance;}),roots.end());
  return roots;
}

inline void Position(const State& s,const double* a,double t,double* x) {
  for(int d=0;d<3;++d)x[d]=s.x[d]+t*s.v[d]+0.5*t*t*a[d];
}
inline std::vector<long double> PlaneRoots(const State& s,const double* a,
    const double* normal,const double* point,double h,double tolerance) {
  std::vector<long double> c(3,0);
  for(int d=0;d<3;++d) {
    c[0]+=normal[d]*(static_cast<long double>(s.x[d])-point[d]);
    c[1]+=normal[d]*static_cast<long double>(s.v[d])*h;
    c[2]+=0.5L*normal[d]*static_cast<long double>(a[d])*h*h;
  }
  return Roots01(c,std::max(1.0e-18L,static_cast<long double>(tolerance/h)));
}
inline std::vector<long double> SphereRoots(const State& s,const double* a,
    const double* centre,double radius,double h,double tolerance) {
  std::vector<long double> c(5,0);c[0]=-static_cast<long double>(radius)*radius;
  for(int d=0;d<3;++d) {
    const long double x=static_cast<long double>(s.x[d])-centre[d],v=s.v[d]*h,b=0.5L*a[d]*h*h;
    c[0]+=x*x;c[1]+=2*x*v;c[2]+=v*v+2*x*b;c[3]+=2*v*b;c[4]+=b*b;
  }
  return Roots01(c,std::max(1.0e-18L,static_cast<long double>(tolerance/h)));
}
}}}}
#endif

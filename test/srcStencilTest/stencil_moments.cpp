/* Exact polynomial moments and compatible identities for ECSIM stencils. */
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "pic.h"
#include "test_harness.h"
#include "test_register.h"

namespace E=PIC::FieldSolver::Electromagnetic::ECSIM::Stencil;
namespace StencilMoments {
using Scalar=std::function<double(double,double,double)>;
using Key=std::tuple<int,int,int>;
using Coefficients=std::map<Key,double>;

static double apply(cStencil source,const Scalar& value,double shift) {
  cStencil::cStencilData stencil;
  source.ExportStencil(&stencil);
  double result=0.0;
  for (int n=0;n<stencil.Length;++n) {
    const auto& item=stencil.Data[n];
    result+=item.a*value(item.i+shift,item.j+shift,item.k+shift);
  }
  return result;
}

static Coefficients coefficients(cStencil source,bool transpose=false) {
  cStencil::cStencilData stencil;
  source.ExportStencil(&stencil);
  Coefficients result;
  for (int n=0;n<stencil.Length;++n) {
    const auto& item=stencil.Data[n];
    const int sign=transpose ? -1 : 1;
    result[{sign*item.i,sign*item.j,sign*item.k}]+=item.a;
  }
  return result;
}

static Coefficients convolve(const Coefficients& a,const Coefficients& b) {
  Coefficients result;
  for (const auto& left:a) for (const auto& right:b) {
    const auto [li,lj,lk]=left.first;
    const auto [ri,rj,rk]=right.first;
    result[{li+ri,lj+rj,lk+rk}]+=left.second*right.second;
  }
  return result;
}

static void add(Coefficients& destination,const Coefficients& source,double scale=1.0) {
  for (const auto& item:source) destination[item.first]+=scale*item.second;
}
static double max_abs(const Coefficients& source) {
  double result=0.0;
  for (const auto& item:source) result=std::max(result,std::abs(item.second));
  return result;
}

struct Builders {
  const char* name;
  void (*curl)(E::cCurlBStencil*,double,double,double);
  void (*grad_div)(E::cGradDivEStencil*,double,double,double);
  void (*laplacian)(E::cLaplacianStencil*,double,double,double);
  void (*curl_curl)(E::cCurlCurlEStencil*,double,double,double);
};
static const Builders builders[]={
  {"second",E::SecondOrder::InitCurlBStencils_edge_based,
   E::SecondOrder::InitGradDivEBStencils_compact,E::SecondOrder::InitLaplacianStencil,
   E::SecondOrder::InitCurlCurlEStencils_compact},
  {"fourth",E::FourthOrder::InitCurlBStencils,E::FourthOrder::InitGradDivEBStencils,
   E::FourthOrder::InitLaplacianStencil,E::FourthOrder::InitCurlCurlEStencils},
  {"sixth",E::SixthOrder::InitCurlBStencils,E::SixthOrder::InitGradDivEBStencils,
   E::SixthOrder::InitLaplacianStencil,E::SixthOrder::InitCurlCurlEStencils},
  {"eighth",E::EighthOrder::InitCurlBStencils,E::EighthOrder::InitGradDivEBStencils,
   E::EighthOrder::InitLaplacianStencil,E::EighthOrder::InitCurlCurlEStencils}
};

static bool near(double actual,double expected,double tolerance=3.0e-10) {
  if (std::isfinite(actual) && std::abs(actual-expected)<=tolerance) return true;
  std::cerr << "moment mismatch: actual=" << actual << " expected=" << expected
            << " tolerance=" << tolerance << "\n";
  return false;
}

int Run(const std::vector<std::string>&) {
  bool pass=true;
  const Scalar one=[](double,double,double){return 1.0;};
  const Scalar x=[](double a,double,double){return a;};
  const Scalar bx=[](double,double y,double z){return y*y+z;};
  const Scalar by=[](double x0,double,double z){return z*z+x0;};
  const Scalar bz=[](double x0,double y,double){return x0*x0+y;};
  const Scalar ex2=[](double x0,double,double){return x0*x0;};
  const Scalar ey2=[](double,double y,double){return y*y;};
  const Scalar ez2=[](double,double,double z){return z*z;};
  const Scalar r2=[](double x0,double y,double z){return x0*x0+y*y+z*z;};
  const Scalar exy=[](double x0,double y,double){return x0*y;};
  const Scalar eyz=[](double,double y,double z){return y*z;};
  const Scalar ezx=[](double x0,double,double z){return z*x0;};

  for (const Builders& build:builders) {
    E::cCurlBStencil curl[3]; E::cGradDivEStencil grad_div[3];
    E::cLaplacianStencil laplacian; E::cCurlCurlEStencil curl_curl[3];
    build.curl(curl,1.0,1.0,1.0); build.grad_div(grad_div,1.0,1.0,1.0);
    build.laplacian(&laplacian,1.0,1.0,1.0);
    build.curl_curl(curl_curl,1.0,1.0,1.0);

    // Differential rows annihilate constants; second derivatives also
    // annihilate affine data.  Curl samples its input at cell centers, hence
    // the +1/2 coordinate shift used only for that staggered operator.
    for (int row=0;row<3;++row) {
      pass&=near(apply(curl[row].Bx,one,0.5)+apply(curl[row].By,one,0.5)+
                 apply(curl[row].Bz,one,0.5),0.0);
      pass&=near(apply(grad_div[row].Ex,one,0.0)+apply(grad_div[row].Ey,one,0.0)+
                 apply(grad_div[row].Ez,one,0.0),0.0);
      pass&=near(apply(curl_curl[row].Ex,one,0.0)+apply(curl_curl[row].Ey,one,0.0)+
                 apply(curl_curl[row].Ez,one,0.0),0.0);
    }
    pass&=near(apply(laplacian.Ex,one,0.0),0.0);
    pass&=near(apply(laplacian.Ex,x,0.0),0.0);

    // Exact quadratic/mixed-polynomial oracles at target coordinate (0,0,0).
    for (int row=0;row<3;++row) {
      pass&=near(apply(curl[row].Bx,bx,0.5)+apply(curl[row].By,by,0.5)+
                 apply(curl[row].Bz,bz,0.5),1.0);
      pass&=near(apply(grad_div[row].Ex,ex2,0.0)+apply(grad_div[row].Ey,ey2,0.0)+
                 apply(grad_div[row].Ez,ez2,0.0),2.0);
      pass&=near(apply(curl_curl[row].Ex,exy,0.0)+apply(curl_curl[row].Ey,eyz,0.0)+
                 apply(curl_curl[row].Ez,ezx,0.0),1.0);
    }
    pass&=near(apply(laplacian.Ex,r2,0.0),6.0);
    pass&=near(apply(laplacian.Ey,ex2,0.0),2.0);
    pass&=near(apply(laplacian.Ez,ey2,0.0),2.0);

    // The reverse-offset operators are the periodic adjoints mapping corners
    // back to centers.  Coefficient-wise commutators independently enforce
    // both div(curl)=0 and curl(grad)=0, without a field that could mask error.
    const Coefficients dx=coefficients(curl[2].By),dy=coefficients(curl[0].Bz),
                       dz=coefficients(curl[1].Bx);
    const Coefficients dxt=coefficients(curl[2].By,true),
                       dyt=coefficients(curl[0].Bz,true),
                       dzt=coefficients(curl[1].Bx,true);
    Coefficients xy,xz,yz;
    add(xy,convolve(dxt,dy)); add(xy,convolve(dyt,dx),-1.0);
    add(xz,convolve(dxt,dz)); add(xz,convolve(dzt,dx),-1.0);
    add(yz,convolve(dyt,dz)); add(yz,convolve(dzt,dy),-1.0);
    const double identity=std::max({max_abs(xy),max_abs(xz),max_abs(yz)});
    if (!std::isfinite(identity) || identity>2.0e-12) {
      std::cerr << build.name << " compatible identity residual=" << identity << "\n";
      pass=false;
    }
  }
  std::cout << "[ECSIM-STENCIL-MOMENTS-IDENTITIES-001] "
            << (pass ? "PASS" : "FAIL") << "\n";
  return pass ? 0 : 1;
}
void ForceLinkAllTests() {}
}

REGISTER_STENCIL_TEST(StencilMoments,"stencil_moments",
  "Exact constant/affine/polynomial moments and compatible curl identities.");

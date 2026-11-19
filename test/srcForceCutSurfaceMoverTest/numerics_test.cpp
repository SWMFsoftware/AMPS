// Deterministic analytical verification of the EXACT production numerical
// primitives. These are prescribed force/geometry fixtures, not observations.
#include "pic_mover_force_cut_surface_math.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
namespace M=PIC::Mover::ForceCutSurface::Math;
void require(bool condition,const char* name){if(!condition)throw std::runtime_error(name);}
int main() try {
  std::cout.precision(17);
  int checks=0;
  const double up[3]={0,0,1},incoming[3]={0,0,-1},leaving[3]={0,0,1},tangent[3]={1,0,0};
  require(M::IncomingVelocity(incoming,up),"actual incoming velocity accepted");++checks;
  require(!M::IncomingVelocity(leaving,up),"leaving velocity refused");++checks;
  require(!M::IncomingVelocity(tangent,up),"exact tangent is not entry");++checks;
  // Independent Decimal evaluation of these exact stored doubles gives
  // +2.8699027721572765e-18; ordinary rounded arithmetic gives zero. This is
  // the archived real adversarial case, not an invented physical observation.
  const double grazing[3]={.10392890309936924,.011313961546403517,.034216624244487133};
  const double grazingNormal[3]={.3241955823240214,-.15651554791520628,-.93295236087508315};
  require(!M::IncomingVelocity(grazing,grazingNormal),"rounded displacement cannot create an entering actual velocity");++checks;
  auto constant=[](const M::State&,double,double* a){a[0]=.3;a[1]=-.2;a[2]=2;};
  M::State initial{{.1,.2,1},{.4,.5,-3}};
  for(int order:{1,2}){
    auto end=M::Step(initial,0,3,order,constant);
    for(int d=0;d<3;++d){double a[3];constant(initial,0,a);require(std::fabs(end.x[d]-(initial.x[d]+3*initial.v[d]+4.5*a[d]))<2e-15,"constant acceleration position");require(std::fabs(end.v[d]-(initial.v[d]+3*a[d]))<2e-15,"constant acceleration velocity");checks+=2;}
  }
  const double acceleration[3]={0,0,2},normal[3]={0,0,1},plane[3]={0,0,0};
  auto roots=M::PlaneRoots(initial,acceleration,normal,plane,3,1e-13);
  require(roots.size()==2,"turn-back must return both roots");++checks;
  require(std::fabs(static_cast<double>(roots[0])*3-(3-std::sqrt(5.0))/2)<1e-13,"first parabolic plane root");++checks;
  require(std::fabs(static_cast<double>(roots[1])*3-(3+std::sqrt(5.0))/2)<1e-13,"second parabolic plane root");++checks;
  const M::State sphereStart{{-2,0,0},{1,0,0}};const double zero[3]={},centre[3]={};
  auto sphere=M::SphereRoots(sphereStart,zero,centre,1,4,1e-13);
  require(sphere.size()==2&&std::fabs(static_cast<double>(sphere[0])*4-1)<1e-12&&std::fabs(static_cast<double>(sphere[1])*4-3)<1e-12,"quartic/degenerate sphere roots");++checks;
  const M::State acceleratedSphere{{-2,0,0},{4,0,0}};const double deceleration[3]={-2,0,0};
  const auto four=M::SphereRoots(acceleratedSphere,deceleration,centre,1,4,1e-13);
  const double exactFour[4]={2-std::sqrt(3.),1,3,2+std::sqrt(3.)};
  require(four.size()==4,"accelerated sphere must expose all four crossings");++checks;
  for(int i=0;i<4;++i){require(std::fabs(static_cast<double>(four[i])*4-exactFour[i])<1e-12,"genuine quartic crossing oracle");++checks;}
  auto rotation=[](const M::State& s,double,double* a){a[0]=s.v[1];a[1]=-s.v[0];a[2]=0;};
  for(int order:{1,2}){
    double errors[2]={};
    for(int level=0;level<2;++level){const int n=64*(1<<level);const double h=1.0/n;M::State s{{0,0,0},{1,0,0}};for(int i=0;i<n;++i)s=M::Step(s,i*h,h,order,rotation);errors[level]=std::hypot(s.v[0]-std::cos(1.0),s.v[1]+std::sin(1.0));}
    const double ratio=errors[0]/errors[1];require(order==1?(ratio>1.8&&ratio<2.2):(ratio>3.8&&ratio<4.2),"velocity-dependent convergence order");++checks;
    std::cout<<"PASS order="<<order<<" refinement_ratio="<<ratio<<" coarse_error="<<errors[0]<<" fine_error="<<errors[1]<<"\n";
  }
  std::cout<<"PASS production numerical primitives checks="<<checks<<"\n";
}
catch(const std::exception& error) {
  std::cerr<<"FAIL production numerical primitive criterion: "<<error.what()<<"\n";
  return 1;
}

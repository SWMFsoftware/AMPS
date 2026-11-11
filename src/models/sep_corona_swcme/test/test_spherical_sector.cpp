#include "sep_corona_swcme/spherical_sector.h"
#include "sep_corona_swcme/cross_ray_transport.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

namespace XR=SEP::CoronaSwcme::CrossRay;
using SEP::CoronalCME::Norm;
using SEP::CoronalCME::Vec3;
constexpr double kPi=3.141592653589793238462643383279502884;
constexpr double kEps=std::numeric_limits<double>::epsilon();

void Require(bool condition,const std::string& message) {
  if(!condition)throw std::runtime_error(message);
}

void VectorClose(Vec3 got,Vec3 exact,double scale) {
  Require(Norm(got-exact)<=128*kEps*scale,
      "independent spherical-sector vector moment mismatch");
}

void CheckClosure(const XR::SectorMeasure& cell) {
  double area=0;
  for(const auto& face:cell.faces)area+=face.physicalAreaM2;
  Require(Norm(cell.areaClosureResidualM2)<=128*kEps*area,
      "static sector violates integral outward-normal closure");
}

void TestAnalyticGeometryAndUniformFlux() {
  // Octant of the shell 1<r<2 m: V=7pi/6 m^3 and each Cartesian centroid
  // component is 45/56 m. Integrating q over an octant gives (pi/4)(1,1,1).
  // The explicit signed face vectors below follow elementary plane/cone/
  // sphere surface integrals, independent of the production angular primitive.
  const XR::SectorBounds bounds{1,2,0,1,0,kPi/2};
  const auto measured=XR::MeasureSphericalSector(bounds);
  Require(measured.ok(),measured.status.message);
  const auto& cell=measured.value;
  Require(std::abs(cell.volumeM3-7*kPi/6)<=128*kEps*cell.volumeM3,
      "octant shell volume differs from independent sphere fraction");
  VectorClose(cell.volumeCentroidHciM,{45.0/56,45.0/56,45.0/56},2);
  const std::array<Vec3,6> expected={Vec3{-kPi/4,-kPi/4,-kPi/4},
      Vec3{kPi,kPi,kPi},Vec3{0,0,-3*kPi/4},Vec3{},
      Vec3{0,-3*kPi/4,0},Vec3{-3*kPi/4,0,0}};
  for(std::size_t i=0;i<expected.size();++i)
    VectorClose(cell.faces[i].outwardAreaVectorM2,expected[i],4*kPi);
  Require(cell.faces[3].physicalAreaM2==0&&
      Norm(cell.faces[3].outwardAreaVectorM2)==0,
      "polar face was floored to a fictitious finite area");
  Require(cell.faces[1].physicalAreaM2>
      Norm(cell.faces[1].outwardAreaVectorM2),
      "curved physical area was confused with its oriented vector measure");
  CheckClosure(cell);
  const XR::SectorBasis rotated{{0,1,0},{-1,0,0},{0,0,1}};
  const auto other=XR::MeasureSphericalSector(bounds,rotated);
  Require(other.ok(),other.status.message);
  for(std::size_t i=0;i<expected.size();++i)
    VectorClose(other.value.faces[i].outwardAreaVectorM2,
        {-expected[i].y,expected[i].x,expected[i].z},4*kPi);
  CheckClosure(other.value);

  // For one uniform physical state the MHD flux tensor is constant. Its
  // integral through a closed cell must vanish. Applying the actual shared
  // HLL kernel along each integrated face vector tests geometry/flux wiring,
  // including all Cartesian momentum and magnetic induction components.
  // This static identity is not moving-mesh geometric conservation or CT.
  const XR::Thermodynamics gas{5.0/3,1.25663706212e-6};
  const XR::Primitive state{2,{2,3,-1},5,{1e-3,2e-3,3e-3}};
  std::array<double,8> sum{},absolute{};
  for(const auto& face:cell.faces) {
    if(face.physicalAreaM2==0)continue;
    const double projectedArea=Norm(face.outwardAreaVectorM2);
    Require(projectedArea>0,"finite face needs flux quadrature rather than a zero normal");
    const auto f=XR::HllFlux(state,state,
        {face.outwardAreaVectorM2/projectedArea,0},gas);
    Require(f.ok(),f.status.message);
    const auto& x=f.value.value;
    const std::array<double,8> values={x.massKgM2S,x.momentumPa.x,
        x.momentumPa.y,x.momentumPa.z,x.energyWPerM2,
        x.inductionTMPerS.x,x.inductionTMPerS.y,x.inductionTMPerS.z};
    for(std::size_t i=0;i<sum.size();++i) {
      sum[i]+=projectedArea*values[i];
      absolute[i]+=projectedArea*std::abs(values[i]);
    }
  }
  for(std::size_t i=0;i<sum.size();++i)
    Require(std::abs(sum[i])<=256*kEps*absolute[i],
        "uniform MHD state acquires a spurious static geometric source");
}

void TestGlobalMeshMeasures() {
  // Explicit reference tessellation, not inferred legacy ray IDs. Shared
  // bounds are numerically identical; phi is periodic, mu caps are collapsed
  // physical boundaries, and two radial layers cover one complete shell.
  constexpr int nMu=7,nPhi=24;
  const double radii[]={1,2,4};
  std::vector<XR::SectorMeasure> mesh;
  const auto index=[](int radial,int mu,int phi) {
    return (radial*nMu+mu)*nPhi+phi;
  };
  double volume=0,outerArea=0,innerArea=0;
  for(int radial=0;radial<2;++radial)for(int mu=0;mu<nMu;++mu)
    for(int phi=0;phi<nPhi;++phi) {
      const auto cell=XR::MeasureSphericalSector({radii[radial],radii[radial+1],
          -1+2.0*mu/nMu,-1+2.0*(mu+1)/nMu,
          2*kPi*phi/nPhi,2*kPi*(phi+1)/nPhi});
      Require(cell.ok(),cell.status.message);CheckClosure(cell.value);
      mesh.push_back(cell.value);volume+=cell.value.volumeM3;
      if(radial==1)outerArea+=cell.value.faces[1].physicalAreaM2;
      if(radial==0)innerArea+=cell.value.faces[0].physicalAreaM2;
    }
  const double reference=4*kPi*(64-1)/3;
  const double budget=8*(mesh.size()+1)*kEps;
  Require(std::abs(volume-reference)<=budget*reference,
      "complete sector tessellation lost or duplicated physical volume");
  Require(std::abs(outerArea-64*kPi)<=budget*64*kPi&&
      std::abs(innerArea-4*kPi)<=budget*4*kPi,
      "complete shell boundary physical area does not match sphere reference");
  int shared=0;
  const auto opposite=[&](const XR::SectorFaceMeasure& a,
      const XR::SectorFaceMeasure& b) {
    VectorClose(a.outwardAreaVectorM2,-1.0*b.outwardAreaVectorM2,
        a.physicalAreaM2+b.physicalAreaM2);
    Require(std::abs(a.physicalAreaM2-b.physicalAreaM2)<=
        128*kEps*(a.physicalAreaM2+b.physicalAreaM2),
        "two owners disagree on shared physical face area");
    ++shared;
  };
  for(int radial=0;radial<2;++radial)for(int mu=0;mu<nMu;++mu)
    for(int phi=0;phi<nPhi;++phi) {
      const auto& cell=mesh[index(radial,mu,phi)];
      if(radial==0)opposite(cell.faces[1],mesh[index(1,mu,phi)].faces[0]);
      if(mu+1<nMu)opposite(cell.faces[3],mesh[index(radial,mu+1,phi)].faces[2]);
      opposite(cell.faces[5],mesh[index(radial,mu,(phi+1)%nPhi)].faces[4]);
    }
  std::cout<<"[CSWC0628-XRAY-MESH-MEASURES] cells="<<mesh.size()
           <<" shared_faces="<<shared<<" volume_m3="<<volume<<'\n';
}

void TestInvalidAndThinGeometry() {
  const double inner=149597870700.0;
  const double outer=std::nextafter(inner,std::numeric_limits<double>::infinity());
  const auto thin=XR::MeasureSphericalSector({inner,outer,-0.2,0.2,0.1,0.3});
  Require(thin.ok()&&thin.value.volumeM3>0,
      "one-ulp heliocentric shell was replaced with zero volume");
  CheckClosure(thin.value);
  for(const XR::SectorBounds bounds:{XR::SectorBounds{-1,2,0,1,0,1},
      XR::SectorBounds{1,1,0,1,0,1},XR::SectorBounds{1,2,-1.1,1,0,1},
      XR::SectorBounds{1,2,0,0,0,1},XR::SectorBounds{1,2,0,1,0,0},
      XR::SectorBounds{1,2,0,1,0,3*kPi}})
    Require(!XR::MeasureSphericalSector(bounds).ok(),
        "invalid sector was silently clipped or floored");
  XR::SectorBasis basis;basis.x={2,0,0};
  Require(!XR::MeasureSphericalSector({1,2,0,1,0,1},basis).ok(),
      "nonorthonormal frame was silently repaired");
  basis={{-1,0,0},{0,1,0},{0,0,1}};
  Require(!XR::MeasureSphericalSector({1,2,0,1,0,1},basis).ok(),
      "left-handed HCI frame was accepted");
}

} // namespace

int main() {
  const std::pair<const char*,void(*)()> tests[]={
      {"CSWC0628-XRAY-GEOMETRY",TestAnalyticGeometryAndUniformFlux},
      {"CSWC0628-XRAY-MESH-CLOSURE",TestGlobalMeshMeasures},
      {"CSWC0628-XRAY-GEOMETRY-REJECTION",TestInvalidAndThinGeometry}};
  int failures=0;
  for(const auto& test:tests) {
    try {test.second();std::cout<<'['<<test.first<<"] PASS static spherical-sector geometry\n";}
    catch(const std::exception& error) {
      ++failures;std::cout<<'['<<test.first<<"] FAIL "<<error.what()<<'\n';
    }
  }
  return failures?1:0;
}

#include "sep_corona_swcme/cross_ray_mesh.h"
#include "sep_coronal_cme/configuration_parser.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>

#include <type_traits>

static_assert(!std::is_default_constructible<SEP::CoronaSwcme::CrossRay::StaticSphericalMesh>::value,
    "invalid empty authority must not bypass its validated factory");

namespace {

namespace XR=SEP::CoronaSwcme::CrossRay;
using SEP::CoronalCME::Norm;
constexpr double kPi=3.141592653589793238462643383279502884;
constexpr double kEps=std::numeric_limits<double>::epsilon();

void Require(bool condition,const std::string& message) {
  if(!condition)throw std::runtime_error(message);
}

XR::StaticMeshInput Input() {
  XR::StaticMeshInput in;
  in.radiusEdgesM={1,2,4};
  for(int j=0;j<=7;++j)in.muEdges.push_back(-1+2.0*j/7);
  for(int k=0;k<=24;++k)in.phiEdgesRad.push_back(2*kPi*k/24);
  in.periodicPhi=true;in.maximumCells=336;
  return in;
}

std::string Event() {
  return SEP::CoronalCME::ComputeContentChecksum("independent-static-mesh-reference-v1");
}

void TestTopologyAndOwnership() {
  const auto mesh=XR::StaticSphericalMesh::Create(Input(),Event(),12.5,7);
  Require(mesh.ok(),mesh.status.message);const auto& m=*mesh.value;
  // Independent counts for a two-layer shell with 7 equal-mu by 24 periodic
  // longitude cells. Each radial sphere has 6*24 ring vertices plus two
  // unique poles; seams/poles must not duplicate topological nodes.
  Require(m.Cells().size()==336&&m.NodesHciM().size()==438&&
      m.Edges().size()==1228&&m.Faces().size()==1128,
      "spherical shell topology has duplicate/missing nodes, edges, faces or cells");
  const long euler=static_cast<long>(m.NodesHciM().size())-
      static_cast<long>(m.Edges().size())+static_cast<long>(m.Faces().size())-
      static_cast<long>(m.Cells().size());
  Require(euler==2,"hollow shell Euler characteristic is not two");
  std::vector<int> incidence(m.Faces().size(),0);
  int collapsed=0;
  for(const auto& cell:m.Cells()) {
    std::map<std::size_t,int> integerBoundary;
    for(const auto& face:cell.faces) {
      if(face.collapsedGeometricFace) {
        Require(face.index==XR::NoMeshIndex,"collapsed polar face was assigned a finite flux");
        ++collapsed;continue;
      }
      Require(face.index<m.Faces().size()&&std::abs(face.outwardOrientation)==1,
          "physical cell face has no unique orientation/owner");
      ++incidence[face.index];
      for(const auto& edge:m.Faces()[face.index].positiveEdgeLoop)
        integerBoundary[edge.index]+=face.outwardOrientation*edge.orientation;
    }
    for(const auto& term:integerBoundary)
      Require(term.second==0,"cell boundary-of-boundary has a magnetic divergence source");
  }
  int boundary=0,interior=0;
  for(std::size_t i=0;i<m.Faces().size();++i) {
    const auto& face=m.Faces()[i];
    const int owners=(face.leftCell!=XR::NoMeshIndex)+(face.rightCell!=XR::NoMeshIndex);
    Require(incidence[i]==owners&&owners>=1&&owners<=2,
        "shared flux would be missed or applied more than once");
    if(owners==1)++boundary;else ++interior;
  }
  Require(boundary==336&&interior==792&&collapsed==96,
      "shell interface/polar boundary inventory differs from analytic topology");
  std::cout<<"[CSWC0628-XRAY-TOPOLOGY-EVIDENCE] cells="<<m.Cells().size()
           <<" nodes="<<m.NodesHciM().size()<<" edges="<<m.Edges().size()
           <<" faces="<<m.Faces().size()<<" euler="<<euler<<'\n';
}

void TestCurvedEdgeParameterization() {
  const auto mesh=XR::StaticSphericalMesh::Create(Input(),Event(),0,1);
  Require(mesh.ok(),mesh.status.message);const auto& m=*mesh.value;
  int radial=0,meridian=0,latitude=0;
  for(std::size_t i=0;i<m.Edges().size();++i) {
    const auto& edge=m.Edges()[i];
    const auto a=m.SampleEdge(i,0),b=m.SampleEdge(i,1),mid=m.SampleEdge(i,0.5);
    Require(a.ok()&&b.ok()&&mid.ok(),"supported curved edge sample failed");
    const double scale=std::max(edge.start.radiusM,edge.end.radiusM);
    Require(Norm(a.value.positionHciM-m.NodesHciM()[edge.firstNode])<=128*kEps*scale&&
        Norm(b.value.positionHciM-m.NodesHciM()[edge.lastNode])<=128*kEps*scale,
        "curved edge endpoints disagree with canonical vertex identity");
    Require(std::abs(Norm(mid.value.derivativeM)-edge.arcLengthM)<=
        128*kEps*edge.arcLengthM,"edge line element is not its independent analytic arc length");
    if(edge.axis==XR::MeshAxis::Radial)++radial;
    else if(edge.axis==XR::MeshAxis::Mu)++meridian;
    else ++latitude;
  }
  Require(radial==292&&meridian==504&&latitude==432,
      "curved edge families have duplicate/missing shared owners");
  Require(!m.SampleEdge(m.Edges().size(),0).ok()&&!m.SampleEdge(0,-0.1).ok()&&
      !m.SampleEdge(0,1.1).ok(),"edge sample escaped its declared support");
}

void TestIdentityAndTransactionalCreation() {
  const auto first=XR::StaticSphericalMesh::Create(Input(),Event(),1,1);
  const auto second=XR::StaticSphericalMesh::Create(Input(),Event(),2,2);
  Require(first.ok()&&second.ok(),"valid static mesh identity fixture failed");
  Require(first.value->Fingerprint()==second.value->Fingerprint()&&
      second.value->EpochS()==2&&second.value->Generation()==2&&
      first.value->EpochS()==1&&first.value->Generation()==1,
      "geometry identity was conflated with immutable epoch/generation ownership");
  auto changed=Input();changed.radiusEdgesM[1]=2.1;
  const auto different=XR::StaticSphericalMesh::Create(changed,Event(),1,1);
  Require(different.ok()&&different.value->Fingerprint()!=first.value->Fingerprint(),
      "changed radial geometry reused the original mesh fingerprint");
  const auto otherEvent=XR::StaticSphericalMesh::Create(Input(),
      SEP::CoronalCME::ComputeContentChecksum("different-event"),1,1);
  Require(otherEvent.ok()&&otherEvent.value->Fingerprint()!=first.value->Fingerprint(),
      "mesh identity did not bind its owning event");
  changed=Input();changed.maximumCells=335;
  Require(!XR::StaticSphericalMesh::Create(changed,Event(),1,1).ok(),
      "explicit machine/resource cell budget was ignored");
  changed=Input();changed.radiusEdgesM[1]=changed.radiusEdgesM[0];
  Require(!XR::StaticSphericalMesh::Create(changed,Event(),1,1).ok(),
      "unordered radial geometry was silently repaired");
  changed=Input();changed.phiEdgesRad.back()-=0.1;
  Require(!XR::StaticSphericalMesh::Create(changed,Event(),1,1).ok(),
      "incomplete longitude interval was silently made periodic");
  Require(!XR::StaticSphericalMesh::Create(Input(),"missing-event",1,1).ok()&&
      !XR::StaticSphericalMesh::Create(Input(),Event(),1,0).ok(),
      "missing owning identity/generation was accepted");
  Require(first.value->Cells().size()==336&&first.value->EpochS()==1&&
      first.value->Generation()==1,
      "rejected private mesh preparation mutated committed geometry");
}

} // namespace

int main() {
  const std::pair<const char*,void(*)()> tests[]={
      {"CSWC0628-XRAY-TOPOLOGY",TestTopologyAndOwnership},
      {"CSWC0628-XRAY-EDGE-GEOMETRY",TestCurvedEdgeParameterization},
      {"CSWC0628-XRAY-MESH-IDENTITY",TestIdentityAndTransactionalCreation}};
  int failures=0;
  for(const auto& test:tests) {
    try {test.second();std::cout<<'['<<test.first<<"] PASS immutable shared mesh\n";}
    catch(const std::exception& error) {
      ++failures;std::cout<<'['<<test.first<<"] FAIL "<<error.what()<<'\n';
    }
  }
  return failures?1:0;
}

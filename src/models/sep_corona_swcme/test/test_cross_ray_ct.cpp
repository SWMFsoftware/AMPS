#include "sep_corona_swcme/cross_ray_ct.h"
#include "sep_coronal_cme/configuration_parser.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#include <type_traits>

static_assert(!std::is_default_constructible<SEP::CoronaSwcme::CrossRay::StaticMagneticFlux>::value,
    "invalid empty authority must not bypass its validated factory");

namespace {

namespace XR=SEP::CoronaSwcme::CrossRay;
using SEP::CoronalCME::Cross;
using SEP::CoronalCME::Dot;
using SEP::CoronalCME::Norm;
using SEP::CoronalCME::Vec3;
constexpr double kPi=3.141592653589793238462643383279502884;
constexpr double kEps=std::numeric_limits<double>::epsilon();

void Require(bool condition,const std::string& message) {
  if(!condition)throw std::runtime_error(message);
}

std::shared_ptr<const XR::StaticSphericalMesh> Mesh() {
  XR::StaticMeshInput in;in.radiusEdgesM={1,2,4};
  for(int j=0;j<=7;++j)in.muEdges.push_back(-1+2.0*j/7);
  for(int k=0;k<=24;++k)in.phiEdgesRad.push_back(2*kPi*k/24);
  in.periodicPhi=true;in.maximumCells=336;
  const auto mesh=XR::StaticSphericalMesh::Create(in,
      SEP::CoronalCME::ComputeContentChecksum("independent-static-CT-reference-v1"),0,1);
  Require(mesh.ok(),mesh.status.message);return mesh.value;
}

std::vector<double> ConstantFieldFlux(const XR::StaticSphericalMesh& mesh,Vec3 b) {
  std::vector<double> values;
  for(const auto& face:mesh.Faces())values.push_back(Dot(b,face.measure.outwardAreaVectorM2));
  return values;
}

std::vector<double> PotentialLines(const XR::StaticSphericalMesh& mesh,Vec3 b) {
  // Independent Stokes reference A=(B cross X)/2, giving curl A=B. Eight-node
  // Gauss--Legendre integration follows the actual curved edge and integrates
  // both position and line element; a straight chord is not substituted for
  // a latitude/meridian arc. The rule is independently checked on monomials.
  const double x[]={0.18343464249564980494,0.52553240991632898582,
      0.79666647741362673959,0.96028985649753623168};
  const double w[]={0.36268378337836198297,0.31370664587788728734,
      0.22238103445337447054,0.10122853629037625915};
  for(int degree=0;degree<=15;++degree) {
    long double sum=0;
    for(int k=0;k<4;++k)for(double sign:{-1.0,1.0})
      sum+=0.5L*w[k]*std::pow(0.5L*(1+sign*x[k]),degree);
    Require(std::abs(sum-1.0L/(degree+1))<8*kEps,
        "independent Gaussian edge rule lost its polynomial identity");
  }
  std::vector<double> result;
  for(std::size_t edge=0;edge<mesh.Edges().size();++edge) {
    long double sum=0;
    for(int k=0;k<4;++k)for(double sign:{-1.0,1.0}) {
      const auto sample=mesh.SampleEdge(edge,0.5*(1+sign*x[k]));
      Require(sample.ok(),sample.status.message);
      sum+=0.5L*w[k]*Dot(0.5*Cross(b,sample.value.positionHciM),sample.value.derivativeM);
    }
    result.push_back(static_cast<double>(sum));
  }
  return result;
}

XR::EdgeEmfInterval Interval(const XR::StaticMagneticFlux& field,
    std::vector<double> values,double dt) {
  XR::EdgeEmfInterval out;
  out.meshFingerprint=field.Mesh()->Fingerprint();
  out.eventFingerprint=field.Mesh()->EventFingerprint();
  out.fieldGeneration=field.Generation();out.startS=field.EpochS();out.endS=out.startS+dt;
  out.valuesVolt=std::move(values);return out;
}

void TestCurlInitialization() {
  const auto mesh=Mesh();const Vec3 b{2e-5,-3e-5,4e-5};
  const auto actual=XR::StaticMagneticFlux::FromEdgePotential(mesh,PotentialLines(*mesh,b));
  Require(actual.ok(),actual.status.message);
  const auto exact=ConstantFieldFlux(*mesh,b);
  for(std::size_t i=0;i<exact.size();++i) {
    const double tolerance=256*kEps*Norm(b)*mesh->Faces()[i].measure.physicalAreaM2+
        actual.value->NumericalFluxBoundsTm2()[i];
    Require(std::abs(actual.value->FaceFluxTm2()[i]-exact[i])<=tolerance,
        "curved-edge vector-potential curl disagrees with analytic constant B");
  }
  // The initializer must retain raw roundoff, not overwrite diagnostic values
  // with a fabricated exact zero or a cleaned field.
  Require(actual.value->RawDivergenceTPerM().size()==mesh->Cells().size(),
      "CT initialization omitted cell divergence readback");
}

void TestFaradayAndDivergence() {
  const auto mesh=Mesh();const Vec3 initialB{0,0,2e-5},curlE{0,0,1e-6};
  const auto initial=XR::StaticMagneticFlux::FromFaceFlux(
      mesh,ConstantFieldFlux(*mesh,initialB));
  Require(initial.ok(),initial.status.message);
  // E=(curlE cross X)/2 has constant curl E. Faraday's exact solution is
  // B(t)=B0-t*curlE. This manufactured Maxwell reference is not a claim of
  // coupled fluid energy, ideal-MHD edge reconstruction or moving-mesh CT.
  const auto emf=PotentialLines(*mesh,curlE);
  const auto next=initial.value->Advance(Interval(*initial.value,emf,0.25));
  Require(next.ok(),next.status.message);
  const auto exact=ConstantFieldFlux(*mesh,initialB-0.25*curlE);
  for(std::size_t i=0;i<exact.size();++i)
    Require(std::abs(next.value->FaceFluxTm2()[i]-exact[i])<=
        next.value->NumericalFluxBoundsTm2()[i]+256*kEps*
            Norm(initialB)*mesh->Faces()[i].measure.physicalAreaM2,
        "CT Faraday sign/line normalization is inconsistent with analytic solution");
  Require(next.value->EpochS()==0.25&&next.value->Generation()==2&&
      initial.value->EpochS()==0&&initial.value->Generation()==1,
      "CT publication changed committed state or lost epoch/generation lineage");

  // Arbitrary deterministic shared EMFs independently test the topological
  // divergence invariant without presupposing a particular fluid solution.
  auto field=next.value;
  for(int step=0;step<20;++step) {
    std::vector<double> values;
    for(std::size_t edge=0;edge<mesh->Edges().size();++edge)
      values.push_back(1e-5*std::sin(0.37*(edge+1)+0.11*step));
    const auto advanced=field->Advance(Interval(*field,std::move(values),0.001));
    Require(advanced.ok(),advanced.status.message);field=advanced.value;
  }
  double maxDivergence=0;
  for(std::size_t c=0;c<mesh->Cells().size();++c) {
    const auto& cell=mesh->Cells()[c];long double bound=0,closed=0,magnitude=0;
    for(const auto& face:cell.faces)if(!face.collapsedGeometricFace) {
      closed+=face.outwardOrientation*static_cast<long double>(field->FaceFluxTm2()[face.index]);
      bound+=field->NumericalFluxBoundsTm2()[face.index];
      magnitude+=std::abs(field->FaceFluxTm2()[face.index]);
    }
    Require(std::abs(closed)<=bound+64*kEps*magnitude,
        "shared CT circulation created a cell monopole beyond arithmetic uncertainty");
    maxDivergence=std::max(maxDivergence,std::abs(field->RawDivergenceTPerM()[c]));
  }
  std::cout<<"[CSWC0628-XRAY-CT-EVIDENCE] cells="<<mesh->Cells().size()
           <<" faces="<<mesh->Faces().size()<<" edges="<<mesh->Edges().size()
           <<" max_raw_div_b_t_m="<<maxDivergence<<'\n';
}

void TestAuthorityAndRejection() {
  const auto mesh=Mesh();auto flux=ConstantFieldFlux(*mesh,{1e-5,2e-5,3e-5});
  const auto initial=XR::StaticMagneticFlux::FromFaceFlux(mesh,flux);
  Require(initial.ok(),initial.status.message);const auto fingerprint=initial.value->Fingerprint();
  flux[0]+=0.01;
  const auto monopole=XR::StaticMagneticFlux::FromFaceFlux(mesh,flux);
  Require(!monopole.ok()&&monopole.status.code==SEP::Core::StatusCode::DataIntegrityFailure,
      "inconsistent initial magnetic field was silently cleaned");
  std::vector<double> zeros(mesh->Edges().size(),0);
  auto correct=Interval(*initial.value,zeros,0.1);
  auto wrong=correct;wrong.meshFingerprint="stale";
  Require(!initial.value->Advance(wrong).ok(),"foreign CT mesh identity was accepted");
  wrong=correct;wrong.eventFingerprint="foreign";
  Require(!initial.value->Advance(wrong).ok(),"foreign CT event identity was accepted");
  wrong=correct;++wrong.fieldGeneration;
  Require(!initial.value->Advance(wrong).ok(),"foreign CT generation was accepted");
  wrong=correct;wrong.startS=-1;
  Require(!initial.value->Advance(wrong).ok(),"stale CT interval epoch was accepted");
  wrong=correct;wrong.endS=wrong.startS;
  Require(!initial.value->Advance(wrong).ok(),"zero CT interval was accepted");
  wrong=correct;wrong.valuesVolt.pop_back();
  Require(!initial.value->Advance(wrong).ok(),"missing shared edge EMF was silently zero-filled");
  wrong=correct;wrong.valuesVolt[0]=std::numeric_limits<double>::quiet_NaN();
  Require(!initial.value->Advance(wrong).ok(),"nonfinite CT EMF was accepted");
  Require(initial.value->Fingerprint()==fingerprint,
      "failed private CT preparation modified committed magnetic state");
  const auto unchanged=initial.value->Advance(correct);
  Require(unchanged.ok()&&unchanged.value->FaceFluxTm2()==initial.value->FaceFluxTm2()&&
      unchanged.value->Fingerprint()!=fingerprint,
      "zero-curl interval lost exact field data or new epoch identity");
}

} // namespace

int main() {
  const std::pair<const char*,void(*)()> tests[]={
      {"CSWC0628-XRAY-CT-INITIALIZATION",TestCurlInitialization},
      {"CSWC0628-XRAY-CT-FARADAY",TestFaradayAndDivergence},
      {"CSWC0628-XRAY-CT-AUTHORITY",TestAuthorityAndRejection}};
  int failures=0;
  for(const auto& test:tests) {
    try {test.second();std::cout<<'['<<test.first<<"] PASS static shared-edge CT kernel\n";}
    catch(const std::exception& error) {
      ++failures;std::cout<<'['<<test.first<<"] FAIL "<<error.what()<<'\n';
    }
  }
  return failures?1:0;
}

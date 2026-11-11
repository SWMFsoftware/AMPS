#include "sep_corona_swcme/cross_ray_reconstruction.h"
#include "sep_coronal_cme/configuration_parser.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
namespace XR=SEP::CoronaSwcme::CrossRay;
using SEP::CoronalCME::Dot;
using SEP::CoronalCME::Norm;
using SEP::CoronalCME::Vec3;
constexpr double kPi=3.141592653589793238462643383279502884;

void Require(bool ok,const std::string& message) {
  if(!ok)throw std::runtime_error(message);
}

std::shared_ptr<const XR::StaticSphericalMesh> Mesh(int n,bool rotate=false) {
  XR::StaticMeshInput in;in.maximumCells=n*n*n;
  for(int i=0;i<=n;++i) {
    in.radiusEdgesM.push_back(1+1.0*i/n);
    in.muEdges.push_back(-0.6+1.2*i/n);
    in.phiEdgesRad.push_back(0.2+1.0*i/n);
  }
  if(rotate)in.basis={{0,0,1},{1,0,0},{0,1,0}};
  const auto mesh=XR::StaticSphericalMesh::Create(in,
      SEP::CoronalCME::ComputeContentChecksum("independent-CT-field-reconstruction-v1"),5,7);
  Require(mesh.ok(),mesh.status.message);return mesh.value;
}

void Constant() {
  for(bool rotated:{false,true}) {
    const auto mesh=Mesh(4,rotated);
    for(double scale:{1e-12,1.0,1e7}) {
      const Vec3 known=scale*Vec3{0.7,-1.1,0.3};
      std::vector<double> values;
      for(const auto& face:mesh->Faces())
        values.push_back(Dot(known,face.measure.outwardAreaVectorM2));
      const auto field=XR::StaticMagneticFlux::FromFaceFlux(mesh,values);
      Require(field.ok(),field.status.message);
      const auto fingerprint=field.value->Fingerprint();
      const auto cells=XR::ReconstructCellField(*field.value);
      Require(cells.ok(),cells.status.message);
      Require(cells.value.size()==mesh->Cells().size(),"reconstruction omitted cells");
      for(const auto& cell:cells.value) {
        Require(Norm(cell.fieldHciT-known)<3e-14*Norm(known),
            "constant Cartesian field is not recovered under SI scaling/rotation");
        Require(cell.maximumNormalFitResidualT<3e-14*Norm(known)&&
            std::isfinite(cell.pivotRatio)&&cell.pivotRatio>=1,
            "constant field has unexplained face-fit residual/conditioning");
      }
      Require(field.value->Fingerprint()==fingerprint&&field.value->FaceFluxTm2()==values,
          "cell reconstruction modified the authoritative magnetic face flux");
    }
  }
}

double AffineError(int n) {
  const auto mesh=Mesh(n);
  // Independent manufactured field B=(y,x,0) T with coordinates measured in
  // metres; div B=0, A=(0,0,(y^2-x^2)/2) T m. Integrating A along each curved
  // edge and applying Stokes supplies exact-consistent face flux. The known
  // volume average of this affine B is B at the analytic volume centroid.
  // Neither the production reconstruction nor its least-squares equations
  // generate expected field values.
  const double x[]={0.18343464249564980494,0.52553240991632898582,
      0.79666647741362673959,0.96028985649753623168};
  const double w[]={0.36268378337836198297,0.31370664587788728734,
      0.22238103445337447054,0.10122853629037625915};
  std::vector<double> potential;
  for(std::size_t edge=0;edge<mesh->Edges().size();++edge) {
    long double line=0;
    for(int k=0;k<4;++k)for(double sign:{-1.0,1.0}) {
      const auto sample=mesh->SampleEdge(edge,0.5*(1+sign*x[k]));
      Require(sample.ok(),sample.status.message);
      const auto& p=sample.value.positionHciM;
      line+=0.25L*w[k]*(static_cast<long double>(p.y)*p.y-
          static_cast<long double>(p.x)*p.x)*sample.value.derivativeM.z;
    }
    potential.push_back(static_cast<double>(line));
  }
  const auto field=XR::StaticMagneticFlux::FromEdgePotential(mesh,potential);
  Require(field.ok(),field.status.message);
  const auto recovered=XR::ReconstructCellField(*field.value);
  Require(recovered.ok(),recovered.status.message);
  long double numerator=0,denominator=0,fit=0;
  for(std::size_t c=0;c<mesh->Cells().size();++c) {
    const auto& cell=mesh->Cells()[c];const auto& xbar=cell.measure.volumeCentroidHciM;
    const Vec3 exact{xbar.y,xbar.x,0};
    numerator+=cell.measure.volumeM3*Norm(recovered.value[c].fieldHciT-exact);
    denominator+=cell.measure.volumeM3*Norm(exact);
    fit=std::max(fit,static_cast<long double>(recovered.value[c].maximumNormalFitResidualT));
  }
  Require(fit>0,"nonuniform field fit residual was erased");
  return static_cast<double>(numerator/denominator);
}

void Convergence() {
  double previous=0;
  for(int n:{4,8,16}) {
    const double error=AffineError(n);
    std::cout<<"[CSWC0628-XRAY-RECONSTRUCTION-EVIDENCE] n="<<n<<" relative_L1="<<error<<'\n';
    Require(std::isfinite(error)&&error>0,"affine field error is missing/nonfinite");
    if(previous) {
      // A constant local representation is at least first-order consistent;
      // smooth symmetric sectors may exhibit second order. This predeclared
      // lower convergence gate is not a production CME error tolerance.
      const double order=std::log(previous/error)/std::log(2.0);
      Require(order>0.8,"affine magnetic reconstruction is not convergent");
    }
    previous=error;
  }
  Require(previous<0.01,"fine affine magnetic reconstruction exceeds 1% relative L1");
}

void Rejection() {
  XR::StaticMeshInput in;in.radiusEdgesM={1,2};in.muEdges={-0.5,0.5};
  in.phiEdgesRad={0,kPi};in.maximumCells=1;
  const auto mesh=XR::StaticSphericalMesh::Create(in,
      SEP::CoronalCME::ComputeContentChecksum("independent-rank-deficient-sector-v1"),0,1);
  Require(mesh.ok(),mesh.status.message);
  const auto field=XR::StaticMagneticFlux::FromFaceFlux(mesh.value,
      std::vector<double>(mesh.value->Faces().size(),0));
  Require(field.ok(),field.status.message);
  const auto before=field.value->Fingerprint();
  const auto reconstructed=XR::ReconstructCellField(*field.value);
  Require(!reconstructed.ok()&&
      reconstructed.status.code==SEP::Core::StatusCode::NumericalFailure,
      "unresolved geometric field component received an invented value");
  Require(field.value->Fingerprint()==before,"failed reconstruction changed committed flux");
  // Zero physical B on a resolvable cell is valid: no magnetic magnitude floor.
  const auto regular=Mesh(2);
  const auto zero=XR::StaticMagneticFlux::FromFaceFlux(regular,
      std::vector<double>(regular->Faces().size(),0));
  Require(zero.ok(),zero.status.message);
  const auto cells=XR::ReconstructCellField(*zero.value);
  Require(cells.ok(),cells.status.message);
  for(const auto& cell:cells.value)Require(Norm(cell.fieldHciT)==0,
      "zero magnetic field was replaced by a reconstruction floor");
}
} // namespace

int main() {
  try {
    Constant();std::cout<<"[CSWC0628-XRAY-RECONSTRUCTION-CONSTANT] PASS independent static field reconstruction\n";
    Convergence();std::cout<<"[CSWC0628-XRAY-RECONSTRUCTION-CONVERGENCE] PASS independent affine-field refinement\n";
    Rejection();std::cout<<"[CSWC0628-XRAY-RECONSTRUCTION-REJECTION] PASS no regularization or field floor\n";
  }catch(const std::exception& e) {
    std::cerr<<"cross-ray reconstruction FAIL: "<<e.what()<<'\n';return 1;
  }
  return 0;
}

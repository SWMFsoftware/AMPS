#include "sep_corona_swcme/cross_ray_reconstruction.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace SEP { namespace CoronaSwcme { namespace CrossRay {

Core::Result<std::vector<ReconstructedMagneticCell>> ReconstructCellField(
    const StaticMagneticFlux& field) {
  using Return=Core::Result<std::vector<ReconstructedMagneticCell>>;
  const auto& mesh=*field.Mesh();
  std::vector<ReconstructedMagneticCell> result;
  result.reserve(mesh.Cells().size());
  for(const auto& cell:mesh.Cells()) {
    long double matrix[3][3]{},rhs[3]{},factor[3][3]{};
    for(const auto& use:cell.faces)if(!use.collapsedGeometricFace) {
      const auto& face=mesh.Faces()[use.index].measure;
      const auto& s=face.outwardAreaVectorM2;
      const long double a=face.physicalAreaM2;
      const long double n[3]={s.x/a,s.y/a,s.z/a};
      const long double b=field.FaceFluxTm2()[use.index]/a;
      // The same canonical orientation is used for S and Phi. Reversing both
      // for the right cell leaves this least-squares equation unchanged.
      for(int i=0;i<3;++i) {
        rhs[i]+=n[i]*b;
        for(int j=0;j<3;++j)matrix[i][j]+=n[i]*n[j];
      }
    }
    const long double scale=std::max({matrix[0][0],matrix[1][1],matrix[2][2]});
    // Normal equations square the condition number. This numerical boundary
    // rejects geometry with squared conditioning unresolved at double input
    // precision rather than regularizing it into a plausible-looking B.
    const long double minimumPivot=128*std::numeric_limits<double>::epsilon()*scale;
    long double smallest=std::numeric_limits<long double>::infinity(),largest=0;
    for(int i=0;i<3;++i)for(int j=0;j<=i;++j) {
      long double value=matrix[i][j];
      for(int k=0;k<j;++k)value-=factor[i][k]*factor[j][k];
      if(i==j) {
        if(!std::isfinite(value)||!(value>minimumPivot))return Return::Failure(
            Core::StatusCode::NumericalFailure,
            "CT cell-field reconstruction has unresolved geometric rank; no regularization applied");
        factor[i][j]=std::sqrt(value);
        smallest=std::min(smallest,value);largest=std::max(largest,value);
      }else factor[i][j]=value/factor[j][j];
    }
    long double temporary[3]{},answer[3]{};
    for(int i=0;i<3;++i) {
      long double value=rhs[i];
      for(int j=0;j<i;++j)value-=factor[i][j]*temporary[j];
      temporary[i]=value/factor[i][i];
    }
    for(int i=2;i>=0;--i) {
      long double value=temporary[i];
      for(int j=i+1;j<3;++j)value-=factor[j][i]*answer[j];
      answer[i]=value/factor[i][i];
    }
    ReconstructedMagneticCell out;
    out.fieldHciT={static_cast<double>(answer[0]),static_cast<double>(answer[1]),
        static_cast<double>(answer[2])};
    out.pivotRatio=static_cast<double>(largest/smallest);
    if(!std::isfinite(out.fieldHciT.x)||!std::isfinite(out.fieldHciT.y)||
        !std::isfinite(out.fieldHciT.z)||!std::isfinite(out.pivotRatio))return Return::Failure(
            Core::StatusCode::NumericalFailure,"CT cell-field reconstruction overflows SI output");
    for(const auto& use:cell.faces)if(!use.collapsedGeometricFace) {
      const auto& face=mesh.Faces()[use.index].measure;
      const long double prediction=answer[0]*face.outwardAreaVectorM2.x+
          answer[1]*face.outwardAreaVectorM2.y+answer[2]*face.outwardAreaVectorM2.z;
      const double residual=static_cast<double>(std::abs(prediction-
          field.FaceFluxTm2()[use.index])/face.physicalAreaM2);
      if(!std::isfinite(residual))return Return::Failure(Core::StatusCode::NumericalFailure,
          "CT reconstruction fit residual is not representable");
      out.maximumNormalFitResidualT=std::max(out.maximumNormalFitResidualT,residual);
    }
    result.push_back(out);
  }
  return Return::Success(std::move(result));
}

} } } // namespace SEP::CoronaSwcme::CrossRay

#include "perpendicular_transport.h"

#include <algorithm>
#include <cmath>

namespace SEP3D { namespace Transport { namespace {
bool Finite(const Core::Vec3& v) {
  return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}
Core::Status Invalid(const char* text) {
  return Core::Status(Core::StatusCode::InvalidInput, text);
}
}  // namespace

Core::Status BuildPerpendicularBasis(const Core::Vec3& b,
                                     PerpendicularBasis* basis) {
  if (basis == nullptr || !Finite(b) || std::fabs(b.Norm() - 1.0) > 1.0e-12)
    return Invalid("perpendicular basis requires a unit magnetic direction");
  // Pick the Cartesian axis least aligned with b.  This maximizes the cross
  // product norm and avoids the loss of precision caused by an almost
  // parallel seed axis.
  const double ax = std::fabs(b.x), ay = std::fabs(b.y), az = std::fabs(b.z);
  const Core::Vec3 seed = ax <= ay && ax <= az ? Core::Vec3{1,0,0}
                         : ay <= az ? Core::Vec3{0,1,0}
                                    : Core::Vec3{0,0,1};
  basis->first = b.Cross(seed).Normalized();
  basis->second = b.Cross(basis->first).Normalized();
  return basis->first.NormSq() > 0.0 && basis->second.NormSq() > 0.0
      ? Core::Status::OK() : Invalid("perpendicular basis is degenerate");
}

Core::Tensor3 AssembleGyrotropicDiffusionTensor(double kp, double kt,
                                                const Core::Vec3& b) {
  Core::Tensor3 tensor;
  const double v[3] = {b.x,b.y,b.z};
  for (int i=0;i<3;++i) for (int j=0;j<3;++j)
    tensor(i,j) = (i==j ? kt : 0.0) + (kp-kt)*v[i]*v[j];
  return tensor;
}

Core::Vec3 GyrotropicTensorItoDrift(double kp, double dkp,
                                    double kt, double dkt,
                                    const Core::Vec3& b,
                                    const Core::Vec3& curvature,
                                    double divb) {
  // With only field-aligned coefficient gradients available, div(K) is
  // exact for closures whose transverse gradients vanish.  Constant kt has
  // dkt=0; constant-ratio kt=a*kp has dkt=a*dkp.
  return b*dkt + b*(dkp-dkt) +
      (curvature + b*divb)*(kp-kt);
}

Core::Status EvaluateGuidingCenterDrift(const GuidingCenterInput& in,
                                        Core::Vec3* drift) {
  if (drift == nullptr) return Invalid("guiding-centre drift output is null");
  *drift = {};
  if (!in.includeGradientB && !in.includeCurvature) return Core::Status::OK();
  if (!Finite(in.bHat) || !Finite(in.gradAbsBTPerM) ||
      !Finite(in.curvaturePerM) || std::fabs(in.bHat.Norm()-1.0)>1.0e-12 ||
      !std::isfinite(in.absBT) || in.absBT<=0.0 ||
      !std::isfinite(in.momentumKgMPerS) || in.momentumKgMPerS<0.0 ||
      !std::isfinite(in.speedMPerS) || in.speedMPerS<0.0 ||
      !std::isfinite(in.chargeC) || in.chargeC==0.0 ||
      !std::isfinite(in.pitchCosine) || std::fabs(in.pitchCosine)>1.0)
    return Invalid("guiding-centre drift input is outside its physical domain");

  // Relativistic first-order guiding-centre velocities.  p*v replaces m*v^2
  // and therefore remains valid at SEP energies.  For Parker transport the
  // isotropic pitch averages are <mu^2>=1/3 and <(1-mu^2)/2>=1/3.
  const double mu2 = in.pitchAveraged ? 1.0/3.0
                                      : in.pitchCosine*in.pitchCosine;
  const double perpendicularWeight = in.pitchAveraged ? 1.0/3.0
                                      : 0.5*(1.0-mu2);
  if (in.includeGradientB)
    *drift += in.bHat.Cross(in.gradAbsBTPerM) *
        (in.momentumKgMPerS*in.speedMPerS*perpendicularWeight /
         (in.chargeC*in.absBT*in.absBT));
  if (in.includeCurvature)
    *drift += in.bHat.Cross(in.curvaturePerM) *
        (in.momentumKgMPerS*in.speedMPerS*mu2 /
         (in.chargeC*in.absBT));
  return Finite(*drift) ? Core::Status::OK()
                        : Invalid("guiding-centre drift is non-finite");
}

} }  // namespace SEP3D::Transport

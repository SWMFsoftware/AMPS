// ============================================================================
// Controlled cross-field diffusion and guiding-centre drift closures.
//
// This host-neutral module is the single implementation used by both Parker
// and focused transport.  All inputs are SI.  The perpendicular closure is
// scalar in the plane normal to b: K = k_perp I +(k_parallel-k_perp) bb.
// A deterministic basis is used for the two perpendicular Wiener processes;
// isotropy makes its orientation physically irrelevant, while determinism is
// essential for MPI/thread/restart reproducibility.
// ============================================================================
#ifndef SEP3D_TRANSPORT_PERPENDICULAR_TRANSPORT_H
#define SEP3D_TRANSPORT_PERPENDICULAR_TRANSPORT_H

#include "../core/sep3d_types.h"

namespace SEP3D { namespace Transport {

struct PerpendicularBasis { Core::Vec3 first; Core::Vec3 second; };

struct GuidingCenterInput {
  Core::Vec3 bHat;
  Core::Vec3 gradAbsBTPerM;
  Core::Vec3 curvaturePerM;
  double absBT = 0.0;
  double momentumKgMPerS = 0.0;
  double speedMPerS = 0.0;
  double pitchCosine = 0.0;
  double chargeC = 0.0;
  bool pitchAveraged = false;
  bool includeGradientB = false;
  bool includeCurvature = false;
};

Core::Status BuildPerpendicularBasis(const Core::Vec3& bHat,
                                     PerpendicularBasis* basis);
Core::Tensor3 AssembleGyrotropicDiffusionTensor(double kappaParallelM2PerS,
                                                double kappaPerpendicularM2PerS,
                                                const Core::Vec3& bHat);
Core::Vec3 GyrotropicTensorItoDrift(
    double kappaParallelM2PerS, double dKappaParallelDsMPerS,
    double kappaPerpendicularM2PerS, double dKappaPerpendicularDsMPerS,
    const Core::Vec3& bHat, const Core::Vec3& curvaturePerM,
    double divBhatPerM);
Core::Status EvaluateGuidingCenterDrift(const GuidingCenterInput& input,
                                        Core::Vec3* driftMPerS);

} }  // namespace SEP3D::Transport
#endif

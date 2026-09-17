// ============================================================================
// Phase-T adapter from srcSEP3D local state to canonical sep_common kernels.
//
// Include this header only in translation units that evaluate scattering
// coefficients.  Provider-only users should include turbulence_models.h so
// they do not acquire a build-time dependency on the sep_common header path.
// The split is important in production: AMPS compiles main_lib.cpp from the
// copied build/main tree with a historic generic rule that may discard local
// CPPFLAGS/CXXFLAGS, whereas turbulence/*.cpp uses srcSEP3D's explicit rule.
//
// No coefficient formula is implemented here.  CoefficientBridge performs
// typed state conversion and forwards every calculation to sep_common, which
// remains the single owner of D_mumu, lambda_parallel, and kappa_parallel.
// ============================================================================

#ifndef SEP3D_COEFFICIENT_BRIDGE_H
#define SEP3D_COEFFICIENT_BRIDGE_H

#include "turbulence_models.h"

// These headers are intentionally confined to the coefficient boundary.  The
// srcSEP3D makefile supplies SEP_COMMON_DIR when it compiles turbulence/*.cpp
// and the standalone coefficient tests.
#include "sep_coefficient_physics.h"
#include "sep_coefficient_registry.h"

namespace SEP3D {
namespace Turbulence {

class CoefficientBridge final {
 public:
  // Convert the complete local 3-D background/turbulence state into the
  // read-only view accepted by sep_common.  Provenance fields are preserved so
  // a coefficient result can be tied to the exact published wave generation.
  static SEP::Transport::CoefficientPhysics::LocalInputView ToSharedInput(
      const TurbulenceSample& turbulence,
      const Background::BackgroundSample& background,
      double heliocentricRadiusM);

  // Evaluate slab Jokipii D_mumu through the canonical shared implementation.
  // Explicit ballistic state maps to a valid zero-scattering result; missing
  // unapproved turbulence remains an error rather than silently becoming zero.
  static SEP::Transport::CoefficientPhysics::PitchAngleResult JokipiiDmumu(
      const TurbulenceSample& turbulence,
      const Background::BackgroundSample& background,
      double heliocentricRadiusM,
      const SEP::Transport::CoefficientPhysics::SpectrumParameters& spectrum,
      const SEP::Transport::CoefficientPhysics::SpeciesProperties& species,
      double speedMPerS, double mu);

  // Unit-preserving wrappers around the shared lambda/kappa/isotropic-D_mumu
  // conversions.  Keeping these calls thin prevents srcSEP and srcSEP3D from
  // evolving numerically different copies of the same transport relation.
  static SEP::Transport::ScalarResult KappaFromMeanFreePath(
      double lambdaM, double speedMPerS);
  static SEP::Transport::ScalarResult MeanFreePathFromKappa(
      double kappaM2PerS, double speedMPerS);
  static SEP::Transport::ScalarResult IsotropicDmumuFromMeanFreePath(
      double lambdaM, double speedMPerS, double mu);
  static SEP::Transport::ScalarResult MeanFreePathFromIsotropicDmumu(
      double dmumuPerS, double speedMPerS, double mu);
};

}  // namespace Turbulence
}  // namespace SEP3D

#endif  // SEP3D_COEFFICIENT_BRIDGE_H

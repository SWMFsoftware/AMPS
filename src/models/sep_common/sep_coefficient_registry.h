#ifndef SEP_COMMON_SEP_COEFFICIENT_REGISTRY_H
#define SEP_COMMON_SEP_COEFFICIENT_REGISTRY_H

#include "sep_background_snapshot.h"
#include "sep_coefficient_physics.h"
#include "sep_transport_common.h"

#include <string>
#include <vector>

namespace SEP {
namespace Transport {
namespace Coefficient {

enum class SourceMode { Prescribed, SelfConsistent, Swmf };
// ParallelDiffusionLibrary obtains kappa_parallel from the model selected in
// the input file's [parallel_diffusion] section and published by
// SEP::ParallelDiffusion (parallel_diffusion/parallel_diffusion.h) through
// its mover-facing ActiveParallelDiffusion pointer.  It is valid only for the
// Parker mover; the enumerator is appended so existing ordinal values (used
// in ConfigurationFingerprint) are unchanged.
enum class SpatialKind { FromPitchAngle, FromMeanFreePath,
                         ParallelDiffusionLibrary };
// IsotropicFromMeanFreePath evaluates the active mean-free-path provider and
// converts it with the isotropic closure D_mumu = v (1 - mu^2) / (2 lambda),
// dD_mumu/dmu = -v mu / lambda.  With this closure the from-dmumu quadrature
// kappa = (v^2/8) Int (1-mu^2)^2 / D_mumu dmu is exactly v lambda / 3, so the
// parker (from-mfp), fte-mfp and fte-dmumu movers share one lambda_parallel.
// New enumerators are appended so existing ordinals (fingerprints) persist.
enum class PitchAngleKind { Configured, Constant, Jokipii1966, Florinskiy,
                            IsotropicFromMeanFreePath };
enum class MeanFreePathKind {
  Qlt,
  Qlt1,
  Tenishev2005,
  Chen2024,
  FromSpatial,
  // lambda = lambda_ref (r/r_ref)^a (R/R_ref)^b with rigidity R = pc/|q| [V];
  // parameters in Configuration::radialRigidityPowerLaw (no defaults).
  RadialRigidityPowerLaw
};
enum class InvalidPolicy { Fail, Ballistic };
enum class TurbulenceAmplitudePolicy { Reject, LimitToMeanField };

struct Descriptor {
  std::string canonicalName;
  std::string quantity;
  std::string units;
  std::string parameterSchema;
};

struct Configuration {
  // Defaults preserve the pre-Step-10 production choices: Dxx is integrated
  // from configured Dmumu, Dmumu uses the configured callback, and fte-mfp uses
  // the Tenishev et al. analytical mean-free-path parameterization.
  SourceMode source = SourceMode::Prescribed;
  SpatialKind spatial = SpatialKind::FromPitchAngle;
  PitchAngleKind pitchAngle = PitchAngleKind::Configured;
  MeanFreePathKind meanFreePath = MeanFreePathKind::Tenishev2005;
  InvalidPolicy invalidPolicy = InvalidPolicy::Fail;

  // WP17 makes the ninety-degree-resonance decision independent of the MFP
  // invalid-value policy.  The default rejects an unresolved Dmumu integral;
  // Ballistic may be selected only for a mover that directly understands an
  // infinite mean free path.
  CoefficientPhysics::ResonanceGapPolicy resonanceGapPolicy =
      CoefficientPhysics::ResonanceGapPolicy::Reject;
  TurbulenceAmplitudePolicy amplitudePolicy =
      TurbulenceAmplitudePolicy::Reject;

  // Named prescribed-model scales replace QLT/QLT1 literals.  They are part of
  // the configuration fingerprint and are ignored by coupled source adapters
  // except for the spectral shape/bounds explicitly shared with those models.
  double prescribedDeltaBOverB = 0.3;
  double constantDmumuPerS = 0.0;
  double correlationLengthAt1AuM = 0.01 * 1.495978707e11;
  CoefficientPhysics::SpectrumParameters spectrum;
  CoefficientPhysics::FlorinskiyParameters florinskiy;
  CoefficientPhysics::SpatialQuadratureConfiguration spatialQuadrature;
  // Used only when meanFreePath == RadialRigidityPowerLaw.  The default
  // object is invalid by construction (all zero), so selecting that provider
  // without stating lambda_ref, r_ref, R_ref and both exponents is rejected
  // by ValidateConfiguration rather than running with hidden constants.
  CoefficientPhysics::RadialRigidityPowerLawParameters radialRigidityPowerLaw;
};

const std::vector<Descriptor>& SpatialRegistry();
const std::vector<Descriptor>& PitchAngleRegistry();
const std::vector<Descriptor>& MeanFreePathRegistry();

const char* SourceName(SourceMode source);
const char* SpatialName(SpatialKind kind);
const char* PitchAngleName(PitchAngleKind kind);
const char* MeanFreePathName(MeanFreePathKind kind);
const char* InvalidPolicyName(InvalidPolicy policy);
const char* ResonanceGapPolicyName(
    CoefficientPhysics::ResonanceGapPolicy policy);
const char* TurbulenceAmplitudePolicyName(TurbulenceAmplitudePolicy policy);

bool ParseSource(const std::string& text, SourceMode* source);
bool ParseSpatial(const std::string& text, SpatialKind* kind);
bool ParsePitchAngle(const std::string& text, PitchAngleKind* kind);
bool ParseMeanFreePath(const std::string& text, MeanFreePathKind* kind);
bool ParseInvalidPolicy(const std::string& text, InvalidPolicy* policy);
bool ParseResonanceGapPolicy(
    const std::string& text, CoefficientPhysics::ResonanceGapPolicy* policy);
bool ParseTurbulenceAmplitudePolicy(
    const std::string& text, TurbulenceAmplitudePolicy* policy);

Status ValidateConfiguration(const Configuration& configuration);
Status ValidateMoverCompatibility(const Configuration& configuration,
                                  const std::string& moverCanonicalName);
// Stable versioned identity of every coefficient/source scale and policy.
// The fingerprint is published with provider provenance so runs that differ
// only in spectrum or quadrature configuration remain distinguishable.
std::string ConfigurationFingerprint(const Configuration& configuration);
Status ValidateSourceAgainstBackground(SourceMode source,
                                       Background::Provider provider,
                                       Background::Ownership ownership);

// These conversions are centralized because they are valid only under the
// named one-dimensional isotropic-scattering closure. mu is dimensionless;
// speed [m/s], lambda [m], kappa [m^2/s], and Dmumu [s^-1].
ScalarResult KappaFromMeanFreePath(double lambdaM, double speedMPerS);
ScalarResult MeanFreePathFromKappa(double kappaM2PerS, double speedMPerS);
ScalarResult IsotropicDmumuFromMeanFreePath(double lambdaM,
                                            double speedMPerS, double mu);
ScalarResult MeanFreePathFromIsotropicDmumu(double dMuMuPerS,
                                            double speedMPerS, double mu);

Configuration& ActiveConfiguration();
Status SetActiveConfiguration(const Configuration& configuration);

}  // namespace Coefficient
}  // namespace Transport
}  // namespace SEP

#endif  // SEP_COMMON_SEP_COEFFICIENT_REGISTRY_H

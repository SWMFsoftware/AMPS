// ============================================================================
// Phase-T prescribed/AWSoM turbulence models.
//
// This header is deliberately independent of sep_common.  It is included by
// the AMPS-facing main_lib.cpp while that file is compiled from build/main,
// where historic Makefile.conf revisions do not propagate application-local
// include flags to their generic main_lib.cpp rule.  Keeping the provider API
// self-contained also reflects the physics boundary: a turbulence provider
// produces local wave state, while coefficient_bridge.h is the opt-in adapter
// that converts that state to the canonical sep_common scattering API.
// ============================================================================

#ifndef SEP3D_TURBULENCE_MODELS_H
#define SEP3D_TURBULENCE_MODELS_H

#include "turbulence_provider.h"

#include <cstddef>
#include <vector>

namespace SEP3D {
namespace Turbulence {

enum class PrescribedSpectrumModel { PowerLaw, Kolmogorov, Kraichnan };

const char* PrescribedSpectrumModelName(PrescribedSpectrumModel model);

// The historic type name is retained as a source-compatibility surface.  Its
// implementation is now a general normalized prescribed power law whose
// named model is explicit and whose slope is validated against that name.
struct PrescribedKolmogorovConfiguration {
  PrescribedSpectrumModel spectrumModel =
      PrescribedSpectrumModel::Kolmogorov;
  double deltaBOverB = 0.3;
  // sigma_c=(deltaB_+^2-deltaB_-^2)/deltaB^2.  sigma_c=0 is balanced;
  // +1 and -1 are purely one-directional limiting states.
  double normalizedCrossHelicity = 0.0;
  double referenceRadiusM = Core::Const::AU;
  double kMinAtReferencePerM = 1.0e-10;
  double kMaxAtReferencePerM = 1.0e-7;
  double kMinRadialExponent = 2.0;
  double kMaxRadialExponent = 2.0;
  double spectralIndex = 5.0 / 3.0;
  double parallelCorrelationLengthAtReferenceM = 0.03 * Core::Const::AU;
  double correlationLengthRadialExponent = 1.0;
  double validityCadenceS = 60.0;
  std::string coordinateFrame = "HCI-like-inertial";
};

class PrescribedKolmogorovProvider final : public TurbulenceProvider {
 public:
  explicit PrescribedKolmogorovProvider(
      const PrescribedKolmogorovConfiguration& configuration);

  const char* CanonicalName() const override {
    switch (configuration_.spectrumModel) {
      case PrescribedSpectrumModel::PowerLaw:
        return "prescribed-power-law";
      case PrescribedSpectrumModel::Kolmogorov:
        return "prescribed-kolmogorov";
      case PrescribedSpectrumModel::Kraichnan:
        return "prescribed-kraichnan";
    }
    return "prescribed-unknown";
  }
  TurbulenceSource Source() const override {
    return TurbulenceSource::PrescribedPowerLaw;
  }
  Core::Status Validate() const override;
  Core::Status Prepare(double epochS) override;
  // Restart-only provenance restoration.  Physics configuration is already
  // frozen; this restores the published generation without replaying every
  // historical cadence.
  Core::Status PrepareGeneration(double epochS, std::uint64_t generation);
  const TurbulenceMetadata* PreparedMetadata() const override;
  TurbulenceSample Evaluate(
      const Core::Vec3& positionM,
      const Background::BackgroundSample& background) const override;
  std::string ResolvedManifest() const override;

 private:
  PrescribedKolmogorovConfiguration configuration_;
  TurbulenceMetadata metadata_;
  std::uint64_t configurationDigest_ = 0;
  bool prepared_ = false;
};

// AWSoM contract used at the coupling boundary.  wPlus is the total Alfvén
// wave energy density travelling along +B (group velocity U+v_A b-hat), and
// wMinus travels against +B.  Both are J/m^3.  For an equipartitioned Alfvén
// wave, magnetic energy is half the total energy, hence deltaB^2=mu0*w.
struct AwsomWaveRecord {
  Core::Vec3 positionM;
  double epochS = 0.0;
  double wPlusJPerM3 = 0.0;
  double wMinusJPerM3 = 0.0;
  bool complete = false;
};

struct AwsomTurbulenceImport {
  double epochS = 0.0;
  double validUntilS = 0.0;
  std::uint64_t generation = 0;
  std::string coordinateFrame = "HCI-like-inertial";
  std::string configurationFingerprint;
  // AWSoM supplies integrated wave energies.  The scattering closure still
  // needs an explicitly declared spectral band and shape; these values are
  // coupling configuration, not quantities inferred from the two energies.
  double kMinPerM = 1.0e-10;
  double kMaxPerM = 1.0e-7;
  double spectralIndex = 5.0 / 3.0;
  double parallelCorrelationLengthM = 0.03 * Core::Const::AU;
  std::vector<AwsomWaveRecord> records;
};

class AwsomTurbulenceProvider final : public TurbulenceProvider {
 public:
  explicit AwsomTurbulenceProvider(
      MissingTurbulencePolicy policy = MissingTurbulencePolicy::Fail);

  // Load is transactional.  Any malformed record rejects the complete
  // candidate and leaves the last successfully loaded generation available.
  Core::Status Load(const AwsomTurbulenceImport& imported);

  const char* CanonicalName() const override { return "swmf-awsom-waves"; }
  TurbulenceSource Source() const override { return TurbulenceSource::SwmfAwsom; }
  Core::Status Validate() const override;
  Core::Status Prepare(double epochS) override;
  const TurbulenceMetadata* PreparedMetadata() const override;
  TurbulenceSample Evaluate(
      const Core::Vec3& positionM,
      const Background::BackgroundSample& background) const override;
  std::string ResolvedManifest() const override;

 private:
  MissingTurbulencePolicy missingPolicy_;
  AwsomTurbulenceImport imported_;
  TurbulenceMetadata metadata_;
  bool loaded_ = false;
  bool prepared_ = false;
};

// Resonances outside a finite prescribed band require an explicit policy.
// PowerLawExtension evaluates the declared spectrum at the physical resonant
// k; it does not clamp the resonance to an edge and pretend it was in-band.
enum class ResonanceRangePolicy { Reject, PowerLawExtension };

struct SpectrumValue {
  Core::Status status;
  double valueT2M = 0.0;  // one-dimensional P(k), so integral P dk is T^2
  double resolvedWaveNumberPerM = 0.0;
  bool extended = false;
};

class NormalizedPowerLawSpectrum {
 public:
  explicit NormalizedPowerLawSpectrum(const TurbulenceSample& sample);

  Core::Status Validate() const;
  SpectrumValue Evaluate(double waveNumberPerM,
                         ResonanceRangePolicy policy) const;
  double AnalyticBandVarianceT2() const;

 private:
  double varianceT2_ = 0.0;
  double kMinPerM_ = 0.0;
  double kMaxPerM_ = 0.0;
  double index_ = 0.0;
  double normalization_ = 0.0;
};

// sep_common-free result used by the AMPS local-state resolver.  The
// implementation delegates every physical coefficient to CoefficientBridge;
// exposing only primitive SI values here preserves the production header
// boundary that older AMPS Makefile.conf revisions require.
struct LocalScatteringCoefficients {
  Core::Status status;
  double kappaParallelM2PerS = 0.0;
  double dMuMuPerS = 0.0;
  double dDmuMuDmuPerS = 0.0;
  std::uint64_t turbulenceGeneration = 0;
};

// Cell-centred stencil used to recover b-hat dot grad(kappa_parallel) at the
// AMPS application boundary.  The coefficient evaluator remains the canonical
// source of kappa; this helper owns only the finite-difference arithmetic and
// its boundary policy.  At an inner/outer boundary exactly one neighbour may
// be unavailable, in which case a first-order one-sided derivative is explicit
// rather than silently substituting zero.
struct ParallelKappaGradientStencil {
  double centerKappaM2PerS = 0.0;
  double stepM = 0.0;
  bool hasMinus = false;
  double minusKappaM2PerS = 0.0;
  bool hasPlus = false;
  double plusKappaM2PerS = 0.0;
};

Core::Status EvaluateParallelKappaGradient(
    const ParallelKappaGradientStencil& stencil,
    double* dKappaParallelDsMPerS);

LocalScatteringCoefficients EvaluateLocalScattering(
    const TurbulenceSample& turbulence,
    const Background::BackgroundSample& background,
    const Core::Vec3& positionM,
    int modelSpecies,
    double speciesMassKg,
    double signedChargeC,
    double momentumKgMPerS,
    double mu);

}  // namespace Turbulence
}  // namespace SEP3D

#endif  // SEP3D_TURBULENCE_MODELS_H

// ============================================================================
// srcSEP3D/turbulence/turbulence_provider.h
//
// Phase-T turbulence authority contract.
//
// Turbulence is deliberately not a field on BackgroundProvider.  The magnetic
// and plasma background may be analytic while the scattering spectrum is
// imported, or vice versa.  Keeping the authorities separate makes that
// choice visible in configuration, provenance, restart checks, and tests.
//
// LAYER: AMPS-independent model code.  This file may use core/ and
// background/ value types but must never depend on AMPS or MPI declarations.
// ============================================================================

#ifndef SEP3D_TURBULENCE_PROVIDER_H
#define SEP3D_TURBULENCE_PROVIDER_H

#include "../background/bg_provider.h"

#include <cstdint>
#include <string>

namespace SEP3D {
namespace Turbulence {

// Missing coupled wave data are never silently converted to zero scattering.
// Ballistic is an explicit run choice and remains distinguishable from valid
// finite turbulence all the way to the coefficient bridge.
enum class MissingTurbulencePolicy { Fail, Ballistic };

enum class TurbulenceSource {
  PrescribedPowerLaw,
  // Source-compatible alias for callers compiled against the original name.
  // The provider is now generalized to named Kolmogorov/Kraichnan or explicit
  // power-law slopes, all represented by the same normalized finite-band law.
  PrescribedKolmogorov = PrescribedPowerLaw,
  SwmfAwsom
};
enum class TurbulenceOwnership { ModelOwned, ImportedReadOnly };

struct TurbulenceMetadata {
  TurbulenceSource source = TurbulenceSource::PrescribedPowerLaw;
  TurbulenceOwnership ownership = TurbulenceOwnership::ModelOwned;
  double epochS = 0.0;
  double validFromS = 0.0;
  double validUntilS = 0.0;
  std::uint64_t generation = 0;
  std::string coordinateFrame;
  std::string providerIdentity;
  std::string configurationFingerprint;
};

// All values below are SI.  The plus/minus fields refer to propagation along
// and against the local magnetic-field direction, respectively.  The
// outward/inward aliases are resolved only after the local sign of B dot r is
// known; this avoids confusing AWSoM's field-aligned labels with heliocentric
// direction when the magnetic polarity reverses.
struct TurbulenceSample {
  Core::Status status;
  bool valid = false;
  bool ballistic = false;

  double deltaB2T2 = 0.0;
  double deltaBPlus2T2 = 0.0;    // wave propagating along +B
  double deltaBMinus2T2 = 0.0;   // wave propagating against +B
  double deltaBOutward2T2 = 0.0;
  double deltaBInward2T2 = 0.0;

  // Total Alfvén-wave energy density in each field-aligned propagation
  // direction [J/m^3].  Under the equipartition convention shared with AWSoM,
  // w_+/-=deltaB_+/-^2/mu0.  Keeping both representations in the typed sample
  // prevents output code from mistaking magnetic variance [T^2] for energy.
  double waveEnergyPlusJPerM3 = 0.0;
  double waveEnergyMinusJPerM3 = 0.0;

  double kMinPerM = 0.0;
  double kMaxPerM = 0.0;
  double spectralIndex = 0.0;
  double parallelCorrelationLengthM = 0.0;

  std::uint64_t generation = 0;
  std::uint64_t configurationDigest = 0;
  std::string sourceIdentity;
};

class TurbulenceProvider {
 public:
  virtual ~TurbulenceProvider() {}

  virtual const char* CanonicalName() const = 0;
  virtual TurbulenceSource Source() const = 0;
  virtual Core::Status Validate() const = 0;
  virtual Core::Status Prepare(double epochS) = 0;
  virtual const TurbulenceMetadata* PreparedMetadata() const = 0;

  // Background is an input, not an owned authority.  Prescribed amplitudes
  // scale with |B|, and imported wave directions require b-hat and position.
  virtual TurbulenceSample Evaluate(
      const Core::Vec3& positionM,
      const Background::BackgroundSample& background) const = 0;

  virtual std::string ResolvedManifest() const = 0;
};

}  // namespace Turbulence
}  // namespace SEP3D

#endif  // SEP3D_TURBULENCE_PROVIDER_H

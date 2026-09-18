#ifndef SEP_UTIL_SEP_PHYSICS_EXTENSIONS_H
#define SEP_UTIL_SEP_PHYSICS_EXTENSIONS_H

#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)

#include <cstddef>
#include <string>
#include <vector>

namespace SEP {
namespace PhysicsExtensions {

// WP47: the stochastic walker measure is part of the physical equation.  A
// density per arc length q=A*n and a density per volume n require different
// Ito drifts on a nonuniform flux tube; naming the measure prevents an implicit
// extra or missing factor of cross-sectional area in sampling.
using Transport::ParkerMeasure;

struct ParkerGeometryInput {
  ParkerMeasure measure = ParkerMeasure::PerArcLength;
  double plasmaAdvectionMPerS = 0.0;
  double kappaParallelM2PerS = 0.0;
  double dKappaDsMPerS = 0.0;
  double dLnAreaDsPerM = 0.0;
};

Transport::ScalarResult ParkerItoDriftMPerS(
    const ParkerGeometryInput& input);
Transport::ScalarResult WalkerToPhysicalDensity(
    ParkerMeasure measure, double walkerDensity,
    double tubeAreaM2);

// WP48: the signed branch convention is +1 along local B and -1 against B.
// signedWaveNumberPerM retains charge and magnetic polarity.  Finite-bin
// weights linearly partition a point resonance between neighboring log-k bin
// centers, avoiding a discontinuous nearest-bin jump on coarse spectra.
enum class ResonanceStatus {
  Resolved,
  NoRoot,
  SingularComoving,
  OutsideBand,
  Invalid
};

struct DynamicResonanceInput {
  double particleSpeedMPerS = 0.0;
  double mu = 0.0;
  double chargeC = 0.0;
  double massKg = 0.0;
  double signedMagneticFieldT = 0.0;
  double alfvenSpeedMPerS = 0.0;
  double speedOfLightMPerS = 299792458.0;
  int harmonic = 1;
  int branch = 1;
  std::vector<double> binCentersPerM;
};

struct ResonanceMetadata {
  Transport::Status status;
  ResonanceStatus resonanceStatus = ResonanceStatus::Invalid;
  int branch = 0;
  int harmonic = 0;
  int polarization = 0;
  double signedWaveNumberPerM = 0.0;
  std::size_t lowerBin = 0;
  std::size_t upperBin = 0;
  double lowerWeight = 0.0;
  double upperWeight = 0.0;
  std::string reason;
};

ResonanceMetadata SolveDynamicResonance(
    const DynamicResonanceInput& input);

// WP49: a C1 Gaussian resonance-broadening term supplies finite scattering
// through mu=0 while preserving D(+/-1)=0.  amplitudePerS and halfWidthMu are
// physical, configured quantities; amplitude zero reproduces the named pure
// magnetostatic-slab baseline exactly.
struct NinetyDegreeClosure {
  double amplitudePerS = 0.0;
  double halfWidthMu = 0.1;
  std::string model = "gaussian-resonance-broadening";
  std::string provenance;
};

struct PitchAngleCoefficient {
  Transport::Status status;
  double dMuMuPerS = 0.0;
  double derivativePerS = 0.0;
  bool insideClosureRegion = false;
};

PitchAngleCoefficient ApplyNinetyDegreeClosure(
    double mu, double slabDMuMuPerS, double slabDerivativePerS,
    const NinetyDegreeClosure& closure);

// WP50: wave-action/energy conversion is kept in one explicit convention.
// frequencyRadPerS is the positive intrinsic frequency magnitude.  The
// geometric update applies an exact control-volume scaling, while background
// work is kept separate from dissipation in the signed ledger.
enum class WaveInvariant { IntegratedEnergy, WaveAction };

struct WaveActionCell {
  double authoritativeValue = 0.0;
  double volumeM3 = 0.0;
  double intrinsicFrequencyRadPerS = 0.0;
};

struct WaveActionUpdate {
  Transport::Status status;
  WaveActionCell cell;
  double backgroundWorkJ = 0.0;
};

Transport::ScalarResult WaveEnergyJ(const WaveActionCell& cell,
                                    WaveInvariant invariant);
WaveActionUpdate ApplyGeometricConservation(
    const WaveActionCell& oldCell, double newVolumeM3,
    double newIntrinsicFrequencyRadPerS, WaveInvariant invariant);

// WP51: spectral energy is integrated per logarithmic k bin [J].  Interface
// flux arrays have N+1 entries [W]; their telescoping finite-volume difference
// is conservative.  Energy leaving the high-k boundary is deposited into the
// electron and ion heat ledgers rather than disappearing as an unnamed sink.
struct SpectralCascadeResult {
  Transport::Status status;
  std::vector<double> plusEnergyJ;
  std::vector<double> minusEnergyJ;
  double electronHeatJ = 0.0;
  double ionHeatJ = 0.0;
  double closureResidualJ = 0.0;
};

SpectralCascadeResult AdvanceConservativeCascade(
    const std::vector<double>& plusEnergyJ,
    const std::vector<double>& minusEnergyJ,
    const std::vector<double>& plusInterfaceFluxW,
    const std::vector<double>& minusInterfaceFluxW,
    double dtS, double electronHeatingFraction);

// WP52: analytic solar-wind profiles are versioned data, not anonymous helper
// branches.  Cubic Hermite knots provide C1 values and analytic derivatives;
// evaluation outside the calibrated domain is explicitly rejected unless a
// named extrapolation policy is selected.
enum class ProfileExtrapolation { Reject, ConstantEndpoint, PowerLawEndpoint };

struct ProfileKnot {
  double radiusM = 0.0;
  double value = 0.0;
  double derivativePerM = 0.0;
};

struct AnalyticProfile {
  std::string id;
  std::string version;
  std::string units;
  std::string provenance;
  ProfileExtrapolation extrapolation = ProfileExtrapolation::Reject;
  std::vector<ProfileKnot> knots;
};

struct ProfileValue {
  Transport::Status status;
  double value = 0.0;
  double derivativePerM = 0.0;
  bool extrapolated = false;
};

Transport::Status ValidateProfile(const AnalyticProfile& profile);
ProfileValue EvaluateProfile(const AnalyticProfile& profile,
                             double radiusM);
std::string ProfileManifest(const AnalyticProfile& profile);

}  // namespace PhysicsExtensions
}  // namespace SEP

#endif  // SEP_UTIL_SEP_PHYSICS_EXTENSIONS_H

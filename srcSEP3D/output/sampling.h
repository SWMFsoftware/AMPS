// ============================================================================
// Phase-O read-only deterministic sampling.
//
// Particle motion and diagnostics never share writable state. AMPS adapters
// first expose immutable observations, then this module sorts them by stable
// physical identity and reduces cell, virtual-spacecraft, field-line, and
// shock products. Sampling can therefore be repeated without advancing a
// random stream, changing a particle, or depending on traversal order.
// ============================================================================

#ifndef SEP3D_OUTPUT_SAMPLING_H
#define SEP3D_OUTPUT_SAMPLING_H

#include "../adapters/particle_ledger.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace SEP3D {
namespace Output {

struct ParticleObservation {
  std::uint64_t stableId = 0;
  std::uint64_t cellId = 0;
  int species = -1;
  Core::Vec3 positionM;
  double momentumKgMPerS = 0.0;
  double restMassKg = 0.0;
  double mu = 0.0;
  double statisticalWeight = 0.0;
};

struct CellDefinition {
  std::uint64_t cellId = 0;
  Core::Vec3 centerM;
  double volumeM3 = 0.0;
};

struct CellMoment {
  std::uint64_t cellId = 0;
  int species = -1;
  double representedParticles = 0.0;
  double numberDensityM3 = 0.0;
  Core::Vec3 weightedFluxM2PerS;
  double kineticEnergyDensityJPerM3 = 0.0;
  double firstPitchMoment = 0.0;  // <mu>, used for dipole anisotropy 3<mu>
};

struct VirtualSpacecraftDefinition {
  std::string name;
  Core::Vec3 positionM;
  double collectionRadiusM = 0.0;
  std::vector<double> kineticEnergyEdgesJ;
  std::vector<int> acceptedSpecies;
  double minimumMu = -1.0;
  double maximumMu = 1.0;
  std::string observerKind = "fixed-cartesian";
  std::string normalization = "differential-intensity";
};

struct VirtualSpacecraftProduct {
  std::string name;
  int species = -1;
  std::vector<double> kineticEnergyEdgesJ;
  std::vector<double> representedParticlesPerJ;
  // Sum-of-weight-squared propagation.  The square root is the Monte-Carlo
  // standard uncertainty in the same differential units as the spectrum.
  std::vector<double> standardUncertaintyPerJ;
  double dipoleAnisotropy = 0.0;  // 3 * sum(w mu) / sum(w)
  std::string observerKind;
  std::string normalization;
  std::uint64_t acceptedMacroparticles = 0;
};

struct FieldLineProjectionDefinition {
  std::string name;
  Core::Vec3 originM;
  Core::Vec3 direction;
  std::vector<double> distanceEdgesM;
};

struct FieldLineProjection {
  std::string name;
  int species = -1;
  std::vector<double> distanceEdgesM;
  std::vector<double> representedParticlesPerM;
};

struct ShockDiagnostic {
  std::uint64_t step = 0;
  int species = -1;
  std::uint64_t injected = 0;
  std::uint64_t escaped = 0;
  std::uint64_t absorbed = 0;
  std::uint64_t failed = 0;
  std::uint64_t crossings = 0;
};

// These counters are restart-critical because output cadence and cumulative
// observation totals must not jump after a resumed run.
struct SamplingState {
  std::uint64_t completedSamplings = 0;
  std::uint64_t observationsProcessed = 0;
  // Pending-window summaries are checkpointed even before publication.  They
  // let restart validation detect a lost or duplicated observer window.
  std::uint64_t pendingWindows = 0;
  std::uint64_t pendingObservations = 0;
  double pendingRepresentedParticles = 0.0;
};

struct SamplingSnapshot {
  Core::Status status;
  std::vector<CellMoment> cellMoments;
  std::vector<VirtualSpacecraftProduct> spacecraft;
  std::vector<FieldLineProjection> fieldLines;
  std::vector<ShockDiagnostic> shocks;
  SamplingState nextState;
};

struct SamplingRequest {
  std::vector<ParticleObservation> particles;
  std::vector<CellDefinition> cells;
  std::vector<VirtualSpacecraftDefinition> spacecraft;
  std::vector<FieldLineProjectionDefinition> fieldLines;
  std::vector<Adapters::LedgerRow> ledgerRows;
  SamplingState previousState;
};

// Finite representation used by the native AMPS Tecplot callback.
//
// AMPS particle moments and srcSEP3D background quantities have different
// availability rules.  An empty particle sample is a valid physical result,
// while a Cartesian AMR cell outside the configured heliocentric shell has no
// physical background state.  Keeping the three flags separate prevents a
// zero-particle cell from being mistaken for a failed background evaluation.
struct TecplotCellPresentation {
  std::vector<double> backgroundValues;
  double backgroundValid = 0.0;
  double particleSamplingWindowValid = 0.0;
  double particleSamplePresent = 0.0;
};

// Fixed six-column turbulence diagnostic appended to every initialization
// Tecplot record.  Cell storage owns only the two directional magnetic
// variances consumed by scattering; this presentation adds their total and
// converts all three to the AWSoM-compatible total Alfvén-wave energy
// convention w=deltaB^2/mu0.  A named aggregate
// `turbulence_wave_energy_density_J_per_m3` makes the requested pre-existing
// turbulence energy visible without requiring a postprocessor to add w+/w-.
struct TurbulenceTecplotPresentation {
  double deltaB2T2 = 0.0;
  double deltaBPlus2T2 = 0.0;
  double deltaBMinus2T2 = 0.0;
  double waveEnergyJPerM3 = 0.0;
  double waveEnergyPlusJPerM3 = 0.0;
  double waveEnergyMinusJPerM3 = 0.0;
};

// The leading comma is intentional: AMPS has already emitted its native
// VARIABLES entries when the application callback appends this fragment.
const char* TurbulenceTecplotVariableList();

// Validate and derive the immutable six-column presentation.  Negative or
// non-finite directional variance is a storage/initialization error and is
// never serialized as a plausible physical value.
Core::Status PrepareTurbulenceTecplotPresentation(
    double deltaBPlus2T2, double deltaBMinus2T2,
    TurbulenceTecplotPresentation* result);

// Interpolate one complete, cell-centred, static state into the temporary
// centre-node object that AMPS creates while writing a vertex-centred Tecplot
// FEBRICK zone.  AMPS interpolates its built-in sampling and DATAFILE slices,
// but application-requested static bytes are not copied automatically.  The
// production callback therefore supplies pointers to the srcSEP3D slice of
// every node in AMPS' interpolation stencil and uses this AMPS-independent
// routine for the component-wise weighted sum.
//
// All fields in the frozen srcSEP3D static layout are doubles.  Treating the
// slice as one state vector guarantees that background primitives, optional
// gradients, and the two directional turbulence variances use the identical
// AMPS stencil and cannot become spatially misregistered in the output.
Core::Status InterpolateStaticCenterState(
    const double* const* stencilValues, const double* coefficients,
    std::size_t stencilSize, std::size_t valueCount, double* result);

// Convert the internal background/sample state to a finite Tecplot record.
// Undefined background values are represented by zeros and backgroundValid=0
// rather than NaN.  Particle absence is represented independently by
// particleSamplePresent=0; it never invalidates the background record.
TecplotCellPresentation PrepareTecplotCellPresentation(
    const std::vector<double>& storedBackgroundValues,
    bool insidePhysicalShell, long int particleSamplingWindowLength,
    double sampledParticleNumber);

SamplingSnapshot Sample(const SamplingRequest& request);

}  // namespace Output
}  // namespace SEP3D

#endif  // SEP3D_OUTPUT_SAMPLING_H

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
};

struct VirtualSpacecraftProduct {
  std::string name;
  int species = -1;
  std::vector<double> kineticEnergyEdgesJ;
  std::vector<double> representedParticlesPerJ;
  double dipoleAnisotropy = 0.0;  // 3 * sum(w mu) / sum(w)
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

SamplingSnapshot Sample(const SamplingRequest& request);

}  // namespace Output
}  // namespace SEP3D

#endif  // SEP3D_OUTPUT_SAMPLING_H

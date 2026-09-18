#ifndef SEP_UTIL_SEP_RUNTIME_CONTRACTS_H
#define SEP_UTIL_SEP_RUNTIME_CONTRACTS_H

#include "sep_transport_common.h"
#include "sep_turbulence_core.h"

#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace SEP {
namespace RuntimeContracts {

// WP42: the production host is allowed to expose turbulence values through PIC
// data arrays, but those arrays are only synchronized views.  This store is the
// single owner of completed-step counters, handoff state, pending sources,
// restart metadata, and the authoritative E+/E- state for every field line.
// The key is an integer field-line identity because it is stable across AMPS
// pointer relocation and can be serialized without host-library addresses.
struct TurbulenceLineRecord {
  int fieldLineId = -1;
  int ownerRank = -1;
  std::uint64_t geometryGeneration = 0;
  std::string configurationFingerprint;
  Turbulence::State state;
};

enum class PendingRemapPolicy {
  // Rejecting unresolved sources is the safe default: applying a transaction
  // to geometry different from the geometry against which it was generated can
  // deposit physical energy in the wrong cell.
  RejectPending,
  // ApplyBeforeRemap is used only at a driver transaction boundary after all
  // worker/rank records for the old generation have been reduced.
  ApplyBeforeRemap
};

class TurbulenceRuntimeStore {
 public:
  Transport::Status Install(const TurbulenceLineRecord& record);
  bool Contains(int fieldLineId) const;
  const TurbulenceLineRecord* Find(int fieldLineId) const;
  TurbulenceLineRecord* FindMutable(int fieldLineId);

  // RefreshBackground updates geometry coefficients that are sampled from a
  // new immutable background snapshot.  It never imports wave energy from host
  // arrays.  A generation change is handled by the conservative core remap;
  // same-generation changes require an identical cell topology.
  Transport::Status RefreshBackground(
      int fieldLineId, std::uint64_t geometryGeneration,
      const std::string& configurationFingerprint,
      const std::vector<Turbulence::CellState>& geometryAndBackground,
      PendingRemapPolicy pendingPolicy,
      Turbulence::EnergyLedger* remapLedger);

  std::size_t size() const { return lines_.size(); }
  void Clear();

  // The checkpoint contains every per-line turbulence checkpoint as a hex
  // token.  Hex encoding is deliberately boring but robust: embedded newlines
  // from the core schema cannot be mistaken for runtime-store delimiters.
  Transport::Status Serialize(std::string* text) const;
  Transport::Status Deserialize(const std::string& text);

 private:
  std::map<int, TurbulenceLineRecord> lines_;
};

// WP43: a reduced exchange is immutable and is already resolved against the
// frozen field-line generation.  particleEnergyChangeJ is the physical change
// of the represented particle population, including statistical weight.  The
// wave branch receives exactly its negative, so the combined energy ledger
// closes before any limiter or boundary exchange is considered.
struct CouplingTransaction {
  std::uint64_t schema = UINT64_C(0x5352435345504301);
  std::uint64_t transactionId = 0;
  std::uint64_t particleId = 0;
  std::uint64_t generation = 0;
  std::uint64_t event = 0;
  std::uint64_t interval = 0;
  int fieldLineId = -1;
  std::size_t cell = 0;
  int branch = 0;                 // +1 for E+, -1 for E-.
  std::size_t spectralBin = 0;
  double particleEnergyChangeJ = 0.0;
  std::string turbulenceIdentity;
};

enum class InvalidTransactionPolicy { RejectBatch, QuarantineInvalid };

struct CouplingBatchResult {
  Transport::Status status;
  double particleEnergyChangeJ = 0.0;
  double waveEnergyChangeJ = 0.0;
  std::uint64_t accepted = 0;
  std::uint64_t quarantined = 0;
  std::vector<std::uint64_t> quarantinedTransactionIds;
};

// Transactions are sorted and summed by physical identity before a staged
// copy is committed.  Thus an invalid record, stale generation, or insufficient
// branch energy cannot leave a half-applied wave update.
CouplingBatchResult ApplyCouplingTransactions(
    const std::vector<CouplingTransaction>& transactions,
    InvalidTransactionPolicy policy,
    TurbulenceRuntimeStore* store);

// WP44: these derivatives are intentionally separate.  In a curved field,
// d(U.b)/ds = bb:grad(U) + U.kappa, so reusing the former as the strain term is
// physically incorrect even though the two are equal on a straight line.
struct VelocityGradientInput {
  double velocityMPerS[3] = {0.0, 0.0, 0.0};
  double bUnit[3] = {1.0, 0.0, 0.0};
  double curvaturePerM[3] = {0.0, 0.0, 0.0};
  // gradient[i][j] is partial U_i / partial x_j [s^-1].
  double gradientPerS[3][3] = {{0.0, 0.0, 0.0},
                               {0.0, 0.0, 0.0},
                               {0.0, 0.0, 0.0}};
};

struct VelocityDerivatives {
  Transport::Status status;
  double divergencePerS = 0.0;
  double fieldAlignedStrainPerS = 0.0;
  double parallelVelocityGradientPerS = 0.0;
  double curvatureContributionPerS = 0.0;
  std::string method;
};

VelocityDerivatives ComputeVelocityDerivatives(
    const VelocityGradientInput& input);
Transport::ScalarResult ContinuityResidualPerS(
    double materialLogDensityDerivativePerS,
    double divergencePerS);

// WP45: all shock providers populate the same validated SI record.  A
// compression ratio at or below unity is a physical no-shock state, not an
// invitation to substitute a minimum shock speed or clamp a DSA index.
enum class ShockQuality { ValidCompressive, NoCompressiveShock, Invalid };

struct ShockState {
  double epochS = 0.0;
  double normal[3] = {1.0, 0.0, 0.0};
  double shockNormalSpeedMPerS = 0.0;
  double upstreamNormalSpeedMPerS = 0.0;
  double downstreamNormalSpeedMPerS = 0.0;
  double upstreamNumberDensityPerM3 = 0.0;
  double downstreamNumberDensityPerM3 = 0.0;
  double compressionRatio = 0.0;
  double alfvenMach = 0.0;
  double sonicMach = 0.0;
  double obliquityRad = 0.0;
  std::string provider;
  std::string provenance;
};

struct ShockValidation {
  Transport::Status status;
  ShockQuality quality = ShockQuality::Invalid;
  double dsaPhaseSpaceIndex = 0.0;
  std::string reason;
};

ShockValidation ValidateShockState(const ShockState& state);
Transport::ScalarResult ProcessedUpstreamParticleCount(
    const ShockState& state, double shockAreaM2, double dtS,
    double speciesInjectionEfficiency);
Transport::ScalarResult CompressionFromMach(double upstreamMach,
                                            double adiabaticIndex);

// WP46: source identity is based only on persistent integers.  Floating-point
// epochs, MPI rank, OpenMP worker, and particle-buffer addresses are excluded,
// so timestep formatting and decomposition cannot change the event sequence.
enum class SourceRandomPurpose : std::uint64_t {
  EventCount = 1,
  Energy = 2,
  PitchAngle = 3,
  Gyrophase = 4,
  Acceptance = 5
};

struct SourceEventKey {
  std::uint64_t schema = 2;
  std::uint64_t campaign = 0;
  std::uint64_t source = 0;
  std::uint64_t event = 0;
  std::uint64_t fieldLine = 0;
  std::uint64_t species = 0;
  std::uint64_t ordinal = 0;
};

std::uint64_t StableSourceIdentity(const SourceEventKey& key);
Transport::KeyedRandomStream SourceRandomStream(
    const SourceEventKey& key, SourceRandomPurpose purpose);

class SourceEventScheduler {
 public:
  explicit SourceEventScheduler(std::uint64_t campaignSeed = 0)
      : campaignSeed_(campaignSeed) {}
  std::uint64_t AllocateEvent() { return nextEvent_++; }
  std::uint64_t campaignSeed() const { return campaignSeed_; }
  std::uint64_t nextEvent() const { return nextEvent_; }
  Transport::Status Serialize(std::string* text) const;
  Transport::Status Deserialize(const std::string& text);

 private:
  std::uint64_t campaignSeed_ = 0;
  std::uint64_t nextEvent_ = 0;
};

}  // namespace RuntimeContracts
}  // namespace SEP

#endif  // SEP_UTIL_SEP_RUNTIME_CONTRACTS_H

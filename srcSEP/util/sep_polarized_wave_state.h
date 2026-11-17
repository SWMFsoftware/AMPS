#ifndef SRCSEP_UTIL_SEP_POLARIZED_WAVE_STATE_H
#define SRCSEP_UTIL_SEP_POLARIZED_WAVE_STATE_H

#include "sep_resonant_collision.h"
#include <array>
#include <cstdint>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace SEP { namespace Turbulence { namespace Polarized {

const std::uint64_t StateSchema=2;
struct CenteredLogGrid {
  std::vector<double> centersPerM;
  double logWidth=0.0;
};
struct Grid {
  std::vector<double> centersPerM,edgesPerM;
  double logWidth=0.0;
};
struct GridResult {Transport::Status status;Grid grid;};
// User declared the existing k values to be bin centers. With the existing
// uniform log width, e_j=k_0*exp((j-1/2)*Delta ln k), including both outer
// half bins. No endpoint truncation or center-as-boundary reinterpretation.
GridResult ResolveCenteredLogGrid(const CenteredLogGrid& grid);
GridResult ResolveCenteredLogGrid(std::size_t bins,double firstCenterPerM,
                                  double lastCenterPerM);
// Guard the retained two-branch native bridge before reading any source or
// changing its destination. A different live grid is never silently remapped.
Transport::Status CopyLegacyBranchSpectrum(const CenteredLogGrid& requested,
    const CenteredLogGrid& stored,const double* source,std::size_t sourceLength,
    std::vector<double>* destination);

struct BranchPolarization {
  bool known=false;
  double leftFraction=0.0;
  std::string provenance;
};
struct InitialCell {
  std::uint64_t fieldLine=0,segment=0;
  double volumeM3=0.0,meanMagneticFieldT=0.0,massDensityKgM3=0.0;
  // These are the existing InitializeWaveEnergyFromPhysicalParameters
  // outputs in magnetic joules. Counterpropagating balance is preserved;
  // it does not supply a circular-polarization distribution.
  double ePlusJ=0.0,eMinusJ=0.0;
  Transport::ResonantCollision::DomainContract domain;
};
struct Initialization {
  CenteredLogGrid grid;
  double powerLawIndex=5.0/3.0;
  BranchPolarization plus,minus;
  std::vector<InitialCell> cells;
  double vacuumPermeabilityHPerM=0.0;
  double epochS=0.0,validUntilS=0.0;
  std::uint64_t geometryGeneration=0;
  std::string backgroundIdentity,configurationFingerprint,initializerProvenance;
};
struct Cell {
  std::uint64_t fieldLine=0,segment=0;
  double volumeM3=0.0,meanMagneticFieldT=0.0,massDensityKgM3=0.0,alfvenSpeedMPerS=0.0;
  // Sole authoritative storage: [plusLeft,plusRight,minusLeft,minusRight].
  std::array<std::vector<double>,4> magneticEnergyJ;
  // Signed rounding relative to the requested integrated branch initialization.
  double initialRoundoffPlusJ=0.0,initialRoundoffMinusJ=0.0;
  Transport::ResonantCollision::DomainContract domain;
};
struct BranchDiagnostics {Transport::Status status;double ePlusJ=0.0,eMinusJ=0.0;};
BranchDiagnostics DeriveBranches(const Cell& cell);

struct Snapshot {
  std::uint64_t schema=StateSchema,revision=1,geometryGeneration=0;
  double epochS=0.0,validUntilS=0.0,vacuumPermeabilityHPerM=0.0;
  Grid grid;
  std::vector<Cell> cells;
  BranchPolarization initialPlus,initialMinus;
  double initialPowerLawIndex=0.0;
  std::string backgroundIdentity,configurationFingerprint,initializerProvenance;
  std::string identity;
  // Structural validity and explicit initialization do not certify WKB,
  // weak turbulence, dispersion limits, or unfinished PIC coupling adapters.
  bool structurallyValid=false,qualifiedForProduction=false;
};
struct SnapshotResult {
  Transport::Status status;
  std::shared_ptr<const Snapshot> snapshot;
};
SnapshotResult PrepareInitialization(const Initialization& initialization);
// Restart retains the four actual spectra, including independently evolved
// helicities. It never calls the power-law initializer or reconstructs a
// helicity from integrated energy. Validation binds the full value/ledger
// state to its versioned identity; preparation makes an immutable deep copy.
Transport::Status ValidateSnapshot(const Snapshot& snapshot);
SnapshotResult PrepareSnapshotForRestore(const Snapshot& decoded);
// Evolution keeps the exact grid, canonical cell namespace, initializer and
// initial projection ledgers. Only a separate physical operator may supply the
// updated cells and target background realization. This metadata never turns
// structural validity into production qualification.
struct EvolutionMetadata {
  double epochS=0.0,validUntilS=0.0;
  std::uint64_t geometryGeneration=0;
  std::string backgroundIdentity,operationProvenance;
};
struct EvolutionOperationLedger {
  // The operation belongs to its own ledger, not initializerProvenance. The
  // before/after identities bind the staged physical values canonically.
  std::string beforeIdentity,afterIdentity,operationProvenance;
};
struct EvolutionCandidateResult {
  Transport::Status status;
  std::shared_ptr<const Snapshot> snapshot;
  EvolutionOperationLedger operation;
};
EvolutionCandidateResult PrepareEvolutionCandidate(const Snapshot& before,
    const std::vector<Cell>& updatedCells,const EvolutionMetadata& metadata);
struct CollisionViewResult {
  Transport::Status status;
  Transport::ResonantCollision::Spectrum spectrum;
};
CollisionViewResult CollisionView(const Snapshot& snapshot,std::size_t cell);

// Publication is intentionally initialization-only until the common atomic
// particle/wave budget transaction owns updates. Failed candidates, unknown
// polarization and accidental reinitialization never replace a realization.
class StateOwner {
 public:
  Transport::Status PublishInitialization(const Initialization& initialization);
  std::shared_ptr<const Snapshot> Acquire()const;
 private:
  mutable std::mutex mutex_;
  std::shared_ptr<const Snapshot> current_;
};

} } }
#endif

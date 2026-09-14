#ifndef SEP_TURBULENCE_PRODUCTION_ADAPTER_H
#define SEP_TURBULENCE_PRODUCTION_ADAPTER_H

#include "util/sep_transport_common.h"

namespace SEP {
namespace Turbulence {
namespace PICAdapter {

struct ProductionLedger {
  double shockSourceJ = 0.0;
  double particleExchangeJ = 0.0;
  double boundaryExchangeJ = 0.0;
  double physicalDissipationJ = 0.0;
  double limiterCorrectionJ = 0.0;
  double closureResidualJ = 0.0;
  std::uint64_t coreSubcycles = 0;
  std::uint64_t limiterActivations = 0;
  std::uint64_t rejectedUpdates = 0;
  std::uint64_t advectionFaceUpdates = 0;
  std::uint64_t reflectionCellUpdates = 0;
  std::uint64_t cascadeCellUpdates = 0;
  std::uint64_t spectralBinUpdates = 0;
};

struct ShockContribution {
  int fieldLine = -1;
  int segment = -1;
  double plusJ = 0.0;
  double minusJ = 0.0;
  std::string provenance;
};

struct ShockSourceDiagnostics {
  std::uint64_t accepted = 0;
  std::uint64_t noIntersection = 0;
  std::uint64_t invalidGeometry = 0;
  std::uint64_t invalidPhysics = 0;
};

// Shock adapters enqueue typed, immutable source records.  The next Advance
// imports them into CellState::pendingShock*, after which only the common core
// mutates wave energy and closes the signed energy ledger.
Transport::Status QueueShockContribution(const ShockContribution& contribution);
void ClearShockContributions();
const ShockSourceDiagnostics& LastShockSourceDiagnostics();
void RecordShockSourceRejection(bool geometryFailure);

// Advance is the only production mutation entry point for locally evolved
// turbulence.  shockRadiusBeforeM/shockRadiusAfterM are SI metres; NaN disables
// the moving-shock source for library calls that do not own a shock trajectory.
// Particle workers must have completed before this function is entered.
Transport::Status Advance(double dtS, double shockRadiusBeforeM,
                          double shockRadiusAfterM);
const ProductionLedger& LastStepLedger();

}  // namespace PICAdapter
}  // namespace Turbulence
}  // namespace SEP

#endif  // SEP_TURBULENCE_PRODUCTION_ADAPTER_H

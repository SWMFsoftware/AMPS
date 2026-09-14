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
};

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

#ifndef SRCSEP_UTIL_SEP_ANALYTIC_SHOCK_EPOCH_H
#define SRCSEP_UTIL_SEP_ANALYTIC_SHOCK_EPOCH_H

#include "sep_shock_source_core.h"

namespace SEP { namespace Shock {

// One immutable source-departure view of an explicitly supplied trajectory.
// Radius [m], speed [m/s] and epoch [s] describe the same physical instant.
// Configuration is retained by value so callers can validate provenance;
// evaluating a view neither advances an owner cache nor reads the PIC clock.
struct SourceEpochView {
  TrajectoryState trajectory;
  double speedMPerS=0;
};
Transport::Status EvaluateSourceAtEpoch(const TrajectoryConfiguration& configuration,
    double departureEpochS,SourceEpochView* output);
Transport::Status ValidateSourceEpochView(const SourceEpochView& view);
// Bind a source-consumption view to one explicitly prepared owner and declared
// departure [s]. A null owner denotes unprepared authority. Unlike the general
// evaluator above, this rejects other epochs and foreign governing knot tables
// even when their instantaneous radius/speed happen to coincide. It reads only
// the supplied values and never repairs either view or advances owner state.
Transport::Status ValidateSourceDepartureView(const SourceEpochView& view,
    const TrajectoryState* preparedOwner,double departureEpochS);

} namespace ParticleSource { namespace ShockWave { namespace Tenishev2005 {

// Native owner operations: every host rank participates before particle reads
// and after their join. A source query itself is read-only and requires the
// trajectory to have been explicitly prepared; it cannot initialize/advance
// state lazily from a mutable clock during particle injection.
Transport::Status PrepareAnalyticTrajectoryAtEpoch(double departureEpochS);
// Source-specific query: requires an active frozen native driver window and
// its prepared owner at exactly that departure. General arbitrary-epoch math
// remains EvaluateSourceAtEpoch(configuration,epoch,...) above.
Transport::Status QueryAnalyticSourceAtEpoch(double departureEpochS,Shock::SourceEpochView* output);
Transport::Status AdvanceAnalyticTrajectoryToEpoch(double completedEpochS);

// Source consumers use the same validated departure view, never cached end
// geometry. These retain the existing global-zero-centered analytic sphere,
// density and compression laws; no SWCME, thermal, frame or new source closure
// is inferred. Unsupported origin/provider contracts reject explicitly.
Transport::Status InjectionLocationAt(const Shock::SourceEpochView& view,int fieldLine,
    int* segment,double* coordinate,double positionM[3]);
Transport::ScalarResult InjectionRateAt(const Shock::SourceEpochView& view);
Transport::ScalarResult CompressionRatioAt(const Shock::SourceEpochView& view);

} } } }
#endif

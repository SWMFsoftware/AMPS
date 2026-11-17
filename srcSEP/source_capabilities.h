#ifndef SRCSEP_SOURCE_CAPABILITIES_H
#define SRCSEP_SOURCE_CAPABILITIES_H
#include "util/sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)

namespace SEP { namespace Sources {
enum class CartesianSource { GalacticCosmicRayBox, ThermalOuterBox,
                              SphericalInnerBoundary, SphericalShock };

// These legacy samplers lack a physical field-line connection and persisted
// collective source-event scheduler. The typed preflight is observational: no
// allocation, random draw, or provider write is permitted. Their dimensional
// destination is srcSEP3D; this declaration does not assert a completed port of
// a particular legacy sampler into that independent application.
Transport::Status CartesianBirthPreflight(CartesianSource source);

// The retained generic helper stores lab-frame velocities and pre-drifts
// newborns before the full PIC mover pass. Available SEP birth bookkeeping
// does not implement the selected relativistic Maxwell-Juttner flux/bridge, so reject this selected
// physical source before even reserving its persistent source-event counter.
Transport::Status MaxwellianFieldLinePreflight();
// Optional solar-wind helpers retain distinct unqualified thermal variance,
// inflow weighting and global-RNG laws. Neither an initial-population label nor
// new birth metadata qualifies those physics; reject before source mutation.
Transport::Status ThermalPopulationPreflight();
Transport::Status ThermalBeginningPreflight();
// A moving-shock source view must describe the exact frozen birth departure
// [s]. The host prepares the analytic owner outside particle reads, injection
// evaluates an immutable explicit-epoch view, and the host later commits end
// geometry for wave work. This guard compares those declared epochs only; it
// never substitutes the mutable PIC counter, advances state, or reserves RNG.
Transport::Status MovingShockDeparturePreflight(double departureEpochS,
                                               double trajectoryQueryEpochS);

// Call before installing/invoking physical source callbacks. Controlled tests
// must explicitly leave unsupported callbacks uninstalled, not bypass this
// configuration check while retaining a live Cartesian producer.
Transport::Status ValidateConfiguredSources();
} }
#endif

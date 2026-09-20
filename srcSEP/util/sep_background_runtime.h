#ifndef SEP_UTIL_BACKGROUND_RUNTIME_H
#define SEP_UTIL_BACKGROUND_RUNTIME_H

#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_background_snapshot.h)

#include <string>

namespace SEP {
namespace Background {

// PIC owns the application clock.  All SEP physics, output, and snapshot code
// reaches that clock through this function; srcSEP does not maintain an elapsed
// time counter of its own.  Returned units are seconds from PIC's configured
// simulation-time origin.
double SimulationTimeSeconds();

// Resolve the background authority selected by the compiled coupling mode and
// frozen run configuration.  Startup code uses this query before any local
// initializer writes field-line plasma or turbulence storage; it must never
// infer authority from whether a legacy array happens to contain zeroes.
Provider ConfiguredProvider();

// Build the fingerprint of background-affecting configuration.  Mover choice
// is intentionally excluded: selecting Parker versus focused transport changes
// the particle operator, never the solar-wind/IMF/shock/turbulence authority.
std::string CurrentConfigurationFingerprint();

// Publish a model-owned analytic or SWCME state.  Repeated calls may advance
// the state epoch but retain the same field-line generation when geometry has
// not been regenerated.  valid_until_seconds is the last clock value for which
// the state may be consumed.
void PublishModelOwnedSnapshot(Provider provider,
                               double epoch_seconds,
                               double valid_until_seconds,
                               const std::string& provenance);

// Inspect the configured provider immediately before PIC::TimeStep().  For an
// SWMF run this publishes a new read-only generation whenever a new coupling
// epoch is observed.  For standalone runs it verifies that the driver-published
// analytic/SWCME state covers the authoritative simulation time.
void PrepareSnapshotForParticleStep();

// Explicitly transfer an imported SWMF generation to a private local-evolution
// copy.  The copy itself is made by the caller; this function records the new
// owner, epoch, provenance, and generation and rejects every non-SWMF source.
void PublishLocalEvolutionHandoff(double epoch_seconds,
                                  double valid_until_seconds,
                                  const std::string& provenance);

}  // namespace Background
}  // namespace SEP

#endif  // SEP_UTIL_BACKGROUND_RUNTIME_H

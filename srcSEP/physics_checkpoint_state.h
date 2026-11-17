#ifndef SRCSEP_PHYSICS_CHECKPOINT_STATE_H
#define SRCSEP_PHYSICS_CHECKPOINT_STATE_H

#include "util/sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)
#include "util/sep_shock_source_core.h"
#include <cstdint>
#include <string>
#include <vector>

namespace SEP { namespace PhysicsCheckpoint {

// Private incoming-wave reservoir caches are physical boundary conditions.
// They cannot be reconstructed from the current last cell after evolution.
// Values are E- [J] in integrated mode, W-(k) [J/m^3/bin] in spectral mode.
struct BoundaryCache {
  int lineCount = -1;
  std::uint64_t components = 0;
  std::vector<double> values;
  std::vector<std::vector<double>> spectra;
  std::vector<char> isSet;
};
struct AnalyticShockCache {
  bool initialized = false;
  double radiusM = 0.0;
  double minimumLineRadiusM = -1.0;
  Shock::TrajectoryState trajectory;
};
struct HostStepCache {
  bool adiabaticInitialized = false;
  bool shockHistoryInitialized = false;
  double previousShockRadiusM = 0.0;
};

// The driver owns this history. Capture/restore must be called only between
// complete native steps; output counters deliberately have a separate contract.
HostStepCache CaptureHostStepCache();
Transport::Status ValidateHostStepCache(const HostStepCache& cache);
void RestoreHostStepCache(const HostStepCache& cache);

}  // namespace PhysicsCheckpoint
namespace AlfvenTurbulence_Kolmogorov {
PhysicsCheckpoint::BoundaryCache CaptureIntegratedBoundaryCache();
Transport::Status ValidateIntegratedBoundaryCache(const PhysicsCheckpoint::BoundaryCache& cache);
void RestoreIntegratedBoundaryCache(PhysicsCheckpoint::BoundaryCache* cache);
namespace WaveNumberResolved {
PhysicsCheckpoint::BoundaryCache CaptureSpectralBoundaryCache();
Transport::Status ValidateSpectralBoundaryCache(const PhysicsCheckpoint::BoundaryCache& cache);
void RestoreSpectralBoundaryCache(PhysicsCheckpoint::BoundaryCache* cache);
}
}
namespace ParticleSource { namespace ShockWave { namespace Tenishev2005 {
Transport::Status CaptureAnalyticShockCache(PhysicsCheckpoint::AnalyticShockCache* cache);
Transport::Status ValidateAnalyticShockCache(const PhysicsCheckpoint::AnalyticShockCache& cache);
void RestoreAnalyticShockCache(PhysicsCheckpoint::AnalyticShockCache* cache);
} } }
}  // namespace SEP
#endif

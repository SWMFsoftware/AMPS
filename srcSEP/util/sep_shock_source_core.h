#ifndef SEP_UTIL_SEP_SHOCK_SOURCE_CORE_H
#define SEP_UTIL_SEP_SHOCK_SOURCE_CORE_H

#include "sep_transport_common.h"

#include <cstddef>
#include <string>
#include <vector>

namespace SEP {
namespace Shock {

// Dependency-free SI representation used by both the analytical trajectory and
// the field-line intersection adapter.  Keeping geometry out of PIC types makes
// the edge cases testable without initializing an AMPS mesh.
struct Vector3 {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

struct SpeedKnot {
  double radiusM = 0.0;
  double speedMPerS = 0.0;
};

struct TrajectoryConfiguration {
  double launchEpochS = 0.0;
  double launchRadiusM = 0.0;
  std::vector<SpeedKnot> knots;
};

struct TrajectoryState {
  TrajectoryConfiguration configuration;
  double epochS = 0.0;
  double radiusM = 0.0;
};

// The trajectory integrates dr/dt=v(r) analytically on every linear speed
// interval.  It is therefore continuous at knots and independent of how a
// caller partitions a global timestep.
Transport::Status ValidateTrajectoryConfiguration(
    const TrajectoryConfiguration& configuration);
Transport::ScalarResult SpeedAtRadius(
    const TrajectoryConfiguration& configuration, double radiusM);
Transport::Status StateAtEpoch(const TrajectoryConfiguration& configuration,
                               double epochS, TrajectoryState* state);
Transport::Status AdvanceToEpoch(TrajectoryState* state, double epochS);
Transport::Status SerializeTrajectory(const TrajectoryState& state,
                                      std::string* text);
Transport::Status DeserializeTrajectory(const std::string& text,
                                        TrajectoryState* state);

enum class CrossingOrientation { Inward, Outward, Tangent };
enum class IntersectionPolicy { First, FirstOutward, All };
enum class IntersectionStatus { Ok, NoIntersection, Ambiguous, InvalidGeometry };

struct Intersection {
  std::size_t segment = 0;
  double fraction = 0.0;
  double arcLengthM = 0.0;
  Vector3 positionM;
  CrossingOrientation orientation = CrossingOrientation::Tangent;
};

struct IntersectionResult {
  IntersectionStatus status = IntersectionStatus::NoIntersection;
  std::string message;
  std::vector<Intersection> intersections;
};

// Intersect an arbitrary polyline with a sphere centered at the origin.  The
// tolerance is a physical distance in metres.  Endpoint duplicates shared by
// adjacent segments are coalesced, while genuine multiple crossings are kept.
IntersectionResult IntersectSphere(const std::vector<Vector3>& vertices,
                                   double radiusM, double toleranceM,
                                   IntersectionPolicy policy);

struct TurbulenceSourceInput {
  double sweptVolumeM3 = 0.0;
  double upstreamMassDensityKgPerM3 = 0.0;
  double shockNormalSpeedMPerS = 0.0;
  double upstreamNormalSpeedMPerS = 0.0;
  double efficiency = 0.0;
  double plusBranchFraction = 0.5;
};

struct TurbulenceSourceEnergy {
  Transport::Status status;
  double plusJ = 0.0;
  double minusJ = 0.0;
  double totalJ = 0.0;
};

// Converts the upstream normal-relative kinetic energy in the swept volume to
// waves using E_wave = eta * (rho V / 2) * u_rel^2.  A zero relative speed is a
// valid zero source; negative or non-finite physical inputs are rejected.
TurbulenceSourceEnergy ComputeTurbulenceSource(
    const TurbulenceSourceInput& input);

}  // namespace Shock
}  // namespace SEP

#endif  // SEP_UTIL_SEP_SHOCK_SOURCE_CORE_H

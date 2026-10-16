#ifndef SRCMOON_MOON_SURFACE_ENVIRONMENT_H
#define SRCMOON_MOON_SURFACE_ENVIRONMENT_H

#include "LunarSurface.h"

#include <array>
#include <cstdint>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace Moon {
namespace SurfaceEnvironment {

struct RayHit {
  bool hit = false;
  double distanceM = std::numeric_limits<double>::infinity();
  std::size_t face = std::numeric_limits<std::size_t>::max();
  std::array<double, 3> positionM{{0.0, 0.0, 0.0}};
  std::array<double, 3> outwardNormal{{0.0, 0.0, 0.0}};
};

// Reusable body-fixed terrain query service.  Its BVH is derived only from the
// production triangulation; D02 illumination products never enter this class,
// preserving them as independent I15/I20 evidence.
class TerrainRayTracer {
 public:
  bool Initialize(const Surface::Triangulation& triangulation,
                  std::string* error);
  bool IntersectRay(const std::array<double, 3>& originM,
                    const std::array<double, 3>& direction,
                    double minimumDistanceM, double maximumDistanceM,
                    std::size_t excludedFace, RayHit* hit,
                    std::string* error) const;
  bool IsSunVisible(const std::array<double, 3>& surfacePositionM,
                    const std::array<double, 3>& outwardNormal,
                    const std::array<double, 3>& unitSunDirection,
                    double rayOriginOffsetM, std::size_t sourceFace,
                    bool* visible, std::string* error) const;
  const Surface::Triangulation& triangulation() const { return mesh_; }

 private:
  struct Node {
    std::array<double, 3> minimum{{0.0, 0.0, 0.0}};
    std::array<double, 3> maximum{{0.0, 0.0, 0.0}};
    std::size_t begin = 0;
    std::size_t end = 0;
    int left = -1;
    int right = -1;
  };

  int BuildNode(std::size_t begin, std::size_t end);
  Surface::Triangulation mesh_;
  std::vector<std::size_t> faceOrder_;
  std::vector<Node> nodes_;
};

// Process-wide production terrain state.  PrepareProductionSurface installs
// the exact triangulation generated from D01 on every MPI rank.  These hooks
// match PIC::Mover's optional moving-surface ABI and therefore exercise the
// normal particle mover rather than a test-only propagator.
bool InitializeProductionTerrain(const Surface::Triangulation& triangulation,
                                 const Runtime::Configuration& configuration,
                                 std::string* error);
bool ProductionTerrainInitialized();
const TerrainRayTracer* ProductionTerrain();
void BeginProductionParticleStep(double ephemerisTimeBeginS);
bool FindProductionMovingSurfaceIntersection(long int particle,
    const double* segmentPositionLsoM,const double* segmentVelocityLsoMPerS,
    const double* segmentAccelerationLsoMPerS2,double elapsedFlightTimeS,
    double maximumFlightTimeS,double& hitFlightTimeS,
    double* hitPositionLsoM,double* outwardNormalLso,
    double* surfaceVelocityLsoMPerS);
bool LastProductionCollisionFrame(double* outwardNormalBodyFixed,
                                  double* lsoToBodyStateTransform36,
                                  double* bodyToLsoStateTransform36);

// Monotonic process-local diagnostic count of geometrically resolved moving
// terrain intersections.  It is incremented inside the same production finder
// invoked by PIC::Mover and is never used to alter motion or physics.  Linked
// I14 snapshots it before/after a mover fixture to prove that the normal mover
// traversed the production LOLA callback rather than a test-only ray tracer.
std::uint64_t ProductionCollisionCount();

}  // namespace SurfaceEnvironment
}  // namespace Moon

#endif  // SRCMOON_MOON_SURFACE_ENVIRONMENT_H

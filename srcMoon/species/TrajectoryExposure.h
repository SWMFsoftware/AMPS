#ifndef SRCMOON_SPECIES_TRAJECTORY_EXPOSURE_H
#define SRCMOON_SPECIES_TRAJECTORY_EXPOSURE_H

#include "PiecewiseHazard.h"
#include <array>
#include <functional>
#include <string>
#include <vector>

namespace Moon { namespace Species { namespace TrajectoryExposure {
using Vector=std::array<double,3>;
// The production mover's accepted quadratic POSITION extension, in metres
// on one frozen Cartesian frame. This does not approximate the velocity
// quadrature or invoke a second copy of any production force law.
struct Arc { Vector positionM{},velocityMS{},accelerationMS2{}; };
struct Sphere { Vector centreM{};double radiusM=0; };
struct Face { Vector verticesM[3]; };
struct PhotonChannel { std::size_t reaction=0;double rateAtReferenceS=0; };
struct ConstantChannel { std::size_t reaction=0;double rateS=0; };
using RayVisibility=std::function<bool(const Vector&,const Vector&,bool*,std::string*)>;
struct Field {
  std::string frame;             // explicit common axes, e.g. MOON_ME_DE421
  double frozenGeometryEpochS=0; // external ET; not a replacement particle clock
  Vector sunPositionM{};
  double referenceDistanceM=0;   // source rate reference distance, explicitly SI
  std::vector<Sphere> opaqueSpheres;
  // Actual production mesh faces on the SAME axes. The visibility callback
  // must use that mesh's production BVH. Both are supplied together or absent.
  std::vector<Face> terrainFaces;
  RayVisibility terrainVisibility;
};
struct Controls {
  double rootTimeToleranceS=0;
  long double absoluteDepthTolerance=0,relativeDepthTolerance=0;
  std::size_t maximumSegments=0,maximumGeometryIntervals=0;
};
struct History {
  std::vector<PiecewiseHazard::Segment> segments;
  long double interpolationDepthErrorBound=0;
  long double shadowBoundaryDepthErrorBound=0;
  std::size_t silhouetteRoots=0;
  std::size_t geometryIntervalsVisited=0;
};
// Enumerate point-Sun silhouette transitions before constructing the smooth
// inverse-square approximation. No endpoint-only shadow sampling, implicit
// time/distance tolerance or empirical rate correction is used. The two
// reported error bounds must fit the declared TOTAL depth budget. Resource,
// geometry and representability failures leave output unchanged.
bool BuildPhotonHistory(const Arc& arc,double durationS,const Field& field,
    const std::vector<PhotonChannel>& photons,
    const std::vector<ConstantChannel>& constants,const Controls& controls,
    History* output,std::string* error);
}}}
#endif

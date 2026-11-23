#ifndef SRCSEP_UTIL_SEP_REDUCED_GEOMETRY_H
#define SRCSEP_UTIL_SEP_REDUCED_GEOMETRY_H
#include "../adapters/reduced_shock_background_adapter.h"
#include <functional>
#include <vector>

namespace SEP { namespace ReducedShock {

// These are numerical work/accuracy controls, in SI; they never modify the
// caller's physical interval, PIC dt, sample count or background cadence.
struct GeometryControls {
  double positionToleranceM=0;
  double maximumTimeSubstepS=0;
  double maximumArcSubstepM=0;
  std::size_t maximumAttempts=0;
};
using VelocityQuery=std::function<Core::Result<Position>(const Position&,double)>;

// Controlled RK4 with step doubling on an endpoint-preserving bisection mesh.
// The local estimate retains its duration-proportional share of the supplied
// global position tolerance. Bisection avoids an artificially short final
// remainder; maximum-step and accuracy subdivisions share the fixed work bound.
// Signed intervals also support the inverse
// flow map needed to query the actual background history of newly entering
// material. This is a mathematical query, never reverse particle stepping.
Core::Result<Position> IntegrateTrajectory(const VelocityQuery&,const Position&,
    double departureS,double arrivalS,const GeometryControls&);
Core::Result<Position> MoveVertex(const IntervalView&,const Position&,
    double departureS,double arrivalS,const GeometryControls&);
Core::Result<Position> TraceMagnetic(const IntervalView&,const Position&,
    const Position& outwardChord,double arcLengthM,double epochS,
    const GeometryControls&);

// The core's magnetic sector, not magnetic polarity, orients increasing line
// index outward. All returned points/state come from the selected authority.
Core::Result<std::vector<Position>> BuildInitialLine(const IntervalView&,
    const Position& initial,double lengthM,std::size_t pointCount,double epochS,
    const GeometryControls&);

struct Connection {
  std::size_t segment=0;
  double fraction=0,arcLengthM=0,surfaceResidualM=0;
  bool tangent=false;
  FrontSample front;
};
// All genuine supported intersections are retained. Interior shared endpoints
// use the next segment at fraction zero, matching native floor(S) ownership.
// A finite spatial tolerance grades numerical surface/line residuals; it is
// not an obliquity floor or a substitute for missing physical front support.
Core::Result<std::vector<Connection>> ConnectLine(const IntervalView&,
    const std::vector<Position>&,double epochS,double spatialToleranceM);
Core::Result<std::vector<Connection>> ConnectMaterialLine(const IntervalView&,
    const std::vector<Position>& departureVertices,double epochS,
    double spatialToleranceM,const GeometryControls&);

} }
#endif

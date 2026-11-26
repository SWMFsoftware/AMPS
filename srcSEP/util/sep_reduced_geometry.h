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
// firstArcM>0 selects geometric arcs a q^k (k=0..N-2) whose sum is lengthM;
// the default zero keeps the historical uniform arc lengthM/(N-1).
Core::Result<std::vector<Position>> BuildInitialLine(const IntervalView&,
    const Position& initial,double lengthM,std::size_t pointCount,double epochS,
    const GeometryControls&,double firstArcM=0.0);
// The arc law is SEP::ReducedShock::GeometricArcs (sep_geometric_spacing.h).

// Follow the core magnetic field line through 'anchor' SUNWARD (opposite the
// core's outward sector orientation) to the sphere |x|=radiusM at epochS and
// return that point. Each step traces an arc equal to the remaining radial
// distance d=|x|-radiusM; since |dr/ds|<=1 on any curve, such a step cannot
// cross the sphere beyond the trace tolerance, and the iteration stops when
// d<=radiusToleranceM. 'trace' controls each magnetic leg (its tolerance is
// the leg's absolute position accuracy). A field line that does not approach
// the Sun (the radius fails to decrease) or does not converge within
// maximumSteps is an error, never a projected or radially substituted point.
Core::Result<Position> TraceSunwardToRadius(const IntervalView&,const Position& anchor,
    double radiusM,double epochS,const GeometryControls& trace,double radiusToleranceM,
    std::size_t maximumSteps=100000);

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
// 'points' may be the contiguous vertex range [firstIndex, firstIndex+size)
// of a line with lineVertexCount vertices (0: the points are the whole line);
// returned segments are then whole-line indices and the shared-endpoint rule
// uses the whole line, so a range that contains every root gives the same
// connections as the whole line.
Core::Result<std::vector<Connection>> ConnectLine(const IntervalView&,
    const std::vector<Position>&,double epochS,double spatialToleranceM,
    std::size_t firstIndex=0,std::size_t lineVertexCount=0);

// Smallest contiguous vertex range [first, first+count) containing every
// segment of the departure line that can meet a sphere about 'center' with a
// radius in [innerRadiusM, outerRadiusM] after each vertex moves by at most
// reachM. A moved chord spans radii within [min(r_i,r_i+1) - 2 reach - l/2,
// max(r_i,r_i+1) + reach] (triangle inequality and convexity of |x|), so any
// segment outside that band is proved disjoint. count==0: no candidate.
struct VertexRange {std::size_t first=0,count=0;};
// Per-segment radial extents of a departure line about 'center', reduced to a
// nondecreasing prefix maximum of max(r_i,r_i+1) and a nondecreasing suffix
// minimum of min(r_i,r_i+1)-l_i/2, so a candidate query is two binary searches.
// The resulting range is a superset of the candidate segments (exact when the
// line's radius is monotone); extra segments only cost advection, never roots.
struct RadialBounds {
  Position center{};
  std::vector<double> prefixMaximumHigh,suffixMinimumLow;
  std::vector<double> arc; // departure-chord arc length at each vertex [m]
  bool finite=false;
};
RadialBounds BuildRadialBounds(const std::vector<Position>& departure,const Position& center);
VertexRange CandidateVertices(const RadialBounds&,double innerRadiusM,double outerRadiusM,double reachM);
VertexRange CandidateVertices(const std::vector<Position>& departure,
    const Position& center,double innerRadiusM,double outerRadiusM,double reachM);

// Connections of the material line at epochS. Only the candidate range for
// the generating sphere (with the core speed bound times the elapsed time as
// reach) is advected; the remaining segments cannot meet the front. Each
// returned arcLengthM is the root's material (departure-chord) coordinate,
// which equals the geometric arc for the departure epoch.
Core::Result<std::vector<Connection>> ConnectMaterialLine(const IntervalView&,
    const std::vector<Position>& departureVertices,double epochS,
    double spatialToleranceM,const GeometryControls&,const RadialBounds* bounds=nullptr);

} }
#endif

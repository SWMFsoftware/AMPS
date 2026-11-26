// ============================================================================
// srcSEP3D/mesh/observer_flux_tube.h
//
// Observer flux tube (mesh.active_region.mode = observer-flux-tube).
//
// WHY THIS MODULE EXISTS:
//   The analytic `parker-tube` corridor follows one constant-speed Parker
//   spiral with a prescribed cone or constant width.  A run that is compared
//   with real observers needs the domain to be the magnetic flux tube that
//   connects those observers to the Sun in the run's own background field
//   (for the reduced shock-front event: PFSS plus the accelerating open-tube
//   wind), wide enough that every observer's collection sphere stays inside
//   it for as long as the observer has to be sampled.  This module traces
//   that tube.  Its lines are used twice:
//     - by the mesh (mesh_model.h, ObserverTubeGeometry): the active-block
//       mask and the distance-graded refinement;
//     - by the shock source (adapters/tube_line_source.h,
//       source.shock.placement = release-surface-on-tube-lines): each kept
//       line is one flux element of the 1-D-equivalent shock source.
//
// FRAME AND COORDINATES:
//   All positions are Sun-centred coordinates [m] in the computation frame
//   (core/computation_frame.h); the caller removes the mesh coordinate
//   origin, which the run configuration currently requires to be zero
//   (runtime/run_configuration.cpp).  The background field must be static in
//   that frame.  For the reduced PFSS/Parker ambient this is the corotating
//   frame, in which field lines do not move and only the observers do; the
//   caller evaluates the field at the frame's reference epoch and refuses an
//   HCI run with a rotating ambient (main_lib.cpp,
//   BuildObserverFluxTubeBeforeMesh).
//
// HOW THE OBSERVERS SET THE TUBE WIDTH:
//   A tube observer (listed in mesh.active_region.observers) has
//     - a trajectory (output/observer_trajectory.h): an astronomical body
//       from a checksummed ephemeris table or from SPICE, a straight HCI line
//       (moving-cartesian: a point x0 and a velocity v), or a point fixed in
//       the computation frame;
//     - a collection radius rho [m];
//     - a stay interval [tube_stay_begin_s, tube_stay_end_s] of model time
//       during which its collection sphere must be inside the simulated tube.
//   The caller samples the stay interval at the cadence
//   mesh.active_region.time_spacing_s (default 0.25 rho / (|Omega| r + |v|),
//   so consecutive collection spheres overlap) and passes the frame
//   positions c_k of the sphere centres as a TubeObserverPath.  The tube is
//   the union of the flux elements met by all those spheres.  A frame-fixed
//   observer therefore gives a tube about 2 (rho + s) wide.  An inertial
//   observer (Earth, SOHO, ACE) moves through the corotating frame at nearly
//   the corotation speed |Omega| r, so a longer stay sweeps its sphere across
//   more field lines and widens the tube along that motion.
//
// CONSTRUCTION (BuildObserverFluxTube):
//  1. Reference sphere r = R0: the radius of the first tube observer at the
//     start of its stay.  A colatitude-longitude lattice of spacing s
//     (mesh.active_region.seed_spacing_m; default half the requested cell
//     size at that observer) covers it.  Ring i spans the colatitudes
//     [i, i + 1] * dTheta, dTheta = pi / ceil(pi R0 / s), and has
//     n_i = round(2 pi R0 sin(theta_i) / s) cells of equal longitude width
//     (theta_i the ring's mid colatitude).  The centre of each cell seeds one
//     field line.  The cell's exact spherical area
//         dA = R0^2 (cos(i dTheta) - cos((i + 1) dTheta)) * 2 pi / n_i
//     gives the line's flux element  Phi = |B_r| dA  [Wb] (B_r = B . r_hat at
//     the cell centre, midpoint rule).  The cells are disjoint, so the kept
//     lines partition the tube's magnetic flux without overlap.
//  2. Candidates.  The field line through every sample c_k is traced to R0,
//     where it lands at q_k.  By flux conservation (div B = 0) a small
//     cross-section of area A at c_k has area A |B(c_k)| / |B(q_k)| at q_k,
//     so transverse lengths scale by m_k = sqrt(|B(c_k)| / |B(q_k)|).  The
//     lattice cells whose centres lie within the great-circle angle
//     alpha_k = (rho m_k + 2 s) / R0 of q_k are candidates.
//  3. Selection.  Each candidate line is traced sunward to the inner radius
//     and outward to the outer radius.  It is kept if, at the radius |c_k| of
//     some sample, it passes within rho + s_k of c_k, where s_k is the line's
//     local spacing there (step 4): the flux element around the line then
//     reaches the collection sphere.
//  4. Local spacing.  The element's cross-section at x is A(x) = Phi / |B(x)|
//     (flux conservation along the element), so its transverse size is
//         s(x) = sqrt(Phi / |B(x)|).
//     Each kept line stores its points and s(x) at every point.
//
// HOW THE MESH USES THE LINES (main_lib.cpp and mesh_model.cpp):
//   - mask: one capsule per chord of every kept line, with radius
//     max(s_i, s_{i+1}) + chord / 2 (+ a roundoff guard).  The s term covers
//     the flux between neighbouring lines and the chord / 2 term covers the
//     curved line between its two points.  A bounding-volume hierarchy (BVH)
//     over the capsule boxes classifies the AMR leaves (mesh_model.cpp,
//     CapsuleBvh);
//   - refinement: the requested cell size is the radial target inside the
//     tube and about twice the distance d to the tube outside it
//     (RequestedCellSizeM), so only the blocks the tube crosses refine;
//   - coverage: after block allocation every sampled collection sphere of
//     every observer (tube observers over their stay intervals, other
//     observers over the whole run) must lie in active blocks, or the run
//     stops (main_lib.cpp, VerifyObserverTubeCoverage).
//
// FIELD-LINE INTEGRATION (TraceFieldLineToRadius):
//   Classical RK4 in arc length s for dx/ds = sign * sector * B / |B|, with
//   step ds = min(trace_step_m, 0.02 r).  The last step is shortened by
//   bisection so that the final point lies on the target sphere.  A closed
//   field line (sector 0), a change of sector, a zero or nonfinite field, a
//   line whose radius is not monotone, or a line that does not reach the
//   sphere within the step budget is an error.  No point is ever projected
//   onto the sphere or substituted to hide such a failure.  (The final
//   point's radial roundoff residual is removed exactly; see the .cpp.)
//
// ASSUMPTIONS AND LIMITS:
//   - Open field lines that are monotone in r between the inner and the
//     outer sphere (the Parker regime of the reduced ambient).  A candidate
//     line that enters closed field or turns back in r stops the build.
//   - The candidate window of step 2 scales the cross-section isotropically
//     (one factor m_k for both transverse directions).  When the observers
//     lie at radii far from R0 and the real cross-section is elongated, the
//     window rests on the 2 s margin.  An observer sphere that ends up
//     uncovered is caught by the coverage check, not silently accepted.
//   - All MPI ranks must trace identical lines.  The caller compares the
//     digest of all line points across ranks.
//
// The module is independent of AMPS, PIC and MPI.
// Documentation: srcSEP3D/CONFIGURATION.md (the mesh.active_region.* and
// observer.ID.tube_stay_* keys) and srcSEP3D/test/event_2013_04_11/README.md,
// section 2.3.
// ============================================================================

#ifndef SEP3D_MESH_OBSERVER_FLUX_TUBE_H
#define SEP3D_MESH_OBSERVER_FLUX_TUBE_H

#include "../core/sep3d_types.h"

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

namespace SEP3D {
namespace Mesh {

// Static background field of the computation frame.
//   argument 1  x [m]: a Sun-centred point of the computation frame;
//   argument 2  out: B(x) [T], in frame components;
//   argument 3  out: the open-field sector at x.  +1 or -1 on open field,
//               chosen so that sector * B points away from the Sun along the
//               line; 0 on closed field.
// Returns a non-OK status when the field is not available at x (for example
// outside the provider's domain).  The tube functions return that status
// unchanged.  The function must be deterministic: every MPI rank traces the
// same lines from it.
using FrameFieldFunction =
    std::function<Core::Status(const Core::Vec3&, Core::Vec3*, int*)>;

// One tube observer, sampled over its stay interval by the caller
// (main_lib.cpp, BuildObserverFluxTubeBeforeMesh).  BuildObserverFluxTube
// requires: at least one sample, equal lengths of timesS and centresM,
// rho > 0, and every sphere strictly inside the radial shell
// (inner radius, outer radius).
struct TubeObserverPath {
  std::string id;                      // observer ID (diagnostics only)
  double collectionRadiusM = 0.0;      // rho [m], radius of the collection sphere
  std::vector<double> timesS;          // model times [s] of the samples of the stay interval
  std::vector<Core::Vec3> centresM;    // Sun-centred frame positions [m] of the sphere centre at those times
};

// Numerical controls of the tube construction (from mesh.active_region.*).
struct ObserverTubeControls {
  double innerRadiusM = 0.0;   // transport inner sphere [m]; every line starts on it
  double outerRadiusM = 0.0;   // transport outer sphere [m]; every line ends on it
  double traceStepM = 0.0;     // upper bound of the arc step [m] (the step is min(this, 0.02 r))
  double seedSpacingM = 0.0;   // lattice spacing s [m] on the reference sphere (> 0 here; the caller resolves 0)
  std::size_t maximumStepsPerLine = 200000;  // RK4 step budget of one traced leg
};

// One kept field line = one flux element of the tube.
struct TracedTubeLine {
  // Points [m] from the inner sphere outward.  The first and last points lie
  // on the inner and outer spheres to roundoff, and the seed point appears
  // once.  The radius increases strictly along the points, because the trace
  // rejects a line that is not monotone in r (PointAtRadius relies on this).
  std::vector<Core::Vec3> pointsM;
  std::vector<double> localSpacingM;   // s(x) = sqrt(Phi/|B(x)|) [m] at each point (same length as pointsM)
  double fluxWb = 0.0;                 // Phi = |B_r| dA [Wb] of the seed's lattice cell
  int sector = 0;                      // +1/-1 open-field sector of the line (never 0)
  Core::Vec3 referencePointM;          // seed: centre of the lattice cell on r = R0 [m]
};

// The traced tube.
struct ObserverFluxTube {
  double referenceRadiusM = 0.0;       // R0 [m]
  double seedSpacingM = 0.0;           // s [m]
  // The kept lines in lattice order (ring index, then longitude index), so
  // every rank holds the same lines in the same order.
  std::vector<TracedTubeLine> lines;
  std::size_t candidateLines = 0;      // lattice cells tested in step 3
  std::size_t observerSamples = 0;     // observer samples over all paths
  // FNV-1a-style hash of the bytes of every line point, in line order.  It
  // identifies the tube; the caller compares it across MPI ranks.
  std::uint64_t digest = 0;
};

// Trace one field line from `start` to the sphere |x| = targetRadiusM, along
// dx/ds = sign * sector * B/|B| (sign = +1 outward, -1 sunward).  The sector
// is taken from the field at `start` and must stay the same along the line.
//   field         the static frame field;
//   start [m]     the first point (always points->front());
//   sign          +1 or -1;
//   targetRadiusM the radius [m] to reach;
//   traceStepM    the arc-step bound [m] (the step is min(this, 0.02 r));
//   maximumSteps  the RK4 step budget;
//   points        out: the line points from `start` to the landing point,
//                 which lies on the sphere to roundoff.  When `start` is
//                 already at or beyond the target radius in the direction
//                 `sign`, the output is just {start}.
// Errors (InvalidInput): an invalid request, closed field, a change of
// sector, a zero or nonfinite |B|, a step that does not move monotonically
// in r in the direction `sign`, or the budget running out before the
// sphere.  A non-OK status of `field` is returned unchanged.  On error the
// contents of *points are unspecified.
Core::Status TraceFieldLineToRadius(const FrameFieldFunction& field,
                                    const Core::Vec3& start, int sign,
                                    double targetRadiusM, double traceStepM,
                                    std::size_t maximumSteps,
                                    std::vector<Core::Vec3>* points);

// Build the observer flux tube (steps 1-4 of the header).
//   paths     the tube observers' sampled stay intervals; the first sample
//             of paths.front() sets R0;
//   controls  inner/outer radii, trace step, seed spacing s, step budget;
//   field     the static frame field;
//   tube      out: the tube.  It is written only on success.
// Errors: InvalidInput for missing or invalid inputs, an observer sphere
// outside the radial shell, a centre line or a candidate line that cannot be
// traced, or a tube with no kept line.  A non-OK field status is returned
// unchanged.  The result is deterministic: it depends only on the inputs, and
// the candidate set is iterated in sorted lattice order.
Core::Status BuildObserverFluxTube(const std::vector<TubeObserverPath>& paths,
                                   const ObserverTubeControls& controls,
                                   const FrameFieldFunction& field,
                                   ObserverFluxTube* tube);

// The point of a traced line at heliocentric radius r [m], and optionally
// the local spacing there.  The two points whose radii bracket r are found
// by bisection; both values are interpolated linearly, with the fraction
// f = (r - r_a) / (r_b - r_a) along the chord.  Returns false, and leaves the
// outputs untouched, when the line has fewer than two points or r lies
// outside [r_first, r_last].  It requires a radius that increases along the
// line, which TraceFieldLineToRadius guarantees.
bool PointAtRadius(const TracedTubeLine& line, double radiusM,
                   Core::Vec3* pointM, double* localSpacingM = nullptr);

}  // namespace Mesh
}  // namespace SEP3D

#endif  // SEP3D_MESH_OBSERVER_FLUX_TUBE_H

// ============================================================================
// srcSEP3D/adapters/tube_line_source.h
//
// 1-D-equivalent shock source on the lines of an observer flux tube
// (source.shock.placement = release-surface-on-tube-lines).
//
// WHY THIS MODULE EXISTS:
//   The 3-D 2013-04-11 test must inject the SAME physical source as the 1-D
//   srcSEP test (srcSEP/test/event_2013_04_11), so that the three transport
//   laws (Parker, focused diffusion, focused scattering) can be compared
//   between the two codes with only the transport geometry changed.  The 1-D
//   source releases protons from one field line of flux Phi at a release
//   sphere upstream of the reduced shock (srcSEP/reduced_shock_particles.cpp,
//   srcSEP/util/sep_reduced_source.cpp).  In 3-D the active domain is the
//   observer flux tube (mesh/observer_flux_tube.h): a bundle of traced field
//   lines, each the axis of one flux element of flux Phi_l.  This module
//   applies the 1-D source to every one of those elements and samples the
//   result into macroparticles.
//
// RATE OF ONE FLUX ELEMENT:
//   Every traced tube line l is one flux element of magnetic flux Phi_l [Wb].
//   At a birth epoch t the line meets the shock at x_sh (a sign change of
//   ShockSurface::Level along the line); if the front is an accepted fast
//   shock there, the line carries the 1-D srcSEP net proton rate
//   (NetProtonRate in srcSEP/util/sep_reduced_source.cpp)
//
//     R_l(t) = eta n1 w1 Phi_l / |B1 . n|      [protons/s],
//
//   eta = injection_efficiency (dimensionless yield), n1 the upstream proton
//   density [m^-3], w1 = V_n - U1.n the shock-frame normal inflow [m/s], B1
//   the upstream field [T] and n the outward unit normal.  Phi_l / |B1.n| is
//   the area [m^2] the element intercepts on the front, so n1 w1 times that
//   area is the proton flux [1/s] entering the shock through the element.
//
// RELEASE:
//   The particle is released where the line crosses the release surface
//   (ShockSurface::ReleaseLevel = 0).  For the heliocentric sphere that is
//   the sphere of radius (1 + f) R_sh, f =
//   source.shock.release_distance_shock_radius_fraction (0.3 in the
//   2013-04-11 decks, i.e. 1.3 R_sh), which is the 1-D release distance
//   L_ref = f R_sh.  Each particle receives
//   * a momentum from the local-compression DSA spectrum of the crossing,
//     dN/dp ~ p^(2-q), q = 3r/(r-1), r the compression ratio at x_sh,
//     between the source energy limits; either sampled directly
//     (constant-statistical-weight) or from the log-uniform proposal
//     g(p) = 1/(p ln(pmax/pmin)) with the exact importance weight
//     w = f_DSA(p)/g(p) (log-uniform-momentum-importance).  Both use the
//     same core DSA law and the same weight helper as the 1-D source;
//   * a pitch-angle cosine from the positive-normal-flux law through the
//     release surface (BuildFocusedEscapeDistribution, the 1-D focused law)
//     for the focused movers, and an isotropic mu for Parker (the 1-D Parker
//     mover likewise stores the isotropic projection of the source);
//   * a transverse position uniform over the element's footprint on the
//     release surface: a square of side s in the plane normal to B, with s
//     the local line spacing sqrt(Phi/|B|), projected along B onto the
//     sphere.  The 1-D source has no transverse extent; this spread only
//     distributes the element's particles over the cross-section the element
//     occupies, so they do not all start on the traced line.
//
// TIME SAMPLING AND WEIGHTS:
//   The N = samples_per_step samples of one step [t0, t0 + dt] are stratified
//   in time: sample j is born at t_j = t0 + (j + u_j) dt / N, u_j uniform on
//   (0,1), exactly as the 1-D source.  At each sample one line is chosen with
//   probability R_l(t_j) / sum_l R_l(t_j); the particle then represents
//
//     W_j = w (sum_l R_l(t_j)) dt / N      physical protons,
//
//   w the importance weight (1 for direct sampling).  The expected number
//   contributed by line l is sum_j P(l) E[w] (sum R) dt/N = sum_j R_l(t_j) dt/N
//   (E[w] = 1 under the proposal), the 1-D stratified estimate of the time
//   integral of R_l for that element.  The expectation of the whole source is
//   therefore the 1-D source summed over the tube's flux elements.
//   srcSEP3D/test/event_2013_04_11/README.md section 2.4 compares the
//   injected protons per weber with the 1-D line.
//
// FRAME:
//   Positions, fields and velocities are in the run's computation frame
//   (core/computation_frame.h).  With the corotating frame used by the
//   2013-04-11 test the reduced ambient is static and the traced lines are
//   fixed curves, so the lines are not advected between births (the 1-D line
//   turns rigidly with the Sun, which is the same set of field lines).  The
//   plasma velocity passed in by `plasma` is the frame velocity
//   U' = U - Omega x x.  The shock provider converts to its own frame
//   internally (adapters/shock_surface.h).
//
// ASSUMPTIONS AND LIMITS:
//   * The lines run from the inner sphere outward with monotonically
//     increasing radius (the open Parker-regime lines the tube admits), and
//     the shock level is monotone along them; the crossing search is a
//     bisection that relies on this.
//   * The release-surface normal and the transverse spread use the
//     heliocentric radial direction x/|x| and re-project onto the sphere
//     |x| = const.  This is exact for the heliocentric-sphere front, whose
//     release surface is a heliocentric sphere; the run configuration
//     requires the coordinate origin at the Sun.
//   * One line is selected per sample, so the number of macroparticles per
//     step is N per species, independent of the number of tube lines.
//
// DETERMINISM AND MPI:
//   All draws are keyed by (campaign seed, step, sample, species, purpose)
//   (transport/keyed_random.h, RandomPurpose::TubeSource*), so every MPI
//   rank generates the same candidates; the caller (main_lib.cpp,
//   InjectTubeLineParticles) lets exactly one rank, the owner of the AMR
//   leaf, allocate each.
//
// The module is independent of AMPS, PIC and MPI.
//
// DOCUMENTATION:
//   srcSEP3D/CONFIGURATION.md (keys source.shock.placement,
//   release_distance_shock_radius_fraction, returning_particles),
//   srcSEP3D/README.md (3-D real-event test) and
//   srcSEP3D/test/event_2013_04_11/README.md section 2.4.
//   Unit tests: test/event_2013_04_11/test_event3d_units.cpp.
// ============================================================================

#ifndef SEP3D_ADAPTERS_TUBE_LINE_SOURCE_H
#define SEP3D_ADAPTERS_TUBE_LINE_SOURCE_H

#include "../mesh/observer_flux_tube.h"
#include "shock_surface.h"
#include "source_runtime.h"

#include <cstdint>
#include <functional>
#include <vector>

namespace SEP3D {
namespace Adapters {

// Frame field and plasma velocity at a frame point and epoch:
//   status = plasma(x [m, frame], t [s, model time], &B [T, frame vector],
//                   &U [m/s, frame velocity U' = U - Omega x x]).
// The caller binds it to the background provider
// (ShockFrontBackgroundAdapter::AmbientAt); a failed status aborts the batch.
using FramePlasmaFunction = std::function<Core::Status(
    const Core::Vec3&, double, Core::Vec3*, Core::Vec3*)>;

// One step's source request for one species.  Every pointer and the plasma
// function must stay valid for the duration of GenerateTubeLineSource.
struct TubeLineSourceRequest {
  const Mesh::ObserverFluxTube* tube = nullptr;   // traced lines (flux elements)
  const ShockSurface* shock = nullptr;            // front, release surface, front state
  FramePlasmaFunction plasma;                     // B and U' at the release point
  int species = 0;                 // AMPS species index (also keys the streams)
  double massKg = 0.0;             // rest mass [kg]
  double injectionEfficiency = 0.0;   // eta in R_l (> 0)
  double releaseFraction = 0.0;    // f: release at (1+f) times the front radius (> 0)
  double minimumEnergyJ = 0.0, maximumEnergyJ = 0.0;   // kinetic-energy band [J], 0 < min < max
  bool logUniformMomentum = false; // true: log-uniform proposal + exact weight
  bool focusedPitchLaw = true;     // false: isotropic (Parker)
  std::uint64_t samplesPerStep = 0;   // N macroparticles per step (> 0)
  std::uint64_t campaignSeed = 0;     // key of all random streams
  std::uint64_t step = 0;          // completed steps before this one
  double stepBeginS = 0.0;         // t0 [s, model time]
  double dtS = 0.0;                // step length [s] (> 0)
};

// Result of one batch.  On a failed status `particles` is incomplete and
// must not be allocated.
struct TubeLineSourceResult {
  Core::Status status;
  std::vector<InjectedParticle> particles;   // all candidates (every rank)
  // Diagnostics of the first sample (j = 0), whose epoch lies in the first
  // time stratum [t0, t0 + dt/N), not exactly at t0:
  double totalRateBeginPerS = 0.0;           // sum R_l at the first sample's epoch [1/s]
  double representedParticles = 0.0;         // sum of particle weights (physical protons)
  std::uint64_t acceptedLines = 0;           // lines with R_l > 0 at the first sample's epoch
  std::uint64_t emptySamples = 0;            // samples with sum R = 0 (no particle produced)
};

// Net rate of one line at epoch t (R_l above, [protons/s]); zero when the
// line does not cross the front (its inner end is not downstream, i.e. the
// front is still below it, or its outer end is not upstream, i.e. the front
// has passed it) or when the front at the crossing is not an accepted shock.
//   line             a traced tube line, points ordered outward
//   shock            the front (frame interface)
//   epochS           model time [s]
//   injectionEfficiency  eta
//   ratePerS         output R_l [1/s]; set to 0 before any test
//   crossingM        optional: the crossing point on the front [m, frame]
//   front            optional: the provider's front state at the crossing
// `crossing` and `front` are written whenever the line crosses the front,
// also when the front there is not accepted.  Fails (InvalidInput) on a null
// output, a line of fewer than two points, a provider error, or an accepted
// front without a physical upstream state (B1 tangent to the front, zero
// density or non-positive inflow): as in 1-D, no obliquity floor is applied.
Core::Status TubeLineNetRate(const Mesh::TracedTubeLine& line,
                             const ShockSurface& shock, double epochS,
                             double injectionEfficiency, double* ratePerS,
                             Core::Vec3* crossingM = nullptr,
                             ShockFrontState* front = nullptr);

// Generate the N candidates of one step (see the header block).  Pure
// function of the request: identical on every rank.  Fails (InvalidInput)
// on an incomplete request or on any of the 1-D source's fail-closed
// conditions: a selected line with no release-surface crossing, a DSA law
// the core rejects, a non-finite importance weight or particle weight, or a
// source band that contains momenta with no outward escape through the
// release surface (focused movers).  The particles are returned in sample
// order, each with its birth phase encoded in remainingFirstStepFraction.
TubeLineSourceResult GenerateTubeLineSource(const TubeLineSourceRequest& request);

}  // namespace Adapters
}  // namespace SEP3D

#endif  // SEP3D_ADAPTERS_TUBE_LINE_SOURCE_H

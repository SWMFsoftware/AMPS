// ============================================================================
// srcSEP3D/output/observer_trajectory.h
//
// Observer trajectories: an observer's Sun-centred position at model time t.
//
// WHY THIS MODULE EXISTS:
//   In a corotating run (core/computation_frame.h) the field lines are
//   static and the real observers move: Earth, SOHO and ACE follow inertial
//   orbits, so in the corotating frame they drift across the Parker spiral
//   at nearly the corotation speed.  Two consumers need the same position
//   function:
//     - the observer flux tube (mesh/observer_flux_tube.h,
//       mesh.active_region.mode = observer-flux-tube), which is traced
//       around the observers' positions over their stay intervals;
//     - the observer sampling (output/observer_runtime.cpp,
//       BuildObserverDefinitions), which places each collection sphere at
//       every publication.
//   With one module for both, the tube is built from exactly the path that
//   is later sampled.
//
// KINDS (observer.ID.kind; srcSEP3D/CONFIGURATION.md):
//   fixed-cartesian, fixed-heliographic, field-connected, spherical-shell:
//       fixed in the computation frame.  In a corotating run that means
//       fixed on a field line, like a 1-D srcSEP observer.  StateHci returns
//       the configured position with zero velocity, and FrameFixed() is true;
//   moving-cartesian:  x_hci(t) = x0 + v t, a straight line with constant
//       velocity in the inertial HCI frame (x0 and v are the observer
//       block's Cartesian position and velocity, t is model time);
//   body:  an ephemeris in the model HCI frame from either
//       - a checksummed CSV table (columns time_s, x_m, y_m, z_m,
//         vx_m_per_s, vy_m_per_s, vz_m_per_s), interpolated by cubic Hermite
//         polynomials with no extrapolation; or
//       - SPICE (output/spice_ephemeris.h), when the build links CSPICE.
//
// FRAMES AND UNITS:
//   SI throughout: m, m/s, and model time t [s], the seconds from the run's
//   time origin.  An ephemeris table is on the same time axis; the
//   2013-04-11 tables state t = 0 = 2013-04-11T07:26:44Z in their headers,
//   and they may begin before 0.  The model HCI frame is the inertial frame
//   of the background provider.  For the 2013-04-11 event, +z is the IAU
//   solar pole and +x the projection of the Sun-Earth vector at t0
//   (test/event_2013_04_11/README.md, section 3).  ObserverFramePosition
//   rotates inertial trajectories (moving-cartesian, body) into the
//   computation frame with Core::ComputationFrame::FromInertialPosition:
//   x_frame = Rz(-theta(t)) x_hci, theta(t) = Omega (t - t_ref).
//   Frame-fixed trajectories are returned unchanged.
//
// The module is independent of AMPS, PIC and MPI.
// Documentation: srcSEP3D/CONFIGURATION.md (observer.ID.kind = body, the
// ephemeris_* and spice_* keys) and srcSEP3D/test/event_2013_04_11/README.md,
// section 3.
// ============================================================================

#ifndef SEP3D_OUTPUT_OBSERVER_TRAJECTORY_H
#define SEP3D_OUTPUT_OBSERVER_TRAJECTORY_H

#include "../core/computation_frame.h"
#include "../runtime/run_configuration.h"

#include <memory>
#include <string>
#include <vector>

namespace SEP3D {
namespace Output {

// The position function of one observer.  Implementations are immutable
// after construction and are shared through ObserverTrajectoryPointer.  The
// table kind is safe for concurrent queries; the SPICE kind serializes its
// CSPICE calls internally (spice_ephemeris.cpp).
class ObserverTrajectory {
 public:
  virtual ~ObserverTrajectory() = default;
  // State at model time t [s]: position [m] and velocity [m/s] in the model
  // HCI frame for an inertial trajectory (FrameFixed() == false).  For a
  // frame-fixed trajectory it is the fixed computation-frame position with
  // zero velocity; use ObserverFramePosition, which handles both cases.
  // Either output may be null.  Errors: a time outside the support (table:
  // NotFound), or a SPICE failure (InvalidInput; ReservedFeature in a build
  // without SPICE).
  virtual Core::Status StateHci(double epochS, Core::Vec3* positionM,
                                Core::Vec3* velocityMPerS) const = 0;
  // Closed support [begin, end] of the trajectory in model time [s].  The
  // unbounded kinds report +-1e300.  The SPICE kind also reports +-1e300,
  // because its coverage is the loaded kernels' and is checked by every
  // query.
  virtual double BeginS() const = 0;
  virtual double EndS() const = 0;
  // true: fixed in the computation frame; false: an inertial trajectory that
  // ObserverFramePosition rotates into the computation frame.
  virtual bool FrameFixed() const = 0;
  // Human-readable kind and source, for logs and diagnostics.
  virtual std::string Description() const = 0;
};

using ObserverTrajectoryPointer = std::shared_ptr<const ObserverTrajectory>;

// One ephemeris sample: model time [s], position [m] and velocity [m/s] in
// the model HCI frame.
struct EphemerisSample {
  double timeS = 0.0;
  Core::Vec3 positionM, velocityMPerS;
};

// Parse an ephemeris CSV.
//   text     the whole file contents (the caller has already checked the
//            SHA-256 of exactly these bytes);
//   origin   the file name, used only in diagnostics ("origin:line: ...");
//   samples  out: the rows in file order.
// Format: empty lines and lines that start with '#' are skipped.  A trailing
// '\r' (CRLF files) is removed.  The first other line must be exactly
// "time_s,x_m,y_m,z_m,vx_m_per_s,vy_m_per_s,vz_m_per_s".  Every following
// line has exactly 7 finite numbers in SI units, model HCI frame.  Times
// strictly increase, and at least two rows are required.
// Errors: InvalidInput naming the file and line.
Core::Status ParseEphemerisTable(const std::string& text,
                                 const std::string& origin,
                                 std::vector<EphemerisSample>* samples);

// Cubic Hermite interpolation between the two samples that bracket t, using
// the tabulated positions and velocities as end values and end slopes.  The
// position is C1 across samples, and the returned velocity is the exact
// derivative of the interpolated position.
//   samples   at least two rows with strictly increasing times;
//   epochS    model time t [s];
//   positionM, velocityMPerS  out (either may be null), model HCI frame.
// Errors: InvalidInput for too few samples or a nonfinite t; NotFound for t
// outside [first time, last time].  There is no extrapolation.
Core::Status InterpolateEphemeris(const std::vector<EphemerisSample>& samples,
                                  double epochS, Core::Vec3* positionM,
                                  Core::Vec3* velocityMPerS);

// Build the trajectory of one configured observer (the configuration has
// already validated the kind-specific keys; runtime/run_configuration.cpp).
//   moving-cartesian  -> straight HCI line x0 + v t;
//   body, table       -> the file at ephemerisResolvedPath, which is read
//                        whole.  Its SHA-256 must equal ephemerisSha256
//                        (lower-case hex), and then it is parsed;
//   body, spice       -> the meta-kernel is loaded and the reference UTC
//                        converted once (PrepareSpiceEphemeris);
//   any other kind    -> fixed at observer.positionM in the computation frame.
// Errors: InvalidInput for an unreadable file, a checksum mismatch, a parse
// error, or a body with no ephemeris source.  SPICE preparation failures
// are returned with the observer ID prefixed (ReservedFeature in a build
// without SPICE).
Core::Status CreateObserverTrajectory(
    const RuntimeModel::ObserverOptions& observer,
    ObserverTrajectoryPointer* output);

// Observer position [m] in the computation frame at model time t [s]:
// inertial trajectories are rotated with frame.FromInertialPosition(x, t);
// frame-fixed trajectories are returned as they are.  The result is
// Sun-centred.  The caller adds the mesh coordinate origin where it needs
// mesh coordinates.  Errors: those of StateHci.
Core::Status ObserverFramePosition(const ObserverTrajectory& trajectory,
                                   const Core::ComputationFrame& frame,
                                   double epochS, Core::Vec3* positionM);

}  // namespace Output
}  // namespace SEP3D

#endif  // SEP3D_OUTPUT_OBSERVER_TRAJECTORY_H

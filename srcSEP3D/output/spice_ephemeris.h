// ============================================================================
// srcSEP3D/output/spice_ephemeris.h
//
// SPICE ephemeris for body observers (observer.ID.kind = body,
// ephemeris = spice).
//
// WHY THIS MODULE EXISTS:
//   A body observer (a planet, a spacecraft) needs its heliocentric
//   trajectory both to build the observer flux tube and to place its
//   collection sphere (output/observer_trajectory.h).  The checksummed
//   ephemeris table works in every build.  This interface lets a deck ask
//   SPICE instead, for any body that the loaded kernels cover, so more
//   observers can be added without preparing tables.
//
// WHAT IT COMPUTES:
//   The geometric state (aberration correction "NONE") of `spice_target`
//   relative to `spice_center` in the SPICE frame `spice_frame`
//   (CSPICE spkezr_c), at the ephemeris time
//
//     et = str2et(spice_reference_utc) + t,   t = model time [s]
//
//   (et in TDB seconds past J2000, as str2et_c returns it).  The state is
//   rotated into the model HCI frame with the deck's
//   `spice_to_model_rotation` R (row-major 3x3, x_model = R x_spice;
//   validated as proper orthonormal by runtime/run_configuration.cpp).  The
//   event preparer writes R, because it defines the model frame.  SPICE
//   kilometres and km/s become metres and m/s.  The kernels are loaded from
//   the meta-kernel `spice_kernels` once per process and path.
//
// BUILD:
//   The implementation is compiled against SpiceUsr.h, with
//   -DSEP3D_WITH_SPICE, only when the AMPS configuration names a SPICE
//   toolkit (Config.pl -spice-path=DIR; srcSEP3D/makefile).  A `nospice`
//   build compiles a refusal that returns StatusCode::ReservedFeature.  A
//   deck that selects SPICE then stops at startup with that diagnosis
//   instead of silently using another ephemeris.
//
// CONCURRENCY: CSPICE is not thread safe.  Every call takes one process-wide
// mutex.  CSPICE errors are switched to RETURN mode and converted into a
// Status, so a SPICE failure never aborts the process.
//
// Documentation: srcSEP3D/CONFIGURATION.md (observer.ID.spice_* keys) and
// srcSEP3D/test/event_2013_04_11/README.md, section 3.
// ============================================================================

#ifndef SEP3D_OUTPUT_SPICE_EPHEMERIS_H
#define SEP3D_OUTPUT_SPICE_EPHEMERIS_H

#include "../core/sep3d_types.h"

#include <array>
#include <string>

namespace SEP3D {
namespace Output {

// One observer's SPICE query, copied from its observer.ID.spice_* keys.
struct SpiceEphemerisRequest {
  // target, center: SPICE body names or IDs (e.g. "EARTH", "SUN");
  // frame: the SPICE reference frame of the returned state;
  // kernels: path of the meta-kernel that furnsh_c loads;
  // referenceUtc: the UTC of model time 0, in any form that str2et_c accepts.
  std::string target, center, frame, kernels, referenceUtc;
  // R, row-major: x_model = R x_spice (the identity by default).
  std::array<double, 9> toModelRotation{{1, 0, 0, 0, 1, 0, 0, 0, 1}};
};

// true when this build links CSPICE.
bool SpiceEphemerisAvailable();

// Load the meta-kernel `request.kernels` (once per process and path; a
// failed load is not remembered, so a later call retries it) and convert
// `request.referenceUtc` to ephemeris time.
//   referenceEtS  out: str2et(referenceUtc) [s, TDB past J2000].
// Errors: InvalidInput carrying the CSPICE short and long messages (furnsh
// or str2et failure); ReservedFeature in a build without SPICE.
Core::Status PrepareSpiceEphemeris(const SpiceEphemerisRequest& request,
                                   double* referenceEtS);

// Geometric state of the target relative to the center, rotated into the
// model HCI frame, at model time t:
//   referenceEtS  from PrepareSpiceEphemeris;
//   epochS        model time t [s]; the query epoch is referenceEtS + t;
//   positionM [m], velocityMPerS [m/s]  out (either may be null).
// Errors: InvalidInput with the CSPICE messages (for example an epoch that
// the loaded SPK kernels do not cover); ReservedFeature without SPICE.
Core::Status SpiceEphemerisState(const SpiceEphemerisRequest& request,
                                 double referenceEtS, double epochS,
                                 Core::Vec3* positionM,
                                 Core::Vec3* velocityMPerS);

}  // namespace Output
}  // namespace SEP3D

#endif  // SEP3D_OUTPUT_SPICE_EPHEMERIS_H

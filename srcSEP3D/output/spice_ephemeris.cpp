// ============================================================================
// srcSEP3D/output/spice_ephemeris.cpp
//
// CSPICE-backed body ephemeris (SEP3D_WITH_SPICE) or its typed refusal (a
// build without SPICE).  Interface, frames and units: spice_ephemeris.h.
// srcSEP3D/makefile adds -DSEP3D_WITH_SPICE and the CSPICE include path for
// this file only, when Makefile.local names a SPICE toolkit.
// ============================================================================

#include "spice_ephemeris.h"

#include <mutex>
#include <set>

#if defined(SEP3D_WITH_SPICE)
#include "SpiceUsr.h"
#endif

namespace SEP3D {
namespace Output {

#if defined(SEP3D_WITH_SPICE)

namespace {

std::mutex gSpiceMutex;           // CSPICE is not thread safe: one lock for every call
std::set<std::string> gLoadedKernels;  // meta-kernel paths already furnished in this process

// After a CSPICE call: if CSPICE signalled an error (failed_c), return its
// short and long messages as InvalidInput, prefixed with `context`, and
// clear the error state (reset_c) so later calls work.  Otherwise OK.
// The long-message buffer holds CSPICE's maximum of 1840 characters plus
// the terminator.  The short one exceeds the 25-character maximum.
Core::Status SpiceFailure(const std::string& context) {
  if (!failed_c()) return Core::Status::OK();
  SpiceChar shortMessage[64] = {0}, longMessage[1841] = {0};
  getmsg_c("SHORT", sizeof(shortMessage), shortMessage);
  getmsg_c("LONG", sizeof(longMessage), longMessage);
  reset_c();
  return Core::Status(Core::StatusCode::InvalidInput,
                      context + ": " + shortMessage + " " + longMessage);
}

void ReturnOnError() {
  // Report CSPICE errors to the caller instead of aborting the process:
  // erract RETURN makes a failing routine return with failed_c() set, and
  // errprt NONE stops CSPICE from printing the message itself (SpiceFailure
  // puts it in the Status).  This runs once per process; callers hold
  // gSpiceMutex.
  static bool once = false;
  if (once) return;
  erract_c("SET", 0, const_cast<SpiceChar*>("RETURN"));
  errprt_c("SET", 0, const_cast<SpiceChar*>("NONE"));
  once = true;
}

}  // namespace

bool SpiceEphemerisAvailable() { return true; }

Core::Status PrepareSpiceEphemeris(const SpiceEphemerisRequest& request,
                                   double* referenceEtS) {
  if (!referenceEtS) return Core::Status(Core::StatusCode::InvalidInput,
                                         "SPICE reference output is null");
  std::lock_guard<std::mutex> lock(gSpiceMutex);
  ReturnOnError();
  // Furnish each meta-kernel path once per process.  Several observers that
  // share one meta-kernel load it once.  A failed load is removed from the
  // set again, so it is not mistaken for a loaded kernel.
  if (gLoadedKernels.insert(request.kernels).second) {
    furnsh_c(request.kernels.c_str());
    const Core::Status loaded = SpiceFailure("furnsh " + request.kernels);
    if (!loaded.ok()) {
      gLoadedKernels.erase(request.kernels);
      return loaded;
    }
  }
  // Ephemeris time of model time 0 (TDB seconds past J2000).  It needs the
  // leap-seconds kernel, which the meta-kernel must list.
  SpiceDouble et = 0.0;
  str2et_c(request.referenceUtc.c_str(), &et);
  const Core::Status converted =
      SpiceFailure("str2et " + request.referenceUtc);
  if (!converted.ok()) return converted;
  *referenceEtS = et;
  return Core::Status::OK();
}

Core::Status SpiceEphemerisState(const SpiceEphemerisRequest& request,
                                 double referenceEtS, double epochS,
                                 Core::Vec3* positionM,
                                 Core::Vec3* velocityMPerS) {
  std::lock_guard<std::mutex> lock(gSpiceMutex);
  ReturnOnError();
  // state = (x, y, z [km], vx, vy, vz [km/s]) of the target relative to the
  // center in `frame`, geometric ("NONE": no light-time or stellar
  // aberration correction).  The one-way light time is returned but not
  // used.
  SpiceDouble state[6] = {0}, lightTime = 0.0;
  spkezr_c(request.target.c_str(), referenceEtS + epochS,
           request.frame.c_str(), "NONE", request.center.c_str(), state,
           &lightTime);
  const Core::Status queried = SpiceFailure("spkezr " + request.target);
  if (!queried.ok()) return queried;
  // x_model = R x_spice with R row-major, then km -> m (and km/s -> m/s).
  // The same constant R rotates the velocity.  That is correct only when
  // spice_frame is an inertial frame fixed relative to the model HCI frame
  // (for example J2000 or ECLIPJ2000).  For a rotating spice_frame the
  // rotation would depend on time; this code does not apply that and does
  // not check the frame.
  const auto& r = request.toModelRotation;
  const auto rotate = [&](const double* v) {
    return Core::Vec3(r[0] * v[0] + r[1] * v[1] + r[2] * v[2],
                      r[3] * v[0] + r[4] * v[1] + r[5] * v[2],
                      r[6] * v[0] + r[7] * v[1] + r[8] * v[2]) * 1.0e3;
  };
  if (positionM) *positionM = rotate(state);
  if (velocityMPerS) *velocityMPerS = rotate(state + 3);
  return Core::Status::OK();
}

#else  // no CSPICE in this build

// Typed refusal: a deck that selects ephemeris = spice stops with
// ReservedFeature and a message naming the rebuild option.  No other
// ephemeris is substituted.

bool SpiceEphemerisAvailable() { return false; }

Core::Status PrepareSpiceEphemeris(const SpiceEphemerisRequest&, double*) {
  return Core::Status(Core::StatusCode::ReservedFeature,
                      "this build does not link SPICE (configured with "
                      "-spice-path=nospice); use ephemeris = table or rebuild "
                      "with -spice-path=<CSPICE directory>");
}

Core::Status SpiceEphemerisState(const SpiceEphemerisRequest&, double, double,
                                 Core::Vec3*, Core::Vec3*) {
  return Core::Status(Core::StatusCode::ReservedFeature,
                      "this build does not link SPICE");
}

#endif

}  // namespace Output
}  // namespace SEP3D

// Reconstruct Cassini sub-spacecraft tracks with the event-era SPICE kernels.
//
// Usage:
//   spice_ground_tracks LSK PCK SPK TAG UTC_REFERENCE T0 T1 DT OUTPUT.csv
//
// Positions are geometric (aberration correction NONE) in IAU_ENCELADUS, as
// in the legacy trajectory program used for the 2014 paper.  Both
// planetocentric and planetographic surface coordinates are written: the clean
// USGS mosaic is planetocentric, while the legacy paper program plotted the
// RECPGR planetographic latitude.  Lengths written to CSV are SI metres.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "SpiceUsr.h"

namespace {

double degrees(double radians) { return radians * dpr_c(); }

double west_from_east(double east_radians) {
  double west = std::fmod(-degrees(east_radians), 360.0);
  return west < 0.0 ? west + 360.0 : west;
}

double zero_to_360(double radians) {
  double angle = std::fmod(degrees(radians), 360.0);
  return angle < 0.0 ? angle + 360.0 : angle;
}

void fail_if_spice(const char* context) {
  if (!failed_c()) return;
  SpiceChar message[1841];
  getmsg_c("LONG", sizeof(message), message);
  std::fprintf(stderr, "SPICE ERROR during %s: %s\n", context, message);
  reset_c();
  std::exit(2);
}

double parse_number(const char* text, const char* name) {
  char* end = nullptr;
  const double value = std::strtod(text, &end);
  if (end == text || *end != '\0' || !std::isfinite(value)) {
    std::fprintf(stderr, "Invalid %s: %s\n", name, text);
    std::exit(2);
  }
  return value;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 10) {
    std::fprintf(stderr,
                 "Usage: %s LSK PCK SPK TAG UTC_REFERENCE T0 T1 DT OUTPUT.csv\n",
                 argv[0]);
    return 2;
  }

  const double t0 = parse_number(argv[6], "T0");
  const double t1 = parse_number(argv[7], "T1");
  const double dt = parse_number(argv[8], "DT");
  if (!(t1 >= t0) || !(dt > 0.0)) {
    std::fprintf(stderr, "Require T1 >= T0 and DT > 0\n");
    return 2;
  }

  SpiceChar return_action[] = "RETURN";
  SpiceChar null_device[] = "NULL";
  erract_c("SET", 0, return_action);
  errdev_c("SET", 0, null_device);
  furnsh_c(argv[1]);
  furnsh_c(argv[2]);
  furnsh_c(argv[3]);
  fail_if_spice("kernel loading");

  SpiceDouble et_reference;
  str2et_c(argv[5], &et_reference);
  fail_if_spice("UTC conversion");

  SpiceInt count;
  SpiceDouble radii[3];
  bodvrd_c("ENCELADUS", "RADII", 3, &count, radii);
  fail_if_spice("Enceladus radii lookup");
  const double equatorial_radius = radii[0];
  const double flattening = (radii[0] - radii[2]) / radii[0];

  FILE* output = std::fopen(argv[9], "wx");
  if (output == nullptr) {
    std::perror(argv[9]);
    return 2;
  }
  std::fprintf(output,
      "flyby,time_from_reference_s,utc,"
      "spacecraft_x_m,spacecraft_y_m,spacecraft_z_m,"
      "subpoint_x_m,subpoint_y_m,subpoint_z_m,altitude_m,"
      "west_lon_planetocentric_deg,lat_planetocentric_deg,"
      "west_lon_planetographic_deg,lat_planetographic_deg\n");

  for (double offset = t0; offset <= t1 + 0.25 * dt; offset += dt) {
    const SpiceDouble et = et_reference + offset;
    SpiceDouble state[6], light_time;
    spkezr_c("CASSINI", et, "IAU_ENCELADUS", "NONE", "ENCELADUS",
             state, &light_time);
    fail_if_spice("Cassini state");

    SpiceDouble subpoint[3], surface_epoch, surface_vector[3];
    subpnt_c("Intercept: ellipsoid", "ENCELADUS", et,
             "IAU_ENCELADUS", "NONE", "CASSINI", subpoint,
             &surface_epoch, surface_vector);
    fail_if_spice("sub-spacecraft point");

    SpiceDouble radius, east_centric, lat_centric;
    reclat_c(subpoint, &radius, &east_centric, &lat_centric);
    // For Enceladus, RECPGR follows the IAU positive-west planetographic
    // convention.  RECLAT always returns positive-east longitude.
    SpiceDouble west_graphic, lat_graphic, altitude;
    recpgr_c("ENCELADUS", subpoint, equatorial_radius, flattening,
             &west_graphic, &lat_graphic, &altitude);
    fail_if_spice("surface coordinate conversion");

    SpiceChar utc[80];
    et2utc_c(et, "ISOC", 3, sizeof(utc), utc);
    fail_if_spice("UTC formatting");
    std::fprintf(output,
        "%s,%.6f,%s,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,"
        "%.9f,%.9f,%.9f,%.9f\n",
        argv[4], offset, utc, 1000.0 * state[0], 1000.0 * state[1],
        1000.0 * state[2], 1000.0 * subpoint[0], 1000.0 * subpoint[1],
        1000.0 * subpoint[2], 1000.0 * vnorm_c(surface_vector),
        west_from_east(east_centric), degrees(lat_centric),
        zero_to_360(west_graphic), degrees(lat_graphic));
  }

  if (std::fclose(output) != 0) {
    std::perror(argv[9]);
    return 2;
  }
  kclear_c();
  return 0;
}

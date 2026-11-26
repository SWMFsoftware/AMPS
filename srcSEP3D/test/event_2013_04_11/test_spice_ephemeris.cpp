// ============================================================================
// srcSEP3D/test/event_2013_04_11/test_spice_ephemeris.cpp
//
// SPICE ephemeris check (E3D07): the SPICE body-observer path
// (ephemeris = spice, output/spice_ephemeris.cpp; README.md section 3)
// reproduces the checked-in Earth ephemeris table ephemeris/earth.csv, which
// prepare_event.py wrote from JPL Horizons DE441 vectors rotated into the
// model frame.
//
// WHY THIS TEST EXISTS:
//   The 2013-04-11 decks read every body observer from a table.  The SPICE
//   path is the interface for other bodies and events; this test checks that
//   it uses the same frame, epoch and units as the table path.
//
// WHAT IT DOES:
//   1. reads the Earth table with the srcSEP3D table parser;
//   2. reads the rotation from the J2000 ecliptic to the model frame,
//      "model_frame_rows_from_ecliptic_j2000" of deck3d_summary.json (written
//      by prepare_event.py; +z the IAU solar pole, +x the projection of the
//      Sun->Earth vector at model time 0);
//   3. asks CSPICE for the geometric (no aberration correction) state of
//      EARTH relative to the SUN in ECLIPJ2000 at
//      et = str2et("2013-04-11T07:26:44") + t, i.e. at model time t of the
//      event (model time 0 = 2013-04-11T07:26:44 UTC), rotated by those rows,
//      at every epoch t of the table;
//   4. passes when the largest position difference is below 10 km and the
//      largest velocity difference below 1e-2 m/s.
//
//   The two sides use different JPL ephemerides (the SPK given to SPICE,
//   e.g. DE430, and Horizons DE441).  They agree only to the difference
//   between the two ephemerides: README.md section 5 records 0.14 km and
//   3e-5 m/s over the 67 samples.  An error of frame, time scale or units
//   would be far larger than the limits.
//
// HOW IT RUNS:
//   run_spice_test.sh builds this file with -DSEP3D_WITH_SPICE against a
//   CSPICE toolkit and runs it with the arguments
//     <meta-kernel> <ephemeris/earth.csv> <deck3d_summary.json>
//   Exit status: 0 pass; 1 failed comparison or unreadable input; 2 usage.
// ============================================================================
#include "output/observer_trajectory.h"
#include "output/spice_ephemeris.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

using namespace SEP3D;

int main(int argc, char** argv) {
  if (argc != 4) {
    std::fprintf(stderr, "usage: %s META_KERNEL EARTH_CSV DECK3D_SUMMARY_JSON\n", argv[0]);
    return 2;
  }
  std::ifstream table(argv[2]);
  std::stringstream text;
  text << table.rdbuf();
  std::vector<Output::EphemerisSample> samples;
  Core::Status status = Output::ParseEphemerisTable(text.str(), argv[2], &samples);
  if (!status.ok()) { std::fprintf(stderr, "%s\n", status.message.c_str()); return 1; }
  // The ecliptic J2000 -> model rotation recorded by the event preparer.
  // The value is a list of three rows; the first nine numbers after the key's
  // opening bracket are its entries in row-major order.  This is not a
  // general JSON parser: it relies on the layout prepare_event.py writes.
  std::ifstream summary(argv[3]);
  std::stringstream json;
  json << summary.rdbuf();
  const std::string s = json.str();
  const auto at = s.find("\"model_frame_rows_from_ecliptic_j2000\"");
  std::vector<double> rows;
  const std::regex number("[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?");
  const auto open = at == std::string::npos ? at : s.find('[', at);
  for (auto it = std::sregex_iterator(s.begin() + static_cast<long>(open == std::string::npos ? s.size() : open), s.end(), number);
       it != std::sregex_iterator() && rows.size() < 9; ++it) {
    rows.push_back(std::stod(it->str()));
  }
  if (at == std::string::npos || rows.size() != 9) {
    std::fprintf(stderr, "no model-frame rotation in %s\n", argv[3]);
    return 1;
  }
  // The request an observer with ephemeris = spice would build from its deck
  // keys (spice_target, spice_center, spice_frame, the meta-kernel,
  // spice_reference_utc = model time 0, spice_to_model_rotation).
  Output::SpiceEphemerisRequest request;
  request.target = "EARTH";
  request.center = "SUN";
  request.frame = "ECLIPJ2000";
  request.kernels = argv[1];
  request.referenceUtc = "2013-04-11T07:26:44";
  for (int i = 0; i < 9; ++i) request.toModelRotation[i] = rows[i];
  // et0: ephemeris seconds (TDB past J2000) of model time 0.
  double et0 = 0.0;
  status = Output::PrepareSpiceEphemeris(request, &et0);
  if (!status.ok()) { std::fprintf(stderr, "%s\n", status.message.c_str()); return 1; }
  // Largest |dx| [m] and |dv| [m/s] over the table epochs (the table rows,
  // not interpolated values, so the comparison involves no interpolation).
  double worst = 0.0, worstSpeed = 0.0;
  for (const auto& sample : samples) {
    Core::Vec3 x, v;
    status = Output::SpiceEphemerisState(request, et0, sample.timeS, &x, &v);
    if (!status.ok()) { std::fprintf(stderr, "%s\n", status.message.c_str()); return 1; }
    worst = std::max(worst, (x - sample.positionM).Norm());
    worstSpeed = std::max(worstSpeed, (v - sample.velocityMPerS).Norm());
  }
  // DE430 (SPICE) and DE441 (Horizons) agree for Earth to well below 10 km;
  // a frame, time-scale or unit error would show as >> 1000 km.
  const bool pass = worst < 1.0e4 && worstSpeed < 1.0e-2;
  std::printf("%s E3D07 SPICE Earth state vs the Horizons table: max |dx| = %.3f km, "
              "max |dv| = %.3e m/s over %zu samples\n",
              pass ? "PASS" : "FAIL", worst / 1e3, worstSpeed, samples.size());
  return pass ? 0 : 1;
}

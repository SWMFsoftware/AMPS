// ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
// Run copy: $COPS_OUTPUT_ROOT/tests/nightfix_20261118/build/getq_probe.cpp (on Voyager: /data001/highresolu_vt/cops_output/tests/nightfix_20261118/build/getq_probe.cpp)
// SHA256 of the run copy = SHA256 of this file without these header lines: 7edd0efedfb04d77a5108f725d412f838cc60a9c4efe6f7faedc36bf86b069bc
// Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
// WHAT: C++ probe (2026-11-18). Prints cSourceModel::GetQ of the rebuilt LiouvilleCOPS objects; PASS if H2O
//   qmin = 0.02/pow(rh,p)*600 bitwise and CO2 qmin/qmax = 0.1.
// RUN: built and run by jobs/nightfix_build.sbatch against ../LiouvilleCOPS/build/*.o and native/build/ref/*.o.
// ==== end of header ====
// H2O night-side fix (2026-11-18): print cSourceModel::GetQ of the rebuilt LiouvilleCOPS objects (Source.cpp) for both
// heliocentric models and three distances; GetQ is private, so the class is opened for this probe only.
// Build: g++ -O2 -std=c++17 -fopenmp -I<CSPICE>/include -I<REF>/src getq_probe.cpp <REF>/build/*.o (no main.o) cspice.a
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#define private public
#include "Source.h"
#undef private
int main() {
  cSourceModel s;
  int bad = 0;
  for (const char* hm : {"shap5", "shap5_1"}) {
    s.HeliocentricModel = hm;
    for (double rh : {1.3, 3.8310817, 3.8334296}) {
      double qmin, qmax, cmin, cmax;
      s.GetQ(Const::H2O, rh, qmin, qmax);
      s.GetQ(Const::CO2, rh, cmin, cmax);
      double p = (std::string(hm) == "shap5") ? 4.2143229 : 4.8;
      double want_qmax = 1.0 / std::pow(rh, p) * 600, want_qmin = 0.02 / std::pow(rh, p) * 600;
      bool ok = qmin == want_qmin && qmax == want_qmax && std::fabs(qmin / qmax - 0.02) < 1e-15 &&
                std::fabs(cmin / cmax - 0.1) < 1e-15;
      bad += !ok;
      std::printf("%-8s rh_AU %.7f  H2O qmin %.17g qmax %.17g qmin/qmax %.17g | CO2 qmin/qmax %.17g  %s\n", hm, rh,
                  qmin, qmax, qmin / qmax, cmin / cmax, ok ? "OK" : "FAIL");
    }
  }
  std::printf("GetQ probe: %s\n", bad ? "FAILED" : "H2O qmin/qmax = 0.02 (bitwise 0.02/pow(rh,p)*600), CO2 0.1: PASSED");
  return bad ? 1 : 0;
}

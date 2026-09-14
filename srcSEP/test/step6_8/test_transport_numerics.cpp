#include "sep_transport_common.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

const double kSpeedOfLightMPerS = 299792458.0;
const double kProtonMassKg = 1.67262192369e-27;

bool Near(double actual, double expected, double tolerance) {
  return std::fabs(actual - expected) <= tolerance;
}

bool Report(const char* id, bool pass, const std::string& detail) {
  std::cout << (pass ? "PASS " : "FAIL ") << id << ": " << detail << '\n';
  return pass;
}

bool Core01() {
  const std::vector<double> speeds = {
      0.0, 1.0e3, 0.2 * kSpeedOfLightMPerS, 0.99 * kSpeedOfLightMPerS};
  bool pass = true;
  for (std::size_t i = 0; i < speeds.size(); ++i) {
    const SEP::Transport::ScalarResult momentum =
        SEP::Transport::MomentumFromSpeed(
            speeds[i], kProtonMassKg, kSpeedOfLightMPerS);
    const SEP::Transport::ScalarResult velocity =
        SEP::Transport::SpeedFromMomentum(
            momentum.value, kProtonMassKg, kSpeedOfLightMPerS);
    pass = pass && momentum.status.ok() && velocity.status.ok() &&
        Near(velocity.value, speeds[i],
             std::max(1.0e-9, 2.0e-14 * kSpeedOfLightMPerS));
  }
  pass = pass && !SEP::Transport::MomentumFromSpeed(
      kSpeedOfLightMPerS, kProtonMassKg, kSpeedOfLightMPerS).status.ok();
  return Report("CORE01", pass,
      "relativistic SI conversions round-trip and reject v>=c");
}

bool Core02() {
  const SEP::Transport::CoordinateAdvance advance =
      SEP::Transport::AdvanceCoordinate(
          2.0, 3.5, 7.0, 0.0, 10.0,
          SEP::Transport::BoundaryPolicy::Absorb);
  return Report("CORE02", advance.status.ok() &&
      Near(advance.positionM, 5.5, 1.0e-15) && !advance.crossedBoundary,
      "field-line arc-length coordinate advances without metric leakage");
}

bool Core03() {
  const SEP::Transport::CoordinateAdvance absorbed =
      SEP::Transport::AdvanceCoordinate(
          9.0, 3.0, 2.0, 0.0, 10.0,
          SEP::Transport::BoundaryPolicy::Absorb);
  const SEP::Transport::CoordinateAdvance reflected =
      SEP::Transport::AdvanceCoordinate(
          9.0, 3.0, 2.0, 0.0, 10.0,
          SEP::Transport::BoundaryPolicy::Reflect);
  return Report("CORE03",
      absorbed.status.code == SEP::Transport::StatusCode::OutOfDomain &&
      reflected.status.ok() && Near(reflected.positionM, 8.0, 1.0e-15) &&
      Near(reflected.directedSpeedMPerS, -2.0, 1.0e-15),
      "absorbing and reflecting boundary policies are explicit");
}

bool Core04() {
  const double gradient_per_m = -2.5e-10;
  const SEP::Transport::ScalarResult focusing =
      SEP::Transport::FocusingLengthFromDlnBds(gradient_per_m);
  const SEP::Transport::ScalarResult uniform =
      SEP::Transport::FocusingLengthFromDlnBds(0.0);
  return Report("CORE04", focusing.status.ok() &&
      Near(focusing.value, -1.0 / gradient_per_m, 1.0e-6) &&
      uniform.status.ok() && std::isinf(uniform.value),
      "focusing length uses the supplied d ln|B|/ds");
}

bool Core05() {
  const double initial_momentum = 4.0e-19;
  const double divergence_per_s = 2.0e-5;
  const double dt_s = 30.0;
  const SEP::Transport::ScalarResult result =
      SEP::Transport::ApplyAdiabaticMomentum(
          initial_momentum, divergence_per_s, dt_s);
  const double expected = initial_momentum *
      std::exp(-divergence_per_s * dt_s / 3.0);
  return Report("CORE05", result.status.ok() &&
      Near(result.value, expected, 1.0e-32),
      "adiabatic momentum uses the exact plasma-frame exponential");
}

bool Core06() {
  SEP::Transport::StepDiagnostics diagnostics;
  const std::vector<SEP::Transport::StepLimit> limits = {
      {"streaming", 0.4}, {"focusing", 0.2}, {"diffusion", 0.3}};
  const SEP::Transport::ScalarResult selected =
      SEP::Transport::SelectSubstep(1.0, limits, 1.0e-6, &diagnostics);
  const SEP::Transport::ScalarResult underflow =
      SEP::Transport::SelectSubstep(
          1.0, {{"diffusion", 1.0e-9}}, 1.0e-6, &diagnostics);
  return Report("CORE06", selected.status.ok() &&
      Near(selected.value, 0.2, 1.0e-15) && diagnostics.selections == 2 &&
      diagnostics.underflowRejects == 1 &&
      underflow.status.code == SEP::Transport::StatusCode::StepUnderflow,
      "composable limits record their cause and reject underflow");
}

bool Core07() {
  SEP::Transport::KeyedRandomStream first(17, 42, 6, 0);
  SEP::Transport::KeyedRandomStream repeat(17, 42, 6, 0);
  SEP::Transport::KeyedRandomStream other(17, 43, 6, 0);
  bool same = true;
  bool differs = false;
  for (int i = 0; i < 16; ++i) {
    const double a = first.UniformOpen01();
    const double b = repeat.UniformOpen01();
    const double c = other.UniformOpen01();
    same = same && a == b && a > 0.0 && a < 1.0;
    differs = differs || a != c;
  }
  return Report("CORE07", same && differs,
      "keyed streams reproduce independently of particle ordering");
}

}  // namespace

int main(int argc, char** argv) {
  // This file now owns only Step 6 common-kernel checks. Steps 7-9 use the
  // registered callbacks in sep_mover_validation.cpp, eliminating the former
  // duplicate Make-only mover implementations.
  if (argc != 2 || std::string(argv[1]) != "step6") {
    std::cerr << "usage: test_transport_numerics step6\n";
    return 2;
  }
  bool pass = true;
  pass = Core01() && pass;
  pass = Core02() && pass;
  pass = Core03() && pass;
  pass = Core04() && pass;
  pass = Core05() && pass;
  pass = Core06() && pass;
  pass = Core07() && pass;
  return pass ? EXIT_SUCCESS : EXIT_FAILURE;
}

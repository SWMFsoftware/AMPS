// ============================================================================
// Phase R1 bounded SWCME integration executable
//
// The supplied historical SWCME test manifest references several translation
// units that were not present in the source archive.  R1 must nevertheless
// provide executable evidence from the relocated canonical directory.  This
// small registry exercises the production 1-D header API and the compiled 3-D
// object from swcme.a.  It is not a substitute for the extended validation
// campaign; it is the fast common-library gate used by srcSEP and srcSEP3D.
// ============================================================================

#include "swcme1d.hpp"
#include "swcme3d.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>

namespace {

struct TestCase {
  const char* id;
  const char* group;
  const char* name;
  bool (*run)(std::string*);
};

bool FinitePositive(double value) {
  return std::isfinite(value) && value > 0.0;
}

bool NearlyEqual(double a, double b, double relative_tolerance) {
  const double scale = std::max(1.0, std::max(std::abs(a), std::abs(b)));
  return std::abs(a - b) <= relative_tolerance * scale;
}

bool TestConstants(std::string* message) {
  const bool ok = swcme3d::AU == swcme1d::AU &&
                  swcme3d::Rs == swcme1d::Rs &&
                  FinitePositive(swcme3d::AU) &&
                  FinitePositive(swcme3d::Rs);
  *message = ok ? "1-D headers and the compiled 3-D library expose identical SI scales"
                : "1-D/3-D AU or solar-radius constants differ";
  return ok;
}

bool TestOneDimensionalState(std::string* message) {
  try {
    const swcme1d::Model model;
    const swcme1d::StepState state = model.prepare_step(0.0);
    const double radius = swcme1d::AU;
    double density = 0.0;
    double velocity = 0.0;
    double br = 0.0;
    double bphi = 0.0;
    double bmag = 0.0;
    double div_v = 0.0;
    const swcme::ModelStatus status =
        model.evaluate_radii_with_B_div_checked(
            state, &radius, &density, &velocity, &br, &bphi, &bmag,
            &div_v, 1);
    const bool ok = status.ok() && FinitePositive(density) &&
                    FinitePositive(velocity) && FinitePositive(bmag) &&
                    std::isfinite(div_v);
    *message = ok ? "canonical 1-D prepared state evaluates a finite 1-AU ambient record"
                  : "canonical 1-D prepared-state evaluation failed";
    return ok;
  } catch (const std::exception& error) {
    *message = std::string("1-D preparation threw: ") + error.what();
    return false;
  }
}

bool TestThreeDimensionalState(std::string* message) {
  try {
    const swcme3d::Params parameters;
    const swcme3d::Model model(parameters);
    const swcme3d::StepState state = model.prepare_step(0.0);
    const double x = swcme3d::AU;
    const double y = 0.0;
    const double z = 0.0;
    double density = 0.0;
    double vx = 0.0, vy = 0.0, vz = 0.0;
    double bx = 0.0, by = 0.0, bz = 0.0;
    double div_v = 0.0;
    const swcme::ModelStatus status =
        model.evaluate_cartesian_with_B_div_checked(
            state, &x, &y, &z, &density, &vx, &vy, &vz,
            &bx, &by, &bz, &div_v, 1);
    const double speed = std::sqrt(vx * vx + vy * vy + vz * vz);
    const double bmag = std::sqrt(bx * bx + by * by + bz * bz);
    const bool ok = status.ok() && FinitePositive(density) &&
                    FinitePositive(speed) && FinitePositive(bmag) &&
                    std::isfinite(div_v);
    *message = ok ? "canonical compiled 3-D state evaluates a finite 1-AU ambient record"
                  : "canonical compiled 3-D prepared-state evaluation failed";
    return ok;
  } catch (const std::exception& error) {
    *message = std::string("3-D preparation threw: ") + error.what();
    return false;
  }
}

bool TestDimensionalAgreement(std::string* message) {
  try {
    const swcme1d::Model model_1d;
    const swcme1d::StepState state_1d = model_1d.prepare_step(0.0);
    const double radius = swcme1d::AU;
    double density_1d = 0.0, velocity_1d = 0.0;
    double br_1d = 0.0, bphi_1d = 0.0, bmag_1d = 0.0, div_1d = 0.0;
    const swcme::ModelStatus status_1d =
        model_1d.evaluate_radii_with_B_div_checked(
            state_1d, &radius, &density_1d, &velocity_1d, &br_1d,
            &bphi_1d, &bmag_1d, &div_1d, 1);

    const swcme3d::Params parameters_3d;
    const swcme3d::Model model_3d(parameters_3d);
    const swcme3d::StepState state_3d = model_3d.prepare_step(0.0);
    const double x = swcme3d::AU, y = 0.0, z = 0.0;
    double density_3d = 0.0;
    double vx = 0.0, vy = 0.0, vz = 0.0;
    double bx = 0.0, by = 0.0, bz = 0.0, div_3d = 0.0;
    const swcme::ModelStatus status_3d =
        model_3d.evaluate_cartesian_with_B_div_checked(
            state_3d, &x, &y, &z, &density_3d, &vx, &vy, &vz,
            &bx, &by, &bz, &div_3d, 1);

    // At +X in the equatorial plane, the local spherical basis is
    // e_r=+X and e_phi=+Y.  The comparison therefore needs no coordinate
    // rotation and directly detects a divergence between the shared common
    // physics used by the two public interfaces.
    const double tolerance = 2.0e-13;
    const bool ok = status_1d.ok() && status_3d.ok() &&
                    NearlyEqual(density_1d, density_3d, tolerance) &&
                    NearlyEqual(velocity_1d, vx, tolerance) &&
                    NearlyEqual(vy, 0.0, tolerance) &&
                    NearlyEqual(vz, 0.0, tolerance) &&
                    NearlyEqual(br_1d, bx, tolerance) &&
                    NearlyEqual(bphi_1d, by, tolerance) &&
                    NearlyEqual(bz, 0.0, tolerance) &&
                    NearlyEqual(div_1d, div_3d, tolerance);
    *message = ok ? "1-D and 3-D public APIs agree for the common equatorial ambient state"
                  : "1-D and 3-D common ambient records differ";
    return ok;
  } catch (const std::exception& error) {
    *message = std::string("cross-dimensional comparison threw: ") + error.what();
    return false;
  }
}

const TestCase kTests[] = {
    {"R1CFG01", "COMMON", "Canonical constant identity", TestConstants},
    {"R1D01", "1D", "One-dimensional prepared-state smoke", TestOneDimensionalState},
    {"R3D01", "3D", "Three-dimensional prepared-state smoke", TestThreeDimensionalState},
    {"R13D01", "1D<->3D", "Common ambient-state identity", TestDimensionalAgreement},
};

const TestCase* Find(const char* id) {
  for (const TestCase& test : kTests) {
    if (std::strcmp(test.id, id) == 0) return &test;
  }
  return nullptr;
}

int Run(const TestCase& test) {
  std::string message;
  const bool passed = test.run(&message);
  std::cout << '[' << test.id << "] " << (passed ? "PASS" : "FAIL")
            << " - " << message << '\n';
  return passed ? 0 : 1;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc == 2 && std::strcmp(argv[1], "--list") == 0) {
    for (const TestCase& test : kTests) {
      std::cout << std::left << std::setw(10) << test.id
                << std::setw(12) << test.group << test.name << '\n';
    }
    return 0;
  }

  if (argc == 3 && std::strcmp(argv[1], "--test") == 0) {
    const TestCase* test = Find(argv[2]);
    if (test == nullptr) {
      std::cerr << "ERROR: unknown test ID " << argv[2] << '\n';
      return 2;
    }
    return Run(*test);
  }

  if ((argc == 1) ||
      (argc == 2 && std::strcmp(argv[1], "--all") == 0)) {
    int failures = 0;
    for (const TestCase& test : kTests) failures += Run(test);
    return failures == 0 ? 0 : 1;
  }

  std::cerr << "usage: " << argv[0]
            << " [--list | --all | --test TEST_ID]\n";
  return 2;
}

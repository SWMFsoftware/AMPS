#include "test_framework.h"

#include "sep_coronal_cme/constants.h"
#include "sep_coronal_cme/mhd_jump_solver.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace SCCMTest { namespace {
using namespace SEP::CoronalCME;

bool Close(double a, double b, double relative = 1.0e-7,
           double absolute = 1.0e-10) {
  return std::abs(a - b) <= absolute +
      relative * std::max(std::abs(a), std::abs(b));
}

MhdPrimitiveState Upstream(double bnAlfven = 0.0,
                           double btAlfven = 0.0,
                           double gamma = 5.0 / 3.0) {
  const double magneticUnit = std::sqrt(
      Constants::kVacuumPermeabilityHPerM);
  // rho=1 and p=1/gamma give c_s=1 in this manufactured SI state.  Magnetic
  // components are specified directly in Alfvén-speed units.
  return {1.0, 1.0 / gamma, {0.0, 0.2, -0.1},
          {bnAlfven * magneticUnit, btAlfven * magneticUnit, 0.0}};
}

CriticalMachTable Table(MachConvention convention = MachConvention::Fast) {
  return {"ek-manufactured-v1", "0123456789abcdef", 5.0 / 3.0,
          convention, {0.0, 1.0},
          {0.0, Constants::kPi / 4.0, Constants::kPi / 2.0},
          {1.5, 2.0, 2.76, 1.1, 1.5, 2.0}};
}

void RH3D01() {
  const double gamma = 5.0 / 3.0, mach = 3.0;
  auto solution = SolveObliqueFastShock(
      Upstream(0.0, 0.0, gamma), {1.0, 0.0, 0.0}, mach, gamma);
  const double expected = (gamma + 1.0) * mach * mach /
      ((gamma - 1.0) * mach * mach + 2.0);
  Require(solution.ok() && Close(solution.value.compressionRatio, expected) &&
          solution.value.residuals.maximum < 1.0e-9,
          "hydrodynamic normal-shock limit is incorrect");
}

void RH3D02() {
  auto solution = SolveObliqueFastShock(
      Upstream(0.5, 0.0), {1.0, 0.0, 0.0}, 3.0, 5.0 / 3.0);
  Require(solution.ok() && Close(solution.value.compressionRatio, 3.0) &&
          solution.value.branch == "regular-parallel-fast",
          "parallel MHD branch did not recover the gas-dynamic compression");
}

void RH3D03() {
  auto upstream = Upstream(0.0, 0.5);
  auto solution = SolveObliqueFastShock(
      upstream, {1.0, 0.0, 0.0}, 3.0, 5.0 / 3.0);
  Require(solution.ok(), solution.status.message);
  Require(Close(solution.value.downstream.magneticFieldT.y /
                upstream.magneticFieldT.y,
                solution.value.compressionRatio) &&
          solution.value.residuals.maximum < 1.0e-9,
          "perpendicular MHD flux-freezing/jump residual failed");
}

void RH3D04() {
  auto solution = SolveObliqueFastShock(
      Upstream(0.4, 0.7), {1.0, 0.0, 0.0}, 3.5, 5.0 / 3.0);
  Require(solution.ok() && solution.value.compressionRatio > 1.0 &&
          solution.value.entropyLogIncrement > 0.0 &&
          solution.value.residuals.maximum < 1.0e-9,
          "general oblique reference state failed the full jump system");
}

void RH3D05() {
  auto upstream = Upstream(0.0, 0.2);
  auto characteristics = EvaluateMhdCharacteristics(
      upstream, {1.0, 0.0, 0.0}, 5.0 / 3.0);
  auto solution = SolveObliqueFastShock(upstream, {1.0, 0.0, 0.0},
      1.01 * characteristics.value.fastSpeedMPerS, 5.0 / 3.0, 1.0e-8);
  // The closer reference guards against confusing a small energy residual
  // with convergence to the physical root.  The compressive solution at
  // M_f=1.001 is distinct from the identity state but has only an O(1e-9)
  // entropy increase, so premature residual-based termination used to fail
  // the independent entropy/characteristic admission checks.
  auto closer = SolveObliqueFastShock(upstream, {1.0, 0.0, 0.0},
      1.001 * characteristics.value.fastSpeedMPerS, 5.0 / 3.0, 1.0e-8);
  Require(solution.ok() && solution.value.compressionRatio > 1.0 &&
          solution.value.compressionRatio < 1.1 && closer.ok() &&
          closer.value.compressionRatio > 1.0 &&
          closer.value.compressionRatio < solution.value.compressionRatio &&
          closer.value.entropyLogIncrement > 0.0 &&
          closer.value.residuals.maximum < 1.0e-8,
          "weak fast shocks did not converge to the physical unit-compression limit");
}

void RH3D06() {
  auto solution = SolveObliqueFastShock(
      Upstream(), {1.0, 0.0, 0.0}, 100.0, 5.0 / 3.0);
  Require(solution.ok() && solution.value.compressionRatio < 4.0 &&
          solution.value.compressionRatio > 3.99,
          "strong hydrodynamic compression did not approach four from below");
}

void RH3D07() {
  auto invalidPressure = Upstream();
  invalidPressure.pressurePa = -1.0;
  Require(!SolveObliqueFastShock(invalidPressure, {1.0, 0.0, 0.0},
                                 3.0, 5.0 / 3.0).ok(),
          "negative upstream pressure was accepted");
  Require(!SolveObliqueFastShock(Upstream(), {1.0, 0.0, 0.0},
                                 0.9, 5.0 / 3.0).ok(),
          "sub-fast candidate was accepted");
  Require(!SolveObliqueFastShock(Upstream(), {1.0, 0.0, 0.0},
                                 -1.0, 5.0 / 3.0).ok(),
          "expansion front was accepted");
}

void RH3D08() {
  auto coldPerpendicular = QueryCriticalMach(Table(), 0.0,
      Constants::kPi / 2.0, 5.0 / 3.0, MachConvention::Fast);
  auto middle = QueryCriticalMach(Table(), 0.5,
      Constants::kPi / 8.0, 5.0 / 3.0, MachConvention::Fast);
  Require(coldPerpendicular.ok() && middle.ok() &&
          Close(coldPerpendicular.value.value, 2.76) &&
          Close(middle.value.value, 1.525),
          "versioned critical-Mach interpolation/reference value failed");
  auto normalTable = Table(MachConvention::NormalAlfven);
  auto perpendicular = QueryCriticalMach(normalTable, 0.2,
      Constants::kPi / 2.0, 5.0 / 3.0,
      MachConvention::NormalAlfven, true);
  Require(perpendicular.ok() &&
          perpendicular.value.validity ==
              CriticalQueryValidity::NormalAlfvenInapplicable &&
          std::isfinite(perpendicular.value.value),
          "exact perpendicular normal-Alfven limit became NaN/false criticality");
}

void RH3D09() {
  auto mono = SolveObliqueFastShock(
      Upstream(0.0, 0.0, 1.4), {1.0, 0.0, 0.0}, 3.0, 1.4);
  auto fiveThirds = SolveObliqueFastShock(
      Upstream(), {1.0, 0.0, 0.0}, 3.0, 5.0 / 3.0);
  Require(mono.ok() && fiveThirds.ok() &&
          !Close(mono.value.compressionRatio,
                 fiveThirds.value.compressionRatio),
          "gamma_ad mutation did not change compression");
  Require(!QueryCriticalMach(Table(), 0.5, 0.2, 1.4,
                             MachConvention::Fast).ok(),
          "critical-Mach table accepted an EOS mismatch");
}

void RH3D10() {
  auto exact = SolveObliqueFastShock(
      Upstream(0.5, 0.0), {1.0, 0.0, 0.0}, 3.0, 5.0 / 3.0);
  double previousDifference = 1.0;
  for (double tangent : {1.0e-3, 1.0e-4, 1.0e-5}) {
    auto near = SolveObliqueFastShock(
        Upstream(0.5, tangent), {1.0, 0.0, 0.0}, 3.0, 5.0 / 3.0);
    Require(near.ok(), near.status.message);
    const double difference = std::abs(near.value.compressionRatio -
                                       exact.value.compressionRatio);
    Require(difference < previousDifference,
            "near-parallel continuation did not converge");
    previousDifference = difference;
  }
  Require(exact.ok() && previousDifference < 1.0e-8,
          "near-parallel scalar limit is not basis invariant");
}

void RH3D11() {
  auto missing = QueryCriticalMach(Table(), 2.0, 0.4,
      5.0 / 3.0, MachConvention::Fast);
  Require(missing.ok() && missing.value.validity ==
          CriticalQueryValidity::BetaOutsideCoverage,
          "beta-domain miss was not typed");
  auto diagnostic = ApplyCriticalityPolicy(
      2.0, missing.value, false, CriticalityMissPolicy::DiagnosticOnlyFast);
  auto excluded = ApplyCriticalityPolicy(
      2.0, missing.value, true,
      CriticalityMissPolicy::ExcludeSourceBudgeted);
  auto fatal = ApplyCriticalityPolicy(
      2.0, missing.value, true, CriticalityMissPolicy::FailPreflight);
  Require(diagnostic.ok() && diagnostic.value.sourceEligible &&
          excluded.ok() && excluded.value.excludedWithoutRenormalization &&
          !fatal.ok(), "criticality miss policy was not applied exactly");

  std::vector<CriticalCoverageSample> samples{
      {0.5, 0.2, 9.0, 90.0, 900.0, false},
      {2.0, 0.2, 1.0, 10.0, 100.0, false}};
  auto ledger = PreflightCriticalCoverage(
      Table(), samples, 5.0 / 3.0, MachConvention::Fast,
      0.11, 0.11, 0.11);
  Require(ledger.ok() && Close(ledger.value.areaFraction, 0.1) &&
          Close(ledger.value.numberFraction, 0.1) &&
          Close(ledger.value.energyFraction, 0.1),
          "criticality exclusion ledger did not close physical measures");
  Require(!PreflightCriticalCoverage(
      Table(), samples, 5.0 / 3.0, MachConvention::Fast,
      0.05, 0.11, 0.11).ok(),
      "over-budget critical-Mach coverage passed preflight");
}

}  // namespace

void RegisterStage6(Registry* tests) {
  (*tests)["RH3D01"] = RH3D01; (*tests)["RH3D02"] = RH3D02;
  (*tests)["RH3D03"] = RH3D03; (*tests)["RH3D04"] = RH3D04;
  (*tests)["RH3D05"] = RH3D05; (*tests)["RH3D06"] = RH3D06;
  (*tests)["RH3D07"] = RH3D07; (*tests)["RH3D08"] = RH3D08;
  (*tests)["RH3D09"] = RH3D09; (*tests)["RH3D10"] = RH3D10;
  (*tests)["RH3D11"] = RH3D11;
}

}  // namespace SCCMTest

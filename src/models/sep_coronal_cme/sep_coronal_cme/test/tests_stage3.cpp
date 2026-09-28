#include "test_framework.h"

#include "sep_coronal_cme/constants.h"
#include "sep_coronal_cme/interface_balance.h"
#include "sep_coronal_cme/source_surface_coupling.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace SCCMTest { namespace {
using namespace SEP::CoronalCME;

bool Close(double a, double b, double relative = 1.0e-9,
           double absolute = 1.0e-12) {
  return std::abs(a - b) <= absolute +
      relative * std::max(std::abs(a), std::abs(b));
}

FiniteShellScs Shell(double outer = 2.0) {
  // The positive monopole dominates the signed l=1 perturbation, making this
  // a legitimate manufactured unsigned boundary without cellwise clipping.
  auto shell = FiniteShellScs::Create(
      1.0, outer, {{0, 0, 4.0, 0.0}, {1, 0, 0.2, 0.0}});
  Require(shell.ok(), shell.status.message);
  return shell.value;
}

void SCS3D01() {
  auto shell = Shell();
  auto inner = shell.Evaluate(1.0, 1.1, 0.4);
  Require(inner.ok(), inner.status.message);
  const double y00 = 1.0 / std::sqrt(4.0 * Constants::kPi);
  const double y10 = std::sqrt(3.0 / (4.0 * Constants::kPi)) * std::cos(1.1);
  Require(Close(inner.value.brT, 4.0 * y00 + 0.2 * y10),
          "finite-shell inner normal boundary was not reproduced");
}

void SCS3D02() {
  auto outer = Shell().Evaluate(2.0, 0.9, 0.7);
  Require(outer.ok(), outer.status.message);
  Require(std::abs(outer.value.bThetaT) < 1.0e-12 &&
          std::abs(outer.value.bPhiT) < 1.0e-12,
          "finite-shell outer boundary is not radial");
}

void SCS3D03() {
  const Vec3 unsignedB{2.0, -1.0, 0.5};
  auto plus = RestoreSector(unsignedB, MagneticSector::Positive, false, false);
  auto minus = RestoreSector(unsignedB, MagneticSector::Negative, false, false);
  Require(plus.ok() && minus.ok() &&
          Close(Norm(plus.value.valueT), Norm(minus.value.valueT)) &&
          Close(Norm(plus.value.valueT + minus.value.valueT), 0.0),
          "discrete polarity restoration did not preserve unsigned field");
}

void SCS3D04() {
  auto balance = EvaluateMagneticInterface(
      {3.0, 1.0, 0.0}, {3.0, 2.0, 0.0}, {1.0, 0.0, 0.0}, 1.0e-14);
  Require(balance.ok(), balance.status.message);
  Require(Close(balance.value.normalJumpT, 0.0) &&
          Norm(balance.value.surfaceCurrentAPerM) > 0.0,
          "sharp PFSS/SCS jump lost its surface current");
}

void SCS3D05() {
  auto first = RestoreSector({1.0, 2.0, 3.0}, MagneticSector::Negative,
                             true, true, InterfaceSide::Minus);
  auto second = RestoreSector({1.0, 2.0, 3.0}, MagneticSector::Negative,
                              true, true, InterfaceSide::Minus);
  Require(first.ok() && second.ok() &&
          first.value.sector == second.value.sector &&
          first.value.side == second.value.side,
          "sector/side identity is decomposition dependent");
}

void SCS3D06() {
  Require(!FiniteShellScs::Create(1.0, 2.0, {{1, 0, 1.0, 0.0}}).ok(),
          "unsigned shell without monopole was accepted");
  Require(!FiniteShellScs::Create(1.0, 2.0,
      {{0, 0, -1.0, 0.0}, {1, 0, 0.1, 0.0}}).ok(),
      "nonpositive unsigned-flux mode was accepted");
  Require(!FiniteShellScs::Create(1.0, 2.0,
      {{0, 0, 1.0, 0.0}, {1, 0, 10.0, 0.0}}).ok(),
      "constrained SCS boundary with a negative lobe was accepted");
}

void SCS3D07() {
  auto thin = ScsHarmonicAttenuation(2, 2.5 / 2.3);
  auto thick = ScsHarmonicAttenuation(2, 2.0);
  Require(thin.ok() && thick.ok() && Close(thin.value, 0.9800, 6.0e-5) &&
          thick.value < thin.value,
          "SCS attenuation does not respond to shell thickness");
}

void SCS3D08() {
  auto a = FiniteShellScs::Create(1.0, 2.0, {{0, 0, 4.0, 0.0}});
  auto b = FiniteShellScs::Create(1.0, 3.0, {{0, 0, 4.0, 0.0}});
  Require(a.ok() && b.ok(), "monopole shell construction failed");
  auto fa = a.value.Evaluate(1.5, 1.0, 0.0);
  auto fb = b.value.Evaluate(1.5, 1.0, 0.0);
  Require(fa.ok() && fb.ok() && Close(fa.value.brT, fb.value.brT) &&
          Close(fa.value.bThetaT, 0.0),
          "monopole gauge changed the physical magnetic field");
}

void SCS3D09() {
  auto diagnostics = EvaluateScsSpectrum(
      {{0, 0, 4.0, 0.0}, {1, 0, 1.0, 0.0}, {2, 1, 0.5, 0.2}}, 1.5);
  Require(diagnostics.ok() &&
          diagnostics.value.outerNonMonopolePowerT2 <
              diagnostics.value.innerNonMonopolePowerT2 &&
          diagnostics.value.outerZonalFraction > 0.0,
          "SCS radialization power diagnostic is inconsistent");
  Require(!ScsHarmonicAttenuation(2, 1.0).ok(),
          "an inert zero-thickness SCS shell was accepted");
}

void HCS3D01() { SCS3D03(); }

void HCS3D02() {
  Require(!RestoreSector({1.0, 0.0, 0.0}, MagneticSector::Positive,
                         true, false).ok(),
          "an unsided ideal-HCS interpolation was accepted");
}

void HCS3D03() {
  Require(!ValidateIdealHcsTransport(
      IdealHcsTransport::CrossSectorOrDrift, false).ok(),
      "cross-sector transport silently used an ideal sign reversal");
  Require(ValidateIdealHcsTransport(
      IdealHcsTransport::FieldAlignedNoCrossing, false).ok(),
      "field-aligned no-crossing ideal-HCS branch was rejected");
}

void CPL3D01() { SCS3D04(); }

void CPL3D02() {
  auto map = IntegrateLongitudeMap(0.2, 1.0, 5.0, 400,
      [](double, double) { return 0.1; },
      [](double, double) { return 0.0; }, 0.1);
  Require(map.ok(), map.status.message);
  auto state = MapParkerState(map.value, 1.0, 5.0,
      Constants::kPi / 2.0, 2.0, 3.0, 4.0, 4.0, 0.1);
  Require(state.ok(), state.status.message);
  Require(Close(state.value.magneticSphericalT.z /
                state.value.magneticSphericalT.x, -0.5),
          "axisymmetric field did not approach the Parker winding angle");
}

void CPL3D03() {
  const double k1 = 0.08;
  auto map = IntegrateLongitudeMap(0.3, 1.0, 3.0, 1000,
      [=](double, double phi) { return k1 * std::sin(phi); },
      [=](double, double phi) { return k1 * std::cos(phi); }, 0.1);
  Require(map.ok(), map.status.message);
  auto state = MapParkerState(map.value, 1.0, 3.0, 1.1,
                              2.0, 4.0, 5.0, 7.0,
                              k1 * std::sin(map.value.mappedLongitudeRad));
  Require(state.ok(), state.status.message);
  const double expected = 1.0 * 1.0 * 4.0 * 5.0 *
      map.value.inverseJacobian;
  Require(Close(state.value.mappedMassFluxKgPerSPerSr, expected, 1.0e-8),
          "non-axisymmetric mass map omitted J_phi");
  Require(!Close(expected, 20.0, 1.0e-5),
          "manufactured map failed to distinguish a missing Jacobian");
}

void CPL3D04() {
  const double q = 0.04;
  auto forward = IntegrateLongitudeMap(0.7, 1.0, 4.0, 1500,
      [=](double, double phi) { return q * std::sin(phi); },
      [=](double, double phi) { return q * std::cos(phi); }, 0.1);
  Require(forward.ok(), forward.status.message);
  auto inverse = IntegrateLongitudeMap(forward.value.mappedLongitudeRad,
      4.0, 1.0, 1500,
      [=](double, double phi) { return q * std::sin(phi); },
      [=](double, double phi) { return q * std::cos(phi); }, 0.1);
  Require(inverse.ok() &&
          std::abs(std::remainder(inverse.value.mappedLongitudeRad - 0.7,
                                  2.0 * Constants::kPi)) < 1.0e-10 &&
          Close(forward.value.forwardJacobian *
                inverse.value.forwardJacobian, 1.0, 1.0e-9),
          "forward/inverse longitude mapping did not close");
}

void CPL3D05() { CPL3D03(); }

void CPL3D06() {
  auto folded = IntegrateLongitudeMap(0.0, 1.0, 5.0, 100,
      [](double, double) { return 0.0; },
      [](double, double) { return 2.0; }, 0.2);
  Require(!folded.ok(), "A_phi fold/floor did not fail transactionally");
}

void CPL3D07() {
  auto inertial = IntegrateLongitudeMap(0.0, 1.0, 3.0, 100,
      [](double, double) { return 0.2; },
      [](double, double) { return 0.0; }, 0.1);
  auto corotating = IntegrateLongitudeMap(0.0, 1.0, 3.0, 100,
      [](double, double) { return 0.0; },
      [](double, double) { return 0.0; }, 0.1);
  Require(inertial.ok() && corotating.ok() &&
          !Close(inertial.value.mappedLongitudeRad,
                 corotating.value.mappedLongitudeRad),
          "rotation-frame winding conventions collapsed to one result");
}

void CPL3D08() {
  Require(ValidateCouplingMode(CouplingMode::NoScsVerification, false,
                               2.0, 2.0, 0.0).ok(),
          "valid no-SCS verification contract rejected");
  Require(!ValidateCouplingMode(CouplingMode::NoScsVerification, true,
                                2.0, 2.0, 0.0).ok(),
          "production selected no-SCS verification mode");
}

void CPL3D09() {
  const Vec3 pfss{1.0, 2.0, 0.0}, scs{1.0, 0.0, 3.0};
  auto left = BlendVectorPotentialFields(2.0, 2.0, 3.0,
      {1.0, 0.0, 0.0}, pfss, scs, {}, {});
  auto right = BlendVectorPotentialFields(3.0, 2.0, 3.0,
      {1.0, 0.0, 0.0}, pfss, scs, {}, {});
  auto middle = BlendVectorPotentialFields(2.5, 2.0, 3.0,
      {1.0, 0.0, 0.0}, pfss, scs, {}, {0.0, 1.0, 0.0});
  Require(left.ok() && right.ok() && middle.ok() &&
          Close(Norm(left.value.magneticFieldT - pfss), 0.0) &&
          Close(Norm(right.value.magneticFieldT - scs), 0.0) &&
          Norm(middle.value.magneticFieldT -
               (0.5 * pfss + 0.5 * scs)) > 0.0,
          "vector-potential transition lost endpoint or curl cross term");
}

void CPL3D11() {
  Require(ValidateSignedPotentialQualification({true, true, true}).ok(),
          "qualified signed Mie potentials rejected");
  Require(!ValidateSignedPotentialQualification({false, true, true}).ok() &&
          !ValidateSignedPotentialQualification({true, false, true}).ok() &&
          !ValidateSignedPotentialQualification({true, true, false}).ok(),
          "invalid signed-potential/gauge trace was accepted");
}

void CPL3D12() {
  auto d9 = EvaluateTransitionDiagnostics(
      {2.0, 1.0, 0.0}, {-2.0, 1.0, 0.0}, {0.0, 1.0, 0.0},
      3.0, 12.0, 0.1, 0.2, 1.0e-12);
  Require(d9.ok(), d9.status.message);
  Require(Close(d9.value.absoluteCrossingFluxWb, 3.0) &&
          d9.value.antipodalityApplicable &&
          d9.value.antipodalityDefectRad > 0.0,
          "D9 one-sided transition measures are incomplete");
}

void PLS3D01() {
  auto sheet = BuildPlasmaSheetNormalization(
      0.0, 4.0, 0.2, PlasmaSheetAuthority::BaseDensity, true);
  Require(sheet.ok() && Close(sheet.value.baseDensityMultiplier, 4.0) &&
          Close(sheet.value.basePressureMultiplier, 4.0),
          "fixed-temperature plasma-sheet EOS was not preserved");
}

void PLS3D02() {
  auto none = BuildPlasmaSheetNormalization(
      0.3, 1.0, 0.2, PlasmaSheetAuthority::BaseDensity, true);
  Require(none.ok() && Close(none.value.contrast, 1.0),
          "zero plasma-sheet enhancement changed the wind authority");
}

void PLS3D03() {
  auto outer = BuildPlasmaSheetNormalization(
      0.0, 3.0, 0.2, PlasmaSheetAuthority::MassPerMagneticFlux, true);
  Require(outer.ok() && Close(outer.value.massLoadingMultiplier, 3.0) &&
          Close(outer.value.baseDensityMultiplier, 1.0),
          "plasma sheet modified two normalization authorities");
  Require(!BuildPlasmaSheetNormalization(
      0.0, 0.9, 0.2, PlasmaSheetAuthority::BaseDensity, true).ok(),
      "plasma-sheet contrast below unity was accepted");
}

void PLS3D04() {
  auto closed = EvaluateVolumeMomentumResidual(
      {1.0, 2.0, 3.0}, {4.0, 5.0, 6.0},
      {2.0, 3.0, 4.0}, {3.0, 4.0, 5.0});
  auto omitted = EvaluateVolumeMomentumResidual(
      {1.0, 2.0, 3.0}, {4.0, 5.0, 6.0},
      {2.0, 3.0, 4.0}, {});
  Require(closed.ok() && omitted.ok() && Close(closed.value.normNPerM3, 0.0) &&
          omitted.value.normNPerM3 > 0.0,
          "smooth plasma-sheet volume force residual was hidden");
}

void OFX3D01() {
  auto lifecycle = BuildOpenFluxLifecycle(
      2.0, 5.0, 5.0 * (1.0 + 1.0e-9), 1.0e-8,
      "construction-a", "qualification-b");
  Require(lifecycle.ok() && Close(lifecycle.value.passAScale, 2.5),
          "two-pass open-flux lifecycle did not rebuild/qualify Pass B");
  Require(!BuildOpenFluxLifecycle(2.0, 5.0, 5.0, 1.0e-8,
      "same", "same").ok(), "construction data reused for qualification");
}

void LOS3D01() {
  auto footprint = EvaluateBudgetRatio(2.0, 10.0, 0.25);
  auto exceeded = EvaluateBudgetRatio(3.0, 10.0, 0.25);
  const auto point = EvaluatePointClearance(true);
  const auto validPoint = EvaluatePointClearance(false);
  Require(footprint.ok() && footprint.value.withinBound && exceeded.ok() &&
          !exceeded.value.withinBound &&
          point.validity == MeasureValidity::RejectedClearance &&
          validPoint.validity == MeasureValidity::InapplicablePointMeasure,
          "finite-footprint/point clearance accounting is not typed");
}

}  // namespace

void RegisterStage3(Registry* tests) {
  (*tests)["SCS3D01"] = SCS3D01; (*tests)["SCS3D02"] = SCS3D02;
  (*tests)["SCS3D03"] = SCS3D03; (*tests)["SCS3D04"] = SCS3D04;
  (*tests)["SCS3D05"] = SCS3D05; (*tests)["SCS3D06"] = SCS3D06;
  (*tests)["SCS3D07"] = SCS3D07; (*tests)["SCS3D08"] = SCS3D08;
  (*tests)["SCS3D09"] = SCS3D09;
  (*tests)["HCS3D01"] = HCS3D01; (*tests)["HCS3D02"] = HCS3D02;
  (*tests)["HCS3D03"] = HCS3D03;
  (*tests)["CPL3D01"] = CPL3D01; (*tests)["CPL3D02"] = CPL3D02;
  (*tests)["CPL3D03"] = CPL3D03; (*tests)["CPL3D04"] = CPL3D04;
  (*tests)["CPL3D05"] = CPL3D05; (*tests)["CPL3D06"] = CPL3D06;
  (*tests)["CPL3D07"] = CPL3D07; (*tests)["CPL3D08"] = CPL3D08;
  (*tests)["CPL3D09"] = CPL3D09; (*tests)["CPL3D11"] = CPL3D11;
  (*tests)["CPL3D12"] = CPL3D12;
  (*tests)["PLS3D01"] = PLS3D01; (*tests)["PLS3D02"] = PLS3D02;
  (*tests)["PLS3D03"] = PLS3D03; (*tests)["PLS3D04"] = PLS3D04;
  (*tests)["OFX3D01"] = OFX3D01; (*tests)["LOS3D01"] = LOS3D01;
}

}  // namespace SCCMTest

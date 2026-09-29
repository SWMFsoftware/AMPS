#include "test_framework.h"

#include "sep_coronal_cme/constants.h"
#include "sep_coronal_cme/particle_source.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace SCCMTest { namespace {
using namespace SEP::CoronalCME;

bool Close(double a, double b, double tolerance = 1.0e-8) {
  return std::abs(a - b) <= tolerance *
      std::max({1.0, std::abs(a), std::abs(b)});
}

void SRC3D01() {
  auto allocation = AllocateLargestRemainder(10, {1.0, 2.0, 3.0});
  Require(allocation.ok() && allocation.value ==
          std::vector<std::uint64_t>({2, 3, 5}),
          "patch-area/rate allocation did not use conservative remainder");
}

void SRC3D02() {
  auto spectrum = MomentumSpectrum::LocalCompressionDsa(1.0, 10.0, 4.0);
  Require(spectrum.ok() && Close(spectrum.value.NumberExponent(), -2.0) &&
          Close(spectrum.value.Cdf(10.0), 1.0) &&
          Close(spectrum.value.InverseCdf(0.5), 20.0 / 11.0),
          "DSA f(p) to dN/dp p^2-Jacobian conversion is incorrect");
}

void SRC3D03() {
  auto allocation = AllocateLargestRemainder(7, {2.0, 5.0});
  const double baseWeight = 10.0;
  const double represented = baseWeight *
      std::accumulate(allocation.value.begin(), allocation.value.end(), 0.0);
  Require(allocation.ok() && Close(represented, 70.0),
          "sample allocation/base-weight closure changed physical rate");
}

void SRC3D04() {
  const double momentum = KeyedUniform01(7, 1, 2, 3, 4, 5, 6);
  const double pitch = KeyedUniform01(7, 2, 2, 3, 4, 5, 6);
  const double repeated = KeyedUniform01(7, 1, 2, 3, 4, 5, 6);
  Require(momentum == repeated && momentum != pitch,
          "keyed streams are not reproducible/independent");
}

void SRC3D05() {
  const double h = 1.0e-5;
  const double leftDerivative =
      (RadialSourceEnvelope(1.0 + h, 1.0, 2.0, true) - 1.0) / h;
  const double rightDerivative =
      (0.0 - RadialSourceEnvelope(2.0 - h, 1.0, 2.0, true)) / h;
  Require(RadialSourceEnvelope(1.5, 1.0, 2.0, true) == 0.5 &&
          RadialSourceEnvelope(2.0, 1.0, 2.0, true) == 0.0 &&
          std::abs(leftDerivative) < 1.0e-8 &&
          std::abs(rightDerivative) < 1.0e-8,
          "smooth radial source taper lacks exact support/endpoint slope");
}

void SRC3D06() {
  Require(EligiblePhysicalSourceRate(
              12.0, SourcePatchTopology::ClosedDiagnoseOnly, false) == 0.0 &&
          EligiblePhysicalSourceRate(
              12.0, SourcePatchTopology::OpenEligible, false) == 12.0,
          "closed diagnose-only patch emitted physical source");
}

void SRC3D07() {
  std::vector<CompiledSpeciesIdentity> species{
      {"proton", 0, "H+", Constants::kProtonMassKg,
       Constants::kElementaryChargeC, 1, true},
      {"alpha", 1, "He++", Constants::kAlphaMassKg,
       2.0 * Constants::kElementaryChargeC, 4, true}};
  Require(ValidateSourceSpecies(species, 2).ok(),
          "all compiled charged species did not bind exactly once");
  species[1].chargeC = 0.0;
  Require(!ValidateSourceSpecies(species, 2).ok(),
          "neutral source species was accepted");
  Require(ValidateSourceBudget({5.0, 100.0, 0.1, 2.0, 100.0, 0.1}).ok() &&
          !ValidateSourceBudget({20.0, 100.0, 0.1, 2.0, 100.0, 0.1}).ok(),
          "number/energy source budgets were not independently enforced");
}

void SRC3D08() {
  const double mass = Constants::kProtonMassKg;
  const double c = Constants::kSpeedOfLightMPerS;
  FourMomentum rest{mass * c * c, {0.0, 0.0, 0.0}};
  auto boosted = BoostFourMomentum(rest, {1.0e6, 2.0e5, 0.0});
  Require(boosted.ok() && Close(FourMomentumInvariant(rest),
                                FourMomentumInvariant(boosted.value), 1.0e-12),
          "source-frame boost did not preserve the four-momentum invariant");
}

void SRC3D09() {
  Require(EligiblePhysicalSourceRate(
              20.0, SourcePatchTopology::TransitionClearance, false) == 0.0,
          "PFSS/SCS transition-clearance patch injected particles");
}

void SRC3D10() {
  auto spectrum = MomentumSpectrum::FixedPowerLaw(1.0, 4.0, -1.0);
  double integral = 0.0;
  const int cells = 10000;
  for (int i = 0; i < cells; ++i) {
    const double p0 = 1.0 + 3.0 * i / cells;
    const double p1 = 1.0 + 3.0 * (i + 1) / cells;
    integral += 0.5 * (spectrum.value.DensityPerMomentum(p0) +
                       spectrum.value.DensityPerMomentum(p1)) * (p1 - p0);
  }
  Require(spectrum.ok() && Close(integral, 1.0, 1.0e-7),
          "direct momentum quadrature did not recover normalized source rate");
}

void SRC3D11() {
  const std::vector<std::pair<double, double>> rate{
      {0.0, 0.0}, {1.0, 2.0}, {2.0, 0.0}};
  const double one = IntegratePiecewiseLinearRate(rate, 0.0, 2.0);
  const double split = IntegratePiecewiseLinearRate(rate, 0.0, 1.0) +
      IntegratePiecewiseLinearRate(rate, 1.0, 2.0);
  Require(Close(one, 2.0) && Close(one, split),
          "event-split time-varying source integration is step dependent");
}

void SRC3D12() {
  Require(ValidateSourcePopulationMeaning(
              SourcePopulationMeaning::NetFirstPassageAtReferenceSurface,
              true).ok() &&
          !ValidateSourcePopulationMeaning(
              SourcePopulationMeaning::GrossShockEmission, true).ok(),
          "production source meaning was silently reinterpreted");
}

CohortKey Key() { return {"proton", 17, 3, 9}; }

void SRC3D13() {
  InjectionCommitRegistry commits;
  Require(commits.CommitOnce(3, 9) && !commits.CommitOnce(3, 9),
          "same physical generation/tick was injected twice");
  CohortLedger ledger;
  Require(ledger.Add(Key(), ReleaseLedgerTerm::CommittedFirstPassageRelease,
                     {10.0, 20.0, 0.0, {}}).ok() &&
          ledger.Add(Key(), ReleaseLedgerTerm::DelayedFrontReturn,
                     {4.0, 8.0, 9.0, {}}).ok(),
          "return ledger rejected a valid immutable cohort");
  auto result = ReduceFiniteHorizonReturn(ledger, Key(), 2.0, 10.0, 5.0);
  Require(result.ok() && Close(result.value.noFrontReturnFraction, 0.6),
          "delayed front return did not debit committed release once");
}

void SRC3D14() {
  std::vector<ReferenceSurfacePoint> points{
      {1, {0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {1.0, 0.0, 0.0}},
      {2, {0.0, 1.0, 0.0}, {2.0, 1.0, 0.0}, {1.0, 0.0, 0.0}}};
  Require(ValidateReferenceSurface(points, 2.0, 1.0e-12).ok(),
          "fixed one-sided reference surface was rejected");
  points[1].referencePointM.x = -2.0;
  Require(!ValidateReferenceSurface(points, 2.0, 1.0e-12).ok(),
          "front-crossing reference support was accepted");
}

void SRC3D15() {
  auto distribution = BuildFocusedEscapeDistribution(
      2.0, 1, 1.0, 0.0, 1000, 1.0e-12);
  auto sample = SampleFocusedMu(distribution.value, 0.5);
  Require(distribution.ok() &&
          distribution.value.status == FocusedEscapeStatus::Admissible &&
          Close(distribution.value.normalizationMPerS, 0.5, 1.0e-6) &&
          sample.ok() && sample.value > 0.0,
          "signed focused positive-flux law/conditional CDF is incorrect");
}

void SRC3D16() {
  auto tangent = BuildFocusedEscapeDistribution(
      10.0, -1, 0.0, 2.0, 100, 1.0e-12);
  auto none = BuildFocusedEscapeDistribution(
      10.0, 1, 0.0, 0.0, 100, 1.0e-12);
  Require(tangent.ok() && none.ok() &&
          tangent.value.status == FocusedEscapeStatus::Admissible &&
          Close(tangent.value.normalizationMPerS, 2.0) &&
          none.value.status == FocusedEscapeStatus::NoFocusedEscape,
          "tangent-field advection/no-escape typing is incorrect");
}

void SRC3D17() {
  DiffusionTensor tensor{4.0, 1.0, 0.0, 2.0, 0.0, 1.0};
  auto direction = ParkerConormalDirection(tensor, {1.0, 0.0, 0.0}, 1.0);
  Require(direction.ok() && Close(direction.value.x, 1.0) &&
          Close(direction.value.y, 0.25) &&
          !ParkerConormalDirection(tensor, {0.0, 0.0, 1.0}, 2.0).ok(),
          "anisotropic conormal or n.kappa.n guard is incorrect");
}

void SRC3D18() {
  auto zero = ShockAdjacentEscapeProbability(2.0, 0.0, 10.0, 5.0);
  auto finite = ShockAdjacentEscapeProbability(2.0, 1.0, 10.0, 5.0);
  Require(zero.ok() && finite.ok() && zero.value == 0.0 &&
          finite.value > 0.0 && finite.value < 1.0,
          "absorbing verification escape probability has wrong support/limit");
}

void SRC3D19() {
  auto constant = IntegrateTransportDepth(
      {{0.0, 2.0, 5.0}, {10.0, 2.0, 5.0}});
  auto variable = IntegrateTransportDepth(
      {{0.0, 1.0, 2.0}, {1.0, 2.0, 2.0}, {2.0, 3.0, 2.0}});
  Require(constant.ok() && variable.ok() && Close(constant.value, 4.0) &&
          Close(variable.value, 2.0),
          "transport-depth diagnostic quadrature is incorrect");
  CohortLedger ledger;
  ledger.Add(Key(), ReleaseLedgerTerm::CommittedFirstPassageRelease,
             {5.0, 5.0, 0.0, {}});
  ledger.Add(Key(), ReleaseLedgerTerm::SurvivingUpstreamInventory,
             {3.0, 3.0, 0.0, {}});
  auto horizon = ReduceFiniteHorizonReturn(ledger, Key(), 9.0, 10.0, 2.0);
  Require(horizon.ok() && horizon.value.rightCensored &&
          Close(horizon.value.survivingInventoryNumber, 3.0),
          "finite-horizon inventory/right censoring was not retained");
}

void LOS3D02() {
  auto pass = EvaluateRepresentedLossCaps(
      100.0, 5.0, 1000.0, 40.0, 0.05, 0.04, true);
  auto fail = EvaluateRepresentedLossCaps(
      100.0, 1.0, 1000.0, 1.0, 0.0, 0.0, true);
  auto inactive = EvaluateRepresentedLossCaps(
      100.0, 0.0, 1000.0, 0.0, 0.1, 0.1, false);
  Require(pass.ok() && pass.value.withinCaps && fail.ok() &&
          !fail.value.withinCaps && !inactive.ok(),
          "represented-number/birth-energy loss caps are incorrect");
}

}  // namespace

void RegisterStage9(Registry* tests) {
  (*tests)["SRC3D01"] = SRC3D01; (*tests)["SRC3D02"] = SRC3D02;
  (*tests)["SRC3D03"] = SRC3D03; (*tests)["SRC3D04"] = SRC3D04;
  (*tests)["SRC3D05"] = SRC3D05; (*tests)["SRC3D06"] = SRC3D06;
  (*tests)["SRC3D07"] = SRC3D07; (*tests)["SRC3D08"] = SRC3D08;
  (*tests)["SRC3D09"] = SRC3D09; (*tests)["SRC3D10"] = SRC3D10;
  (*tests)["SRC3D11"] = SRC3D11; (*tests)["SRC3D12"] = SRC3D12;
  (*tests)["SRC3D13"] = SRC3D13; (*tests)["SRC3D14"] = SRC3D14;
  (*tests)["SRC3D15"] = SRC3D15; (*tests)["SRC3D16"] = SRC3D16;
  (*tests)["SRC3D17"] = SRC3D17; (*tests)["SRC3D18"] = SRC3D18;
  (*tests)["SRC3D19"] = SRC3D19; (*tests)["LOS3D02"] = LOS3D02;
}

}  // namespace SCCMTest

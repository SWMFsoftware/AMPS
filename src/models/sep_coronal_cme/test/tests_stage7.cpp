#include "test_framework.h"

#include "sep_coronal_cme/constants.h"
#include "sep_coronal_cme/shock_provider.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

namespace SCCMTest { namespace {
using namespace SEP::CoronalCME;

bool Close(double a, double b, double tolerance = 1.0e-9) {
  return std::abs(a - b) <= tolerance *
      std::max({1.0, std::abs(a), std::abs(b)});
}

MhdPrimitiveState State() {
  return {1.0, 0.6, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
}

PatchInterfaceEvidence Evidence(
    InterfacePolicy policy = InterfacePolicy::BoundedApproximation) {
  return {policy, 7, true, true, true, true, "manufactured-interface-v1"};
}

ShockPatchInput Patch(std::uint64_t id, double speed, double area = 1.0,
                      bool clearance = false) {
  ShockPatchInput patch;
  patch.stableId = id;
  patch.areaM2 = area;
  patch.outwardNormal = {1.0, 0.0, 0.0};
  patch.upstream = State();
  patch.shockNormalSpeedMPerS = speed;
  patch.incidentNumberRatePerS = 10.0 * area;
  patch.incidentKineticEnergyRateW = 100.0 * area;
  patch.intersectsTransitionClearance = clearance;
  patch.interface = Evidence();
  return patch;
}

ShockPreparationOptions Options() {
  ShockPreparationOptions options;
  options.initialGate = InitialShockGate::None;
  return options;
}

void SHK3D05() {
  TransactionalShockProvider provider;
  auto result = provider.Prepare(0.0, 7, {Patch(1, 3.0), Patch(2, 0.5)},
                                 Options());
  Require(result.ok() && result.value->patches[0].sourceActive &&
          !result.value->patches[1].sourceActive &&
          Close(result.value->measures.fastAreaM2, 1.0),
          "mixed fast/sub-fast surface classification failed");
}

void SHK3D06() {
  TransactionalShockProvider provider;
  auto first = provider.Prepare(0.0, 7, {Patch(1, 3.0)}, Options());
  auto bad = Patch(2, 3.0);
  bad.areaM2 = -1.0;
  auto failed = provider.Prepare(1.0, 7, {bad}, Options());
  Require(first.ok() && !failed.ok() && provider.CurrentGeneration() == 1 &&
          provider.PreparedSurface() == first.value,
          "failed update changed the published generation");
}

void SHK3D07() {
  TransactionalShockProvider provider;
  auto result = provider.Prepare(0.0, 7, {Patch(1, 0.5)}, Options());
  Require(result.ok() && result.value->measures.fastAreaM2 == 0.0 &&
          result.value->measures.sourceActiveAreaM2 == 0.0,
          "gate none rejected a valid zero-source initial dome");
}

void SHK3D08() {
  auto crossing = LocateFirstFastCrossing(
      0.0, 10.0, [](double time) { return 0.5 + 0.1 * time; }, 1.0e-10);
  Require(crossing.ok() && Close(crossing.value, 5.0, 1.0e-9),
          "M_f=1 crossing was not independently root located");
}

void SHK3D09() {
  auto subfast = Patch(1, 0.5);
  TransactionalShockProvider noneProvider, anyProvider, fractionProvider;
  auto none = noneProvider.Prepare(0.0, 7, {subfast}, Options());
  auto anyOptions = Options();
  anyOptions.initialGate = InitialShockGate::AnyFastPatch;
  auto any = anyProvider.Prepare(0.0, 7, {subfast}, anyOptions);
  auto fractionOptions = Options();
  fractionOptions.initialGate = InitialShockGate::MinimumFastAreaFraction;
  fractionOptions.minimumFastAreaFraction = 0.6;
  auto fraction = fractionProvider.Prepare(
      0.0, 7, {Patch(2, 3.0), Patch(3, 0.5)}, fractionOptions);
  Require(none.ok() && !any.ok() && !fraction.ok(),
          "initial shock gate selectors did not preserve their semantics");
}

void SHK3D10() {
  auto parent = Patch(42, 3.0, 10.0);
  auto split = SplitShockPatch(
      parent, {0.25, 0.75},
      {Evidence(InterfacePolicy::BoundedApproximation),
       Evidence(InterfacePolicy::StationaryTangentialDiscontinuity)});
  Require(split.ok() && split.value[0].parentId == 42 &&
          split.value[0].stableId != split.value[1].stableId &&
          Close(split.value[0].areaM2 + split.value[1].areaM2, 10.0) &&
          Close(split.value[0].incidentNumberRatePerS +
                split.value[1].incidentNumberRatePerS, 100.0),
          "one-sided patch split did not close geometry/source measures");
}

void SHK3D11() {
  auto terminated = Patch(1, 3.0);
  terminated.sourceTerminated = true;
  TransactionalShockProvider provider;
  auto result = provider.Prepare(20.0, 7, {terminated}, Options());
  Require(result.ok() && result.value->patches[0].fast &&
          result.value->patches[0].sourceTerminated &&
          !result.value->patches[0].sourceActive,
          "terminated source was reactivated by a later fast geometry");
}

void SHK3D12() {
  TransactionalShockProvider provider;
  auto first = provider.Prepare(0.0, 7, {Patch(1, 3.0)}, Options());
  std::shared_ptr<const ShockSurfaceSnapshot> retained = first.value;
  auto second = provider.Prepare(1.0, 7, {Patch(2, 3.0)}, Options());
  Require(second.ok() && retained->generation == 1 &&
          retained->patches[0].stableId == 1 &&
          provider.PreparedSurface()->generation == 2,
          "immutable prepared-surface handle did not retain ownership");
}

void SHK3D13() {
  const auto law = [](double time) { return 0.8 + 0.04 * time; };
  auto coarse = LocateFirstFastCrossing(0.0, 20.0, law, 1.0e-8);
  auto fine = LocateFirstFastCrossing(4.0, 6.0, law, 1.0e-10);
  Require(coarse.ok() && fine.ok() && Close(coarse.value, fine.value, 1.0e-7),
          "event location depends on the outer update cadence");
}

void SHK3D14() {
  TransactionalShockProvider provider;
  auto options = Options();
  options.maximumTransitionAreaFraction = 0.26;
  options.maximumTransitionNumberFraction = 0.26;
  options.maximumTransitionEnergyFraction = 0.26;
  auto result = provider.Prepare(
      0.0, 7, {Patch(1, 3.0, 3.0, false), Patch(2, 3.0, 1.0, true)},
      options);
  Require(result.ok() && Close(result.value->measures.counterfactualAreaM2, 4.0) &&
          Close(result.value->measures.excludedAreaM2, 1.0) &&
          Close(result.value->measures.sourceActiveAreaM2, 3.0),
          "transition mask did not preserve counterfactual/excluded ledgers");
  options.maximumTransitionAreaFraction = 0.2;
  auto rejected = provider.Prepare(
      1.0, 7, {Patch(1, 3.0, 3.0, false), Patch(2, 3.0, 1.0, true)},
      options);
  Require(!rejected.ok() && provider.CurrentGeneration() == 1,
          "over-budget history was not rejected transactionally");
}

void SNAP3D09() {
  TransactionalShockProvider provider;
  auto first = provider.Prepare(0.0, 7, {Patch(1, 3.0)}, Options());
  auto diagnostic = Patch(2, 3.0);
  diagnostic.interface = Evidence(InterfacePolicy::DiagnosticKinematic);
  auto production = Options();
  production.productionIntent = true;
  auto failed = provider.Prepare(1.0, 7, {diagnostic}, production);
  auto bounded = Patch(3, 3.0);
  auto passed = provider.Prepare(2.0, 7, {bounded}, production);
  Require(first.ok() && !failed.ok() && passed.ok() &&
          passed.value->generation == 2,
          "interface provenance did not guard transactional publication");
}

}  // namespace

void RegisterStage7(Registry* tests) {
  (*tests)["SHK3D05"] = SHK3D05; (*tests)["SHK3D06"] = SHK3D06;
  (*tests)["SHK3D07"] = SHK3D07; (*tests)["SHK3D08"] = SHK3D08;
  (*tests)["SHK3D09"] = SHK3D09; (*tests)["SHK3D10"] = SHK3D10;
  (*tests)["SHK3D11"] = SHK3D11; (*tests)["SHK3D12"] = SHK3D12;
  (*tests)["SHK3D13"] = SHK3D13; (*tests)["SHK3D14"] = SHK3D14;
  (*tests)["SNAP3D09"] = SNAP3D09;
}

}  // namespace SCCMTest

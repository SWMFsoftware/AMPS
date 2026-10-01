// Portable regression checks for the production native boundary evaluator and
// the shared corridor/solar geometry.  This is deliberately not an AMPS mock:
// it cannot establish native allocation, MPI exchange or particle stepping.
#include "../validation/coronal_cme_application_test.h"
#include "../mesh/mesh_model.h"
#include "../runtime/configuration_io.h"

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
namespace V = SEP3D::Validation;
namespace M = SEP3D::Mesh;
namespace R = SEP3D::RuntimeModel;
namespace C = SEP3D::Core;
unsigned checks = 0;
void Require(bool condition, const std::string& description) {
  ++checks;
  if (!condition) throw std::runtime_error(description);
}
V::NativeApplicationState FullDomain() {
  V::NativeApplicationState state;
  state.inputSchemaVersion = 4;
  state.backgroundAuthority = "analytic-parker";
  state.shockAuthority = "swcme";
  state.solarBoundaryRegistered = true;
  state.activeMaskInstalled = true;
  state.activeRegionAllocationVerified = true;
  state.activeRegionMode = "full-domain";
  state.plannedActiveLeaves = 100;
  state.globalAllocatedBlocks = 100;
  return state;
}
V::NativeTestResult Check(const V::NativeApplicationState& state) {
  std::vector<V::NativeTestDescriptor> selected;
  const auto status = V::SelectCoronalCmeNativeTests(false, {"SCCM3D04"}, &selected);
  if (!status.ok()) throw std::runtime_error(status.message);
  return V::EvaluateCoronalCmeNativeTests(state, selected).front();
}
void ExpectFailure(V::NativeApplicationState state, const std::string& reason) {
  const auto result = Check(state);
  Require(result.status == V::NativeTestStatus::Fail &&
              result.message.find(reason) != std::string::npos,
          "invalid native boundary state was accepted or misdiagnosed: " + reason);
}
// Exercise the production parser and selector without allocating AMPS state.
// Expected membership is read from descriptors so this check grows with the
// suite. It does not assert a fixed SCCM ID range or silently drop new cases.
void SuiteSelectionChecks() {
  std::vector<V::NativeTestDescriptor> selected;
  auto status = V::SelectCoronalCmeNativeTests(false, {}, &selected, "sep-corona");
  Require(status.ok() && !selected.empty(), "coupled suite selection failed");
  std::vector<std::string> expected, actual;
  for (const auto& test : V::CoronalCmeNativeTests())
    if (test.suite == "sep-corona") expected.push_back(test.id);
  for (const auto& test : selected) actual.push_back(test.id);
  Require(actual == expected, "coupled suite omitted or reordered registered tests");
  const auto results = V::EvaluateCoronalCmeNativeTests(FullDomain(), selected);
  Require(results.size() == expected.size(), "suite evaluator omitted a test");
  actual.clear();
  for (const auto& result : results) actual.push_back(result.id);
  Require(actual == expected, "suite evaluator lost test identity/order");
  status = V::SelectCoronalCmeNativeTests(true, {}, &selected);
  Require(status.ok() && selected.size() == V::CoronalCmeNativeTests().size(),
          "legacy all-tests selector no longer covers the whole registry");
  Require(!V::SelectCoronalCmeNativeTests(false, {}, &selected, "typo").ok(),
          "unknown suite was silently accepted");
  Require(!V::SelectCoronalCmeNativeTests(true, {}, &selected, "sep-corona").ok(),
          "conflicting all/suite selectors were accepted");
  Require(!V::SelectCoronalCmeNativeTests(false, {"SCCM3D04"}, &selected,
              "sep-corona").ok(), "conflicting suite/ID selectors were accepted");
  auto parse = [](std::vector<std::string> args, R::StandaloneCommandLine* cli) {
    std::vector<char*> argv;
    for (auto& arg : args) argv.push_back(&arg[0]);
    return R::ParseStandaloneCommandLine(static_cast<int>(argv.size()), argv.data(), cli);
  };
  R::StandaloneCommandLine cli;
  const std::vector<std::string> base = {"amps", "--test-suite", "sep-corona",
      "--test-input", "run.in", "--test-steps", "0", "--expect-mpi-ranks", "4"};
  Require(parse(base, &cli).ok() && cli.testSuite == "sep-corona" &&
              cli.tests.empty() && !cli.allTests && cli.testSteps == 0 &&
              cli.expectedMpiRanks == 4 && cli.testInputPath == "run.in",
          "suite CLI did not normalize the documented initialization command");
  auto mixedCase = base; mixedCase[2] = "SEP-CORONA";
  Require(parse(mixedCase, &cli).ok() && cli.testSuite == "sep-corona",
          "suite name did not normalize case");
  for (const auto& extra : std::vector<std::vector<std::string>>{
           {"--all-tests"}, {"--test", "SCCM3D04"}, {"--list-tests"},
           {"--dry-run"}, {"--initialization-only"},
           {"--test-suite", "sep-corona"}}) {
    auto args = base; args.insert(args.end(), extra.begin(), extra.end());
    Require(!parse(args, &cli).ok(), "suite CLI accepted a conflicting selector/mode");
  }
  Require(!parse({"amps", "--test-suite", "unknown", "--test-input", "run.in"}, &cli).ok(),
          "suite CLI accepted an unknown suite");
  Require(!parse({"amps", "--test-suite"}, &cli).ok(),
          "suite CLI accepted a missing name");
  Require(!parse({"amps", "--test-suite", "sep-corona"}, &cli).ok(),
          "suite CLI accepted a missing input deck");
  Require(parse({"amps", "--all-tests", "--test-input", "run.in"}, &cli).ok() && cli.allTests,
          "legacy all-tests CLI regressed");
  Require(parse({"amps", "--input", "run.in"}, &cli).ok() && cli.testSuite.empty(),
          "production CLI regressed");
}
void NativeEvidenceChecks() {
  auto state = FullDomain();
  Require(Check(state).status == V::NativeTestStatus::Pass,
          "full-domain identity plan with no pruning must pass");
  // Full-domain mode may still remove leaves wholly inside the physical Sun.
  state.plannedInactiveLeaves = state.plannedSolarInteriorLeaves = 8;
  state.activeRegionPruningApplied = true;
  Require(Check(state).status == V::NativeTestStatus::Pass,
          "full-domain solar-interior exclusion must pass");
  state.activeRegionMode = "parker-tube";
  state.plannedInactiveLeaves = 900;
  Require(Check(state).status == V::NativeTestStatus::Pass,
          "pruned corridor composed with solar exclusion must pass");
  // A corridor wider than the mesh is valid, but provides no pruning benefit.
  state = FullDomain(); state.activeRegionMode = "parker-tube";
  Require(Check(state).status == V::NativeTestStatus::Pass,
          "valid wide corridor must not require an artificial inactive leaf");
  state = FullDomain(); state.activeMaskInstalled = false;
  ExpectFailure(state, "plan was not installed");
  state = FullDomain(); state.solarBoundaryRegistered = false;
  ExpectFailure(state, "solar absorption boundary");
  state = FullDomain(); state.activeRegionAllocationVerified = false;
  ExpectFailure(state, "allocation was not verified");
  state = FullDomain(); state.activeRegionMode = "unknown";
  ExpectFailure(state, "mode is missing or unsupported");
  state = FullDomain(); state.plannedActiveLeaves = 0;
  ExpectFailure(state, "no physical active leaves");
  state = FullDomain(); state.globalAllocatedBlocks = 99;
  ExpectFailure(state, "allocated owner-block count");
  state = FullDomain(); state.activeRegionPruningApplied = true;
  ExpectFailure(state, "pruning flag disagrees");
  state = FullDomain(); state.plannedInactiveLeaves = 1;
  ExpectFailure(state, "pruning flag disagrees");
  state = FullDomain(); state.activeRegionMode = "parker-tube";
  state.plannedInactiveLeaves = 1; state.plannedSolarInteriorLeaves = 2;
  state.activeRegionPruningApplied = true;
  ExpectFailure(state, "solar-interior leaves are not all excluded");
  state = FullDomain(); state.plannedInactiveLeaves = 9;
  state.plannedSolarInteriorLeaves = 8; state.activeRegionPruningApplied = true;
  ExpectFailure(state, "removed leaves outside the solid Sun");
}
void CorridorAndSolarGeometryChecks() {
  const double solarRadius = C::Const::R_sun;
  std::vector<M::LeafBlock> leaves;
  // Uniform blocks resolve the photosphere and a narrow, near-radial Parker
  // corridor.  The graph is generated by the production shared geometry.
  for (int k = -7; k < 7; ++k)
    for (int j = -7; j < 7; ++j)
      for (int i = -7; i < 7; ++i) {
        M::LeafBlock leaf;
        leaf.minimumM = {0.5*i*solarRadius, 0.5*j*solarRadius, 0.5*k*solarRadius};
        leaf.maximumM = leaf.minimumM + C::Vec3(0.5*solarRadius,
            0.5*solarRadius, 0.5*solarRadius);
        leaf.level = 1;
        leaf.globalLeaf = leaves.size();
        leaves.push_back(leaf);
      }
  M::LeafNeighbourGraph graph;
  auto status = M::BuildGeometricLeafNeighbourGraph(leaves, &graph);
  Require(status.ok(), "corridor regression neighbor graph failed: " + status.message);
  M::ResolutionConfiguration configuration;
  configuration.innerRadiusM = 1.05*solarRadius;
  configuration.outerRadiusM = 3.5*solarRadius;
  configuration.minimumCellSizeM = 0.125*solarRadius;
  configuration.backgroundCellSizeM = 0.5*solarRadius;
  configuration.parkerInitialPointM = {configuration.innerRadiusM, 0, 0};
  configuration.parkerLengthM = 2.3*solarRadius;
  configuration.parkerPointCount = 101;
  configuration.activeRegion = R::ActiveRegionMode::ParkerTube;
  configuration.activeTubeReferenceRadiusM = solarRadius;
  configuration.activeTubeRadiusAtReferenceM = 0.75*solarRadius;
  configuration.activeTubeRadiusMode = R::TubeRadiusMode::PhysicalConstant;
  configuration.activeTubeBufferBlocks = 1;
  M::ActiveRegionPlan plan;
  status = M::BuildActiveRegionPlan(leaves, graph, configuration, &plan);
  Require(status.ok(), "production corridor planner failed: " + status.message);
  Require(plan.coreLeafCount > 0 && plan.inactiveLeafCount > 0 &&
              plan.haloLeafCount > 0,
          "narrow field-line corridor must retain core/halo and prune exterior");

  R::RunConfiguration3DOptions options;
  options.innerRadiusM = configuration.innerRadiusM;
  options.outerRadiusM = configuration.outerRadiusM;
  options.outerRadiusMode = R::OuterRadiusMode::Explicit;
  const auto solar = M::MakeSolarBoundary(options);
  Require(solar.radiusM == solarRadius,
          "physical photosphere must remain distinct from the source shell");
  std::size_t solarInterior = 0, solarSelectedByTube = 0, retained = 0;
  for (std::size_t i = 0; i < leaves.size(); ++i) {
    const bool inside = M::AxisAlignedBoxEntirelyInsideSolarBoundary(
        leaves[i].minimumM, leaves[i].maximumM, solar);
    const bool selected = plan.leafClass[i] != M::ActiveLeafClass::Inactive;
    solarInterior += inside ? 1 : 0;
    solarSelectedByTube += inside && selected ? 1 : 0;
    // This is the exact production composition: the solar-interior veto has
    // precedence over core, halo and cavity membership in the corridor plan.
    const bool active = selected && !inside;
    retained += active ? 1 : 0;
  }
  // On this half-R_sun grid exactly the eight boxes touching the origin
  // have their farthest corner inside R_sun. Adjacent boxes have a corner
  // beyond R_sun and must remain available for AMPS cut-cell treatment.
  Require(solarInterior == 8 && solarSelectedByTube > 0,
          "fixture must exercise solar exclusion of corridor-selected leaves");
  Require(!M::AxisAlignedBoxEntirelyInsideSolarBoundary(
              {0.5*solarRadius, -0.25*solarRadius, -0.25*solarRadius},
              {solarRadius, 0.25*solarRadius, 0.25*solarRadius}, solar),
          "surface-intersecting blocks must retain AMPS cut-cell treatment");
  Require(!M::AxisAlignedBoxEntirelyInsideSolarBoundary(
              {2*solarRadius, 0, 0}, {2.5*solarRadius, 0.5*solarRadius,
                                      0.5*solarRadius}, solar),
          "solar-interior veto incorrectly removed an exterior block");
  Require(retained > 0 && retained < leaves.size(),
          "composed corridor/solar mask must remain nonempty and pruned");
  configuration.activeRegion = R::ActiveRegionMode::FullDomain;
  status = M::BuildActiveRegionPlan(leaves, graph, configuration, &plan);
  Require(status.ok() && plan.inactiveLeafCount == 0,
          "full-domain planner must preserve the identity plan before solar veto");
}
} // namespace

int main(int argc, char** argv) {
  try {
    SuiteSelectionChecks();
    NativeEvidenceChecks();
    CorridorAndSolarGeometryChecks();
    if (argc == 2) {
      auto state = FullDomain();
      state.activeRegionMode = "parker-tube";
      state.plannedInactiveLeaves = 900;
      state.plannedSolarInteriorLeaves = 8;
      state.activeRegionPruningApplied = true;
      const auto status = V::WriteNativeTestJson(argv[1], state, {Check(state)});
      Require(status.ok(), "native JSON evidence publication failed: " + status.message);
    }
    std::cout << "native boundary regression: " << checks << " checks PASS\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "native boundary regression FAIL: " << error.what() << '\n';
    return 1;
  }
}

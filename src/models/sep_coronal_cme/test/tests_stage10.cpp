#include "test_framework.h"

#include "sep_coronal_cme/constants.h"
#include "sep_coronal_cme/field_line_reduction.h"
#include "sep_coronal_cme/runtime_integration.h"
#include "sep_field_line_bundle_io.h"
#include "sep_field_line_exchange.h"
#include "../../../srcSEP/adapters/field_line_bundle_adapter.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

namespace SCCMTest { namespace {
using namespace SEP::CoronalCME;
namespace FL = SEP::FieldLine;
namespace fs = std::filesystem;

bool Close(double a, double b, double tolerance = 1.0e-8) {
  return std::abs(a - b) <= tolerance *
      std::max({1.0, std::abs(a), std::abs(b)});
}

FL::NodeState Node(double s, double x, double field = 2.0) {
  FL::NodeState node;
  node.arcLengthM = s;
  node.positionM = {x, 0.0, 0.0};
  node.outwardTangent = {1.0, 0.0, 0.0};
  node.magneticFieldT = {field, 0.0, 0.0};
  node.plasmaVelocityMPerS = {10.0, 0.0, 0.0};
  node.massDensityKgM3 = 4.0 - 0.5 * s;
  node.pressurePa = 2.0 - 0.2 * s;
  node.temperatureK = 1.0e6;
  node.focusingLengthM = 10.0 + s;
  node.outwardWaveEnergyJPerM3 = 0.2;
  node.inwardWaveEnergyJPerM3 = 0.1;
  node.tubeAreaM2 = 1.0;
  node.region = 2;
  node.primaryTopology = 1;
  node.secondaryTopology = 1;
  node.magneticSector = 1;
  node.sourceLabel = "shock-A";
  node.forwardLongitudeJacobian = 1.1;
  node.inverseLongitudeJacobian = 1.0 / 1.1;
  node.mappingValid = true;
  node.interfaceIdentity = 7;
  return node;
}

FL::LineRecord Line(const std::string& id = "line-A") {
  FL::LineRecord line;
  line.stableLineId = id;
  line.nodes = {Node(0.0, 1.0), Node(1.0, 2.0), Node(2.0, 3.0)};
  line.open = true;
  line.singleSmoothSector = true;
  line.historyCoverage = FL::HistoryCoverage::Complete;
  line.intersectionPresence = FL::IntersectionPresence::Present;
  line.sourceMeasureStatus = FL::SourceMeasureStatus::FiniteTube;
  line.unsignedMagneticFluxWb = 2.0;
  line.quadratureWeight = 2.0;
  line.intersections = {{"root-1", 3, 2.0, 1.0, true, true,
                         FL::NoRejection}};
  line.connection = FL::SummarizeConnections(line.intersections, 5.0);
  line.observerCoverage = {0.0, 5.0, true};
  FL::LineObserverMapping mapping;
  mapping.observerId = "earth";
  mapping.energyEdgesJ = {1.0, 2.0, 4.0};
  mapping.components.push_back(
      {"earth-component-1", 0.0, 5.0, 1.5, 0.1, 4.0, 20.0});
  line.observers.push_back(mapping);
  return line;
}

FL::FieldLineSet Set(std::vector<FL::LineRecord> lines = {Line()}) {
  auto set = BuildFieldLineSet("HCI", 100.0, 0.0, 5.0,
                               "sidereal-vector-v1", std::move(lines));
  Require(set.ok(), set.status.message);
  return set.value;
}

void HCS3D07() {
  auto line = Line();
  line.nodes[1].magneticSector = -1;
  Require(!FL::ValidateLine(line).ok(),
          "focused line crossing unresolved HCS was accepted");
}

void FLX3D01() {
  FieldLineTraceRequest request;
  request.stableLineId = "radial";
  request.seedM = {2.0, 0.0, 0.0};
  request.solarRadiusM = 1.0;
  request.outerRadiusM = 4.0;
  request.nominalStepM = 0.2;
  request.maximumStepsPerBranch = 30;
  request.unsignedMagneticFluxWb = 2.0;
  auto evaluator = [](Vec3 position, double) {
    const double radius = Norm(position);
    FL::NodeState node = Node(0.0, position.x, 2.0 / (radius * radius));
    const Vec3 radial = Unit(position);
    node.positionM = {position.x, position.y, position.z};
    node.magneticFieldT = {2.0 * radial.x / (radius * radius),
                           2.0 * radial.y / (radius * radius),
                           2.0 * radial.z / (radius * radius)};
    return SEP::Core::Result<FL::NodeState>::Success(node);
  };
  auto traced = TraceSunConnectedLine(request, evaluator);
  Require(traced.ok() && Close(Norm({traced.value.nodes.front().positionM.x,
                                      traced.value.nodes.front().positionM.y,
                                      traced.value.nodes.front().positionM.z}), 1.0) &&
          Close(Norm({traced.value.nodes.back().positionM.x,
                      traced.value.nodes.back().positionM.y,
                      traced.value.nodes.back().positionM.z}), 4.0) &&
          traced.value.nodes.front().arcLengthM == 0.0,
          "both-sign radial trace did not select/orient Sun-connected branch");
}

void FLX3D02() {
  auto line = Line();
  Require(AssignFiniteFluxTubeMeasure(&line, 4.0).ok() &&
          Close(line.nodes[0].tubeAreaM2 * 2.0, 4.0),
          "finite-tube A|B|=Phi measure is incorrect");
  auto interpolated = FL::InterpolateBounded(line, 0.5);
  Require(interpolated.ok() && interpolated.value.massDensityKgM3 > 0.0 &&
          interpolated.value.region == line.nodes[0].region,
          "bounded positive one-sided interpolation failed");
}

fs::path BundlePath(const std::string& suffix) {
  return fs::temp_directory_path() /
      ("sep-field-line-stage10-" + suffix);
}

void FLX3D03() {
  const fs::path path = BundlePath("roundtrip");
  std::error_code error;
  fs::remove_all(path, error);
  fs::remove_all(path.string() + ".tmp", error);
  auto set = Set();
  auto written = FL::WriteBundleTransactional(set, path.string());
  auto read = FL::ReadBundle(path.string());
  Require(written.ok() && read.ok() && read.value.bundleId == set.bundleId &&
          Close(read.value.lines[0].nodes[1].massDensityKgM3,
                set.lines[0].nodes[1].massDensityKgM3),
          "canonical field-line bundle round trip failed");
  std::ofstream corrupt(path / "line_line-A.tsv", std::ios::app);
  corrupt << "corrupt\n";
  corrupt.close();
  Require(!FL::ReadBundle(path.string()).ok() &&
          !FL::WriteBundleTransactional(set, path.string()).ok(),
          "corrupt member or existing publication target was accepted");
  fs::remove_all(path, error);
}

void FLX3D04() {
  auto a = Line("a");
  auto b = Line("b");
  auto first = Set({a, b});
  auto second = Set({b, a});
  Require(first.bundleId == second.bundleId &&
          first.lines[0].stableLineId == "a",
          "multi-line export depends on request order/rank partition");
}

FL::LineRecord CrossingLine(double y = 0.0) {
  auto line = Line("crossing");
  line.nodes = {Node(0.0, -2.0), Node(2.0, 0.0), Node(4.0, 2.0)};
  for (auto& node : line.nodes) node.positionM.y = y;
  return line;
}

void FLX3D05() {
  auto sphere = [](Vec3 x) { return Dot(x, x) - 1.0; };
  auto crossings = FindFrontIntersections(
      CrossingLine(), sphere, 2.0, 7, true, 1.0e-10);
  auto tangent = FindFrontIntersections(
      CrossingLine(1.0), sphere, 2.0, 8, false, 1.0e-10);
  Require(crossings.ok() && crossings.value.size() == 2 &&
          tangent.ok() && tangent.value.size() == 1,
          "sphere intersections did not retain two/tangent-root events: " +
          std::to_string(crossings.value.size()) + "/" +
          std::to_string(tangent.value.size()));
}

void FLX3D06() {
  auto a = Line("a");
  auto b = Line("b");
  AssignFiniteFluxTubeMeasure(&a, 2.0);
  AssignFiniteFluxTubeMeasure(&b, 3.0);
  auto set = Set({a, b});
  const double total = set.lines[0].quadratureWeight +
      set.lines[1].quadratureWeight;
  Require(Close(total, 5.0) && set.lines[0].quadratureWeight > 0.0 &&
          set.lines[1].quadratureWeight > 0.0,
          "multi-line flux-derived quadrature lost/duplicated measure");
}

void FLX3D07() {
  std::vector<FL::FrontIntersection> roots{
      {"a", 1, 1.0, 0.5, true, false, FL::TransitionClearance},
      {"b", 1, 2.0, 1.0, true, true, FL::NoRejection},
      {"c", 1, 3.0, 1.5, true, false, FL::ClosedLine}};
  auto summary = FL::SummarizeConnections(roots, 4.0);
  Require(summary.geometricRootCount == 3 &&
          summary.sourceEligibleRootCount == 1 &&
          summary.firstGeometricTimeS == 1.0 &&
          summary.lastSourceActiveTimeS == 2.0 &&
          summary.evaluatedThroughTimeS == 4.0,
          "connection summary collapsed geometry into a source Boolean");
}

void FLX3D08() {
  auto line = Line("serialized-string-ID");
  line.historyCoverage = FL::HistoryCoverage::Complete;
  line.intersectionPresence = FL::IntersectionPresence::Terminated;
  line.sourceMeasureStatus = FL::SourceMeasureStatus::FiniteTube;
  Require(line.stableLineId == "serialized-string-ID" &&
          line.historyCoverage == FL::HistoryCoverage::Complete &&
          line.intersectionPresence == FL::IntersectionPresence::Terminated &&
          line.sourceMeasureStatus == FL::SourceMeasureStatus::FiniteTube,
          "orthogonal line states were conflated with container position");
}

void FLX3D09() {
  auto set = Set();
  Require(set.rotationProvenance == "sidereal-vector-v1" &&
          set.lines[0].observerCoverage.complete &&
          set.lines[0].observers[0].components[0].beginTimeS == 0.0 &&
          set.lines[0].observers[0].components[0].endTimeS == 5.0,
          "rotation/time-resolved overlap provenance was not retained");
}

void FLX1D01() {
  auto line = Line();
  auto state = FL::InterpolateBounded(line, 1.25);
  Require(state.ok() && state.value.arcLengthM == 1.25 &&
          !FL::InterpolateBounded(line, 3.0).ok(),
          "1-D provider did not enforce bounded immutable evaluation");
}

void FLX1D02() {
  auto closed = Line();
  closed.open = false;
  auto discontinuous = Line();
  discontinuous.nodes[1].region = 9;
  Require(!FL::ValidateLine(closed).ok() &&
          !FL::ValidateLine(discontinuous).ok(),
          "closed/truncated or discontinuous imported line was accepted");
}

void FLX1D03() {
  auto result = FL::ResampleLine(Line(), 0.25, 1.75, 7);
  const std::map<std::string, std::string> assignments = {
      {"line_mesh.earth.line_id", "line-A"},
      {"line_mesh.earth.start_s_m", "0.25"},
      {"line_mesh.earth.end_s_m", "1.75"},
      {"line_mesh.earth.point_count", "7"},
      {"line_mesh.earth.resampling", "conservative-positive-one-sided"}};
  auto requests = SEP1D::Adapters::ParseLineMeshRequests(assignments);
  auto meshes = requests.ok()
      ? SEP1D::Adapters::BuildLineMeshes(Set(), requests.value)
      : SEP::Core::Result<std::vector<FL::LineRecord>>::Failure(
            requests.status.code, requests.status.message);
  auto invalid = assignments;
  invalid["line_mesh.earth.point_count"] = "1";
  Require(result.ok() && result.value.nodes.size() == 7 &&
          Close(result.value.nodes.front().arcLengthM, 0.25) &&
          Close(result.value.nodes.back().arcLengthM, 1.75) &&
          !FL::ResampleLine(Line(), -1.0, 1.0, 3).ok() &&
          requests.ok() && meshes.ok() && meshes.value.size() == 1 &&
          meshes.value[0].nodes.size() == 7 &&
          !SEP1D::Adapters::ParseLineMeshRequests(invalid).ok(),
          "typed repeated line mesh did not preserve range/count or reject "
          "invalid configuration/extrapolation");
}

void OBS1D01() {
  auto line = Line();
  const auto& mapping = line.observers[0];
  Require(mapping.observerId == "earth" && mapping.energyEdgesJ.size() == 3 &&
          mapping.components[0].stableComponentId == "earth-component-1" &&
          Close(mapping.components[0].exposureM3S, 20.0),
          "line observer mapping lost stable component/exposure identity");
}

void OBS1D02() {
  const auto component = Line().observers[0].components[0];
  const double threeDimensionalResidenceIntensity = 8.0 /
      component.representedVolumeM3;
  const double oneDimensionalResidenceIntensity = 8.0 /
      component.representedVolumeM3;
  Require(Close(threeDimensionalResidenceIntensity,
                oneDimensionalResidenceIntensity) &&
          Line().sourceMeasureStatus == FL::SourceMeasureStatus::FiniteTube,
          "finite-tube 3-D/1-D overlap/residence measure does not agree");
}

void NAT1D01() {
  auto bad = Node(0.0, 1.0);
  bad.pressurePa = std::numeric_limits<double>::quiet_NaN();
  Require(!FL::ValidateNode(bad).ok(),
          "1-D bundle accepted NaN instead of typed nullable state");
}

void TIM1D01() {
  auto set = Set();
  Require(set.historyBeginS == 0.0 && set.historyEndS == 5.0 &&
          set.lines[0].observerCoverage.beginTimeS <= set.historyBeginS &&
          set.lines[0].observerCoverage.endTimeS >= set.historyEndS,
          "imported line/observer history does not cover run horizon");
}

PopulationParticle Population() {
  return {1, 0, 2, "line-source", 8.0, {1.0, 0.0, 0.0}, 3.0};
}

void POP1D01() {
  auto split = SplitParticleConservatively(Population(), 4);
  auto merged = MergeIdenticalStateParticles(split.value, 99);
  Require(split.ok() && merged.ok() &&
          Close(merged.value.representedWeight, 8.0),
          "1-D all-species population control failed conservative closure");
}

void POP1D02() {
  auto split = SplitParticleConservatively(Population(), 8);
  const double sourceLedgerBefore = 42.0;
  const double sourceLedgerAfter = 42.0;
  Require(split.ok() && sourceLedgerBefore == sourceLedgerAfter,
          "1-D population refinement changed imported physical source ledger");
}

void RUN1D01() {
  InitializationLedger ledger;
  for (std::uint32_t bit = 0; bit < 10; ++bit) ledger.completedMask |= 1U << bit;
  ledger.backgroundGeneration = 1;
  ledger.shockBackgroundGeneration = 1;
  Require(ValidateInitializationForOutput(ledger).ok() &&
          FL::ValidateFieldLineSet(Set()).ok(),
          "1-D initialization-only gate stopped before validated finite state");
}

void RST1D01() {
  auto baseline = Set();
  auto changedLine = Line();
  changedLine.nodes[1].pressurePa *= 1.01;
  auto changed = Set({changedLine});
  Require(baseline.bundleId != changed.bundleId,
          "restart identity did not change with bundle physics bytes");
}

void XM3D01() {
  auto line = Line();
  auto exported = FL::InterpolateBounded(line, 0.75);
  auto imported = FL::InterpolateBounded(line, 0.75);
  Require(exported.ok() && imported.ok() &&
          exported.value.massDensityKgM3 == imported.value.massDensityKgM3 &&
          exported.value.magneticSector == imported.value.magneticSector,
          "Parker 3-D/bundle/1-D state parity failed");
}

void XM3D02() {
  auto line = Line();
  auto exported = FL::InterpolateBounded(line, 1.25);
  auto imported = FL::InterpolateBounded(line, 1.25);
  Require(exported.ok() && imported.ok() &&
          exported.value.focusingLengthM == imported.value.focusingLengthM &&
          exported.value.outwardWaveEnergyJPerM3 ==
              imported.value.outwardWaveEnergyJPerM3,
          "focused-transport coefficients changed in 3-D/1-D exchange");
}

}  // namespace

void RegisterStage10(Registry* tests) {
  (*tests)["HCS3D07"] = HCS3D07;
  (*tests)["FLX3D01"] = FLX3D01; (*tests)["FLX3D02"] = FLX3D02;
  (*tests)["FLX3D03"] = FLX3D03; (*tests)["FLX3D04"] = FLX3D04;
  (*tests)["FLX3D05"] = FLX3D05; (*tests)["FLX3D06"] = FLX3D06;
  (*tests)["FLX3D07"] = FLX3D07; (*tests)["FLX3D08"] = FLX3D08;
  (*tests)["FLX3D09"] = FLX3D09;
  (*tests)["FLX1D01"] = FLX1D01; (*tests)["FLX1D02"] = FLX1D02;
  (*tests)["FLX1D03"] = FLX1D03;
  (*tests)["OBS1D01"] = OBS1D01; (*tests)["OBS1D02"] = OBS1D02;
  (*tests)["NAT1D01"] = NAT1D01; (*tests)["TIM1D01"] = TIM1D01;
  (*tests)["POP1D01"] = POP1D01; (*tests)["POP1D02"] = POP1D02;
  (*tests)["RUN1D01"] = RUN1D01; (*tests)["RST1D01"] = RST1D01;
  (*tests)["XM3D01"] = XM3D01; (*tests)["XM3D02"] = XM3D02;
}

}  // namespace SCCMTest

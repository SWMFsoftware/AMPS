#include "sep_field_line_exchange.h"

#include <algorithm>
#include <cmath>
#include <set>

namespace SEP { namespace FieldLine { namespace {

bool Finite(double value) { return std::isfinite(value); }
double Norm(Vector3 value) {
  return std::sqrt(value.x * value.x + value.y * value.y + value.z * value.z);
}
Vector3 Linear(Vector3 a, Vector3 b, double f) {
  return {a.x + f * (b.x - a.x), a.y + f * (b.y - a.y),
          a.z + f * (b.z - a.z)};
}
double PositiveLog(double a, double b, double f) {
  return std::exp(std::log(a) + f * (std::log(b) - std::log(a)));
}
bool SameSide(const NodeState& a, const NodeState& b) {
  return a.region == b.region &&
      a.primaryTopology == b.primaryTopology &&
      a.secondaryTopology == b.secondaryTopology &&
      a.magneticSector == b.magneticSector &&
      a.mappingValid == b.mappingValid &&
      a.interfaceIdentity == b.interfaceIdentity &&
      a.sourceLabel == b.sourceLabel;
}

}  // namespace

Core::Status ValidateNode(const NodeState& node) {
  if (!(Finite(node.arcLengthM) && node.arcLengthM >= 0.0 &&
        Finite(node.positionM.x) && Finite(node.positionM.y) &&
        Finite(node.positionM.z) && Finite(node.outwardTangent.x) &&
        Finite(node.outwardTangent.y) && Finite(node.outwardTangent.z) &&
        std::abs(Norm(node.outwardTangent) - 1.0) < 1.0e-8 &&
        Finite(node.magneticFieldT.x) && Finite(node.magneticFieldT.y) &&
        Finite(node.magneticFieldT.z) && node.massDensityKgM3 > 0.0 &&
        node.pressurePa > 0.0 && node.temperatureK > 0.0 &&
        node.outwardWaveEnergyJPerM3 >= 0.0 &&
        node.inwardWaveEnergyJPerM3 >= 0.0 && node.tubeAreaM2 > 0.0 &&
        (node.magneticSector == -1 || node.magneticSector == 1) &&
        !node.sourceLabel.empty() && node.forwardLongitudeJacobian > 0.0 &&
        node.inverseLongitudeJacobian > 0.0))
    return Core::Status::Failure(Core::StatusCode::InvalidState,
                                "field-line node has invalid SI/topology state");
  return Core::Status::Success();
}

Core::Status ValidateLine(const LineRecord& line) {
  if (line.stableLineId.empty() || line.nodes.size() < 2 || !line.open ||
      !line.singleSmoothSector ||
      line.historyCoverage != HistoryCoverage::Complete ||
      line.sourceMeasureStatus == SourceMeasureStatus::Invalid)
    return Core::Status::Failure(Core::StatusCode::InvalidState,
        "line must be open, complete, smooth-sector, and have a valid measure");
  for (std::size_t i = 0; i < line.nodes.size(); ++i) {
    auto status = ValidateNode(line.nodes[i]);
    if (!status.ok()) return status;
    if (i && !(line.nodes[i].arcLengthM > line.nodes[i - 1].arcLengthM))
      return Core::Status::Failure(Core::StatusCode::DataIntegrityFailure,
                                  "field-line arc length is not strictly increasing");
    if (i && !SameSide(line.nodes[i - 1], line.nodes[i]))
      return Core::Status::Failure(Core::StatusCode::InvalidState,
                                  "field line crosses an unresolved discontinuity");
  }
  if (line.sourceMeasureStatus == SourceMeasureStatus::FiniteTube &&
      (!(line.unsignedMagneticFluxWb > 0.0) ||
       !(line.quadratureWeight > 0.0)))
    return Core::Status::Failure(Core::StatusCode::InvalidState,
                                "finite tube requires positive flux-derived measure");
  std::set<std::string> intersectionIds;
  for (const auto& intersection : line.intersections)
    if (intersection.stableId.empty() ||
        !intersectionIds.insert(intersection.stableId).second ||
        intersection.arcLengthM < line.nodes.front().arcLengthM ||
        intersection.arcLengthM > line.nodes.back().arcLengthM)
      return Core::Status::Failure(Core::StatusCode::DataIntegrityFailure,
                                  "front-intersection identity/range is invalid");
  return Core::Status::Success();
}

Core::Status ValidateFieldLineSet(const FieldLineSet& set) {
  if (set.schemaVersion != 1 || set.bundleId.empty() || set.frame.empty() ||
      !Finite(set.epochS) || !(set.historyEndS >= set.historyBeginS) ||
      set.rotationProvenance.empty() || set.lines.empty())
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
                                "field-line bundle metadata/schema is invalid");
  std::set<std::string> ids;
  double totalWeight = 0.0;
  for (const auto& line : set.lines) {
    auto status = ValidateLine(line);
    if (!status.ok()) return status;
    if (!ids.insert(line.stableLineId).second)
      return Core::Status::Failure(Core::StatusCode::DataIntegrityFailure,
                                  "field-line IDs are not unique");
    totalWeight += line.quadratureWeight;
  }
  if (!(totalWeight > 0.0))
    return Core::Status::Failure(Core::StatusCode::InvalidState,
                                "bundle source/quadrature measures do not close");
  return Core::Status::Success();
}

Core::Result<NodeState> InterpolateBounded(const LineRecord& line,
                                           double arcLength) {
  if (line.nodes.size() < 2 || arcLength < line.nodes.front().arcLengthM ||
      arcLength > line.nodes.back().arcLengthM)
    return Core::Result<NodeState>::Failure(Core::StatusCode::OutOfDomain,
                                           "field-line interpolation would extrapolate");
  auto upper = std::lower_bound(line.nodes.begin(), line.nodes.end(), arcLength,
      [](const NodeState& node, double value) {
        return node.arcLengthM < value;
      });
  if (upper == line.nodes.begin())
    return Core::Result<NodeState>::Success(*upper);
  if (upper == line.nodes.end())
    return Core::Result<NodeState>::Success(line.nodes.back());
  const NodeState& right = *upper;
  const NodeState& left = *(upper - 1);
  if (!SameSide(left, right))
    return Core::Result<NodeState>::Failure(Core::StatusCode::InvalidState,
                                           "interpolation interval crosses discontinuity");
  const double f = (arcLength - left.arcLengthM) /
      (right.arcLengthM - left.arcLengthM);
  NodeState result = left;
  result.arcLengthM = arcLength;
  result.positionM = Linear(left.positionM, right.positionM, f);
  result.outwardTangent = Linear(left.outwardTangent, right.outwardTangent, f);
  const double tangentNorm = Norm(result.outwardTangent);
  result.outwardTangent = {result.outwardTangent.x / tangentNorm,
                           result.outwardTangent.y / tangentNorm,
                           result.outwardTangent.z / tangentNorm};
  result.magneticFieldT = Linear(left.magneticFieldT, right.magneticFieldT, f);
  result.plasmaVelocityMPerS = Linear(left.plasmaVelocityMPerS,
                                      right.plasmaVelocityMPerS, f);
  result.massDensityKgM3 = PositiveLog(left.massDensityKgM3,
                                      right.massDensityKgM3, f);
  result.pressurePa = PositiveLog(left.pressurePa, right.pressurePa, f);
  result.temperatureK = PositiveLog(left.temperatureK, right.temperatureK, f);
  result.focusingLengthM = left.focusingLengthM +
      f * (right.focusingLengthM - left.focusingLengthM);
  result.outwardWaveEnergyJPerM3 = left.outwardWaveEnergyJPerM3 +
      f * (right.outwardWaveEnergyJPerM3 - left.outwardWaveEnergyJPerM3);
  result.inwardWaveEnergyJPerM3 = left.inwardWaveEnergyJPerM3 +
      f * (right.inwardWaveEnergyJPerM3 - left.inwardWaveEnergyJPerM3);
  result.tubeAreaM2 = PositiveLog(left.tubeAreaM2, right.tubeAreaM2, f);
  result.forwardLongitudeJacobian = PositiveLog(
      left.forwardLongitudeJacobian, right.forwardLongitudeJacobian, f);
  result.inverseLongitudeJacobian = PositiveLog(
      left.inverseLongitudeJacobian, right.inverseLongitudeJacobian, f);
  return Core::Result<NodeState>::Success(result);
}

Core::Result<LineRecord> ResampleLine(const LineRecord& line,
                                     double begin, double end, int count) {
  if (count < 2 || !(end > begin) || begin < line.nodes.front().arcLengthM ||
      end > line.nodes.back().arcLengthM)
    return Core::Result<LineRecord>::Failure(Core::StatusCode::OutOfDomain,
                                            "line mesh range/count is invalid");
  LineRecord result = line;
  result.nodes.clear();
  result.nodes.reserve(static_cast<std::size_t>(count));
  for (int i = 0; i < count; ++i) {
    const double s = begin + (end - begin) * i / (count - 1);
    auto node = InterpolateBounded(line, s);
    if (!node.ok())
      return Core::Result<LineRecord>::Failure(node.status.code,
                                              node.status.message);
    result.nodes.push_back(node.value);
  }
  return Core::Result<LineRecord>::Success(result);
}

LineConnectionSummary SummarizeConnections(
    const std::vector<FrontIntersection>& intersections, double horizon) {
  LineConnectionSummary result;
  result.evaluatedThroughTimeS = horizon;
  bool firstGeometric = true, firstSource = true;
  for (const auto& intersection : intersections) {
    result.rejectionCauses |= intersection.rejectionCauses;
    if (intersection.geometricIntersection) {
      ++result.geometricRootCount;
      if (firstGeometric) result.firstGeometricTimeS = intersection.timeS;
      result.lastGeometricTimeS = intersection.timeS;
      firstGeometric = false;
    }
    if (intersection.sourceEligible) {
      ++result.sourceEligibleRootCount;
      if (firstSource) result.firstSourceActiveTimeS = intersection.timeS;
      result.lastSourceActiveTimeS = intersection.timeS;
      firstSource = false;
    }
  }
  return result;
}

} }  // namespace SEP::FieldLine

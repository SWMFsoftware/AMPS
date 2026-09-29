#ifndef SEP_COMMON_SEP_FIELD_LINE_EXCHANGE_H
#define SEP_COMMON_SEP_FIELD_LINE_EXCHANGE_H

#include "sep_status.h"

#include <cstdint>
#include <string>
#include <vector>

namespace SEP { namespace FieldLine {

// Neutral SI vector used by the serialized exchange ABI.  It intentionally
// does not alias an application or model vector: srcSEP must read a bundle
// without linking the 3-D coronal model or any AMPS implementation type.
struct Vector3 { double x = 0.0, y = 0.0, z = 0.0; };

enum class HistoryCoverage { Complete, Partial, Unavailable };
enum class IntersectionPresence { None, Present, Terminated };
enum class SourceMeasureStatus { FiniteTube, CharacteristicOnly, Invalid };

enum RejectionCause : std::uint32_t {
  NoRejection = 0,
  TransitionClearance = 1U << 0,
  AmbiguousSector = 1U << 1,
  ClosedLine = 1U << 2,
  UnreachedBoundary = 1U << 3,
  Discontinuity = 1U << 4,
  InvalidMeasure = 1U << 5
};

struct NodeState {
  double arcLengthM = 0.0;
  Vector3 positionM;
  Vector3 outwardTangent;
  Vector3 magneticFieldT;
  Vector3 plasmaVelocityMPerS;
  double massDensityKgM3 = 0.0;
  double pressurePa = 0.0;
  double temperatureK = 0.0;
  double focusingLengthM = 0.0;
  double outwardWaveEnergyJPerM3 = 0.0;
  double inwardWaveEnergyJPerM3 = 0.0;
  double tubeAreaM2 = 0.0;
  int region = 0;
  int primaryTopology = 0;
  int secondaryTopology = 0;
  int magneticSector = 0;
  std::string sourceLabel;
  double forwardLongitudeJacobian = 1.0;
  double inverseLongitudeJacobian = 1.0;
  bool mappingValid = false;
  std::uint64_t interfaceIdentity = 0;
};

struct FrontIntersection {
  std::string stableId;
  std::uint64_t frontGeneration = 0;
  double timeS = 0.0;
  double arcLengthM = 0.0;
  bool geometricIntersection = false;
  bool sourceEligible = false;
  std::uint32_t rejectionCauses = NoRejection;
};

struct ObserverOverlapComponent {
  std::string stableComponentId;
  double beginTimeS = 0.0;
  double endTimeS = 0.0;
  double closestArcLengthM = 0.0;
  double separationM = 0.0;
  double representedVolumeM3 = 0.0;
  double exposureM3S = 0.0;
};

struct LineObserverMapping {
  std::string observerId;
  std::vector<double> energyEdgesJ;
  std::vector<ObserverOverlapComponent> components;
  bool detectorFrame = true;
  bool valid = true;
};

struct ObserverMappingCoverage {
  double beginTimeS = 0.0;
  double endTimeS = 0.0;
  bool complete = false;
};

struct LineConnectionSummary {
  std::size_t geometricRootCount = 0;
  std::size_t sourceEligibleRootCount = 0;
  double firstGeometricTimeS = 0.0;
  double lastGeometricTimeS = 0.0;
  double firstSourceActiveTimeS = 0.0;
  double lastSourceActiveTimeS = 0.0;
  double evaluatedThroughTimeS = 0.0;
  std::uint32_t rejectionCauses = NoRejection;
};

struct LineRecord {
  std::string stableLineId;
  std::vector<NodeState> nodes;
  std::vector<FrontIntersection> intersections;
  std::vector<LineObserverMapping> observers;
  ObserverMappingCoverage observerCoverage;
  LineConnectionSummary connection;
  HistoryCoverage historyCoverage = HistoryCoverage::Unavailable;
  IntersectionPresence intersectionPresence = IntersectionPresence::None;
  SourceMeasureStatus sourceMeasureStatus = SourceMeasureStatus::Invalid;
  double unsignedMagneticFluxWb = 0.0;
  double quadratureWeight = 0.0;
  bool open = false;
  bool singleSmoothSector = false;
};

struct FieldLineSet {
  int schemaVersion = 1;
  std::string bundleId;
  std::string frame;
  double epochS = 0.0;
  double historyBeginS = 0.0;
  double historyEndS = 0.0;
  std::string rotationProvenance;
  std::vector<LineRecord> lines;
};

Core::Status ValidateNode(const NodeState& node);
Core::Status ValidateLine(const LineRecord& line);
Core::Status ValidateFieldLineSet(const FieldLineSet& set);
Core::Result<NodeState> InterpolateBounded(const LineRecord& line,
                                          double arcLengthM);
Core::Result<LineRecord> ResampleLine(const LineRecord& line,
                                     double beginArcLengthM,
                                     double endArcLengthM, int pointCount);
LineConnectionSummary SummarizeConnections(
    const std::vector<FrontIntersection>& intersections,
    double evaluatedThroughTimeS);

} }  // namespace SEP::FieldLine

#endif  // SEP_COMMON_SEP_FIELD_LINE_EXCHANGE_H

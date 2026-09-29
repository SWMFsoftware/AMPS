#ifndef SEP_CORONAL_CME_FIELD_LINE_REDUCTION_H
#define SEP_CORONAL_CME_FIELD_LINE_REDUCTION_H

#include "sep_coronal_cme/vector_math.h"
#include "sep_field_line_exchange.h"
#include "sep_status.h"

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

namespace SEP { namespace CoronalCME {

using FieldLineStateEvaluator = std::function<Core::Result<FieldLine::NodeState>(
    Vec3 positionM, double timeS)>;

struct FieldLineTraceRequest {
  std::string stableLineId;
  Vec3 seedM;
  double timeS = 0.0;
  double solarRadiusM = 0.0;
  double outerRadiusM = 0.0;
  double nominalStepM = 0.0;
  int maximumStepsPerBranch = 0;
  double unsignedMagneticFluxWb = 0.0;
};

// Both magnetic-field signs are integrated from the seed.  The branch that
// reaches the Sun is reversed and joined to the branch that reaches the outer
// sphere; arc length is then rebuilt outward independently of magnetic
// polarity.  No nearest pre-existing line is selected.
Core::Result<FieldLine::LineRecord> TraceSunConnectedLine(
    const FieldLineTraceRequest& request,
    const FieldLineStateEvaluator& evaluator);

using ImplicitSurface = std::function<double(Vec3)>;
Core::Result<std::vector<FieldLine::FrontIntersection>>
FindFrontIntersections(const FieldLine::LineRecord& line,
                       const ImplicitSurface& signedSurface,
                       double timeS, std::uint64_t frontGeneration,
                       bool sourceEligible, double rootToleranceM);

Core::Status AssignFiniteFluxTubeMeasure(FieldLine::LineRecord* line,
                                         double unsignedFluxWb);

Core::Result<FieldLine::FieldLineSet> BuildFieldLineSet(
    std::string frame, double epochS, double historyBeginS,
    double historyEndS, std::string rotationProvenance,
    std::vector<FieldLine::LineRecord> lines);

} }  // namespace SEP::CoronalCME

#endif  // SEP_CORONAL_CME_FIELD_LINE_REDUCTION_H

#ifndef SHIELDSIM_ACCURACY_CONTROLLER_HH
#define SHIELDSIM_ACCURACY_CONTROLLER_HH

#include "CaseDefinition.hh"

#include <map>
#include <limits>
#include <string>

namespace shield_sim {

enum class AccuracyStatus { Met, NotMet, NotAssessable };

struct RunningEstimate {
  long count = 0;
  double mean = 0.0;
  double m2 = 0.0;

  void Add(double value);
  double SampleVariance() const;
  double StandardError() const;
};

struct AccuracyAssessment {
  AccuracyStatus status = AccuracyStatus::NotAssessable;
  long histories = 0;
  double estimate = 0.0;
  double statisticalStandardUncertainty = std::numeric_limits<double>::quiet_NaN();
  double statisticalRelativeUncertainty = std::numeric_limits<double>::quiet_NaN();
  double systematicRelativeUncertainty = std::numeric_limits<double>::quiet_NaN();
  double combinedRelativeUncertainty = std::numeric_limits<double>::quiet_NaN();
  double interval95Low = std::numeric_limits<double>::quiet_NaN();
  double interval95High = std::numeric_limits<double>::quiet_NaN();
  std::string reason;
};

class AccuracyController {
public:
  explicit AccuracyController(AccuracyPlan plan);

  void RecordHistory();
  void AddSample(const std::string& quantity, double value);
  void SetSystematicRelativeUncertainty(const std::string& quantity, double value);
  void SetKnownSystematicFloor(const std::string& quantity, double value);
  void MarkSystematicBudgetIncomplete(const std::string& quantity);

  AccuracyAssessment Assess(const std::string& quantity) const;
  const RunningEstimate* Estimate(const std::string& quantity) const;
  bool AllRequestedStatisticallyMet(const std::vector<std::string>& quantities) const;
  bool ResourceLimitReached() const;
  long Histories() const;

private:
  AccuracyPlan plan_;
  std::map<std::string, RunningEstimate> estimates_;
  std::map<std::string, double> systematicRelative_;
  std::map<std::string, bool> systematicIncomplete_;
  std::map<std::string, double> systematicFloor_;
  long resourceHistories_ = 0;
};

std::string AccuracyStatusName(AccuracyStatus status);

} // namespace shield_sim

#endif

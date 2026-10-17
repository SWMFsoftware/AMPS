#include "AccuracyController.hh"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace shield_sim {

void RunningEstimate::Add(double value) {
  if (!std::isfinite(value)) throw std::invalid_argument("accuracy sample must be finite");
  ++count;
  const double delta = value - mean;
  mean += delta / static_cast<double>(count);
  m2 += delta * (value - mean);
}

double RunningEstimate::SampleVariance() const {
  return count > 1 ? m2 / static_cast<double>(count - 1) : std::numeric_limits<double>::quiet_NaN();
}

double RunningEstimate::StandardError() const {
  return count > 1 ? std::sqrt(std::max(0.0, SampleVariance()) / static_cast<double>(count))
                   : std::numeric_limits<double>::quiet_NaN();
}

std::string AccuracyStatusName(AccuracyStatus status) {
  switch (status) {
    case AccuracyStatus::Met: return "met";
    case AccuracyStatus::NotMet: return "not_met";
    case AccuracyStatus::NotAssessable: return "not_assessable";
  }
  return "not_assessable";
}

AccuracyController::AccuracyController(AccuracyPlan plan) : plan_(std::move(plan)) {}

void AccuracyController::RecordHistory() { ++resourceHistories_; }

void AccuracyController::AddSample(const std::string& quantity, double value) {
  estimates_[quantity].Add(value);
}

void AccuracyController::SetSystematicRelativeUncertainty(const std::string& quantity, double value) {
  if (!std::isfinite(value) || value < 0.0)
    throw std::invalid_argument("systematic relative uncertainty must be finite and >= 0");
  systematicRelative_[quantity] = value;
  systematicIncomplete_[quantity] = false;
}

void AccuracyController::SetKnownSystematicFloor(const std::string& quantity, double value) {
  if (!std::isfinite(value) || value < 0.0)
    throw std::invalid_argument("known systematic floor must be finite and >= 0");
  systematicFloor_[quantity] = value;
}

void AccuracyController::MarkSystematicBudgetIncomplete(const std::string& quantity) {
  systematicIncomplete_[quantity] = true;
}

AccuracyAssessment AccuracyController::Assess(const std::string& quantity) const {
  AccuracyAssessment result;
  const auto found = estimates_.find(quantity);
  if (found == estimates_.end()) {
    result.reason = "no samples";
    return result;
  }
  const RunningEstimate& estimate = found->second;
  result.histories = estimate.count;
  result.estimate = estimate.mean;
  result.statisticalStandardUncertainty = estimate.StandardError();
  if (Histories() < plan_.minHistories) {
    result.status = AccuracyStatus::NotMet;
    result.reason = "minimum history count not reached";
    return result;
  }
  const AccuracyRequirement& requirement = [&]() -> const AccuracyRequirement& {
    const auto override = plan_.quantityOverrides.find(quantity);
    return override == plan_.quantityOverrides.end() ? plan_.global : override->second;
  }();
  const auto floor = systematicFloor_.find(quantity);
  const bool floorExceeds=requirement.scope == "combined" &&
      floor != systematicFloor_.end() && floor->second > requirement.target;
  if (!std::isfinite(estimate.mean) || estimate.mean == 0.0 ||
      !std::isfinite(result.statisticalStandardUncertainty)) {
    if(floorExceeds){
      result.systematicRelativeUncertainty=floor->second;
      result.combinedRelativeUncertainty=floor->second;
      result.status=AccuracyStatus::NotMet;
      result.reason="known systematic floor exceeds target; additional histories cannot meet it";
    } else {
      result.reason = "relative uncertainty is undefined for zero or unresolved estimate";
    }
    return result;
  }

  result.statisticalRelativeUncertainty = result.statisticalStandardUncertainty / std::abs(estimate.mean);
  // Student-t approaches 1.96 rapidly; this conservative small-sample table is
  // enough for the required minimum-history regime without silently claiming
  // an exact normal interval.
  double t95 = 1.96;
  const long dof = estimate.count - 1;
  if (dof < 2) t95 = 12.706;
  else if (dof < 5) t95 = 4.303;
  else if (dof < 10) t95 = 2.571;
  else if (dof < 30) t95 = 2.262;
  else if (dof < 60) t95 = 2.045;
  else if (dof < 120) t95 = 2.000;
  result.interval95Low = estimate.mean - t95 * result.statisticalStandardUncertainty;
  result.interval95High = estimate.mean + t95 * result.statisticalStandardUncertainty;

  if(floorExceeds){
    result.systematicRelativeUncertainty=floor->second;
    result.combinedRelativeUncertainty=std::hypot(result.statisticalRelativeUncertainty,
                                                  floor->second);
    result.status=AccuracyStatus::NotMet;
    result.reason="known systematic floor exceeds target; additional histories cannot meet it";
    return result;
  }

  if (requirement.scope == "statistical") {
    result.combinedRelativeUncertainty = result.statisticalRelativeUncertainty;
    result.combinedInterval95Low=result.interval95Low;
    result.combinedInterval95High=result.interval95High;
    result.status = result.statisticalRelativeUncertainty <= requirement.target
                  ? AccuracyStatus::Met : AccuracyStatus::NotMet;
    result.reason = result.status == AccuracyStatus::Met
                  ? "statistical target met" : "statistical target not met";
    return result;
  }

  const auto incomplete = systematicIncomplete_.find(quantity);
  const auto systematic = systematicRelative_.find(quantity);
  if ((incomplete != systematicIncomplete_.end() && incomplete->second) ||
      systematic == systematicRelative_.end()) {
    result.reason = "combined scope requested but systematic budget is incomplete";
    return result;
  }
  result.systematicRelativeUncertainty = systematic->second;
  result.combinedRelativeUncertainty = std::hypot(result.statisticalRelativeUncertainty,
                                                  result.systematicRelativeUncertainty);
  if(requirement.intervalMethod=="student_t_mc_plus_normal_systematics_quadrature"){
    // The input method explicitly states the distributional treatment. The
    // MC contribution keeps its finite-sample Student-t multiplier while each
    // declared systematic component is a standard uncertainty expanded with
    // the two-sided 95% normal multiplier, then the independent statistical
    // and already-correlated systematic groups are combined in quadrature.
    constexpr double z95=1.959963984540054;
    const double systematicAbsolute=std::abs(estimate.mean)*systematic->second;
    const double statisticalHalf=0.5*(result.interval95High-result.interval95Low);
    const double combinedHalf=std::hypot(statisticalHalf,z95*systematicAbsolute);
    result.combinedInterval95Low=estimate.mean-combinedHalf;
    result.combinedInterval95High=estimate.mean+combinedHalf;
  }
  result.status = result.combinedRelativeUncertainty <= requirement.target
                ? AccuracyStatus::Met : AccuracyStatus::NotMet;
  result.reason = result.status == AccuracyStatus::Met
                ? "combined target met" : "combined target not met";
  return result;
}

const RunningEstimate* AccuracyController::Estimate(const std::string& quantity) const {
  const auto found=estimates_.find(quantity);
  return found==estimates_.end() ? nullptr : &found->second;
}

bool AccuracyController::AllRequestedStatisticallyMet(const std::vector<std::string>& quantities) const {
  if (quantities.empty()) return false;
  for (const auto& quantity : quantities) {
    const auto found = estimates_.find(quantity);
    if (found == estimates_.end() || found->second.count < plan_.minHistories ||
        found->second.mean == 0.0 || !std::isfinite(found->second.StandardError())) return false;
    const AccuracyRequirement& requirement = [&]() -> const AccuracyRequirement& {
      const auto override = plan_.quantityOverrides.find(quantity);
      return override == plan_.quantityOverrides.end() ? plan_.global : override->second;
    }();
    if (found->second.StandardError() / std::abs(found->second.mean) > requirement.target) return false;
  }
  return true;
}

long AccuracyController::Histories() const {
  if(resourceHistories_>0) return resourceHistories_;
  long maximum = 0;
  for (const auto& item : estimates_) maximum = std::max(maximum, item.second.count);
  return maximum;
}

bool AccuracyController::ResourceLimitReached() const { return Histories() >= plan_.maxHistories; }

} // namespace shield_sim
